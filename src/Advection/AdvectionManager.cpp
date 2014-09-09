#include "AdvectionManager.h"
#include "ShapeFunction.h"
#include "TracerSkeleton.h"

namespace lasm {

AdvectionManager::AdvectionManager() {
    domain = NULL;
    mesh = NULL;
    timeManager = NULL;
    regrid = NULL;
    cellTree = NULL;
    numMixedTracer = 0;
    numVoidCell = 0;
    filamentLimit = 100;
    strictFilamentLimit = 5;
    radialMixing = 1;
    lateralMixing = 1000;
    strictLateralMixing = 10;
    shrinkFactor = 0.05;
    disorderDegreeLimit = 1.05;
    restartFileName = "N/A";
    isMassFixed = true;
    REPORT_ONLINE;
}

AdvectionManager::~AdvectionManager() {
    if (regrid != NULL) {
        delete regrid;
    }
    if (cellTree != NULL) {
        delete cellTree;
    }
    REPORT_OFFLINE;
}

void AdvectionManager::init(const Domain &domain, const Mesh &mesh,
                            const ConfigManager &configManager,
                            const TimeManager &timeManager) {
    this->domain = &domain;
    this->mesh = &mesh;
    this->timeManager = &timeManager;
    TimeLevelIndex<2> timeIdx;
    assert(timeIdx.get() == 0);
    // get parameters from configuration manager
    if (configManager.hasKey("lasm", "filament_limit")) {
        configManager.getValue("lasm", "filament_limit", filamentLimit);
    }
    if (configManager.hasKey("lasm", "strict_filament_limit")) {
        configManager.getValue("lasm", "strict_filament_limit", strictFilamentLimit);
    }
    if (configManager.hasKey("lasm", "radial_mixing")) {
        configManager.getValue("lasm", "radial_mixing", radialMixing);
    }
    if (configManager.hasKey("lasm", "lateral_mixing")) {
        configManager.getValue("lasm", "lateral_mixing", lateralMixing);
    }
    if (configManager.hasKey("lasm", "strict_lateral_mixing")) {
        configManager.getValue("lasm", "strict_lateral_mixing", strictLateralMixing);
    }
    if (configManager.hasKey("lasm", "disorder_degree_limit")) {
        configManager.getValue("lasm", "disorder_degree_limit", disorderDegreeLimit);
    }
    if (configManager.hasKey("lasm", "shrink_factor")) {
        configManager.getValue("lasm", "shrink_factor", shrinkFactor);
    }
    if (configManager.hasKey("lasm", "is_mass_fixed")) {
        configManager.getValue("lasm", "is_mass_fixed", isMassFixed);
    }
    if (configManager.hasKey("lasm", "restart_file")) {
        configManager.getValue("lasm", "restart_file", restartFileName);
    }
    // initialize tracer manager
    tracerManager.init(domain, mesh, configManager);
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        tracerManager.tracers[t]->actualFilamentLimit = filamentLimit;
    }
    // initialize shape function (or kernel function)
    ShapeFunction::init(domain);
    // initialize regrid object
    if (regrid == NULL) {
        regrid = new Regrid(mesh);
    }
    // initialize mesh adpator
    meshAdaptor.init(domain, mesh);
    // initialize tree structure of mesh grids
#if defined USE_SPHERE_DOMAIN
    cellCoords.reshape(3, mesh.getTotalNumGrid(CENTER, domain.getNumDim()));
#elif defined USE_CARTESIAN_DOMAIN
    cellCoords.reshape(domain.getNumDim(), mesh.getTotalNumGrid(CENTER, domain.getNumDim()));
#endif
    for (int i = 0; i < mesh.getTotalNumGrid(CENTER, domain.getNumDim()); ++i) {
        cellCoords.col(i) = meshAdaptor.getCoord(i).getCartCoord();
    }
    cellTree = new Tree(cellCoords, cellCoordsMap);
    cellCoords = cellTree->Dataset();
    // connect tracers and mesh grids
    embedTracersIntoMesh(timeIdx);
    connectTracersAndMesh(timeIdx);
}

void AdvectionManager::registerTracer(const string &name, const string &units,
                                      const string &brief) {
    tracerManager.registerTracer(name, units, brief);
    meshAdaptor.registerTracer(name, units, brief);
    for (int l = 0; l < 2; ++l) {
        totalMass.getLevel(l).push_back(0);
    }
}

void AdvectionManager::input(const TimeLevelIndex<2> &timeIdx, double *q) {
    // copy the input tracer density onto internal mesh grids
    int l = 0;
    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER, domain->getNumDim()); ++i) {
            meshAdaptor.getDensity(timeIdx, s, i) = q[l];
            meshAdaptor.getMass(timeIdx, s, i) = q[l]*meshAdaptor.getVolume(i);
            l++;
        }
    }
    // transfer the tracer mass from cells to tracers
    remapMeshToTracers(timeIdx);
    diagnose(timeIdx);
    // TODO: The mass on cells could be different after remapping from tracers.
    remapTracersToMesh(timeIdx);
    diagnose(timeIdx);
}

void AdvectionManager::restart(const TimeLevelIndex<2> &timeIdx) {
    
}

void AdvectionManager::output(const TimeLevelIndex<2> &timeIdx, int ncId) {
    // append global attributes
    nc_redef(ncId);
    nc_put_att(ncId, NC_GLOBAL, "filament_limit", NC_DOUBLE, 1, &filamentLimit);
    nc_put_att(ncId, NC_GLOBAL, "strict_filament_limit", NC_DOUBLE, 1, &strictFilamentLimit);
    nc_put_att(ncId, NC_GLOBAL, "radial_mixing", NC_DOUBLE, 1, &radialMixing);
    nc_put_att(ncId, NC_GLOBAL, "lateral_mixing", NC_DOUBLE, 1, &lateralMixing);
    nc_put_att(ncId, NC_GLOBAL, "strict_lateral_mixing", NC_DOUBLE, 1, &strictLateralMixing);
    nc_put_att(ncId, NC_GLOBAL, "disorder_degree_limit", NC_DOUBLE, 1, &disorderDegreeLimit);
    nc_put_att(ncId, NC_GLOBAL, "shrink_factor", NC_DOUBLE, 1, &shrinkFactor);
    nc_enddef(ncId);
    // append tracer data
    tracerManager.output(timeIdx, ncId);
}

void AdvectionManager::diagnose(const TimeLevelIndex<2> &timeIdx) {
    calcTotalMass(timeIdx);
    // print total mass for each species
    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
        REPORT_NOTICE("Total mass of \"" <<
                      tracerManager.getSpeciesInfo(s).getName() << "\" is " <<
                      setw(30) << setprecision(20) <<
                      totalMass.getLevel(timeIdx)[s] << ".");
    }
}

void AdvectionManager::advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                               const VelocityField &velocity) {
    static clock_t time1, time2;

    time1 = clock();
    integrate_RK4(dt, newTimeIdx, velocity);
    time2 = clock();
    REPORT_NOTICE("integrate_RK uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    embedTracersIntoMesh(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("embedTracersIntoMesh uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    connectTracersAndMesh(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("connectTracersAndMesh uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");
    
    time1 = clock();
    remapTracersToMesh(newTimeIdx, &velocity);
    time2 = clock();
    REPORT_NOTICE("remapTracersToMesh uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    checkTracerShapes(newTimeIdx, velocity);
    time2 = clock();
    REPORT_NOTICE("checkTracerShapes uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    mixTracers(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("mixTracers uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");
}

// -----------------------------------------------------------------------------
// private member functions

void AdvectionManager::calcTotalMass(const TimeLevelIndex<2> &timeIdx) {
    // print total mass for each species
    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
        totalMass.getLevel(timeIdx)[s] = 0;
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER, domain->getNumDim()); ++i) {
            totalMass.getLevel(timeIdx)[s] += meshAdaptor.getMass(timeIdx, s, i);
        }
    }
}

/*
    Trajectory equation:

                                        dx
                                        -- = v,
                                        dt

    Deformation equation (H is not computed from this equation):

                                      dH
                                      -- = ∇v H.
                                      dt
 */

void AdvectionManager::integrate_RK4(double dt,
                                     const TimeLevelIndex<2> &newTimeIdx,
                                     const VelocityField &velocity) {
    TimeLevelIndex<2> oldTimeIdx = newTimeIdx-1;
    TimeLevelIndex<2> halfTimeIdx = newTimeIdx-0.5;
    double dt05 = 0.5*dt;
    const auto &divergence = velocity.getDivergence();
#pragma omp parallel for
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        Velocity v1(domain->getNumDim());
        Velocity v2(domain->getNumDim());
        Velocity v3(domain->getNumDim());
        Velocity v4(domain->getNumDim());
        Velocity v(domain->getNumDim());
        double div;
        vec rho(tracerManager.tracers.size());
        double k1_rho[tracerManager.getNumSpecies()];
        double k2_rho[tracerManager.getNumSpecies()];
        double k3_rho[tracerManager.getNumSpecies()];
        double k4_rho[tracerManager.getNumSpecies()];
        // update location and deformation matrix of tracer
        SpaceCoord &x0 = tracer->getX(oldTimeIdx);
        SpaceCoord &x1 = tracer->getX(newTimeIdx);
        MeshIndex &idx0 = tracer->getMeshIndex(oldTimeIdx);
        MeshIndex &idx1 = tracer->getMeshIndex(newTimeIdx);
#ifdef USE_RLL_MESH
        // TODO: Should we hide the following codes? Because they are
        //       related to sphere domain and RLL mesh.
        if (idx0.isOnPole()) {
            idx0.setMoveOnPole(true);
            idx1.setMoveOnPole(true);
            x0.transformToPS(*domain);
        } else {
            idx0.setMoveOnPole(false);
            idx1.setMoveOnPole(false);
        }
#endif
        // stage 1
        regrid->run(BILINEAR, oldTimeIdx, velocity, x0, v1, &idx0);
        regrid->run(BILINEAR, oldTimeIdx, divergence, x0, div, &idx0);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            k1_rho[s] = -tracer->getDensity(s)*div;
            rho[s] = tracer->getDensity(s)+dt05*k1_rho[s];
        }
        mesh->move(x0, dt05, v1, idx0, x1);
        idx1.locate(*mesh, x1);
        // stage 2
        regrid->run(BILINEAR, halfTimeIdx, velocity, x1, v2, &idx1);
        regrid->run(BILINEAR, halfTimeIdx, divergence, x1, div, &idx1);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            k2_rho[s] = -rho[s]*div;
            rho[s] = tracer->getDensity(s)+dt05*k2_rho[s];
        }
        mesh->move(x0, dt05, v2, idx0, x1);
        idx1.locate(*mesh, x1);
        // stage 3
        regrid->run(BILINEAR, halfTimeIdx, velocity, x1, v3, &idx1);
        regrid->run(BILINEAR, newTimeIdx, divergence, x1, div, &idx1);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            k3_rho[s] = -rho[s]*div;
            rho[s] = tracer->getDensity(s)+dt*k3_rho[s];
        }
        mesh->move(x0, dt, v3, idx0, x1);
        idx1.locate(*mesh, x1);
        // stage 4
        regrid->run(BILINEAR, newTimeIdx, velocity, x1, v4, &idx1);
        regrid->run(BILINEAR, newTimeIdx, divergence, x1, div, &idx1);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            k4_rho[s] = -rho[s]*div;
            tracer->getDensity(s) += dt*
                (k1_rho[s]+2.0*k2_rho[s]+2.0*k3_rho[s]+k4_rho[s])/6.0;
        }
        v = (v1+v2*2.0+v3*2.0+v4)/6.0;
        mesh->move(x0, dt, v, idx0, x1);
        idx1.locate(*mesh, x1);
#ifdef USE_SPHERE_DOMAIN
        x1.transformToCart(*domain);
#endif
        // update skeleton points of tracer
        TracerSkeleton &s = tracer->getSkeleton();
        vector<SpaceCoord*> &x0s = s.getSpaceCoords(oldTimeIdx);
        vector<SpaceCoord*> &x1s = s.getSpaceCoords(newTimeIdx);
        vector<MeshIndex*> &idx0s = s.getMeshIdxs(oldTimeIdx);
        vector<MeshIndex*> &idx1s = s.getMeshIdxs(newTimeIdx);
        for (int i = 0; i < x0s.size(); ++i) {
#ifdef USE_RLL_MESH
            if (idx0s[i]->isOnPole()) {
                idx0s[i]->setMoveOnPole(true);
                idx1s[i]->setMoveOnPole(true);
                x0s[i]->transformToPS(*domain);
            } else {
                idx0s[i]->setMoveOnPole(false);
                idx1s[i]->setMoveOnPole(false);
            }
#endif
            // stage 1
            regrid->run(BILINEAR, oldTimeIdx, velocity, *x0s[i], v1, idx0s[i]);
            mesh->move(*x0s[i], dt05, v1, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(*mesh, *x1s[i]);
            // stage 2
            regrid->run(BILINEAR, halfTimeIdx, velocity, *x1s[i], v2, idx1s[i]);
            mesh->move(*x0s[i], dt05, v2, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(*mesh, *x1s[i]);
            // stage 3
            regrid->run(BILINEAR, halfTimeIdx, velocity, *x1s[i], v3, idx1s[i]);
            mesh->move(*x0s[i], dt, v3, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(*mesh, *x1s[i]);
            // stage 4
            regrid->run(BILINEAR, newTimeIdx, velocity, *x1s[i], v4, idx1s[i]);
            v = (v1+v2*2.0+v3*2.0+v4)/6.0;
            mesh->move(*x0s[i], dt, v, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(*mesh, *x1s[i]);
#ifdef USE_SPHERE_DOMAIN
            x1s[i]->transformToCart(*domain);
#endif
        }
        tracer->updateDeformMatrix(*domain, *mesh, newTimeIdx);
    }
}

void AdvectionManager::embedTracersIntoMesh(const TimeLevelIndex<2> &timeIdx) {
    meshAdaptor.resetContainedTracers();
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        int i = tracer->getMeshIndex(timeIdx).getIndex(*mesh, CENTER);
        meshAdaptor.containTracer(i, tracer);
    }
}

void AdvectionManager::connectTracerAndMesh(const TimeLevelIndex<2> &timeIdx,
                                            Tracer *tracer) {
    // call mlpack::range::RangeSearch to find out the neighbor cells of tracers
    // and set the data structures for both cells and tracers for remapping
    Searcher a(cellTree, NULL, cellCoords,
               tracer->getX(timeIdx).getCartCoord(), true);
    double longAxisSize = tracer->getShapeSize(timeIdx)[0];
    mlpack::math::Range r(0.0, longAxisSize);
    vector<vector<size_t> > neighbors;
    vector<vector<double> > distances;
    a.Search(r, neighbors, distances);
    BodyCoord y(domain->getNumDim());
    for (int i = 0; i < neighbors[0].size(); ++i) {
        int j = cellCoordsMap[neighbors[0][i]];
        // calculate the tracer shape function for the cell
        tracer->getBodyCoord(*domain, timeIdx, meshAdaptor.getCoord(j), y);
        double w = tracer->getShapeFunction(timeIdx, y);
        if (w > 0.0) {
#pragma omp critical
            meshAdaptor.connectTracer(j, tracer, w);
            tracer->connectCell(j);
        }
    }
    if (tracer->getNumConnectedCell() == 0) {
        // tracer has not connect with any cells, so connect with its host cell
#pragma omp critical
        meshAdaptor.connectTracer(tracer->getHostCell(), tracer,
                                  ShapeFunction::getMaxValue());
        tracer->connectCell(tracer->getHostCell());
    }
}

void AdvectionManager::connectTracersAndMesh(const TimeLevelIndex<2> &timeIdx) {
    meshAdaptor.resetConnectedTracers();
#pragma omp parallel for
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        tracer->resetConnectedCells();
        connectTracerAndMesh(timeIdx, tracer);
    }
}

void AdvectionManager::checkTracerShapes(const TimeLevelIndex<2> &timeIdx,
                                         const VelocityField &velocity) {
    numMixedTracer = 0;
    const double cosThetaBound = cos(20*RAD);

    double tmp1 = -1.0e15, tmp2 = 1.0e15;

//#define PRINT_DISORDER_TRACERS
#ifdef PRINT_DISORDER_TRACERS
    vector<int> disorder;
#endif
    // calculate tracer actual filament limit
#pragma omp parallel for
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        const vec &S = tracer->getS();
        double filament = S[0]/S[1];
        if (filament < strictFilamentLimit ||
            tracer->actualFilamentLimit == strictFilamentLimit) continue;
        const vector<int> &connectedCells = tracer->getConnectedCells();
        if (connectedCells.size() <= 1) {
            CHECK_POINT;
            continue;
        }
        vector<Tracer*> tracers;
        for (int i = 0; i < tracer->getNumConnectedCell(); ++i) {
            for (int j = 0; j < meshAdaptor.getNumConnectedTracer(connectedCells[i]); ++j) {
                Tracer *tracer1 = meshAdaptor.getConnectedTracers(connectedCells[i])[j];
                bool alreadyAdded = false;
                if (tracer1 != tracer) {
                    for (int k = 0; k < tracers.size(); ++k) {
                        if (tracer1 == tracers[k]) {
                            alreadyAdded = true;
                            break;
                        }
                    }
                    if (!alreadyAdded) {
                        tracers.push_back(tracer1);
                    }
                }
            }
        }
        if (tracers.size() > 1) {
            vec x0(2);
#ifdef USE_SPHERE_DOMAIN
            domain->project(geomtk::SphereDomain::STEREOGRAPHIC,
                            tracer->getX(timeIdx),
                            tracer->getLongAxisVertexSpaceCoord(), x0);
#else
            x0 = tracer->getLongAxisVertexSpaceCoord()()-tracer->getX(timeIdx)();
#endif
            vec x1(2), x2(2);
            vec cosThetas(tracers.size());
            for (int i = 0; i < tracers.size(); ++i) {
                Tracer *tracer1 = tracers[i];
                if (tracer1 == tracer) continue;
#ifdef USE_SPHERE_DOMAIN
                domain->project(geomtk::SphereDomain::STEREOGRAPHIC,
                                tracer->getX(timeIdx),
                                tracer1->getX(timeIdx), x1);
                domain->project(geomtk::SphereDomain::STEREOGRAPHIC,
                                tracer->getX(timeIdx),
                                tracer1->getLongAxisVertexSpaceCoord(), x2);
#else
                x1 = tracer1->getX(timeIdx)()-tracer->getX(timeIdx)();
                x2 = tracer1->getLongAxisVertexSpaceCoord()()-tracer->getX(timeIdx)();
#endif
                cosThetas[i] = fabs(norm_dot(x0, x2-x1));
            }
            double disorderDegree = mean(cosThetas)/(min(cosThetas)+1.0e-15);
            if (tmp1 < disorderDegree) tmp1 = disorderDegree;
            if (tmp2 > disorderDegree) tmp2 = disorderDegree;
            if (disorderDegree > disorderDegreeLimit &&
                min(cosThetas) < cosThetaBound) {
                std::ofstream file("disorder_tracers.txt");
                tracer->dump(timeIdx, *domain, meshAdaptor, file, 0);
                int idx = 1;
                for (int i = 0; i < tracers.size(); ++i) {
                    tracers[i]->dump(timeIdx, *domain, meshAdaptor, file, idx++);
                }
                file << "vertices = new((/" << tracers.size() << ",2/), double)" << endl;
                for (int m = 0; m < domain->getNumDim(); ++m) {
                    file << "vertices(:," << m << ") = (/";
                    for (int i = 0; i < tracers.size(); ++i) {
                        file << tracers[i]->getLongAxisVertexSpaceCoord()(m)/RAD;
                        if (i != tracers.size()-1) {
                            file << ",";
                        } else {
                            file << "/)" << endl;
                        }
                    }
                }
                file.close();
                for (int i = 0; i < tracers.size(); ++i) {
                    tracers[i]->actualFilamentLimit = strictFilamentLimit;
                    tracers[i]->actualLateralMixing = strictLateralMixing;
                }
                tracer->actualFilamentLimit = strictFilamentLimit;
                tracer->actualLateralMixing = strictLateralMixing;
            }
        }
    }
    // check the degree of filamentation
#pragma omp parallel for
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
#ifdef PRINT_DISORDER_TRACERS
        if (tracer->actualFilamentLimit < filamentLimit) {
#pragma omp critical
            disorder.push_back(tracer->getID());
        }
#endif
        const vec &S = tracer->getS();
        double filament = S[0]/S[1];
        if (filament > tracer->actualFilamentLimit) {
#pragma omp critical
            recordTracer(Tracer::NEED_MIXING, tracer);
        } else {
            tracer->actualFilamentLimit = filamentLimit;
            tracer->actualLateralMixing = lateralMixing;
        }
    }
#ifdef PRINT_DISORDER_TRACERS
    for (int i = 0; i < disorder.size(); ++i) {
        if (i != disorder.size()-1) {
            cout << disorder[i] << ",";
        } else {
            cout << disorder[i] << endl;
        }
    }
#endif
    // cout << setw(20) << setprecision(10) << tmp1;
    // cout << setw(20) << setprecision(10) << tmp2 << endl;
}

void AdvectionManager::mixTracers(const TimeLevelIndex<2> &timeIdx) {
    for (int t = 0; t < numMixedTracer; ++t) {
        Tracer *tracer = mixedTracers[t];
        // get surrounding tracers
        int cell = tracer->getHostCell();
        const vector<Tracer*> &surroundTracers = meshAdaptor.getConnectedTracers(cell);
#ifndef NDEBUG
//#define CHECK_MIX_TRACER
#ifdef CHECK_MIX_TRACER
        std::ofstream file("mixed_tracers.txt");
        int idx = 1;
        tracer->dump(timeIdx, *domain, file, 0);
#endif
        assert(surroundTracers.size() > 1);
        vec totalMass1(tracerManager.getNumSpecies(), arma::fill::zeros);
        vec totalMass2(tracerManager.getNumSpecies(), arma::fill::zeros);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            totalMass1[s] += tracer->getMass(s);
        }
#endif
        // calcuate mixing weights
        vec x0(2);
#ifdef USE_SPHERE_DOMAIN
        domain->project(geomtk::SphereDomain::STEREOGRAPHIC,
                        tracer->getX(timeIdx),
                        tracer->getLongAxisVertexSpaceCoord(), x0);
#else
        x0 = tracer->getLongAxisVertexSpaceCoord()()-tracer->getX(timeIdx)();
#endif
        double n0 = norm(x0, 2);
        vec x1(2);
        vec weights(meshAdaptor.getNumConnectedTracer(cell), arma::fill::zeros);
        for (int i = 0; i < meshAdaptor.getNumConnectedTracer(cell); ++i) {
            Tracer *tracer1 = surroundTracers[i];
            if (tracer1 == tracer) continue;
#ifdef USE_SPHERE_DOMAIN
            domain->project(geomtk::SphereDomain::STEREOGRAPHIC,
                            tracer->getX(timeIdx),
                            tracer1->getX(timeIdx), x1);
#else
            x1 = tracer1->getX(timeIdx)()-tracer->getX(timeIdx)();
#endif
            x1 /= n0;
            double cosTheta = norm_dot(x0, x1);
            double sinTheta = sqrt(1-cosTheta*cosTheta);
            double n1 = norm(x1, 2);
            double d1 = n1*cosTheta;
            double d2 = n1*sinTheta;
            weights[i] = exp(-(radialMixing*d1*d1+tracer->actualLateralMixing*d2*d2));
            if (weights[i] < 1.0e-6) {
                weights[i] = 0;
            }
#ifndef NDEBUG
            if (weights[i] != 0) {
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    totalMass1[s] += tracer1->getMass(s);
                }
            }
            assert(weights[i] <= 1);
#ifdef CHECK_MIX_TRACER
            tracer1->dump(timeIdx, *domain, file, idx++);
#endif
#endif
        }
        double sumWeights = sum(weights);
        if (sumWeights < 1.0e-14) {
            continue;
        }
        weights /= sumWeights;
#ifndef NDEBUG
#ifdef CHECK_MIX_TRACER
        file << "weights = (/";
        for (int i = 0; i < meshAdaptor.getNumConnectedTracer(cell); ++i) {
            if (i == meshAdaptor.getNumConnectedTracer(cell)-1) {
                file << weights[i] << "/)" << endl;
            } else {
                file << weights[i] << ",";
            }
        }
#endif
#endif
        // distribute mass to surrounding tracers
        vec &S = tracer->getS();
        double oldVolume = tracer->getDetH(timeIdx);
        double newVolume = (1-shrinkFactor)*oldVolume;
        for (int i = 0; i < meshAdaptor.getNumConnectedTracer(cell); ++i) {
            Tracer *tracer1 = surroundTracers[i];
            if (weights[i] == 0) continue;
            double volume1 = tracer1->getDetH(timeIdx);
            double volume2 = shrinkFactor*oldVolume*weights[i]+volume1;
            assert(volume2 > volume1);
            tracer1->getS() *= pow(volume2/volume1, 1.0/domain->getNumDim());
            tracer1->updateDeformMatrix(*domain, timeIdx);
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                double m1 = tracer->getMass(s)*shrinkFactor;
                double &m2 = tracer1->getMass(s);
                m2 += m1*weights[i];
#ifndef NDEBUG
                totalMass2[s] += m2;
                double rho1 = tracer->getDensity(s);
                double rho2 = tracer1->getDensity(s);
#endif
                tracer1->calcDensity(timeIdx, s);
#ifndef NDEBUG
                double rho3 = tracer1->getDensity(s);
                assert(((rho1 <= rho3 || (rho1-rho3)/rho1 <= 1.0e-12) &&
                        (rho3 <= rho2 || (rho3-rho2)/rho2 <= 1.0e-12)) ||
                       ((rho2 <= rho3 || (rho2-rho3)/rho2 <= 1.0e-12) &&
                        (rho3 <= rho1 || (rho3-rho1)/rho1 <= 1.0e-12)));
#endif
            }
        }
        // change tracer shape and mass
        S[0] = newVolume/S[1];
        tracer->updateDeformMatrix(*domain, timeIdx);
        tracer->resetSkeleton(*domain, *mesh, timeIdx);
        tracer->setType(Tracer::GOOD_SHAPE);
        tracer->actualFilamentLimit = filamentLimit;
        tracer->actualLateralMixing = lateralMixing;
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            tracer->calcMass(timeIdx, s);
#ifndef NDEBUG
            totalMass2[s] += tracer->getMass(s);
#endif
        }
#ifndef NDEBUG
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            assert(fabs(totalMass1[s]-totalMass2[s])/totalMass1[s] < 1.0e-12);
        }
#ifdef CHECK_MIX_TRACER
        tracer->dump(timeIdx, *domain, file, idx);
        file.close();
#endif
#endif
    }
}
    
void AdvectionManager::fillVoidCells(const TimeLevelIndex<2> &timeIdx) {
    // NOTE: The occurance of void cells should be rare.
#pragma omp parallel for
    for (int c = 0; c < numVoidCell; ++c) {
        int cell = voidCells[c];
        Searcher a(cellTree, NULL, cellCoords,
                   meshAdaptor.getCoord(cell).getCartCoord(), true);
#ifdef USE_SPHERE_DOMAIN
        double searchRadius = 1*RAD*domain->getRadius();
#else
        double searchRadius = 0.1*domain->getAxisSpan(0);
#endif
        while (true) {
            mlpack::math::Range r(0.0, searchRadius);
            vector<vector<size_t> > neighbors;
            vector<vector<double> > distances;
            a.Search(r, neighbors, distances);
            if (neighbors[0].size() != 0) {
                vector<int> ngbCells;
                for (int i = 0; i < neighbors[0].size(); ++i) {
                    int ngbCell = cellCoordsMap[neighbors[0][i]];
                    // check if the cell is not a void one
                    if (find(voidCells.begin(), voidCells.end(), ngbCell)
                        == voidCells.end()) {
                        ngbCells.push_back(ngbCell);
                    }
                }
                if (ngbCells.size() == 0) {
                    searchRadius *= 2;
                    continue;
                }
                vec weights(ngbCells.size());
                for (int i = 0; i < ngbCells.size(); ++i) {
                    double d = domain->calcDistance(meshAdaptor.getCoord(cell),
                                                    meshAdaptor.getCoord(ngbCells[i]));
                    weights[i] = 1/d;
                }
                weights /= sum(weights);
                for (int i = 0; i < ngbCells.size(); ++i) {
                    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                        meshAdaptor.getDensity(timeIdx, s, cell) += meshAdaptor.getDensity(timeIdx, s, ngbCells[i])*weights[i];
                    }
                }
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    meshAdaptor.getMass(timeIdx, s, cell) = meshAdaptor.getDensity(timeIdx, s, cell)*meshAdaptor.getVolume(cell);
                }
                break;
            } else {
                searchRadius *= 2;
            }
        }
    }
    
}

#define REMAP_DENSITY
//#define REMAP_MASS

void AdvectionManager::remapMeshToTracers(const TimeLevelIndex<2> &timeIdx) {
#pragma omp parallel for
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        tracer->resetSpecies();
        const vector<int> &cells = tracer->getConnectedCells();
        double totalWeight = 0;
#if defined REMAP_DENSITY
        for (int i = 0; i < tracer->getNumConnectedCell(); ++i) {
            double weight = meshAdaptor.getRemapWeight(cells[i], tracer);
            totalWeight += weight;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                tracer->getDensity(s) += meshAdaptor.getDensity(timeIdx, s, cells[i])*weight;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            tracer->getDensity(s) /= totalWeight;
            tracer->calcMass(timeIdx, s);
        }
#elif defined REMAP_MASS
        for (int i = 0; i < tracer->getNumConnectedCell(); ++i) {
            double weight = cells[i]->getRemapWeight(tracer);
            totalWeight += weight;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                tracer->getSpeciesMass(s) += meshAdaptor.getMass(timeIdx, s, cells[i])*weight;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            tracer->getSpeciesMass(s) /= totalWeight;
            tracer->getSpeciesDensity(s) = tracer->getSpeciesMass(s)/tracer->getDetH(timeIdx);
        }
#endif
    }
}

void AdvectionManager::remapTracersToMesh(const TimeLevelIndex<2> &timeIdx,
                                          const VelocityField *velocity) {
    numVoidCell = 0;
    meshAdaptor.resetSpecies(timeIdx);
#pragma omp parallel for
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, domain->getNumDim()); ++i) {
        if (meshAdaptor.getNumConnectedTracer(i) == 0) {
            if (numVoidCell == voidCells.size()) {
                voidCells.push_back(i);
            } else {
                voidCells[numVoidCell] = i;
            }
            numVoidCell++;
            continue;
        };
        const vector<Tracer*> &tracers = meshAdaptor.getConnectedTracers(i);
        double W = 0;
#if defined REMAP_DENSITY
        for (int j = 0; j < meshAdaptor.getNumConnectedTracer(i); ++j) {
            double w = meshAdaptor.getRemapWeight(i, tracers[j]);
            W += w;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                meshAdaptor.getDensity(timeIdx, s, i) += tracers[j]->getDensity(s)*w;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            meshAdaptor.getDensity(timeIdx, s, i) /= W;
            meshAdaptor.getMass(timeIdx, s, i) = meshAdaptor.getDensity(timeIdx, s, i)*meshAdaptor.getVolume(i);
        }
#elif defined REMAP_MASS
        for (int j = 0; j < cell.getNumConnectedTracer(); ++j) {
            double w = cell.getRemapWeight(tracers[j]);
            W += w;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                meshAdaptor.getMass(timeIdx, s, i) += tracers[j]->getSpeciesMass(s)*w;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            meshAdaptor.getMass(timeIdx, s, i) /= W;
            meshAdaptor.getDensity(timeIdx, s, i) = meshAdaptor.getMass(timeIdx, s, i)/meshAdaptor.getVolume(i);
        }
#endif
    }
    fillVoidCells(timeIdx);
#ifdef REMAP_DENSITY
    if (isMassFixed) {
        correctTotalMassOnMesh(timeIdx);
    }
#endif
}

void AdvectionManager::correctTotalMassOnMesh(const TimeLevelIndex<2> &timeIdx) {
    double expectedTotalMass[tracerManager.getNumSpecies()];
    if (timeIdx.isCurrentIndex()) {
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            expectedTotalMass[s] = totalMass.getLevel(timeIdx)[s];
        }
    } else {
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            expectedTotalMass[s] = totalMass.getLevel(timeIdx-1)[s];
        }
    }
    calcTotalMass(timeIdx);
    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
        double fixer = expectedTotalMass[s]/totalMass.getLevel(timeIdx)[s];
        double biasPercent = (totalMass.getLevel(timeIdx)[s]-
                              expectedTotalMass[s])/expectedTotalMass[s];
        REPORT_NOTICE("Mass conservation bias percentage is " <<
                      std::fixed << setw(10) << setprecision(4) <<
                      biasPercent*100 << "%.");
#pragma omp parallel for
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER, domain->getNumDim()); ++i) {
            meshAdaptor.getDensity(timeIdx, s, i) *= fixer;
            meshAdaptor.getMass(timeIdx, s, i) =
                meshAdaptor.getDensity(timeIdx, s, i)*meshAdaptor.getVolume(i);
        }
    }
}

void AdvectionManager::recordTracer(Tracer::TracerType type, Tracer *tracer) {
    tracer->setType(type);
    switch (type) {
        case Tracer::NEED_MIXING:
            if (numMixedTracer == mixedTracers.size()) {
                mixedTracers.push_back(tracer);
            } else {
                mixedTracers[numMixedTracer] = tracer;
            }
            numMixedTracer++;
            break;
        default:
            REPORT_ERROR("Unexpected branch!");
            break;
    }
}

}
