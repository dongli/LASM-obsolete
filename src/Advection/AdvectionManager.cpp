#include "AdvectionManager.h"
#include "ShapeFunction.h"
#include "TracerSkeleton.h"
#include "AdvectionTestCase.h"

namespace lasm {

AdvectionManager::AdvectionManager() {
    domain = NULL;
    mesh = NULL;
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
                            const ConfigManager &configManager) {
    this->domain = &domain;
    this->mesh = &mesh;
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
    cellCoords.reshape(3, mesh.totalNumGrid(CENTER, domain.numDim()));
#elif defined USE_CARTESIAN_DOMAIN
    cellCoords.reshape(domain.numDim(), mesh.totalNumGrid(CENTER, domain.numDim()));
#endif
    for (int i = 0; i < mesh.totalNumGrid(CENTER, domain.numDim()); ++i) {
        cellCoords.col(i) = meshAdaptor.coord(i).cartCoord();
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
        totalMass.level(l).push_back(0);
    }
}

void AdvectionManager::input(const TimeLevelIndex<2> &timeIdx, double *q) {
    // copy the input tracer density onto internal mesh grids
    int l = 0;
    for (int s = 0; s < tracerManager.numSpecies(); ++s) {
        for (int i = 0; i < mesh->totalNumGrid(CENTER, domain->numDim()); ++i) {
            meshAdaptor.density(timeIdx, s, i) = q[l];
            meshAdaptor.mass(timeIdx, s, i) = q[l]*meshAdaptor.volume(i);
            l++;
        }
    }
    // transfer the tracer mass from cells to tracers
    remapMeshToTracers(timeIdx);
    diagnose(timeIdx);
    // TODO: The mass on cells could be different after remapping from tracers.
    remapTracersToMesh(timeIdx);
}

void AdvectionManager::restart(const string &fileName) {
    TimeLevelIndex<2> timeIdx;
    tracerManager.input(fileName);
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        tracerManager.tracers[t]->actualFilamentLimit = filamentLimit;
    }
    // connect tracers and mesh grids
    embedTracersIntoMesh(timeIdx);
    connectTracersAndMesh(timeIdx);
    meshAdaptor.input(fileName, tracerManager);
    diagnose(timeIdx);
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
    for (int s = 0; s < tracerManager.numSpecies(); ++s) {
        REPORT_NOTICE("Total mass of \"" <<
                      tracerManager.speciesInfo(s).name() << "\" is " <<
                      setw(30) << setprecision(20) <<
                      totalMass.level(timeIdx)[s] << ".");
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
    remapTracersToMesh(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("remapTracersToMesh uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    checkTracerShapes(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("checkTracerShapes uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    mixTracers(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("mixTracers uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");
}

void AdvectionManager::advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                               const AdvectionTestCase &testCase) {
    integrate_RK4(dt, newTimeIdx, testCase);
    embedTracersIntoMesh(newTimeIdx);
    connectTracersAndMesh(newTimeIdx);
    remapTracersToMesh(newTimeIdx);
//    checkTracerShapes(newTimeIdx);
//    mixTracers(newTimeIdx);
}

// -----------------------------------------------------------------------------
// private member functions

void AdvectionManager::calcTotalMass(const TimeLevelIndex<2> &timeIdx) {
    // print total mass for each species
    for (int s = 0; s < tracerManager.numSpecies(); ++s) {
        totalMass.level(timeIdx)[s] = 0;
        for (int i = 0; i < mesh->totalNumGrid(CENTER, domain->numDim()); ++i) {
            totalMass.level(timeIdx)[s] += meshAdaptor.mass(timeIdx, s, i);
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
    const auto &divergence = velocity.divergence();
#pragma omp parallel for
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        Velocity v1(domain->numDim());
        Velocity v2(domain->numDim());
        Velocity v3(domain->numDim());
        Velocity v4(domain->numDim());
        Velocity v(domain->numDim());
        double div;
        vec rho(tracerManager.tracers.size());
        double k1_rho[tracerManager.numSpecies()];
        double k2_rho[tracerManager.numSpecies()];
        double k3_rho[tracerManager.numSpecies()];
        double k4_rho[tracerManager.numSpecies()];
        // Update the centroid and deformation matrix of the tracer.
        SpaceCoord &x0 = tracer->x(oldTimeIdx);
        SpaceCoord &x1 = tracer->x(newTimeIdx);
        MeshIndex &idx0 = tracer->meshIndex(oldTimeIdx);
        MeshIndex &idx1 = tracer->meshIndex(newTimeIdx);
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
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            k1_rho[s] = -tracer->density(s)*div;
            rho[s] = tracer->density(s)+dt05*k1_rho[s];
        }
        mesh->move(x0, dt05, v1, idx0, x1);
        idx1.locate(*mesh, x1);
        // stage 2
        regrid->run(BILINEAR, halfTimeIdx, velocity, x1, v2, &idx1);
        regrid->run(BILINEAR, halfTimeIdx, divergence, x1, div, &idx1);
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            k2_rho[s] = -rho[s]*div;
            rho[s] = tracer->density(s)+dt05*k2_rho[s];
        }
        mesh->move(x0, dt05, v2, idx0, x1);
        idx1.locate(*mesh, x1);
        // stage 3
        regrid->run(BILINEAR, halfTimeIdx, velocity, x1, v3, &idx1);
        regrid->run(BILINEAR, halfTimeIdx, divergence, x1, div, &idx1);
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            k3_rho[s] = -rho[s]*div;
            rho[s] = tracer->density(s)+dt*k3_rho[s];
        }
        mesh->move(x0, dt, v3, idx0, x1);
        idx1.locate(*mesh, x1);
        // stage 4
        regrid->run(BILINEAR, newTimeIdx, velocity, x1, v4, &idx1);
        regrid->run(BILINEAR, newTimeIdx, divergence, x1, div, &idx1);
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            k4_rho[s] = -rho[s]*div;
            tracer->density(s) += dt*
                (k1_rho[s]+2.0*k2_rho[s]+2.0*k3_rho[s]+k4_rho[s])/6.0;
        }
        v = (v1+v2*2.0+v3*2.0+v4)/6.0;
        mesh->move(x0, dt, v, idx0, x1);
        idx1.locate(*mesh, x1);
#ifdef USE_SPHERE_DOMAIN
        x1.transformToCart(*domain);
#endif
        // Update the skeleton points of the tracer.
        TracerSkeleton &s = tracer->skeleton();
        vector<SpaceCoord*> &x0s = s.spaceCoords(oldTimeIdx);
        vector<SpaceCoord*> &x1s = s.spaceCoords(newTimeIdx);
        vector<MeshIndex*> &idx0s = s.meshIndices(oldTimeIdx);
        vector<MeshIndex*> &idx1s = s.meshIndices(newTimeIdx);
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

void AdvectionManager::integrate_RK4(double dt,
                                     const TimeLevelIndex<2> &newTimeIdx,
                                     const AdvectionTestCase &testCase) {
    TimeLevelIndex<2> oldTimeIdx = newTimeIdx-1;
    double dt05 = 0.5*dt;
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        Velocity v1(domain->numDim());
        Velocity v2(domain->numDim());
        Velocity v3(domain->numDim());
        Velocity v4(domain->numDim());
        Velocity v(domain->numDim());
        double div;
        vec rho(tracerManager.tracers.size());
        double k1_rho[tracerManager.numSpecies()];
        double k2_rho[tracerManager.numSpecies()];
        double k3_rho[tracerManager.numSpecies()];
        double k4_rho[tracerManager.numSpecies()];
        // Update the centroid and deformation matrix of the tracer.
        SpaceCoord &x0 = tracer->x(oldTimeIdx);
        SpaceCoord &x1 = tracer->x(newTimeIdx);
        MeshIndex &idx0 = tracer->meshIndex(oldTimeIdx);
        MeshIndex &idx1 = tracer->meshIndex(newTimeIdx);
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
        testCase.evalVelocity(0, x0, idx0.isMoveOnPole(), v1);
        testCase.evalDivergence(0, x0, div);
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            k1_rho[s] = -tracer->density(s)*div;
            rho[s] = tracer->density(s)+dt05*k1_rho[s];
        }
        mesh->move(x0, dt05, v1, idx0, x1);
        idx1.locate(*mesh, x1);
        // stage 2
        testCase.evalVelocity(dt05, x1, idx0.isMoveOnPole(), v2);
        testCase.evalDivergence(dt05, x1, div);
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            k2_rho[s] = -rho[s]*div;
            rho[s] = tracer->density(s)+dt05*k2_rho[s];
        }
        mesh->move(x0, dt05, v2, idx0, x1);
        idx1.locate(*mesh, x1);
        // stage 3
        testCase.evalVelocity(dt05, x1, idx0.isMoveOnPole(), v3);
        testCase.evalDivergence(dt05, x1, div);
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            k3_rho[s] = -rho[s]*div;
            rho[s] = tracer->density(s)+dt*k3_rho[s];
        }
        mesh->move(x0, dt, v3, idx0, x1);
        idx1.locate(*mesh, x1);
        // stage 4
        testCase.evalVelocity(dt, x1, idx0.isMoveOnPole(), v4);
        testCase.evalDivergence(dt, x1, div);
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            k4_rho[s] = -rho[s]*div;
            tracer->density(s) += dt*
            (k1_rho[s]+2.0*k2_rho[s]+2.0*k3_rho[s]+k4_rho[s])/6.0;
        }
        v = (v1+v2*2.0+v3*2.0+v4)/6.0;
        mesh->move(x0, dt, v, idx0, x1);
        idx1.locate(*mesh, x1);
#ifdef USE_SPHERE_DOMAIN
        x1.transformToCart(*domain);
#endif
        // Update the skeleton points of the tracer.
        TracerSkeleton &s = tracer->skeleton();
        vector<SpaceCoord*> &x0s = s.spaceCoords(oldTimeIdx);
        vector<SpaceCoord*> &x1s = s.spaceCoords(newTimeIdx);
        vector<MeshIndex*> &idx0s = s.meshIndices(oldTimeIdx);
        vector<MeshIndex*> &idx1s = s.meshIndices(newTimeIdx);
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
            testCase.evalVelocity(0, *x0s[i], idx0s[i]->isMoveOnPole(), v1);
            mesh->move(*x0s[i], dt05, v1, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(*mesh, *x1s[i]);
            // stage 2
            testCase.evalVelocity(dt05, *x1s[i], idx0s[i]->isMoveOnPole(), v2);
            mesh->move(*x0s[i], dt05, v2, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(*mesh, *x1s[i]);
            // stage 3
            testCase.evalVelocity(dt05, *x1s[i], idx0s[i]->isMoveOnPole(), v3);
            mesh->move(*x0s[i], dt, v3, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(*mesh, *x1s[i]);
            // stage 4
            testCase.evalVelocity(dt, *x1s[i], idx0s[i]->isMoveOnPole(), v4);
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
        int i = tracer->meshIndex(timeIdx).getIndex(*mesh, CENTER);
        meshAdaptor.containTracer(i, tracer);
    }
}

void AdvectionManager::connectTracerAndMesh(const TimeLevelIndex<2> &timeIdx,
                                            Tracer *tracer) {
    // call mlpack::range::RangeSearch to find out the neighbor cells of tracers
    // and set the data structures for both cells and tracers for remapping
    Searcher a(cellTree, NULL, cellCoords,
               tracer->x(timeIdx).cartCoord(), true);
    double longAxisSize = tracer->shapeSize(timeIdx)[0];
    mlpack::math::Range r(0.0, longAxisSize);
    vector<vector<size_t> > neighbors;
    vector<vector<double> > distances;
    a.Search(r, neighbors, distances);
    BodyCoord y(domain->numDim());
    for (int i = 0; i < neighbors[0].size(); ++i) {
        int j = cellCoordsMap[neighbors[0][i]];
        // calculate the tracer shape function for the cell
        tracer->calcBodyCoord(*domain, timeIdx, meshAdaptor.coord(j), y);
        double w = tracer->shapeFunction(timeIdx, y);
        if (w > 0.0) {
#pragma omp critical
            {
                meshAdaptor.connectTracer(j, tracer, w);
                tracer->connectCell(j);
            }
        }
    }
    if (tracer->numConnectedCell() == 0) {
        // tracer has not connect with any cells, so connect with its host cell
#pragma omp critical
        {
            meshAdaptor.connectTracer(tracer->hostCellIndex(), tracer,
                                      ShapeFunction::maxValue());
            tracer->connectCell(tracer->hostCellIndex());
        }
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

void AdvectionManager::checkTracerShapes(const TimeLevelIndex<2> &timeIdx) {
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
        const vec &S = tracer->S();
        double filament = S[0]/S[1];
        if (filament < strictFilamentLimit ||
            tracer->actualFilamentLimit == strictFilamentLimit) continue;
        const vector<int> &connectedCells = tracer->connectedCellIndices();
        if (connectedCells.size() <= 1 ||
            tracer->detH(timeIdx) < 1.5*meshAdaptor.volume(tracer->hostCellIndex())) {
            continue;
        }
        vector<Tracer*> tracers;
        for (int i = 0; i < tracer->numConnectedCell(); ++i) {
            for (int j = 0; j < meshAdaptor.numConnectedTracer(connectedCells[i]); ++j) {
                Tracer *tracer1 = meshAdaptor.connectedTracers(connectedCells[i])[j];
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
                            tracer->x(timeIdx),
                            tracer->longAxisVertexSpaceCoord(), x0);
#else
            x0 = tracer->longAxisVertexSpaceCoord()()-tracer->x(timeIdx)();
#endif
            vec x1(2), x2(2);
            vec cosThetas(tracers.size());
            for (int i = 0; i < tracers.size(); ++i) {
                Tracer *tracer1 = tracers[i];
                if (tracer1 == tracer) continue;
#ifdef USE_SPHERE_DOMAIN
                domain->project(geomtk::SphereDomain::STEREOGRAPHIC,
                                tracer->x(timeIdx),
                                tracer1->x(timeIdx), x1);
                domain->project(geomtk::SphereDomain::STEREOGRAPHIC,
                                tracer->x(timeIdx),
                                tracer1->longAxisVertexSpaceCoord(), x2);
#else
                x1 = tracer1->x(timeIdx)()-tracer->x(timeIdx)();
                x2 = tracer1->longAxisVertexSpaceCoord()()-tracer->x(timeIdx)();
#endif
                cosThetas[i] = fabs(norm_dot(x0, x2-x1));
            }
            double disorderDegree = mean(cosThetas)/(min(cosThetas)+1.0e-15);
            if (tmp1 < disorderDegree) tmp1 = disorderDegree;
            if (tmp2 > disorderDegree) tmp2 = disorderDegree;
            if (disorderDegree > disorderDegreeLimit &&
                min(cosThetas) < cosThetaBound) {
                // std::ofstream file("disorder_tracers.txt");
                // tracer->dump(timeIdx, *domain, meshAdaptor, file, 0);
                // int idx = 1;
                // for (int i = 0; i < tracers.size(); ++i) {
                //     tracers[i]->dump(timeIdx, *domain, meshAdaptor, file, idx++);
                // }
                // file << "vertices = new((/" << tracers.size() << ",2/), double)" << endl;
                // for (int m = 0; m < domain->numDim(); ++m) {
                //     file << "vertices(:," << m << ") = (/";
                //     for (int i = 0; i < tracers.size(); ++i) {
                //         file << tracers[i]->longAxisVertexSpaceCoord()(m)/RAD;
                //         if (i != tracers.size()-1) {
                //             file << ",";
                //         } else {
                //             file << "/)" << endl;
                //         }
                //     }
                // }
                // file.close();
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
            disorder.push_back(tracer->ID());
        }
#endif
        const vec &S = tracer->S();
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
        int cell = tracer->hostCellIndex();
        const vector<Tracer*> &surroundTracers = meshAdaptor.connectedTracers(cell);
#ifndef NDEBUG
//#define CHECK_MIX_TRACER
#ifdef CHECK_MIX_TRACER
        std::ofstream file("mixed_tracers.txt");
        int idx = 1;
        tracer->dump(timeIdx, *domain, file, 0);
#endif
        assert(surroundTracers.size() > 1);
        vec totalMass1(tracerManager.numSpecies(), arma::fill::zeros);
        vec totalMass2(tracerManager.numSpecies(), arma::fill::zeros);
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            totalMass1[s] += tracer->mass(s);
        }
#endif
        // calcuate mixing weights
        vec x0(2);
#ifdef USE_SPHERE_DOMAIN
        domain->project(geomtk::SphereDomain::STEREOGRAPHIC,
                        tracer->x(timeIdx),
                        tracer->longAxisVertexSpaceCoord(), x0);
#else
        x0 = tracer->longAxisVertexSpaceCoord()()-tracer->x(timeIdx)();
#endif
        double n0 = norm(x0, 2);
        vec x1(2);
        vec weights(meshAdaptor.numConnectedTracer(cell), arma::fill::zeros);
        for (int i = 0; i < meshAdaptor.numConnectedTracer(cell); ++i) {
            Tracer *tracer1 = surroundTracers[i];
            if (tracer1 == tracer) continue;
#ifdef USE_SPHERE_DOMAIN
            domain->project(geomtk::SphereDomain::STEREOGRAPHIC,
                            tracer->x(timeIdx),
                            tracer1->x(timeIdx), x1);
#else
            x1 = tracer1->x(timeIdx)()-tracer->x(timeIdx)();
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
                for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                    totalMass1[s] += tracer1->mass(s);
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
        for (int i = 0; i < meshAdaptor.numConnectedTracer(cell); ++i) {
            if (i == meshAdaptor.numConnectedTracer(cell)-1) {
                file << weights[i] << "/)" << endl;
            } else {
                file << weights[i] << ",";
            }
        }
#endif
#endif
        // distribute mass to surrounding tracers
        vec &S = tracer->S();
        double oldVolume = tracer->detH(timeIdx);
        double newVolume = (1-shrinkFactor)*oldVolume;
        for (int i = 0; i < meshAdaptor.numConnectedTracer(cell); ++i) {
            Tracer *tracer1 = surroundTracers[i];
            if (weights[i] == 0) continue;
            double volume1 = tracer1->detH(timeIdx);
            double volume2 = shrinkFactor*oldVolume*weights[i]+volume1;
            assert(volume2 > volume1);
            tracer1->S() *= pow(volume2/volume1, 1.0/domain->numDim());
            tracer1->updateDeformMatrix(*domain, timeIdx);
            for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                double m1 = tracer->mass(s)*shrinkFactor;
                double &m2 = tracer1->mass(s);
                m2 += m1*weights[i];
#ifndef NDEBUG
                totalMass2[s] += m2;
                double rho1 = tracer->density(s);
                double rho2 = tracer1->density(s);
#endif
                tracer1->calcDensity(timeIdx, s);
#ifndef NDEBUG
                double rho3 = tracer1->density(s);
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
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            tracer->calcMass(timeIdx, s);
#ifndef NDEBUG
            totalMass2[s] += tracer->mass(s);
#endif
        }
#ifndef NDEBUG
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
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
                   meshAdaptor.coord(cell).cartCoord(), true);
#ifdef USE_SPHERE_DOMAIN
        double searchRadius = 1*RAD*domain->radius();
#else
        double searchRadius = 0.1*domain->axisSpan(0);
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
                    double d = domain->calcDistance(meshAdaptor.coord(cell),
                                                    meshAdaptor.coord(ngbCells[i]));
                    weights[i] = 1/d;
                }
                weights /= sum(weights);
                for (int i = 0; i < ngbCells.size(); ++i) {
                    for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                        meshAdaptor.density(timeIdx, s, cell) += meshAdaptor.density(timeIdx, s, ngbCells[i])*weights[i];
                    }
                }
                for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                    meshAdaptor.mass(timeIdx, s, cell) = meshAdaptor.density(timeIdx, s, cell)*meshAdaptor.volume(cell);
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
        const vector<int> &cells = tracer->connectedCellIndices();
        double totalWeight = 0;
#if defined REMAP_DENSITY
        for (int i = 0; i < tracer->numConnectedCell(); ++i) {
            double weight = meshAdaptor.remapWeight(cells[i], tracer);
            totalWeight += weight;
            for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                tracer->density(s) += meshAdaptor.density(timeIdx, s, cells[i])*weight;
            }
        }
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            tracer->density(s) /= totalWeight;
            tracer->calcMass(timeIdx, s);
        }
#elif defined REMAP_MASS
        for (int i = 0; i < tracer->numConnectedCell(); ++i) {
            double weight = cells[i]->remapWeight(tracer);
            totalWeight += weight;
            for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                tracer->mass(s) += meshAdaptor.mass(timeIdx, s, cells[i])*weight;
            }
        }
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            tracer->mass(s) /= totalWeight;
            tracer->density(s) = tracer->mass(s)/tracer->detH(timeIdx);
        }
#endif
    }
}

void AdvectionManager::remapTracersToMesh(const TimeLevelIndex<2> &timeIdx) {
    numVoidCell = 0;
    meshAdaptor.resetSpecies(timeIdx);
#pragma omp parallel for
    for (int i = 0; i < mesh->totalNumGrid(CENTER); ++i) {
        if (meshAdaptor.numConnectedTracer(i) == 0) {
            if (numVoidCell == voidCells.size()) {
                voidCells.push_back(i);
            } else {
                voidCells[numVoidCell] = i;
            }
            numVoidCell++;
            continue;
        };
        const vector<Tracer*> &tracers = meshAdaptor.connectedTracers(i);
        double totalWeight = 0;
#if defined REMAP_DENSITY
        for (int j = 0; j < meshAdaptor.numConnectedTracer(i); ++j) {
            double weight = meshAdaptor.remapWeight(i, tracers[j]);
            totalWeight += weight;
            for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                meshAdaptor.density(timeIdx, s, i) += tracers[j]->density(s)*weight;
            }
        }
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            meshAdaptor.density(timeIdx, s, i) /= totalWeight;
            meshAdaptor.mass(timeIdx, s, i) = meshAdaptor.density(timeIdx, s, i)*meshAdaptor.volume(i);
        }
#elif defined REMAP_MASS
        for (int j = 0; j < cell.numConnectedTracer(); ++j) {
            double weight = cell.remapWeight(tracers[j]);
            totalWeight += w;
            for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                meshAdaptor.mass(timeIdx, s, i) += tracers[j]->mass(s)*weight;
            }
        }
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            meshAdaptor.mass(timeIdx, s, i) /= totalWeight;
            meshAdaptor.density(timeIdx, s, i) = meshAdaptor.mass(timeIdx, s, i)/meshAdaptor.volume(i);
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
    double expectedTotalMass[tracerManager.numSpecies()];
    if (timeIdx.isCurrentIndex()) {
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            expectedTotalMass[s] = totalMass.level(timeIdx)[s];
        }
    } else {
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            expectedTotalMass[s] = totalMass.level(timeIdx-1)[s];
        }
    }
    calcTotalMass(timeIdx);
    for (int s = 0; s < tracerManager.numSpecies(); ++s) {
        double fixer = expectedTotalMass[s]/totalMass.level(timeIdx)[s];
        double biasPercent = (totalMass.level(timeIdx)[s]-
                              expectedTotalMass[s])/expectedTotalMass[s];
        REPORT_NOTICE("Mass conservation bias percentage is " <<
                      std::fixed << setw(10) << setprecision(4) <<
                      biasPercent*100 << "%.");
#pragma omp parallel for
        for (int i = 0; i < mesh->totalNumGrid(CENTER, domain->numDim()); ++i) {
            meshAdaptor.mass(timeIdx, s, i) *= fixer;
            meshAdaptor.density(timeIdx, s, i) =
                meshAdaptor.mass(timeIdx, s, i)/meshAdaptor.volume(i);
        }
    }
    diagnose(timeIdx); // NOTE: Do not delete this line!
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
