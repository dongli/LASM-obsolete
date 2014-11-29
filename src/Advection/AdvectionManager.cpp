#include "AdvectionManager.h"
#include "ShapeFunction.h"
#include "TracerSkeleton.h"
#include "AdvectionTestCase.h"

namespace lasm {

AdvectionManager::AdvectionManager() {
    domain = NULL;
    mesh = NULL;
    regrid = NULL;
    filter = NULL;
    cellTree = NULL;
    numVoidCell = 0;
    filamentLimit = 100;
    radialMixing = 1;
    lateralMixing = 1000;
    restoreFactor = 0.1;
    isMassFixed = true;
    REPORT_ONLINE;
}

AdvectionManager::~AdvectionManager() {
    if (regrid != NULL) {
        delete regrid;
    }
    if (filter != NULL) {
        delete filter;
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
    // Get parameters from configuration manager.
    if (configManager.hasKey("lasm", "filament_limit")) {
        configManager.getValue("lasm", "filament_limit", filamentLimit);
    }
    if (configManager.hasKey("lasm", "radial_mixing")) {
        configManager.getValue("lasm", "radial_mixing", radialMixing);
    }
    if (configManager.hasKey("lasm", "lateral_mixing")) {
        configManager.getValue("lasm", "lateral_mixing", lateralMixing);
    }
    if (configManager.hasKey("lasm", "restore_factor")) {
        configManager.getValue("lasm", "restore_factor", restoreFactor);
    }
    if (configManager.hasKey("lasm", "is_mass_fixed")) {
        configManager.getValue("lasm", "is_mass_fixed", isMassFixed);
    }
    // Initialize tracer manager.
    tracerManager.init(domain, mesh, configManager);
    // Initialize shape function (or kernel function).
    ShapeFunction::init(domain);
    // Initialize regrid object.
    if (regrid == NULL) {
        regrid = new Regrid(mesh);
    }
    // Initialize filter object.
    if (filter == NULL) {
        filter = new Filter(mesh, geomtk::NINE_POINT_SMOOTHING);
    }
    // Initialize mesh adpator.
    meshAdaptor.init(domain, mesh);
    // Initialize tree structure of mesh grids.
#if defined LASM_SPHERE_DOMAIN
    cellCoords.reshape(3, mesh.totalNumGrid(CENTER, domain.numDim()));
#elif defined LASM_CARTESIAN_DOMAIN
    cellCoords.reshape(domain.numDim(), mesh.totalNumGrid(CENTER, domain.numDim()));
#endif
    for (int i = 0; i < mesh.totalNumGrid(CENTER, domain.numDim()); ++i) {
        cellCoords.col(i) = meshAdaptor.coord(i).cartCoord();
    }
    cellTree = new Tree(cellCoords, cellCoordsMap);
    cellCoords = cellTree->Dataset();
    // Connect tracers and mesh grids.
    embedTracersIntoMesh(timeIdx);
    connectTracersAndMesh(timeIdx);
}

void AdvectionManager::registerTracer(const string &name, const string &units,
                                      const string &brief, bool smooth) {
    tracerManager.registerTracer(name, units, brief, smooth);
    meshAdaptor.registerTracer(name, units, brief);
    for (int l = 0; l < 2; ++l) {
        totalMass.level(l).push_back(0);
    }
}

void AdvectionManager::input(const TimeLevelIndex<2> &timeIdx, double *q) {
    // Copy input tracer density onto internal mesh grids.
    int l = 0;
    for (int s = 0; s < tracerManager.numSpecies(); ++s) {
        for (int i = 0; i < mesh->totalNumGrid(CENTER, domain->numDim()); ++i) {
            meshAdaptor.density(timeIdx, s, i) = q[l];
            meshAdaptor.mass(timeIdx, s, i) = q[l]*meshAdaptor.volume(i);
            l++;
        }
    }
    // Transfer tracer mass from cells to tracers.
    remapMeshToTracers(timeIdx);
    diagnose(timeIdx);
    // NOTE: The mass on cells could be different after remapping from tracers.
    remapTracersToMesh(timeIdx);
}

void AdvectionManager::restart(const string &fileName) {
    TimeLevelIndex<2> timeIdx;
    tracerManager.input(fileName);
    // Connect tracers and mesh grids.
    embedTracersIntoMesh(timeIdx);
    connectTracersAndMesh(timeIdx);
    meshAdaptor.input(fileName, tracerManager);
    diagnose(timeIdx);
}

void AdvectionManager::output(const TimeLevelIndex<2> &timeIdx, int ncId) {
    // Append global attributes.
    nc_redef(ncId);
    nc_put_att(ncId, NC_GLOBAL, "filament_limit", NC_DOUBLE, 1, &filamentLimit);
    nc_put_att(ncId, NC_GLOBAL, "radial_mixing", NC_DOUBLE, 1, &radialMixing);
    nc_put_att(ncId, NC_GLOBAL, "lateral_mixing", NC_DOUBLE, 1, &lateralMixing);
    nc_put_att(ncId, NC_GLOBAL, "restore_factor", NC_DOUBLE, 1, &restoreFactor);
    nc_enddef(ncId);
    // Append tracer data.
    tracerManager.output(timeIdx, ncId);
}

void AdvectionManager::diagnose(const TimeLevelIndex<2> &timeIdx) {
    calcTotalMass(timeIdx);
    // Print total mass for each species.
    for (int s = 0; s < tracerManager.numSpecies(); ++s) {
        REPORT_NOTICE("Total mass of \"" <<
                      tracerManager.speciesInfo(s).name() << "\" is " <<
                      setw(30) << setprecision(20) <<
                      totalMass.level(timeIdx)[s] << ".");
    }
}

void AdvectionManager::advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                               const VelocityField &velocity) {
//#if defined LASM_TENDENCY_ON_MESH
//    PRINT_USED_TIME(remapTendencyFromMesh(newTimeIdx-1));
//#endif
    PRINT_USED_TIME(integrate_RK4(dt, newTimeIdx, velocity));
    PRINT_USED_TIME(embedTracersIntoMesh(newTimeIdx));
    PRINT_USED_TIME(connectTracersAndMesh(newTimeIdx));
    PRINT_USED_TIME(remapTracersToMesh(newTimeIdx));
    PRINT_USED_TIME(mixTracers(newTimeIdx));
}

void AdvectionManager::advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                               const AdvectionTestCase &testCase) {
    PRINT_USED_TIME(integrate_RK4(dt, newTimeIdx, testCase));
    PRINT_USED_TIME(embedTracersIntoMesh(newTimeIdx));
    PRINT_USED_TIME(connectTracersAndMesh(newTimeIdx));
    PRINT_USED_TIME(remapTracersToMesh(newTimeIdx));
    PRINT_USED_TIME(mixTracers(newTimeIdx));
}

void AdvectionManager::calcTotalMass(const TimeLevelIndex<2> &timeIdx) {
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
        // Add tracer density tendencies recorded.
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            tracer->density(s) += tracer->tendency(s)*dt;
            tracer->mass(s) = tracer->density(s)*tracer->detH(oldTimeIdx);
        }
        Velocity v1(domain->numDim());
        Velocity v2(domain->numDim());
        Velocity v3(domain->numDim());
        Velocity v4(domain->numDim());
        Velocity v(domain->numDim());
        double div;
        double rho[tracerManager.numSpecies()];
        double k1_rho[tracerManager.numSpecies()];
        double k2_rho[tracerManager.numSpecies()];
        double k3_rho[tracerManager.numSpecies()];
        double k4_rho[tracerManager.numSpecies()];
        // Update the centroid and deformation matrix of the tracer.
        SpaceCoord &x0 = tracer->x(oldTimeIdx);
        SpaceCoord &x1 = tracer->x(newTimeIdx);
        MeshIndex &idx0 = tracer->meshIndex(oldTimeIdx);
        MeshIndex &idx1 = tracer->meshIndex(newTimeIdx);
#ifdef LASM_RLL_MESH
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
#ifdef LASM_SPHERE_DOMAIN
        x1.transformToCart(*domain);
#endif
        // Update the skeleton points of the tracer.
        TracerSkeleton &s = tracer->skeleton();
        vector<SpaceCoord*> &x0s = s.spaceCoords(oldTimeIdx);
        vector<SpaceCoord*> &x1s = s.spaceCoords(newTimeIdx);
        vector<MeshIndex*> &idx0s = s.meshIndices(oldTimeIdx);
        vector<MeshIndex*> &idx1s = s.meshIndices(newTimeIdx);
        for (int i = 0; i < x0s.size(); ++i) {
#ifdef LASM_RLL_MESH
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
#ifdef LASM_SPHERE_DOMAIN
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
#pragma omp parallel for
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        // Add tracer density tendencies recorded.
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            tracer->density(s) += tracer->tendency(s)*dt;
            tracer->mass(s) = tracer->density(s)*tracer->detH(oldTimeIdx);
        }
        Velocity v1(domain->numDim());
        Velocity v2(domain->numDim());
        Velocity v3(domain->numDim());
        Velocity v4(domain->numDim());
        Velocity v(domain->numDim());
        double div;
        double rho[tracerManager.numSpecies()];
        double k1_rho[tracerManager.numSpecies()];
        double k2_rho[tracerManager.numSpecies()];
        double k3_rho[tracerManager.numSpecies()];
        double k4_rho[tracerManager.numSpecies()];
        // Update the centroid and deformation matrix of the tracer.
        SpaceCoord &x0 = tracer->x(oldTimeIdx);
        SpaceCoord &x1 = tracer->x(newTimeIdx);
        MeshIndex &idx0 = tracer->meshIndex(oldTimeIdx);
        MeshIndex &idx1 = tracer->meshIndex(newTimeIdx);
#ifdef LASM_RLL_MESH
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
#ifdef LASM_SPHERE_DOMAIN
        x1.transformToCart(*domain);
#endif
        // Update the skeleton points of the tracer.
        TracerSkeleton &s = tracer->skeleton();
        vector<SpaceCoord*> &x0s = s.spaceCoords(oldTimeIdx);
        vector<SpaceCoord*> &x1s = s.spaceCoords(newTimeIdx);
        vector<MeshIndex*> &idx0s = s.meshIndices(oldTimeIdx);
        vector<MeshIndex*> &idx1s = s.meshIndices(newTimeIdx);
        for (int i = 0; i < x0s.size(); ++i) {
#ifdef LASM_RLL_MESH
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
#ifdef LASM_SPHERE_DOMAIN
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
//#pragma omp critical
            meshAdaptor.connectTracer(j, tracer, w);
            tracer->connectCell(j);
        }
    }
    if (tracer->numConnectedCell() == 0) {
        // tracer has not connect with any cells, so connect with its host cell
//#pragma omp critical
        meshAdaptor.connectTracer(tracer->hostCellIndex(), tracer,
                                  ShapeFunction::maxValue());
        tracer->connectCell(tracer->hostCellIndex());
    }
}

void AdvectionManager::connectTracersAndMesh(const TimeLevelIndex<2> &timeIdx) {
    meshAdaptor.resetConnectedTracers();
//#pragma omp parallel for
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        tracer->resetConnectedCells();
        connectTracerAndMesh(timeIdx, tracer);
    }
}

void AdvectionManager::mixTracers(const TimeLevelIndex<2> &timeIdx) {
    int numMixedTracer = 0;
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        TracerSkeleton &s = tracer->skeleton();
        // Check the validity of the linear assumption.
        vector<SpaceCoord*> &x = s.spaceCoords(timeIdx);
        const vector<BodyCoord*> &y = s.bodyCoords();
        SpaceCoord X(2);
        double bias = -1.0e+33;
        for (int i = 0; i < 4; ++i) {
            tracer->calcSpaceCoord(*domain, timeIdx, *(y[i]), X);
            double bias0 = domain->calcDistance(*(x[i]), X)/tracer->shapeSize(timeIdx)[0];
            if (bias0 > bias) bias = bias0;
        }
        // Set the bias limit based on the filament of the tracer and its volume
        // compared with the neighbor tracers.
        const double minBiasLimit = 0.1;
        const double maxBiasLimit = 0.5;
        vector<Tracer*> neighborTracers = getNeighborTracers(tracer);
        double meanVolume = tracer->detH(timeIdx);
        for (int i = 0; i < neighborTracers.size(); ++i) {
            meanVolume += neighborTracers[i]->detH(timeIdx);
        }
        meanVolume /= neighborTracers.size()+1;
        double ratio = tracer->detH(timeIdx)/meanVolume;
        double biasLimit = transitionFunction(1, maxBiasLimit, 5, minBiasLimit,
                                              ratio*tracer->filament());
        // Check tracer bias.
        tracer->linearDegeneration() = bias;
        if (bias < biasLimit && tracer->filament() < filamentLimit) {
            continue;
        }
        // Calcuate the mixing weights.
        vec x0(2);
#ifdef LASM_SPHERE_DOMAIN
        domain->project(geomtk::SphereDomain::STEREOGRAPHIC, tracer->x(timeIdx),
                        tracer->longAxisVertexSpaceCoord(), x0);
#else
        x0 = tracer->longAxisVertexSpaceCoord()()-tracer->x(timeIdx)();
#endif
        double n0 = norm(x0, 2);
        vec x1(2);
        vec weights(neighborTracers.size(), arma::fill::zeros);
        for (int i = 0; i < neighborTracers.size(); ++i) {
#ifdef LASM_SPHERE_DOMAIN
            domain->project(geomtk::SphereDomain::STEREOGRAPHIC, tracer->x(timeIdx),
                            neighborTracers[i]->x(timeIdx), x1);
#else
            x1 = neighborTracers[i]->x(timeIdx)()-tracer->x(timeIdx)();
#endif
            x1 /= n0;
            double cosTheta = norm_dot(x0, x1);
            double sinTheta = sqrt(1-cosTheta*cosTheta);
            double n1 = norm(x1, 2);
            double d1 = n1*cosTheta;
            double d2 = n1*sinTheta;
            weights[i] = exp(-(radialMixing*d1*d1+lateralMixing*d2*d2));
            if (weights[i] < 1.0e-6) {
                weights[i] = 0;
            }
        }
        double sumWeights = sum(weights);
        if (sumWeights < 1.0e-14) {
            continue;
        }
        weights /= sumWeights;
        // Caclulate total mass.
        double totalMass[tracerManager.numSpecies()];
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            totalMass[s] = tracer->mass(s);
            for (int i = 0; i < neighborTracers.size(); ++i) {
                if (weights[i] == 0) continue;
                totalMass[s] += neighborTracers[i]->mass(s);
            }
        }
        // Caclulate weighted total volume.
        double weightedTotalVolume = tracer->detH(timeIdx);
        for (int i = 0; i < neighborTracers.size(); ++i) {
            if (weights[i] == 0) continue;
            weightedTotalVolume += neighborTracers[i]->detH(timeIdx)*weights[i];
        }
        // Calculate weighted mean tracer density.
        double weightedMeanDensity[tracerManager.numSpecies()];
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            weightedMeanDensity[s] = tracer->mass(s);
            for (int i = 0; i < neighborTracers.size(); ++i) {
                if (weights[i] == 0) continue;
                weightedMeanDensity[s] += neighborTracers[i]->mass(s)*weights[i];
            }
            weightedMeanDensity[s] /= weightedTotalVolume;
        }
        // Restore tracer density to mean density.
        double newTotalMass[tracerManager.numSpecies()];
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            tracer->density(s) += restoreFactor*(weightedMeanDensity[s]-tracer->density(s));
            tracer->mass(s) = tracer->density(s)*tracer->detH(timeIdx);
            newTotalMass[s] = tracer->mass(s);
        }
        for (int i = 0; i < neighborTracers.size(); ++i) {
            if (weights[i] == 0) continue;
            double c = restoreFactor*weights[i];
            for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                neighborTracers[i]->density(s) += c*(weightedMeanDensity[s]-neighborTracers[i]->density(s));
                neighborTracers[i]->mass(s) = neighborTracers[i]->density(s)*neighborTracers[i]->detH(timeIdx);
                newTotalMass[s] += neighborTracers[i]->mass(s);
            }
        }
        // Fix mass inconservation due to floating point inaccuracy.
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            if (totalMass[s] == 0) {
                if (newTotalMass[s] != 0) {
                    REPORT_ERROR("totalMass[" << s << "] is zero, but " <<
                                 "newTotalMass[" << s << "] is not!");
                }
                continue;
            }
            double fixer = totalMass[s]/newTotalMass[s];
            if (fixer == 1) continue;
            tracer->mass(s) *= fixer;
            tracer->density(s) = tracer->mass(s)/tracer->detH(timeIdx);
            for (int i = 0; i < neighborTracers.size(); ++i) {
                neighborTracers[i]->mass(s) *= fixer;
                neighborTracers[i]->density(s) = neighborTracers[i]->mass(s)/neighborTracers[i]->detH(timeIdx);
            }
        }
        // Change problematic tracer shape (make tracer more uniform).
        const double maxUniformFactor = 0.5;
        double a, b;
        double uniformFactor = transitionFunction(1, 1, 5, maxUniformFactor,
                                                  tracer->filament());
        a = pow(pow(uniformFactor, domain->numDim()-1), 1.0/domain->numDim());
        b = 1/a;
        tracer->S()[0] *= a;
        for (int i = 1; i < domain->numDim(); ++i) {
            tracer->S()[i] *= b;
        }
        tracer->updateDeformMatrix(*domain, timeIdx);
        tracer->resetSkeleton(*domain, *mesh, timeIdx);
        numMixedTracer++;
//        cout << "+++++++++++++++++++++++++++++++++++++++++" << endl;
//        cout << setw(8) << tracer->ID();
//        cout << setw(20) << setprecision(5) << tracer->filament();
//        cout << setw(20) << setprecision(5) << uniformFactor;
//        cout << setw(20) << setprecision(5) << a;
//        cout << setw(20) << setprecision(5) << b;
//        cout << endl;
    }
    if (numMixedTracer != 0) {
        REPORT_NOTICE(numMixedTracer << " tracers are mixed.");
    }
}
    
void AdvectionManager::fillVoidCells(const TimeLevelIndex<2> &timeIdx) {
    // NOTE: The occurance of void cells should be rare.
//#pragma omp parallel for
    for (int c = 0; c < numVoidCell; ++c) {
        int cell = voidCells[c];
        Searcher a(cellTree, NULL, cellCoords,
                   meshAdaptor.coord(cell).cartCoord(), true);
#ifdef LASM_SPHERE_DOMAIN
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

#if defined REMAP_MASS
void AdvectionManager::remapMeshToTracers(const TimeLevelIndex<2> &timeIdx) {
    tracerManager.resetSpecies();
    for (int i = 0; i < mesh->totalNumGrid(CENTER); ++i) {
        const vector<Tracer*> &tracers = meshAdaptor.connectedTracers(i);
        double totalWeight = 0;
        for (int j = 0; j < meshAdaptor.numConnectedTracer(i); ++j) {
            totalWeight += meshAdaptor.remapWeight(i, tracers[j])*tracers[j]->detH(timeIdx);
        }
        for (int j = 0; j < meshAdaptor.numConnectedTracer(i); ++j) {
            double weight = meshAdaptor.remapWeight(i, tracers[j])*tracers[j]->detH(timeIdx)/totalWeight;
            for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                tracers[j]->mass(s) += meshAdaptor.mass(timeIdx, s, i)*weight;
            }
        }
    }
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            tracer->density(s) = tracer->mass(s)/tracer->detH(timeIdx);
        }
    }
}

void AdvectionManager::remapTracersToMesh(const TimeLevelIndex<2> &timeIdx) {
    meshAdaptor.resetSpecies(timeIdx);
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        const vector<int> &cellIndices = tracer->connectedCellIndices();
        double totalWeight = 0;
        for (int i = 0; i < tracer->numConnectedCell(); ++i) {
            int j = cellIndices[i];
            totalWeight += meshAdaptor.remapWeight(j, tracer)*meshAdaptor.volume(j);
        }
        for (int i = 0; i < tracer->numConnectedCell(); ++i) {
            int j = cellIndices[i];
            double weight = meshAdaptor.remapWeight(j, tracer)*meshAdaptor.volume(j)/totalWeight;
            for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                meshAdaptor.mass(timeIdx, s, j) += tracer->mass(s)*weight;
            }
        }
    }
    for (int i = 0; i < mesh->totalNumGrid(CENTER); ++i) {
        if (meshAdaptor.numConnectedTracer(i) == 0) {
            REPORT_ERROR("Encounter void cell! Handle it.");
        }
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            meshAdaptor.density(timeIdx, s, i) = meshAdaptor.mass(timeIdx, s, i)/meshAdaptor.volume(i);
        }
    }
}
#elif defined REMAP_DENSITY
void AdvectionManager::remapTendencyFromMesh(const TimeLevelIndex<2> &timeIdx) {
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            tracer->tendency(s) = 0;
        }
        const vector<int> &cellIndices = tracer->connectedCellIndices();
        double totalWeight = 0;
        for (int i = 0; i < tracer->numConnectedCell(); ++i) {
            int j = cellIndices[i];
            totalWeight += meshAdaptor.remapWeight(j, tracer);
        }
        for (int i = 0; i < tracer->numConnectedCell(); ++i) {
            int j = cellIndices[i];
            double weight = meshAdaptor.remapWeight(j, tracer)/totalWeight;
            for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                tracer->tendency(s) += meshAdaptor.tendency(s, j)*weight;
            }
        }
    }
}

void AdvectionManager::remapMeshToTracers(const TimeLevelIndex<2> &timeIdx) {
    tracerManager.resetSpecies();
    for (int t = 0; t < tracerManager.tracers.size(); ++t) {
        Tracer *tracer = tracerManager.tracers[t];
        const vector<int> &cellIndices = tracer->connectedCellIndices();
        double totalWeight = 0;
        for (int i = 0; i < tracer->numConnectedCell(); ++i) {
            int j = cellIndices[i];
            totalWeight += meshAdaptor.remapWeight(j, tracer);
        }
        for (int i = 0; i < tracer->numConnectedCell(); ++i) {
            int j = cellIndices[i];
            double weight = meshAdaptor.remapWeight(j, tracer)/totalWeight;
            for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                tracer->density(s) += meshAdaptor.density(timeIdx, s, j)*weight;
            }
        }
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            tracer->mass(s) = tracer->density(s)*tracer->detH(timeIdx);
        }
    }
}

void AdvectionManager::remapTracersToMesh(const TimeLevelIndex<2> &timeIdx) {
    numVoidCell = 0;
    meshAdaptor.resetSpecies(timeIdx);
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
        for (int j = 0; j < meshAdaptor.numConnectedTracer(i); ++j) {
            totalWeight += meshAdaptor.remapWeight(i, tracers[j]);
        }
        for (int j = 0; j < meshAdaptor.numConnectedTracer(i); ++j) {
            double weight = meshAdaptor.remapWeight(i, tracers[j])/totalWeight;
            for (int s = 0; s < tracerManager.numSpecies(); ++s) {
                meshAdaptor.density(timeIdx, s, i) += tracers[j]->density(s)*weight;
            }
        }
        for (int s = 0; s < tracerManager.numSpecies(); ++s) {
            meshAdaptor.mass(timeIdx, s, i) = meshAdaptor.density(timeIdx, s, i)*meshAdaptor.volume(i);
        }
    }
    fillVoidCells(timeIdx);
    for (int s = 0; s < tracerManager.numSpecies(); ++s) {
        if (!tracerManager.speciesInfos[s]->smooth()) continue;
        meshAdaptor.density()[s]->applyBndCond(timeIdx);
        filter->run(timeIdx, *meshAdaptor.density()[s]);
        for (int i = 0; i < mesh->totalNumGrid(CENTER); ++i) {
            meshAdaptor.mass(timeIdx, s, i) = meshAdaptor.density(timeIdx, s, i)*meshAdaptor.volume(i);
        }
    }
    if (isMassFixed) {
        correctTotalMassOnMesh(timeIdx);
    }
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
//#pragma omp parallel for
        for (int i = 0; i < mesh->totalNumGrid(CENTER, domain->numDim()); ++i) {
            meshAdaptor.mass(timeIdx, s, i) *= fixer;
            meshAdaptor.density(timeIdx, s, i) =
                meshAdaptor.mass(timeIdx, s, i)/meshAdaptor.volume(i);
        }
    }
    diagnose(timeIdx); // NOTE: Do not delete this line!
}
#endif

vector<Tracer*> AdvectionManager::getNeighborTracers(Tracer *tracer) const {
    const vector<int> &connectedCellIndices = tracer->connectedCellIndices();
    vector<Tracer*> neighborTracers;
    for (int i = 0; i < tracer->numConnectedCell(); ++i) {
        for (int j = 0; j < meshAdaptor.numConnectedTracer(connectedCellIndices[i]); ++j) {
            Tracer *tracer1 = meshAdaptor.connectedTracers(connectedCellIndices[i])[j];
            bool alreadyAdded = false;
            if (tracer1 != tracer) {
                for (int k = 0; k < neighborTracers.size(); ++k) {
                    if (tracer1 == neighborTracers[k]) {
                        alreadyAdded = true;
                        break;
                    }
                }
                if (!alreadyAdded) {
                    neighborTracers.push_back(tracer1);
                }
            }
        }
    }
    return neighborTracers;
}

}
