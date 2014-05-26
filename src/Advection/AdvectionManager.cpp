#include "AdvectionManager.h"
#include "ShapeFunction.h"
#include "TracerSkeleton.h"

namespace lasm {

AdvectionManager::AdvectionManager() {
    domain = NULL;
    mesh = NULL;
    regrid = NULL;
    numLongTracer = 0;
    numUnresolvedTracer = 0;
    numVoidCell = 0;
    alpha = 0.85;
    beta1 = 1;
    beta2 = 10;
    isMassFixed = true;
    REPORT_ONLINE;
}

AdvectionManager::~AdvectionManager() {
    if (regrid != NULL) {
        delete regrid;
    }
    delete cellTree;
    REPORT_OFFLINE;
}

void AdvectionManager::init(const Domain &domain, const Mesh &mesh,
                            const geomtk::ConfigManager &configManager) {
    this->domain = &domain;
    this->mesh = &mesh;
    TimeLevelIndex<2> initTimeIdx;
#ifndef NDEBUG
    assert(initTimeIdx.get() == 0);
#endif
    if (configManager.hasKey("lasm", "alpha")) {
        configManager.getValue("lasm", "alpha", alpha);
    }
    if (configManager.hasKey("lasm", "beta1")) {
        configManager.getValue("lasm", "beta1", beta1);
    }
    if (configManager.hasKey("lasm", "beta2")) {
        configManager.getValue("lasm", "beta2", beta2);
    }
    if (configManager.hasKey("lasm", "is_mass_fixed")) {
        configManager.getValue("lasm", "is_mass_fixed", isMassFixed);
    }
    // -------------------------------------------------------------------------
    // initialize tracer manager
    tracerManager.init(domain, mesh, configManager);
    // initialize shape function (or kernel function)
    ShapeFunction::init(domain);
    // initialize regrid object
    if (regrid == NULL) {
        regrid = new Regrid(mesh);
    }
    // initialize tracer mesh cells
    tracerMeshCells.create("", "", "mesh cells for storing tracers", mesh, CENTER);
    for (int i = 0; i < mesh.getTotalNumGrid(CENTER); ++i) {
        SpaceCoord x(domain.getNumDim());
        mesh.getGridCoord(i, CENTER, x);
        x.transformToCart(domain);
        double volume = mesh.getCellVolume(i);
        for (int l = 0; l < 2; ++l) {
            tracerMeshCells(initTimeIdx+l, i).setCoord(x);
            tracerMeshCells(initTimeIdx+l, i).setVolume(volume);
            tracerMeshCells(initTimeIdx+l, i).setID(i);
        }
    }
    // initialize tree structure of mesh grids
    cellCoords.reshape(3, mesh.getTotalNumGrid(CENTER));
    for (int i = 0; i < mesh.getTotalNumGrid(CENTER); ++i) {
        cellCoords.col(i) = tracerMeshCells(initTimeIdx, i).getCoord().getCartCoord();
    }
    cellTree = new Tree(cellCoords, cellCoordsMap);
    cellCoords = cellTree->Dataset();
    // -------------------------------------------------------------------------
    connectTracersAndMesh(initTimeIdx);
}

void AdvectionManager::registerTracer(const string &name, const string &units,
                                      const string &brief) {
    tracerManager.registerTracer(name, units, brief);
    TimeLevelIndex<2> initTimeIdx;
    for (int l = 0; l < 2; ++l) {
        totalMass.getLevel(l).push_back(0);
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            tracerMeshCells(initTimeIdx+l, i).addSpecies();
        }
    }
}

void AdvectionManager::input(const TimeLevelIndex<2> &timeIdx,
                             vector<ScalarField*> &q) {
#ifndef NDEBUG
    assert(q.size() == tracerManager.getNumSpecies());
#endif
    // -------------------------------------------------------------------------
    // copy the input tracer density onto internal mesh grids
    for (int s = 0; s < q.size(); ++s) {
#ifndef NDEBUG
        assert(q[s]->getMesh().getTotalNumGrid(CENTER) ==
               tracerMeshCells.getMesh().getTotalNumGrid(CENTER));
#endif
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            TracerMeshCell &cell = tracerMeshCells(timeIdx, i);
            cell.getSpeciesDensity(s) = (*q[s])(timeIdx, i);
            cell.getSpeciesMass(s) = (*q[s])(timeIdx, i)*cell.getVolume();
        }
        calcTotalMass(timeIdx);
    }
    // -------------------------------------------------------------------------
    // transfer the tracer mass from cells to tracers
    remapMeshToTracers(timeIdx);
    diagnose(timeIdx);
    // TODO: The mass on cells could be different after remapping from tracers.
    remapTracersToMesh(timeIdx);
    diagnose(timeIdx);
}

void AdvectionManager::output(const string &fileName,
                              const TimeLevelIndex<2> &oldTimeIdx) {
    tracerManager.output(fileName, oldTimeIdx);
    // output the tracer density on the mesh
    int ncId, lonDimId, latDimId;
    int lonVarId, latVarId;
    int dimIds[domain->getNumDim()];
    int qVarIds[tracerManager.getNumSpecies()], volVarId;
    vec lon, lat;
    
    if (nc_open(fileName.c_str(), NC_WRITE, &ncId) != NC_NOERR) {
        REPORT_ERROR("Failed to open " << fileName << "!");
    }
    nc_redef(ncId);
    // -------------------------------------------------------------------------
    if (nc_put_att(ncId, NC_GLOBAL, "alpha", NC_DOUBLE, 1, &alpha) != NC_NOERR){
        REPORT_ERROR("Failed to put global attribute \"alpha\"!");
    }
    if (nc_put_att(ncId, NC_GLOBAL, "beta1", NC_DOUBLE, 1, &beta1) != NC_NOERR){
        REPORT_ERROR("Failed to put global attribute \"beta1\"!");
    }
    if (nc_put_att(ncId, NC_GLOBAL, "beta2", NC_DOUBLE, 1, &beta2) != NC_NOERR){
        REPORT_ERROR("Failed to put global attribute \"beta2\"!");
    }
    // -------------------------------------------------------------------------
    // define dimensions
    // =========================================================================
    // longitude dimension
    if (nc_def_dim(ncId, "lon", mesh->getNumGrid(0, FULL), &lonDimId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define dimension lon!");
    }
    if (nc_def_var(ncId, "lon", NC_DOUBLE, 1, &lonDimId, &lonVarId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define coordinate variable lon!");
    }
    nc_put_att(ncId, lonVarId, "long_name", NC_CHAR, 9, "longitude");
    nc_put_att(ncId, lonVarId, "units", NC_CHAR, 12, "degrees_east");
    // =========================================================================
    // latitude dimension
    if (nc_def_dim(ncId, "lat", mesh->getNumGrid(1, FULL), &latDimId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define dimension lat!");
    }
    if (nc_def_var(ncId, "lat", NC_DOUBLE, 1, &latDimId, &latVarId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define coordinate variable lat!");
    }
    nc_put_att(ncId, latVarId, "long_name", NC_CHAR, 8, "latitude");
    nc_put_att(ncId, latVarId, "units", NC_CHAR, 13, "degrees_north");
    // =========================================================================
    dimIds[0] = latDimId;
    dimIds[1] = lonDimId;
    for (int i = 0; i < tracerManager.getNumSpecies(); ++i) {
        const TracerSpeciesInfo &speciesInfo = tracerManager.getSpeciesInfo(i);
        if (nc_def_var(ncId, speciesInfo.getName().c_str(), NC_DOUBLE,
                       domain->getNumDim(), dimIds, &qVarIds[i])
            != NC_NOERR) {
            REPORT_ERROR("Failed to define variable " <<
                         speciesInfo.getName().c_str() << "!");
        }
        nc_put_att(ncId, qVarIds[i], "long_name", NC_CHAR,
                   speciesInfo.getBrief().length(), speciesInfo.getBrief().c_str());
        nc_put_att(ncId, qVarIds[i], "units", NC_CHAR,
                   speciesInfo.getUnits().length(), speciesInfo.getUnits().c_str());
    }
    if (nc_def_var(ncId, "volume", NC_DOUBLE, domain->getNumDim(), dimIds, &volVarId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define variable volume!");
    }
    nc_put_att(ncId, volVarId, "long_name", NC_CHAR, 16, "mesh cell volume");
    nc_enddef(ncId);
    // -------------------------------------------------------------------------
    // put variables
    // =========================================================================
    lon = mesh->getGridCoords(0, FULL)/RAD;
    lat = mesh->getGridCoords(1, FULL)/RAD;
    nc_put_var(ncId, lonVarId, lon.memptr());
    nc_put_var(ncId, latVarId, lat.memptr());
    double *x  = new double[mesh->getTotalNumGrid(CENTER)];
    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            x[i] = tracerMeshCells(oldTimeIdx, i).getSpeciesDensity(s);
        }
        nc_put_var(ncId, qVarIds[s], x);
    }
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
        x[i] = tracerMeshCells(oldTimeIdx, i).getVolume();
    }
    nc_put_var(ncId, volVarId, x);
    delete [] x;
    nc_close(ncId);
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
                               const geomtk::RLLVelocityField &velocity) {
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
    splitTracers(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("splitTracers uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    mergeTracers(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("mergeTracers uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    remapTracersToMesh(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("remapTracersToMesh uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");
}

// -----------------------------------------------------------------------------
// private member functions

void AdvectionManager::calcTotalMass(const TimeLevelIndex<2> &timeIdx) {
    // print total mass for each species
    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
        totalMass.getLevel(timeIdx)[s] = 0;
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            totalMass.getLevel(timeIdx)[s] +=
                tracerMeshCells(timeIdx, i).getSpeciesDensity(s)*
                tracerMeshCells(timeIdx, i).getVolume();
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
    numLongTracer = 0;
    list<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        Velocity v1(domain->getNumDim());
        Velocity v2(domain->getNumDim());
        Velocity v3(domain->getNumDim());
        Velocity v4(domain->getNumDim());
        Velocity v(domain->getNumDim());
        const VelocityField::FieldType &divergence = velocity.getDivergence();
        double div;
        vec rho(tracerManager.tracers.size());
        double k1_rho[tracerManager.getNumSpecies()];
        double k2_rho[tracerManager.getNumSpecies()];
        double k3_rho[tracerManager.getNumSpecies()];
        double k4_rho[tracerManager.getNumSpecies()];
        // ---------------------------------------------------------------------
        // update location and deformation matrix of tracer
        SpaceCoord &x0 = (*tracer)->getX(oldTimeIdx);
        SpaceCoord &x1 = (*tracer)->getX(newTimeIdx);
        MeshIndex &idx0 = (*tracer)->getMeshIndex(oldTimeIdx);
        MeshIndex &idx1 = (*tracer)->getMeshIndex(newTimeIdx);
        // TODO: Should we hide the following codes? Because they are
        //       related to sphere domain.
        if (idx0.isOnPole()) {
            idx0.setMoveOnPole(true);
            idx1.setMoveOnPole(true);
            x0.transformToPS(*domain);
        } else {
            idx0.setMoveOnPole(false);
            idx1.setMoveOnPole(false);
        }
        // =====================================================================
        // stage 1
        regrid->run(BILINEAR, oldTimeIdx, velocity, x0, v1, &idx0);
        regrid->run(BILINEAR, oldTimeIdx, divergence, x0, div, &idx0);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            k1_rho[s] = -(*tracer)->getSpeciesDensity(s)*div;
            rho[s] = (*tracer)->getSpeciesDensity(s)+dt05*k1_rho[s];
        }
        mesh->move(x0, dt05, v1, idx0, x1);
        idx1.locate(*mesh, x1);
        // =====================================================================
        // stage 2
        regrid->run(BILINEAR, halfTimeIdx, velocity, x1, v2, &idx1);
        regrid->run(BILINEAR, halfTimeIdx, divergence, x1, div, &idx1);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            k2_rho[s] = -rho[s]*div;
            rho[s] = (*tracer)->getSpeciesDensity(s)+dt05*k2_rho[s];
        }
        mesh->move(x0, dt05, v2, idx0, x1);
        idx1.locate(*mesh, x1);
        // =====================================================================
        // stage 3
        regrid->run(BILINEAR, halfTimeIdx, velocity, x1, v3, &idx1);
        regrid->run(BILINEAR, newTimeIdx, divergence, x1, div, &idx1);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            k3_rho[s] = -rho[s]*div;
            rho[s] = (*tracer)->getSpeciesDensity(s)+dt*k3_rho[s];
        }
        mesh->move(x0, dt, v3, idx0, x1);
        idx1.locate(*mesh, x1);
        // =====================================================================
        // stage 4
        regrid->run(BILINEAR, newTimeIdx, velocity, x1, v4, &idx1);
        regrid->run(BILINEAR, newTimeIdx, divergence, x1, div, &idx1);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            k4_rho[s] = -rho[s]*div;
            (*tracer)->getSpeciesDensity(s) += dt*
                (k1_rho[s]+2.0*k2_rho[s]+2.0*k3_rho[s]+k4_rho[s])/6.0;
        }
        v = (v1+v2*2.0+v3*2.0+v4)/6.0;
        mesh->move(x0, dt, v, idx0, x1);
        idx1.locate(*mesh, x1);
        x1.transformToCart(*domain);
        // ---------------------------------------------------------------------
        // update skeleton points of tracer
        TracerSkeleton &s = (*tracer)->getSkeleton();
        vector<SpaceCoord*> &x0s = s.getSpaceCoords(oldTimeIdx);
        vector<SpaceCoord*> &x1s = s.getSpaceCoords(newTimeIdx);
        vector<MeshIndex*> &idx0s = s.getMeshIdxs(oldTimeIdx);
        vector<MeshIndex*> &idx1s = s.getMeshIdxs(newTimeIdx);
        for (int i = 0; i < x0s.size(); ++i) {
            if (idx0s[i]->isOnPole()) {
                idx0s[i]->setMoveOnPole(true);
                idx1s[i]->setMoveOnPole(true);
                x0s[i]->transformToPS(*domain);
            } else {
                idx0s[i]->setMoveOnPole(false);
                idx1s[i]->setMoveOnPole(false);
            }
            // =================================================================
            // stage 1
            regrid->run(BILINEAR, oldTimeIdx, velocity, *x0s[i], v1, idx0s[i]);
            mesh->move(*x0s[i], dt05, v1, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(*mesh, *x1s[i]);
            // =================================================================
            // stage 2
            regrid->run(BILINEAR, halfTimeIdx, velocity, *x1s[i], v2, idx1s[i]);
            mesh->move(*x0s[i], dt05, v2, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(*mesh, *x1s[i]);
            // =================================================================
            // stage 3
            regrid->run(BILINEAR, halfTimeIdx, velocity, *x1s[i], v3, idx1s[i]);
            mesh->move(*x0s[i], dt, v3, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(*mesh, *x1s[i]);
            // =================================================================
            // stage 4
            regrid->run(BILINEAR, newTimeIdx, velocity, *x1s[i], v4, idx1s[i]);
            v = (v1+v2*2.0+v3*2.0+v4)/6.0;
            mesh->move(*x0s[i], dt, v, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(*mesh, *x1s[i]);
            x1s[i]->transformToCart(*domain);
        }
        // ---------------------------------------------------------------------
        (*tracer)->updateDeformMatrix(*domain, *mesh, newTimeIdx);
        if ((*tracer)->getBadType() == Tracer::EXTREME_FILAMENTATION) {
            if (numLongTracer == longTracers.size()) {
                longTracers.push_back(tracer);
            } else {
                longTracers[numLongTracer] = tracer;
            }
            numLongTracer++;
        }
    }
}

void AdvectionManager::embedTracersIntoMesh(const TimeLevelIndex<2> &timeIdx) {
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
        tracerMeshCells(timeIdx, i).resetContainedTracers();
    }
    list<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        int i = mesh->wrapIndex((*tracer)->getMeshIndex(timeIdx), CENTER);
        TracerMeshCell &cell = tracerMeshCells(timeIdx, i);
        cell.contain(*tracer);
    }
}

void AdvectionManager::connectTracerAndMesh(const TimeLevelIndex<2> &timeIdx,
                                            list<Tracer*>::iterator &tracer) {
    // call mlpack::range::RangeSearch to find out the neighbor cells of tracers
    // and set the data structures for both cells and tracers for remapping
    Searcher a(cellTree, NULL, cellCoords,
               (*tracer)->getX(timeIdx).getCartCoord(), true);
    double longAxisSize = (*tracer)->getShapeSize(timeIdx)[0];
    mlpack::math::Range r(0.0, longAxisSize);
    vector<vector<size_t> > neighbors;
    vector<vector<double> > distances;
    a.Search(r, neighbors, distances);
    BodyCoord y(domain->getNumDim());
    for (int i = 0; i < neighbors[0].size(); ++i) {
        int cellIdx = cellCoordsMap[neighbors[0][i]];
        TracerMeshCell *cell = &tracerMeshCells(timeIdx, cellIdx);
        // calculate the tracer shape function for the cell
        (*tracer)->getBodyCoord(*domain, timeIdx, cell->getCoord(), y);
        double f = (*tracer)->getShapeFunction(timeIdx, y);
        if (f > 0.0) {
            cell->connect(*tracer, f);
            (*tracer)->connect(cell, f);
        }
    }
    if ((*tracer)->getConnectedCells().size() == 0) {
        // tracer has not connect with any cells, so connect with its host cell
        double weight = ShapeFunction::getMaxValue();
        (*tracer)->getHostCell()->connect(*tracer, weight);
        (*tracer)->connect((*tracer)->getHostCell(), weight);
// The following codes are commented to ensure the split small tracers are
// merged by other tracers, and the total tracer number is not increased.
//        if ((*tracer)->getDetH(timeIdx) < (*tracer)->getHostCell()->getVolume()) {
//            // tracer is not resolved by the mesh
//            (*tracer)->setBadType(Tracer::NOT_RESOLVED);
//            if (numUnresolvedTracer == unresolvedTracers.size()) {
//                unresolvedTracers.push_back(tracer);
//            } else {
//                unresolvedTracers[numUnresolvedTracer] = tracer;
//            }
//            numUnresolvedTracer++;
//        } else {
//#ifndef NDEBUG
//            (*tracer)->dump(timeIdx, *domain);
//            CHECK_POINT;
//#endif
//        }
//    } else if ((*tracer)->getBadType() == Tracer::NOT_RESOLVED) {
//#ifndef NDEBUG
//        (*tracer)->dump(timeIdx, *domain);
//        CHECK_POINT;
//#endif
//        (*tracer)->setBadType(Tracer::GOOD_SHAPE);
    }
}

void AdvectionManager::connectTracersAndMesh(const TimeLevelIndex<2> &timeIdx) {
    numUnresolvedTracer = 0;
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
        tracerMeshCells(timeIdx, i).resetConnectedTracers();
    }
    list<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        (*tracer)->resetConnectedCells();
        connectTracerAndMesh(timeIdx, tracer);
    }
// Use the following code to show the connected tracers with the given cell.
//    TracerMeshCell &cell = tracerMeshCells(timeIdx, 19978);
//    std::ofstream file("connected_tracers.txt");
//    for (int i = 0; i < cell.getNumConnectedTracer(); ++i) {
//        cout << i+1 << ": " << cell.getConnectedTracers()[i]->getID() << endl;
//        cell.getConnectedTracers()[i]->outputShape(timeIdx, *domain, file, i);
//    }
//    file.close();
//    CHECK_POINT;
}

void AdvectionManager::splitTracers(const TimeLevelIndex<2> &timeIdx) {
    for (int t = 0; t < numLongTracer; ++t) {
        list<Tracer*>::iterator tracer = longTracers[t];
//#define CHECK_SPLIT_TRACER
#ifdef CHECK_SPLIT_TRACER
        std::ofstream file("split_parcel.txt");
        (*tracer)->dump(timeIdx, *domain, file, 0);
        int idx = 1;
#endif
        // ---------------------------------------------------------------------
        const SpaceCoord &c0 = (*tracer)->getX(timeIdx);
        Tracer *tracer1 = new Tracer(domain->getNumDim());
        Tracer *tracer2 = new Tracer(domain->getNumDim());
        const mat &H = (*tracer)->getH(timeIdx);
        const mat &invH = (*tracer)->getInvH(timeIdx);
        const mat &V = (*tracer)->getV();
        vec S;
        // ---------------------------------------------------------------------
        // add two new tracers along major axis of the needle tracer
        S = (*tracer)->getS(); S[0] *= 1-alpha; S /= sqrt(2.0);
        SpaceCoord x(domain->getNumDim()), x0(domain->getNumDim());
        BodyCoord y(domain->getNumDim()), y0(domain->getNumDim());
        y0() = invH*H*V.col(0);
        // first tracer
        SpaceCoord &c1 = tracer1->getX(timeIdx);
        y() = y0()*alpha;
        (*tracer)->getSpaceCoord(*domain, timeIdx, y, c1);
        c1.transformToCart(*domain);
        tracer1->getMeshIndex(timeIdx).locate(*mesh, tracer1->getX(timeIdx));
        tracerMeshCells(timeIdx, mesh->wrapIndex(tracer1->getMeshIndex(timeIdx),
                                                 CENTER)).contain(tracer1);
        tracer1->resetDeformMatrix(*domain, *mesh, timeIdx, c0, S);
        tracer1->setID(tracerManager.tracers.back()->getID()+1);
        tracerManager.tracers.push_back(tracer1);
        connectTracerAndMesh(timeIdx, --tracerManager.tracers.end());
        // merge the small tracer at the same step
        tracer1->setBadType(Tracer::NOT_RESOLVED);
        if (numUnresolvedTracer == unresolvedTracers.size()) {
            unresolvedTracers.push_back(--tracerManager.tracers.end());
        } else {
            unresolvedTracers[numUnresolvedTracer] = --tracerManager.tracers.end();
        }
        numUnresolvedTracer++;
        tracer1->fatherID = (*tracer)->getID(); // do not merged by the original tracer
        // second tracer
        SpaceCoord &c2 = tracer2->getX(timeIdx);
        y() *= -1;
        (*tracer)->getSpaceCoord(*domain, timeIdx, y, c2);
        c2.transformToCart(*domain);
        tracer2->getMeshIndex(timeIdx).locate(*mesh, tracer2->getX(timeIdx));
        tracerMeshCells(timeIdx, mesh->wrapIndex(tracer2->getMeshIndex(timeIdx),
                                                 CENTER)).contain(tracer2);
        tracer2->resetDeformMatrix(*domain, *mesh, timeIdx, c0, S);
        tracer2->setID(tracerManager.tracers.back()->getID()+1);
        tracerManager.tracers.push_back(tracer2);
        connectTracerAndMesh(timeIdx, --tracerManager.tracers.end());
        // merge the small tracer at the same step
        tracer2->setBadType(Tracer::NOT_RESOLVED);
        if (numUnresolvedTracer == unresolvedTracers.size()) {
            unresolvedTracers.push_back(--tracerManager.tracers.end());
        } else {
            unresolvedTracers[numUnresolvedTracer] = --tracerManager.tracers.end();
        }
        numUnresolvedTracer++;
        tracer2->fatherID = (*tracer)->getID();
        // ---------------------------------------------------------------------
        // share part of tracer density and mass to the new tracers
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            tracer1->addSpecies();
            tracer1->getSpeciesDensity(s) = (*tracer)->getSpeciesDensity(s);
            tracer1->calcSpeciesMass(timeIdx, s);
#ifndef NDEBUG
            assert(fabs(tracer1->getSpeciesMass(s)-tracer1->getSpeciesDensity(s)*tracer1->getDetH(timeIdx)) < 1.0e-12);
#endif
            tracer2->addSpecies();
            tracer2->getSpeciesDensity(s) = (*tracer)->getSpeciesDensity(s);
            tracer2->calcSpeciesMass(timeIdx, s);
#ifndef NDEBUG
            assert(fabs(tracer2->getSpeciesMass(s)-tracer2->getSpeciesDensity(s)*tracer2->getDetH(timeIdx)) < 1.0e-12);
#endif
        }
        // ---------------------------------------------------------------------
        // shrink the old parcel
        S = (*tracer)->getS(); S[0] *= alpha;
        (*tracer)->resetDeformMatrix(*domain, *mesh, timeIdx, c1, S);
        for (int i = 0; i < (*tracer)->getNumConnectedCell(); ++i) {
            TracerMeshCell *cell = (*tracer)->getConnectedCells()[i];
            cell->disconnect(*tracer);
        }
        (*tracer)->resetConnectedCells();
        connectTracerAndMesh(timeIdx, tracer);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
#ifndef NDEBUG
            double m = (*tracer)->getSpeciesMass(s);
#endif
            (*tracer)->calcSpeciesMass(timeIdx, s);
#ifndef NDEBUG
            double m0 = (*tracer)->getSpeciesMass(s);
            double m1 = tracer1->getSpeciesMass(s);
            double m2 = tracer2->getSpeciesMass(s);
            if (m != 0) {
                assert(fabs(m-m0-m1-m2)/m < 1.0e-12);
            }
#endif
        }
        (*tracer)->setBadType(Tracer::GOOD_SHAPE);
        // ---------------------------------------------------------------------
#ifndef NDEBUG
        REPORT_NOTICE("Long tracer " << (*tracer)->getID() << " is split into "
                      << tracer1->getID() << " and " << tracer2->getID() << ".");
#endif
#ifdef CHECK_SPLIT_TRACER
        (*tracer)->dump(timeIdx, *domain, file, idx++);
        tracer1->dump(timeIdx, *domain, file, idx++);
        tracer2->dump(timeIdx, *domain, file, idx++);
        file.close();
#endif
    }
}

void AdvectionManager::mergeTracers(const TimeLevelIndex<2> &timeIdx) {
    for (int t = 0; t < numUnresolvedTracer; ++t) {
        list<Tracer*>::iterator tracer = unresolvedTracers[t];
        TracerMeshCell *cell = (*tracer)->getHostCell();
        vector<Tracer*> &tracers = cell->getConnectedTracers();
        // ---------------------------------------------------------------------
        // merge unresolved tracer
//#define CHECK_MERGE_TRACER
#ifdef CHECK_MERGE_TRACER
        std::ofstream file("merge_parcel.txt");
        (*tracer)->dump(timeIdx, *domain, file, 0);
        int idx = 1;
        arma::uvec tags(cell->getNumConnectedTracer(), arma::fill::zeros);
#endif
        const mat &H = (*tracer)->getH(timeIdx);
        const mat &invH = (*tracer)->getInvH(timeIdx);
        const mat &V = (*tracer)->getV();
        BodyCoord y(domain->getNumDim()), y0(domain->getNumDim());
        y0() = invH*H*V.col(0);
        vec weights(cell->getNumConnectedTracer(), arma::fill::zeros);
#ifndef NDEBUG
        vec totalMass1(tracerManager.getNumSpecies(), arma::fill::zeros);
        vec totalMass2(tracerManager.getNumSpecies(), arma::fill::zeros);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            totalMass1[s] += (*tracer)->getSpeciesMass(s);
        }
#endif
        vec d1(cell->getNumConnectedTracer(), arma::fill::zeros);
        vec d2(cell->getNumConnectedTracer(), arma::fill::zeros);
        for (int i = 0; i < cell->getNumConnectedTracer(); ++i) {
            if (tracers[i]->getID() != (*tracer)->getID() &&
                tracers[i]->getID() != (*tracer)->fatherID &&
                tracers[i]->getBadType() == Tracer::GOOD_SHAPE) {
                (*tracer)->getBodyCoord(*domain, timeIdx, tracers[i]->getX(timeIdx), y);
                double cosTheta = norm_dot(y(), y0());
                double sinTheta = sqrt(1-cosTheta*cosTheta);
                double n = norm(y(), 2);
                d1[i] = n*cosTheta;
                d2[i] = n*sinTheta;
            }
        }
        double scale = 1;
        while (true) {
            weights.zeros();
            for (int i = 0; i < cell->getNumConnectedTracer(); ++i) {
                if (tracers[i]->getID() != (*tracer)->getID() &&
                    tracers[i]->getID() != (*tracer)->fatherID &&
                    tracers[i]->getBadType() == Tracer::GOOD_SHAPE) {
                    weights[i] = exp(-scale*(beta1*d1[i]*d1[i]+beta2*d2[i]*d2[i]));
                }
            }
            double sumWeights = sum(weights);
            if (sumWeights < 1.0e-15) {
                scale *= 0.99;
            } else {
                weights /= sumWeights;
                break;
            }
        }
#ifdef CHECK_MERGE_TRACER
        for (int i = 0; i < cell->getNumConnectedTracer(); ++i) {
            if (tracers[i]->getID() != (*tracer)->getID() &&
                tracers[i]->getID() != (*tracer)->fatherID &&
                tracers[i]->getBadType() == Tracer::GOOD_SHAPE) {
                tracers[i]->dump(timeIdx, *domain, file, idx++);
                tags[i] = 1;
            }
        }
        file << "weights = (/";
        arma::uvec tmp = find(tags == 1, 1, "last");
        for (int i = 0; i < cell->getNumConnectedTracer(); ++i) {
            if (tags[i] == 1) {
                if (i == tmp[0]) {
                    file << weights[i] << "/)" << endl;
                } else {
                    file << weights[i] << ",";
                }
            }
        }
        file.close();
#endif
        double volume1 = (*tracer)->getDetH(timeIdx);
        for (int i = 0; i < cell->getNumConnectedTracer(); ++i) {
            if (tracers[i]->getID() != (*tracer)->getID() &&
                tracers[i]->getID() != (*tracer)->fatherID &&
                tracers[i]->getBadType() == Tracer::GOOD_SHAPE) {
                double volume2 = tracers[i]->getDetH(timeIdx);
                double volume = volume1*weights[i]+volume2;
                const mat &U = tracers[i]->getU();
                const mat &V = tracers[i]->getV();
                vec &S = tracers[i]->getS();
                S *= pow(volume/volume2, 1.0/domain->getNumDim());
                tracers[i]->getH(timeIdx) = U*diagmat(S)*V.t();
                tracers[i]->getDetH(timeIdx) = volume;
                tracers[i]->getInvH(timeIdx) = inv(tracers[i]->getH(timeIdx));
                tracers[i]->resetSkeleton(*domain, *mesh, timeIdx);
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    double m1 = (*tracer)->getSpeciesMass(s);
                    double &m2 = tracers[i]->getSpeciesMass(s);
#ifndef NDEBUG
                    totalMass1[s] += m2;
                    double rho1 = (*tracer)->getSpeciesDensity(s);
                    double rho2 = tracers[i]->getSpeciesDensity(s);
#endif
                    m2 = m1*weights[i]+m2;
                    tracers[i]->calcSpeciesDensity(timeIdx, s);
#ifndef NDEBUG
                    totalMass2[s] += m2;
                    double rho3 = tracers[i]->getSpeciesDensity(s);
                    assert(((rho1 <= rho3 || (rho1-rho3)/rho1 <= 1.0e-12) &&
                            (rho3 <= rho2 || (rho3-rho2)/rho2 <= 1.0e-12)) ||
                           ((rho2 <= rho3 || (rho2-rho3)/rho2 <= 1.0e-12) &&
                            (rho3 <= rho1 || (rho3-rho1)/rho1 <= 1.0e-12)));
#endif
                }
            }
        }
#ifndef NDEBUG
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            assert(fabs(totalMass1[s]-totalMass2[s])/totalMass1[s] < 1.0e-12);
        }
#endif
        // ---------------------------------------------------------------------
        // remove unsolved tracer
        for (int i = 0; i < (*tracer)->getNumConnectedCell(); ++i) {
            TracerMeshCell *cell = (*tracer)->getConnectedCells()[i];
            cell->disconnect(*tracer);
        }
        (*tracer)->getHostCell()->discontain(*tracer);
        tracerManager.tracers.erase(tracer);
    }
}
    
void AdvectionManager::handleVoidCells(const TimeLevelIndex<2> &timeIdx) {
    // NOTE: The occurance of void cells should be rare.
    for (int c = 0; c < numVoidCell; ++c) {
        TracerMeshCell *cell = voidCells[c];
        Searcher a(cellTree, NULL, cellCoords,
                   cell->getCoord().getCartCoord(), true);
        double searchRadius = 1*RAD*domain->getRadius();
        while (true) {
            mlpack::math::Range r(0.0, searchRadius);
            vector<vector<size_t> > neighbors;
            vector<vector<double> > distances;
            a.Search(r, neighbors, distances);
            if (neighbors[0].size() != 0) {
                vector<TracerMeshCell*> ngbCells;
                for (int i = 0; i < neighbors[0].size(); ++i) {
                    int cellIdx = cellCoordsMap[neighbors[0][i]];
                    TracerMeshCell *ngbCell = &tracerMeshCells(timeIdx, cellIdx);
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
                    double d = domain->calcDistance(cell->getCoord(),
                                                    ngbCells[i]->getCoord());
                    weights[i] = 1/d;
                }
                weights /= sum(weights);
                for (int i = 0; i < ngbCells.size(); ++i) {
                    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                        cell->getSpeciesDensity(s) += ngbCells[i]->getSpeciesDensity(s)*weights[i];
                    }
                }
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    cell->getSpeciesMass(s) = cell->getSpeciesDensity(s)*cell->getVolume();
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
    list<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        (*tracer)->resetSpecies();
        const vector<TracerMeshCell*> &cells = (*tracer)->getConnectedCells();
        double totalWeight = 0;
#if defined REMAP_DENSITY
        for (int i = 0; i < (*tracer)->getNumConnectedCell(); ++i) {
            double weight = cells[i]->getRemapWeight(*tracer);
            totalWeight += weight;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                (*tracer)->getSpeciesDensity(s) += cells[i]->getSpeciesDensity(s)*weight;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            (*tracer)->getSpeciesDensity(s) /= totalWeight;
            (*tracer)->calcSpeciesMass(timeIdx, s);
        }
#elif defined REMAP_MASS
        for (int i = 0; i < (*tracer)->getNumConnectedCell(); ++i) {
            double weight = cells[i]->getRemapWeight(*tracer);
            totalWeight += weight;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                (*tracer)->getSpeciesMass(s) += cells[i]->getSpeciesMass(s)*weight;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            (*tracer)->getSpeciesMass(s) /= totalWeight;
            (*tracer)->getSpeciesDensity(s) = (*tracer)->getSpeciesMass(s)/(*tracer)->getDetH(timeIdx);
        }
#endif
    }
}

void AdvectionManager::remapTracersToMesh(const TimeLevelIndex<2> &timeIdx) {
    numVoidCell = 0;
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
        TracerMeshCell &cell = tracerMeshCells(timeIdx, i);
        cell.resetSpecies();
        if (cell.getNumConnectedTracer() == 0) {
            if (numVoidCell == voidCells.size()) {
                voidCells.push_back(&cell);
            } else {
                voidCells[numVoidCell] = &cell;
            }
            numVoidCell++;
            continue;
        };
        vector<Tracer*> &tracers = cell.getConnectedTracers();
        double totalWeight = 0;
#if defined REMAP_DENSITY
        for (int j = 0; j < cell.getNumConnectedTracer(); ++j) {
            double weight = cell.getRemapWeight(tracers[j]);
            totalWeight += weight;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                cell.getSpeciesDensity(s) += tracers[j]->getSpeciesDensity(s)*weight;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            cell.getSpeciesDensity(s) /= totalWeight;
            cell.getSpeciesMass(s) = cell.getSpeciesDensity(s)*cell.getVolume();
        }
#elif defined REMAP_MASS
        for (int j = 0; j < cell.getNumConnectedTracer(); ++j) {
            double weight = cell.getRemapWeight(tracers[j]);
            totalWeight += weight;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                cell.getSpeciesMass(s) += tracers[j]->getSpeciesMass(s)*weight;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            cell.getSpeciesMass(s) /= totalWeight;
            cell.getSpeciesDensity(s) = cell.getSpeciesMass(s)/cell.getVolume();
        }
#endif
    }
    handleVoidCells(timeIdx);
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
#ifndef NDEBUG
        double biasPercent = (totalMass.getLevel(timeIdx)[s]-expectedTotalMass[s])/expectedTotalMass[s];
        REPORT_NOTICE("Mass conservation bias percentage is " << std::fixed << setw(10) << setprecision(4) << biasPercent*100 << "%.");
#endif
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            tracerMeshCells(timeIdx, i).getSpeciesDensity(s) *= fixer;
        }
        totalMass.getLevel(timeIdx)[s] = expectedTotalMass[s];
    }
}

}
