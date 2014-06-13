#include "AdvectionManager.h"
#include "ShapeFunction.h"
#include "TracerSkeleton.h"

namespace lasm {

double transitionFunction(double a0, double a1, double b0, double b1, double x) {
    if (x <= b0) {
        return a0;
    } else if (x > b0 && x <= b1) {
        double t = (x-b0)/(b1-b0);
        return (a1-a0)*(4-3*t)*pow(t, 3)+a0;
    } else {
        return a1;
    }
}

AdvectionManager::AdvectionManager() {
    domain = NULL;
    mesh = NULL;
    timeManager = NULL;
    outputFileFormat = NULL;
    regrid = NULL;
    cellTree = NULL;
    numFilamentTracer = 0;
    numVoidCell = 0;
    filamentLimit = 2;
    radialMixing = 1;
    lateralMixing = 10;
    isMassFixed = true;
    REPORT_ONLINE;
}

AdvectionManager::~AdvectionManager() {
    if (outputFileFormat != NULL) {
        delete outputFileFormat;
    }
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
    TimeLevelIndex<2> initTimeIdx;
#ifndef NDEBUG
    assert(initTimeIdx.get() == 0);
#endif
    if (configManager.hasKey("lasm", "filament_limit")) {
        configManager.getValue("lasm", "filament_limit", filamentLimit);
    }
    if (configManager.hasKey("lasm", "radial_mixing")) {
        configManager.getValue("lasm", "radial_mixing", radialMixing);
    }
    if (configManager.hasKey("lasm", "lateral_mixing")) {
        configManager.getValue("lasm", "lateral_mixing", lateralMixing);
    }
    if (configManager.hasKey("lasm", "is_mass_fixed")) {
        configManager.getValue("lasm", "is_mass_fixed", isMassFixed);
    }
    string outputFilePrefix;
    if (configManager.hasKey("lasm", "output_prefix")) {
        configManager.getValue("lasm", "output_prefix", outputFilePrefix);
    } else {
        outputFilePrefix = "lasm";
    }
    outputFileFormat = new StampString(outputFilePrefix+".", ".nc");
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

void AdvectionManager::output(const TimeLevelIndex<2> &oldTimeIdx) {
    string fileName = outputFileFormat->run("%5.5d", timeManager->getNumStep());
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
    if (nc_put_att(ncId, NC_GLOBAL, "filament_limit", NC_DOUBLE, 1,
                   &filamentLimit) != NC_NOERR){
        REPORT_ERROR("Failed to put global dimension \"filament_limit\"!");
    }
    if (nc_put_att(ncId, NC_GLOBAL, "radial_mixing", NC_DOUBLE, 1,
                   &radialMixing) != NC_NOERR){
        REPORT_ERROR("Failed to put global dimension \"radial_mixing\"!");
    }
    if (nc_put_att(ncId, NC_GLOBAL, "lateral_mixing", NC_DOUBLE, 1,
                   &lateralMixing) != NC_NOERR){
        REPORT_ERROR("Failed to put global dimension \"lateral_mixing\"!");
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
    remapTracersToMesh(newTimeIdx, &velocity);
    time2 = clock();
    REPORT_NOTICE("remapTracersToMesh uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    checkTracerShapes(newTimeIdx, velocity);
    time2 = clock();
    REPORT_NOTICE("checkTracerShapes uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    handleFilamentTracers(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("handleFilamentTracers uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");
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
    const VelocityField::FieldType &divergence = velocity.getDivergence();
    list<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
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
        double w = (*tracer)->getShapeFunction(timeIdx, y);
        if (w > 0.0) {
            double d = domain->calcDistance((*tracer)->getX(timeIdx),
                                            cell->getCoord());
            cell->connect(*tracer, w, d);
            (*tracer)->connect(cell);
        }
    }
    if ((*tracer)->getConnectedCells().size() == 0) {
        // tracer has not connect with any cells, so connect with its host cell
        double w = ShapeFunction::getMaxValue();
        double d = domain->calcDistance((*tracer)->getX(timeIdx),
                                        (*tracer)->getHostCell()->getCoord());
        (*tracer)->getHostCell()->connect(*tracer, w, d);
        (*tracer)->connect((*tracer)->getHostCell());
    }
}

void AdvectionManager::connectTracersAndMesh(const TimeLevelIndex<2> &timeIdx) {
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
        tracerMeshCells(timeIdx, i).resetConnectedTracers();
    }
    list<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        (*tracer)->resetConnectedCells();
        connectTracerAndMesh(timeIdx, tracer);
    }
}

void AdvectionManager::checkTracerShapes(const TimeLevelIndex<2> &timeIdx,
                                         const VelocityField &velocity) {
    const VelocityField::FieldType &vorticity = velocity.getVorticity()[0];
    numFilamentTracer = 0;
    list<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        // check the degree of filamentation
        const vec &S = (*tracer)->getS();
        double vor = fabs(vorticity(timeIdx, (*tracer)->getHostCell()->getID()));
        double actualLimit = transitionFunction(filamentLimit, 5, 1.0e-9, 1.0e-5, vor);
        if (S[0]/S[1] > actualLimit) {
            recordTracer(Tracer::FILAMENT, tracer);
        }
    }
}

void AdvectionManager::handleFilamentTracers(const TimeLevelIndex<2> &timeIdx) {
    const double alpha = 0.01;
    for (int t = 0; t < numFilamentTracer; ++t) {
        list<Tracer*>::iterator tracer = filamentTracers[t];
        // get surrounding tracers
        TracerMeshCell *cell = (*tracer)->getHostCell();
        const vector<Tracer*> &surroundTracers = cell->getConnectedTracers();
#ifndef NDEBUG
#ifdef CHECK_FILAMENT_TRACER
        std::ofstream file("filament_parcel.txt");
        (*tracer)->dump(timeIdx, *domain, file, 0);
        int idx = 1;
#endif
        assert(surroundTracers.size() > 1);
        vec totalMass1(tracerManager.getNumSpecies(), arma::fill::zeros);
        vec totalMass2(tracerManager.getNumSpecies(), arma::fill::zeros);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            totalMass1[s] += (*tracer)->getSpeciesMass(s);
        }
#endif
        // calcuate mixing weights
        mat &H = (*tracer)->getH(timeIdx);
        mat &invH = (*tracer)->getInvH(timeIdx);
        const mat &U = (*tracer)->getU();
        const mat &V = (*tracer)->getV();
        BodyCoord y(domain->getNumDim()), y0(domain->getNumDim());
        y0() = invH*H*V.col(0);
        vec weights(cell->getNumConnectedTracer(), arma::fill::zeros);
        for (int i = 0; i < cell->getNumConnectedTracer(); ++i) {
            Tracer *tracer1 = surroundTracers[i];
            if (tracer1 == *tracer) continue;
            (*tracer)->getBodyCoord(*domain, timeIdx, tracer1->getX(timeIdx), y);
            double cosTheta = norm_dot(y(), y0());
            double sinTheta = sqrt(1-cosTheta*cosTheta);
            double n = norm(y(), 2);
            double d1 = n*cosTheta;
            double d2 = n*sinTheta;
            weights[i] = exp(-(radialMixing*d1*d1+lateralMixing*d2*d2));
            if (weights[i] < 1.0e-10) {
                weights[i] = 0;
            }
#ifndef NDEBUG
            if (weights[i] != 0) {
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    totalMass1[s] += tracer1->getSpeciesMass(s);
                }
            }
            assert(weights[i] <= 1);
#ifdef CHECK_FILAMENT_TRACER
            tracer1->dump(timeIdx, *domain, file, idx++);
#endif
#endif

        }
        double sumWeights = sum(weights);
        if (sumWeights < 1.0e-14) {
            continue;
        }
        weights /= sumWeights;
#if (defined NDEBUG && defined CHECK_FILAMENT_TRACER)
        file << "weights = (/";
        for (int i = 0; i < cell->getNumConnectedTracer(); ++i) {
            if (i == cell->getNumConnectedTracer()-1) {
                file << weights[i] << "/)" << endl;
            } else {
                file << weights[i] << ",";
            }
        }
        file.close();
#endif
        // distribute mass to surrounding tracers
        double oldVolume = (*tracer)->getDetH(timeIdx);
        double newVolume = oldVolume*(1-alpha);
        for (int i = 0; i < cell->getNumConnectedTracer(); ++i) {
            Tracer *tracer1 = surroundTracers[i];
            if (weights[i] == 0) continue;
            double volume1 = tracer1->getDetH(timeIdx);
            double volume2 = alpha*oldVolume*weights[i]+volume1;
            assert(volume2 > volume1);
            mat &H1 = tracer1->getH(timeIdx);
            const mat &U1 = tracer1->getU();
            const mat &V1 = tracer1->getV();
            vec &S1 = tracer1->getS();
            S1 *= pow(volume2/volume1, 1.0/domain->getNumDim());
            H1 = U1*diagmat(S1)*V1.t();
            tracer1->getDetH(timeIdx) = volume2;
            tracer1->getInvH(timeIdx) = inv(H1);
            tracer1->resetSkeleton(*domain, *mesh, timeIdx);
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                double m1 = (*tracer)->getSpeciesMass(s)*alpha;
                double &m2 = tracer1->getSpeciesMass(s);
                m2 += m1*weights[i];
#ifndef NDEBUG
                totalMass2[s] += m2;
                double rho1 = (*tracer)->getSpeciesDensity(s);
                double rho2 = tracer1->getSpeciesDensity(s);
#endif
                tracer1->calcSpeciesDensity(timeIdx, s);
#ifndef NDEBUG
                double rho3 = tracer1->getSpeciesDensity(s);
                assert(((rho1 <= rho3 || (rho1-rho3)/rho1 <= 1.0e-12) &&
                        (rho3 <= rho2 || (rho3-rho2)/rho2 <= 1.0e-12)) ||
                       ((rho2 <= rho3 || (rho2-rho3)/rho2 <= 1.0e-12) &&
                        (rho3 <= rho1 || (rho3-rho1)/rho1 <= 1.0e-12)));
#endif
            }
        }
        // change tracer shape and mass
        vec &S = (*tracer)->getS();
        S[0] = newVolume/S[1];
        H = U*diagmat(S)*V.t();
        (*tracer)->getDetH(timeIdx) = newVolume;
        invH = inv(H);
        (*tracer)->resetSkeleton(*domain, *mesh, timeIdx);
        (*tracer)->setType(Tracer::GOOD_SHAPE);
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            (*tracer)->calcSpeciesMass(timeIdx, s);
#ifndef NDEBUG
            totalMass2[s] += (*tracer)->getSpeciesMass(s);
#endif
        }
#ifndef NDEBUG
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            assert(fabs(totalMass1[s]-totalMass2[s])/totalMass1[s] < 1.0e-12);
        }
#endif
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

void AdvectionManager::remapTracersToMesh(const TimeLevelIndex<2> &timeIdx,
                                          const VelocityField *velocity) {
    const VelocityField::FieldType *vorticity = NULL;
    if (velocity != NULL) {
        vorticity = &velocity->getVorticity()[0];
    }
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
        const vector<Tracer*> &tracers = cell.getConnectedTracers();
        bool useIDW = false;

        // TEST
        if (vorticity != NULL) {
            int idx[3];
            mesh->unwrapIndex(cell.getID(), idx, CENTER);
            vec::fixed<9> vors; vors.zeros();
            if (idx[1] != 0) {
                vors[0] = fabs((*vorticity)(timeIdx, idx[0]-1, idx[1]-1));
                vors[1] = fabs((*vorticity)(timeIdx, idx[0],   idx[1]-1));
                vors[2] = fabs((*vorticity)(timeIdx, idx[0]+1, idx[1]-1));
            }
            vors[3] = fabs((*vorticity)(timeIdx, idx[0]-1, idx[1]));
            vors[4] = fabs((*vorticity)(timeIdx, idx[0],   idx[1]));
            vors[5] = fabs((*vorticity)(timeIdx, idx[0]+1, idx[1]));
            if (idx[1] != mesh->getNumGrid(1, FULL)-1) {
                vors[6] = fabs((*vorticity)(timeIdx, idx[0]-1, idx[1]+1));
                vors[7] = fabs((*vorticity)(timeIdx, idx[0],   idx[1]+1));
                vors[8] = fabs((*vorticity)(timeIdx, idx[0]+1, idx[1]+1));
            }
            if (any(vors > 1.0e-5) && fabs((*vorticity)(timeIdx, cell.getID())) < 1.0e-5) {
                useIDW = true;
            }
        }
        
        if (vorticity != NULL && fabs((*vorticity)(timeIdx, cell.getID())) > 1.0e-5) {
            // turn on inverse distance weight remapping mode
            useIDW = true;
        }
        double W = 0;
#if defined REMAP_DENSITY
        for (int j = 0; j < cell.getNumConnectedTracer(); ++j) {
            double w;
            if (useIDW) {
                w = 1/cell.getRemapDistance(tracers[j]);
            } else {
                w = cell.getRemapWeight(tracers[j]);
            }
            W += w;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                cell.getSpeciesDensity(s) += tracers[j]->getSpeciesDensity(s)*w;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            cell.getSpeciesDensity(s) /= W;
            cell.getSpeciesMass(s) = cell.getSpeciesDensity(s)*cell.getVolume();
        }
#elif defined REMAP_MASS
        for (int j = 0; j < cell.getNumConnectedTracer(); ++j) {
            double w = cell.getRemapWeight(tracers[j]);
            W += w;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                cell.getSpeciesMass(s) += tracers[j]->getSpeciesMass(s)*w;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            cell.getSpeciesMass(s) /= W;
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

void AdvectionManager::recordTracer(Tracer::TracerType type,
                                    list<Tracer*>::iterator &tracer) {
    (*tracer)->setType(type);
    switch (type) {
        case Tracer::FILAMENT:
            if (numFilamentTracer == filamentTracers.size()) {
                filamentTracers.push_back(tracer);
            } else {
                filamentTracers[numFilamentTracer] = tracer;
            }
            numFilamentTracer++;
            break;
        default:
            REPORT_ERROR("Unexpected branch!");
            break;
    }
}

}
