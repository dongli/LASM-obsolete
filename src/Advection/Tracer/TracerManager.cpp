#include "TracerManager.h"
#include "TracerSkeleton.h"

namespace lasm
{

TracerManager::TracerManager() {
    scale0 = 1.5;
    REPORT_ONLINE;
}

TracerManager::~TracerManager() {
    for (int t = 0; t < tracers.size(); ++t) {
        delete tracers[t];
    }
    REPORT_OFFLINE;
}

#if defined USE_CARTESIAN_DOMAIN
void TracerManager::init(const Domain &domain, const Mesh &mesh,
                         const ConfigManager &configManager) {
    this->domain = &domain;
    this->mesh = &mesh;
    if (configManager.hasKey("lasm", "scale0")) {
        configManager.getValue("lasm", "scale0", scale0);
    }
    int numParcel = 0;
    int numParcelX, numParcelY;
    configManager.getValue("lasm", "num_parcel_x", numParcelX);
    configManager.getValue("lasm", "num_parcel_y", numParcelY);
    numParcel = numParcelX*numParcelY;
    tracers.resize(numParcel);
    int id = 0;
    for (int t = 0; t < tracers.size(); ++t) {
        tracers[t] = new Tracer(domain.numDim());
        tracers[t]->setID(id++);
    }
    TimeLevelIndex<2> timeIdx;
    for (int t = 0; t < tracers.size(); ++t) {
        Tracer *tracer = tracers[t];
        // set coordinate
        SpaceCoord &x0 = tracer->x(timeIdx);
        vec h(domain.numDim());
        double dx = domain.axisSpan(0)/numParcelX;
        double dy = domain.axisSpan(1)/numParcelY;
        int l = 0;
        for (int j = 0; j < numParcelY; ++j) {
            for (int i = 0; i < numParcelX; ++i) {
                if (l == tracer->getID()) {
                    double y = domain.axisStart(1)+dy*0.5+dy*j;
                    double x = domain.axisStart(0)+dx*0.5+dx*i;
                    x0.setCoord(x, y);
                    l = -1;
                    break;
                }
                l++;
            }
            if (l == -1) break;
        }
        h(0) = dx;
        h(1) = dy;
        h *= scale0;
        MeshIndex &idx0 = tracer->meshIndex(timeIdx);
        idx0.locate(mesh, x0);
        tracer->skeleton().init(domain, mesh, h.max());
        tracer->HtimeIdx).eye();
        tracer->updateDeformMatrix(domain, mesh, timeIdx);
    }
}
#elif defined USE_SPHERE_DOMAIN
void TracerManager::init(const Domain &domain, const Mesh &mesh,
                         const geomtk::ConfigManager &configManager) {
    this->domain = &domain;
    this->mesh = &mesh;
    if (configManager.hasKey("lasm", "scale0")) {
        configManager.getValue("lasm", "scale0", scale0);
    }
    TimeLevelIndex<2> timeIdx;
    int numParcel = 0;
    int numParcelX, numParcelY;
    configManager.getValue("lasm", "num_parcel_x", numParcelX);
    configManager.getValue("lasm", "num_parcel_y", numParcelY);
//#define USE_FULL_LAT_LON
#ifdef USE_FULL_LAT_LON
    numParcel = numParcelX*numParcelY;
#else
    // Note: Use reduced lat-lon mesh as the initial distribution of tracers.
    if (numParcelX%2 != 0) {
        REPORT_ERROR("Tracer number (now is " << numParcelX <<
                     ") along longitude axis should be even!");
    }
    double dlat = domain.axisSpan(1)/numParcelY;
    double shiftLat = 45.0*RAD;
    double minNumTracerX = 4;
    double cosLat0, cosLat1 = cos(shiftLat);
    for (int j = 0; j < numParcelY; ++j) {
        double lat = domain.axisStart(1)+dlat*0.5+dlat*j;
        int numParcelX_;
        if (fabs(lat) < shiftLat) {
            numParcelX_ = numParcelX;
        } else {
            if (j == 0) cosLat0 = fabs(cos(lat));
            double ratio = (fabs(cos(lat))-cosLat0)/(cosLat1-cosLat0);
            numParcelX_ = ratio*numParcelX+(1-ratio)*minNumTracerX;
            if (numParcelX_%2 != 0) numParcelX_++;
        }
        numParcel += numParcelX_;
    }
#endif
    tracers.resize(numParcel);
    int id = 0;
    for (int t = 0; t < tracers.size(); ++t) {
        tracers[t] = new Tracer(domain.numDim());
        tracers[t]->ID() = id++;
    }
    for (int t = 0; t < tracers.size(); ++t) {
        Tracer *tracer = tracers[t];
        // set coordinate
        SpaceCoord &x0 = tracer->x(timeIdx);
        vec h(domain.numDim());
#ifdef USE_FULL_LAT_LON
        double dlon = domain.axisSpan(0)/numParcelX;
        double dlat = domain.axisSpan(1)/numParcelY;
        int l = 0;
        for (int j = 0; j < numParcelY; ++j) {
            for (int i = 0; i < numParcelX; ++i) {
                if (l == tracer->ID()) {
                    double lat = domain.axisStart(1)+dlat*0.5+dlat*j;
                    double lon = domain.axisStart(0)+dlon*0.5+dlon*i;
                    x0.setCoord(lon, lat);
                    l = -1;
                    break;
                }
                l++;
            }
            if (l == -1) break;
        }
#else
        // Note: Use reduced lat-lon mesh as the initial distribution of tracers.
        double dlat = domain.axisSpan(1)/numParcelY;
        double shiftLat = 45.0*RAD;
        double minNumTracerX = 4;
        double cosLat0, cosLat1 = cos(shiftLat);
        double dlon;
        int l = 0;
        for (int j = 0; j < numParcelY; ++j) {
            double lat = domain.axisStart(1)+dlat*0.5+dlat*j;
            int numParcelX_;
            if (fabs(lat) < shiftLat) {
                numParcelX_ = numParcelX;
            } else {
                if (j == 0) cosLat0 = fabs(cos(lat));
                double ratio = (fabs(cos(lat))-cosLat0)/(cosLat1-cosLat0);
                numParcelX_ = ratio*numParcelX+(1-ratio)*minNumTracerX;
                if (numParcelX_%2 != 0) numParcelX_++;
            }
            for (int i = 0; i < numParcelX_; ++i) {
                if (l == tracer->ID()) {
                    dlon = domain.axisSpan(0)/numParcelX_;
                    double lon = domain.axisStart(0)+dlon*0.5+dlon*i;
                    x0.setCoord(lon, lat);
                    l = -1;
                    break;
                }
                l++;
            }
            if (l == -1) break;
        }
#endif
        h(0) = dlon*domain.radius()*x0.cosLat();
        h(1) = dlat*domain.radius();
        // When h is small, there may be small spots on the grid densities.
        // The selection of h can be tricky. When h *= 1.5, the errors of solid
        // rotation test is oscillatory with time.
        h *= scale0;
        // set tracer skeleton
        tracer->skeleton().init(domain, mesh, h.max());
        // set deformation matrix
        tracer->H(timeIdx).eye();
        tracer->updateDeformMatrix(domain, mesh, timeIdx);
    }
    // Do the rest initialization jobs.
    for (int t = 0; t < tracers.size(); ++t) {
        Tracer *tracer = tracers[t];
        tracer->x(timeIdx).transformToCart(domain);
        tracer->meshIndex(timeIdx).locate(mesh, tracer->x(timeIdx));
        // TODO: This may be unnecessary.
        // when tracer is on Pole, transform its coordinate to PS for later use
        if (tracer->meshIndex(timeIdx).isOnPole()) {
            tracer->x(timeIdx).transformToPS(domain);
        }
    }
    REPORT_NOTICE(tracers.size() << " tracers are initialized.");
}
#endif

void TracerManager::registerTracer(const string &name, const string &units,
                                   const string &brief, bool smooth) {
    speciesInfos.push_back(new TracerSpeciesInfo(name, units, brief, smooth));
    for (int t = 0; t < tracers.size(); ++t) {
        tracers[t]->addSpecies();
    }
    REPORT_NOTICE("\"" << name << "\" is registered.");
}

int TracerManager::speciesIndex(const string &name) const {
    for (int i = 0; i < speciesInfos.size(); ++i) {
        if (speciesInfos[i]->name() == name) {
            return i;
        }
    }
    REPORT_ERROR("Unregistered tracer species \"" << name << "\"!");
}
    
int TracerManager::numSpecies() const {
    return speciesInfos.size();
}
    
const TracerSpeciesInfo& TracerManager::speciesInfo(int speciesIdx) const {
    return *speciesInfos[speciesIdx];
}

void TracerManager::resetSpecies() {
    for (int t = 0; t < tracers.size(); ++t) {
        tracers[t]->resetSpecies();
    }
}

void TracerManager::input(const string &fileName) {
    TimeLevelIndex<2> timeIdx;
    int numTracerDimId, numSkel1DimId, numDimDimId, numSpeciesDimId;
    int idVarId, cVarId, rhoVarId, mVarId, s1VarId;
    int ncId, ret;
    size_t numTracer, numDim, _numSpecies, numSkel1;
    int l;
    int *intData;
    double *doubleData;

    ret = nc_open(fileName.c_str(), NC_NOWRITE, &ncId);
    CHECK_NC_OPEN(ret, fileName);

    ret = nc_inq_dimid(ncId, "num_tracer", &numTracerDimId);
    CHECK_NC_INQ_DIMID(ret, fileName, "num_tracer");

    ret = nc_inq_dimlen(ncId, numTracerDimId, &numTracer);
    CHECK_NC_INQ_DIMLEN(ret, fileName, "num_tracer");

    ret = nc_inq_dimid(ncId, "num_dim", &numDimDimId);
    CHECK_NC_INQ_DIMID(ret, fileName, "num_dim");

    ret = nc_inq_dimlen(ncId, numDimDimId, &numDim);
    CHECK_NC_INQ_DIMLEN(ret, fileName, "num_dim");

    ret = nc_inq_dimid(ncId, "num_species", &numSpeciesDimId);
    CHECK_NC_INQ_DIMID(ret, fileName, "num_species");

    ret = nc_inq_dimlen(ncId, numSpeciesDimId, &_numSpecies);
    CHECK_NC_INQ_DIMLEN(ret, fileName, "num_species");

    ret = nc_inq_dimid(ncId, "num_skel1", &numSkel1DimId);
    CHECK_NC_INQ_DIMID(ret, fileName, "num_skel1");

    ret = nc_inq_dimlen(ncId, numSkel1DimId, &numSkel1);
    CHECK_NC_INQ_DIMLEN(ret, fileName, "num_skel1");

    ret = nc_inq_varid(ncId, "id", &idVarId);
    CHECK_NC_INQ_VARID(ret, fileName, "id");

    ret = nc_inq_varid(ncId, "c", &cVarId);
    CHECK_NC_INQ_VARID(ret, fileName, "c");

    ret = nc_inq_varid(ncId, "rho", &rhoVarId);
    CHECK_NC_INQ_VARID(ret, fileName, "rho");

    ret = nc_inq_varid(ncId, "m", &mVarId);
    CHECK_NC_INQ_VARID(ret, fileName, "m");

    ret = nc_inq_varid(ncId, "s1", &s1VarId);
    CHECK_NC_INQ_VARID(ret, fileName, "s1");

    for (int t = 0; t < tracers.size(); ++t) {
        delete tracers[t];
    }
    tracers.resize(numTracer);
    for (int t = 0; t < numTracer; ++t) {
        tracers[t] = new Tracer(numDim);
    }

    intData = new int[numTracer];
    ret = nc_get_var(ncId, idVarId, intData);
    CHECK_NC_GET_VAR(ret, fileName, "id");
    for (int t = 0; t < numTracer; ++t) {
        tracers[t]->ID() = intData[t];
    }
    delete [] intData;

    doubleData = new double[numTracer*numDim];
    ret = nc_get_var(ncId, cVarId, doubleData);
    CHECK_NC_GET_VAR(ret, fileName, "c");
    l = 0;
    if (numDim == 2) {
        for (int t = 0; t < numTracer; ++t) {
            tracers[t]->x(timeIdx).setCoord(doubleData[l], doubleData[l+1]);
            l += 2;
        }
    } else if (numDim == 3) {
        for (int t = 0; t < numTracer; ++t) {
            tracers[t]->x(timeIdx).setCoord(doubleData[l], doubleData[l+1],
                                               doubleData[l+2]);
            l += 3;
        }
    }
    delete [] doubleData;

    assert(numSpecies() == _numSpecies);
    doubleData = new double[numTracer*_numSpecies];
    ret = nc_get_var(ncId, rhoVarId, doubleData);
    CHECK_NC_GET_VAR(ret, fileName, "rho");
    l = 0;
    for (int t = 0; t < numTracer; ++t) {
        for (int s = 0; s < _numSpecies; ++s) {
            tracers[t]->addSpecies();
        }
        for (int s = 0; s < _numSpecies; ++s) {
            tracers[t]->density(s) = doubleData[l++];
        }
    }
    ret = nc_get_var(ncId, mVarId, doubleData);
    CHECK_NC_GET_VAR(ret, fileName, "m");
    l = 0;
    for (int t = 0; t < numTracer; ++t) {
        for (int s = 0; s < _numSpecies; ++s) {
            tracers[t]->mass(s) = doubleData[l++];
        }
    }
    delete [] doubleData;

    if (numDim == 2) {
        doubleData = new double[numTracer*4*numDim];
        ret = nc_get_var(ncId, s1VarId, doubleData);
        CHECK_NC_GET_VAR(ret, fileName, "s1");
        l = 0;
        for (int t = 0; t < numTracer; ++t) {
            TracerSkeleton &s = tracers[t]->skeleton();
            vector<SpaceCoord*> &xs = s.spaceCoords(timeIdx);
            for (int i = 0; i < xs.size(); ++i) {
                xs[i]->setCoord(doubleData[l], doubleData[l+1]);
                l += 2;
            }
        }
    } else {
        REPORT_ERROR("Under construction!");
    }
    delete [] doubleData;
    // Do the rest initialization jobs.
    for (int t = 0; t < tracers.size(); ++t) {
        Tracer *tracer = tracers[t];
        tracer->x(timeIdx).transformToCart(*domain);
        tracer->meshIndex(timeIdx).locate(*mesh, tracer->x(timeIdx));
        // TODO: This may be unnecessary.
        // when tracer is on Pole, transform its coordinate to PS for later use
        if (tracer->meshIndex(timeIdx).isOnPole()) {
            tracer->x(timeIdx).transformToPS(*domain);
        }
        TracerSkeleton &s = tracer->skeleton();
        vector<SpaceCoord*> &xs = s.spaceCoords(timeIdx);
        vector<MeshIndex*> &idxs = s.meshIndices(timeIdx);
        for (int i = 0; i < xs.size(); ++i) {
            xs[i]->transformToCart(*domain);
            idxs[i]->locate(*mesh, *xs[i]);
        }
        tracer->updateDeformMatrix(*domain, *mesh, timeIdx);
    }

    ret = nc_close(ncId);
    CHECK_NC_CLOSE(ret, fileName);
}

void TracerManager::output(const TimeLevelIndex<2> &timeIdx, int ncId) {
    int numTracerDimId, numSkel1DimId, numDimDimId, numSpeciesDimId;
    int idVarId;
    int cDimIds[2], cVarId;
    int rhoDimIds[2], rhoVarId, mVarId;
    int sDimIds[3], s1VarId;
#define OUTPUT_TRACER_SHAPE
#ifdef OUTPUT_TRACER_SHAPE
    int numSkel2DimId, s2VarId, numSkel2 = 40;
#endif
    char str[100];
    int l;
    int *intData;
    double *doubleData;

    nc_redef(ncId);

    nc_def_dim(ncId, "num_tracer", tracers.size(), &numTracerDimId);
    nc_def_dim(ncId, "num_dim", domain->numDim(), &numDimDimId);
    nc_def_dim(ncId, "num_species", numSpecies(), &numSpeciesDimId);

    if (domain->numDim() == 2) {
        // only output skeleton in 2D domain, since in 3D it could be messy.
        nc_def_dim(ncId, "num_skel1", 4, &numSkel1DimId);
#ifdef OUTPUT_TRACER_SHAPE
        nc_def_dim(ncId, "num_skel2", numSkel2, &numSkel2DimId);
#endif
    } else {
        REPORT_ERROR("Under construction!");
    }

    time_t curr_time;
    time(&curr_time);
    struct tm *timeinfo;
    timeinfo = gmtime(&curr_time);
    sprintf(str, "UTC %4.2d-%2.2d-%2.2d",
            timeinfo->tm_year+1900,
            timeinfo->tm_mon+1,
            timeinfo->tm_mday);
    nc_put_att(ncId, NC_GLOBAL, "create_date", NC_CHAR, strlen(str), str);
    nc_put_att(ncId, NC_GLOBAL, "scale0", NC_DOUBLE, 1, &scale0);
    
    nc_def_var(ncId, "id", NC_INT, 1, &numTracerDimId, &idVarId);
    sprintf(str, "tracer identifier");
    nc_put_att(ncId, idVarId, "long_name", NC_CHAR, strlen(str), str);

    cDimIds[0] = numTracerDimId;
    cDimIds[1] = numDimDimId;
    nc_def_var(ncId, "c", NC_DOUBLE, 2, cDimIds, &cVarId);
    sprintf(str, "tracer centroid coordinates on %s", domain->brief().c_str());
    nc_put_att(ncId, cVarId, "long_name", NC_CHAR, strlen(str), str);
    
    rhoDimIds[0] = numTracerDimId;
    rhoDimIds[1] = numSpeciesDimId;
    nc_def_var(ncId, "rho", NC_DOUBLE, 2, rhoDimIds, &rhoVarId);
    sprintf(str, "tracer species denisty");
    nc_put_att(ncId, rhoVarId, "long_name", NC_CHAR, strlen(str), str);

    nc_def_var(ncId, "m", NC_DOUBLE, 2, rhoDimIds, &mVarId);
    sprintf(str, "tracer species mass");
    nc_put_att(ncId, mVarId, "long_name", NC_CHAR, strlen(str), str);

    if (domain->numDim() == 2) {
        sDimIds[0] = numTracerDimId;
        sDimIds[1] = numSkel1DimId;
        sDimIds[2] = numDimDimId;
        nc_def_var(ncId, "s1", NC_DOUBLE, 3, sDimIds, &s1VarId);
        sprintf(str, "tracer actual skeleton");
        nc_put_att(ncId, s1VarId, "long_name", NC_CHAR, strlen(str), str);
#ifdef OUTPUT_TRACER_SHAPE
        sDimIds[1] = numSkel2DimId;
        nc_def_var(ncId, "s2", NC_DOUBLE, 3, sDimIds, &s2VarId);
        sprintf(str, "tracer fitted skeleton");
        nc_put_att(ncId, s2VarId, "long_name", NC_CHAR, strlen(str), str);
#endif
    }
    
    nc_enddef(ncId);

    intData = new int[tracers.size()];
    l = 0;
    for (int t = 0; t < tracers.size(); ++t) {
        intData[l++] = tracers[t]->ID();
    }
    nc_put_var(ncId, idVarId, intData);
    delete [] intData;
    
    doubleData = new double[tracers.size()*domain->numDim()];
    l = 0;
    for (int t = 0; t < tracers.size(); ++t) {
        for (int m = 0; m < domain->numDim(); ++m) {
            doubleData[l++] = tracers[t]->x(timeIdx)(m);
        }
    }
    nc_put_var(ncId, cVarId, doubleData);
    delete [] doubleData;
    
    doubleData = new double[tracers.size()*numSpecies()];
    l = 0;
    for (int t = 0; t < tracers.size(); ++t) {
        for (int s = 0; s < numSpecies(); ++s) {
            doubleData[l++] = tracers[t]->density(s);
        }
    }
    nc_put_var(ncId, rhoVarId, doubleData);
    l = 0;
    for (int t = 0; t < tracers.size(); ++t) {
        for (int s = 0; s < numSpecies(); ++s) {
            doubleData[l++] = tracers[t]->mass(s);
        }
    }
    nc_put_var(ncId, mVarId, doubleData);
    delete [] doubleData;

    if (domain->numDim() == 2) {
        doubleData = new double[tracers.size()*4*domain->numDim()];
        l = 0;
        for (int t = 0; t < tracers.size(); ++t) {
            TracerSkeleton &s = tracers[t]->skeleton();
            vector<SpaceCoord*> &xs = s.spaceCoords(timeIdx);
            for (int i = 0; i < xs.size(); ++i) {
                for (int m = 0; m < domain->numDim(); ++m) {
                    doubleData[l++] = (*xs[i])(m);
                }
            }
        }
        nc_put_var(ncId, s1VarId, doubleData);
        delete [] doubleData;
#ifdef OUTPUT_TRACER_SHAPE
        double dtheta = PI2/numSkel2;
        BodyCoord y(2);
        SpaceCoord x(2);
        doubleData = new double[tracers.size()*numSkel2*domain->numDim()];
        l = 0;
        for (int t = 0; t < tracers.size(); ++t) {
            for (int i = 0; i < numSkel2; ++i) {
                double theta = i*dtheta;
                y(0) = cos(theta);
                y(1) = sin(theta);
                tracers[t]->calcSpaceCoord(*domain, timeIdx, y, x);
                for (int m = 0; m < domain->numDim(); ++m) {
                    doubleData[l++] = x(m);
                }
            }
        }
        nc_put_var(ncId, s2VarId, doubleData);
        delete [] doubleData;
#endif
    } else {
        REPORT_ERROR("Under construction!");
    }
}

}
