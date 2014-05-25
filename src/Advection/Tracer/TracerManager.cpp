#include "TracerManager.h"
#include "TracerSkeleton.h"

namespace lasm
{

TracerManager::TracerManager() {
    scale0 = 1.5;
    REPORT_ONLINE;
}

TracerManager::~TracerManager() {
    REPORT_OFFLINE;
}

void TracerManager::init(const Domain &domain, const Mesh &mesh,
                         const geomtk::ConfigManager &configManager) {
    this->domain = &domain;
    int numParcel = 0;
    int numParcelX, numParcelY;
    configManager.getValue("lasm", "num_parcel_x", numParcelX);
    configManager.getValue("lasm", "num_parcel_y", numParcelY);
    if (configManager.hasKey("lasm", "gammas")) {
        configManager.getValue("lasm", "gammas", Tracer::gammas);
    }
    if (configManager.hasKey("lasm", "scale0")) {
        configManager.getValue("lasm", "scale0", scale0);
    }
    // calculate the total number of tracers
//#define USE_FULL_LAT_LON
#ifdef LASM_USE_SPHERE_DOMAIN
#ifdef USE_FULL_LAT_LON
    numParcel = numParcelX*numParcelY;
#else
    // Note: Use reduced lat-lon mesh as the initial distribution of tracers.
    if (numParcelX%2 != 0) {
        REPORT_ERROR("Tracer number (now is " << numParcelX <<
                     ") along longitude axis should be even!");
    }
    double dlat = domain.getAxisSpan(1)/numParcelY;
    double shiftLat = 45.0*RAD;
    double minNumTracerX = 4;
    double cosLat0, cosLat1 = cos(shiftLat);
    for (int j = 0; j < numParcelY; ++j) {
        double lat = domain.getAxisStart(1)+dlat*0.5+dlat*j;
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
#else
    REPORT_ERROR("Under construction!");
#endif
    tracers.resize(numParcel);
    int ID = 0;
    list<Tracer*>::iterator tracer = tracers.begin();
    for (; tracer != tracers.end(); ++tracer) {
        *tracer = new Tracer(domain.getNumDim());
        (*tracer)->setID(ID++);
    }
    TimeLevelIndex<2> initTimeIdx;
    for (tracer = tracers.begin(); tracer != tracers.end(); ++tracer) {
        // set coordinate
        SpaceCoord &x0 = (*tracer)->getX(initTimeIdx);
        vec h(domain.getNumDim());
#ifdef LASM_USE_SPHERE_DOMAIN
#ifdef USE_FULL_LAT_LON
        double dlon = domain.getAxisSpan(0)/numParcelX;
        double dlat = domain.getAxisSpan(1)/numParcelY;
        int l = 0;
        for (int j = 0; j < numParcelY; ++j) {
            for (int i = 0; i < numParcelX; ++i) {
                if (l == (*tracer)->getID()) {
                    double lat = domain.getAxisStart(1)+dlat*0.5+dlat*j;
                    double lon = domain.getAxisStart(0)+dlon*0.5+dlon*i;
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
        double dlat = domain.getAxisSpan(1)/numParcelY;
        double shiftLat = 45.0*RAD;
        double minNumTracerX = 4;
        double cosLat0, cosLat1 = cos(shiftLat);
        double dlon;
        int l = 0;
        for (int j = 0; j < numParcelY; ++j) {
            double lat = domain.getAxisStart(1)+dlat*0.5+dlat*j;
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
                if (l == (*tracer)->getID()) {
                    dlon = domain.getAxisSpan(0)/numParcelX_;
                    double lon = domain.getAxisStart(0)+dlon*0.5+dlon*i;
                    x0.setCoord(lon, lat);
                    l = -1;
                    break;
                }
                l++;
            }
            if (l == -1) break;
        }
#endif
        x0.transformToCart(domain);
        h(0) = dlon*domain.getRadius()*x0.getCosLat();
        h(1) = dlat*domain.getRadius();
        // When h is small, there may be small spots on the grid densities.
        // The selection of h can be tricky. When h *= 1.5, the errors of solid
        // rotation test is oscillatory with time.
        h *= scale0;
#else
        REPORT_ERROR("Under construction!");
#endif
        // set mesh index
        MeshIndex &idx0 = (*tracer)->getMeshIndex(initTimeIdx);
        idx0.locate(mesh, x0);
        // TODO: This may be unnecessary.
        // when tracer is on Pole, transform its coordinate to PS for later use
        if (idx0.isOnPole()) {
            x0.transformToPS(domain);
        }
        // set tracer skeleton
        (*tracer)->getSkeleton().init(domain, mesh, h.max());
        // set deformation matrix
        (*tracer)->getH(initTimeIdx).eye();
        (*tracer)->updateDeformMatrix(domain, mesh, initTimeIdx);
    }
    // -------------------------------------------------------------------------
    REPORT_NOTICE(numParcel << " tracers are initialized.");
}

void TracerManager::registerTracer(const string &name, const string &units,
                                   const string &brief) {
    speciesInfos.push_back(new TracerSpeciesInfo(name, units, brief));
    list<Tracer*>::iterator tracer;
    for (tracer = tracers.begin(); tracer != tracers.end(); ++tracer) {
        (*tracer)->addSpecies();
    }
    REPORT_NOTICE("\"" << name << "\" is registered.");
}

int TracerManager::getSpeciesIndex(const string &name) const {
    for (int i = 0; i < speciesInfos.size(); ++i) {
        if (speciesInfos[i]->getName() == name) {
            return i;
        }
    }
    REPORT_ERROR("Unregistered tracer species \"" << name << "\"!");
}
    
int TracerManager::getNumSpecies() const {
    return speciesInfos.size();
}
    
const TracerSpeciesInfo& TracerManager::getSpeciesInfo(int speciesIdx) const {
    return *speciesInfos[speciesIdx];
}

void TracerManager::output(const string &fileName,
                           const TimeLevelIndex<2> &oldTimeIdx) {
    int ncId;
    int numTracerDimId, numSkel1DimId, numDimDimId, numSpeciesDimId;
    int idVarId;
    int cDimIds[2], cVarId;
    int hDimIds[3], hVarId;
    int rhoDimIds[2], rhoVarId;
    int sDimIds[3], s1VarId;
#define OUTPUT_TRACER_SHAPE
#ifdef OUTPUT_TRACER_SHAPE
    int numSkel2DimId, s2VarId, numSkel2 = 40;
#endif
    char str[100];
    int l;
    int *intData;
    double *doubleData;
    list<Tracer*>::iterator tracer;

    if (nc_create(fileName.c_str(), NC_CLOBBER, &ncId) != NC_NOERR) {
        REPORT_ERROR("Failed to open \"" << fileName << "\"!");
    }

    if (nc_def_dim(ncId, "num_tracer", tracers.size(), &numTracerDimId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define dimension \"num_tracer\"!");
    }

    if (domain->getNumDim() == 2) {
        // only output skeleton in 2D domain, since in 3D it could be messy.
        if (nc_def_dim(ncId, "num_skel1", 4, &numSkel1DimId)
            != NC_NOERR) {
            REPORT_ERROR("Failed to define dimension \"num_skel1\"!");
        }
#ifdef OUTPUT_TRACER_SHAPE
        if (nc_def_dim(ncId, "num_skel2", numSkel2, &numSkel2DimId)
             != NC_NOERR) {
            REPORT_ERROR("Failed to define dimension \"num_skel2\"!");
        }
#endif
    }

    if (nc_def_dim(ncId, "num_dim", domain->getNumDim(), &numDimDimId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define dimension \"num_dim\"!");
    }
    
    if (nc_def_dim(ncId, "num_species", getNumSpecies(), &numSpeciesDimId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define dimension \"num_species\"!");
    }

    time_t curr_time;
    time(&curr_time);
    struct tm *timeinfo;
    timeinfo = gmtime(&curr_time);
    sprintf(str, "UTC %4.2d-%2.2d-%2.2d",
            timeinfo->tm_year+1900,
            timeinfo->tm_mon+1,
            timeinfo->tm_mday);
    if (nc_put_att(ncId, NC_GLOBAL, "create_date", NC_CHAR, strlen(str), str)
        != NC_NOERR) {
        REPORT_ERROR("Failed to put attribute in \"" << fileName << "\"!");
    }

    if (nc_put_att(ncId, NC_GLOBAL, "scale0", NC_DOUBLE, 1, &scale0)
        != NC_NOERR) {
        REPORT_ERROR("Failed to put global dimension \"scale0\"!");
    }

    if (nc_put_att(ncId, NC_GLOBAL, "gammas", NC_DOUBLE, 1, &Tracer::gammas)
        != NC_NOERR) {
        REPORT_ERROR("Failed to put global dimension \"gammas\"!");
    }
    
    if (nc_def_var(ncId, "id", NC_INT, 1, &numTracerDimId, &idVarId) != NC_NOERR) {
        REPORT_ERROR("Failed to define variable \"id\"!");
    }

    cDimIds[0] = numTracerDimId;
    cDimIds[1] = numDimDimId;
    if (nc_def_var(ncId, "c", NC_DOUBLE, 2, cDimIds, &cVarId) != NC_NOERR) {
        REPORT_ERROR("Failed to define variable \"c\"!");
    }
    sprintf(str, "tracer centroid coordinates on %s", domain->getBrief().c_str());
    if (nc_put_att(ncId, cVarId, "long_name", NC_CHAR, strlen(str), str)
        != NC_NOERR) {
        REPORT_ERROR("Failed to put attribute in \"" << fileName << "\"!");
    }

    hDimIds[0] = numTracerDimId;
    hDimIds[1] = numDimDimId;
    hDimIds[2] = numDimDimId;
    if (nc_def_var(ncId, "h", NC_DOUBLE, 3, hDimIds, &hVarId) != NC_NOERR) {
        REPORT_ERROR("Failed to define variable \"h\"!");
    }
    sprintf(str, "tracer deformation matrix");
    if (nc_put_att(ncId, hVarId, "long_name", NC_CHAR, strlen(str), str)
        != NC_NOERR) {
        REPORT_ERROR("Failed to put attribute in \"" << fileName << "\"!");
    }
    
    rhoDimIds[0] = numTracerDimId;
    rhoDimIds[1] = numSpeciesDimId;
    if (nc_def_var(ncId, "rho", NC_DOUBLE, 2, rhoDimIds, &rhoVarId) != NC_NOERR) {
        REPORT_ERROR("Failed to define variable \"rho\"!");
    }
    sprintf(str, "tracer species denisty");
    if (nc_put_att(ncId, rhoVarId, "long_name", NC_CHAR, strlen(str), str)
        != NC_NOERR) {
        REPORT_ERROR("Failed to put attribute in \"" << fileName << "\"!");
    }

    if (domain->getNumDim() == 2) {
        sDimIds[0] = numTracerDimId;
        sDimIds[1] = numSkel1DimId;
        sDimIds[2] = numDimDimId;
        if (nc_def_var(ncId, "s1", NC_DOUBLE, 3, sDimIds, &s1VarId)
            != NC_NOERR) {
            REPORT_ERROR("Failed to define variable \"s1\"!");
        }
        sprintf(str, "tracer actual skeleton");
        if (nc_put_att(ncId, s1VarId, "long_name", NC_CHAR, strlen(str), str)
            != NC_NOERR) {
            REPORT_ERROR("Failed to put attribute in \"" << fileName << "\"!");
        }
#ifdef OUTPUT_TRACER_SHAPE
        sDimIds[1] = numSkel2DimId;
        if (nc_def_var(ncId, "s2", NC_DOUBLE, 3, sDimIds, &s2VarId)
            != NC_NOERR) {
            REPORT_ERROR("Failed to define variable \"s2\"!");
        }
        sprintf(str, "tracer fitted skeleton");
        if (nc_put_att(ncId, s2VarId, "long_name", NC_CHAR, strlen(str), str)
            != NC_NOERR) {
            REPORT_ERROR("Failed to put attribute in \"" << fileName << "\"!");
        }
#endif
    }
    
    if (nc_enddef(ncId) != NC_NOERR) {

    }
    // -------------------------------------------------------------------------
    intData = new int[tracers.size()];
    l = 0;
    for (tracer = tracers.begin(); tracer != tracers.end(); ++tracer) {
        intData[l++] = (*tracer)->getID();
    }
    if (nc_put_var(ncId, idVarId, intData) != NC_NOERR) {
        REPORT_ERROR("Failed to put variable in \"" << fileName << "\"!");
    }
    delete [] intData;
    
    doubleData = new double[tracers.size()*domain->getNumDim()];
    l = 0;
    for (tracer = tracers.begin(); tracer != tracers.end(); ++tracer) {
        for (int m = 0; m < domain->getNumDim(); ++m) {
            doubleData[l++] = (*tracer)->getX(oldTimeIdx)(m);
        }
    }
    if (nc_put_var(ncId, cVarId, doubleData) != NC_NOERR) {
        REPORT_ERROR("Failed to put variable in \"" << fileName << "\"!");
    }
    delete [] doubleData;

    doubleData = new double[tracers.size()*domain->getNumDim()*domain->getNumDim()];
    l = 0;
    for (tracer = tracers.begin(); tracer != tracers.end(); ++tracer) {
        for (int m1 = 0; m1 < domain->getNumDim(); ++m1) {
            for (int m2 = 0; m2 < domain->getNumDim(); ++m2) {
                doubleData[l++] = (*tracer)->getH(oldTimeIdx)(m1, m2);
            }
        }
    }
    if (nc_put_var(ncId, hVarId, doubleData) != NC_NOERR) {
        REPORT_ERROR("Failed to put variable in \"" << fileName << "\"!");
    }
    delete [] doubleData;
    
    doubleData = new double[tracers.size()*getNumSpecies()];
    l = 0;
    for (tracer = tracers.begin(); tracer != tracers.end(); ++tracer) {
        for (int s = 0; s < getNumSpecies(); ++s) {
            doubleData[l++] = (*tracer)->getSpeciesDensity(s);
        }
    }
    if (nc_put_var(ncId, rhoVarId, doubleData) != NC_NOERR) {
        REPORT_ERROR("Failed to put variable in \"" << fileName << "\"!");
    }
    delete [] doubleData;

    if (domain->getNumDim() == 2) {
        doubleData = new double[tracers.size()*4*domain->getNumDim()];
        l = 0;
        for (tracer = tracers.begin(); tracer != tracers.end(); ++tracer) {
            TracerSkeleton &s = (*tracer)->getSkeleton();
            vector<SpaceCoord*> &xs = s.getSpaceCoords(oldTimeIdx);
            for (int i = 0; i < xs.size(); ++i) {
                for (int m = 0; m < domain->getNumDim(); ++m) {
                    doubleData[l++] = (*xs[i])(m);
                }
            }
        }
        if (nc_put_var(ncId, s1VarId, doubleData) != NC_NOERR) {
            REPORT_ERROR("Failed to put variable in \"" << fileName << "\"!");
        }
        delete [] doubleData;
#ifdef OUTPUT_TRACER_SHAPE
        double dtheta = PI2/numSkel2;
        BodyCoord y(2);
        SpaceCoord x(2);
        doubleData = new double[tracers.size()*numSkel2*domain->getNumDim()];
        l = 0;
        for (tracer = tracers.begin(); tracer != tracers.end(); ++tracer) {
            for (int i = 0; i < numSkel2; ++i) {
                double theta = i*dtheta;
                y(0) = cos(theta);
                y(1) = sin(theta);
                (*tracer)->getSpaceCoord(*domain, oldTimeIdx, y, x);
                for (int m = 0; m < domain->getNumDim(); ++m) {
                    doubleData[l++] = x(m);
                }
            }
        }
        if (nc_put_var(ncId, s2VarId, doubleData) != NC_NOERR) {
            REPORT_ERROR("Failed to put variable in \"" << fileName << "\"!");
        }
        delete [] doubleData;
#endif
    }
    // -------------------------------------------------------------------------
    if (nc_close(ncId) != NC_NOERR) {
        REPORT_ERROR("Failed to close file \"" << fileName << "\"!");
    }
    REPORT_NOTICE("File \"" << fileName << "\" is outputted.");
}

}
