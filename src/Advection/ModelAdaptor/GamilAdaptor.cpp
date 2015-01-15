#ifndef __LASM_GamilAdaptor__
#define __LASM_GamilAdaptor__

#include "lasm.h"

/**
 *  This is the LASM adaptor for GAMIL. Note LASM is run in 2.5-D mode.
 *
 *  2014-11-19
 */

static geomtk::SphereDomain domain(2);
static geomtk::RLLMesh mesh(domain);
static geomtk::ConfigManager configManager;
static std::vector<geomtk::RLLVelocityField> velocityFields;
static std::vector<lasm::AdvectionManager> advectionManagers;
static geomtk::TimeLevelIndex<2> oldTimeIdx, newTimeIdx;
static double dt;

static const int FULL = geomtk::RLLStagger::GridType::FULL;
static const int HALF = geomtk::RLLStagger::GridType::HALF;
static const int CENTER = geomtk::RLLStagger::Location::CENTER;

// For debug output only:
static geomtk::TimeManager timeManager;
static geomtk::IOManager<geomtk::RLLDataFile> io;
static int debugFileIdx;

extern "C" {

void debug_output() {
    io.create(debugFileIdx);
    io.output<double, 2>(debugFileIdx, oldTimeIdx,
                         {&velocityFields[25](0),
                          &velocityFields[25](1)});
                          //&velocityFields[25].divergence(),
                          //advectionManagers[25].density()[0]});
    //advectionManagers[25].output(oldTimeIdx, io.file(debugFileIdx).fileID);
    io.close(debugFileIdx);
}

void lasm_init_cpp_(double *earthRadius, double *stepSize,
                    int *numLon, double *lon,
                    int *numLat, double *lat, int *numLev,
                    int *startYear, int *startMonth,
                    int *startDay, int *startSecond,
                    int *stopYear, int *stopMonth,
                    int *stopDay, int *stopSecond) {
    REPORT_NOTICE("Initialize LASM.");
    // Configure LASM.
    configManager.addKeyValue("lasm", "num_parcel_x", *numLon);
    configManager.addKeyValue("lasm", "num_parcel_y", *numLat);
    //
    domain.radius() = *earthRadius;
    dt = *stepSize;
    // Initialize lat-lon mesh.
    arma::vec _lon(*numLon);
    for (int i = 0; i < *numLon; ++i) {
        _lon[i] = lon[i];
    }
    mesh.setGridCoordComps(0, *numLon, _lon);
    arma::vec _lat(*numLat);
    for (int j = 0; j < *numLat; ++j) {
        _lat[j] = lat[j];
    }
    mesh.setGridCoordComps(1, *numLat, _lat);
    // Initialize velocity fields on each level.
    velocityFields.resize(*numLev);
    for (int k = 0; k < *numLev; ++k) {
        velocityFields[k].create(mesh, true, true);
    }
    // Initialize advection managers on each level.
    advectionManagers.resize(*numLev);
    for (int k = 0; k < *numLev; ++k) {
        advectionManagers[k].init(domain, mesh, configManager);
        advectionManagers[k].registerTracer("Q", "?", "water vapor");
    }
    // Debug output.
    geomtk::Time startTime, stopTime;
    char str[30];
    sprintf(str, "%4.4d-%2.2d-%2.2d %5.5d", *startYear, *startMonth, *startDay, *startSecond);
    startTime = str;
    sprintf(str, "%4.4d-%2.2d-%2.2d %5.5d", *stopYear, *stopMonth, *stopDay, *stopSecond);
    stopTime = str;
    timeManager.init(startTime, stopTime, dt);
    io.init(timeManager);
    std::string debugFilePattern = "lasm.gamil.128x60.%5s.nc";
    debugFileIdx = io.registerOutputFile(mesh, debugFilePattern, geomtk::TimeStepUnit::STEP, 1);
    io.file(debugFileIdx).registerField("double", geomtk::RLLSpaceDimensions::FULL_DIMENSION,
                                        {&velocityFields[25](0), &velocityFields[25](1)});
}

void lasm_first_step_cpp_(double *u, double *v, double *q) {
    // Copy external velocity arrays into internal velocity objects.
    int nx = mesh.numGrid(0, FULL);
    int ny = mesh.numGrid(1, FULL);
    int nz = advectionManagers.size();
    int ns = advectionManagers[0].numSpecies();
    for (int k = 0; k < nz; ++k) {
        for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
            for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
                int l = k*(ny+2)*(nx+2)+(j+1)*(nx+2)+i;
                velocityFields[k](0)(oldTimeIdx, i, ny-1-j) = u[l];
            }
        }
    }
    for (int k = 0; k < nz; ++k) {
        for (int j = mesh.js(HALF); j <= mesh.je(HALF); ++j) {
            for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                int l = k*(ny+2)*(nx+2)+(j+1)*(nx+2)+i;
                velocityFields[k](1)(oldTimeIdx, i, ny-2-j) = -v[l];
            }
        }
    }
    for (int k = 0; k < nz; ++k) {
        velocityFields[k].applyBndCond(oldTimeIdx);
    }
    // Copy external tracer density arrays into internal tracer density objects.
    double *q_ = new double[ns*mesh.totalNumGrid(CENTER, 2)];
    for (int k = 0; k < nz; ++k) {
        for (int s = 0; s < ns; ++s) {
            for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
                for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                    int l1 = s*ny*nx+(ny-1-j)*nx+i-1;
                    int l2 = s*nz*(ny+2)*(nx+2)+k*(ny+2)*(nx+2)+(j+1)*(nx+2)+i;
                    q_[l1] = q[l2];
                }
            }
        }
        advectionManagers[k].input(oldTimeIdx, q_);
    }
    delete [] q_;
    // Copy out internal tracer density after remapping.
    for (int k = 0; k < nz; ++k) {
        for (int s = 0; s < ns; ++s) {
            for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
                for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                    int l = s*nz*(ny+2)*(nx+2)+k*(ny+2)*(nx+2)+(j+1)*(nx+2)+i;
                    q[l] = (*advectionManagers[k].density()[s])(newTimeIdx, i, ny-1-j);
                }
            }
        }
    }
    // Debug output.
    debug_output();
}

void lasm_advance_cpp_(double *u, double *v, double *q) {
    newTimeIdx = oldTimeIdx+1;
    // Copy external velocity arrays into internal velocity objects.
    int nx = mesh.numGrid(0, FULL);
    int ny = mesh.numGrid(1, FULL);
    int nz = advectionManagers.size();
    int ns = advectionManagers[0].numSpecies();
    for (int k = 0; k < nz; ++k) {
        for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
            for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
                int l = k*(ny+2)*(nx+2)+(j+1)*(nx+2)+i;
                velocityFields[k](0)(newTimeIdx, i, ny-1-j) = u[l];
            }
        }
    }
    for (int k = 0; k < nz; ++k) {
        for (int j = mesh.js(HALF); j <= mesh.je(HALF); ++j) {
            for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                int l = k*(ny+2)*(nx+2)+(j+1)*(nx+2)+i;
                velocityFields[k](1)(newTimeIdx, i, ny-2-j) = -v[l];
            }
        }
    }
    for (int k = 0; k < nz; ++k) {
        velocityFields[k].applyBndCond(newTimeIdx, true);
    }
    // Calculate external tracer density tendencies. Input 'q' should include
    // those tendencies.
    for (int k = 0; k < nz; ++k) {
        for (int s = 0; s < ns; ++s) {
            for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
                for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                    int l = s*nz*(ny+2)*(nx+2)+k*(ny+2)*(nx+2)+(j+1)*(nx+2)+i;
                    (*advectionManagers[k].tendency()[s])(i, ny-1-j) =
                        (q[l]-(*advectionManagers[k].density()[s])(oldTimeIdx, i, ny-1-j))/dt;
                }
            }
        }
    }
    // Advance advection managers on each level.
    std::cerr << timeManager.currTime().s() << std::endl;
    for (int k = 0; k < nz; ++k) {
        advectionManagers[k].advance(dt, newTimeIdx, velocityFields[k]);
    }
    // Copy out the advected tracer density.
    for (int k = 0; k < nz; ++k) {
        for (int s = 0; s < ns; ++s) {
            for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
                for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                    int l = s*nz*(ny+2)*(nx+2)+k*(ny+2)*(nx+2)+(j+1)*(nx+2)+i;
                    q[l] = (*advectionManagers[k].density()[s])(newTimeIdx, i, ny-1-j);
                }
            }
        }
    }
    // Shift time index.
    oldTimeIdx.shift();
    // Debug output.
    timeManager.advance();
    debug_output();
}

void lasm_final_cpp_() {
    
}

}

#endif // __LASM_GamilAdaptor__
