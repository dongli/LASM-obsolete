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

extern "C" {

void lasm_init_cpp_(double *earthRadius, double *stepSize, int *numLon, double *lon,
                    int *numLat, double *lat, int *numLev) {
    //configManager.parse(...);
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
}

void lasm_advance_cpp_(double *u, double *v, double *q) {
    // Copy external velocity arrays into internal velocity objects.
    int nx = mesh.numGrid(0, geomtk::RLLStagger::GridType::FULL);
    int ny = mesh.numGrid(1, geomtk::RLLStagger::GridType::FULL);
    int nz = advectionManagers.size();
    int ns = advectionManagers[0].numSpecies();
    for (int k = 0; k < nz; ++k) {
        for (int j = mesh.js(geomtk::RLLStagger::GridType::FULL)-1; j <= mesh.je(geomtk::RLLStagger::GridType::FULL)+1; ++j) {
            for (int i = mesh.is(geomtk::RLLStagger::GridType::HALF)-1; i <= mesh.ie(geomtk::RLLStagger::GridType::HALF)+1; ++i) {
                int l = k*(ny+2)*(nx+2)+j*(nx+2)+i;
                velocityFields[k](0)(newTimeIdx, i, j) = u[l];
            }
        }
    }
    for (int k = 0; k < nz; ++k) {
        for (int j = mesh.js(geomtk::RLLStagger::GridType::HALF)-1; j <= mesh.je(geomtk::RLLStagger::GridType::HALF)+1; ++j) {
            for (int i = mesh.is(geomtk::RLLStagger::GridType::FULL)-1; i <= mesh.ie(geomtk::RLLStagger::GridType::FULL)+1; ++i) {
                int l = k*(ny+2)*(nx+2)+j*(nx+2)+i;
                velocityFields[k](1)(newTimeIdx, i, j) = v[l];
            }
        }
    }
    // Calculate external tracer density tendencies. Input 'q' should include
    // those tendencies.
    for (int s = 0; s < ns; ++s) {
        for (int k = 0; k < nz; ++k) {
            for (int j = mesh.js(geomtk::RLLStagger::GridType::HALF)-1; j <= mesh.je(geomtk::RLLStagger::GridType::HALF)+1; ++j) {
                for (int i = mesh.is(geomtk::RLLStagger::GridType::FULL)-1; i <= mesh.ie(geomtk::RLLStagger::GridType::FULL)+1; ++i) {
                    int l = s*nz*(ny+2)*(nx+2)+k*(ny+2)*(nx+2)+j*(nx+2)+i;
                    (*advectionManagers[k].tendency()[s])(i, j) =
                        (q[l]-(*advectionManagers[k].density()[s])(oldTimeIdx, i, j))/dt;
                }
            }
        }
    }
    // Advance advection managers on each level.
    for (int k = 0; k < nz; ++k) {
        advectionManagers[k].advance(dt, newTimeIdx, velocityFields[k]);
    }
}

void lasm_final_cpp_() {
    
}

}

#endif // __LASM_GamilAdaptor__
