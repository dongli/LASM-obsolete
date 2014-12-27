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
static int numLev;

extern "C" {

void lasm_init_cpp_(double earthRadius, double stepSize,
                    int numLon, double *lonFull, double *lonHalf,
                    int numLat, double *latFull, double *latHalf,
                    int _numLev) {
    //configManager.parse(...);
    domain.radius() = earthRadius;
    dt = stepSize;
    // Initialize lat-lon mesh.
    arma::vec _lonFull(numLon), _lonHalf(numLon);
    for (int i = 0; i < numLon; ++i) {
        _lonFull[i] = lonFull[i];
        _lonHalf[i] = lonHalf[i];
    }
    mesh.setGridCoordComps(0, numLon, _lonFull, _lonHalf);
    arma::vec _latFull(numLat), _latHalf(numLat-1);
    for (int j = 0; j < numLat; ++j) {
        _latFull[j] = latFull[j];
    }
    for (int j = 0; j < numLat-1; ++j) {
        _latHalf[j] = latHalf[j];
    }
    mesh.setGridCoordComps(1, numLat, _latFull, _latHalf);
    // Initialize velocity fields on each level.
    velocityFields.resize(numLev);
    for (int k = 0; k < numLev; ++k) {
        velocityFields[k].create(mesh, true, true);
    }
    // Initialize advection managers on each level.
    numLev = _numLev;
    advectionManagers.resize(numLev);
    for (int k = 0; k < numLev; ++k) {
        advectionManagers[k].init(domain, mesh, configManager);
        advectionManagers[k].registerTracer("Q", "?", "water vapor");
    }
}

void lasm_advance_cpp_(double *u, double *v, double *q) {
    // Copy external velocity arrays into internal velocity objects.
    int l = 0;
    for (int k = 0; k < numLev; ++k) {
        for (int j = mesh.js(geomtk::RLLStagger::GridType::FULL); j <= mesh.je(geomtk::RLLStagger::GridType::FULL); ++j) {
            for (int i = mesh.is(geomtk::RLLStagger::GridType::HALF); i <= mesh.ie(geomtk::RLLStagger::GridType::HALF); ++i) {
                velocityFields[k](0)(newTimeIdx, i, j) = u[l++];
            }
        }
    }
    l = 0;
    for (int k = 0; k < numLev; ++k) {
        for (int j = mesh.js(geomtk::RLLStagger::GridType::HALF); j <= mesh.je(geomtk::RLLStagger::GridType::HALF); ++j) {
            for (int i = mesh.is(geomtk::RLLStagger::GridType::FULL); i <= mesh.ie(geomtk::RLLStagger::GridType::FULL); ++i) {
                velocityFields[k](1)(newTimeIdx, i, j) = v[l++];
            }
        }
    }
    // Calculate external tracer density tendencies. Input 'q' should include
    // those tendencies.
    l = 0;
    for (int k = 0; k < numLev; ++k) {
        for (int i = 0; i < mesh.totalNumGrid(geomtk::RLLStagger::Location::CENTER); ++i) {
            for (int s = 0; s < advectionManagers[k].numSpecies(); ++s) {
                (*advectionManagers[k].tendency()[s])(i) =
                    (q[l++]-(*advectionManagers[k].density()[s])(oldTimeIdx, i))/dt;
            }
        }
    }
    // Advance advection managers on each level.
    for (int k = 0; k < numLev; ++k) {
        advectionManagers[k].advance(dt, newTimeIdx, velocityFields[k]);
    }
}

void lasm_final_cpp_() {
    
}

}

#endif // __LASM_GamilAdaptor__
