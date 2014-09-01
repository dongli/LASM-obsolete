#include "SolidRotationTestCase.h"

namespace lasm {

SolidRotationTestCase::SolidRotationTestCase() {
    REPORT_ONLINE;
}

SolidRotationTestCase::~SolidRotationTestCase() {
    delete mesh;
    delete domain;
    delete axisPole;
    delete c0;
    delete cr0;
    REPORT_OFFLINE;
}

void SolidRotationTestCase::init(const ConfigManager &configManager,
                                 TimeManager &timeManager) {
    // initialize domain
    domain = new geomtk::SphereDomain(2);
    domain->setRadius(1.0);
    // initialize mesh
    mesh = new geomtk::RLLMesh(*domain);
    int numLon = 240, numLat = 121;
    if (configManager.hasKey("test_case", "num_lon")) {
        configManager.getValue("test_case", "num_lon", numLon);
    }
    if (configManager.hasKey("test_case", "num_lat")) {
        configManager.getValue("test_case", "num_lat", numLat);
    }
    mesh->init(numLon, numLat);
    // initialize velocity
    velocity.create(*mesh, true, HAS_HALF_LEVEL);
    // set parameters
    angleSpeed = PI2/12/86400;
    U0 = domain->getRadius()*angleSpeed;
    alpha = M_PI_2;
    axisPole = new SpaceCoord(2);
    c0 = new SpaceCoord(2);
    cr0 = new SpaceCoord(2);
    axisPole->setCoord(M_PI,  M_PI_2-alpha);
    c0->setCoord(M_PI_2, alpha);
    domain->rotate(*axisPole, *c0, *cr0);
    R = domain->getRadius()/3;
    H0 = 1000;
    AdvectionTestCase::init(configManager, timeManager);
}

Time SolidRotationTestCase::getStartTime() const {
    Time time;
    return time;
}

Time SolidRotationTestCase::getEndTime() const {
    Time time;
    return time+12*TimeUnit::DAYS;
}

double SolidRotationTestCase::getStepSize() const {
    return 30*TimeUnit::MINUTES;
}

void SolidRotationTestCase::advance(double time,
                                    const TimeLevelIndex<2> &timeIdx) {
    double sinAlpha = sin(alpha), cosAlpha = cos(alpha);
    for (int i = 0; i < mesh->getTotalNumGrid(velocity(0).getStaggerLocation(), 2); ++i) {
        const SpaceCoord &x = mesh->getGridCoord(velocity(0).getStaggerLocation(), i);
        double cosLat = x.getCosLat();
        double sinLat = x.getSinLat();
        double cosLon = x.getCosLon();
        velocity(0)(timeIdx, i) = U0*(cosLat*cosAlpha+sinLat*cosLon*sinAlpha);
    }
    for (int i = 0; i < mesh->getTotalNumGrid(velocity(1).getStaggerLocation(), 2); ++i) {
        const SpaceCoord &x = mesh->getGridCoord(velocity(1).getStaggerLocation(), i);
        double sinLon = x.getSinLon();
        velocity(1)(timeIdx, i) = -U0*sinLon*sinAlpha;
    }
    if (timeIdx.isCurrentIndex()) {
        velocity.applyBndCond(timeIdx);
    } else {
        velocity.applyBndCond(timeIdx, UPDATE_HALF_LEVEL);
    }
}

void SolidRotationTestCase::calcInitCond(AdvectionManager &advectionManager) {
    TimeLevelIndex<2> timeIdx;
    advectionManager.registerTracer("q0", "N/A", "background tracer");
    advectionManager.registerTracer("q1", "N/A", "cosine hill tracer");
    densities = &advectionManager.getDensities();
    AdvectionTestCase::registerDefaultOutput();
    double q[2*mesh->getTotalNumGrid(CENTER, 2)];
    int l = 0;
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, 2); ++i) {
        q[l++] = 1.0;
    }
    calcSolution(0, timeIdx, *(*densities)[1]);
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, 2); ++i) {
        q[l++] = (*(*densities)[1])(timeIdx, i);
    }
    // propagate initial conditions to advection manager
    advectionManager.input(timeIdx, q);
}

void SolidRotationTestCase::calcSolution(double dt,
                                         const TimeLevelIndex<2> &timeIdx,
                                         AdvectionManager &advectionManager) {
    calcSolution(dt, timeIdx, *(*densities)[1]);
    advectionManager.remapMeshToTracers(timeIdx);
    REPORT_NOTICE("Overwrite tracers with the true solution.");
}

void SolidRotationTestCase::calcSolution(double dt,
                                         const TimeLevelIndex<2> &timeIdx,
                                         ScalarField &q) {
    cr0->setCoordComp(0, (*cr0)(0)+angleSpeed*dt);
    domain->rotateBack(*axisPole, *c0, *cr0);
    
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh->getGridCoord(CENTER, i);
        double d = domain->calcDistance(*c0, x);
        if (d < R) {
            q(timeIdx, i) = H0*(1+cos(M_PI*d/R))/2;
        } else {
            q(timeIdx, i) = 0;
        }
    }
}

}
