#include "SolidRotationTestCase.h"

#ifdef USE_SPHERE_DOMAIN

namespace lasm {

SolidRotationTestCase::SolidRotationTestCase() {
    REPORT_ONLINE;
}

SolidRotationTestCase::~SolidRotationTestCase() {
    delete _mesh;
    delete _domain;
    delete axisPole;
    delete c0;
    delete cr0;
    REPORT_OFFLINE;
}

void SolidRotationTestCase::init(const ConfigManager &configManager,
                                 TimeManager &timeManager) {
    // initialize domain
    _domain = new geomtk::SphereDomain(2);
    _domain->radius() = 1;
    // initialize mesh
    _mesh = new Mesh(*_domain);
    int numLon = 240, numLat = 121;
    if (configManager.hasKey("test_case", "num_lon")) {
        configManager.getValue("test_case", "num_lon", numLon);
    }
    if (configManager.hasKey("test_case", "num_lat")) {
        configManager.getValue("test_case", "num_lat", numLat);
    }
    _mesh->init(numLon, numLat);
    // initialize velocity
    velocity.create(mesh(), true, HAS_HALF_LEVEL);
    // set parameters
    angleSpeed = PI2/12/86400;
    U0 = domain().radius()*angleSpeed;
    alpha = M_PI_2;
    axisPole = new SpaceCoord(2);
    c0 = new SpaceCoord(2);
    cr0 = new SpaceCoord(2);
    axisPole->setCoord(M_PI,  M_PI_2-alpha);
    c0->setCoord(M_PI_2, alpha);
    domain().rotate(*axisPole, *c0, *cr0);
    R = domain().radius()/3;
    H0 = 1000;
    AdvectionTestCase::init(configManager, timeManager);
}

Time SolidRotationTestCase::startTime() const {
    Time time;
    return time;
}

Time SolidRotationTestCase::endTime() const {
    Time time;
    return time+12*TimeUnit::DAYS;
}

double SolidRotationTestCase::stepSize() const {
    return 30*TimeUnit::MINUTES;
}

void SolidRotationTestCase::advance(double time,
                                    const TimeLevelIndex<2> &timeIdx) {
    double sinAlpha = sin(alpha), cosAlpha = cos(alpha);
    for (int i = 0; i < mesh().totalNumGrid(velocity(0).staggerLocation(), 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(velocity(0).staggerLocation(), i);
        double cosLat = x.cosLat();
        double sinLat = x.sinLat();
        double cosLon = x.cosLon();
        velocity(0)(timeIdx, i) = U0*(cosLat*cosAlpha+sinLat*cosLon*sinAlpha);
    }
    for (int i = 0; i < mesh().totalNumGrid(velocity(1).staggerLocation(), 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(velocity(1).staggerLocation(), i);
        double sinLon = x.sinLon();
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
    density = &advectionManager.density();
    AdvectionTestCase::registerDefaultOutput();
    double q[2*mesh().totalNumGrid(CENTER, 2)];
    int l = 0;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        q[l++] = 1.0;
    }
    calcSolution(0, timeIdx, *(*density)[1]);
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        q[l++] = (*(*density)[1])(timeIdx, i);
    }
    // propagate initial conditions to advection manager
    advectionManager.input(timeIdx, q);
}

void SolidRotationTestCase::calcSolution(double dt,
                                         const TimeLevelIndex<2> &timeIdx,
                                         AdvectionManager &advectionManager) {
    calcSolution(dt, timeIdx, *(*density)[1]);
    advectionManager.remapMeshToTracers(timeIdx);
    REPORT_NOTICE("Overwrite tracers with the true solution.");
}

void SolidRotationTestCase::calcSolution(double dt,
                                         const TimeLevelIndex<2> &timeIdx,
                                         ScalarField &q) {
    cr0->setCoordComp(0, (*cr0)(0)+angleSpeed*dt);
    domain().rotateBack(*axisPole, *c0, *cr0);
    
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double d = domain().calcDistance(*c0, x);
        if (d < R) {
            q(timeIdx, i) = H0*(1+cos(M_PI*d/R))/2;
        } else {
            q(timeIdx, i) = 0;
        }
    }
}

} // lasm

#endif
