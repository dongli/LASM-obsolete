#include "CartesianRotationTestCase.h"

#ifdef LASM_CARTESIAN_DOMAIN

namespace lasm {

CartesianRotationTestCase::CartesianRotationTestCase() {
    angleSpeed = 1.0*RAD;
    _stepSize = 1;
    REPORT_ONLINE;
}

CartesianRotationTestCase::~CartesianRotationTestCase() {
    REPORT_OFFLINE;
}

void CartesianRotationTestCase::init(const ConfigManager &configManager,
                                     TimeManager &timeManager) {
    _domain = new geomtk::CartesianDomain(2);
    _domain->setAxis(0, "x", "x axis", "m", -1, geomtk::OPEN, 1, geomtk::OPEN);
    _domain->setAxis(1, "y", "y axis", "m", -1, geomtk::OPEN, 1, geomtk::OPEN);
    _mesh = new geomtk::CartesianMesh(*_domain);
    int nx = 100, ny = 100;
    if (configManager.hasKey("test_case", "num_x")) {
        configManager.getValue("test_case", "num_x", nx);
    }
    if (configManager.hasKey("test_case", "num_y")) {
        configManager.getValue("test_case", "num_y", ny);
    }
    _mesh->init(nx, ny);
    velocity.create(*_mesh, true, HAS_HALF_LEVEL);
    AdvectionTestCase::init(configManager, timeManager);
}

Time CartesianRotationTestCase::startTime() const {
    Time res;
    return res;
}

Time CartesianRotationTestCase::endTime() const {
    Time res = startTime()+360;
    return res;
}

void CartesianRotationTestCase::advance(double time,
                                        const TimeLevelIndex<2> &timeIdx) {
    for (int i = 0; i < mesh().totalNumGrid(velocity(0).staggerLocation(), 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(velocity(0).staggerLocation(), i);
        double theta = atan2(x(1), x(0));
        velocity(0)(timeIdx, i) = -sin(theta)*norm(x())*angleSpeed;
    }
    for (int i = 0; i < mesh().totalNumGrid(velocity(1).staggerLocation(), 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(velocity(1).staggerLocation(), i);
        double theta = atan2(x(1), x(0));
        velocity(1)(timeIdx, i) = cos(theta)*norm(x())*angleSpeed;
    }
}

void CartesianRotationTestCase::calcInitCond(AdvectionManager &advectionManager) {
    TimeLevelIndex<2> timeIdx;
    advectionManager.registerTracer("q0", "N/A", "background tracer");
    advectionManager.registerTracer("q1", "N/A", "cone tracer");
    advectionManager.registerTracer("q2", "N/A", "slotted cyliner tracer");
    density = &advectionManager.density();
    AdvectionTestCase::registerDefaultOutput();
    double q[3*mesh().totalNumGrid(CENTER, 2)];
    int l = 0;
    SpaceCoord x0(2); x0(0) = 0.5, x0(1) = 0.0;
    double R = 0.2;
    // Background tracer
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        q[l++] = 1.0;
    }
    // Cone tracer
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double d = domain().calcDistance(x0, x);
        if (d <= R) {
            q[l++] = R-d;
        } else {
            q[l++] = 0.0;
        }
    }
    // Slotted cylinder tracer
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double d = domain->calcDistance(x0, x);
        if (d <= R && fabs(x(0)-x0(0)) >= R/6.0) {
            q[l++] = 1.0;
        } else if (d <= R && fabs(x(0)-x0(0)) < R/6.0 && x(1)-x0(1) > 5.0/12.0*R) {
            q[l++] = 1.0;
        } else {
            q[l++] = 0.1;
        }
    }
    advectionManager.input(timeIdx, q);
}

} //lasm

#endif