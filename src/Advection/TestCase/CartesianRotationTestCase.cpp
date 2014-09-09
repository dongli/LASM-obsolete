#include "CartesianRotationTestCase.h"

#ifdef USE_CARTESIAN_DOMAIN

namespace lasm {

CartesianRotationTestCase::CartesianRotationTestCase() {
    angleSpeed = 1.0*RAD;
    REPORT_ONLINE;
}

CartesianRotationTestCase::~CartesianRotationTestCase() {
    REPORT_OFFLINE;
}

void CartesianRotationTestCase::init(const ConfigManager &configManager,
                                     TimeManager &timeManager) {
    domain = new geomtk::CartesianDomain(2);
    domain->setAxis(0, "x", "x axis", "m", -1, geomtk::OPEN, 1, geomtk::OPEN);
    domain->setAxis(1, "y", "y axis", "m", -1, geomtk::OPEN, 1, geomtk::OPEN);
    mesh = new geomtk::CartesianMesh(*domain);
    int nx = 100, ny = 100;
    if (configManager.hasKey("test_case", "num_x")) {
        configManager.getValue("test_case", "num_x", nx);
    }
    if (configManager.hasKey("test_case", "num_y")) {
        configManager.getValue("test_case", "num_y", ny);
    }
    mesh->init(nx, ny);
    velocity.create(*mesh, true, HAS_HALF_LEVEL);
    AdvectionTestCase::init(configManager, timeManager);
}

Time CartesianRotationTestCase::getStartTime() const {
    Time res;
    return res;
}

Time CartesianRotationTestCase::getEndTime() const {
    Time res = getStartTime()+360;
    return res;
}

double CartesianRotationTestCase::getStepSize() const {
    return 1.0;
}

void CartesianRotationTestCase::advance(double time,
                                        const TimeLevelIndex<2> &timeIdx) {
    for (int i = 0; i < mesh->getTotalNumGrid(velocity(0).getStaggerLocation(), 2); ++i) {
        const SpaceCoord &x = mesh->getGridCoord(velocity(0).getStaggerLocation(), i);
        double theta = atan2(x(1), x(0));
        velocity(0)(timeIdx, i) = -sin(theta)*norm(x())*angleSpeed;
    }
    for (int i = 0; i < mesh->getTotalNumGrid(velocity(1).getStaggerLocation(), 2); ++i) {
        const SpaceCoord &x = mesh->getGridCoord(velocity(1).getStaggerLocation(), i);
        double theta = atan2(x(1), x(0));
        velocity(1)(timeIdx, i) = cos(theta)*norm(x())*angleSpeed;
    }
}

void CartesianRotationTestCase::calcInitCond(AdvectionManager &advectionManager) {
    TimeLevelIndex<2> timeIdx;
    advectionManager.registerTracer("q0", "N/A", "background tracer");
    advectionManager.registerTracer("q1", "N/A", "cone tracer");
    advectionManager.registerTracer("q2", "N/A", "slotted cyliner tracer");
    densities = &advectionManager.getDensities();
    AdvectionTestCase::registerDefaultOutput();
    double q[3*mesh->getTotalNumGrid(CENTER, 2)];
    int l = 0;
    SpaceCoord x0(2); x0(0) = 0.5, x0(1) = 0.0;
    double R = 0.2;
    // Background tracer
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, 2); ++i) {
        q[l++] = 1.0;
    }
    // Cone tracer
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh->getGridCoord(CENTER, i);
        double d = domain->calcDistance(x0, x);
        if (d <= R) {
            q[l++] = R-d;
        } else {
            q[l++] = 0.0;
        }
    }
    // Slotted cylinder tracer
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh->getGridCoord(CENTER, i);
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