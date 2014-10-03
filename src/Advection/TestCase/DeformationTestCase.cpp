#include "DeformationTestCase.h"

#ifdef USE_SPHERE_DOMAIN

namespace lasm {

DeformationTestCase::DeformationTestCase() {
    subcase = "case4";
    period = 5;
    REPORT_ONLINE;
}

DeformationTestCase::~DeformationTestCase() {
    delete mesh;
    delete domain;
    REPORT_OFFLINE;
}

void DeformationTestCase::init(const ConfigManager &configManager,
                               TimeManager &timeManager) {
    // initialize domain
    domain = new geomtk::SphereDomain(2);
    domain->setRadius(1);
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

    AdvectionTestCase::init(configManager, timeManager);
}

Time DeformationTestCase::getStartTime() const {
    Time time;
    return time;
}

Time DeformationTestCase::getEndTime() const {
    Time time = getStartTime()+period;
    return time;
}

double DeformationTestCase::getStepSize() const {
    return period/600.0;
}

void DeformationTestCase::advance(double time,
                                  const TimeLevelIndex<2> &timeIdx) {
    double cosT = cos(M_PI*time/period);
    double k, R = domain->getRadius();
    // advance velocity
    if (subcase == "case1") {
        k = 2.4;
        for (int i = 0; i < mesh->getTotalNumGrid(velocity(0).getStaggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh->getGridCoord(velocity(0).getStaggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(0)(timeIdx, i) = k*pow(sin(lon*0.5), 2.0)*sin(lat*2.0)*cosT;
        }
        for (int i = 0; i < mesh->getTotalNumGrid(velocity(1).getStaggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh->getGridCoord(velocity(1).getStaggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(1)(timeIdx, i) = k*0.5*sin(lon)*cos(lat)*cosT;
        }
    } else if (subcase == "case2") {
        k = 2.0;
        for (int i = 0; i < mesh->getTotalNumGrid(velocity(0).getStaggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh->getGridCoord(velocity(0).getStaggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(0)(timeIdx, i) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT;
        }
        for (int i = 0; i < mesh->getTotalNumGrid(velocity(1).getStaggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh->getGridCoord(velocity(1).getStaggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(1)(timeIdx, i) = k*sin(lon*2.0)*cos(lat)*cosT;
        }
    } else if (subcase == "case3") {
        k = 5.0*R/period;
        for (int i = 0; i < mesh->getTotalNumGrid(velocity(0).getStaggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh->getGridCoord(velocity(0).getStaggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(0)(timeIdx, i) = -k*pow(sin(lon), 2.0)*sin(lat*2.0)*pow(cos(lat), 2.0)*cosT;
        }
        for (int i = 0; i < mesh->getTotalNumGrid(velocity(1).getStaggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh->getGridCoord(velocity(1).getStaggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(1)(timeIdx, i) = k*0.5*sin(lon)*pow(cos(lat), 3.0)*cosT;
        }
    } else if (subcase == "case4") {
        k = 10.0*R/period;
        double c1 = PI2*time/period;
        double c2 = PI2*R/period;
        for (int i = 0; i < mesh->getTotalNumGrid(velocity(0).getStaggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh->getGridCoord(velocity(0).getStaggerLocation(), i);
            double lon = x(0)-c1, lat = x(1);
            velocity(0)(timeIdx, i) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT+c2*cos(lat);
        }
        for (int i = 0; i < mesh->getTotalNumGrid(velocity(1).getStaggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh->getGridCoord(velocity(1).getStaggerLocation(), i);
            double lon = x(0)-c1, lat = x(1);
            velocity(1)(timeIdx, i) = k*sin(lon*2.0)*cos(lat)*cosT;
        }
    }
    if (timeIdx.isCurrentIndex()) {
        velocity.applyBndCond(timeIdx);
    } else {
        velocity.applyBndCond(timeIdx, UPDATE_HALF_LEVEL);
    }
}

void DeformationTestCase::calcInitCond(AdvectionManager &advectionManager) {
    advectionManager.registerTracer("q0", "N/A", "background tracer");
    advectionManager.registerTracer("q1", "N/A", "cosine hills tracer");
    advectionManager.registerTracer("q2", "N/A", "q1 correlated tracer");
    advectionManager.registerTracer("q3", "N/A", "slotted cylinders tracer");
    advectionManager.registerTracer("q4", "N/A", "Gaussian hills tracer");
    densities = &advectionManager.getDensities();
    AdvectionTestCase::registerDefaultOutput();
    SpaceCoord c0(2), c1(2);
    c0.setCoord(M_PI*5.0/6.0, 0.0); c0.transformToCart(*domain);
    c1.setCoord(M_PI*7.0/6.0, 0.0); c1.transformToCart(*domain);
    double hmax, r, g, a, b, c;
    double *q = new double[5*mesh->getTotalNumGrid(CENTER, 2)];
    int l = 0;
    // background tracer
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, 2); ++i) {
        q[l++] = 1.0;
    }
    // cosine hills tracer
    hmax = 1, r = domain->getRadius()*0.5, g = 0.1, c = 0.9;
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh->getGridCoord(CENTER, i);
        double r0 = domain->calcDistance(x, c0);
        double r1 = domain->calcDistance(x, c1);
        if (r0 < r) {
            q[l++] = g+c*hmax*0.5*(1+cos(M_PI*r0/r));
        } else if (r1 < r) {
            q[l++] = g+c*hmax*0.5*(1+cos(M_PI*r1/r));
        } else {
            q[l++] = g;
        }
    }
    // tracer correlated to cosine hills tracer
    a = -0.8, b = 0.9;
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, 2); ++i) {
        q[l] = a*pow(q[l-mesh->getTotalNumGrid(CENTER, 2)], 2)+b; l++;
    }
    // slotted cylinders tracer
    b = 0.1, c = 1.0, r = 0.5;
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh->getGridCoord(CENTER, i);
        double r0 = domain->calcDistance(x, c0);
        double r1 = domain->calcDistance(x, c1);
        if ((r0 <= r && fabs(x(0)-c0(0)) >= r/6.0) ||
            (r1 <= r && fabs(x(0)-c1(0)) >= r/6.0))
            q[l++] = c;
        else if (r0 <= r && fabs(x(0)-c0(0)) < r/6.0 &&
                 x(1)-c0(1) < -5.0/12.0*r)
            q[l++] = c;
        else if (r1 <= r && fabs(x(0)-c1(0)) < r/6.0 &&
                 x(1)-c1(1) > 5.0/12.0*r)
            q[l++] = c;
        else
            q[l++] = b;
    }
    // Gaussian hills tracer
    hmax = 0.95, b = 5.0;
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh->getGridCoord(CENTER, i);
        vec d0 = x.getCartCoord()-c0.getCartCoord();
        vec d1 = x.getCartCoord()-c1.getCartCoord();
        q[l++] = hmax*(exp(-b*dot(d0, d0))+exp(-b*dot(d1, d1)));
    }
    // propagate initial conditions to advection manager
    TimeLevelIndex<2> timeIdx;
    advectionManager.input(timeIdx, q);
}

} // lasm

#endif
