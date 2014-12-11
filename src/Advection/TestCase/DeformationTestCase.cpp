#include "DeformationTestCase.h"

#ifdef LASM_SPHERE_DOMAIN

namespace lasm {

DeformationTestCase::DeformationTestCase() {
    subcase = "case4";
    period = 5;
    _stepSize = period/600;
    REPORT_ONLINE;
}

DeformationTestCase::~DeformationTestCase() {
    delete _mesh;
    delete _domain;
    REPORT_OFFLINE;
}

void DeformationTestCase::init(const ConfigManager &configManager,
                               TimeManager &timeManager) {
    // Initialize domain.
    _domain = new geomtk::SphereDomain(2);
    _domain->radius() = 1;
    // Initialize mesh.
    _mesh = new geomtk::RLLMesh(*_domain);
    int numLon = 240, numLat = 121;
    if (configManager.hasKey("test_case", "num_lon")) {
        configManager.getValue("test_case", "num_lon", numLon);
    }
    if (configManager.hasKey("test_case", "num_lat")) {
        configManager.getValue("test_case", "num_lat", numLat);
    }
    _mesh->init(numLon, numLat);
    // Call super class initialization.
    AdvectionTestCase::init(configManager, timeManager);
    // Initialize velocity.
    if (!useAnalyticalVelocity) {
        velocity.create(mesh(), true, HAS_HALF_LEVEL);
    }
}

Time DeformationTestCase::startTime() const {
    Time time;
    return time;
}

Time DeformationTestCase::endTime() const {
    Time time = startTime()+period;
    return time;
}

void DeformationTestCase::advanceDynamics(double time,
                                          const TimeLevelIndex<2> &timeIdx) {
    if (useAnalyticalVelocity) return;
    double cosT = cos(M_PI*time/period);
    double k, R = domain().radius();
    // advance velocity
    if (subcase == "case1") {
        k = 2.4;
        for (int i = 0; i < mesh().totalNumGrid(velocity(0).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh().gridCoord(velocity(0).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(0)(timeIdx, i) = k*pow(sin(lon*0.5), 2.0)*sin(lat*2.0)*cosT;
        }
        for (int i = 0; i < mesh().totalNumGrid(velocity(1).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh().gridCoord(velocity(1).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(1)(timeIdx, i) = k*0.5*sin(lon)*cos(lat)*cosT;
        }
    } else if (subcase == "case2") {
        k = 2.0;
        for (int i = 0; i < mesh().totalNumGrid(velocity(0).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh().gridCoord(velocity(0).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(0)(timeIdx, i) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT;
        }
        for (int i = 0; i < mesh().totalNumGrid(velocity(1).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh().gridCoord(velocity(1).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(1)(timeIdx, i) = k*sin(lon*2.0)*cos(lat)*cosT;
        }
    } else if (subcase == "case3") {
        k = 5.0*R/period;
        for (int i = 0; i < mesh().totalNumGrid(velocity(0).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh().gridCoord(velocity(0).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(0)(timeIdx, i) = -k*pow(sin(lon), 2.0)*sin(lat*2.0)*pow(cos(lat), 2.0)*cosT;
        }
        for (int i = 0; i < mesh().totalNumGrid(velocity(1).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh().gridCoord(velocity(1).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocity(1)(timeIdx, i) = k*0.5*sin(lon)*pow(cos(lat), 3.0)*cosT;
        }
    } else if (subcase == "case4") {
        k = 10.0*R/period;
        double c1 = PI2*time/period;
        double c2 = PI2*R/period;
        for (int i = 0; i < mesh().totalNumGrid(velocity(0).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh().gridCoord(velocity(0).staggerLocation(), i);
            double lon = x(0)-c1, lat = x(1);
            velocity(0)(timeIdx, i) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT+c2*cos(lat);
        }
        for (int i = 0; i < mesh().totalNumGrid(velocity(1).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = mesh().gridCoord(velocity(1).staggerLocation(), i);
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
    // Register tracer species.
    advectionManager.registerTracer("q0", "N/A", "background tracer");
    advectionManager.registerTracer("q1", "N/A", "cosine hills tracer");
    advectionManager.registerTracer("q2", "N/A", "q1 correlated tracer");
    advectionManager.registerTracer("q3", "N/A", "slotted cylinders tracer");
    advectionManager.registerTracer("q4", "N/A", "Gaussian hills tracer");
    advectionManager.registerTracer("q5", "N/A", "slotted cylinders tracer");
    advectionManager.registerTracer("q6", "N/A", "displaced slotted cylinders tracer");
    advectionManager.registerTracer("q7", "N/A", "slotted cylinders residual tracer");
    density = &advectionManager.density();
    AdvectionTestCase::registerDefaultOutput();
    // Set initial conditions for each tracer species.
    SpaceCoord c0(2), c1(2);
    c0.setCoord(M_PI*5.0/6.0, 0.0); c0.transformToCart(domain());
    c1.setCoord(M_PI*7.0/6.0, 0.0); c1.transformToCart(domain());
    double hmax, r, g, a, b, c;
    double *q = new double[8*mesh().totalNumGrid(CENTER, 2)];
    int l = 0;
    // - background tracer
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        q[l++] = 1.0;
    }
    // - cosine hills tracer
    hmax = 1, r = domain().radius()*0.5, g = 0.1, c = 0.9;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double r0 = domain().calcDistance(x, c0);
        double r1 = domain().calcDistance(x, c1);
        if (r0 < r) {
            q[l++] = g+c*hmax*0.5*(1+cos(M_PI*r0/r));
        } else if (r1 < r) {
            q[l++] = g+c*hmax*0.5*(1+cos(M_PI*r1/r));
        } else {
            q[l++] = g;
        }
    }
    // - tracer correlated to cosine hills tracer
    a = -0.8, b = 0.9;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        q[l] = a*pow(q[l-mesh().totalNumGrid(CENTER, 2)], 2)+b; l++;
    }
    // - slotted cylinders tracer
    b = 0.1, c = 1.0, r = 0.5;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double r0 = domain().calcDistance(x, c0);
        double r1 = domain().calcDistance(x, c1);
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
    // - Gaussian hills tracer
    hmax = 0.95, b = 5.0;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        vec d0 = x.cartCoord()-c0.cartCoord();
        vec d1 = x.cartCoord()-c1.cartCoord();
        q[l++] = hmax*(exp(-b*dot(d0, d0))+exp(-b*dot(d1, d1)));
    }
    // - Another slotted cylinders tracer
    c0.setCoord(M_PI*3.0/4.0, 0.0); c0.transformToCart(domain());
    c1.setCoord(M_PI*5.0/4.0, 0.0); c1.transformToCart(domain());
    b = 0.1, c = 1.0/3.0, r = 0.5;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double r0 = domain().calcDistance(x, c0);
        double r1 = domain().calcDistance(x, c1);
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
    // - Displaced slotted cylinders tracer
    c0.setCoord(M_PI*3.0/4.0, M_PI/18.0); c0.transformToCart(domain());
    c1.setCoord(M_PI*5.0/4.0, -M_PI/18.0); c1.transformToCart(domain());
    b = 0.1, c = 2.0/3.0, r = 0.5;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double r0 = domain().calcDistance(x, c0);
        double r1 = domain().calcDistance(x, c1);
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
    // - Slotted cylinders residual tracer
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        q[l] = 1.0-q[l-2*mesh().totalNumGrid(CENTER, 2)]-q[l-mesh().totalNumGrid(CENTER, 2)];
        l++;
    }
    // Propagate initial conditions to advection manager.
    TimeLevelIndex<2> timeIdx;
    advectionManager.input(timeIdx, q);
    delete [] q;
}

void DeformationTestCase::evalVelocity(double dt, const SpaceCoord &x,
                                       bool isMoveOnPole, Velocity &v) const {
    double time = timeManager->seconds()+dt;
    double cosT = cos(M_PI*time/period);
    double k, R = domain().radius();
    if (subcase == "case4") {
        k = 10.0*R/period;
        double c1 = PI2*time/period;
        double c2 = PI2*R/period;
        double lon = x(0)-c1, lat = x(1);
        v(0) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT+c2*cos(lat);
        v(1) = k*sin(lon*2.0)*cos(lat)*cosT;
        if (isMoveOnPole) {
            v.transformToPS(x);
        }
    } else {
        REPORT_ERROR("Under construction!");
    }
}

void DeformationTestCase::evalDivergence(double dt, const SpaceCoord &x,
                                         double &div) const {
    double time = timeManager->seconds()+dt;
    if (subcase == "case4") {
        div = 0.0;
    } else {
        REPORT_ERROR("Under construction!");
    }
}

} // lasm

#endif
