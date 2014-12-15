#include "BarotropicTestCase.h"

#ifdef LASM_SPHERE_DOMAIN

namespace lasm {

BarotropicTestCase::BarotropicTestCase() {
    subcase = "case1";
    REPORT_ONLINE;
}

BarotropicTestCase::~BarotropicTestCase() {
    delete _domain;
    delete _mesh;
    REPORT_OFFLINE;
}

void BarotropicTestCase::init(const ConfigManager &configManager,
                              TimeManager &timeManager) {
    int numLon = 160, numLat = 81;
    if (configManager.hasKey("test_case", "num_lon")) {
        configManager.getValue("test_case", "num_lon", numLon);
    }
    if (configManager.hasKey("test_case", "num_lat")) {
        configManager.getValue("test_case", "num_lat", numLat);
    }
    if (configManager.hasKey("test_case", "subcase")) {
        configManager.getValue("test_case", "subcase", subcase);
    }
    model.init(timeManager, numLon, numLat);
    _domain = &model.domain();
    _mesh = &model.mesh();
    velocity.create(mesh(), false, HAS_HALF_LEVEL);
    AdvectionTestCase::init(configManager, timeManager);
    // Append output fields.
    io.file(outputFileIdx).registerField("double", FULL_DIMENSION,
        {&model.geopotentialDepth(), &model.surfaceGeopotential()});
    // Set time step size.
    if (mesh().numGrid(0, FULL) == 160 && mesh().numGrid(1, FULL) == 81) {
        _stepSize = 1*TimeUnit::MINUTES;
    } else if (mesh().numGrid(0, FULL) == 240 &&
               mesh().numGrid(1, FULL) == 121) {
        _stepSize = 20*TimeUnit::SECONDS;
    } else {
        REPORT_ERROR("Unspecified time step size!");
    }
}

Time BarotropicTestCase::startTime() const {
    Time time;
    return time;
}

Time BarotropicTestCase::endTime() const {
    Time time;
    return time+10*TimeUnit::DAYS;
}

void BarotropicTestCase::calcInitCond(AdvectionManager &advectionManager) {
    TimeLevelIndex<2> timeIdx;
    advectionManager.registerTracer("q0", "N/A", "background tracer");
    advectionManager.registerTracer("q1", "N/A", "geopotential depth", true);
    advectionManager.registerTracer("q2", "N/A", "step tracer");
    density = &advectionManager.density();
    AdvectionTestCase::registerDefaultOutput(advectionManager);
    // set initial condition for barotropic model
    if (restartFileName == "N/A") {
        if (subcase == "case1") {
            testCase.calcInitCond(model);
        } else if (subcase == "case2") {
            int fileIdx = io.registerInputFile(*_mesh, "ic.barotropic.case2.nc");
            io.file(fileIdx).registerField("double", FULL_DIMENSION,
                {&model.zonalWind(), &model.meridionalWind(), &model.geopotentialDepth()});
            io.open(fileIdx);
            io.input<double, 2>(fileIdx, timeIdx, {&model.zonalWind(),
                                &model.meridionalWind(), &model.geopotentialDepth()});
            io.close(fileIdx);
            io.removeFile(fileIdx);
            model.zonalWind().applyBndCond(timeIdx);
            model.meridionalWind().applyBndCond(timeIdx);
            model.geopotentialDepth().applyBndCond(timeIdx);
        } else {
            REPORT_ERROR("Unknown subcase \"" << subcase << "\"!");
        }
        // set initial condition for tracers
        double *q = new double[3*mesh().totalNumGrid(CENTER, 2)];
        int l = 0;
        // background tracer
        for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
            q[l++] = 1.0;
        }
        // geopotential tracer
        const ScalarField &gd = model.geopotentialDepth();
        for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
            q[l++] = gd(i);
        }
        // step tracer
        for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
            const SpaceCoord &x = mesh().gridCoord(CENTER, i);
            if (x(0) > 160*RAD && x(0) < 200*RAD &&
                x(1) > 10*RAD  && x(1) < 40*RAD) {
                q[l++] = 1;
            } else {
                q[l++] = 0.1;
            }
        }
        // propagate initial conditions to advection manager
        advectionManager.input(timeIdx, q);
        delete [] q;
    } else {
        model.input(restartFileName);
        advectionManager.restart(restartFileName);
    }
}
    
void BarotropicTestCase::output(const TimeLevelIndex<2> &timeIdx,
                                AdvectionManager &advectionManager) {
    AdvectionTestCase::startOutput(timeIdx);
    io.output<double, 2>(outputFileIdx, timeIdx, {&model.geopotentialDepth()});
    io.output<double>(outputFileIdx, {&model.surfaceGeopotential()});
    AdvectionTestCase::finishOutput(timeIdx, advectionManager);
}

void BarotropicTestCase::advanceDynamics(double time,
                                         const TimeLevelIndex<2> &timeIdx) {
    if (!timeIdx.isCurrentIndex()) {
        model.integrate(timeIdx-1, stepSize());
    }
    int lonGridType, latGridType;
    lonGridType = velocity(0).gridType(0);
    latGridType = velocity(0).gridType(1);
    for (int j = mesh().js(latGridType); j <= mesh().je(latGridType); ++j) {
        for (int i = mesh().is(lonGridType); i <= mesh().ie(lonGridType); ++i) {
            velocity(0)(timeIdx, i, j) = model.zonalWind()(timeIdx, i, j);
        }
    }
    lonGridType = velocity(1).gridType(0);
    latGridType = velocity(1).gridType(1);
    for (int j = mesh().js(latGridType); j <= mesh().je(latGridType); ++j) {
        for (int i = mesh().is(lonGridType); i <= mesh().ie(lonGridType); ++i) {
            velocity(1)(timeIdx, i, j) = model.meridionalWind()(timeIdx, i, j);
        }
    }
    if (timeIdx.isCurrentIndex()) {
        velocity.applyBndCond(timeIdx);
    } else {
        velocity.applyBndCond(timeIdx, UPDATE_HALF_LEVEL);
    }
}

} // lasm

#endif
