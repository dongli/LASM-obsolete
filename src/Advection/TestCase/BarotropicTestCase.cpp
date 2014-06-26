#include "BarotropicTestCase.h"

namespace lasm {

BarotropicTestCase::BarotropicTestCase() {
    subcase = "case1";
    REPORT_ONLINE;
}

BarotropicTestCase::~BarotropicTestCase() {
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
    model.init(timeManager, numLon, numLat);
    domain = &model.getDomain();
    mesh = &model.getMesh();
    velocity.create(model.getMesh(), false, HAS_HALF_LEVEL);
    AdvectionTestCase::init(configManager, timeManager);
    // append output fields
    io.file(outputFileIdx).registerOutputField<double, 2, FULL_DIMENSION>
        (1, &model.getGeopotentialDepth());
    io.file(outputFileIdx).registerOutputField<double, 1, FULL_DIMENSION>
        (1, &model.getSurfaceGeopotential());
}

Time BarotropicTestCase::getStartTime() const {
    Time time;
    return time;
}

Time BarotropicTestCase::getEndTime() const {
    Time time;
    return time+1*TimeUnit::DAYS;
}

double BarotropicTestCase::getStepSize() const {
    if (model.getMesh().getNumGrid(0, FULL) == 160 &&
        model.getMesh().getNumGrid(1, FULL) == 81) {
        return 1*TimeUnit::MINUTES;
    } else if (model.getMesh().getNumGrid(0, FULL) == 240 &&
               model.getMesh().getNumGrid(1, FULL) == 121) {
        return 20*TimeUnit::SECONDS;
    } else {
        REPORT_ERROR("Unspecified time step size!");
    }
}

void BarotropicTestCase::calcInitCond(AdvectionManager &advectionManager) {
    TimeLevelIndex<2> timeIdx;
    // set initial condition for barotropic model
    if (restartFileName == "N/A") {
        if (subcase == "case1") {
            testCase.calcInitCond(model);
        } else if (subcase == "case2") {
            int fileIdx = io.registerInputFile(model.getMesh(), "erai_input.nc");
            io.file(fileIdx).registerInputField<double, 2, FULL_DIMENSION>
                (3, &model.getZonalWind(), &model.getMeridionalWind(), &model.getGeopotentialDepth());
            io.open(fileIdx);
            io.input<double, 2>(fileIdx, timeIdx, 3, &model.getZonalWind(),
                                &model.getMeridionalWind(), &model.getGeopotentialDepth());
            io.close(fileIdx);
            io.removeFile(fileIdx);
            model.getZonalWind().applyBndCond(timeIdx);
            model.getMeridionalWind().applyBndCond(timeIdx);
            model.getGeopotentialDepth().applyBndCond(timeIdx);
        } else {
            REPORT_ERROR("Unknown subcase \"" << subcase << "\"!");
        }
    } else {
        model.input(restartFileName);
    }
    // set initial condition for tracers
    advectionManager.registerTracer("q0", "N/A", "background tracer");
    advectionManager.registerTracer("q1", "N/A", "geopotential depth");
    advectionManager.registerTracer("q2", "N/A", "step tracer");
    meshedDensities = &advectionManager.getMeshedDensities();
    AdvectionTestCase::registerDefaultOutput();
    double q[3*model.getMesh().getTotalNumGrid(CENTER)];
    int l = 0;
    // background tracer
    for (int i = 0; i < model.getMesh().getTotalNumGrid(CENTER); ++i) {
        q[l++] = 1.0;
    }
    // geopotential tracer
    const ScalarField &gd = model.getGeopotentialDepth();
    for (int i = 0; i < model.getMesh().getTotalNumGrid(CENTER); ++i) {
        q[l++] = gd(i);
    }
    // step tracer
    SpaceCoord x(2);
    for (int i = 0; i < model.getMesh().getTotalNumGrid(CENTER); ++i) {
        mesh->getGridCoord(i, CENTER, x);
        if (x(0) > 160*RAD && x(0) < 200*RAD &&
            x(1) > 10*RAD  && x(1) < 40*RAD) {
            q[l++] = 1;
        } else {
            q[l++] = 0.1;
        }
    }
    // propagate initial conditions to advection manager
    advectionManager.input(timeIdx, q);
}
    
void BarotropicTestCase::output(const TimeLevelIndex<2> &timeIdx,
                                AdvectionManager &advectionManager) {
    AdvectionTestCase::startOutput(timeIdx);
    io.output<double, 2>(outputFileIdx, timeIdx, 1, &model.getGeopotentialDepth());
    io.output<double, 1>(outputFileIdx, 1, &model.getSurfaceGeopotential());
    AdvectionTestCase::finishOutput(timeIdx, advectionManager);
}

void BarotropicTestCase::advance(double time,
                                 const TimeLevelIndex<2> &timeIdx) {
    if (timeIdx.isCurrentIndex()) {
        model.integrate(timeIdx, getStepSize());
    } else {
        model.integrate(timeIdx-1, getStepSize());
    }
    for (int j = 0; j < model.getMesh().getNumGrid(1, velocity(0).getGridType(1)); ++j) {
        for (int i = 0; i < model.getMesh().getNumGrid(0, velocity(0).getGridType(0)); ++i) {
            velocity(0)(timeIdx, i, j) = model.getZonalWind()(timeIdx, i, j);
        }
    }
    for (int j = 0; j < model.getMesh().getNumGrid(1, velocity(1).getGridType(1)); ++j) {
        for (int i = 0; i < model.getMesh().getNumGrid(0, velocity(1).getGridType(0)); ++i) {
            velocity(1)(timeIdx, i, j) = model.getMeridionalWind()(timeIdx, i, j);
        }
    }
    if (timeIdx.isCurrentIndex()) {
        velocity.applyBndCond(timeIdx);
    } else {
        velocity.applyBndCond(timeIdx, UPDATE_HALF_LEVEL);
    }
}
    
}
