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
                              const TimeManager &timeManager) {
    AdvectionTestCase::init(configManager, timeManager);
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
    model.init(numLon, numLat);
    velocity.create(model.getMesh(), false, HAS_HALF_LEVEL);

    io.init(timeManager);
    fileIdx = io.registerOutputFile(model.getMesh(), "barotropic-testcase",
                                    geomtk::IOFrequencyUnit::STEPS, 10);
    io.file(fileIdx).registerOutputField<double, 2, FULL_DIMENSION>
        (5, &model.getGeopotentialDepth(), &velocity(0), &velocity(1),
         &velocity.getDivergence(), &velocity.getVorticity()[0]);
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
    SpaceCoord x(2);
    TimeLevelIndex<2> timeIdx;
    // -------------------------------------------------------------------------
    // set initial condition for barotropic model
    if (subcase == "case1") {
        testCase.calcInitCond(model);
    } else if (subcase == "case2") {
        int inputFileIdx = io.registerInputFile(model.getMesh(), "erai_input.nc");
        io.file(inputFileIdx).registerInputField<double, 2, FULL_DIMENSION>
            (3, &model.getZonalWind(), &model.getMeridionalWind(), &model.getGeopotentialDepth());
        io.open(inputFileIdx);
        io.input<double, 2>(inputFileIdx, timeIdx, 3, &model.getZonalWind(),
                            &model.getMeridionalWind(), &model.getGeopotentialDepth());
        io.close(inputFileIdx);
        model.getZonalWind().applyBndCond(timeIdx);
        model.getMeridionalWind().applyBndCond(timeIdx);
        model.getGeopotentialDepth().applyBndCond(timeIdx);
    }
    const ScalarField &gd = model.getGeopotentialDepth();
    // -------------------------------------------------------------------------
    // set initial condition for tracers
    ScalarField *q0; // reference tracer
    ScalarField *q1; // continuous tracer
    ScalarField *q2; // discontinuous tracer
    // reference tracer
    q.push_back(new ScalarField); q0 = q.back();
    q0->create("", "", "", model.getMesh(), CENTER);
    for (int i = 0; i < model.getMesh().getTotalNumGrid(CENTER); ++i) {
        (*q0)(timeIdx, i) = 1.0;
    }
    // continuous tracer
    q.push_back(new ScalarField); q1 = q.back();
    q1->create("", "", "", model.getMesh(), CENTER);
    for (int i = 0; i < model.getMesh().getTotalNumGrid(CENTER); ++i) {
        (*q1)(timeIdx, i) = gd(i);
    }
    // discontinous tracer
    q.push_back(new ScalarField); q2 = q.back();
    q2->create("", "", "", model.getMesh(), CENTER);
    for (int i = 0; i < model.getMesh().getTotalNumGrid(CENTER); ++i) {
        model.getMesh().getGridCoord(i, CENTER, x);
        if (x(0) > 160*RAD && x(0) < 200*RAD &&
            x(1) > 10*RAD  && x(1) < 40*RAD) {
            (*q2)(timeIdx, i) = 1;
        } else {
            (*q2)(timeIdx, i) = 0.1;
        }
    }
    // -------------------------------------------------------------------------
    AdvectionTestCase::calcInitCond(advectionManager);
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
    io.create(fileIdx);
    io.output<double, 2>(fileIdx, timeIdx, 5, &model.getGeopotentialDepth(),
                         &velocity(0), &velocity(1), &velocity.getDivergence(),
                         &velocity.getVorticity()[0]);
    io.close(fileIdx);
}
    
}
