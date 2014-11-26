#include "AdvectionTestCase.h"

namespace lasm {

AdvectionTestCase::AdvectionTestCase() {
    _domain = NULL;
    _mesh = NULL;
    density = NULL;
    restartFileName = "N/A";
    useAnalyticalVelocity = false;
}

AdvectionTestCase::~AdvectionTestCase() {
}

void AdvectionTestCase::init(const ConfigManager &configManager,
                             TimeManager &timeManager) {
    this->timeManager = &timeManager;
    io.init(timeManager);
    string pattern, unit;
    double freq;
    configManager.getValue("test_case", "output_pattern", pattern);
    configManager.getValue("test_case", "output_freq_unit", unit);
    configManager.getValue("test_case", "output_freq", freq);
    if (configManager.hasKey("test_case", "restart_file")) {
        configManager.getValue("test_case", "restart_file", restartFileName);
    }
    if (configManager.hasKey("test_case", "subcase")) {
        configManager.getValue("test_case", "subcase", subcase);
    }
    if (configManager.hasKey("test_case", "use_analytical_velocity")) {
        configManager.getValue("test_case", "use_analytical_velocity", useAnalyticalVelocity);
    }
    if (unit == "steps") {
        outputFileIdx = io.registerOutputFile(*_mesh, pattern,
                                              TimeStepUnit::STEP, freq);
    } else if (unit == "minutes") {
        outputFileIdx = io.registerOutputFile(*_mesh, pattern,
                                              TimeStepUnit::MINUTE, freq);
    } else if (unit == "hours") {
        outputFileIdx = io.registerOutputFile(*_mesh, pattern,
                                              TimeStepUnit::HOUR, freq);
    } else if (unit == "days") {
        outputFileIdx = io.registerOutputFile(*_mesh, pattern,
                                              TimeStepUnit::DAY, freq);
    } else {
        REPORT_ERROR("Unknown output frequency unit \"" << unit << "\"!");
    }
    volumes.create("vol", "m2", "grid cell volume", mesh(), CENTER, 2);
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        volumes(i) = mesh().cellVolume(i);
    }
}

void AdvectionTestCase::registerDefaultOutput() {
    if (!useAnalyticalVelocity) {
        io.file(outputFileIdx).registerField("double", FULL_DIMENSION,
                                             {&velocity(0), &velocity(1),
                                              &velocity.divergence()});
    }
    io.file(outputFileIdx).registerField("double", FULL_DIMENSION, {&volumes});
    for (int s = 0; s < density->size(); ++s) {
        io.file(outputFileIdx).registerField("double", FULL_DIMENSION,
                                             {(*density)[s]});
    }
}

void AdvectionTestCase::startOutput(const TimeLevelIndex<2> &timeIdx) {
    io.create(outputFileIdx);
    if (!useAnalyticalVelocity) {
        io.output<double, 2>(outputFileIdx, timeIdx,
                             {&velocity(0), &velocity(1),
                              &velocity.divergence()});
    }
    io.output<double>(outputFileIdx, {&volumes});
    for (int s = 0; s < density->size(); ++s) {
        io.output<double, 2>(outputFileIdx, timeIdx, {(*density)[s]});
    }
}

void AdvectionTestCase::finishOutput(const TimeLevelIndex<2> &timeIdx,
                                     AdvectionManager &advectionManager) {
    if (io.isFileActive(outputFileIdx)) {
        advectionManager.output(timeIdx, io.file(outputFileIdx).fileID);
        REPORT_NOTICE("File \"" << io.file(outputFileIdx).fileName << "\" is outputted.");
    }
    io.close(outputFileIdx);
}

void AdvectionTestCase::output(const TimeLevelIndex<2> &timeIdx,
                               AdvectionManager &advectionManager) {
    startOutput(timeIdx);
    finishOutput(timeIdx, advectionManager);
}

void AdvectionTestCase::calcSolution(double dt,
                                     const TimeLevelIndex<2> &timeIdx,
                                     AdvectionManager &advectionManager) {
    REPORT_ERROR("calcSolution is not available!");
}

void AdvectionTestCase::calcSolution(double dt,
                                     const TimeLevelIndex<2> &timeIdx,
                                     ScalarField &q) {
    REPORT_ERROR("calcSolution is not available!");
}

void AdvectionTestCase::evalVelocity(double dt, const SpaceCoord &x,
                                     bool isMoveOnPole, Velocity &v) const {
    REPORT_ERROR("calcSolution is not available!");
}

void AdvectionTestCase::evalDivergence(double dt, const SpaceCoord &x,
                                       double &div) const {
    REPORT_ERROR("calcSolution is not available!");
}

void AdvectionTestCase::advanceDynamics(double time,
                                        const TimeLevelIndex<2> &timeIdx) {
    REPORT_ERROR("advance is not available!");
}

void AdvectionTestCase::advancePhysics(double time,
                                       const TimeLevelIndex<2> &timeIdx,
                                       lasm::AdvectionManager &advectionManager) {
    REPORT_ERROR("advance is not available!");
}

}
