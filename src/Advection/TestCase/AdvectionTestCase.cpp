#include "AdvectionTestCase.h"

namespace lasm {

AdvectionTestCase::AdvectionTestCase() {
    domain = NULL;
    mesh = NULL;
    densities = NULL;
    restartFileName = "N/A";
}

AdvectionTestCase::~AdvectionTestCase() {
}

void AdvectionTestCase::init(const ConfigManager &configManager,
                             TimeManager &timeManager) {
    this->timeManager = &timeManager;
    io.init(timeManager);
    string pattern, unit;
    int freq;
    configManager.getValue("test_case", "output_pattern", pattern);
    configManager.getValue("test_case", "output_freq_unit", unit);
    configManager.getValue("test_case", "output_freq", freq);
    if (configManager.hasKey("test_case", "restart_file")) {
        configManager.getValue("test_case", "restart_file", restartFileName);
    }
    if (configManager.hasKey("test_case", "subcase")) {
        configManager.getValue("test_case", "subcase", subcase);
    }
    if (unit == "steps") {
        outputFileIdx = io.registerOutputFile(*mesh, pattern,
                                              TimeStepUnit::STEP, freq);
    } else if (unit == "minutes") {
        outputFileIdx = io.registerOutputFile(*mesh, pattern,
                                              TimeStepUnit::MINUTE, freq);
    } else if (unit == "hours") {
        outputFileIdx = io.registerOutputFile(*mesh, pattern,
                                              TimeStepUnit::HOUR, freq);
    } else if (unit == "days") {
        outputFileIdx = io.registerOutputFile(*mesh, pattern,
                                              TimeStepUnit::DAY, freq);
    } else {
        REPORT_ERROR("Under construction!");
    }
}

void AdvectionTestCase::registerDefaultOutput() {
    io.file(outputFileIdx).registerField("double", FULL_DIMENSION,
                                         {&velocity(0), &velocity(1),
                                          &velocity.getDivergence()});
    for (int s = 0; s < densities->size(); ++s) {
        io.file(outputFileIdx).registerField("double", FULL_DIMENSION,
                                             {(*densities)[s]});
    }
}

void AdvectionTestCase::startOutput(const TimeLevelIndex<2> &timeIdx) {
    io.create(outputFileIdx);
    io.output<double, 2>(outputFileIdx, timeIdx,
                         {&velocity(0), &velocity(1),
                          &velocity.getDivergence()});
    for (int s = 0; s < densities->size(); ++s) {
        io.output<double, 2>(outputFileIdx, timeIdx, {(*densities)[s]});
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
}

void AdvectionTestCase::calcSolution(double dt,
                                     const TimeLevelIndex<2> &timeIdx,
                                     ScalarField &q) {
    REPORT_ERROR("calcSolution is not available!");
}

}
