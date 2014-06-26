#include "AdvectionTestCase.h"

namespace lasm {

AdvectionTestCase::AdvectionTestCase() {
    domain = NULL;
    mesh = NULL;
    meshedDensities = NULL;
    restartFileName = "N/A";
}

AdvectionTestCase::~AdvectionTestCase() {
}

void AdvectionTestCase::init(const ConfigManager &configManager,
                             TimeManager &timeManager) {
    this->timeManager = &timeManager;
    io.init(timeManager);
    string prefix, unit;
    int freq;
    configManager.getValue("test_case", "output_prefix", prefix);
    configManager.getValue("test_case", "output_freq_unit", unit);
    configManager.getValue("test_case", "output_freq", freq);
    if (configManager.hasKey("test_case", "restart_file")) {
        configManager.getValue("test_case", "restart_file", restartFileName);
    }
    if (configManager.hasKey("test_case", "subcase")) {
        configManager.getValue("test_case", "subcase", subcase);
    }
    if (unit == "steps") {
        outputFileIdx = io.registerOutputFile(*mesh, prefix,
                                              IOFrequencyUnit::STEPS, freq);
    } else if (unit == "minutes") {
        outputFileIdx = io.registerOutputFile(*mesh, prefix,
                                              IOFrequencyUnit::MINUTES, freq);
    } else if (unit == "hours") {
        outputFileIdx = io.registerOutputFile(*mesh, prefix,
                                              IOFrequencyUnit::HOURS, freq);
    } else if (unit == "days") {
        outputFileIdx = io.registerOutputFile(*mesh, prefix,
                                              IOFrequencyUnit::DAYS, freq);
    } else {
        REPORT_ERROR("Under construction!");
    }
}

void AdvectionTestCase::registerDefaultOutput() {
    io.file(outputFileIdx).registerOutputField<double, 2, FULL_DIMENSION>
        (4, &velocity(0), &velocity(1), &velocity.getDivergence(),
         &velocity.getVorticity()[0]);
    for (int s = 0; s < meshedDensities->size(); ++s) {
        io.file(outputFileIdx).registerOutputField<double, 2, FULL_DIMENSION>
            (1, (*meshedDensities)[s]);
    }
}

void AdvectionTestCase::startOutput(const TimeLevelIndex<2> &timeIdx) {
    io.create(outputFileIdx);
    io.output<double, 2>(outputFileIdx, timeIdx, 4, &velocity(0), &velocity(1),
                         &velocity.getDivergence(), &velocity.getVorticity()[0]);
    for (int s = 0; s < meshedDensities->size(); ++s) {
        io.output<double, 2>(outputFileIdx, timeIdx, 1, (*meshedDensities)[s]);
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
