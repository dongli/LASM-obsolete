#include "lasm.h"

int main(int argc, const char *argv[])
{
    if (argc != 2) {
        REPORT_ERROR("Wrong usage! A configuration file path is needed.");
    }
    geomtk::ConfigManager configManager;
    lasm::AdvectionTestCase *testCase;
    lasm::AdvectionManager advectionManager;
    geomtk::TimeManager timeManager;
    geomtk::TimeLevelIndex<2> oldTimeIdx, newTimeIdx;
    // Parse configuration.
    configManager.parse(argv[1]);
    // Choose test case.
    std::string testCaseName, subcaseName = "";
    configManager.getValue("test_case", "case_name", testCaseName);
    if (testCaseName == "terminator_chemistry") {
        testCase = new lasm::TerminatorChemistryTestCase();
    } else {
        REPORT_ERROR("Unknown test case!");
    }
    // Initialization.
    testCase->init(configManager, timeManager);
    timeManager.init(testCase->startTime(), testCase->endTime(),
                     testCase->stepSize());
    advectionManager.init(testCase->domain(), testCase->mesh(), configManager);
    // Calculate initial conditions.
    testCase->calcInitCond(advectionManager);
    testCase->advanceDynamics(timeManager.seconds(), oldTimeIdx);
    testCase->output(oldTimeIdx, advectionManager);
    // Integration loop.
    while (!timeManager.isFinished()) {
        newTimeIdx = oldTimeIdx+1;
        double time = timeManager.seconds()+timeManager.stepSize();
        testCase->advancePhysics(time, oldTimeIdx, advectionManager);
        for (int i = 0; i < testCase->numSubcycledStep(); ++i) {
            testCase->advanceDynamics(time, newTimeIdx);
            advectionManager.advance(testCase->subcycledStepSize(), newTimeIdx,
                                     testCase->velocityField());
        }
        timeManager.advance();
        oldTimeIdx.shift();
        testCase->output(oldTimeIdx, advectionManager);
    }
}
