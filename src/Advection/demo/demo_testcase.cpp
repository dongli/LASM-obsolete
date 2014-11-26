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
    geomtk::TimeLevelIndex<2> oldTimeIdx;
    // Parse configuration.
    configManager.parse(argv[1]);
    // Choose test case.
    bool isTrueSolution = false;
    std::string testCaseName, subcaseName = "";
    configManager.getValue("test_case", "case_name", testCaseName);
#if defined USE_CARTESIAN_DOMAIN
    if (testCaseName == "cartesian_rotation") {
        testCase = new lasm::CartesianRotationTestCase();
    }
#elif defined USE_SPHERE_DOMAIN
    if (testCaseName == "rotation") {
        testCase = new lasm::SolidRotationTestCase();
        if (configManager.hasKey("test_case", "is_true_solution")) {
            configManager.getValue("test_case", "is_true_solution", isTrueSolution);
        }
    } else if (testCaseName == "deform") {
        testCase = new lasm::DeformationTestCase();
    } else if (testCaseName == "barotropic") {
        testCase = new lasm::BarotropicTestCase();
    }
#endif
    else {
        REPORT_ERROR("Unknown test_case \"" << testCaseName << "\"!");
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
        geomtk::TimeLevelIndex<2> newTimeIdx = oldTimeIdx+1;
        double time = timeManager.seconds()+timeManager.stepSize();
        testCase->advanceDynamics(time, newTimeIdx);
        if (!testCase->isUseAnalyticalVelocity()) {
            advectionManager.advance(timeManager.stepSize(), newTimeIdx,
                                     testCase->velocityField());
        } else {
            advectionManager.advance(timeManager.stepSize(), newTimeIdx,
                                     *testCase);
        }
        if (isTrueSolution) {
            testCase->calcSolution(timeManager.stepSize(), newTimeIdx,
                                   advectionManager);
        }
        timeManager.advance();
        oldTimeIdx.shift();
        testCase->output(oldTimeIdx, advectionManager);
    }
    delete testCase;
    return 0;
}
