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
    timeManager.init(testCase->getStartTime(), testCase->getEndTime(),
                     testCase->getStepSize());
    advectionManager.init(testCase->getDomain(), testCase->getMesh(), configManager);

    testCase->calcInitCond(advectionManager);
    testCase->advance(timeManager.getSeconds(), oldTimeIdx);
    testCase->output(oldTimeIdx, advectionManager);
    // Integration loop.
    while (!timeManager.isFinished()) {
        geomtk::TimeLevelIndex<2> newTimeIdx = oldTimeIdx+1;
        double time = timeManager.getSeconds()+timeManager.getStepSize();
        testCase->advance(time, newTimeIdx);
        if (!testCase->isUseAnalyticalVelocity()) {
            advectionManager.advance(timeManager.getStepSize(), newTimeIdx,
                                     testCase->getVelocityField());
        } else {
            advectionManager.advance(timeManager.getStepSize(), newTimeIdx,
                                     *testCase);
        }
        if (isTrueSolution) {
            testCase->calcSolution(timeManager.getStepSize(),
                                   newTimeIdx, advectionManager);
        }
        timeManager.advance();
        oldTimeIdx.shift();
        testCase->output(oldTimeIdx, advectionManager);
    }
    delete testCase;
    return 0;
}
