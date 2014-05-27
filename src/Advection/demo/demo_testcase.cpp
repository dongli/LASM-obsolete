#include "lasm.h"

int main(int argc, const char *argv[])
{
    // -------------------------------------------------------------------------
    // check arguments
    if (argc != 2) {
        REPORT_ERROR("Wrong usage! A configuration file path is needed.");
    }
    // -------------------------------------------------------------------------
    geomtk::ConfigManager configManager;
    lasm::AdvectionTestCase *testCase;
    lasm::AdvectionManager advectionManager;
    geomtk::TimeManager timeManager;
    geomtk::TimeLevelIndex<2> oldTimeIdx;
    // -------------------------------------------------------------------------
    // parse configuration
    configManager.parse(argv[1]);
    // -------------------------------------------------------------------------
    // choose test case
    bool isTrueSolution = false;
    std::string testCaseName, subcaseName = "";
    configManager.getValue("lasm", "test_case", testCaseName);
    if (testCaseName == "rotation") {
        testCase = new lasm::SolidRotationTestCase();
        if (configManager.hasKey("lasm", "is_true_solution")) {
            configManager.getValue("lasm", "is_true_solution", isTrueSolution);
        }
    } else if (testCaseName == "deform") {
        testCase = new lasm::DeformationTestCase();
        if (configManager.hasKey("lasm", "subcase")) {
            configManager.getValue("lasm", "subcase", subcaseName);
            testCase->selectSubcase(subcaseName);
        }
    } else if (testCaseName == "barotropic") {
        testCase = new lasm::BarotropicTestCase();
    } else {
        REPORT_ERROR("Unknown test_case \"" << testCaseName << "\"!");
    }
    // -------------------------------------------------------------------------
    // initialization
    timeManager.init(testCase->getStartTime(), testCase->getEndTime(),
                     testCase->getStepSize());
    testCase->init(configManager, timeManager);
    advectionManager.init(testCase->getDomain(), testCase->getMesh(),
                          configManager, timeManager);
    
    testCase->calcInitCond(advectionManager);
    advectionManager.output(oldTimeIdx);
    testCase->advance(timeManager.getSeconds(), oldTimeIdx);
    // -------------------------------------------------------------------------
    // integration loop
    while (!timeManager.isFinished()) {
        geomtk::TimeLevelIndex<2> newTimeIdx = oldTimeIdx+1;
        double time = timeManager.getSeconds()+timeManager.getStepSize();
        testCase->advance(time, newTimeIdx);
        advectionManager.advance(timeManager.getStepSize(), newTimeIdx,
                                 testCase->getVelocityField());
        if (isTrueSolution) {
            testCase->calcSolution(timeManager.getStepSize(),
                                   newTimeIdx, advectionManager);
        }
        timeManager.advance();
        oldTimeIdx.shift();
        advectionManager.output(oldTimeIdx);
    }
    delete testCase;
    return 0;
}
