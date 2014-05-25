#include "AdvectionTestCase.h"

namespace lasm {

AdvectionTestCase::AdvectionTestCase() {
}

AdvectionTestCase::~AdvectionTestCase() {
    for (int i = 0; i < q.size(); ++i) {
        delete q[i];
    }
}

void AdvectionTestCase::init(const geomtk::TimeManager &timeManager) {
    this->timeManager = &timeManager;
}

void AdvectionTestCase::outputVelocity(const string &fileName,
                                       const TimeLevelIndex<2> &oldTimeIdx) const {
}

void AdvectionTestCase::selectSubcase(const string &subcase) {
    this->subcase = subcase;
}

void AdvectionTestCase::calcInitCond(AdvectionManager &advectionManager) {
    assert(q.size() > 0);
    geomtk::StampString name("q", ""), brief("test tracer ", "");
    for (int i = 0; i < q.size(); ++i) {
        advectionManager.registerTracer(name.run("%d", i), "N/A",
                                        brief.run("%d", i));
    }
    TimeLevelIndex<2> initTimeIdx;
    advectionManager.input(initTimeIdx, q);
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
