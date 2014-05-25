#ifndef __LASM_SolidRotationTestCase__
#define __LASM_SolidRotationTestCase__

#include "AdvectionTestCase.h"

namespace lasm {

class SolidRotationTestCase : public AdvectionTestCase {
protected:
    double angleSpeed, U0, alpha;
    SpaceCoord *axisPole, *c0, *cr0;
    double R, H0;
public:
    SolidRotationTestCase();
    virtual ~SolidRotationTestCase();

    virtual void init(const geomtk::TimeManager &timeManager);

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    void calcInitCond(AdvectionManager &advectionManager);
    void calcSolution(double dt, const TimeLevelIndex<2> &timeIdx,
                      AdvectionManager &advectionManager);
    void advance(double time, const TimeLevelIndex<2> &timeIdx);
protected:
    void calcSolution(double dt, const TimeLevelIndex<2> &timeIdx,
                      ScalarField &q);
};

}

#endif // __LASM_SolidRotationTestCase__
