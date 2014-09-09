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

    virtual void init(const ConfigManager &configManager,
                      TimeManager &timeManager);

    virtual Time getStartTime() const;
    virtual Time getEndTime() const;
    virtual double getStepSize() const;

    virtual void calcInitCond(AdvectionManager &advectionManager);

    void calcSolution(double dt, const TimeLevelIndex<2> &timeIdx,
                      AdvectionManager &advectionManager);

    virtual void advance(double time, const TimeLevelIndex<2> &timeIdx);
protected:
    void calcSolution(double dt, const TimeLevelIndex<2> &timeIdx,
                      ScalarField &q);
};

}

#endif // __LASM_SolidRotationTestCase__
