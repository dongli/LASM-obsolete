#ifndef __LASM_DeformationTestCase__
#define __LASM_DeformationTestCase__

#include "AdvectionTestCase.h"

namespace lasm {

class DeformationTestCase : public AdvectionTestCase {
protected:
    double period;
public:
    DeformationTestCase();
    ~DeformationTestCase();

    virtual void init(const ConfigManager &configManager,
                      TimeManager &timeManager);

    virtual Time getStartTime() const;
    virtual Time getEndTime() const;
    virtual double getStepSize() const;

    virtual void calcInitCond(AdvectionManager &advectionManager);

    virtual void advance(double time, const TimeLevelIndex<2> &timeIdx);

    virtual void evalVelocity(double dt, const SpaceCoord &x,
                              bool isMoveOnPole, Velocity &v) const;

    virtual void evalDivergence(double dt, const SpaceCoord &x,
                                double &div) const;
};

}

#endif // __LASM_DeformationTestCase__
