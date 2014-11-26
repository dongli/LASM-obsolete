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

    virtual Time startTime() const;

    virtual Time endTime() const;

    virtual void calcInitCond(AdvectionManager &advectionManager);

    virtual void advanceDynamics(double time, const TimeLevelIndex<2> &timeIdx);

    virtual void evalVelocity(double dt, const SpaceCoord &x,
                              bool isMoveOnPole, Velocity &v) const;

    virtual void evalDivergence(double dt, const SpaceCoord &x,
                                double &div) const;
};

}

#endif // __LASM_DeformationTestCase__
