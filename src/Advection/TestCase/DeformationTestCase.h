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

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    virtual void calcInitCond(AdvectionManager &advectionManager);

    virtual void advance(double time, const TimeLevelIndex<2> &timeIdx);
};

}

#endif // __LASM_DeformationTestCase__
