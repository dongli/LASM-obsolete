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

    virtual void init(const geomtk::TimeManager &timeManager);

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    void calcInitCond(AdvectionManager &advectionManager);
    void advance(double time, const TimeLevelIndex<2> &timeIdx);
};

}

#endif // __LASM_DeformationTestCase__
