#ifndef __Geomtk_CartesianRotationTestCase__
#define __Geomtk_CartesianRotationTestCase__

#include "lasm_commons.h"
#include "AdvectionTestCase.h"

namespace lasm {

class CartesianRotationTestCase : public AdvectionTestCase {
    double angleSpeed;
public:
    CartesianRotationTestCase();
    virtual ~CartesianRotationTestCase();

    virtual void init(const ConfigManager &configManager,
                      TimeManager &timeManager);

    virtual Time getStartTime() const;
    virtual Time getEndTime() const;
    virtual double getStepSize() const;

    virtual void calcInitCond(AdvectionManager &advectionManager);

    virtual void advance(double time, const TimeLevelIndex<2> &timeIdx);
};

}

#endif