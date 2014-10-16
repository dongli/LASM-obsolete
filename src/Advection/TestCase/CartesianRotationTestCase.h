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

    virtual Time startTime() const;

    virtual Time endTime() const;

    virtual double stepSize() const;

    virtual void calcInitCond(AdvectionManager &advectionManager);

    virtual void advance(double time, const TimeLevelIndex<2> &timeIdx);
};

}

#endif