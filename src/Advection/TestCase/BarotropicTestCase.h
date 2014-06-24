#ifndef __LASM_BarotropicTestCase__
#define __LASM_BarotropicTestCase__

#include "AdvectionTestCase.h"
#include "barotropic_model.h"

namespace lasm {

class BarotropicTestCase : public AdvectionTestCase {
protected:
    barotropic_model::ToyTestCase testCase;
    barotropic_model::BarotropicModel_A_ImplicitMidpoint model;
public:
    BarotropicTestCase();
    virtual ~BarotropicTestCase();

    virtual void init(const ConfigManager &configManager,
                      const TimeManager &timeManager);

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    virtual const Domain& getDomain() const { return model.getDomain(); }
    virtual const Mesh& getMesh() const { return model.getMesh(); }
    
    void calcInitCond(AdvectionManager &advectionManager);
    void advance(double time, const TimeLevelIndex<2> &timeIdx);
};

}

#endif // __LASM_BarotropicTestCase__
