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
                      TimeManager &timeManager);

    virtual Time startTime() const;

    virtual Time endTime() const;
    
    virtual void calcInitCond(AdvectionManager &advectionManager);
    
    virtual void output(const TimeLevelIndex<2> &timeIdx,
                        AdvectionManager &advectionManager);
    
    virtual void advanceDynamics(double time, const TimeLevelIndex<2> &timeIdx);
};

}

#endif // __LASM_BarotropicTestCase__
