#ifndef __LASM_TerminatorChemistryTestCase__
#define __LASM_TerminatorChemistryTestCase__


#include "DeformationTestCase.h"

namespace lasm {

class TerminatorChemistryTestCase : public DeformationTestCase {
protected:
    const double k1CenterLon = 300*RAD;
    const double k1CenterLat = 20*RAD;
    const double cly0 = 4e-6;
public:
    TerminatorChemistryTestCase();
    ~TerminatorChemistryTestCase();

    virtual void init(const ConfigManager &configManager,
                      TimeManager &timeManager);

    virtual Time startTime() const;

    virtual Time endTime() const;

    virtual void calcInitCond(AdvectionManager &advectionManager);

    virtual void advancePhysics(double time, const TimeLevelIndex<2> &timeIdx,
                                AdvectionManager &advectionManager);
protected:
    /**
     *  Calculate solar photolysis rate and recombination rate.
     *
     *  @param lon longitude.
     *  @param lat latitude.
     *  @param k1  solar photolysis rate.
     *  @param k2  recombination rate.
     */
    void calcChemicalReactionRate(double lon, double lat, double &k1, double &k2);

    void calcTendency(double lon, double lat, double cl, double cl2,
                      double &dCl, double &dCl2);
};

} // lasm

#endif // __LASM_TerminatorChemistryTestCase__
