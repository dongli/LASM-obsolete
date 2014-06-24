#ifndef __LASM_AdvectionTestCase_h__
#define __LASM_AdvectionTestCase_h__

#include "lasm_commons.h"
#include "AdvectionManager.h"

namespace lasm {

class AdvectionTestCase {
protected:
    Domain *domain;
    Mesh *mesh;
    const geomtk::TimeManager *timeManager;
    geomtk::IOManager<geomtk::RLLDataFile> io;
    int fileIdx;
    string subcase;
    VelocityField velocity;
    vector<ScalarField*> q;
public:
    AdvectionTestCase();
    virtual ~AdvectionTestCase();

    virtual void init(const ConfigManager &configManager,
                      const TimeManager &timeManager);

    virtual const Domain& getDomain() const { return *domain; }
    virtual const Mesh& getMesh() const { return *mesh; }
    virtual const VelocityField& getVelocityField() const { return velocity; }
    
    /**
     *  Output velocity field.
     *
     *  @param fileName   the output file name.
     *  @param oldTimeIdx the old time level index.
     */
    void outputVelocity(const string &fileName,
                        const TimeLevelIndex<2> &oldTimeIdx) const;

    /**
     *  Return the start time of the test case.
     *
     *  @return A Time object.
     */
    virtual Time getStartTime() const = 0;

    /**
     *  Return the end time of the test case.
     *
     *  @return A Time object.
     */
    virtual Time getEndTime() const = 0;

    /**
     *  Return the time step size of the test case.
     *
     *  @return The step size in seconds.
     */
    virtual double getStepSize() const = 0;

    /**
     *  Select the subcase to run.
     *
     *  @param subcaseName the name of the subcase.
     */
    virtual void selectSubcase(const string &subcaseName);

    /**
     *  Calculate initial condition and set tracers. This base function just
     *  transfer the settings of derived ones to tracers, so it will be called
     *  at the end of derived ones.
     */
    virtual void calcInitCond(AdvectionManager &advectionManager);

    /**
     *  Advance the test case one time step.
     *
     *  @param time    the time in seconds.
     *  @param timeIdx the time level index.
     */
    virtual void advance(double time, const TimeLevelIndex<2> &timeIdx) = 0;

    /**
     *  Calculate the solution of the test case if any, and reset the tracers
     *  for latter outputting.
     *
     *  @param dt               the time step size in seconds.
     *  @param timeIdx          the time level index.
     *  @param advectionManager the advection manager.
     */
    virtual void calcSolution(double dt, const TimeLevelIndex<2> &timeIdx,
                              AdvectionManager &advectionManager);
protected:
    /**
     *  Calculate the solution of the test case if any.
     *
     *  @param time    the time in seconds.
     *  @param timeIdx the time level index.
     *  @param q       the output solution.
     */
    virtual void calcSolution(double time, const TimeLevelIndex<2> &timeIdx,
                              ScalarField &q);
};

}

#endif // __LASM_AdvectionTestCase_h__
