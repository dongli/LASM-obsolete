#ifndef __LASM_AdvectionTestCase_h__
#define __LASM_AdvectionTestCase_h__

#include "lasm_commons.h"
#include "AdvectionManager.h"

namespace lasm {

class AdvectionTestCase {
protected:
    Domain *domain;
    Mesh *mesh;
    TimeManager *timeManager;
    IOManager io;
    int outputFileIdx;
    string restartFileName;
    string subcase;
    VelocityField velocity;
    vector<ScalarField*> *densities;
public:
    AdvectionTestCase();
    virtual ~AdvectionTestCase();

    virtual void init(const ConfigManager &configManager,
                      TimeManager &timeManager);

    virtual const Domain& getDomain() const { return *domain; }
    virtual const Mesh& getMesh() const { return *mesh; }
    virtual const VelocityField& getVelocityField() const { return velocity; }

    /**
     *  Register default output fields including velocity field and meshed
     *  tracer density fields.
     */
    void registerDefaultOutput();

    /**
     *  Output data.
     *
     *  @param timeIdx          the time level index.
     *  @param advectionManager the advection manager.
     */
    virtual void output(const TimeLevelIndex<2> &timeIdx,
                        AdvectionManager &advectionManager);

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
     *  Calculate initial condition and set tracers.
     *
     *  @param advectionManager the advection manager.
     */
    virtual void calcInitCond(AdvectionManager &advectionManager) = 0;

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

    void startOutput(const TimeLevelIndex<2> &timeIdx);
    
    void finishOutput(const TimeLevelIndex<2> &timeIdx,
                      AdvectionManager &advectionManager);
};

}

#endif // __LASM_AdvectionTestCase_h__
