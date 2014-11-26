#ifndef __LASM_AdvectionTestCase_h__
#define __LASM_AdvectionTestCase_h__

#include "lasm_commons.h"
#include "AdvectionManager.h"

namespace lasm {

class AdvectionTestCase {
protected:
    Domain *_domain;
    Mesh *_mesh;
    TimeManager *timeManager;
    IOManager io;
    int outputFileIdx;
    string restartFileName;
    string subcase;
    VelocityField velocity;
    SingleScalarField volumes;
    vector<ScalarField*> *density;
    bool useAnalyticalVelocity;
    double _stepSize;
    double _subcycledStepSize;
    int _numSubcycledStep;
public:
    AdvectionTestCase();
    virtual ~AdvectionTestCase();

    virtual void init(const ConfigManager &configManager,
                      TimeManager &timeManager);

    virtual const Domain& domain() const { return *_domain; }

    virtual const Mesh& mesh() const { return *_mesh; }

    virtual const VelocityField& velocityField() const { return velocity; }

    bool isUseAnalyticalVelocity() const { return useAnalyticalVelocity; }

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
    virtual Time startTime() const = 0;

    /**
     *  Return the end time of the test case.
     *
     *  @return A Time object.
     */
    virtual Time endTime() const = 0;

    /**
     *  Return the time step size of the test case.
     *
     *  @return The step size in seconds.
     */
    double stepSize() const { return _stepSize; }

    double subcycledStepSize() const { return _subcycledStepSize; }

    int numSubcycledStep() const { return _numSubcycledStep; }

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
    virtual void advanceDynamics(double time, const TimeLevelIndex<2> &timeIdx);

    virtual void advancePhysics(double time, const TimeLevelIndex<2> &timeIdx,
                                AdvectionManager &advectionManager);

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

    virtual void evalVelocity(double dt, const SpaceCoord &x,
                              bool isMoveOnPole, Velocity &v) const;

    virtual void evalDivergence(double dt, const SpaceCoord &x,
                                double &div) const;
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
