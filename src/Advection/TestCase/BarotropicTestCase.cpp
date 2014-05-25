#include "BarotropicTestCase.h"

namespace lasm {

BarotropicTestCase::BarotropicTestCase() {
    REPORT_ONLINE;
}

BarotropicTestCase::~BarotropicTestCase() {
    REPORT_OFFLINE;
}

void BarotropicTestCase::init(const geomtk::TimeManager &timeManager) {
    AdvectionTestCase::init(timeManager);
    model.init(160, 81);
    velocity.create(model.getMesh(), false, HAS_HALF_LEVEL);

    io.init(timeManager);
    fileIdx = io.registerOutputFile(model.getMesh(), "barotropic-output",
                                    geomtk::IOFrequencyUnit::STEPS, 1);
    io.file(fileIdx).registerOutputField<double, 1, barotropic_model::FULL_DIMENSION>(1, &model.getGeopotentialHeight());
}

Time BarotropicTestCase::getStartTime() const {
    Time time;
    return time;
}

Time BarotropicTestCase::getEndTime() const {
    Time time;
    return time+2*TimeUnit::DAYS;
}

double BarotropicTestCase::getStepSize() const {
    return 1*TimeUnit::MINUTES;
}

void BarotropicTestCase::calcInitCond(AdvectionManager &advectionManager) {
    SpaceCoord x(2);
    // -------------------------------------------------------------------------
    // set initial condition for barotropic model
    x.setCoord(120*RAD, 35*RAD);
    testCase.addPeak(x, 1500*barotropic_model::G, model.getDomain().getRadius()*0.5);
    testCase.calcInitCond(model);
    const SingleScalarField &gh = model.getGeopotentialHeight();
    // -------------------------------------------------------------------------
    // set initial condition for tracers
    ScalarField *q0; // reference tracer
    ScalarField *q1; // continuous tracer
    ScalarField *q2; // discontinuous tracer
    TimeLevelIndex<2> timeIdx;
    // reference tracer
    q.push_back(new ScalarField); q0 = q.back();
    q0->create("", "", "", model.getMesh(), CENTER);
    for (int i = 0; i < model.getMesh().getTotalNumGrid(CENTER); ++i) {
        (*q0)(timeIdx, i) = 1.0;
    }
    // continuous tracer
    q.push_back(new ScalarField); q1 = q.back();
    q1->create("", "", "", model.getMesh(), CENTER);
    for (int i = 0; i < model.getMesh().getTotalNumGrid(CENTER); ++i) {
        (*q1)(timeIdx, i) = gh(i);
    }
    // discontinous tracer
    q.push_back(new ScalarField); q2 = q.back();
    q2->create("", "", "", model.getMesh(), CENTER);
    for (int i = 0; i < model.getMesh().getTotalNumGrid(CENTER); ++i) {
        model.getMesh().getGridCoord(i, CENTER, x);
        if (x(0) > 140*RAD && x(0) < 160*RAD &&
            x(1) > 10*RAD  && x(1) < 40*RAD) {
            (*q2)(timeIdx, i) = 1;
        } else {
            (*q2)(timeIdx, i) = 0.1;
        }
    }
    // -------------------------------------------------------------------------
    AdvectionTestCase::calcInitCond(advectionManager);
}

void BarotropicTestCase::advance(double time,
                                 const TimeLevelIndex<2> &timeIdx) {
    if (timeIdx.isCurrentIndex()) {
        model.integrate(timeIdx, getStepSize());
    } else {
        model.integrate(timeIdx-1, getStepSize());
    }
    io.create(fileIdx);
    io.output<double, 1>(fileIdx, 1, &model.getGeopotentialHeight());
    io.close(fileIdx);
    for (int j = 0; j < model.getMesh().getNumGrid(1, velocity(0).getGridType(1)); ++j) {
        for (int i = 0; i < model.getMesh().getNumGrid(0, velocity(0).getGridType(0)); ++i) {
            velocity(0)(timeIdx, i, j) = model.getZonalWind()(timeIdx, i, j);
        }
    }
    for (int j = 0; j < model.getMesh().getNumGrid(1, velocity(1).getGridType(1)); ++j) {
        for (int i = 0; i < model.getMesh().getNumGrid(0, velocity(1).getGridType(0)); ++i) {
            velocity(1)(timeIdx, i, j) = model.getMeridionalWind()(timeIdx, i, j);
        }
    }
    if (timeIdx.isCurrentIndex()) {
        velocity.applyBndCond(timeIdx);
    } else {
        velocity.applyBndCond(timeIdx, UPDATE_HALF_LEVEL);
    }
}
    
}
