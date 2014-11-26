#include "TerminatorChemistryTestCase.h"

namespace lasm {

TerminatorChemistryTestCase::TerminatorChemistryTestCase() {
    period = 12*TimeUnit::DAYS;
    _stepSize = 1800*TimeUnit::SECONDS;
    _subcycledStepSize = 1800*TimeUnit::SECONDS;
    _numSubcycledStep = _stepSize/_subcycledStepSize;
    REPORT_ONLINE;
}

TerminatorChemistryTestCase::~TerminatorChemistryTestCase() {
    delete _domain;
    delete _mesh;
    REPORT_OFFLINE;
}

void TerminatorChemistryTestCase::init(const ConfigManager &configManager,
                                       TimeManager &timeManager) {
    // Initialize domain.
    _domain = new geomtk::SphereDomain(2);
    _domain->radius() = 6.3172e6;
    // Initialize mesh.
    _mesh = new geomtk::RLLMesh(*_domain);
    int numLon = 360, numLat = 181;
    if (configManager.hasKey("test_case", "num_lon")) {
        configManager.getValue("test_case", "num_lon", numLon);
    }
    if (configManager.hasKey("test_case", "num_lat")) {
        configManager.getValue("test_case", "num_lat", numLat);
    }
    _mesh->init(numLon, numLat);
    // Call super class initialization.
    AdvectionTestCase::init(configManager, timeManager);
    // Initialize velocity.
    if (!useAnalyticalVelocity) {
        velocity.create(mesh(), true, HAS_HALF_LEVEL);
    }
}

Time TerminatorChemistryTestCase::startTime() const {
    Time time;
    return time;
}

Time TerminatorChemistryTestCase::endTime() const {
    Time time = startTime()+period;
    return time;
}

void TerminatorChemistryTestCase::advancePhysics(double time,
                                                 const TimeLevelIndex<2> &timeIdx,
                                                 AdvectionManager &advectionManager) {
    // Calculate tracer density tendencies.
#if defined LASM_EVALUATE_TENDENCY_ON_MESH
    vector<SingleScalarField*> &tendency = advectionManager.tendency();
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double cl = (*(*density)[1])(timeIdx, i)/(*(*density)[0])(timeIdx, i);
        double cl2 = (*(*density)[2])(timeIdx, i)/(*(*density)[0])(timeIdx, i);
        calcTendency(x(0), x(1), cl, cl2,
                     (*tendency[1])(i), (*tendency[2])(i));
        (*tendency[1])(i) *= (*(*density)[0])(timeIdx, i)/numSubcycledStep();
        (*tendency[2])(i) *= (*(*density)[0])(timeIdx, i)/numSubcycledStep();
    }
#elif defined LASM_EVALUATE_TENDENCY_ON_PARCEL
    if (numSubcycledStep() != 1) {
        REPORT_ERROR("When evaluating tendency on parcels, the subcycling should be turned off!");
    }
    vector<Tracer*> &tracers = advectionManager.tracers();
    for (int i = 0; i < tracers.size(); ++i) {
        const SpaceCoord &x = tracers[i]->x(timeIdx);
        double cl = tracers[i]->density(1)/tracers[i]->density(0);
        double cl2 = tracers[i]->density(2)/tracers[i]->density(0);
        double dCl, dCl2;
        calcTendency(x(0), x(1), cl, cl2,
                     tracers[i]->tendency(1), tracers[i]->tendency(2));
        tracers[i]->tendency(1) *= tracers[i]->density(0);
        tracers[i]->tendency(2) *= tracers[i]->density(0);
    }
#endif
}

void TerminatorChemistryTestCase::calcInitCond(AdvectionManager &advectionManager) {
    // Register tracer species.
    advectionManager.registerTracer("q0", "N/A", "air");
    advectionManager.registerTracer("q1", "N/A", "Cl");
    advectionManager.registerTracer("q2", "N/A", "Cl2");
    density = &advectionManager.density();
    AdvectionTestCase::registerDefaultOutput();
    // Calculate initial conditions for each tracer species.
    double *q = new double[3*mesh().totalNumGrid(CENTER, 2)];
    int l = 0;
    // - air
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        q[l++] = 1.0;
    }
    // - Cl
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double k1, k2;
        calcChemicalReactionRate(x(0), x(1), k1, k2);
        double r = k1/(4*k2);
        double det = sqrt(r*r+2*r*cly0);
        q[l++] = det-r;
    }
    // - Cl2
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double k1, k2;
        calcChemicalReactionRate(x(0), x(1), k1, k2);
        double r = k1/(4*k2);
        double det = sqrt(r*r+2*r*cly0);
        q[l++] = cly0*0.5-(det-r)*0.5;
    }
    // Propagate initial conditions to advection manager.
    TimeLevelIndex<2> timeIdx;
    advectionManager.input(timeIdx, q);
    delete [] q;
}

void TerminatorChemistryTestCase::calcChemicalReactionRate(double lon, double lat,
                                                           double &k1, double &k2) {
    k1 = fmax(0, sin(lat)*sin(k1CenterLat)+cos(lat)*cos(k1CenterLat)*cos(lon-k1CenterLon));
    k2 = 1;
}

void TerminatorChemistryTestCase::calcTendency(double lon, double lat, double cl,
                                               double cl2, double &dCl, double &dCl2) {
    double k1, k2;

    calcChemicalReactionRate(lon, lat, k1, k2);

    double r = k1/(4*k2);
    double cly = cl+2*cl2;

    double det = sqrt(r*r+2*r*cly);
    double expdt = exp(-4*k2*det*stepSize());

    double el;
    if (fabs(det*k2*stepSize()) > 1e-16) {
        el = (1-expdt)/det/stepSize();
    } else {
        el = 4*k2;
    }

    dCl = -el*(cl-det+r)*(cl+det+r)/(1+expdt+stepSize()*el*(cl+r));
    dCl2 = -dCl*0.5;
}

}
