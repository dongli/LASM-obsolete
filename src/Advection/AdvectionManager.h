#ifndef __LASM_AdvectionManager__
#define __LASM_AdvectionManager__

#include "TracerManager.h"
#include "MeshAdaptor.h"

namespace lasm {

class AdvectionTestCase;

/**
 *  This class specifies the manager of linear advection.
 */
class AdvectionManager {
protected:
    const Domain *domain;
    const Mesh *mesh;
    TracerManager tracerManager;
    MeshAdaptor meshAdaptor;
    Regrid *regrid;
    /**
     *  Interparcel mixing control parameters
     */
    double filamentLimit;       //>! control the tracer filament degree
    double strictFilamentLimit; //>! strict value of filamentLimit
    double radialMixing;        //>! control the radial mixing degree
    double lateralMixing;       //>! control the lateral mixing degree
    double strictLateralMixing; //>! strict value of lateralMixing
    double shrinkFactor;        //>! control the mixed parcel shrinking degree
    double disorderDegreeLimit;
    bool isMassFixed;           //>! control whether use mass fixer
    TimeLevels<vector<double>, 2> totalMass;
    /**
     *  Neighbor cell searching variables
     */
    Tree *cellTree;                 //>! tree data structure for mesh cells for
                                    //>! avoiding rebuild of tree each time
    mat cellCoords;                 //>! collection of cell space coordinates
    vector<size_t> cellCoordsMap;   //>! mapping for cells since tree building
                                    //>! will modify the order of cells
    // some array recording objects need to be processed
    int numMixedTracer;
    vector<Tracer*> mixedTracers;
    int numVoidCell;
    vector<int> voidCells;
public:
    AdvectionManager();
    ~AdvectionManager();

    void init(const Domain &domain, const Mesh &mesh,
              const ConfigManager &configManager);

    vector<ScalarField*>& density() { return meshAdaptor.density(); }

    /**
     *  Register a tracer species.
     *
     *  @param name  the name of the tracer.
     *  @param units the units of the tracer.
     *  @param brief the brief about the tracer.
     */
    void registerTracer(const string &name, const string &units,
                        const string &brief);

    /**
     *  Input one tracer species from given scalar field.
     *
     *  @param timeIdx the time level index.
     *  @param q       the input scalar field.
     */
    void input(const TimeLevelIndex<2> &timeIdx, double *q);

    void restart(const string &fileName);

    /**
     *  Output tracers on old time level into netCDF file.
     *
     *  @param timeIdx the time level index.
     *  @param ncId    the netCDF output file ID.
     */
    void output(const TimeLevelIndex<2> &timeIdx, int ncId);

    /**
     *  Diagnose total tracer mass.
     *
     *  @param timeIdx the time level index.
     */
    void diagnose(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Advect tracers one time step forward, and remap tracers onto mesh cells.
     *
     *  @param dt         the time step size.
     *  @param oldTimeIdx the new time level index.
     *  @param velocity   the velocity field.
     */
    void advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                 const VelocityField &velocity);

    void advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                 const AdvectionTestCase &testCase);

    /**
     *  Calculate the total mass of every tracer species, and store the results
     *  into tracerMeshCells.
     *
     *  @param timeIdx the time level index.
     */
    void calcTotalMass(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Integrate the advection equations by using 4th-order Runge-Kutta method.
     *
     *  @param dt         the time step size.
     *  @param newTimeIdx the new time level index.
     *  @param velocity   the velocity field.
     */
    void integrate_RK4(double dt, const TimeLevelIndex<2> &newTimeIdx,
                       const VelocityField &velocity);

    void integrate_RK4(double dt, const TimeLevelIndex<2> &newTimeIdx,
                       const AdvectionTestCase &testCase);

    /**
     *  Find out in which cell the tracer is.
     *
     *  @param timeIdx the time level index.
     */
    void embedTracersIntoMesh(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Connect one tracer to the mesh grids.
     *
     *  @param timeIdx the time step size.
     *  @param tracer  the tracer iterator.
     */
    void connectTracerAndMesh(const TimeLevelIndex<2> &timeIdx, Tracer *tracer);

    /**
     *  Prepare the bidirectional remapping between tracers and mesh. Find out
     *  the mesh cells that a tracer will affected and calculate the weights.
     *
     *  @param timeIdx the time level index.
     */
    void connectTracersAndMesh(const TimeLevelIndex<2> &timeIdx);

    void checkTracerShapes(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Mix tracer with its surrounding tracers.
     *
     *  @param timeIdx the time level index.
     */
    void mixTracers(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Fill void cells by tracer densities on their neighbors.
     *
     *  @param timeIdx the time level index.
     */
    void fillVoidCells(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Remap the tracer mass from mesh cells to tracers.
     *
     *  @param timeIdx the time level index.
     */
    void remapMeshToTracers(const TimeLevelIndex<2> &timeIdx);
    
    /**
     *  Remap the tracer mass from tracers to mesh cells.
     *
     *  @param timeIdx the time level index.
     */
    void remapTracersToMesh(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Correct the total mass when remapping tracer density.
     *
     *  @param timeIdx the time level index.
     */
    void correctTotalMassOnMesh(const TimeLevelIndex<2> &timeIdx);

    void recordTracer(Tracer::TracerType type, Tracer *tracer);
};
}

#endif
