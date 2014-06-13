#ifndef __LASM_AdvectionManager__
#define __LASM_AdvectionManager__

#include "TracerManager.h"
#include "TracerMeshCell.h"

namespace lasm {

/**
 *  This class specifies the manager of linear advection.
 */
class AdvectionManager {
protected:
    const Domain *domain;
    const Mesh *mesh;
    const TimeManager *timeManager;
    TracerManager tracerManager;
    Field<TracerMeshCell> tracerMeshCells;
    Regrid *regrid;
    // -------------------------------------------------------------------------
    // key parameters
    double filamentLimit;   //>! control the tracer filament degree
    double radialMixing;    //>! control the radial mixing degree
    double lateralMixing;   //>! control the lateral mixing degree
    bool isMassFixed;       //>! control whether use mass fixer
    // -------------------------------------------------------------------------
    StampString *outputFileFormat;
    TimeLevels<vector<double>, 2> totalMass;
    // -------------------------------------------------------------------------
    // range search parameters
    Tree *cellTree;                 //>! tree data structure for mesh cells for
                                    //>! avoiding rebuild of tree each time
    mat cellCoords;                 //>! collection of cell space coordinates
    vector<size_t> cellCoordsMap;   //>! mapping for cells since tree building
                                    //>! will modify the order of cells
    // -------------------------------------------------------------------------
    // some array recording objects need to be processed
    int numFilamentTracer;
    vector<list<Tracer*>::iterator> filamentTracers;
    int numVoidCell;
    vector<TracerMeshCell*> voidCells;
public:
    AdvectionManager();
    ~AdvectionManager();

    /**
     *  Initialize advection manager.
     *
     *  @param domain        the spatial domain.
     *  @param mesh          the mesh.
     *  @param configManager the configuration manager.
     */
    void init(const Domain &domain, const Mesh &mesh,
              const ConfigManager &configManager,
              const TimeManager &timeManager);


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
    void input(const TimeLevelIndex<2> &timeIdx, vector<ScalarField*> &q);

    /**
     *  Output tracers on old time level into netCDF file.
     *
     *  @param timeIdx the time level index.
     */
    void output(const TimeLevelIndex<2> &timeIdx);

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
     *  @param V          the velocity field.
     */
    void advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                 const VelocityField &V);
private:
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
     *  @param oldTimeIdx the old time level index.
     *  @param velocity   the velocity field.
     */
    void integrate_RK4(double dt, const TimeLevelIndex<2> &oldTimeIdx,
                       const VelocityField &velocity);

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
    void connectTracerAndMesh(const TimeLevelIndex<2> &timeIdx,
                              list<Tracer*>::iterator &tracer);

    /**
     *  Prepare the bidirectional remapping between tracers and mesh. Find out
     *  the mesh cells that a tracer will affected and calculate the weights.
     *
     *  @param timeIdx the time level index.
     */
    void connectTracersAndMesh(const TimeLevelIndex<2> &timeIdx);

    void checkTracerShapes(const TimeLevelIndex<2> &timeIdx,
                           const VelocityField &velocity);

    /**
     *  Handle the filament tracers.
     *
     *  @param timeIdx the time level index.
     */
    void handleFilamentTracers(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Handle void cells.
     *
     *  @param timeIdx the time level index.
     */
    void handleVoidCells(const TimeLevelIndex<2> &timeIdx);

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
    void remapTracersToMesh(const TimeLevelIndex<2> &timeIdx,
                            const VelocityField *velocity = NULL);

    /**
     *  Correct the total mass when remapping tracer density.
     *
     *  @param timeIdx the time level index.
     */
    void correctTotalMassOnMesh(const TimeLevelIndex<2> &timeIdx);

    void recordTracer(Tracer::TracerType type, list<Tracer*>::iterator &tracer);
};
}

#endif
