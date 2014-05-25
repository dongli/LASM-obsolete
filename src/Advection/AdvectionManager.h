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
    TracerManager tracerManager;
    Field<TracerMeshCell> tracerMeshCells;
    Regrid *regrid; //>! used to interpolate velocity onto tracers
    TimeLevels<vector<double>, 2> totalMass; //>! tracer species total mass
    // key parameters
    double alpha;     //>! control the location of new split parcels
    double beta1;     //>! control the merging weight along the major axis
    double beta2;     //>! control the merging weight vertical to the major axis
    bool isMassFixed; //>! control whether use mass fixer
private:
    // range search parameters
    Tree *cellTree;                 //>! tree data structure for mesh cells for
                                    //>! avoiding rebuild of tree each time
    mat cellCoords;                 //>! collection of cell space coordinates
    vector<size_t> cellCoordsMap;   //>! mapping for cells since tree building
                                    //>! will modify the order of cells
    // some array recording objects need to be processed
    int numLongTracer;
    vector<list<Tracer*>::iterator> longTracers;
    int numUnresolvedTracer;
    vector<list<Tracer*>::iterator> unresolvedTracers;
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
              const geomtk::ConfigManager &configManager);

    const Field<TracerMeshCell>& getTracerMeshCells() const {
        return tracerMeshCells;
    }

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
     *  @param fileName   the file name.
     *  @param newTimeIdx the new time level index.
     */
    void output(const string &fileName, const TimeLevelIndex<2> &newTimeIdx);

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
     *  @param V          the velocity field.
     */
    void integrate_RK4(double dt, const TimeLevelIndex<2> &oldTimeIdx,
                       const VelocityField &V);

    void embedTracersIntoMesh(const TimeLevelIndex<2> &timeIdx);

    void connectTracerAndMesh(const TimeLevelIndex<2> &timeIdx,
                              list<Tracer*>::iterator &tracer);

    /**
     *  Prepare the bidirectional remapping between tracers and mesh. Find out
     *  the mesh cells that a tracer will affected and calculate the weights.
     *
     *  @param timeIdx the time level index.
     */
    void connectTracersAndMesh(const TimeLevelIndex<2> &timeIdx);

    void splitTracers(const TimeLevelIndex<2> &timeIdx);
    
    void mergeTracers(const TimeLevelIndex<2> &timeIdx);

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
    void remapTracersToMesh(const TimeLevelIndex<2> &timeIdx);

    void correctTotalMassOnMesh(const TimeLevelIndex<2> &timeIdx);
};
}

#endif
