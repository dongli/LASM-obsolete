#ifndef __LASM_MeshAdaptor_h__
#define __LASM_MeshAdaptor_h__

#include "lasm_commons.h"

namespace lasm {

class Tracer;

class MeshAdaptor
{
private:
    const Domain *domain;
    const Mesh *mesh;
    vector<ScalarField*> densities;
    vector<ScalarField*> masses;
    vector<double> volumes;
    vector<int> numConnectedTracer;
    vector<vector<Tracer*> > connectedTracers;
    vector<vector<double> > remapWeights;
    vector<int> numContainedTracer;
    vector<vector<Tracer*> > containedTracers;
public:
    MeshAdaptor();
    ~MeshAdaptor();

    void init(const Domain &domain, const Mesh &mesh);

    /**
     *  Register tracer by increasing several arrays.
     *
     *  @param name  the tracer name.
     *  @param units the tracer density units.
     *  @param brief the tracer brief.
     */
    void registerTracer(const string &name, const string &units, const string &brief);

    /**
     *  Reset all the densities and masses to zero.
     *
     *  @param timeIdx the time level index.
     */
    void resetSpecies(const TimeLevelIndex<2> &timeIdx);
    
    /**
     *  Reset the number of connected tracers of all the cells to zero.
     */
    void resetConnectedTracers();

    /**
     *  Reset the number of contained tracers of all the cells to zero.
     */
    void resetContainedTracers();

    /**
     *  Return the meshed tracer density arrays.
     *
     *  @return The meshed tracer densities.
     */
    vector<ScalarField*>& getDensities() { return densities; }

    /**
     *  Return the meshed tracer mass arrays.
     *
     *  @return The meshed tracer masses.
     */
    vector<ScalarField*>& getMasses() { return masses; }

    /**
     *  Return the center grid coordinate of a cell.
     *
     *  @param i the cell index.
     *
     *  @return The center grid coordinate.
     */
    const SpaceCoord& getCoord(int i) const { return mesh->getGridCoord(CENTER, i); }

    /**
     *  Return the volume of a cell.
     *
     *  @param i the cell index.
     *
     *  @return The cell volume.
     */
    double getVolume(int i) const { return volumes[i]; }

    /**
     *  Return the number of connected tracers of a cell.
     *
     *  @param i the cell index.
     *
     *  @return The connected tracer number.
     */
    int getNumConnectedTracer(int i) const { return numConnectedTracer[i]; }

    /**
     *  Return the connected tracer pointer array of a cell. Note the array size
     *  may be larger than the actual connected tracer number.
     *
     *  @param i the cell index.
     *
     *  @return The connected tracer pointer array.
     */
    const vector<Tracer*>& getConnectedTracers(int i) const { return connectedTracers[i]; }

    /**
     *  Return the remapping weight of a cell with a tracer.
     *
     *  @param i      the cell index.
     *  @param tracer the tracer pointer.
     *
     *  @return The remapping weight.
     */
    double getRemapWeight(int i, Tracer *tracer) const;

    /**
     *  Return the number of contained tracers of a cell.
     *
     *  @param i the cell index.
     *
     *  @return The contained tracer number.
     */
    int getNumContainedTracer(int i) const { return numContainedTracer[i]; }

    /**
     *  Return the contained tracer pointer array of a cell. Note the array size
     *  may be larger than the actual contained tracer number.
     *
     *  @param i the cell index.
     *
     *  @return The contained tracer pointer array.
     */
    const vector<Tracer*>& getContainedTracers(int i) const { return containedTracers[i]; }

    /**
     *  Return the tracer density of a species in a cell.
     *
     *  @param timeIdx the time level index.
     *  @param s       the tracer species index.
     *  @param i       the cell index.
     *
     *  @return The tracer density.
     */
    double& getDensity(const TimeLevelIndex<2> &timeIdx, int s, int i) {
        return (*densities[s])(timeIdx, i);
    }

    /**
     *  Return the tracer mass of a species in a cell.
     *
     *  @param timeIdx the time level index.
     *  @param s       the tracer species index.
     *  @param i       the cell index.
     *
     *  @return The tracer mass.
     */
    double& getMass(const TimeLevelIndex<2> &timeIdx, int s, int i) {
        return (*masses[s])(timeIdx, i);
    }

    /**
     *  Connect a tracer with a cell.
     *
     *  @param i      the cell index
     *  @param tracer the tracer pointer.
     *  @param weight the remapping weight.
     */
    void connectTracer(int i, Tracer *tracer, double weight);

    /**
     *  Contain a tracer into a cell.
     *
     *  @param i      the cell index.
     *  @param tracer the tracer pointer.
     */
    void containTracer(int i, Tracer *tracer);

    /**
     *  Check if a tracer is connected with a cell.
     *
     *  @param i      the cell index.
     *  @param tracer the tracer pointer.
     *
     *  @return The boolean result.
     */
    bool isTracerConnected(int i, Tracer *tracer);
};

} // lasm

#endif // __LASM_MeshAdaptor_h__