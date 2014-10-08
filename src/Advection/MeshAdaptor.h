#ifndef __LASM_MeshAdaptor_h__
#define __LASM_MeshAdaptor_h__

#include "lasm_commons.h"

namespace lasm {

class Tracer;
class TracerManager;

class MeshAdaptor
{
private:
    const Domain *domain;
    const Mesh *mesh;
    vector<ScalarField*> _density;
    vector<ScalarField*> _mass;
    vector<int> _numConnectedTracer;
    vector<vector<Tracer*> > _connectedTracers;
    vector<vector<double> > _remapWeights;
    vector<int> _numContainedTracer;
    vector<vector<Tracer*> > _containedTracers;
public:
    MeshAdaptor();
    ~MeshAdaptor();

    void init(const Domain &domain, const Mesh &mesh);

    void input(const string &fileName, const TracerManager &tracerManager);

    void registerTracer(const string &name, const string &units, const string &brief);

    void resetSpecies(const TimeLevelIndex<2> &timeIdx);

    void resetConnectedTracers();

    void resetContainedTracers();

    const SpaceCoord& coord(int i) const { return mesh->getGridCoord(CENTER, i); }

    double volume(int i) const { return mesh->getCellVolume(i); }

    vector<ScalarField*>& density() { return _density; }

    double& density(const TimeLevelIndex<2> &timeIdx, int s, int i) {
        return (*_density[s])(timeIdx, i);
    }

    double& mass(const TimeLevelIndex<2> &timeIdx, int s, int i) {
        return (*_mass[s])(timeIdx, i);
    }

    void connectTracer(int i, Tracer *tracer, double weight);

    int numConnectedTracer(int i) const { return _numConnectedTracer[i]; }

    const vector<Tracer*>& connectedTracers(int i) const {
        return _connectedTracers[i];
    }

    void containTracer(int i, Tracer *tracer);

    int numContainedTracer(int i) const { return _numContainedTracer[i]; }

    const vector<Tracer*>& containedTracers(int i) const {
        return _containedTracers[i];
    }

    double remapWeight(int i, Tracer *tracer) const;

    bool isTracerConnected(int i, Tracer *tracer);
};

} // lasm

#endif // __LASM_MeshAdaptor_h__