#ifndef __LASM_Tracer__
#define __LASM_Tracer__

#include "lasm_commons.h"
#include "Parcel.h"
#include "MeshAdaptor.h"

namespace lasm {

class TracerSkeleton;
class TracerMeshCell;

/**
 *  This class describes the tracer that is used to be advected by external wind
 *  field. It derives from parcel.
 */
class Tracer : public Parcel {
public:
    enum TracerType {
        GOOD_SHAPE, BAD_SHAPE
    };
protected:
    vec _density;  //>! species density array
    vec _mass;     //>! species mass array
    TracerSkeleton *_skeleton;
    double _filament;
    double _linearDegeneration;
    TracerType type;
    BodyCoord *vy;  //>! long axis vertex body coordinate
    SpaceCoord *vx; //>! long axis vertex space coordinate
    /**
     *  Remapping parameters
     */
    int _numConnectedCell;
    vector<int> _connectedCellIndices;

    int _hostCellIndex;
public:
    Tracer(int numDim);
    virtual ~Tracer();
    
    /**
     *  Add a species.
     */
    void addSpecies() {
        _density.resize(_density.size()+1);
        _density[_density.size()-1] = 0;
        _mass.resize(_mass.size()+1);
        _mass[_mass.size()-1] = 0;
    }

    double& density(int s) { return _density[s]; }

    double density(int s) const { return _density[s]; }

    double& mass(int s) { return _mass[s]; }

    double mass(int s) const { return _mass[s]; }

    void resetSpecies() {
        _density.zeros();
        _mass.zeros();
    }

    double filament() const { return _filament; }

    double linearDegeneration() const { return _linearDegeneration; }

    double& linearDegeneration() { return _linearDegeneration; }

    Tracer& operator=(const Tracer &other);

    const BodyCoord& longAxisVertexBodyCoord() const { return *vy; }

    SpaceCoord& longAxisVertexSpaceCoord() { return *vx; }

    /**
     *  Get the skeleton of tracer.
     *
     *  @return The skeleton object.
     */
    TracerSkeleton& skeleton() { return *_skeleton; }

    void setType(TracerType type) { this->type = type; }

    TracerType getType() const { return type; }

    /**
     *  Reset the connected cells to empty for later updating.
     */
    void resetConnectedCells();

    /**
     *  Connect the tracer with the given mesh cell.
     *
     *  @param i the cell index.
     */
    void connectCell(int i);

    /**
     *  Get the connected mesh cells.
     *
     *  @return The connected mesh cell index list.
     */
    const vector<int>& connectedCellIndices() const {
        return _connectedCellIndices;
    }

    /**
     *  Get the number of connected cells. It should be used instead of the size
     *  method of the connected cell array.
     *
     *  @return The connected cell number.
     */
    int numConnectedCell() const { return _numConnectedCell; }

    int& hostCellIndex() { return _hostCellIndex; }

    int hostCellIndex() const { return _hostCellIndex; }

    /**
     *  Update deformation matrix from tracer skeleton.
     *
     *  @param domain  the spatial domain.
     *  @param mesh    the model mesh.
     *  @param timeIdx the time level index.
     */
    void updateDeformMatrix(const Domain &domain,
                            const Mesh &mesh,
                            const TimeLevelIndex<2> &timeIdx);

    void updateDeformMatrix(const Domain &domain,
                            const TimeLevelIndex<2> &timeIdx);

    /**
     *  Reset tracer skeleton points from deformation matrix
     *
     *  @param domain  the spatial domain.
     *  @param mesh    the model mesh.
     *  @param timeIdx the time level index.
     */
    void resetSkeleton(const Domain &domain, const Mesh &mesh,
                       const TimeLevelIndex<2> &timeIdx);

    /**
     *  Dump the centroid, shape, neighbor cells and neighbor tracers into file
     *  using by external NCL script.
     *
     *  @param timeIdx     the time level index.
     *  @param domain      the spatial domain.
     *  @param meshAdaptor the mesh adaptor.
     *  @param file        the text file in class std::ofstream.
     *  @param idx         the index to distinguish tracer in the file.
     */
    void dump(const TimeLevelIndex<2> &timeIdx, const Domain &domain,
              const MeshAdaptor &meshAdaptor, ofstream &file, int idx);

    /**
     *  Dump the centroid, shape, neighbor cells and neighbor tracers into file
     *  using by external NCL script.
     *
     *  @param timeIdx     the time level index.
     *  @param domain      the spatial domain.
     *  @param meshAdaptor the mesh adaptor.
     */
    void dump(const TimeLevelIndex<2> &timeIdx, const Domain &domain,
              const MeshAdaptor &meshAdaptor);
};

}

#endif // __LASM_Tracer__
