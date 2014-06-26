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
        GOOD_SHAPE, NEED_MIXING
    };
    double actualFilamentLimit;
    double actualLateralMixing;
protected:
    vec density;        //>! species density array
    vec mass;           //>! species mass array
    TracerSkeleton *skeleton;
    TracerType type;
    BodyCoord *vy;  //>! long axis vertex body coordinate
    SpaceCoord *vx; //>! long axis vertex space coordinate
    /**
     *  Remapping parameters
     */
    int numConnectedCell;
    vector<int> connectedCellIdxs;

    int hostCellIdx;
public:
    Tracer(int numDim);
    virtual ~Tracer();
    
    /**
     *  Add a species.
     */
    void addSpecies() {
        density.resize(density.size()+1);
        density[density.size()-1] = 0;
        mass.resize(mass.size()+1);
        mass[mass.size()-1] = 0;
    }

    /**
     *  Calculate species density from mass.
     *
     *  @param timeIdx the time level index.
     *  @param s       the species index.
     */
    void calcDensity(const TimeLevelIndex<2> &timeIdx, int s) {
        density[s] = mass[s]/detH.getLevel(timeIdx);
    }

    double& getDensity(int s) { return density[s]; }

    double getDensity(int s) const { return density[s]; }

    /**
     *  Calculate species mass from density.
     *
     *  @param timeIdx the time level index.
     *  @param s       the species index.
     */
    void calcMass(const TimeLevelIndex<2> &timeIdx, int s) {
        mass[s] = density[s]*detH.getLevel(timeIdx);
    }

    double& getMass(int s) { return mass[s]; }

    double getMass(int s) const { return mass[s]; }

    void resetSpecies() {
        density.zeros();
        mass.zeros();
    }

    Tracer& operator=(const Tracer &other);

    const BodyCoord& getLongAxisVertexBodyCoord() const { return *vy; }

    SpaceCoord& getLongAxisVertexSpaceCoord() { return *vx; }

    /**
     *  Get the mesh index of tracer.
     *
     *  @param timeIdx the time level index.
     *
     *  @return The mesh index.
     */
    MeshIndex& getMeshIndex(const TimeLevelIndex<2> &timeIdx) {
        return *idx.getLevel(timeIdx);
    }

    /**
     *  Get the skeleton of tracer.
     *
     *  @return The skeleton object.
     */
    TracerSkeleton& getSkeleton() { return *skeleton; }

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
    const vector<int>& getConnectedCells() const { return connectedCellIdxs; }

    /**
     *  Get the number of connected cells. It should be used instead of the size
     *  method of the connected cell array.
     *
     *  @return The connected cell number.
     */
    int getNumConnectedCell() const { return numConnectedCell; }

    /**
     *  Set the host cell where the tracer is.
     *
     *  @param cell the host cell index.
     */
    void setHostCell(int i) { hostCellIdx = i; }

    /**
     *  Get the host cell where the tracer is.
     *
     *  @return The host cell.
     */
    int getHostCell() const { return hostCellIdx; }

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
