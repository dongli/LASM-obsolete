#ifndef __LASM_Tracer__
#define __LASM_Tracer__

#include "lasm_commons.h"
#include "Parcel.h"

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
        GOOD_SHAPE, EXTREME_FILAMENTATION, NOT_RESOLVED, POOR_APPROXIMATION
    };
    Tracer *father;
protected:
    vec density;        //>! species density array
    vec mass;           //>! species mass array
    TracerSkeleton *skeleton;
    TracerType type;
    /**
     *  Remapping parameters
     */
    int numConnectedCell;
    vector<TracerMeshCell*> connectedCells;

    TracerMeshCell *hostCell;
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
    void calcSpeciesDensity(const TimeLevelIndex<2> &timeIdx, int s) {
        density[s] = mass[s]/detH.getLevel(timeIdx);
    }

    double& getSpeciesDensity(int speciesIdx) {
#ifndef NDEBUG
        if (speciesIdx >= density.size()) {
            REPORT_ERROR("Species index " << speciesIdx << " exceeds range [0," <<
                         density.size()-1 << "]!");
        }
#endif
        return density[speciesIdx];
    }

    double getSpeciesDensity(int speciesIdx) const {
#ifndef NDEBUG
        if (speciesIdx >= density.size()) {
            REPORT_ERROR("Species index " << speciesIdx << " exceeds range [0," <<
                         density.size()-1 << "]!");
        }
#endif
        return density[speciesIdx];
    }

    /**
     *  Calculate species mass from density.
     *
     *  @param timeIdx the time level index.
     *  @param s       the species index.
     */
    void calcSpeciesMass(const TimeLevelIndex<2> &timeIdx, int s) {
        mass[s] = density[s]*detH.getLevel(timeIdx);
    }

    double& getSpeciesMass(int speciesIdx) {
        return mass[speciesIdx];
    }

    double getSpeciesMass(int speciesIdx) const {
        return mass[speciesIdx];
    }

    void resetSpecies() {
        for (int s = 0; s < density.size(); ++s) {
            density[s] = 0.0;
        }
    }

    Tracer& operator=(const Tracer &other);

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
     *  @param cell   the mesh cell to be connected.
     *  @param weight the remapping weight.
     */
    void connect(TracerMeshCell *cell, double weight);

    /**
     *  Get the connected mesh cells.
     *
     *  @return The connected mesh cell list.
     */
    vector<TracerMeshCell*>& getConnectedCells() { return connectedCells; }

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
     *  @param cell the host cell.
     */
    void setHostCell(TracerMeshCell *cell) { hostCell = cell; }

    /**
     *  Get the host cell where the tracer is.
     *
     *  @return The host cell.
     */
    TracerMeshCell* getHostCell() const { return hostCell; }

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

    /**
     *  Reset deformation matrix from given major axis vertex and scale matrix.
     *
     *  @param domain  the spatial domain.
     *  @param mesh    the model mesh.
     *  @param timeIdx the time level index.
     *  @param x       the major axis vertex spatial coordinate.
     *  @param S       the scale matrix.
     */
    void resetDeformMatrix(const Domain &domain,
                           const Mesh &mesh,
                           const TimeLevelIndex<2> &timeIdx,
                           const SpaceCoord &x, const vec &S);

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
     *  @param timeIdx the time level index.
     *  @param domain  the spatial domain.
     *  @param file    the text file in class std::ofstream.
     *  @param idx     the index to distinguish tracer in the file.
     */
    void dump(const TimeLevelIndex<2> &timeIdx, const Domain &domain,
              std::ofstream &file, int idx);

    /**
     *  Dump the centroid, shape, neighbor cells and neighbor tracers into file
     *  using by external NCL script.
     *
     *  @param timeIdx the time level index.
     *  @param domain  the spatial domain.
     */
    void dump(const TimeLevelIndex<2> &timeIdx, const Domain &domain);
};

}

#endif // __LASM_Tracer__
