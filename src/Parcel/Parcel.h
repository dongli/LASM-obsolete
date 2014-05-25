#ifndef __LASM_Parcel__
#define __LASM_Parcel__

#include "lasm_commons.h"

namespace lasm {

/**
 *  This class specifies the parcel, which is the computational unit of LASM.
 */
class Parcel {
protected:
    int ID;
    TimeLevels<SpaceCoord*, 2> q;           //>! centroid coordinate
    TimeLevels<mat*, 2> H;                  //>! deformation matrix
    TimeLevels<double, 2> detH;             //>! determinant of H
    TimeLevels<mat*, 2> invH;               //>! inversion of H
    TimeLevels<vec::fixed<2>, 2> shapeSize; //>! parcel shape size
    TimeLevels<MeshIndex*, 2> idx;          //>! parcel mesh index
    
    mat U, V;   //>! SVD decomposed matrices
    vec S;      //>! SVD decomposed matrices
public:
    Parcel(int numDim);
    virtual ~Parcel();

    Parcel& operator=(const Parcel &other);
    
    /**
     *  Set the ID of tracer.
     *
     *  @param ID the given ID.
     */
    virtual void setID(int ID) { this->ID = ID; }

    /**
     *  Get the ID of tracer.
     *
     *  @return The ID.
     */
    virtual int getID() const { return ID; }

    /**
     *  Get the spatial coordinate of tracer.
     *
     *  @param timeIdx the time level index.
     *
     *  @return The spatial coordinate.
     */
    SpaceCoord& getX(const TimeLevelIndex<2> &timeIdx) {
        return *(q.getLevel(timeIdx));
    }

    /**
     *  Get the deformation matrix of tracer.
     *
     *  @param timeIdx the time level index.
     *
     *  @return The deformation matrix.
     */
    mat& getH(const TimeLevelIndex<2> &timeIdx) {
        return *(H.getLevel(timeIdx));
    }

    double getDetH(const TimeLevelIndex<2> &timeIdx) const {
        return detH.getLevel(timeIdx);
    }
    
    double& getDetH(const TimeLevelIndex<2> &timeIdx) {
        return detH.getLevel(timeIdx);
    }

    /**
     *  Get the inverse deformation matrix of tracer.
     *
     *  @param timeIdx the time level index.
     *
     *  @return The inverse deformation matrix.
     */
    mat& getInvH(const TimeLevelIndex<2> &timeIdx) {
        return *(invH.getLevel(timeIdx));
    }

    /**
     *  Transform body coordinate into spatial coordinate.
     *
     *  @param domain  the domain used for differencing coordinates.
     *  @param timeIdx the time level index.
     *  @param y       the body coordinate.
     *  @param x       the spatial coordinate.
     */
    void getSpaceCoord(const Domain &domain,
                       const TimeLevelIndex<2> &timeIdx,
                       const BodyCoord &y, SpaceCoord &x);
    
    /**
     *  Transform spatial coordinate into body coordinate.
     *
     *  @param domain  the domain used for differencing coordinates.
     *  @param timeIdx the time level index.
     *  @param x       the spatial coordinate.
     *  @param y       the body coordinate.
     */
    void getBodyCoord(const Domain &domain,
                      const TimeLevelIndex<2> &timeIdx,
                      const SpaceCoord &x, BodyCoord &y);

    /**
     *  Get the shape function value for a given body coordinate.
     *
     *  @param timeIdx the time level index.
     *  @param y       the body coordinate.
     *
     *  @return The shape function value.
     */
    double getShapeFunction(const TimeLevelIndex<2> &timeIdx,
                            const BodyCoord &y);

    mat& getU() { return U; }
    
    mat& getV() { return V; }
    
    vec& getS() { return S; }
    
    void updateShapeSize(const Domain &domain,
                         const TimeLevelIndex<2> &timeIdx);

    const vec::fixed<2>& getShapeSize(const TimeLevelIndex<2> &timeIdx) const {
        return shapeSize.getLevel(timeIdx);
    }
};

}

#endif // __LASM_Parcel__
