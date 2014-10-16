#ifndef __LASM_Parcel__
#define __LASM_Parcel__

#include "lasm_commons.h"

namespace lasm {

/**
 *  This class specifies the parcel, which is the computational unit of LASM.
 */
class Parcel {
protected:
    int _ID;
    TimeLevels<SpaceCoord*, 2> _q;           //>! centroid coordinate
    TimeLevels<mat*, 2> _H;                  //>! deformation matrix
    TimeLevels<double, 2> _detH;             //>! determinant of H
    TimeLevels<mat*, 2> _invH;               //>! inversion of H
    TimeLevels<vec::fixed<2>, 2> _shapeSize; //>! parcel shape size
    TimeLevels<MeshIndex*, 2> _idx;           //>! parcel mesh index
    
    mat _U, _V; //>! SVD decomposed matrices
    vec _S;     //>! SVD decomposed matrices
public:
    Parcel(int numDim);
    virtual ~Parcel();

    Parcel& operator=(const Parcel &other);

    virtual int& ID() { return _ID; }

    virtual int ID() const { return _ID; }

    SpaceCoord& x(const TimeLevelIndex<2> &timeIdx) {
        return *(_q.level(timeIdx));
    }

    MeshIndex& meshIndex(const TimeLevelIndex<2> &timeIdx) {
        return *_idx.level(timeIdx);
    }

    mat& H(const TimeLevelIndex<2> &timeIdx) {
        return *(_H.level(timeIdx));
    }

    mat& U() { return _U; }

    mat& V() { return _V; }

    vec& S() { return _S; }

    double detH(const TimeLevelIndex<2> &timeIdx) const {
        return _detH.level(timeIdx);
    }
    
    double& detH(const TimeLevelIndex<2> &timeIdx) {
        return _detH.level(timeIdx);
    }

    mat& invH(const TimeLevelIndex<2> &timeIdx) {
        return *(_invH.level(timeIdx));
    }

    /**
     *  Transform body coordinate into spatial coordinate.
     *
     *  @param domain  the domain used for differencing coordinates.
     *  @param timeIdx the time level index.
     *  @param y       the body coordinate.
     *  @param x       the spatial coordinate.
     */
    void calcSpaceCoord(const Domain &domain,
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
    void calcBodyCoord(const Domain &domain,
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
    double shapeFunction(const TimeLevelIndex<2> &timeIdx, const BodyCoord &y);
    
    void updateShapeSize(const Domain &domain,
                         const TimeLevelIndex<2> &timeIdx);

    const vec::fixed<2>& shapeSize(const TimeLevelIndex<2> &timeIdx) const {
        return _shapeSize.level(timeIdx);
    }
};

}

#endif // __LASM_Parcel__
