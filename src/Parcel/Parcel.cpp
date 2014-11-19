#include "Parcel.h"
#include "ShapeFunction.h"

namespace lasm {

Parcel::Parcel(int numDim) {
    for (int l = 0; l < _q.numLevel(); ++l) {
        _q.level(l) = new SpaceCoord(numDim);
        _H.level(l) = new mat(numDim, numDim);
        _invH.level(l) = new mat(numDim, numDim);
        _idx.level(l) = new MeshIndex(numDim);
    }
}

Parcel::~Parcel() {
    for (int l = 0; l < _q.numLevel(); ++l) {
        delete _q.level(l);
        delete _H.level(l);
        delete _invH.level(l);
        delete _idx.level(l);
    }
}
    
Parcel& Parcel::operator=(const Parcel &other) {
    if (this != &other) {
        for (int l = 0; l < _q.numLevel(); ++l) {
            *(_q.level(l)) = *(other._q.level(l));
            *(_H.level(l)) = *(other._H.level(l));
            _detH.level(l) = other._detH.level(l);
            *(_invH.level(l)) = *(other._invH.level(l));
            _shapeSize.level(l) = other._shapeSize.level(l);
            *(_idx.level(l)) = *(other._idx.level(l));
        }
    }
    return *this;
}
    
void Parcel::calcSpaceCoord(const Domain &domain,
                            const TimeLevelIndex<2> &timeIdx,
                            const BodyCoord &y, SpaceCoord &x) {
#if defined USE_SPHERE_DOMAIN
    // In sphere domain, we calculate deformation matrix stuffs on local
    // stereographic projection of tracer centroid.
    x() = (*_H.level(timeIdx))*y();
    domain.projectBack(geomtk::SphereDomain::STEREOGRAPHIC,
                       *_q.level(timeIdx), x, x());
#elif defined USE_CARTESIAN_DOMAIN
    x() = (*_q.level(timeIdx))()+(*_H.level(timeIdx))*y();
    domain.constrain(x);
#endif
}

void Parcel::calcBodyCoord(const Domain &domain,
                           const TimeLevelIndex<2> &timeIdx,
                           const SpaceCoord &x, BodyCoord &y) {
#ifdef USE_SPHERE_DOMAIN
    domain.project(geomtk::SphereDomain::STEREOGRAPHIC,
                   *_q.level(timeIdx), x, y());
    y() = (*_invH.level(timeIdx))*y();
#elif defined USE_CARTESIAN_DOMAIN
    y() = (*_invH.level(timeIdx))*domain.diffCoord(x, (*_q.level(timeIdx)));
#endif
}
    
double Parcel::shapeFunction(const TimeLevelIndex<2> &timeIdx,
                             const BodyCoord &y) {
    double f;
    ShapeFunction::evalFunc(y, f);
    return f;
}

void Parcel::updateShapeSize(const Domain &domain,
                             const TimeLevelIndex<2> &timeIdx) {
    BodyCoord y(domain.numDim());
    SpaceCoord x(domain.numDim());
    for (int m = 0; m < domain.numDim(); ++m) {
        y() = (*_invH.level(timeIdx))*(*_H.level(timeIdx))*_V.col(m);
        calcSpaceCoord(domain, timeIdx, y, x);
        _shapeSize.level(timeIdx)[m] = domain.calcDistance(x, *_q.level(timeIdx));
    }
}

}
