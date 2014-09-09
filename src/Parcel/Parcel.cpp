#include "Parcel.h"
#include "ShapeFunction.h"

namespace lasm {

Parcel::Parcel(int numDim) {
    for (int l = 0; l < q.getNumLevel(); ++l) {
        q.getLevel(l) = new SpaceCoord(numDim);
        H.getLevel(l) = new mat(numDim, numDim);
        invH.getLevel(l) = new mat(numDim, numDim);
        idx.getLevel(l) = new MeshIndex(numDim);
    }
}

Parcel::~Parcel() {
    for (int l = 0; l < q.getNumLevel(); ++l) {
        delete q.getLevel(l);
        delete H.getLevel(l);
        delete invH.getLevel(l);
        delete idx.getLevel(l);
    }
}
    
Parcel& Parcel::operator=(const Parcel &other) {
    if (this != &other) {
        for (int l = 0; l < q.getNumLevel(); ++l) {
            *(q.getLevel(l)) = *(other.q.getLevel(l));
            *(H.getLevel(l)) = *(other.H.getLevel(l));
            detH.getLevel(l) = other.detH.getLevel(l);
            *(invH.getLevel(l)) = *(other.invH.getLevel(l));
            shapeSize.getLevel(l) = other.shapeSize.getLevel(l);
            *(idx.getLevel(l)) = *(other.idx.getLevel(l));
        }
    }
    return *this;
}
    
void Parcel::getSpaceCoord(const Domain &domain,
                           const TimeLevelIndex<2> &timeIdx,
                           const BodyCoord &y, SpaceCoord &x) {
#if defined USE_SPHERE_DOMAIN
    // In sphere domain, we calculate deformation matrix stuffs on local
    // stereographic projection of tracer centroid.
    x() = (*H.getLevel(timeIdx))*y();
    domain.projectBack(geomtk::SphereDomain::STEREOGRAPHIC,
                       *q.getLevel(timeIdx), x, x());
#elif defined USE_CARTESIAN_DOMAIN
    x() = (*q.getLevel(timeIdx))()+(*H.getLevel(timeIdx))*y();
    domain.constrain(x);
#endif
}

void Parcel::getBodyCoord(const Domain &domain,
                          const TimeLevelIndex<2> &timeIdx,
                          const SpaceCoord &x, BodyCoord &y) {
#ifdef USE_SPHERE_DOMAIN
    domain.project(geomtk::SphereDomain::STEREOGRAPHIC,
                   *q.getLevel(timeIdx), x, y());
    y() = (*invH.getLevel(timeIdx))*y();
#elif defined USE_CARTESIAN_DOMAIN
    y() = (*invH.getLevel(timeIdx))*domain.diffCoord(x, (*q.getLevel(timeIdx)));
#endif
}
    
double Parcel::getShapeFunction(const TimeLevelIndex<2> &timeIdx,
                                const BodyCoord &y) {
    double f;
    ShapeFunction::evalFunc(y, f);
//    f /= detH.getLevel(timeIdx);
    return f;
}

void Parcel::updateShapeSize(const Domain &domain,
                             const TimeLevelIndex<2> &timeIdx) {
    BodyCoord y(domain.getNumDim());
    SpaceCoord x(domain.getNumDim());
    for (int m = 0; m < domain.getNumDim(); ++m) {
        y() = (*invH.getLevel(timeIdx))*(*H.getLevel(timeIdx))*V.col(m);
        getSpaceCoord(domain, timeIdx, y, x);
        shapeSize.getLevel(timeIdx)[m] = domain.calcDistance(x, *q.getLevel(timeIdx));
    }
}

}
