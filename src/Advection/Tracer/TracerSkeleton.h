#ifndef __LASM_TracerSkeleton__
#define __LASM_TracerSkeleton__

#include "lasm_commons.h"
#include "Tracer.h"

namespace lasm {

class TracerSkeleton {
protected:
    Tracer *host;
    TimeLevels<vector<SpaceCoord*>, 2> x; //>! spatial coordinates
    vector<BodyCoord*> y;                 //>! fixed body coordinates
    TimeLevels<vector<MeshIndex*>, 2> idx;
    // Note: In sphere domain, xl are local stereographic projection coordinates,
    //       and in normal Cartesian domain, they should be the same with x.
    TimeLevels<vector<vec>, 2> xl;              //>! local coordinates
public:
    TracerSkeleton(Tracer *host, int numDim);
    virtual ~TracerSkeleton();

    TracerSkeleton& operator=(const TracerSkeleton &other);

    void init(const Domain &domain, const Mesh &mesh, double size);

    void updateLocalCoord(const Domain &domain,
                          const TimeLevelIndex<2> &timeIdx);

    vector<SpaceCoord*>& spaceCoords(const TimeLevelIndex<2> &timeIdx) {
        return x.level(timeIdx);
    }

    const vector<BodyCoord*>& bodyCoords() const { return y; }

    vector<MeshIndex*>& meshIndices(const TimeLevelIndex<2> &timeIdx) {
        return idx.level(timeIdx);
    }

    const vector<vec>& localCoords(const TimeLevelIndex<2> &timeIdx) const {
        return xl.level(timeIdx);
    }
};

}

#endif // __LASM_TracerSkeleton__
