#include "TracerSkeleton.h"

namespace lasm {

TracerSkeleton::TracerSkeleton(Tracer *host, int numDim) {
    this->host = host;
    y.resize(numDim*2);
    for (int l = 0; l < x.numLevel(); ++l) {
        x.level(l).resize(y.size());
        idx.level(l).resize(y.size());
        xl.level(l).resize(y.size());
        for (int i = 0; i < x.level(l).size(); ++i) {
            x.level(l)[i] = new SpaceCoord(numDim);
            idx.level(l)[i] = new MeshIndex(numDim);
            xl.level(l)[i].set_size(numDim);
        }
    }
    for (int i = 0; i < y.size(); ++i) {
        y[i] = new BodyCoord(numDim);
    }
    double d = 1;
    if (numDim == 2) {
        (*y[0])() <<  -d << 0.0 << arma::endr;
        (*y[1])() << 0.0 <<  -d << arma::endr;
        (*y[2])() <<   d << 0.0 << arma::endr;
        (*y[3])() << 0.0 <<   d << arma::endr;
    } else if (numDim == 3) {
        (*y[0])() <<  -d << 0.0 << 0.0 << arma::endr;
        (*y[1])() << 0.0 <<  -d << 0.0 << arma::endr;
        (*y[2])() <<   d << 0.0 << 0.0 << arma::endr;
        (*y[3])() << 0.0 <<   d << 0.0 << arma::endr;
        (*y[4])() << 0.0 << 0.0 <<  -d << arma::endr;
        (*y[5])() << 0.0 << 0.0 <<   d << arma::endr;
    }
}

TracerSkeleton::~TracerSkeleton() {
    for (int l = 0; l < x.numLevel(); ++l) {
        for (int i = 0; i < x.level(l).size(); ++i) {
            delete x.level(l)[i];
            delete idx.level(l)[i];
        }
    }
    for (int i = 0; i < y.size(); ++i) {
        delete y[i];
    }
}
    
TracerSkeleton& TracerSkeleton::operator=(const TracerSkeleton &other) {
    if (this != &other) {
        for (int l = 0; l < x.numLevel(); ++l) {
            for (int i = 0; i < x.level(l).size(); ++i) {
                *(x.level(l)[i]) = *(other.x.level(l)[i]);
                *(idx.level(l)[i]) = *(other.idx.level(l)[i]);
            }
        }
        for (int i = 0; i < y.size(); ++i) {
            *y[i] = *other.y[i];
        }
    }
    return *this;
}
    
/*
                  (0,1)
                    o
                    |
                    |
                    |
     (-1,0) o-------x------o (1,0)
                    |
                    |
                    |
                    o
                  (0,-1)
*/

void TracerSkeleton::init(const Domain &domain, const Mesh &mesh,
                          double size) {
    // set the body and initial spatial coordinates of skeleton points
    TimeLevelIndex<2> initTimeIdx;
    const SpaceCoord &x0 = host->x(initTimeIdx);
    double dtheta = PI2/4;
#if defined USE_SPHERE_DOMAIN
    double lon, lat = M_PI_2-size/domain.radius();
    SpaceCoord xr(domain.numDim());
    for (int i = 0; i < 4; ++i) {
        lon = i*dtheta;
        xr.setCoord(lon, lat);
        domain.rotateBack(x0, *x.level(initTimeIdx)[i], xr);
        x.level(initTimeIdx)[i]->transformToCart(domain); // TODO: Do we need this?
    }
#elif defined USE_CARTESIAN_DOMAIN
    for (int i = 0; i < 4; ++i) {
        double theta = i*dtheta;
        (*x.level(initTimeIdx)[i])(0) = size*cos(theta)+x0(0);
        (*x.level(initTimeIdx)[i])(1) = size*sin(theta)+x0(1);
        domain.constrain(*x.level(initTimeIdx)[i]);
    }
#endif
    if (domain.numDim() == 3) {
        REPORT_ERROR("Under construction!");
        
    }
    for (int i = 0; i < y.size(); ++i) {
        idx.level(initTimeIdx)[i]->locate(mesh, *x.level(initTimeIdx)[i]);
    }
}
    
void TracerSkeleton::updateLocalCoord(const Domain &domain,
                                      const TimeLevelIndex<2> &timeIdx) {
    const SpaceCoord &x0 = host->x(timeIdx);
#ifdef USE_SPHERE_DOMAIN
    for (int i = 0; i < x.level(timeIdx).size(); ++i) {
        domain.project(geomtk::SphereDomain::STEREOGRAPHIC, x0,
                       *x.level(timeIdx)[i], xl.level(timeIdx)[i]);
    }
#else
    for (int i = 0; i < x.level(timeIdx).size(); ++i) {
        xl.level(timeIdx)[i] = domain.diffCoord((*x.level(timeIdx)[i]), x0);
    }
#endif
}
    
}
