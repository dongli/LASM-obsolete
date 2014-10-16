#include "ShapeFunction.h"

namespace lasm {

double ShapeFunction::J;
vec ShapeFunction::nodes;
vec ShapeFunction::weights;
const Domain* ShapeFunction::domain;
double ShapeFunction::_maxValue;

void ShapeFunction::init(const Domain &domain_) {
    J = 1.0/12.0;
    nodes.resize(5);
    weights.resize(5);

    nodes(0) = -2.0/3.0; weights(0) =  41.0/1280.0;
    nodes(1) = -1.0/3.0; weights(1) = 316.0/1280.0;
    nodes(2) =      0.0; weights(2) = 566.0/1280.0;
    nodes(3) =  1.0/3.0; weights(3) = 316.0/1280.0;
    nodes(4) =  2.0/3.0; weights(4) =  41.0/1280.0;
    
    _maxValue = pow(4.0/3.0, domain_.numDim());

    domain = &domain_;
}

void ShapeFunction::evalFunc(const BodyCoord &y, double &f) {
    f = 1.0;
    for (int i = 0; i < domain->numDim(); ++i) {
        if (-1.0 <= y(i) && y(i) <= -0.5) {
            f *= 2.0*pow(1.0+y(i), 3.0);
        } else if (-0.5 <= y(i) && y(i) <= 0.0) {
            f *= 1.0-6.0*pow(y(i), 2.0)*(1.0+y(i));
        } else if (0.0 <= y(i) && y(i) <= 0.5) {
            f *= 1.0-6.0*pow(y(i), 2.0)*(1.0-y(i));
        } else if (0.5 <= y(i) && y(i) <= 1.0) {
            f *= 2.0*pow(1.0-y(i), 3.0);
        } else {
            f = 0.0;
            return;
        }
        f *= 4.0/3.0;
    }
}

void ShapeFunction::evalDerv(const BodyCoord& y, vec &d) {
#ifndef NDEBUG
    assert(d.size() == domain->numDim());
#endif
    d.ones();
    for (int i = 0; i < domain->numDim(); ++i) {
        for (int j = 0; j < domain->numDim(); ++j) {
            if (i == j) {
                if (-1.0 <= y(j) && y(j) <= -0.5) {
                    d(i) *= 6.0*pow(1.0+y(j), 2.0);
                } else if (-0.5 <= y(j) && y(j) <= 0.0) {
                    d(i) *= -12.0*y(j)-18.0*pow(y(j), 2.0);
                } else if (0.0 <= y(j) && y(j) <= 0.5) {
                    d(i) *= -12.0*y(j)+18.0*pow(y(j), 2.0);
                } else if (0.5 <= y(j) && y(j) <= 1.0) {
                    d(i) *= -6.0*pow(1.0-y(j), 2.0);
                } else {
                    d(i) = 0.0;
                    continue;
                }
            } else {
                if (-1.0 <= y(j) && y(j) <= -0.5) {
                    d(i) *= 2.0*pow(1.0+y(j), 3.0);
                } else if (-0.5 <= y(j) && y(j) <= 0.0) {
                    d(i) *= 1.0-6.0*pow(y(j), 2.0)*(1.0+y(j));
                } else if (0.0 <= y(j) && y(j) <= 0.5) {
                    d(i) *= 1.0-6.0*pow(y(j), 2.0)*(1.0-y(j));
                } else if (0.5 <= y(j) && y(j) <= 1.0) {
                    d(i) *= 2.0*pow(1.0-y(j), 3.0);
                } else {
                    d(i) = 0.0;
                    continue;
                }
            }
            d(i) *= 4.0/3.0;
        }
    }
}

}
