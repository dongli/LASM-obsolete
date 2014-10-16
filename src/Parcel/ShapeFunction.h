#ifndef __LASM_ShapeFunction__
#define __LASM_ShapeFunction__

#include "lasm_commons.h"

namespace lasm {

class ShapeFunction {
public:
    static double J;
    // one-dimensional quadrature nodes and weights
    static vec nodes;
    static vec weights;
private:
    static const Domain *domain;
    static double _maxValue;
public:
    static void init(const Domain &domain);
    static double maxValue() { return _maxValue; }
    static void evalFunc(const BodyCoord& y, double &f);
    static void evalDerv(const BodyCoord& y, vec &d);

private:
    ShapeFunction() {}
    virtual ~ShapeFunction() {}
};

}

#endif
