#ifndef __lasm_commons__
#define __lasm_commons__

#include "geomtk.h"
#include <mlpack/methods/range_search/range_search.hpp>

#include <iostream>
#include <iomanip>
#include <list>
#include <string>
#include <ctime>
#include <fstream>

// Control macros.
//#define LASM_CARTESIAN_DOMAIN
#define LASM_SPHERE_DOMAIN
#define LASM_RLL_MESH
#define LASM_TENDENCY_ON_MESH
//#define LASM_TENDENCY_ON_PARCEL

namespace lasm {

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::string;
using std::vector;
using std::list;
using std::find;
using std::ofstream;

using arma::vec;
using arma::mat;
using arma::cube;
using arma::sp_mat;
using arma::svd;

// shortcuts for MLPACK classes
typedef mlpack::tree::BinarySpaceTree<
mlpack::bound::HRectBound<2>,
mlpack::range::RangeSearchStat> Tree;
typedef mlpack::metric::EuclideanDistance Metric;
typedef mlpack::range::RangeSearch<Metric, Tree> Searcher;

#if defined LASM_CARTESIAN_DOMAIN
// ############################################################
// Cartesian objects
typedef geomtk::CartesianDomain Domain;
typedef geomtk::SpaceCoord SpaceCoord;
typedef geomtk::Velocity Velocity;
typedef geomtk::CartesianMesh Mesh;
typedef geomtk::CartesianMeshIndex MeshIndex;
template <typename T, int N = 2>
using Field = geomtk::CartesianField<T, N>;
typedef geomtk::CartesianField<double, 2> ScalarField;
typedef geomtk::CartesianField<double, 1> SingleScalarField;
typedef geomtk::CartesianVelocityField VelocityField;
typedef geomtk::CartesianRegrid Regrid;
typedef geomtk::IOManager<geomtk::CartesianDataFile> IOManager;
#elif defined LASM_SPHERE_DOMAIN
// ############################################################
// Sphere objects
typedef geomtk::SphereDomain Domain;
#if defined LASM_RLL_MESH
// ============================================================
// Regular latitude-longitude mesh objects
typedef geomtk::SphereCoord SpaceCoord;
typedef geomtk::SphereVelocity Velocity;
typedef geomtk::RLLMesh Mesh;
typedef geomtk::RLLMeshIndex MeshIndex;
template <typename T, int N = 2>
using Field = geomtk::RLLField<T, N>;
typedef geomtk::RLLField<double, 2> ScalarField;
typedef geomtk::RLLField<double, 1> SingleScalarField;
typedef geomtk::RLLVelocityField VelocityField;
typedef geomtk::RLLRegrid Regrid;
typedef geomtk::IOManager<geomtk::RLLDataFile> IOManager;
typedef geomtk::RLLFilter<Mesh> Filter;
#elif defined USE_CUBIC_SPHERE_MESH
// ============================================================
// Cubic sphere mesh objects
// To be implemented.
#endif
#endif

typedef geomtk::BodyCoord BodyCoord;
typedef geomtk::TimeManager TimeManager;
typedef geomtk::ConfigManager ConfigManager;
typedef geomtk::StampString StampString;

#ifdef LASM_SPHERE_DOMAIN
using geomtk::RAD;
#else
const double RAD = 1;
#endif
using geomtk::PI2;
using geomtk::BILINEAR;
using geomtk::Time;
using geomtk::TimeLevels;
using geomtk::TimeLevelIndex;
using geomtk::SystemTools;
using geomtk::TimeUnit;
using geomtk::TimeStepUnit;

const int FULL = geomtk::RLLStagger::GridType::FULL;
const int HALF = geomtk::RLLStagger::GridType::HALF;
const int CENTER = geomtk::RLLStagger::Location::CENTER;
const int FULL_DIMENSION = geomtk::RLLSpaceDimensions::FULL_DIMENSION;

// Transition function.
static double transitionFunction(double x0, double y0,
                                 double x1, double y1,
                                 double x) {
    if (x < x0) {
        return y0;
    } else if (x > x1) {
        return y1;
    } else {
        double dx = x1-x0;
        double dy = y1-y0;
        double t = (x-x0)/dx;
        return dy*(4-3*t)*pow(t, 3)+y0;
    }
}

}

#endif
