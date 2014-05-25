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

namespace lasm {

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::string;
using std::vector;
using std::list;
using std::find;

using arma::vec;
using arma::mat;
using arma::cube;
using arma::sp_mat;
using arma::svd;

using geomtk::RAD;
using geomtk::PI2;
using geomtk::BILINEAR;
using geomtk::Time;
using geomtk::TimeLevels;
using geomtk::TimeLevelIndex;
using geomtk::SystemTools;
using geomtk::TimeUnit;

const int FULL = geomtk::RLLStagger::GridType::FULL;
const int HALF = geomtk::RLLStagger::GridType::HALF;
const int CENTER = geomtk::RLLStagger::Location::CENTER;

// shortcuts for MLPACK classes
typedef mlpack::tree::BinarySpaceTree<
mlpack::bound::HRectBound<2>,
mlpack::range::RangeSearchStat> Tree;
typedef mlpack::metric::EuclideanDistance Metric;
typedef mlpack::range::RangeSearch<Metric, Tree> Searcher;

typedef geomtk::SphereDomain Domain;
typedef geomtk::SphereCoord SpaceCoord;
typedef geomtk::BodyCoord BodyCoord;
typedef geomtk::SphereVelocity Velocity;
typedef geomtk::RLLMesh Mesh;
typedef geomtk::RLLMeshIndex MeshIndex;
template <typename T, int N = 2>
using Field = geomtk::RLLField<T, N>;
typedef geomtk::NumericRLLField<double, 2> ScalarField;
typedef geomtk::NumericRLLField<double, 1> SingleScalarField;
typedef geomtk::RLLVelocityField VelocityField;
typedef geomtk::RLLRegrid Regrid;

#define LASM_USE_SPHERE_DOMAIN

}

#endif
