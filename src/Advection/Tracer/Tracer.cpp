#include "Tracer.h"
#include "TracerSkeleton.h"
#include "TracerMeshCell.h"

namespace lasm {

Tracer::Tracer(int numDim) : Parcel(numDim) {
    skeleton = new TracerSkeleton(this, numDim);
    type = GOOD_SHAPE;
    numConnectedCell = 0;
    father = NULL;
}

Tracer::~Tracer() {
    delete skeleton;
}

Tracer& Tracer::operator=(const Tracer &other) {
    Parcel::operator=(other);
    if (this != &other) {
        for (int s = 0; s < density.size(); ++s) {
            density[s] = other.density[s];
        }
        *skeleton = *(other.skeleton);
        // TODO: Handle connected cells.
    }
    return *this;
}

void Tracer::resetConnectedCells() {
    numConnectedCell = 0;
}

void Tracer::connect(TracerMeshCell *cell, double weight) {
    if (numConnectedCell == connectedCells.size()) {
        connectedCells.push_back(cell);
    } else {
        connectedCells[numConnectedCell] = cell;
    }
    numConnectedCell++;
}
    
void Tracer::updateDeformMatrix(const Domain &domain,
                                const Mesh &mesh,
                                const TimeLevelIndex<2> &timeIdx) {
    skeleton->updateLocalCoord(domain, timeIdx);
    const vector<vec> &xl = skeleton->getLocalCoords(timeIdx);
    const vector<BodyCoord*> &y = skeleton->getBodyCoords();
    mat &H0 = *H.getLevel(timeIdx);
    if (domain.getNumDim() == 2) {
        // elements of four matrices H_12, H_14, H_32, H_34
        double h11_1 = xl[0][0]/(*y[0])(0);
        double h21_1 = xl[0][1]/(*y[0])(0);
        double h11_3 = xl[2][0]/(*y[2])(0);
        double h21_3 = xl[2][1]/(*y[2])(0);
        double h12_2 = xl[1][0]/(*y[1])(1);
        double h22_2 = xl[1][1]/(*y[1])(1);
        double h12_4 = xl[3][0]/(*y[3])(1);
        double h22_4 = xl[3][1]/(*y[3])(1);
        // the final matrix is the combination of the four
        H0(0, 0) = (h11_1+h11_3)*0.5;
        H0(0, 1) = (h12_2+h12_4)*0.5;
        H0(1, 0) = (h21_1+h21_3)*0.5;
        H0(1, 1) = (h22_2+h22_4)*0.5;
        if (!svd(U, S, V, H0)) {
            REPORT_ERROR("Encounter error with arma::svd!");
        }
#ifndef NDEBUG
        assert(S[0] >= S[1]);
#endif
        // reset the determinant of H to parcel volume
        // NOTE: Parcel volume is represented by the first tracer.
        if (density.size() != 0) {
            int s;
            for (s = 0; s < mass.size(); ++s) {
                if (mass[s] != 0) {
                    break;
                }
            }
            if (s != mass.size()) {
                double volume = mass[s]/density[s];
                double tmp = sqrt(volume/(S[0]*S[1]));
                S[0] *= tmp; S[1] *= tmp;
                H0 = U*diagmat(S)*V.t();
                resetSkeleton(domain, mesh, timeIdx);
            }
        }
        detH.getLevel(timeIdx) = arma::prod(S);
        *invH.getLevel(timeIdx) = inv(H0);
    } else if (domain.getNumDim() == 3) {
        REPORT_ERROR("Under construction!");
    }
    updateShapeSize(domain, timeIdx);
}

void Tracer::resetDeformMatrix(const Domain &domain,
                               const Mesh &mesh,
                               const TimeLevelIndex<2> &timeIdx,
                               const SpaceCoord &x, const vec &S) {
    vec xs1(domain.getNumDim()), xs2(domain.getNumDim());
    mat &H0 = *H.getLevel(timeIdx);
#ifdef LASM_USE_SPHERE_DOMAIN
    domain.project(geomtk::SphereDomain::STEREOGRAPHIC,
                   *q.getLevel(timeIdx), x, xs1);
    // x should be one vertex of the major axis
    xs1 *= S[0]/norm(xs1);
    if (domain.getNumDim() == 2) {
        // xs1 -> (1, 0)  xs2 -> (0, 1)
        xs2[0] = S[1]/S[0];
        xs2[1] = -xs1[0]*xs2[0]/xs1[1];
        xs2 *= S[1]/norm(xs2);
        double cosTheta = xs1[0]/norm(xs1);
        double sinTheta = xs1[1]/norm(xs1);
        if (-sinTheta*xs2[0]+cosTheta*xs2[1] < 0) {
            xs2 *= -1;
        }
        H0(0, 0) = xs1[0]; H0(0, 1) = xs2[0];
        H0(1, 0) = xs1[1]; H0(1, 1) = xs2[1];
        if (!svd(U, this->S, V, H0)) {
            REPORT_ERROR("Encounter error with arma::svd!");
        }
        detH.getLevel(timeIdx) = arma::prod(S);
        *invH.getLevel(timeIdx) = inv(H0);
        resetSkeleton(domain, mesh, timeIdx);
    } else if (domain.getNumDim() == 3) {
        REPORT_ERROR("Under construction!");
    }
#else
    REPORT_ERROR("Under construction!");
#endif
    updateShapeSize(domain, timeIdx);
}

void Tracer::resetSkeleton(const Domain &domain, const Mesh &mesh,
                           const TimeLevelIndex<2> &timeIdx) {
    // reset skeleton points
    const vector<BodyCoord*> &y = skeleton->getBodyCoords();
    vector<SpaceCoord*> &x = skeleton->getSpaceCoords(timeIdx);
    vector<MeshIndex*> &meshIdx = skeleton->getMeshIdxs(timeIdx);
    for (int i = 0; i < x.size(); ++i) {
        getSpaceCoord(domain, timeIdx, *y[i], *x[i]);
        meshIdx[i]->reset();
        meshIdx[i]->locate(mesh, *x[i]);
    }
}

void Tracer::dump(const TimeLevelIndex<2> &timeIdx, const Domain &domain,
                  std::ofstream &file, int idx) {
    // -------------------------------------------------------------------------
    // tracer centroid
    file << "centroid" << idx << " = (/" << getX(timeIdx)(0)/RAD << "," <<
        getX(timeIdx)(1)/RAD << "/)" << endl;
    // -------------------------------------------------------------------------
    // tracer shape
    int n = 100;
    file << "shape" << idx << " = new((/" << n << ",2/), double)" << endl;
    double dtheta = PI2/(n-1);
    BodyCoord y(2);
    SpaceCoord x(2);
    vector<vec::fixed<2> > shape(n);
    for (int i = 0; i < shape.size(); ++i) {
        double theta = i*dtheta;
        y(0) = cos(theta);
        y(1) = sin(theta);
        getSpaceCoord(domain, timeIdx, y, x);
        shape[i] = x()/RAD;
    }
    for (int m = 0; m < 2; ++m) {
        file << "shape" << idx << "(:," << m << ") = (/";
        for (int i = 0; i < shape.size(); ++i) {
            if (i != shape.size()-1) {
                file << shape[i][m] << ",";
            } else {
                file << shape[i][m] << "/)" << endl;
            }
        }
    }
    // -------------------------------------------------------------------------
    // connected cells
    if (numConnectedCell != 0) {
        file << "ngb_cells" << idx << " = new((/" << numConnectedCell
            << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_cells" << idx << "(:," << m << ") = (/";
            for (int i = 0; i < numConnectedCell; ++i) {
                TracerMeshCell *cell = connectedCells[i];
                if (i != numConnectedCell-1) {
                    file << cell->getCoord()(m) << ",";
                } else {
                    file << cell->getCoord()(m) << "/)" << endl;
                }
            }
        }
    }
    // -------------------------------------------------------------------------
    // neighbor tracer locations
    int numNeighborTracer = 0;
    for (int i = 0; i < numConnectedCell; ++i) {
        numNeighborTracer += connectedCells[i]->getNumContainedTracer();
    }
    if (numNeighborTracer != 0) {
        file << "ngb_tracers" << idx << " = new((/" << numNeighborTracer <<
            ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_tracers" << idx << "(:," << m << ") = (/";
            int k = 0;
            for (int i = 0; i < numConnectedCell; ++i) {
                TracerMeshCell *cell = connectedCells[i];
                vector<Tracer*> &tracers = cell->getContainedTracers();
                for (int j = 0; j < cell->getNumContainedTracer(); ++j) {
                    if (k != numNeighborTracer-1) {
                        file << tracers[j]->getX(timeIdx)(m) << ",";
                    } else {
                        file << tracers[j]->getX(timeIdx)(m) << "/)" << endl;
                    }
                    k++;
                }
            }
        }
    }
}

void Tracer::dump(const TimeLevelIndex<2> &timeIdx, const Domain &domain) {
    std::ofstream file; file.open("tracer_dump.txt");
    // -------------------------------------------------------------------------
    // centroid location
    file << "centroid = (/" << getX(timeIdx)(0) << "," << getX(timeIdx)(1) << "/)" << endl;
    // -------------------------------------------------------------------------
    // neighbor cell locations
    if (numConnectedCell != 0) {
        file << "ngb_cells = new((/" << numConnectedCell << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_cells(:," << m << ") = (/";
            for (int i = 0; i < numConnectedCell; ++i) {
                TracerMeshCell *cell = connectedCells[i];
                if (i != numConnectedCell-1) {
                    file << cell->getCoord()(m) << ",";
                } else {
                    file << cell->getCoord()(m) << "/)" << endl;
                }
            }
        }
    }
    // -------------------------------------------------------------------------
    // neighbor tracer locations
    int numNeighborTracer = 0;
    for (int i = 0; i < numConnectedCell; ++i) {
        numNeighborTracer += connectedCells[i]->getNumContainedTracer();
    }
    if (numNeighborTracer != 0) {
        file << "ngb_tracers = new((/" << numNeighborTracer << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_tracers(:," << m << ") = (/";
            int k = 0;
            for (int i = 0; i < numConnectedCell; ++i) {
                TracerMeshCell *cell = connectedCells[i];
                vector<Tracer*> &tracers = cell->getContainedTracers();
                for (int j = 0; j < cell->getNumContainedTracer(); ++j) {
                    if (k != numNeighborTracer-1) {
                        file << tracers[j]->getX(timeIdx)(m) << ",";
                    } else {
                        file << tracers[j]->getX(timeIdx)(m) << "/)" << endl;
                    }
                    k++;
                }
            }
        }
    }
    // -------------------------------------------------------------------------
    // tracer shape
    int n = 100;
    file << "shape = new((/" << n << ",2/), double)" << endl;
    double dtheta = PI2/(n-1);
    BodyCoord y(2);
    SpaceCoord x(2);
    vector<vec::fixed<2> > shape(n);
    for (int i = 0; i < shape.size(); ++i) {
        double theta = i*dtheta;
        y(0) = cos(theta);
        y(1) = sin(theta);
        getSpaceCoord(domain, timeIdx, y, x);
        shape[i] = x();
    }
    for (int m = 0; m < 2; ++m) {
        file << "shape(:," << m << ") = (/";
        for (int i = 0; i < shape.size(); ++i) {
            if (i != shape.size()-1) {
                file << shape[i][m] << ",";
            } else {
                file << shape[i][m] << "/)" << endl;
            }
        }
    }
    file.close();
}

}
