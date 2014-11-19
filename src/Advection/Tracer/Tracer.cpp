#include "Tracer.h"
#include "TracerSkeleton.h"

namespace lasm {

Tracer::Tracer(int numDim) : Parcel(numDim) {
    _linearDegeneration = -999.0;
    _skeleton = new TracerSkeleton(this, numDim);
    type = GOOD_SHAPE;
    vy = new BodyCoord(numDim);
    vx = new SpaceCoord(numDim);
    _numConnectedCell = 0;
}

Tracer::~Tracer() {
    delete _skeleton;
    delete vy;
    delete vx;
}

Tracer& Tracer::operator=(const Tracer &other) {
    Parcel::operator=(other);
    if (this != &other) {
        for (int s = 0; s < _density.size(); ++s) {
            _density[s] = other._density[s];
            _mass[s] = other._mass[s];
        }
        *_skeleton = *(other._skeleton);
        // TODO: Handle connected cells.
    }
    return *this;
}

void Tracer::resetConnectedCells() {
    _numConnectedCell = 0;
}

void Tracer::connectCell(int i) {
    if (_numConnectedCell == _connectedCellIndices.size()) {
        _connectedCellIndices.push_back(i);
    } else {
        _connectedCellIndices[_numConnectedCell] = i;
    }
    _numConnectedCell++;
}
    
void Tracer::updateDeformMatrix(const Domain &domain,
                                const Mesh &mesh,
                                const TimeLevelIndex<2> &timeIdx) {
    _skeleton->updateLocalCoord(domain, timeIdx);
    const vector<vec> &xl = _skeleton->localCoords(timeIdx);
    const vector<BodyCoord*> &y = _skeleton->bodyCoords();
    mat &H0 = *_H.level(timeIdx);
    if (domain.numDim() == 2) {
        // Calculate the elements of four matrices H1, H2, H3, H4.
        double h11_1 = xl[0][0]/(*y[0])(0);
        double h21_1 = xl[0][1]/(*y[0])(0);
        double h11_3 = xl[2][0]/(*y[2])(0);
        double h21_3 = xl[2][1]/(*y[2])(0);
        double h12_2 = xl[1][0]/(*y[1])(1);
        double h22_2 = xl[1][1]/(*y[1])(1);
        double h12_4 = xl[3][0]/(*y[3])(1);
        double h22_4 = xl[3][1]/(*y[3])(1);
        // The final matrix is the combination of the four.
        H0(0, 0) = (h11_1+h11_3)*0.5;
        H0(0, 1) = (h12_2+h12_4)*0.5;
        H0(1, 0) = (h21_1+h21_3)*0.5;
        H0(1, 1) = (h22_2+h22_4)*0.5;
        if (!svd(_U, _S, _V, H0)) {
            REPORT_ERROR("Failed to do SVD on a matrix!");
        }
#ifndef NDEBUG
        assert(_S[0] >= _S[1]);
#endif
        // reset the determinant of H to parcel volume
        // NOTE: Parcel volume is represented by the first tracer.
        if (_density.size() != 0) {
            int s;
            for (s = 0; s < _mass.size(); ++s) {
                if (_mass[s] != 0) {
                    break;
                }
            }
            if (s != _mass.size()) {
                double volume = _mass[s]/_density[s];
                double tmp = sqrt(volume/(_S[0]*_S[1]));
                _S[0] *= tmp; _S[1] *= tmp;
                H0 = _U*diagmat(_S)*_V.t();
            }
        }
        _detH.level(timeIdx) = arma::prod(_S);
        *_invH.level(timeIdx) = inv(H0);
        (*vy)() = *_invH.level(timeIdx)*H0*_V.col(0);
        calcSpaceCoord(domain, timeIdx, *vy, *vx);
        _filament = _S[0]/_S[1];
    } else if (domain.numDim() == 3) {
        REPORT_ERROR("Under construction!");
    }
    updateShapeSize(domain, timeIdx);
}

void Tracer::updateDeformMatrix(const Domain &domain,
                                const TimeLevelIndex<2> &timeIdx) {
    mat &H0 = *_H.level(timeIdx);
    H0 = _U*diagmat(_S)*_V.t();
    _detH.level(timeIdx) = arma::prod(_S);
    *_invH.level(timeIdx) = inv(H0);
    (*vy)() = *_invH.level(timeIdx)*H0*_V.col(0);
    calcSpaceCoord(domain, timeIdx, *vy, *vx);
    _filament = _S[0]/_S[1];
}

void Tracer::resetSkeleton(const Domain &domain, const Mesh &mesh,
                           const TimeLevelIndex<2> &timeIdx) {
    // reset skeleton points
    const vector<BodyCoord*> &y = _skeleton->bodyCoords();
    vector<SpaceCoord*> &x = _skeleton->spaceCoords(timeIdx);
    vector<MeshIndex*> &meshIdx = _skeleton->meshIndices(timeIdx);
    for (int i = 0; i < x.size(); ++i) {
        calcSpaceCoord(domain, timeIdx, *y[i], *x[i]);
        meshIdx[i]->reset();
        meshIdx[i]->locate(mesh, *x[i]);
    }
}

void Tracer::dump(const TimeLevelIndex<2> &timeIdx, const Domain &domain,
                  const MeshAdaptor &meshAdaptor, ofstream &file, int idx) {
    // tracer centroid
    file << "centroid" << idx << " = (/" << x(timeIdx)(0)/RAD << ",";
    file << x(timeIdx)(1)/RAD << "/)" << endl;
    // tracer skeleton points
    file << "skel_points" << idx << " = new((/4,2/), double)" << endl;
    for (int m = 0; m < 2; ++m) {
        file << "skel_points" << idx << "(:," << m << ") = (/";
        for (int i = 0; i < _skeleton->spaceCoords(timeIdx).size(); ++i) {
            file << (*_skeleton->spaceCoords(timeIdx)[i])(m)/RAD;
            if (i != 3) {
                file << ",";
            } else {
                file << "/)" << endl;
            }
        }
    }
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
        calcSpaceCoord(domain, timeIdx, y, x);
        shape[i] = x();
    }
    for (int m = 0; m < 2; ++m) {
        file << "shape" << idx << "(:," << m << ") = (/";
        for (int i = 0; i < shape.size(); ++i) {
            if (i != shape.size()-1) {
                file << shape[i][m]/RAD << ",";
            } else {
                file << shape[i][m]/RAD << "/)" << endl;
            }
        }
    }
    // connected cells
    if (_numConnectedCell != 0) {
        file << "ngb_cells" << idx << " = new((/" << _numConnectedCell
            << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_cells" << idx << "(:," << m << ") = (/";
            for (int i = 0; i < _numConnectedCell; ++i) {
                file << meshAdaptor.coord(_connectedCellIndices[i])(m)/RAD;
                if (i != _numConnectedCell-1) {
                    file << ",";
                } else {
                    file << "/)" << endl;
                }
            }
        }
    }
    // neighbor tracer locations
    int numNeighborTracer = 0;
    for (int i = 0; i < _numConnectedCell; ++i) {
        int j = _connectedCellIndices[i];
        numNeighborTracer += meshAdaptor.numContainedTracer(j);
    }
    if (numNeighborTracer != 0) {
        file << "ngb_tracers" << idx << " = new((/" << numNeighborTracer << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_tracers" << idx << "(:," << m << ") = (/";
            int k = 0;
            for (int i = 0; i < _numConnectedCell; ++i) {
                const vector<Tracer*> &tracers = meshAdaptor.containedTracers(_connectedCellIndices[i]);
                for (int j = 0; j < meshAdaptor.numContainedTracer(_connectedCellIndices[i]); ++j) {
                    file << tracers[j]->x(timeIdx)(m)/RAD;
                    if (k != numNeighborTracer-1) {
                        file << ",";
                    } else {
                        file << "/)" << endl;
                    }
                    k++;
                }
            }
        }
    }
}

void Tracer::dump(const TimeLevelIndex<2> &timeIdx, const Domain &domain,
                  const MeshAdaptor &meshAdaptor) {
    std::ofstream file; file.open("tracer_dump.txt");
    // tracer centroid
    file << "centroid = (/" << x(timeIdx)(0)/RAD << ",";
    file << x(timeIdx)(1)/RAD << "/)" << endl;
    // tracer skeleton points
    file << "skel_points = new((/4,2/), double)" << endl;
    for (int m = 0; m < 2; ++m) {
        file << "skel_points(:," << m << ") = (/";
        for (int i = 0; i < _skeleton->spaceCoords(timeIdx).size(); ++i) {
            file << (*_skeleton->spaceCoords(timeIdx)[i])(m)/RAD;
            if (i != 3) {
                file << ",";
            } else {
                file << "/)" << endl;
            }
        }
    }
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
        calcSpaceCoord(domain, timeIdx, y, x);
        shape[i] = x();
    }
    for (int m = 0; m < 2; ++m) {
        file << "shape(:," << m << ") = (/";
        for (int i = 0; i < shape.size(); ++i) {
            if (i != shape.size()-1) {
                file << shape[i][m]/RAD << ",";
            } else {
                file << shape[i][m]/RAD << "/)" << endl;
            }
        }
    }
    // neighbor cell locations
    if (_numConnectedCell != 0) {
        file << "ngb_cells = new((/" << _numConnectedCell << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_cells(:," << m << ") = (/";
            for (int i = 0; i < _numConnectedCell; ++i) {
                int j = _connectedCellIndices[i];
                file << meshAdaptor.coord(j)(m)/RAD;
                if (i != _numConnectedCell-1) {
                    file << ",";
                } else {
                    file << "/)" << endl;
                }
            }
        }
    }
    // neighbor tracer locations
    int numNeighborTracer = 0;
    for (int i = 0; i < _numConnectedCell; ++i) {
        int j = _connectedCellIndices[i];
        numNeighborTracer += meshAdaptor.numContainedTracer(j);
    }
    if (numNeighborTracer != 0) {
        file << "ngb_tracers = new((/" << numNeighborTracer << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_tracers(:," << m << ") = (/";
            int k = 0;
            for (int i = 0; i < _numConnectedCell; ++i) {
                const vector<Tracer*> &tracers = meshAdaptor.containedTracers(_connectedCellIndices[i]);
                for (int j = 0; j < meshAdaptor.numContainedTracer(_connectedCellIndices[i]); ++j) {
                    file << tracers[j]->x(timeIdx)(m)/RAD;
                    if (k != numNeighborTracer-1) {
                        file << ",";
                    } else {
                        file << "/)" << endl;
                    }
                    k++;
                }
            }
        }
    }
    file.close();
}

}
