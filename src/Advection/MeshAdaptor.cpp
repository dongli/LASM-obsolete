#include "MeshAdaptor.h"
#include "Tracer.h"
#include "TracerManager.h"

namespace lasm {
    
MeshAdaptor::MeshAdaptor() {
    domain = NULL;
    mesh = NULL;
    REPORT_ONLINE;
}
    
MeshAdaptor::~MeshAdaptor() {
    if (mesh != NULL) {
        for (int i = 0; i < _density.size(); ++i) {
            delete _density[i];
            delete _mass[i];
        }
    }
    REPORT_OFFLINE;
}

void MeshAdaptor::init(const Domain &domain, const Mesh &mesh) {
    this->domain = &domain;
    this->mesh = &mesh;
    _numConnectedTracer.resize(mesh.totalNumGrid(CENTER));
    _connectedTracers.resize(mesh.totalNumGrid(CENTER));
    _remapWeights.resize(mesh.totalNumGrid(CENTER));
    _numContainedTracer.resize(mesh.totalNumGrid(CENTER));
    _containedTracers.resize(mesh.totalNumGrid(CENTER));
}

void MeshAdaptor::input(const string &fileName,
                        const TracerManager &tracerManager) {
    TimeLevelIndex<2> timeIdx;
    int ncId, ret;
    double *data = new double[mesh->totalNumGrid(CENTER)];
    ret = nc_open(fileName.c_str(), NC_NOWRITE, &ncId);
    CHECK_NC_OPEN(ret, fileName);
    for (int s = 0; s < tracerManager.numSpecies(); ++s) {
        int varId;
        const string &varName = tracerManager.speciesInfo(s).name();
        ret = nc_inq_varid(ncId, varName.c_str(), &varId);
        CHECK_NC_INQ_VARID(ret, fileName, varName);
        ret = nc_get_var(ncId, varId, data);
        for (int i = 0; i < mesh->totalNumGrid(CENTER); ++i) {
            (*_density[s])(timeIdx, i) = data[i];
            (*_mass[s])(timeIdx, i) = data[i]*volume(i);
        }
    }
    delete [] data;
    ret = nc_close(ncId);
    CHECK_NC_CLOSE(ret, fileName);
}

void MeshAdaptor::registerTracer(const string &name, const string &units,
                                 const string &brief) {
    _density.push_back(new ScalarField);
    _density.back()->create(name, units, brief, *mesh, CENTER, domain->numDim());
    _mass.push_back(new ScalarField);
    _mass.back()->create(name, units, brief, *mesh, CENTER, domain->numDim());
}

void MeshAdaptor::resetSpecies(const TimeLevelIndex<2> &timeIdx) {
    for (int s = 0; s < _density.size(); ++s) {
        for (int i = 0; i < mesh->totalNumGrid(CENTER, domain->numDim()); ++i) {
            (*_density[s])(timeIdx, i) = 0;
            (*_mass[s])(timeIdx, i) = 0;
        }
    }
}
    
void MeshAdaptor::resetContainedTracers() {
    for (int i = 0; i < mesh->totalNumGrid(CENTER, domain->numDim()); ++i) {
        _numContainedTracer[i] = 0;
    }
}

void MeshAdaptor::resetConnectedTracers() {
    for (int i = 0; i < mesh->totalNumGrid(CENTER, domain->numDim()); ++i) {
        _numConnectedTracer[i] = 0;
    }
}
    
double MeshAdaptor::remapWeight(int i, Tracer *tracer) const {
    for (int j = 0; j < _numConnectedTracer[i]; ++j) {
        if (_connectedTracers[i][j] == tracer) {
            return _remapWeights[i][j];
        }
    }
    REPORT_ERROR("Tracer is not connected!");
}

void MeshAdaptor::containTracer(int i, Tracer *tracer) {
#ifndef NDEBUG
    for (int j = 0; j < _numContainedTracer[i]; ++j) {
        if (_containedTracers[i][j] == tracer) {
            REPORT_ERROR("Tracer (ID = " << tracer->ID() <<
                         ") has already been contained!");
        }
    }
#endif
    if (_numContainedTracer[i] == _containedTracers[i].size()) {
        _containedTracers[i].push_back(tracer);
    } else {
        _containedTracers[i][_numContainedTracer[i]] = tracer;
    }
    tracer->hostCellIndex() = i;
    _numContainedTracer[i]++;
}
    
void MeshAdaptor::connectTracer(int i, Tracer *tracer, double weight) {
#ifndef NDEBUG
    if (isTracerConnected(i, tracer)) {
        REPORT_ERROR("Tracer (ID = " << tracer->ID() <<
                     ") has already been connected!");
    }
#endif
    if (_numConnectedTracer[i] == _connectedTracers[i].size()) {
        _connectedTracers[i].push_back(tracer);
        _remapWeights[i].push_back(weight);
    } else {
        _connectedTracers[i][_numConnectedTracer[i]] = tracer;
        _remapWeights[i][_numConnectedTracer[i]] = weight;
    }
    _numConnectedTracer[i]++;
}
    
bool MeshAdaptor::isTracerConnected(int i, Tracer *tracer) {
    for (int j = 0; j < _numConnectedTracer[i]; ++j) {
        if (_connectedTracers[i][j] == tracer) {
            return true;
        }
    }
    return false;
}

} // lasm
