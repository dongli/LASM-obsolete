#include "MeshAdaptor.h"
#include "Tracer.h"

namespace lasm {
    
MeshAdaptor::MeshAdaptor() {
    domain = NULL;
    mesh = NULL;
    REPORT_ONLINE;
}
    
MeshAdaptor::~MeshAdaptor() {
    if (mesh != NULL) {
        for (int i = 0; i < coords.size(); ++i) {
            delete coords[i];
        }
        for (int i = 0; i < densities.size(); ++i) {
            delete densities[i];
            delete masses[i];
        }
    }
    REPORT_OFFLINE;
}

void MeshAdaptor::init(const Domain &domain, const Mesh &mesh) {
    this->domain = &domain;
    this->mesh = &mesh;
    // TODO: 'coords' should be in geomtk::Mesh classes.
    coords.resize(mesh.getTotalNumGrid(CENTER, domain.getNumDim()));
    volumes.resize(mesh.getTotalNumGrid(CENTER, domain.getNumDim()));
    numConnectedTracer.resize(mesh.getTotalNumGrid(CENTER, domain.getNumDim()));
    connectedTracers.resize(mesh.getTotalNumGrid(CENTER, domain.getNumDim()));
    remapWeights.resize(mesh.getTotalNumGrid(CENTER, domain.getNumDim()));
    numContainedTracer.resize(mesh.getTotalNumGrid(CENTER, domain.getNumDim()));
    containedTracers.resize(mesh.getTotalNumGrid(CENTER, domain.getNumDim()));
    for (int i = 0; i < mesh.getTotalNumGrid(CENTER, domain.getNumDim()); ++i) {
        coords[i] = new SpaceCoord(domain.getNumDim());
        mesh.getGridCoord(i, CENTER, *coords[i]);
        coords[i]->transformToCart(domain);
        volumes[i] = mesh.getCellVolume(i);
    }
}

void MeshAdaptor::registerTracer(const string &name, const string &units,
                                 const string &brief) {
    densities.push_back(new ScalarField);
    densities.back()->create(name, units, brief, *mesh, CENTER, domain->getNumDim());
    masses.push_back(new ScalarField);
    masses.back()->create(name, units, brief, *mesh, CENTER, domain->getNumDim());
}

void MeshAdaptor::resetSpecies(const TimeLevelIndex<2> &timeIdx) {
    for (int s = 0; s < densities.size(); ++s) {
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER, domain->getNumDim()); ++i) {
            (*densities[s])(timeIdx, i) = 0;
            (*masses[s])(timeIdx, i) = 0;
        }
    }
}
    
void MeshAdaptor::resetContainedTracers() {
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, domain->getNumDim()); ++i) {
        numContainedTracer[i] = 0;
    }
}

void MeshAdaptor::resetConnectedTracers() {
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER, domain->getNumDim()); ++i) {
        numConnectedTracer[i] = 0;
    }
}
    
double MeshAdaptor::getRemapWeight(int i, Tracer *tracer) const {
    for (int j = 0; j < numConnectedTracer[i]; ++j) {
        if (connectedTracers[i][j] == tracer) {
            return remapWeights[i][j];
        }
    }
    REPORT_ERROR("Tracer is not connected!");
}

void MeshAdaptor::containTracer(int i, Tracer *tracer) {
#ifndef NDEBUG
    for (int j = 0; j < numContainedTracer[i]; ++j) {
        if (containedTracers[i][j] == tracer) {
            REPORT_ERROR("Tracer (ID = " << tracer->getID() <<
                         ") has already been contained!");
        }
    }
#endif
    if (numContainedTracer[i] == containedTracers[i].size()) {
        containedTracers[i].push_back(tracer);
    } else {
        containedTracers[i][numContainedTracer[i]] = tracer;
    }
    tracer->setHostCell(i);
    numContainedTracer[i]++;
}
    
void MeshAdaptor::connectTracer(int i, Tracer *tracer, double weight) {
#ifndef NDEBUG
    if (isTracerConnected(i, tracer)) {
        REPORT_ERROR("Tracer (ID = " << tracer->getID() <<
                     ") has already been connected!");
    }
#endif
    if (numConnectedTracer[i] == connectedTracers[i].size()) {
        connectedTracers[i].push_back(tracer);
        remapWeights[i].push_back(weight);
    } else {
        connectedTracers[i][numConnectedTracer[i]] = tracer;
        remapWeights[i][numConnectedTracer[i]] = weight;
    }
    numConnectedTracer[i]++;
}
    
bool MeshAdaptor::isTracerConnected(int i, Tracer *tracer) {
    for (int j = 0; j < numConnectedTracer[i]; ++j) {
        if (connectedTracers[i][j] == tracer) {
            return true;
        }
    }
    return false;
}

} // lasm
