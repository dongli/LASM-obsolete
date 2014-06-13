#include "TracerMeshCell.h"

namespace lasm {

TracerMeshCell::TracerMeshCell() {
    x = NULL;
    numConnectedTracer = 0;
    numContainedTracer = 0;
}

TracerMeshCell::~TracerMeshCell() {
    if (x != NULL) {
        delete x;
    }
}

void TracerMeshCell::resetConnectedTracers() {
    numConnectedTracer = 0;
}

void TracerMeshCell::connect(Tracer *tracer, double weight, double distance) {
#ifndef NDEBUG
    for (int i = 0; i < numConnectedTracer; ++i) {
        if (connectedTracers[i] == tracer) {
            REPORT_ERROR("Tracer (ID = " << tracer->getID() <<
                         ") has already been connected!");
        }
    }
#endif
    if (numConnectedTracer == connectedTracers.size()) {
        connectedTracers.push_back(tracer);
        remapWeights.push_back(weight);
        remapDistances.push_back(distance);
    } else {
        connectedTracers[numConnectedTracer] = tracer;
        remapWeights[numConnectedTracer] = weight;
        remapDistances[numConnectedTracer] = distance;
    }
    numConnectedTracer++;
}

void TracerMeshCell::disconnect(Tracer *tracer) {
#ifndef NDEBUG
    int i = 0;
    for (; i < numConnectedTracer; ++i) {
        if (connectedTracers[i] == tracer) {
            break;
        }
    }
    assert(i != numConnectedTracer);
#endif
    for (int i = 0; i < numConnectedTracer; ++i) {
        if (connectedTracers[i] == tracer) {
            for (int j = i+1; j < numConnectedTracer; ++j) {
                connectedTracers[j-1] = connectedTracers[j];
                remapWeights[j-1] = remapWeights[j];
            }
            numConnectedTracer--;
            break;
        }
    }
}

double TracerMeshCell::getRemapWeight(Tracer *tracer) const {
    for (int i = 0; i < numConnectedTracer; ++i) {
        if (connectedTracers[i] == tracer) {
            return remapWeights[i];
        }
    }
    REPORT_ERROR("Tracer is not connected!");
}

double TracerMeshCell::getRemapDistance(Tracer *tracer) const {
    for (int i = 0; i < numConnectedTracer; ++i) {
        if (connectedTracers[i] == tracer) {
            return remapDistances[i];
        }
    }
    REPORT_ERROR("Tracer is not connected!");
}

void TracerMeshCell::resetContainedTracers() {
    numContainedTracer = 0;
}

void TracerMeshCell::contain(Tracer *tracer) {
#ifndef NDEBUG
    for (int i = 0; i < numContainedTracer; ++i) {
        if (containedTracers[i] == tracer) {
            REPORT_ERROR("Tracer (ID = " << tracer->getID() <<
                         ") has already been contained!");
        }
    }
#endif
    if (numContainedTracer == containedTracers.size()) {
        containedTracers.push_back(tracer);
    } else {
        containedTracers[numContainedTracer] = tracer;
    }
    tracer->setHostCell(this);
    numContainedTracer++;
}

void TracerMeshCell::discontain(Tracer *tracer) {
#ifndef NDEBUG
    int i = 0;
    for (; i < numContainedTracer; ++i) {
        if (containedTracers[i] == tracer) {
            break;
        }
    }
    assert(i != numContainedTracer);
#endif
    for (int i = 0; i < numContainedTracer; ++i) {
        if (containedTracers[i] == tracer) {
            for (int j = i+1; j < numContainedTracer; ++j) {
                containedTracers[j-1] = containedTracers[j];
            }
            numContainedTracer--;
            break;
        }
    }
}

}
