#include "TracerSpeciesInfo.h"

namespace lasm {

TracerSpeciesInfo::TracerSpeciesInfo() {
    name = "N/A";
    units = "N/A";
    brief = "N/A";
}

TracerSpeciesInfo::TracerSpeciesInfo(const string &name, const string &units,
                                     const string &brief) {
    this->name = name;
    this->units = units;
    this->brief = brief;
}

TracerSpeciesInfo::~TracerSpeciesInfo() {
}

}
