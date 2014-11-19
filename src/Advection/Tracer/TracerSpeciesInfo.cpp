#include "TracerSpeciesInfo.h"

namespace lasm {

TracerSpeciesInfo::TracerSpeciesInfo() {
    _name = "N/A";
    _units = "N/A";
    _brief = "N/A";
}

TracerSpeciesInfo::TracerSpeciesInfo(const string &name, const string &units,
                                     const string &brief, bool smooth) {
    _name = name;
    _units = units;
    _brief = brief;
    _smooth = smooth;
}

TracerSpeciesInfo::~TracerSpeciesInfo() {
}

}
