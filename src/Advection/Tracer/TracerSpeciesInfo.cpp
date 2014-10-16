#include "TracerSpeciesInfo.h"

namespace lasm {

TracerSpeciesInfo::TracerSpeciesInfo() {
    _name = "N/A";
    _units = "N/A";
    _brief = "N/A";
}

TracerSpeciesInfo::TracerSpeciesInfo(const string &name, const string &units,
                                     const string &brief) {
    _name = name;
    _units = units;
    _brief = brief;
}

TracerSpeciesInfo::~TracerSpeciesInfo() {
}

}
