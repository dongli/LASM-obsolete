#ifndef __LASM_TracerSpeciesInfo__
#define __LASM_TracerSpeciesInfo__

#include "lasm_commons.h"

namespace lasm {

class TracerSpeciesInfo {
    string name;
    string units;
    string brief;
public:
    TracerSpeciesInfo();
    TracerSpeciesInfo(const string &name, const string &units,
                      const string &brief);
    ~TracerSpeciesInfo();

    const string &getName() const { return name; }
    const string &getUnits() const { return units; }
    const string &getBrief() const { return brief; }
};

}

#endif // __LASM_TracerSpeciesInfo__
