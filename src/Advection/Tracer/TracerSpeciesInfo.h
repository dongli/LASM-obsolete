#ifndef __LASM_TracerSpeciesInfo__
#define __LASM_TracerSpeciesInfo__

#include "lasm_commons.h"

namespace lasm {

class TracerSpeciesInfo {
    string _name;
    string _units;
    string _brief;
    bool _smooth;
public:
    TracerSpeciesInfo();
    TracerSpeciesInfo(const string &name, const string &units,
                      const string &brief, bool smooth = false);
    ~TracerSpeciesInfo();

    const string &name() const { return _name; }

    const string &units() const { return _units; }

    const string &brief() const { return _brief; }

    bool smooth() const { return _smooth; }
};

}

#endif // __LASM_TracerSpeciesInfo__
