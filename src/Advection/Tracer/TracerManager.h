#ifndef __LASM_TracerManager__
#define __LASM_TracerManager__

#include "Tracer.h"
#include "TracerSpeciesInfo.h"

namespace lasm {

/**
 *  This class stores the tracer objects and manages the initialization,
 *  registration and output.
 */
class TracerManager {
    friend class AdvectionManager;
protected:
    const Domain *domain;
    vector<Tracer*> tracers;
    vector<TracerSpeciesInfo*> speciesInfos;
    double scale0; //>! initial parcel size scale (relative to grid cell)
public:
    TracerManager();
    virtual ~TracerManager();

    /**
     *  Initialize tracer manager.
     *
     *  @param domain     the space domain.
     *  @param mesh       the mesh where flow is defined.
     *  @param configManager the configuration manager.
     */
    void init(const Domain &domain, const Mesh &mesh,
              const geomtk::ConfigManager &configManager);

    void registerTracer(const string &name, const string &units,
                        const string &brief);

    /**
     *  Return the tracer species index based on the given name.
     *
     *  @param name the tracer species name.
     *
     *  @return The species index.
     */
    int getSpeciesIndex(const string &name) const;

    /**
     *  Get the number of tracer species.
     *
     *  @return The species number.
     */
    int getNumSpecies() const;
    
    const TracerSpeciesInfo& getSpeciesInfo(int speciesIdx) const;

    /**
     *  Output tracers on old time level into netCDF file.
     *
     *  @param fileName the output netCDF file name.
     *  @param timeIdx  the time level index.
     */
    void output(const string &fileName, const TimeLevelIndex<2> &timeIdx);
};

}

#endif // __LASM_TracerManager__
