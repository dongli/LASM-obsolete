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
    const Mesh *mesh;
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
                        const string &brief, bool smooth = false);

    /**
     *  Return the tracer species index based on the given name.
     *
     *  @param name the tracer species name.
     *
     *  @return The species index.
     */
    int speciesIndex(const string &name) const;

    /**
     *  Get the number of tracer species.
     *
     *  @return The species number.
     */
    int numSpecies() const;
    
    const TracerSpeciesInfo& speciesInfo(int speciesIdx) const;

    void resetSpecies();

    /**
     *  Input tracers on the old time level from a netCDF file.
     *
     *  @param fileName the input netCDF file name.
     */
    void input(const string &fileName);

    /**
     *  Output tracers on the old time level into a netCDF file.
     *
     *  @param timeIdx  the time level index.
     *  @param ncId     the output netCDF file ID.
     */
    void output(const TimeLevelIndex<2> &timeIdx, int ncId);
};

}

#endif // __LASM_TracerManager__
