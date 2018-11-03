#ifndef FpmcInterface_Fpmc_h
#define FpmcInterface_Fpmc_h

/** \class Fpmc
 *
 * Generates Fpmc HepMC events
 *
 ***************************************/

#include "Fpmc/Fpmc.h"
#include "HepMC/Version.h"

#ifndef HEPMC_VERSION_CODE // HepMC v2
# define HEPMC_VERSION2
# include "HepMC/IO_HERWIG.h"
#endif

namespace HepMC { class GenEvent; }
namespace fpmc
{
  class HepMCWrapper : public Fpmc
  {
    public:
      HepMCWrapper( const char* );
      ~HepMCWrapper();

      /// Retrieve the last event generated
      const HepMC::GenEvent& event();
      /// Write the last event generated onto a file
      void write( const char* );

    private:
#ifdef HEPMC_VERSION2
      HepMC::IO_HERWIG conv_;
#endif
      /// Last event generated
      std::unique_ptr<HepMC::GenEvent> hepMCEvt_;

      /// HepMC verbosity
      bool hepMCVerbosity_;
  };
}

#endif
