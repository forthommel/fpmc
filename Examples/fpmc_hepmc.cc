#include "Fpmc/HepMCWrapper.h"
#include "Fpmc/ArgsParser.h"

#include <iostream>

#ifdef HEPMC_VERSION2
#include "HepMC/IO_GenEvent.h"
#else
#include "HepMC/WriterAscii.h"
#include "HepMC/LHEF.h"
#include "HepMC/GenParticle.h" //FIXME
#endif

using namespace std;

int main( int argc, char* argv[] )
{
  fpmc::ArgsParser args( argc, argv, { "cfg", "nevents" }, { "fileout", "type" } );
  const auto required_params = args.required_parameters(), optional_params = args.optional_parameters();
  //----------------------------

  // Required parameters
  string datacard = required_params.at( "cfg" );
  unsigned int maxEvents = stoi( required_params.at( "nevents" ) );

  // Optional parameters
  const string outputFileName = ( optional_params.count( "fileout" ) == 0 )
    ? "fpmc.hepmc"
    : optional_params.at( "fileout" );
  const string outputType = ( optional_params.count( "type" ) == 0 )
    ? "hepmc"
    : optional_params.at( "type" );

  cout << "=========================================================" << endl
       << "FPMC (Wrapper) will initialize with parameters: " << endl
       << "  Datacard:    " << datacard << endl
       << "  N events:    " << maxEvents << endl
       << "  Output file: " << outputFileName << endl
       << "  Output type: " << outputType << endl
       << "=========================================================" << endl;

  fpmc::HepMCWrapper generator( datacard.c_str() );

  if ( outputType == "hepmc" ) {
#ifdef HEPMC_VERSION2
    std::ofstream output_file( outputFileName, ios::out );
    HepMC::IO_GenEvent output( output_file );
#else
    HepMC::WriterAscii output( outputFileName );
#endif
    for ( unsigned int evt = 0; evt < maxEvents; ++evt ) {
      cout << "[FPMC Wrapper] Processing event " << ( evt+1 ) << endl;
#ifdef HEPMC_VERSION2
      output.write_event( &generator.event() );
#else
      output.write_event( generator.event() );
#endif
    }
  }
  else if ( outputType == "lhef" ) {
#ifdef HEPMC_VERSION2
    throw runtime_error( "LHEF only works for HepMC v3+!" );
#else
    LHEF::Writer output( outputFileName );
    for ( unsigned int evt = 0; evt < maxEvents; ++evt ) {
      const auto& ev = generator.event();
      cout << ev.particles().size() << endl;
      for ( const auto& parts : ev.particles() ) {
        cout << parts->pid() << endl;
      }
    }
#endif
  }
  else
    throw runtime_error( "Unrecognised output mode: \""+outputType+"\"!" );
  //cout << "Cross section: " << generator.crossSection() << " +/- " << generator.crossSection() / sqrt( maxEvents ) << " pb"<< endl;
  return 0;
}
