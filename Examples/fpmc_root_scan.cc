#include "Fpmc/Fpmc.h"
#include "Fpmc/ArgsParser.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TGraph2D.h>

#include <iostream>
#include <fstream>

using namespace std;

int main( int argc, char* argv[] )
{
  fpmc::Fpmc gen;
  gen.parameters().setSqrtS( 13.e3 );
  //gen.parameters().setProcessId( 16065 ); // spin-2 even resonance
  //gen.parameters().setProcessId( 16059 ); // anomalous gg>gg through QGC
  gen.parameters().setProcessId( 16016 ); // anomalous gg>gg through QGC
  gen.parameters().add( "aaanom", 2 );
  gen.parameters().setExotic( true );

  gen.parameters().setProcessType( fpmc::ProcessType::ExclusiveProcess );
  gen.parameters().setInteractionType( fpmc::InteractionType::QED );
  gen.parameters().setIntermediateFlux( fpmc::Flux::PhotonPhotonBudnevCoherent );
  gen.parameters().setPtRange( 75. );
  gen.parameters().setEtaRange( -2.5, 2.5 );
//  gen.parameters().setMassRange( 300., 2000. );

  gen.parameters().add( "iprint", 0 );
  gen.parameters().dump();

  unique_ptr<TFile> file( TFile::Open( "output.root", "recreate" ) );

  //gen.parameters().add( "aaa2", 0. ); //FIXME

  TGraph2D g_scan;
  g_scan.SetTitle( "#sigma(#zeta_{1},#zeta_{2}) (pb);#zeta_{1};#zeta_{2}" );
  std::ofstream out( "scan.dat" );

  const double min_z1 = -2.5e-12, max_z1 = 2.5e-12;
  const double min_z2 = -2.5e-12, max_z2 = 2.5e-12;
  const unsigned short num_z1 = 21, num_z2 = 21;
  for ( unsigned short i = 0; i < num_z1; ++i ) {
    const double z1 = min_z1+( max_z1-min_z1 )*i / ( num_z1-1 );
    for ( unsigned short j = 0; j < num_z2; ++j ) {
      const double z2 = min_z2+( max_z2-min_z2 )*j / ( num_z2-1 );

      gen.parameters().add( "a1a", z1 );
      gen.parameters().add( "a2a", z2 );
      hwpram_.IPRINT = 0;

      gen.initialise();
      const double xsec = gen.crossSection();
      g_scan.SetPoint( g_scan.GetN(), z1, z2, xsec );
      out << z1 << "\t" << z2 << "\t" << xsec << "\n";
    }
  }
  g_scan.Write();

  file->Write();
  out.close();

  return 0;
}

