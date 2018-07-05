#include "Fpmc/Fpmc.h"
#include "Fpmc/ArgsParser.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TGraph2D.h>

#include <iostream>

using namespace std;

int main( int argc, char* argv[] )
{
  fpmc::Fpmc gen;
  gen.parameters().dump();
  gen.parameters().setSqrtS( 13.e3 );
  gen.parameters().setProcessId( 16065 ); // spin-0 even resonance
  gen.parameters().setExotic( true );

  gen.parameters().setProcessType( fpmc::ProcessType::ExclusiveProcess );
  gen.parameters().setInteractionType( fpmc::InteractionType::QED );
  gen.parameters().setIntermediateFlux( fpmc::Flux::PhotonPhotonBudnevCoherent );
  gen.parameters().setPtRange( 75. );
  gen.parameters().setEtaRange( -2.5, 2.5 );
  gen.parameters().setMassRange( 300., 2000. );

  unique_ptr<TFile> file( TFile::Open( "output.root", "recreate" ) );

  gen.parameters().add( "aaa2", 0. ); //FIXME

  TGraph2D g_scan;

  const double min_lmass = 1., max_lmass = 4.;
  const double min_lfinv = -3., max_lfinv = 1.;
  const unsigned short num_mass = 25, num_finv = 25;
  for ( unsigned short i = 0; i < num_mass; ++i ) {
    const double lmass = min_lmass + ( max_lmass-min_lmass )*i / ( num_mass-1 ), mass = pow( 10., lmass );
    for ( unsigned short j = 0; j < num_finv; ++j ) {
      const double lfinv = min_lfinv + ( max_lfinv-min_lfinv )*j / ( num_finv-1 ), finv = pow( 10., lfinv );

      gen.parameters().add( "aam", mass );
      gen.parameters().add( "aaf0", 1./finv );
      gen.parameters().add( "aaw", pow( mass, 3 )*0.25*M_1_PI*finv*finv ); //FIXME

      gen.initialise();
      const double xsec = gen.crossSection();
      g_scan.SetPoint( g_scan.GetN(), lmass, lfinv, xsec );
    }
  }
  {
    TCanvas c;
    g_scan.Draw( "colz" );
    //c.SetLogx();
    //c.SetLogy();
    c.Write();
  }

  file->Write();

  return 0;
}

