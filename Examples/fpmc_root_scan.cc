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
  gen.parameters().setProcessId( 16065 ); // spin-0 even resonance
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
  ofstream of( "scan.dat" );

  gen.parameters().add( "aaa2", 0. ); //FIXME

  TGraph2D g_scan_linlin, g_scan_loglog;
  g_scan_linlin.SetName( "linlin" );
  g_scan_loglog.SetName( "loglog" );
  g_scan_linlin.SetTitle( ";m_{a} (GeV);f_{0}^{-1} (GeV^{-1});#sigma(m_{a},f_{0}^{-1})" );
  g_scan_loglog.SetTitle( ";log_{10}(m_{a}/GeV);log_{10}(f_{0}^{-1}/GeV^{-1});#sigma(m_{a},f_{0}^{-1})" );

  const double min_mass = 250., max_mass = 2000.; // in GeV
  const double min_finv = 1.e0, max_finv = 1.e-3; // in 1/GeV
  const double min_lmass = log10( min_mass ), max_lmass = log10( max_mass ), min_lfinv = log10( min_finv ), max_lfinv = log10( max_finv );
  const unsigned short num_mass = 25, num_finv = 25;
  for ( unsigned short i = 0; i < num_mass; ++i ) {
    const double lmass = min_lmass + ( max_lmass-min_lmass )*i / ( num_mass-1 ), mass = pow( 10., lmass );
    for ( unsigned short j = 0; j < num_finv; ++j ) {
      const double lfinv = min_lfinv + ( max_lfinv-min_lfinv )*j / ( num_finv-1 ), finv = pow( 10., lfinv );

      gen.parameters().add( "aam", mass );
      gen.parameters().add( "aaf0", 1./finv );
      gen.parameters().add( "aaw", mass*mass*mass*0.25*M_1_PI*finv*finv ); //FIXME
      hwpram_.IPRINT = 0;

      gen.initialise();
      const double xsec = gen.crossSection();
      g_scan_loglog.SetPoint( g_scan_loglog.GetN(), lmass, lfinv, xsec );
      g_scan_linlin.SetPoint( g_scan_linlin.GetN(), mass, finv, xsec );
      of << lmass << "\t" << lfinv << "\t" << xsec << "\n";
    }
  }
  g_scan_linlin.Write();
  g_scan_loglog.Write();

  file->Write();
  of.close();

  return 0;
}

