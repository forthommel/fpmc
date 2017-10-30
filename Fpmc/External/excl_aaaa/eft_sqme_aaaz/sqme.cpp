#include <iostream>
#include <fstream>
#include <math.h>
#include "helicity_amplitudes.h"

using namespace std;

namespace eft_aaaz
{
  void
  me_SM( void (*me)( double, double, double&, double&, int ), double s,double t, double& re, double& im, int exclude_loops )
  {
    // This routine computes the complex SM amplitude
    // The first argument can be any of the helicity amplitudes Mpppp,Mppmm,Mpmpm,Mpmmp,Mpppm

    double d_re, d_im;

    re = im = 0.;

    for ( unsigned int i = 0; i < 9; ++i ) {
      me(s/(4*SM_masses[i]*SM_masses[i]),t/(4*SM_masses[i]*SM_masses[i]), d_re, d_im, exclude_loops);
      re += d_re * SM_weight[i];
      im += d_im * SM_weight[i];
    }

    // Add also the W contribution 
    const double norm = 0.25/mW/mW;

    if      ( me == Mpppp_fermion ) Mpppp_vector( s*norm, t*norm, d_re, d_im, exclude_loops );
    else if ( me == Mppmm_fermion ) Mppmm_vector( s*norm, t*norm, d_re, d_im, exclude_loops );
    else if ( me == Mpmpm_fermion ) Mpmpm_vector( s*norm, t*norm, d_re, d_im, exclude_loops );
    else if ( me == Mpmmp_fermion ) Mpmmp_vector( s*norm, t*norm, d_re, d_im, exclude_loops );
    else if ( me == Mpppm_fermion ) Mpppm_vector( s*norm, t*norm, d_re, d_im, exclude_loops );

    re += d_re; 
    im += d_im;

    re *= 8*alpha_em*alpha_em;
    im *= 8*alpha_em*alpha_em;
    // the factor of 8 is needed because of the conventions in
    // Costantini, DeTollis, Pistoni
  }

  // Computes the  squared matrix element and the SM interference from free zeta_1, zeta_2
  double
  sqme( double s, double t, int exclude_loops_SM, double z1, double z3 )
  {
    double zeta1 = z1; // are in GeV^-4 units
    double zeta3 = z3; //

    if ( s < 0 || t > 0 || t < -s ) {
      cerr << "Invalid domain. Valid range is s>=0 and -s<=t<=0" << endl;
      return 0;
    }

    const double mZ2 = mZ*mZ;
    //Symmetry factor; returns squared matrix element
    return 4.
            *( 3.*zeta1*zeta1+3.*zeta3*zeta3-2.*zeta1*zeta3 )
            *( pow( s*s+t*t+s*t, 2 )
               - 2.*mZ2*( s+t )*( 3.*( s*s+t*t )+s*t )
               + mZ2*mZ2*( 3.*( s*s+t*t )+2.*s*t ) )
         + 16./3.*zeta1*zeta3*mZ*mZ*( mZ2-s-t )*s*t;
  }
}  //namespace eft_aaaa

