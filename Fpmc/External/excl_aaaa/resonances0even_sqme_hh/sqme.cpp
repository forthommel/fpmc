#include <iostream>
#include <fstream>
#include <math.h>
#include "helicity_amplitudes.h"

using namespace std;

namespace resonances0even_hh
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

    re *= -8*alpha_em*alpha_em;
    im *= -8*alpha_em*alpha_em;
    // the factor of - 8 is needed because of the conventions in
    // Costantini, DeTollis, Pistoni
  }

  /////////////////////////////// Z Z final state from OZ operator

  // Computes the exotic spin 0 even  squared matrix element and the SM interference
  double
  sqme( double s, double t, int exclude_loops_SM, int exclude_loops_EX, double m, double c, double ch, double w_c, double aa )
  {
    double mass = m;
    double f0 = c;
    double f0H = ch;
    double w_const = w_c;
    double a2 = aa;

    if ( s < 0 || t > 0 || t < -s ) {
      cerr << "Invalid domain. Valid range is s>=0 and -s<=t<=0" << endl;
      return 0;
    }

    // read in input data

    double re_ex, im_ex, value = 0.;

    // Mpppp:
    // the exotic matrix element:
    MhhOHpppp_spin0even(s, t, mass, f0, f0H, w_const,a2, re_ex, im_ex);
    value += re_ex*(re_ex) + im_ex*(im_ex) ;

    // repeat for the other helicities

    // Mppmm:
    MhhOHppmm_spin0even(s, t, mass, f0, f0H, w_const,a2, re_ex, im_ex);
    value += re_ex*(re_ex) + im_ex*(im_ex) ;

    // Mpmmp:
    MhhOHpmmp_spin0even(s, t, mass, f0, f0H, w_const,a2, re_ex, im_ex);
    value += re_ex*(re_ex) + im_ex*(im_ex) ;

    // Mpmpm:
    MhhOHpmpm_spin0even(s, t, mass, f0, f0H, w_const,a2, re_ex, im_ex);
    value += re_ex*(re_ex) + im_ex*(im_ex) ;

    // Mpppm
    MhhOHpppm_spin0even(s, t, mass, f0, f0H, w_const,a2, re_ex, im_ex);
    value += re_ex*(re_ex) + im_ex*(im_ex) ;

    // Mpp00
    MhhOHpp00_spin0even(s, t, mass, f0, f0H, w_const,a2, re_ex, im_ex);
    value += ( re_ex*re_ex + im_ex*im_ex );

    return 0.5*value;
  }
}
