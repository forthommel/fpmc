#include<iostream>
#include<fstream>
#include<math.h>
#include"helicity_amplitudes.h"

using namespace std;

namespace sm_aaaa
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

  // compute the SM squared matrix element, including leptons, quarks and the W boson
  double
  sqme( double s, double t, int exclude_loops )
  {
    double re, im, value = 0.;

    if ( s < 0 || t > 0 || t < -s ) {
      cerr << "Invalid domain. Valid range is s>=0 and -s<=t<=0" << endl;
      return 0;
    }

    me_SM( Mpppm_fermion, s, t, re, im, exclude_loops );
    value += 4*( re*re+im*im );
   
    me_SM( Mppmm_fermion, s, t, re, im, exclude_loops );
    value += re*re+im*im;

    me_SM( Mpppp_fermion, s, t, re, im, exclude_loops );
    value += re*re+im*im;

    me_SM( Mpmmp_fermion, s, t, re, im, exclude_loops );
    value += re*re+im*im;

    me_SM( Mpmpm_fermion, s, t, re, im, exclude_loops );
    value += re*re+im*im;

    return 0.5 * value;  
  }
} //namespace sm_aaaa
