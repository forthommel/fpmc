#include <iostream>
#include <fstream>
#include <math.h>
#include "helicity_amplitudes.h"

using namespace std;

namespace bsmf_aaaa
{
  void
  me_SM( void (*me)( double, double, double&, double&, int ), double s,double t, double& re, double& im, int exclude_loops )
  {
    // This routine computes the complex SM amplitude
    // The first argument can be any of the helicity amplitudes Mpppp,Mppmm,Mpmpm,Mpmmp,Mpppm

    double d_re, d_im;

    re = im = 0;

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

  // Computes the exotic fermion squared matrix element and the SM interference
  double
  sqme( double s, double t, int exclude_loops_SM, int exclude_loops_EX, double m, double q, double n )
  {
    double mass = m;
    double charge = q;
    double multiplicity = n;

    if ( s < 0 || t > 0 || t < -s ) {
      cerr << "Invalid domain. Valid range is s>=0 and -s<=t<=0" << endl;
      return 0;
    }
    //M.S now parameters
    /*

    // read in input data

    ifstream data;
  
    data.open ("./data");
  
    data >> mass;
    data >> charge;
    data >> multiplicity;

    data.close();
    */

    double weight_exotic = multiplicity*pow(charge,4);

    // reduced Mandelstam variables
    double sred = s/(4 * mass * mass);
    double tred = t/(4 * mass * mass);   
  
    double re_ex, im_ex, re_SM, im_SM, value = 0.;

    // Mpppp:
    // the exotic matrix element:
    Mpppp_fermion(sred,tred, re_ex, im_ex, exclude_loops_EX);
    re_ex *= 8*alpha_em*alpha_em*weight_exotic;
    im_ex *= 8*alpha_em*alpha_em*weight_exotic;
    // the factor of 8 is needed because of the conventions in
    // Costantini, DeTollis, Pistoni
    // the SM matrix element:
    me_SM(Mpppp_fermion,s,t, re_SM, im_SM, exclude_loops_SM);
    value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;

    // repeat for the other helicities

    // Mppmm:
    Mppmm_fermion(sred,tred, re_ex, im_ex, exclude_loops_EX);
    re_ex *= 8*alpha_em*alpha_em*weight_exotic;
    im_ex *= 8*alpha_em*alpha_em*weight_exotic;
    me_SM(Mppmm_fermion,s,t, re_SM, im_SM, exclude_loops_SM);
    value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;

    // Mpmmp:
    Mpmmp_fermion(sred,tred, re_ex, im_ex, exclude_loops_EX);
    re_ex *= 8*alpha_em*alpha_em*weight_exotic;
    im_ex *= 8*alpha_em*alpha_em*weight_exotic;
    me_SM(Mpmmp_fermion,s,t, re_SM, im_SM, exclude_loops_SM);
    value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;

    // Mpmpm:
    Mpmpm_fermion(sred,tred, re_ex, im_ex, exclude_loops_EX);
    re_ex *= 8*alpha_em*alpha_em*weight_exotic;
    im_ex *= 8*alpha_em*alpha_em*weight_exotic;
    me_SM(Mpmpm_fermion,s,t, re_SM, im_SM, exclude_loops_SM);
    value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;

    // Mpppm
    Mpppm_fermion(sred,tred, re_ex, im_ex, exclude_loops_EX);
    re_ex *= 8*alpha_em*alpha_em*weight_exotic;
    im_ex *= 8*alpha_em*alpha_em*weight_exotic;
    me_SM(Mpppm_fermion,s,t, re_SM, im_SM, exclude_loops_SM);
    value += 4* (  re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) );
 
    return 0.5*value;
  }
}  //namespace bsmf_aaaa

