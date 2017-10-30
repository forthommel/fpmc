#include <math.h>

//---M.Saimpert 01/2014--matthias.saimpert@cern.ch--------//
//--------Coding of the exclusive aa->aa process----------//
//--formulas from G. von Gersdorff (gersdorff@gmail.com)--// 
//--formulas from S. Fichet  sylvain.fichet@gmail.com--)--// 
//modification of the comphep external module used for aaww/aazz//

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif
  // routines called by fpmc 

  // SM
  void sm_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops);
  // Exotic fermions
  void bsmf_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& q, double& n );
  // Exotic vectors
  void bsmv_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& q, double& n );
  // Spin-0-even neutral resonances
  void resonances0even_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& w_c, double& aa );
  // Spin-0-even neutral resonances
  void resonances0even_sqme_zz_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& _cz, double& w_c, double& aa );
  // Spin-0-even neutral resonances
  void resonances0evenoh_sqme_zz_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& _cz, double& w_c, double& aa );
  // Spin-0-even neutral resonances
  void resonances0even_sqme_ww_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& _cw, double& w_c, double& aa );
  // Spin-0-even neutral resonances
  void resonances0even_sqme_hh_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& _cw, double& w_c, double& aa );
  // Spin-0-even neutral resonances
  void resonances0even_sqme_gluglu_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& _cw, double& w_c, double& aa );
  // Spin-0-even neutral resonances
  void resonances0even_sqme_az_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& _czg, double& w_c, double& aa );
  // Spin-2 neutral resonances
  void resonances2_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& w_c, double& aa );
  // EFT limit
  void eft_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, double& zeta1, double& zeta2, double& cutoff );
  void eft_sqme_aaaz_c_( double&, double&, double&, int&, double&, double&, double& );

#ifdef __cplusplus
}
#endif

//////////////////////////////////////////////////////////////////// 
// SM aaaa including:
//  fermion+W loop (exclude_loops=0), no fermions (1), no W (2)
//////////////////////////////////////////////////////////////////// 
namespace sm_aaaa { extern double sqme( double, double, int ); }
void
sm_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops )
{
  mp2 = sm_aaaa::sqme(s, t, exclude_loops);
}
#define sm_sqme_aaaa_c__ sm_sqme_aaaa_c_ //wrapper for g77

//////////////////////////////////////////////////////////////////// 
// BSM aaaa with exotic fermions of mass m charge q multiplicity n
// including interference with SM (SM fermion loops excluded is default)
//////////////////////////////////////////////////////////////////// 
namespace bsmf_aaaa { extern double sqme( double, double, int, int, double, double, double ); }
void
bsmf_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& q, double& n )
{
  mp2 = bsmf_aaaa::sqme( s, t, exclude_loops_SM, exclude_loops_EX, m, q, n );
}
#define bsmf_sqme_aaaa_c__ bsmf_sqme_aaaa_c_ //wrapper for g77

//////////////////////////////////////////////////////////////////// 
// BSM aaaa with exotic vector of mass m charge q multiplicity n
// including interferences with SM
//////////////////////////////////////////////////////////////////// 
namespace bsmv_aaaa { extern double sqme( double, double, int, int, double, double, double ); }
void
bsmv_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& q, double& n )
{
  mp2 = bsmv_aaaa::sqme(s, t, exclude_loops_SM, exclude_loops_EX, m, q, n);
}
#define bsmv_sqme_aaaa_c__ bsmv_sqme_aaaa_c_ //wrapper for g77

//////////////////////////////////////////////////////////////////// 
// BSM aaaa with neutral resonance of mass m coupling c width w
// including interferences with SM
//////////////////////////////////////////////////////////////////// 
namespace resonances0even_aaaa { extern double sqme( double, double, int, int, double, double, double, double ); }
void
resonances0even_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& w_c, double& aa )
{
  mp2 = resonances0even_aaaa::sqme(s, t, exclude_loops_SM, exclude_loops_EX, m, c, w_c, aa);
}
#define resonances0even_sqme_aaaa_c__ resonances0even_sqme_aaaa_c_ //wrapper for g77
 
//////////////////////////////////////////////////////////////////// 
// BSM aaaa with neutral resonance of mass m coupling c width w
// including interferences with SM
//////////////////////////////////////////////////////////////////// 
namespace resonances0even_zz { extern double sqme( double, double, int, int, double, double, double, double, double ); }
void
resonances0even_sqme_zz_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& cz, double& w_c, double& aa )
{
  mp2 = resonances0even_zz::sqme(s, t, exclude_loops_SM, exclude_loops_EX, m, c, cz, w_c, aa);
}
#define resonances0even_sqme_zz_c__ resonances0even_sqme_zz_c_ //wrapper for g77

namespace resonances0evenoh_zz { extern double sqme( double, double, int, int, double, double, double, double, double ); }
void
resonances0evenoh_sqme_zz_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& cz, double& w_c, double& aa )
{
  mp2 = resonances0evenoh_zz::sqme( s, t, exclude_loops_SM, exclude_loops_EX, m, c, cz, w_c, aa );
}
#define resonances0evenoh_sqme_zz_c__ resonances0evenoh_sqme_zz_c_ //wrapper for g77

////////////////////////////////////////////////////////////////////
// BSM aaaa with neutral resonance of mass m coupling c width w
// including interferences with SM
////////////////////////////////////////////////////////////////////
namespace resonances0even_ww { extern double sqme( double, double, int, int, double, double, double, double, double ); }
void
resonances0even_sqme_ww_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& cw, double& w_c, double& aa )
{
  mp2 = resonances0even_ww::sqme(s, t, exclude_loops_SM, exclude_loops_EX, m, c, cw, w_c, aa);
}
#define resonances0even_sqme_ww_c__ resonances0even_sqme_ww_c_ //wrapper for g77

////////////////////////////////////////////////////////////////////
// BSM aahh with neutral resonance of mass m coupling c width w
// including interferences with SM
////////////////////////////////////////////////////////////////////
namespace resonances0even_hh { extern double sqme( double, double, int, int, double, double, double, double, double ); }
void
resonances0even_sqme_hh_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& cw, double& w_c, double& aa )
{
  mp2 = resonances0even_hh::sqme( s, t, exclude_loops_SM, exclude_loops_EX, m, c, cw, w_c, aa );
}
#define resonances0even_sqme_hh_c__ resonances0even_sqme_hh_c_ //wrapper for g77

////////////////////////////////////////////////////////////////////
// BSM gluon-gluon with neutral resonance of mass m coupling c width w
// including interferences with SM
////////////////////////////////////////////////////////////////////
namespace resonances0even_gluglu { extern double sqme( double, double, int, int, double, double, double, double, double ); }
void
resonances0even_sqme_gluglu_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& cw, double& w_c, double& aa )
{
  mp2 = resonances0even_gluglu::sqme( s, t, exclude_loops_SM, exclude_loops_EX, m, c, cw, w_c, aa );
}
#define resonances0even_sqme_gluglu_c__ resonances0even_sqme_gluglu_c_ //wrapper for g77

////////////////////////////////////////////////////////////////////
// BSM aaaa with neutral resonance of mass m coupling c width w
// including interferences with SM
////////////////////////////////////////////////////////////////////
namespace resonances0even_az { extern double sqme( double, double, int, int, double, double, double, double, double ); }
void
resonances0even_sqme_az_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& czg, double& w_c, double& aa )
{
  mp2 = resonances0even_az::sqme( s, t, exclude_loops_SM, exclude_loops_EX, m, c, czg, w_c, aa );
}
#define resonances0even_sqme_az_c__ resonances0even_sqme_az_c_ //wrapper for g77

//////////////////////////////////////////////////////////////////// 
// BSM aaaa with neutral resonance of mass m coupling c width w
// including interferences with SM
//////////////////////////////////////////////////////////////////// 
namespace resonances2_aaaa { extern double sqme( double, double, int, int, double, double, double, double ); }
void
resonances2_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, int& exclude_loops_EX, double& m, double& c, double& w_c, double& aa )
{
  mp2 = resonances2_aaaa::sqme( s, t, exclude_loops_SM, exclude_loops_EX, m, c, w_c, aa );
}
#define resonances2_sqme_aaaa_c__ resonances2_sqme_aaaa_c_ //wrapper for g77

/////////////////////////////////////////////////////////////////// 
// Anomalous aaaa in the EFT limit parametrized by z1 and z2
// including interferences with SM
//  cutoff_scale means the scale cutoff - if < 0, no formfactor used
//////////////////////////////////////////////////////////////////// 
namespace eft_aaaa { extern double sqme( double, double, int, double, double ); }
void
eft_sqme_aaaa_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, double& zeta1, double& zeta2, double& cutoff )
{
  // apply cutoff  for anomalous couplings 
  //R.S. Gupta arXiv:1111.3354 [hep-ph] 
  double fact = 1.;
  if ( cutoff > 0 ) fact = 1/( 1+pow( s/ cutoff / cutoff, 2 ) );
  mp2 = eft_aaaa::sqme( s, t, exclude_loops_SM, zeta1*fact, zeta2*fact );
}

/////////////////////////////////////////////////////////////////// 
// Anomalous aaaz in the EFT limit parametrized by z1 and z2
// including interferences with SM
//  cutoff_scale means the scale cutoff - if < 0, no formfactor used
//////////////////////////////////////////////////////////////////// 
namespace eft_aaaz { extern double sqme( double, double, int, double, double ); }
void
eft_sqme_aaaz_c_( double& mp2, double& s, double& t, int& exclude_loops_SM, double& zeta1, double& zeta2, double& cutoff )
{
  // apply cutoff  for anomalous couplings 
  //R.S. Gupta arXiv:1111.3354 [hep-ph] 
  double fact = 1.;
  if (cutoff > 0 ) fact = 1/( 1+pow( s/ cutoff / cutoff, 2 ) );
  mp2 = eft_aaaz::sqme( s, t, exclude_loops_SM, zeta1*fact, zeta2*fact );
}

#define eft_sqme_aaaa_c__ eft_sqme_aaaa_c_ //wrapper for g77

