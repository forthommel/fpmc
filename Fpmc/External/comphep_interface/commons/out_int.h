/*$Log */
#ifndef __OUT_INT__
#define __OUT_INT__

#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef int (DNN) (double *);
typedef double (FNN) (void);
extern double sqrMom (char *, double *);
extern double computer_eps;
extern double Fmax;
extern int *calcCoef;
extern double DP[];

#endif
