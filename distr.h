/*

  Statistical distributions
  C-source by Claus EKstrøm 1996-1999
  Based on StatUnit source by Tue Tjur


*/

#ifndef _DISTR_H_
#define _DISTR_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PI
#define PI 3.141592654
#endif


// Computes the logarithm of the gamma function at f/2
double lnGamma(long f);

// Returns the right tail probability in the gamma
// distribution with lambda = f/2.
double pGamma(long f, double y);

//  Returns the right tail probability in the chi square distribution
//  with f degrees of freedom.
double pChi2(long f, double y);

// Standardized norm.
// BAD HACK!
double pNorm(double y);

// Returns the left tail probability of the beta distribution
// with paramters lambda1=f1/2 and lambda2=f2/2. Use only f1+f2<40.
// Accuracy around +/- 1E-16 .
double pBeta0(long f1, long f2, double y);

// Returns the LEFT tail probability in the beta distribution
// with paramters lambda1=f1/2 and lambda2=f2/2. Use only
// f1 and f2 < 1E6.
double pBeta(long f1, long f2, double y);

// Returns the right tail probability in the T distribution.
// Use only  f < 1E6.
double pTdistr(long f, double y);

double pGammaInv(long f, double p);
double pChi2Inv(long f, double p);
double pNormalInv(double p);

#ifdef __cplusplus
}
#endif

#endif
