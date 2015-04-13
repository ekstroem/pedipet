/*

  Statistical distributions
  C-source by Claus Ekstrøm 1996-1999
  Based on StatUnit source by Tue Tjur


*/

#include <stdio.h>
#include <math.h>

#include "distr.h"
#include "rmath.h"


// Computes the logarithm of the gamma function at f/2
double lnGamma(long f)
{
  double sum, y;
  long k;

  y = (double) f/2.0;

  if (f>600)
  {
    sum = log(2*PI)/2.0  + (y-0.5)*log(y) -y + 1/(12.0*y) - 1/360.0/y/y/y;
    return sum;
  }
  else
  {
    k = f;
    sum = 0;

    while (k>2)
    {
      k   -= 2;
      sum += log(0.5*k);
    }
    if (k==1)
      sum += 0.5*log(PI);
    return sum;
  }
}

// Returns the right tail probability in the gamma
// distribution with lambda = f/2.

double pGamma(long f, double y)
{

  double term, sum;
  long k;

  if (y<=0)
    return 1;

  // CHECK if long=> double
  /*  if (y<f/2.0 || y<42)
  {
  */
    term = (f/2.0)*log(y) - y - lnGamma(f+2);
    if (term>-1000)
      term = exp(term);
    else
      term = 0;

    sum = 0.0;
    k = 0;

    while ( (f+k)*term>(f+k-2*y)*1E-20)
    {
      sum += term;
      term *= 2*y/((double)(f+k+2));
      k += 2;
    }
    return fabs(1-sum);
    /*
  }
  else
  {
    term = (0.5*f-1)*log(y)-y-lnGamma(f);
    if (term >-1000)
      term = exp(term);
    else
      term = 0;

    sum = 0.0;
    k = 0;

    while (term*y > (2*y-f+k)*0.5E-20 && f-k>1)
    {
      sum += term;
      k += 2;
      term *= 0.5*(f-k)/y;
    }
    return fabs(sum);

  }
    */
}

//  Returns the right tail probability in the chi square distribution
//  with f degrees of freedom.
double pChi2(long f, double y)
{
  return (pGamma(f,0.5*y));
}

//  Returns the right tail probability in the standardized normal distribution
double pNorm(double y)
{
  if (y>0)
    return (pGamma(1,0.5*y*y)/2.0);
  else
    return (1-(pGamma(1,0.5*y*y))/2.0);
}


// Returns the left tail probability of the beta distribution
// with paramters lambda1=f1/2 and lambda2=f2/2. Use only f1+f2<40.
// Accuracy around +/- 1E-16 .
double pBeta0(long f1, long f2, double y)
{
  double sum, term;
  long k;

  sum = 0.0;
  k = 0;

  term = exp(lnGamma(f1+f2)-lnGamma(f2) -lnGamma(f1+2)+f1*log(y)*.5);

  while (k<f2 || fabs(term)>1E-20)
  {
    sum += term;
    k += 2;
    term = -term*y*(f2-k)*(f1+k-2)/(double)k/(f1+k);
  }
  return sum;
}


// Returns the LEFT tail probability in the beta distribution
// with paramters lambda1=f1/2 and lambda2=f2/2. Use only
// f1 and f2 < 1E6.
double pBeta(long f1, long f2, double y)
{
  double sum, term;
  long k;
  int intch;

  if (f1==f2 && y==0.5)
    return 0.5;
  else
  {
    if (y<=0)
      return 0;
    else if (y>=1)
      return 1;
  }

  intch = 0;

  if (y>1-y)
  {
    intch = 1;
    k = f1;
    f1 = f2;
    f2 = k;
    y = 1-y;
  }

  if (f1+f2<41)
    sum = pBeta0(f1,f2,y);
  else
  {
    term = (f2/2.0-1)*log(1-y) + (f1/2.0)*log(y) + lnGamma(f1+f2) - lnGamma(f1+2) - lnGamma(f2);
    if (term>-1000)
      term = exp(term);
    else
      term = 0;

    if (term<1E-35 && (y<f1/(double)(f1+f2)))
      sum = 0.0;
    else
      if (term<1E-35 && (y>f1/(double)(f1+f2)))
        sum = 1.0;
      else
      {
        k = 0;
        sum = 0.0;

        while (fabs(term)>1E-25 || y*(f2-k) > (1-y)*(f1+k))
        {
          sum += term;
          k += 2;
          term *= y*(f2-k)/((1-y)*(double)(f1+k));
        }
      }
  }

  if (intch)
    sum = 1- sum;
  return fabs(sum);
}


// Returns the right tail probability in the T distribution.
// Use only  f < 1E6.
double pTdistr(long f, double y)
{
  double p;

  if (y==0.0)
    return 0.5;
  else
  {
    p = (double) f/(y*y+f);
    p = pBeta(f,1,p);
    p /=2.0;

    if (y<0)
      p = 1-p;
    return p;
  }
}

double pNormalInv(double p)
{
  // Should really call the source code from R using:
  return(qnorm5(p, 0.0, 1.0, 0, 0));
  // but currently I stick to the old Tjur code:

  /*
  double pp,y,a,b,y0;
  int sign;
  y= 0;  
  y0=1;
  pp=0.5;

  if (p=pp)
    return(0);
  if (p<pp)
    sign = -1;
  else 
    sign = 1;
  



  // Klytløsning
  if (p < 1E-9) {
    return(5.997807);
  }
  while (y0>1E-5) {
    y0=y;
    // log normal
    a=-.5*log(2*PI)-.5*y*y;
    b =y;    
    if (fabs(b)<1E-2)
      y+=(pp-p)*exp(-a);
    else {
      y+=log(fabs(1+b*(pp-p)*exp(-a)))/b;
    }
    pp=pNorm(y); 

    if (pp-p)

    //    printf("-y:   %f \n",y0-y);
    //    fflush(stdout);


    y0=fabs(y-y0);
  }

  //    printf("YYYY:   %f     %f\n",y,pp);
  */
  //  return(y);
}


double pGammaInv(long f, double p)
{
  double pp, y, y0, a, b, a0;
  y = 0;

  a0 = -lnGamma(f);

  if (f == 1) {
    y = pNormalInv(p/2);
    y = 0.5*y*y;
  }
  else {
    if (f>100) {
      y = sqrt(2*f-1)+pNormalInv(p);
      y = 0.25*y*y;      
    }
    else {
      y = .5*f;
      y0 = 1;
      pp = pGamma(f,y);
      while (y0>1E-7) {
	y0 = y;
	a = a0 + (.5*f-1)*log(y)-y;
	b = (.5*f-1)/y-1;
	if (fabs(b*(pp-p)*exp(-a)) < 1E-5 ) {
	  y += (pp-p)*exp(-a);	  
	}
	else {
	  y += log(1+b*(pp-p)*exp(-a))/b;
	}
        pp = pGamma(f,y);
        y0 = fabs(y-y0);
      }
    }

  }
  return y;
}

double pChi2Inv(long f, double p) 
{
  double y;

  y=pGammaInv(f,p);
  return(2*y);
}
