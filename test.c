#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "rmath.h"
#include "distr.h"



main() {

  printf("%15.13f\n", pNormalInv(.025));
  printf("%15.13f\n", pGammaInv(1, 0.000000001));
  printf("%15.13f\n", pChi2Inv(1,  0.0000000000000000000001));
}
