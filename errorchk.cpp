/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * Functions to check for errors
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA.
 */


#include <cmath>
#include <iostream>
#include <cstring>
#include <cstdio>
#include "errorchk.h"
#include "matrix.h"

//#define artikel
//#define fivealleles

#define newoutput

#ifdef artikel
#include <cstring2>
#endif

int IBS(int a1, int a2, int b1, int b2) {
  if ((a1==b1 && a2==b2) || (a1==b2 && a2==b1))
    return 2;
  if ((a1!=b1 && a2!=b2) && (a1!=b2 && a2!=b1))
    return 0;
  return 1;
}


void VerifyIndependence(individual *indlist, int all) {
  individual *ind,*ind2;
  markerlist *mkr, *mkr2;
  int i, j, k, l, n;
  double p2, p1, pval, ntests;
#ifdef artikel
  int *a, *b, *c, *d;
#endif

  ntests = 0;

  // Start by calculating alle the necessary mean IBS-scores

  double *meanval, *varval, *P2, *P1;

  meanval = new double[numberofmarkers()+1];
  varval = new double[numberofmarkers()+1];
  P2 = new double[numberofmarkers()+1];
  P1 = new double[numberofmarkers()+1];

#ifdef newoutput
  printf("Testing independence between individuals\n");
  printf("Comparing ");
  if (all)
    printf("all individuals\n");
  else 
    printf("only founders\n");
  printf("\nIndividuals                               t        p\n---------------------------------------------------------\n");

#endif

#ifdef artikel
  a = new int[numberofmarkers()+1];
  b = new int[numberofmarkers()+1];
  c = new int[numberofmarkers()+1];
  d = new int[numberofmarkers()+1];

  CopyAlleleFreq(1);
#endif


#ifdef artikel
  for (i = 1; i<= numberofmarkers() ; i++) {
    meanval[i] = 1.25;
    varval[i]  = 0.4375;
  }
#else  
  for (i = 1; i<= numberofmarkers() ; i++) {
    CopyAlleleFreq(i);
    meanval[i] = 0;
    p2 = 0;
    p1 = 0;
    for (j=1; j<= NumberOfAlleles(i); j++) {
      p2 += SQR(SQR(allelefreq[j]));
      for (k=1; k<= NumberOfAlleles(i); k++) {
	if (k==j) 
	  continue;
	p2 += 2*SQR(allelefreq[j])*SQR(allelefreq[k]);
	p1 += 4*SQR(allelefreq[j])*allelefreq[j]*allelefreq[k];
	for (l=1; l<= NumberOfAlleles(i); l++) {
	  if (l==j || l==k)
	    continue;
	  p1 += 4*SQR(allelefreq[j])*allelefreq[k]*allelefreq[l];
	}
      }
    }
    P2[i] = p2;
    P1[i] = p1;

    meanval[i] = 2*p2 + p1;
    varval[i] = 4*p2 + p1 - SQR(meanval[i]);
    //    printf("Mean IBS for marker %d is %f  (%f)\n", i, meanval[i], varval[i]);
  }
#endif

  // First look at founders only

  for (ind = indlist; ind; ind = ind->next) {
    if ((!all && !founder(ind)) || (!IsGenotyped(ind))) {
      continue;
    }
    for (ind2 = ind->next; ind2; ind2 = ind2->next) {
      if ((!all && !founder(ind2)) || !IsGenotyped(ind2)) {
	continue;
      }

      // Look at this pair
      double score = 0;
      n = 0;

      for (i = 1;i<=numberofmarkers() ; i++) {
	if (istypedformarker(ind, i) && istypedformarker(ind2, i)) {
	  n++;
	  mkr = (markerlist *) markernumber(ind,i);
	  mkr2 = (markerlist *) markernumber(ind2,i);

	  score += ((IBS(mkr->allele1, mkr->allele2, mkr2->allele1, mkr2->allele2) - meanval[i]) / sqrt(varval[i]));

#ifdef artikel
	  if (!strcmpl(ind->id, "6") && (!strcmpl(ind2->id, "10"))) {
	    a[i] = IBS(mkr->allele1, mkr->allele2, mkr2->allele1, mkr2->allele2);
	  }
	  if (!strcmpl(ind->id, "6") && (!strcmpl(ind2->id, "15"))) {
	    b[i] = IBS(mkr->allele1, mkr->allele2, mkr2->allele1, mkr2->allele2);
	  }
	  if (!strcmpl(ind->id, "10") && (!strcmpl(ind2->id, "11"))) {
	    c[i] = IBS(mkr->allele1, mkr->allele2, mkr2->allele1, mkr2->allele2);
	  }
	  if (!strcmpl(ind->id, "11") && (!strcmpl(ind2->id, "15"))) {
	    d[i] = IBS(mkr->allele1, mkr->allele2, mkr2->allele1, mkr2->allele2);
	  }
#endif
	}
      }
      ntests += 1;
      score /= sqrt(n);
      pval = pNorm(score);
      TextColour(COL_GREEN);
      if (pval < 0.05) 
	TextColour(COL_YELLOW);
      if (pval < 0.01) 
	TextColour(COL_RED);
#ifdef newoutput
      if (pval < 0.05) {
	printf("%-12s and %-12s           %6.3f (p = %5.3f)", ind->id, ind2->id, score, (score<0 ? 1 : pNorm(score)));
	
	// Check for inbreeding
	if (!strcmpl(ind->pedigree, ind2->pedigree)) {
	  printf("  I");	  	
	}
	printf("\n");
      }
#else
      printf("Indep.-test for %-12s and %-12s is %f (p = %5.3f)\n", ind->id, ind2->id, score, (score<0 ? 1 : pNorm(score)));
#endif
    }
  }

  TextColour(COL_WHITE);

#ifdef newoutput
  printf("\nA total of %-.0f pairwise tests were made\n", ntests);
  printf("Tests were based on %d markers\n", numberofmarkers());
#endif

#ifdef artikel
    MATRIX res(4,1);

    res = 0;

    for (i = 1;i<=numberofmarkers() ; i++) {
      res += Paper10Function(a[i],b[i],c[i],d[i]);
    }    
    res *=  1.0/(sqrt(numberofmarkers()));
    printf("Samlet test:\n");
    res.Print();
    
#endif


  delete[] meanval;
  delete[] varval;
  delete[] P2;
  delete[] P1;

#ifdef artikel
  delete[] a;
  delete[] b;
  delete[] c;
  delete[] d;
#endif
}


MATRIX Paper10Function(int a, int b, int c, int d) {
  int i;
  MATRIX chol(4,4), variance(4,4), mean(4,1), res(4,1);

  // Creates the vectors and matrices
  //
  // S_1 = 1,3
  // S_2 = 1,4
  // S_3 = 2,3
  // S_4 = 2,4
  //


#ifdef fivealleles
  mean = .6560;
  for (i =1; i<=4; i++) {
    variance(i,i) = .8 - SQR(.6560);
  }
  variance(1,2) = 0.016384;
  variance(2,1) = 0.016384;
  variance(1,3) = 0.016384;
  variance(3,1) = 0.016384;
  variance(2,4) = 0.016384;
  variance(4,2) = 0.016384;
  variance(3,4) = 0.016384;
  variance(4,3) = 0.016384;
#else

  mean = 1.25;
  for (i =1; i<=4; i++) {
    variance(i,i) = 0.4375;
  }
  variance(1,2) = 0.0625;
  variance(2,1) = 0.0625;
  variance(1,3) = 0.0625;
  variance(3,1) = 0.0625;
  variance(2,4) = 0.0625;
  variance(4,2) = 0.0625;
  variance(3,4) = 0.0625;
  variance(4,3) = 0.0625;
#endif

  // Cholesky(variance).Print();

  //  chol = Cholesky(Inverse(variance));
  chol = GaussJordan(Cholesky(variance));

  //  chol.Print();

  res(1,1) = a;
  res(2,1) = b;
  res(3,1) = c;
  res(4,1) = d;

  return(chol*(res - mean));
}


#undef artikel
