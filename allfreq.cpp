/*
 * PEDIPET
 * 
 * allfreq.cpp
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
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

/*

  Likelihood calculation of general pedigrees
  for single markers

*/

/*

Possible options for speed optimizations:
  1) Remove calculation of combinations from each step in
     calculating the likelihood
  2) Improve/correctly implement Langes step B
  3) Remove pair of untyped parents with only a single typed son
     Remember to update the offspring to no parents
  4) Remove untyped, non-founder individual without any offspring

*/

#include <cmath>
#include <cstring>
#include <cstdio>
#include "allfreq.h"
#include "time.h"


// Define EM convergence criteria for log likelihoods
#define EMCONVCRIT  0.001
#define EMLAG       3
#define UPSCALE     1000

#define XLINKED

// This holds the number of non-founder combination that is allowed
// before reducing genotypes for each founder combination
// This does not currently work
#define REDUCE_PEDIGREE

#ifdef REDUCE_PEDIGREE
#define NFCOMBINATIONS   63
#endif

// #define debug

// Show time spend in function
// #define SHOWTIME

//
// Calculates the expected allele distribution for this pedigree
//
MATRIX ExpectedAlleleDistribution(individual *indlist, MATRIX theta, double *loglike) 
{
  MATRIX newtheta(theta.Rows(),1), currenttheta(theta.Rows(),1);
  individual *ind, *newlist;
  unsigned long int combinations, i, j, nfcomb, fcomb, newcomb, orignfcomb;
  double checkcomb, total, res, ll;
  markerlist *mkr, *fmkr, *mmkr;
  int al1, al2;

  // Calculate the number of combinations
  // required to evaluate
  combinations = 1;
  fcomb = 1;
  nfcomb = 1;
  checkcomb    = 0;
  for (ind = indlist; ind; ind = ind->next) {    
    combinations *= ind->tmpint2;    
    checkcomb += log(ind->tmpint2);
    if (founder(ind)) {
      //      printf("%s is a founder   (%d)\n",ind->id, ind->tmpint2);
      fcomb *= ind->tmpint2;
    }
    else 
      nfcomb *= ind->tmpint2;
  }

  if (fabs(log(combinations)-checkcomb) > .001) {
    printf("ERROR: Too many combinations (%ld) for me to handle in pedigree %s\n", combinations, indlist->pedigree);
    //    halt(1);
  }

  newlist = indlist;

  newtheta = 0;
  total = 0;
  res = 0;
  orignfcomb = nfcomb;

#ifdef debug
  printf("Combinations for pedigree %s : %ld (founders: %ld, non-founder: %ld)\n", indlist->pedigree, combinations, fcomb, nfcomb);
#endif



  // Go through all possible founder combinations
  for (i = 0; i < fcomb; i++) {

    newlist = indlist;
    // Number of combination in new list
    // This changes below
    newcomb = 1;

#ifdef REDUCE_PEDIGREE
    // Should only reduce genotypes if there are many non-founder combinations
    // Remember to also change the check below

    if (orignfcomb>NFCOMBINATIONS) {
#ifdef debug      
      //      printf("Reducing pedigree\n");
#endif
      newlist = CopyIndividualList(indlist);
      // Removes the founders genotypes
      for (ind = newlist; ind; ind = ind->next) {      
	if (founder(ind)) {
	  mkr = markernumber(ind, ind->tmpint1);
	  al1 = mkr->allele1;
	  al2 = mkr->allele2;
	  freelist(ind->marker);
	  ind->marker = 0;
	
	  mkr = (markerlist *) cmalloc (sizeof (markerlist));
	  memset (mkr,0,sizeof (markerlist));
	  mkr->allele1 = al1;
	  mkr->allele2 = al2;
	  
	  addlist(&ind->marker, mkr);
	  ind->tmpint1 = 1;
	  ind->tmpint2 = 1;
	}
      }
      ReduceGenotypes(newlist);

      nfcomb = 1;
      for (ind = newlist; ind; ind = ind->next) {      
	//	printf("%s has %d options\n", ind->id, listlen(ind->marker));
	newcomb *= listlen(ind->marker);
	ind->tmpint2 = listlen(ind->marker);  // Maximum number of combinations
	if (!founder(ind)) {
	  //	  printf("%s has %d options\n", ind->id, listlen(ind->marker));
	  nfcomb *= listlen(ind->marker);
	}
      }
    }
      
#endif
    // Go through all possible combinations
    // and calculate the likelihood of each
    // combination. These are used for calulating
    // the expected scores in the EM-algorithm

    // Illegal pedigree
    if (newcomb > 0) {

    for (j = 0; j < nfcomb; j++) {
      currenttheta = 0;

        
      res = 0;

      // Calculate the likelihood of this combination
      for (ind = newlist; ind; ind = ind->next) {      
	mkr = markernumber(ind, ind->tmpint1);
	if (founder(ind)) {
	  //	  printf("F\n");
	  res += log(theta(mkr->allele1,1))+ log(theta(mkr->allele2,1));

	  //	  printf("RES : %f  ===> %f   (%s)\n", res, total, ind->id);

	  currenttheta(mkr->allele1,1) += 1;
	  currenttheta(mkr->allele2,1) += 1;

	  // If the founder is heterozygote, then we have removed one of the 
	  // combinations and should therefore multiply by 2
	  // HOWEVER, that is not true for parental founder for sex-linked markers
	  if (options->Index[O_XLINKED] && ind->sex == S_MALE) {
	    // Do nothing
	  }
	  else {	   
	    if (mkr->allele1 != mkr->allele2)
	      res += log(2);
	  }
	}
	else {
	  // Non-founder
	  //	  printf("O\n");	  
	  fmkr = markernumber(ind->father,ind->father->tmpint1);
	  mmkr = markernumber(ind->mother,ind->mother->tmpint1);


	  //	  printf("RES : %f  ===> %f   (%s)\n", res, total, ind->id);

	  //	  printf("OO     : %d %d    %d %d    %d %d %s\n", mkr->allele1, mkr->allele2, fmkr->allele1, fmkr->allele2, mmkr->allele1, mmkr->allele2, ind->father->id);

	  //	  printf("OO+\n");	  	  
	  // Check Mendel
	  if (mkr->allele1 != fmkr->allele1 && mkr->allele1 != fmkr->allele2) {
	    res = 1;
	    break;
	  }
	  //	  printf("OO-\n");	  
	  if (mkr->allele2 != mmkr->allele1 && mkr->allele2!= mmkr->allele2) {
	    res = 1;
	    break;
	  }

	  //	  printf("OO0\n");	  

	  if (heterozygous(ind->father,ind->father->tmpint1))
	    res += log(.5);
	  if (heterozygous(ind->mother,ind->mother->tmpint1))
	    res += log(.5);
	}
	//	printf("u\n");
      }

      //      printf("Hopsa\n");

      // add to theta
      if (res<1) {
	newtheta += currenttheta*(exp(res)*UPSCALE); 
	total += UPSCALE*exp(res);
	
      }
      //      printf("res : %f\n", res);
      //      printf("Total : %f\n", total);

      // Update the array for the non-founders
      for (ind = newlist; ind; ind = ind->next) {
	// At the end for this person
	if (!founder(ind)) {
	  if (ind->tmpint1 == ind->tmpint2) {
	    ind->tmpint1 = 1;
	  }
	  else {
	    ind->tmpint1++;
	    break;
	  }
	}
      }
    }


    }

#ifdef REDUCE_PEDIGREE
    // Should be the same check as above
    if (orignfcomb>NFCOMBINATIONS)
      FreeIndividualList(newlist);
#endif

    // Update the array
    for (ind = indlist; ind; ind = ind->next) {
      if (founder(ind)) {
	// At the end for this person
	if (ind->tmpint1 == ind->tmpint2) {
	  ind->tmpint1 = 1;
	}
	else {
	  ind->tmpint1++;
	  break;
	}
      }
    }
  }
  newtheta *= 1.0/total;

  ll = log(total)-log(UPSCALE);
  *loglike = ll;

  for (ind = indlist; ind; ind = ind->next) {
    ind->tmpint1 = 1;
  }

  return newtheta;
}


void MLAlleleFrequencyEstimation(individual *listofpersons, int markernum)
{
  individual *gelist, *ind, *ind2, *indlist, **pOptList;
  namelist *pedigrees, *ped;
  markerlist *mkr, *mkr2, *mkr3;
  int classtype, i, j, nPedigrees, nFounders, nAlleles, removedone;
  double deltaloglike, loglike, startloglike;
  double *lag;
  FreqList *fl;

#ifdef SHOWTIME
  GetTime(1);
#endif

  lag = new double[EMLAG];

  CopyAlleleFreq(markernum);
  nAlleles = NumberOfAlleles(markernum);

#ifdef XLINKED
  if (options->Index[O_XLINKED]) {
    // Reduce the allele frequency count
    // This is just an ad hoc start value
    // Assuming that half the founder are men

    for (j=1; j<=NumberOfAlleles(markernum); j++) {
      allelefreq[j] *= .75;
    }    
    nAlleles++;
    allelefreq[nAlleles] = .25;
  }
#endif

  MATRIX theta(nAlleles,1), newtheta(nAlleles,1);

  // Make a copy of the individual list
  indlist = CopyIndividualList(listofpersons);

  // Make a list of the pedigrees to include
  pedigrees =  MakePedigreeList(indlist);
  nPedigrees = listlen(pedigrees);

  pOptList = new (individual *[nPedigrees]);

  // Start by preparing the data
  nFounders = 0;

  // Remove non-founders that have no offspring and aren't genotyped 
  for (ind = indlist; ind; ) {
    ind2 = ind->next;
    if (ind->offspring == NULL && !istypedformarker(ind, markernum)) {
      RemoveIndividual(indlist, ind);      
    }
    ind = ind2;
  }

  for (ped = pedigrees, i=0; ped; ped = ped->next, i++) {
    // Assume all founders are genotyped
    classtype = 1; 
    gelist = GenotypeElimination(ped->name,markernum);

    for (ind = gelist; ind; ind = ind->next) {
      // Start by reducing heterozygous genotypes for founders
      if (founder(ind)) {
	nFounders++;
        for (mkr = ind->marker; mkr; mkr = mkr->next) {
          for (mkr2 = mkr->next; mkr2; ) {
	    mkr3 = mkr2->next;
            if (mkr->allele1 == mkr2->allele2 && mkr->allele2 == mkr2->allele1) {
              removelist(&ind->marker, mkr2);
	      break;
            }
	    mkr2 = mkr3;
          }
        }
        // How many combination are left?
        if (listlen(ind->marker)>1) {
	  // More than one: all founders not genotyped
          classtype = 0;
        }
      }

      // Initializes the tmpint variables to speed up things later
      ind->tmpint1 = 1;                     // Current combination
      ind->tmpint2 = listlen(ind->marker);  // Maximum number of combinations
    }

    // Remove all non-founders for pedigrees
    // where all founder are typed (or deduced)
    // Remember to remove from the bottom up!
    if (classtype == 1) {
      removedone = 1;
      while (removedone) {
	removedone = 0;
	for (ind = gelist; ind; ) {
	  ind2 = ind->next;
	  if (!founder(ind) && !ind->offspring) {
	    removedone = 1;
	    RemoveIndividual(gelist, ind);
	    // Fix things if we removed the first person
	    if (gelist == ind)
	      gelist = ind2;
	  }        
	  else {
	    ind->tmpint1 = 1;                     // Current combination
	    ind->tmpint2 = listlen(ind->marker);  // Maximum number of combinations
	  }
	  ind = ind2;
	}
      }
    }
    pOptList[i] = gelist;
  }
  // Now do an EM-algorithm for finding the maximum
  // Initialize by using the allele frequencies
  for (i = 1; i<= nAlleles; i++) {
    theta(i,1) = allelefreq[i];
  }
  startloglike = 0;
  j = 0;
  while(1) {
    j++;
    printf("Iteration : %3d",j);

    loglike = 0;
    deltaloglike = 0;
  
    newtheta = 0;
    for (i = 0; i < nPedigrees; i++) {
      newtheta += ExpectedAlleleDistribution(pOptList[i], theta, &loglike);
#ifdef debug      
      printf("  -> Pedigree %d (%f)\n",i+1, loglike);
#endif      

      deltaloglike += loglike;
    }
    newtheta *= (1.0/(2.0*nFounders));
    theta = newtheta;

    if (j == 1) {
      startloglike = deltaloglike;
      for (i = 0; i<EMLAG; i++)
	lag[i] = deltaloglike + 2*EMCONVCRIT;
    }

    deltaloglike -= startloglike;

    printf(" : Delta log-like: %f\n", deltaloglike);

    // Change in ll since start
    if (fabs(lag[EMLAG-1]-deltaloglike)<EMCONVCRIT)
      break;

    // Fix the lags
    for (i = EMLAG-1; i; i--)
      lag[i]= lag[i-1];
    lag[0] = deltaloglike;
  }

  printf("ML estimates of allele frequencies:\n");
#ifdef XLINKED
  if (options->Index[O_XLINKED]) {
    newtheta *= 1.0/(1.0-newtheta(newtheta.Rows(),1));
    MATRIX save = newtheta;
    newtheta.Resize(newtheta.Rows()-1,1);
    for (i = 1; i <= newtheta.Rows(); i++) {
      newtheta(i,1) = save(i,1);
    }

    long fndr=0;
    long mfounders=0;

    for (ind = listofpersons; ind; ind = ind->next) {
      if (founder(ind)) {
	fndr++;
	if (ind->sex == S_MALE)
	  mfounders++;
      }
    }
    // Verify pseudo allele frequency
    //    printf("%f  === %f\n", (double) mfounders/(double) (2.0*fndr), save(save.Rows(),1));
  }  
#endif  
  Transpose(newtheta).Print();



  // Converged. Now set new allele frequencies
  fl = FrequencyNumber(markernum);
  fl->num_alleles = (int) newtheta.Rows();

  for (j = 1 ; j<=fl->num_alleles; j++) {
    fl->frequency[j] = newtheta(j,1);
  }

  delete[] lag;

  for (i = 0; i<nPedigrees; i++)
    FreeIndividualList(pOptList[i]);
  delete[] pOptList;
  FreeNameList(pedigrees);
  FreeIndividualList(indlist);

#ifdef SHOWTIME
  printf("Time spend (in s) : %d\n", GetTime());
#endif
}


