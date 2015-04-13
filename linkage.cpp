/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * Functions for two and threepoint linkage analysis
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

/********************************************************************
 *
 * Twopoint linkage analysis
 *
 ********************************************************************/


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>


#include "linkage.h"

//#define debug

// Lousyfix introduces a simple speedup, that will halve the number of combinations to search through
#define lousyfix  

extern "C" int founder(individual *ind);
extern "C" void adduniqueidlist (IDlist **l1, individual *p1);

// This is _really_ bad programming

individual *INDLIST;
int MARKER1, MARKER2, MARKER3;
IDlist *ID1, *ID2, *ID3;

double LogTwopointLikelihood(individual *indlist1, individual *indlist2, int markernum, int markernum2, double theta)
{
  individual *ind, *ind2;
  long int combin1, combin2, c1, c2;
  int fhettype, mhettype;
  double checkcomb, res, tempres;
  markerlist *mkr, *mkr2, *fmkr, *mmkr, *fmkr2, *mmkr2;

  // Determine the number of combinations required for each marker
  combin1 = 1;
  combin2 = 1;
  checkcomb = 0;
  for (ind = indlist1; ind; ind = ind->next) {
    ind->tmpint1 = 1;
    ind->tmpint2 = listlen(ind->marker);
    combin1 *= listlen(ind->marker);
    checkcomb += log(listlen(ind->marker));
  }
  // Check for integer overflow
  if (fabs(log(combin1)-checkcomb) > .001) {
    printf("ERROR: Too many combinations for me to handle in pedigree %s (%f %f)\n", indlist1->pedigree, (float) combin1, exp(checkcomb));
  }
  // Do the same for marker 2
  checkcomb = 0;
  for (ind = indlist2; ind; ind = ind->next) {
    // Remember that this a different list so we can assign the same variable
    // here as above
    ind->tmpint1 = 1;
    ind->tmpint2 = listlen(ind->marker);
    combin2 *= listlen(ind->marker);
    checkcomb += log(listlen(ind->marker));
  }
  // Check for integer overflow
  if (fabs(log(combin2)-checkcomb) > .001) {
    printf("ERROR: Too many combinations for me to handle in pedigree %s (%f %f)\n", indlist2->pedigree, (float) combin2, exp(checkcomb));
  }

  res = 0;

  //  printf("Combinations: %ld  %ld\n", combin1, combin2);

  for (c1 = 0; c1 < combin1; c1++) {
    for (c2 = 0; c2 < combin2; c2++) {
      tempres = 0;

      for (ind=indlist1; ind; ind=ind->next) {
	ind2 = FindIndListMember(indlist2, ind->id);

	mkr = markernumber(ind, ind->tmpint1);
	mkr2 = markernumber(ind2, ind2->tmpint1);

	// Calculate the likelihood for this combination
	if (founder(ind)) {
	  CopyAlleleFreq(markernum);
	  tempres += log(allelefreq[mkr->allele1])+log(allelefreq[mkr->allele2]);
	  CopyAlleleFreq(markernum2);
	  tempres += log(allelefreq[mkr2->allele1])+log(allelefreq[mkr2->allele2]);
	}
	else {
	  // Get the current genotypes of the parents
	  fmkr = markernumber(ind->father, ind->father->tmpint1);
	  mmkr = markernumber(ind->mother, ind->mother->tmpint1);
	  fmkr2= markernumber(ind2->father, ind2->father->tmpint1);
	  mmkr2= markernumber(ind2->mother, ind2->mother->tmpint1);

#ifdef debugs
	  printf("Person: %s\n", ind->id);
      	  printf("Father   : %d | %d\n", fmkr->allele1, fmkr->allele2);
	  printf("Father   : %d | %d\n", fmkr2->allele1, fmkr2->allele2);
	  printf("Mother   : %d | %d\n", mmkr->allele1, mmkr->allele2);
	  printf("Mother   : %d | %d\n", mmkr2->allele1, mmkr2->allele2);
	  printf("Offspring: %d | %d\n", mkr->allele1, mkr->allele2);
	  printf("Offspring: %d | %d\n", mkr2->allele1, mkr2->allele2);
	  printf("%d  <-> %d\n", c1, c2);
#endif

	  // Check to see if combination legal
	  if ((fmkr->allele1 != mkr->allele1 && fmkr->allele2 != mkr->allele1) || 
              (mmkr->allele1 != mkr->allele2 && mmkr->allele2 != mkr->allele2) || 
              (fmkr2->allele1 != mkr2->allele1 && fmkr2->allele2 != mkr2->allele1) || 
              (mmkr2->allele1 != mkr2->allele2 && mmkr2->allele2 != mkr2->allele2)) {
	    tempres = 0;
	    break;
	  }

	  // If the father is heterozygous at either locus then multiply by .5
	  fhettype = 1*(fmkr->allele1 != fmkr->allele2) + 10*(fmkr2->allele1 != fmkr2->allele2);
	  mhettype = 1*(mmkr->allele1 != mmkr->allele2) + 10*(mmkr2->allele1 != mmkr2->allele2);
	  
	  if (fhettype>0)
	    tempres += log(.5);
	  if (mhettype>0)
	    tempres += log(.5);

	  // Father is double heterozygous)
	  if (fhettype==11) {
	    if ((mkr->allele1 == fmkr->allele1 && mkr2->allele1 == fmkr2->allele2) || (mkr->allele1 == fmkr->allele2 && mkr2->allele1 == fmkr2->allele1)) {
	      // A recombination must have occured
	      tempres +=log(theta);
	    }
	    else 
	      tempres +=log(1-theta);
	  }

	  // Mother is Doub. Het.
   	  if (mhettype==11) {
	    if ((mkr->allele2 == mmkr->allele1 && mkr2->allele2 == mmkr2->allele2) || (mkr->allele2 == mmkr->allele2 && mkr2->allele2 == mmkr2->allele1)) {
	      // A recombination must have occured
	      tempres +=log(theta);
	    }
	    else 
	      tempres +=log(1-theta);
	    
	  }
	  //	  printf("Tempres: %f\n", tempres);	  
	}
      }
      res += exp(tempres);

	  //	  printf("Res    : %f\n", res);


	
      // Update the array for the second marker
      for (ind = indlist2; ind; ind = ind->next) {
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

    // Update the array for the first marker
    for (ind = indlist1; ind; ind = ind->next) {
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
  
  //  printf("RES    : %f\n", log(res));  
  return log(res);
}


double TwopointLinkage(individual *indlist, int marker1, int marker2, double theta)
{
  individual *ind, *ind2;
  namelist *Pedigrees, *ped;
  double res;

  res = 0;

  // Same marker?
  if (marker1 == marker2)
    return 0;

  Pedigrees = MakePedigreeList(indlist);  
  for (ped = Pedigrees; ped; ped = ped->next) {
    ind = GenotypeElimination(ped->name,marker1);
    ind2 = GenotypeElimination(ped->name,marker2);
    res += LogTwopointLikelihood(ind, ind2, marker1, marker2, theta);
    FreeIndividualList(ind);
    FreeIndividualList(ind2);
  }

  FreeNameList(Pedigrees);

  return res;
}


double TwopointLinkageAnalysis(individual *indlist, int marker1, int marker2, double theta)
{
  int i;

  double points[3], value[3];
  double lpoint, rpoint, lvalue, rvalue;

  points[0] = 0.0001;
  points[2] =.5;
  if (theta>points[0] && theta < .5)
    points[1] = theta;
  else
    points[1] = .25;

  for (i = 0; i<3; i++) {
    value[i] = TwopointLinkage(indlist, marker1, marker2, points[i]);
    //    printf("%f\n", value[i]);
  }
  
  i = 1;
  // Start bracketing the minimum
  while (value[1]<value[2] || value[1]<value[0] && i<10) {
    points[1] = 0.05*(double) i;
    value[1] = TwopointLinkage(indlist, marker1, marker2, points[1]);
    //    printf("%f\n", value[1]);
    i++;
  }

  if (value[0]>=value[1] && value[0]>=value[2])
    return 0;

  if (value[2]>=value[1] && value[2]>=value[0])
    return .5;

  while ((points[2]-points[0])>LINKAGE_PRECISION) {
    lpoint = .5*(points[1]+points[0]);
    lvalue = TwopointLinkage(indlist, marker1, marker2, lpoint);

    //    printf("%f  %f   ::: %f   %f\n", lpoint, lvalue, value[0], value[1]);

    if (lvalue<value[1]) {
      points[0] = lpoint;
      value[0] = lvalue;
    }
    else {
      points[2] = points[1];
      points[1] = lpoint;
      value[2] = value[1];
      value[1] = lvalue;
    }      

    rpoint = .5*(points[2]+points[1]);
    rvalue = TwopointLinkage(indlist, marker1, marker2, rpoint);

    //    printf("%f  %f   ::: %f   %f\n", rpoint, rvalue, value[1], value[2]);
    if (rvalue<value[1]) {
      points[2] = rpoint;
      value[2] = rvalue;
    }
    else {
      points[0] = points[1];
      points[1] = rpoint;
      value[0] = value[1];
      value[1] = rvalue;
    }      

    /*
    for (i = 0; i<3; i++)
      printf("%5.3f  ", points[i]);
    printf("\n");
    for (i = 0; i<3; i++)
      printf("%5.3f  ", value[i]);
    printf("\n");
    */
  }

  return points[1];
}


/************************

 * Here comes all the three point analyses
 * Not pretty but useful for my thesis

 ************************/


double LogThreepointLikelihood(individual *indlist1, individual *indlist2, individual *indlist3, int markernum, int markernum2, int markernum3, double gamma01, double gamma10, double gamma11)
{
  individual *ind, *ind2, *ind3;
  long int combin1, combin2, combin3, c1, c2, c3;
  int fhettype, mhettype;
  double checkcomb, res, tempres;
  markerlist *mkr, *mkr2, *mkr3, *fmkr, *mmkr, *fmkr2, *mmkr2, *fmkr3, *mmkr3;

  // Determine the number of combinations required for each marker
  combin1 = 1;
  combin2 = 1;
  combin3 = 1;
  checkcomb = 0;
  for (ind = indlist1; ind; ind = ind->next) {
    ind->tmpint1 = 1;
    ind->tmpint2 = listlen(ind->marker);
    combin1 *= listlen(ind->marker);
    checkcomb += log(listlen(ind->marker));
  }
  // Check for integer overflow
  if (fabs(log(combin1)-checkcomb) > .001) {
    printf("ERROR: Too many combinations for me to handle in pedigree %s (%f %f)\n", indlist1->pedigree, (double) combin1, exp(checkcomb));
  }
  // Do the same for marker 2
  checkcomb = 0;
  for (ind = indlist2; ind; ind = ind->next) {
    // Remember that this a different list so we can assign the same variable
    // here as above
    ind->tmpint1 = 1;
    ind->tmpint2 = listlen(ind->marker);
    combin2 *= listlen(ind->marker);
    checkcomb += log(listlen(ind->marker));
  }
  // Check for integer overflow
  if (fabs(log(combin2)-checkcomb) > .001) {
    printf("ERROR: Too many combinations for me to handle in pedigree %s (%f %f)\n", indlist2->pedigree, (float) combin2, exp(checkcomb));
  }

  checkcomb = 0;
  for (ind = indlist3; ind; ind = ind->next) {
    // Remember that this a different list so we can assign the same variable
    // here as above
    ind->tmpint1 = 1;
    ind->tmpint2 = listlen(ind->marker);
    combin3 *= listlen(ind->marker);
    checkcomb += log(listlen(ind->marker));
    //    printf("%19s : %d\n", ind->id, listlen(ind-marker);
  }
  // Check for integer overflow
  if (fabs(log(combin3)-checkcomb) > .001) {
    printf("ERROR: Too many combinations for me to handle in pedigree %s (%f %f)\n", indlist2->pedigree, (float) combin2, exp(checkcomb));
  }

  res = 0;

  //  printf("Total combinations: %ld  %ld  %ld\n", combin1, combin2, combin3);

  for (c1 = 0; c1 < combin1; c1++) {
    for (c2 = 0; c2 < combin2; c2++) {
      for (c3 = 0; c3 < combin3; c3++) {
	tempres = 0;

	for (ind=indlist1; ind; ind=ind->next) {
	  ind2 = FindIndListMember(indlist2, ind->id);
	  ind3 = FindIndListMember(indlist3, ind->id);

	  mkr = markernumber(ind, ind->tmpint1);
	  mkr2 = markernumber(ind2, ind2->tmpint1);
	  mkr3 = markernumber(ind3, ind3->tmpint1);

	  // Calculate the likelihood for this combination
	  if (founder(ind)) {
	    CopyAlleleFreq(markernum);
	    tempres += log(allelefreq[mkr->allele1])+log(allelefreq[mkr->allele2]);
	    CopyAlleleFreq(markernum2);
	    tempres += log(allelefreq[mkr2->allele1])+log(allelefreq[mkr2->allele2]);
	    CopyAlleleFreq(markernum3);
	    tempres += log(allelefreq[mkr3->allele1])+log(allelefreq[mkr3->allele2]);

#ifdef lousyfix
	    // Check if the founder could

	    //	    printf("===> %s   %d\n", ind->id, ind->tmpint4);

	    if (ind->tmpint4 == 1 && mkr->allele1 != mkr->allele2) {
	      tempres += log(2);
	    }
	    else if (ind->tmpint4 == 2 && mkr2->allele1 != mkr2->allele2) {
	      tempres += log(2);
	    }
	    else if (ind->tmpint4 == 3 && mkr3->allele1 != mkr3->allele2) {
	      tempres += log(2);
	    }

	    /*
	    if (mkr2->allele1 != mkr2->allele2) {
	      tempres += log(2);
	    }
	    if (mkr3->allele1 != mkr3->allele2) {
	      tempres += log(2);
	    }
	    */
#endif	    

	  }
	  else {
	    // Get the current genotypes of the parents
	    fmkr = markernumber(ind->father, ind->father->tmpint1);
	    mmkr = markernumber(ind->mother, ind->mother->tmpint1);
	    fmkr2= markernumber(ind2->father, ind2->father->tmpint1);
	    mmkr2= markernumber(ind2->mother, ind2->mother->tmpint1);
	    fmkr3= markernumber(ind3->father, ind3->father->tmpint1);
	    mmkr3= markernumber(ind3->mother, ind3->mother->tmpint1);
	    
#ifdef debugs
	    printf("Person: %s\n", ind->id);
	    printf("Father   : %d | %d\n", fmkr->allele1, fmkr->allele2);
	    printf("Father   : %d | %d\n", fmkr2->allele1, fmkr2->allele2);
	    printf("Father   : %d | %d\n", fmkr3->allele1, fmkr3->allele2);
	    printf("Mother   : %d | %d\n", mmkr->allele1, mmkr->allele2);
	    printf("Mother   : %d | %d\n", mmkr2->allele1, mmkr2->allele2);
	    printf("Mother   : %d | %d\n", mmkr3->allele1, mmkr3->allele2);
	    printf("Offspring: %d | %d\n", mkr->allele1, mkr->allele2);
	    printf("Offspring: %d | %d\n", mkr2->allele1, mkr2->allele2);
	    printf("Offspring: %d | %d\n", mkr3->allele1, mkr3->allele2);
#endif
	    
	    
	    // Check to see if combination legal
	    if ((fmkr->allele1 != mkr->allele1 && fmkr->allele2 != mkr->allele1) || 
		(mmkr->allele1 != mkr->allele2 && mmkr->allele2 != mkr->allele2) || 
		(fmkr2->allele1 != mkr2->allele1 && fmkr2->allele2 != mkr2->allele1) || 
		(mmkr2->allele1 != mkr2->allele2 && mmkr2->allele2 != mkr2->allele2) ||
		(fmkr3->allele1 != mkr3->allele1 && fmkr3->allele2 != mkr3->allele1) || 
		(mmkr3->allele1 != mkr3->allele2 && mmkr3->allele2 != mkr3->allele2)) {
	      tempres = 0;
	      break;
	    }

	    // If the father is heterozygous at either locus then multiply by .5
	    // (transmission probability)
	    fhettype = 1*(fmkr->allele1 != fmkr->allele2) + 10*(fmkr2->allele1 != fmkr2->allele2) + 100*(fmkr3->allele1 != fmkr3->allele2);
	    mhettype = 1*(mmkr->allele1 != mmkr->allele2) + 10*(mmkr2->allele1 != mmkr2->allele2) + 100*(mmkr3->allele1 != mmkr3->allele2);
	    
	    if (fhettype>0)
	      tempres += log(.5);
	    if (mhettype>0)
	      tempres += log(.5);

#ifdef reduceends		
	    if (ind->tmpint4==10 && mkr->allele1 != mkr->allele2) {
	      if (fhettype>0)
		tempres += log(.5);
	      if (mhettype>0)
		tempres += log(.5);
	    }
#endif

	    
	    //CLAUS: HERTIL
	    
	    switch (fhettype) {
	    case 0:
	    case 1:  // Het only at first marker
	    case 10: // Het only at second marker
	    case 100:// Het only at third marker
	      break; 
	    case 11: 
	      if ((mkr->allele1 == fmkr->allele1 && mkr2->allele1 == fmkr2->allele2) || (mkr->allele1 == fmkr->allele2 && mkr2->allele1 == fmkr2->allele1)) {
		// A recombination must have occured
		tempres +=log(gamma10 + gamma11);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(1-gamma10-gamma11);
		}
#endif
	      }
	      else {
		tempres +=log(1-gamma10-gamma11);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(gamma10 + gamma11);
		}
#endif
	      }

	    break;
	    case 101:
	      if ((mkr->allele1 == fmkr->allele1 && mkr3->allele1 == fmkr3->allele2) || (mkr->allele1 == fmkr->allele2 && mkr3->allele1 == fmkr3->allele1)) {
		// A recombination must have occured
		tempres +=log(gamma10 + gamma01);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(gamma01 + gamma10);
		}
#endif
	      }
	      else {
		tempres +=log(1-gamma10-gamma01);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(1-gamma10-gamma01);
		}
#endif

	      }

	      break;
	    case 110:
	      if ((mkr2->allele1 == fmkr2->allele1 && mkr3->allele1 == fmkr3->allele2) || (mkr2->allele1 == fmkr2->allele2 && mkr3->allele1 == fmkr3->allele1)) {
		// A recombination must have occured
		tempres +=log(gamma01 + gamma11);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(gamma01 + gamma11);
		}
#endif
	      }
	      else {
		tempres +=log(1-gamma01-gamma11);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(1-gamma01-gamma11);
		}
#endif

	      }

	      break;
	    case 111: // Father is triple heterozygous
	      if ((mkr->allele1 == fmkr->allele1 && mkr2->allele1 == fmkr2->allele1 && mkr3->allele1 == fmkr3->allele1) 
		  || (mkr->allele1 == fmkr->allele2 && mkr2->allele1 == fmkr2->allele2 && mkr3->allele1 == fmkr3->allele2)) {
		// No recombinations
		tempres += log(1-gamma11-gamma01-gamma10);

#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres += log(gamma10);
		}
#endif

	      }
	      else if ((mkr->allele1 == fmkr->allele1 && mkr2->allele1 == fmkr2->allele1 && mkr3->allele1 == fmkr3->allele2) 
		       || (mkr->allele1 == fmkr->allele2 && mkr2->allele1 == fmkr2->allele2 && mkr3->allele1 == fmkr3->allele1)) {
		// Recombination between 2 and 3
		tempres += log(gamma01);

#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres += log(gamma11);
		}
#endif

	      }
	      else if ((mkr->allele1 == fmkr->allele1 && mkr2->allele1 == fmkr2->allele2 && mkr3->allele1 == fmkr3->allele2) 
		       || (mkr->allele1 == fmkr->allele2 && mkr2->allele1 == fmkr2->allele1 && mkr3->allele1 == fmkr3->allele1)) {
		// Recombination between 1 and 2
		tempres += log(gamma10);		
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres += log(1-gamma11-gamma01-gamma10);
		}
#endif

	      }
	      else if ((mkr->allele1 == fmkr->allele1 && mkr2->allele1 == fmkr2->allele2 && mkr3->allele1 == fmkr3->allele1) 
		       || (mkr->allele1 == fmkr->allele2 && mkr2->allele1 == fmkr2->allele1 && mkr3->allele1 == fmkr3->allele2)) {
		// Recombination between 1 and 2 and 3
		tempres += log(gamma11);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres += log(gamma01);
		}
#endif

	      }
	      else {
		printf("ERROR: Programming error for three point linkage analysis. I should not be able to reach this point.\n");
	      }    
	      break;
	    }



	    switch (mhettype) {
	    case 0:
	    case 1:  // Het only at first marker
	    case 10: // Het only at second marker
	    case 100:// Het only at third marker
	      break; 
	    case 11: 
	      if ((mkr->allele2 == mmkr->allele1 && mkr2->allele2 == mmkr2->allele2) || (mkr->allele2 == mmkr->allele2 && mkr2->allele2 == mmkr2->allele1)) {
		// A recombination must have occured
		tempres +=log(gamma10 + gamma11);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(1-gamma10-gamma11);
		}
#endif
	      }
	      else {
		tempres +=log(1-gamma10-gamma11);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(gamma10 + gamma11);
		}
#endif
	      }


	    break;
	    case 101:
	      if ((mkr->allele2 == mmkr->allele1 && mkr3->allele2 == mmkr3->allele2) || (mkr->allele2 == mmkr->allele2 && mkr3->allele2 == mmkr3->allele1)) {
		// A recombination must have occured
		tempres +=log(gamma10 + gamma01);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(gamma01 + gamma10);
		}
#endif
	      }
	      else {
		tempres +=log(1-gamma10-gamma01);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(1-gamma10-gamma01);
		}
#endif
	      }
	      break;
	    case 110:
	      if ((mkr2->allele2 == mmkr2->allele1 && mkr3->allele2 == mmkr3->allele2) || (mkr2->allele2 == mmkr2->allele2 && mkr3->allele2 == mmkr3->allele1)) {
		// A recombination must have occured
		tempres +=log(gamma01 + gamma11);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(gamma01 + gamma11);
		}
#endif
	      }
	      else {
		tempres +=log(1-gamma01-gamma11);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres +=log(1-gamma01-gamma11);
		}
#endif
	      }

	      break;
	    case 111: // Father is triple heterozygous
	      if ((mkr->allele2 == mmkr->allele1 && mkr2->allele2 == mmkr2->allele1 && mkr3->allele2 == mmkr3->allele1) 
		  || (mkr->allele2 == mmkr->allele2 && mkr2->allele2 == mmkr2->allele2 && mkr3->allele2 == mmkr3->allele2)) {
		// No recombinations
		tempres += log(1-gamma11-gamma01-gamma10);
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres += log(gamma10);
		}
#endif

	      }
	      else if ((mkr->allele2 == mmkr->allele1 && mkr2->allele2 == mmkr2->allele1 && mkr3->allele2 == mmkr3->allele2) 
		       || (mkr->allele2 == mmkr->allele2 && mkr2->allele2 == mmkr2->allele2 && mkr3->allele2 == mmkr3->allele1)) {
		// Recombination between 2 and 3
		tempres += log(gamma01);	
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres += log(gamma11);
		}
#endif
	
	      }
	      else if ((mkr->allele2 == mmkr->allele1 && mkr2->allele2 == mmkr2->allele2 && mkr3->allele2 == mmkr3->allele2) 
		       || (mkr->allele2 == mmkr->allele2 && mkr2->allele2 == mmkr2->allele1 && mkr3->allele2 == mmkr3->allele1)) {
		// Recombination between 2 and 3
		tempres += log(gamma10);	
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres += log(1-gamma11-gamma01-gamma10);
		}
#endif
	
	      }
	      else if ((mkr->allele2 == mmkr->allele1 && mkr2->allele2 == mmkr2->allele2 && mkr3->allele2 == mmkr3->allele1) 
		       || (mkr->allele2 == mmkr->allele2 && mkr2->allele2 == mmkr2->allele1 && mkr3->allele2 == mmkr3->allele2)) {
		// Recombination between 2 and 3
		tempres += log(gamma11);			       
#ifdef reduceends		
		if (ind->tmpint4==10) {
		  tempres += log(gamma01);
		}
#endif

	      }
	      else {
		printf("ERROR: Programming error for three point linkage analysis. I should not be able to reach this point.\n");
	      }    
	      break;
	    }



	  }  // end of else

	  //	  printf("Tempres: %f\n", tempres);	  
	}

	res += exp(tempres);

	for (ind = indlist3; ind; ind = ind->next) {
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
	
      // Update the array for the second marker
      for (ind = indlist2; ind; ind = ind->next) {
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
    
    // Update the array for the first marker
    for (ind = indlist1; ind; ind = ind->next) {
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

  //  printf("RES    : %f\n", log(res));  
  return log(res);
}





double ThreepointLinkage2(IDlist *list1, IDlist *list2, IDlist *list3, int marker1, int marker2, int marker3, double gamma01, double gamma10, double gamma11)
{
  double res;

  IDlist *id1, *id2, *id3;

  res = 0;

  if (gamma01 < 0)
    return -1000000;
  if (gamma10 < 0)
    return -1000000;
  if (gamma11 < 0)
    return -1000000;

  if (gamma01 +gamma10 > .5)
    return -1000000;
  if (gamma10 +gamma11 > .5)
    return -1000000;
  if (gamma01 +gamma11 > .5)
    return -1000000;


  for (id1=list1, id2=list2, id3=list3; id1; id1=id1->next, id2=id2->next, id3=id3->next) {    
    res += LogThreepointLikelihood(id1->ind, id2->ind, id3->ind, marker1, marker2, marker3, gamma01, gamma10, gamma11);
  }




  return res;
}




double ThreepointLinkage(individual *indlist, int marker1, int marker2, int marker3, double gamma01, double gamma10, double gamma11)
{
  individual *ind, *ind2, *ind3, *ind4, *ind5;
  markerlist *marker, *mkr, *mkr2;
  namelist *Pedigrees, *ped;
  double res;


#ifdef reduceends
  int okcode, allok;
#endif

  res = 0;

  // Same marker?
  if (marker1 == marker2)
    return 0;
  if (marker1 == marker3)
    return 0;
  if (marker2 == marker3)
    return 0;

  if (gamma01 < 0)
    return -1000000;
  if (gamma10 < 0)
    return -1000000;
  if (gamma11 < 0)
    return -1000000;

  if (gamma01 +gamma10 > .5)
    return -1000000;
  if (gamma10 +gamma11 > .5)
    return -1000000;
  if (gamma01 +gamma11 > .5)
    return -1000000;

  Pedigrees = MakePedigreeList(indlist);  
  for (ped = Pedigrees; ped; ped = ped->next) {
    ind = GenotypeElimination(ped->name,marker1);
    ind2 = GenotypeElimination(ped->name,marker2);
    ind3 = GenotypeElimination(ped->name,marker3);

#ifdef debug    
    printf("Ped: %s\n", ped->name);
#endif

#ifdef lousyfix
    // Include special fixed for my case in order to finish my paper
    // Ideas:
    // 1) Homozygous founders need only have one and not both copies of the heterozygous genotypes (multiply the likelihood by two afterwards)
    for (ind4 = ind; ind4; ind4=ind4->next) {
      ind4->tmpint4=0;
      if (founder(ind4)) {
	ind4->tmpint2 = 0;
	// Go through genotype list
	for (marker = ind4->marker; marker; marker = marker->next) {
	  if (marker->allele1 == marker->allele2) {
	    continue;
	  }
	  for (mkr = marker->next; mkr; ) {
	    mkr2 = mkr->next;
	    if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
	      removelist(&ind4->marker, mkr);     
	      ind4->tmpint4 = 1;
	    }
	    mkr = mkr2;
	  }	  
	}
      }
#ifdef reduceends	
      else if (listlen(ind->offspring) == 0) {
	// It is an offspring/pedigree end
	
	// Check to see if _all_ possible heterogygous genotypes occur unordered
	allok = 1;
	for (marker = ind4->marker; marker; marker = marker->next) {	  
	  // Check to see if heterozygous
	  if (marker->allele1 != marker->allele2) {
	    okcode = 0;
	    for (mkr = marker->next; mkr; mkr=mkr->next) {	  
	      //	      if (marker
	      if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
		okcode = 1;
		break;		
	      }
	    }
	    allok *= okcode;
	  }
	}
	// If allok = 1 then all heterozygous genotypes occur in two copies
	// Then we can reduce them
	if (allok == 1) {
	  for (marker = ind4->marker; marker; marker = marker->next) {
	    if (marker->allele1 == marker->allele2) {
	      continue;
	    }
	    for (mkr = marker->next; mkr; ) {
	      mkr2 = mkr->next;
	      if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
		removelist(&ind4->marker, mkr);     
		ind4->tmpint4 = 10;
	      }
	      mkr = mkr2;
	    }	  
	  }
	}	
      }
#endif
    }

    // Go through marker 2
    // If a founder on the previous marker was only homoz. the reduce this one.
    for (ind4 = ind2; ind4; ind4=ind4->next) {
      if (founder(ind4)) {

	// Get a pointer to the same individual for marker 1
	ind5 = FindIndListMember(ind, ind4->id);

	// Was the person only homoz at the previous marker
	if (ind5->tmpint4 == 0) {

	  // Yes so reduce this one if possible

	  for (marker = ind4->marker; marker; marker = marker->next) {
	    if (marker->allele1 == marker->allele2) {
	      continue;
	    }
	    for (mkr = marker->next; mkr; ) {
	      mkr2 = mkr->next;
	      if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
		removelist(&ind4->marker, mkr);     
		ind5->tmpint4 = 2;
	      }
	      mkr = mkr2;
	    }	 
	  } 
	  
	}

      }
    }

    // Go through marker 3
    // If a founder on the previous marker was only homoz. the reduce this one.
    for (ind4 = ind3; ind4; ind4=ind4->next) {
      if (founder(ind4)) {

	// Get a pointer to the same individual for marker 1
	ind5 = FindIndListMember(ind, ind4->id);

	// Was the person only homoz at the previous marker
	if (ind5->tmpint4 == 0) {

	  // Yes so reduce this one if possible

	  for (marker = ind4->marker; marker; marker = marker->next) {
	    if (marker->allele1 == marker->allele2) {
	      continue;
	    }
	    for (mkr = marker->next; mkr; ) {
	      mkr2 = mkr->next;
	      if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
		removelist(&ind4->marker, mkr);     
		ind5->tmpint4 = 3;
	      }
	      mkr = mkr2;
	    }	 
	  } 
	  
	}

      }
    }
#endif

    res += LogThreepointLikelihood(ind, ind2, ind3, marker1, marker2, marker3, gamma01, gamma10, gamma11);
    FreeIndividualList(ind);
    FreeIndividualList(ind2);
    FreeIndividualList(ind3);
  }

  FreeNameList(Pedigrees);


  return res;
}

double ThreepointLinkageHaldane(individual *indlist, int marker1, int marker2, int marker3, double gamma01, double gamma11)
{
  individual *ind, *ind2, *ind3, *ind4, *ind5;
  markerlist *marker, *mkr, *mkr2;
  namelist *Pedigrees, *ped;
  double res, gamma10;

  res = 0;

  gamma10 = gamma11*(1-gamma11-gamma01)/(gamma01+gamma11);

  // Same marker?
  if (marker1 == marker2)
    return 0;
  if (marker1 == marker3)
    return 0;
  if (marker2 == marker3)
    return 0;

  if (gamma01 < 0)
    return -1000000;
  if (gamma10 < 0)
    return -1000000;
  if (gamma11 < 0)
    return -1000000;

  if (gamma01 +gamma10 > .5)
    return -1000000;
  if (gamma10 +gamma11 > .5)
    return -1000000;
  if (gamma01 +gamma11 > .5)
    return -1000000;


  Pedigrees = MakePedigreeList(indlist);  
  for (ped = Pedigrees; ped; ped = ped->next) {
    ind = GenotypeElimination(ped->name,marker1);
    ind2 = GenotypeElimination(ped->name,marker2);
    ind3 = GenotypeElimination(ped->name,marker3);


#ifdef lousyfix
    // Include special fixed for my case in order to finish my paper
    // Ideas:
    // 1) Homozygous founders need only have one and not both copies of the heterozygous genotypes (multiply the likelihood by two afterwards)
    for (ind4 = ind; ind4; ind4=ind4->next) {
      if (founder(ind4)) {
	ind4->tmpint2 = 0;
	// Go through genotype list
	for (marker = ind4->marker; marker; marker = marker->next) {
	  if (marker->allele1 == marker->allele2) {
	    continue;
	  }
	  for (mkr = marker->next; mkr; ) {
	    mkr2 = mkr->next;
	    if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
	      removelist(&ind4->marker, mkr);     
	      ind4->tmpint4 = 1;
	    }
	    mkr = mkr2;
	  }	  
	}
      }
    }

    // Go through marker 2
    // If a founder on the previous marker was only homoz. the reduce this one.
    for (ind4 = ind2; ind4; ind4=ind4->next) {
      if (founder(ind4)) {

	// Get a pointer to the same individual for marker 1
	ind5 = FindIndListMember(ind, ind4->id);

	// Was the person only homoz at the previous marker
	if (ind5->tmpint4 == 0) {

	  // Yes so reduce this one if possible

	  for (marker = ind4->marker; marker; marker = marker->next) {
	    if (marker->allele1 == marker->allele2) {
	      continue;
	    }
	    for (mkr = marker->next; mkr; ) {
	      mkr2 = mkr->next;
	      if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
		removelist(&ind4->marker, mkr);     
		ind5->tmpint4 = 2;
	      }
	      mkr = mkr2;
	    }	 
	  } 
	  
	}

      }
    }


    // Go through marker 3
    // If a founder on the previous marker was only homoz. the reduce this one.
    for (ind4 = ind3; ind4; ind4=ind4->next) {
      if (founder(ind4)) {

	// Get a pointer to the same individual for marker 1
	ind5 = FindIndListMember(ind, ind4->id);

	// Was the person only homoz at the previous marker
	if (ind5->tmpint4 == 0) {

	  // Yes so reduce this one if possible

	  for (marker = ind4->marker; marker; marker = marker->next) {
	    if (marker->allele1 == marker->allele2) {
	      continue;
	    }
	    for (mkr = marker->next; mkr; ) {
	      mkr2 = mkr->next;
	      if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
		removelist(&ind4->marker, mkr);     
		ind5->tmpint4 = 3;
	      }
	      mkr = mkr2;
	    }	 
	  } 
	  
	}

      }
    }

#endif

    res += LogThreepointLikelihood(ind, ind2, ind3, marker1, marker2, marker3, gamma01, gamma10, gamma11);
    FreeIndividualList(ind);
    FreeIndividualList(ind2);
    FreeIndividualList(ind3);
  }

  FreeNameList(Pedigrees);

  return res;
}


double ThreepointLinkageKosambi(individual *indlist, int marker1, int marker2, int marker3, double gamma01, double gamma11)
{
  individual *ind, *ind2, *ind3, *ind4, *ind5;
  markerlist *marker, *mkr, *mkr2;
  namelist *Pedigrees, *ped;
  double res, gamma10;

  res = 0;

  gamma10 = (-(gamma11+gamma01)*(gamma11+gamma01) +sqrt(SQR(gamma11*gamma11+gamma01*gamma01) + 2*gamma11*(gamma11+gamma01)) )/(2*(gamma11+gamma01));

  // Same marker?
  if (marker1 == marker2)
    return 0;
  if (marker1 == marker3)
    return 0;
  if (marker2 == marker3)
    return 0;

  if (gamma01 < 0)
    return -1000000;
  if (gamma10 < 0)
    return -1000000;
  if (gamma11 < 0)
    return -1000000;

  if (gamma01 +gamma10 > .5)
    return -1000000;
  if (gamma10 +gamma11 > .5)
    return -1000000;
  if (gamma01 +gamma11 > .5)
    return -1000000;


  Pedigrees = MakePedigreeList(indlist);  
  for (ped = Pedigrees; ped; ped = ped->next) {
    ind = GenotypeElimination(ped->name,marker1);
    ind2 = GenotypeElimination(ped->name,marker2);
    ind3 = GenotypeElimination(ped->name,marker3);


#ifdef lousyfix
    // Include special fixed for my case in order to finish my paper
    // Ideas:
    // 1) Homozygous founders need only have one and not both copies of the heterozygous genotypes (multiply the likelihood by two afterwards)
    for (ind4 = ind; ind4; ind4=ind4->next) {
      if (founder(ind4)) {
	ind4->tmpint2 = 0;
	// Go through genotype list
	for (marker = ind4->marker; marker; marker = marker->next) {
	  if (marker->allele1 == marker->allele2) {
	    continue;
	  }
	  for (mkr = marker->next; mkr; ) {
	    mkr2 = mkr->next;
	    if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
	      removelist(&ind4->marker, mkr);     
	      ind4->tmpint4 = 1;
	    }
	    mkr = mkr2;
	  }	  
	}
      }
    }

    // Go through marker 2
    // If a founder on the previous marker was only homoz. the reduce this one.
    for (ind4 = ind2; ind4; ind4=ind4->next) {
      if (founder(ind4)) {

	// Get a pointer to the same individual for marker 1
	ind5 = FindIndListMember(ind, ind4->id);

	// Was the person only homoz at the previous marker
	if (ind5->tmpint4 == 0) {

	  // Yes so reduce this one if possible

	  for (marker = ind4->marker; marker; marker = marker->next) {
	    if (marker->allele1 == marker->allele2) {
	      continue;
	    }
	    for (mkr = marker->next; mkr; ) {
	      mkr2 = mkr->next;
	      if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
		removelist(&ind4->marker, mkr);     
		ind5->tmpint4 = 2;
	      }
	      mkr = mkr2;
	    }	 
	  } 
	  
	}

      }
    }


    // Go through marker 3
    // If a founder on the previous marker was only homoz. the reduce this one.
    for (ind4 = ind3; ind4; ind4=ind4->next) {
      if (founder(ind4)) {

	// Get a pointer to the same individual for marker 1
	ind5 = FindIndListMember(ind, ind4->id);

	// Was the person only homoz at the previous marker
	if (ind5->tmpint4 == 0) {

	  // Yes so reduce this one if possible

	  for (marker = ind4->marker; marker; marker = marker->next) {
	    if (marker->allele1 == marker->allele2) {
	      continue;
	    }
	    for (mkr = marker->next; mkr; ) {
	      mkr2 = mkr->next;
	      if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
		removelist(&ind4->marker, mkr);     
		ind5->tmpint4 = 3;
	      }
	      mkr = mkr2;
	    }	 
	  } 
	  
	}

      }
    }

#endif

    res += LogThreepointLikelihood(ind, ind2, ind3, marker1, marker2, marker3, gamma01, gamma10, gamma11);
    FreeIndividualList(ind);
    FreeIndividualList(ind2);
    FreeIndividualList(ind3);
  }

  FreeNameList(Pedigrees);

  return res;
}



double wrapper(double x[]) {
  return (-ThreepointLinkage(INDLIST, MARKER1, MARKER2, MARKER3, x[0], x[1], x[2]));
}

double wrapperny(double x[]) {
  return (-ThreepointLinkage2(ID1, ID2, ID3, MARKER1, MARKER2, MARKER3, x[0], x[1], x[2]));
}


double wrapper2(double x[]) {
  return (-ThreepointLinkageHaldane(INDLIST, MARKER1, MARKER2, MARKER3, x[0], x[1]));
}

double wrapper3(double x[]) {
  return (-ThreepointLinkageKosambi(INDLIST, MARKER1, MARKER2, MARKER3, x[0], x[1]));
}


int PrepareDataForThreePointLinkageAnalysis(individual *indlist, int marker1, int marker2, int marker3)
{
  individual *ind, *ind2, *ind3, *ind4, *ind5;
  namelist *Pedigrees, *ped;
  markerlist *marker, *mkr, *mkr2;
  IDlist *list1, *list2, *list3;

  // Check whether we have specified the same marker twice
  if (marker1 == marker2)
    return -1;
  if (marker1 == marker3)
    return -1;
  if (marker2 == marker3)
    return -1;

  list1 = NULL;
  list2 = NULL;
  list3 = NULL;

  // Go through all pedigrees
  Pedigrees = MakePedigreeList(indlist);  
  for (ped = Pedigrees; ped; ped = ped->next) {
    ind = GenotypeElimination(ped->name,marker1);
    ind2 = GenotypeElimination(ped->name,marker2);
    ind3 = GenotypeElimination(ped->name,marker3);

#ifdef lousyfix
    // Include special fixed for my case in order to finish my paper
    // Ideas:
    // 1) Homozygous founders need only have one and not both copies of the heterozygous genotypes (multiply the likelihood by two afterwards)
    for (ind4 = ind; ind4; ind4=ind4->next) {
      if (founder(ind4)) {
	ind4->tmpint2 = 0;
	// Go through genotype list
	for (marker = ind4->marker; marker; marker = marker->next) {
	  if (marker->allele1 == marker->allele2) {
	    continue;
	  }
	  for (mkr = marker->next; mkr; ) {
	    mkr2 = mkr->next;
	    if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
	      removelist(&ind4->marker, mkr);     
	      ind4->tmpint4 = 1;
	    }
	    mkr = mkr2;
	  }	  
	}
      }
    }

    // Go through marker 2
    // If a founder on the previous marker was only homoz. the reduce this one.
    for (ind4 = ind2; ind4; ind4=ind4->next) {
      if (founder(ind4)) {

	// Get a pointer to the same individual for marker 1
	ind5 = FindIndListMember(ind, ind4->id);

	// Was the person only homoz at the previous marker
	if (ind5->tmpint4 == 0) {

	  // Yes so reduce this one if possible

	  for (marker = ind4->marker; marker; marker = marker->next) {
	    if (marker->allele1 == marker->allele2) {
	      continue;
	    }
	    for (mkr = marker->next; mkr; ) {
	      mkr2 = mkr->next;
	      if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
		removelist(&ind4->marker, mkr);     
		ind5->tmpint4 = 2;
	      }
	      mkr = mkr2;
	    }	 
	  } 
	  
	}

      }
    }


    // Go through marker 3
    // If a founder on the previous marker was only homoz. the reduce this one.
    for (ind4 = ind3; ind4; ind4=ind4->next) {
      if (founder(ind4)) {

	// Get a pointer to the same individual for marker 1
	ind5 = FindIndListMember(ind, ind4->id);

	// Was the person only homoz at the previous marker
	if (ind5->tmpint4 == 0) {

	  // Yes so reduce this one if possible

	  for (marker = ind4->marker; marker; marker = marker->next) {
	    if (marker->allele1 == marker->allele2) {
	      continue;
	    }
	    for (mkr = marker->next; mkr; ) {
	      mkr2 = mkr->next;
	      if (mkr->allele1 == marker->allele2 && mkr->allele2 == marker->allele1) {
		removelist(&ind4->marker, mkr);     
		ind5->tmpint4 = 3;
	      }
	      mkr = mkr2;
	    }	 
	  } 
	  
	}

      }
    }

#endif

    adduniqueidlist(&list1, ind);
    adduniqueidlist(&list2, ind2);
    adduniqueidlist(&list3, ind3);

  }
  FreeNameList(Pedigrees);

  ID1 = list1;
  ID2 = list2;
  ID3 = list3;

  return 0;
}

void CleanUpThreePointLinkageData(IDlist *list1, IDlist *list2, IDlist *list3) 
{
  IDlist *id;

  for (id = list1; id; id=id->next) {
    FreeIndividualList(id->ind);
  }
  for (id = list2; id; id=id->next) {
    FreeIndividualList(id->ind);
  }
  for (id = list3; id; id=id->next) {
    FreeIndividualList(id->ind);
  }

  freelist(list1);
  freelist(list2);
  freelist(list3);
}


double ThreepointLinkageAnalysis(individual *indlist, int marker1, int marker2, int marker3, double gamma01, double gamma10, double gamma11)
{
  int i, oldmap;
  double kosambi, haldane;
  
  // gamma_10, gamma_01, gamma_11
  double start[] = {.10 , .09 , .01 };
  double min, minh, mink;

  INDLIST = indlist;
  MARKER1 = marker1;
  MARKER2 = marker2;
  MARKER3 = marker3;

  if (PrepareDataForThreePointLinkageAnalysis(indlist, marker1, marker2, marker3) == -1) {
    printf("WARNING: Ending three point linkage analysis because less than three markers were used\n");
    return 0;
  }


  printf("===> %f\n", ThreepointLinkage(INDLIST, MARKER1, MARKER2, MARKER3, .1, .09, .01));

  min=simplex(wrapper,start,3,5.0e-3,.1);
  for (i=0;i<3;i++) {
    printf("%f\n",start[i]);
  }
  printf("Loglike (Full model): %f\n", min);

  double start2[] = {.1 , .01 };
  minh=simplex(wrapper2,start2,2,3.0e-3,.1);
  for (i=0;i<2;i++) {
    printf("%f\n",start2[i]);
  }
  printf("Loglike (Haldane   ): %f\n", minh);
  haldane = minh;

  double start3[] = {.1 , .01 };
  mink=simplex(wrapper3,start3,2,3.0e-3,.1);
  for (i=0;i<2;i++) {
    printf("%f\n",start3[i]);
  }
  printf("Loglike (Kosambi   ): %f\n", mink);

  kosambi = mink;

  if (haldane < kosambi) {
    printf("BEST MODEL: Haldane\n");
  }
  else 
    printf("BEST MODEL: Kosambi\n");


  oldmap = usedmap;
  usedmap = M_HALDANE;
  // Verify that the map is expanded

  if (MapFunction(start[0]+start[2]) + MapFunction(start[1]+start[2])>MapFunction(start[0]+start[1])) {
    // Need
    printf("MAPEXPAND Haldane   %f\n", .5*(pChi2(1, 2*(minh-min))));
    
  }
  else {
    printf("MAPEXPAND Haldane   1.000\n");
  }
  usedmap = M_KOSAMBI;
  // Verify that the map is expanded
  if (MapFunction(start[0]+start[2]) + MapFunction(start[1]+start[2])>MapFunction(start[0]+start[1])) {
    // Need
    printf("MAPEXPAND Kosambi   %f\n", .5*(pChi2(1,2*(mink-min))));
    
  }
  else {
    printf("MAPEXPAND Kosambi   1.000\n");
  }


  // Flush buffers
  fflush(0);

  CleanUpThreePointLinkageData(ID1, ID2, ID3);

  usedmap = oldmap;

  return min;
}





  // This function is stolen from somewhere

double rosen(double x[])
{
        return (100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])+(1.0-x[0])*(1.0-x[0]));
}


/* 
Program: rosen.c
Author : Michael F. Hutt
hutt@ieee.org
http://home.earthlink.net/~mfhutt
11/3/97

An implementation of the Nelder-Mead simplex method applied to
Rosenbrock's function.

Jan. 6, 1999 
Modified to conform to the algorithm presented
in Margaret H. Wright's paper on Direct Search Methods.

*/

double simplex(double (*func)(double[]), double start[],int n, double EPSILON,
                        double scale)
{

        int vs;        /* vertex with smallest value */
        int vh;        /* vertex with next smallest value */
        int vg;        /* vertex with largest value */

        int i,j,m,row;
        int k;          /* track the number of function evaluations */
        int itr;                /* track the number of iterations */

        double **v;          /* holds vertices of simplex */
        double pn,qn;        /* values used to create initial simplex */
        double *f;           /* value of function at each vertex */
        double fr;           /* value of function at reflection point */
        double fe;           /* value of function at expansion point */
        double fc;           /* value of function at contraction point */
        double *vr;         /* reflection - coordinates */
        double *ve;         /* expansion - coordinates */
        double *vc;         /* contraction - coordinates */
        double *vm;         /* centroid - coordinates */
        double min;

        double fsum,favg,s,cent;

        /* dynamically allocate arrays */

        /* allocate the rows of the arrays */
        v =  (double **) malloc ((n+1) * sizeof(double *));
        f =  (double *) malloc ((n+1) * sizeof(double));
        vr = (double *) malloc (n * sizeof(double));
        ve = (double *) malloc (n * sizeof(double));  
        vc = (double *) malloc (n * sizeof(double));  
        vm = (double *) malloc (n * sizeof(double));  

        /* allocate the columns of the arrays */
        for (i=0;i<=n;i++) {
                v[i] = (double *) malloc (n * sizeof(double));
        }


	/* create the initial simplex */
        /* assume one of the vertices is 0,0 */
           
        pn = scale*(sqrt(n+1)-1+n)/(n*sqrt(2));
        qn = scale*(sqrt(n+1)-1)/(n*sqrt(2));

        for (i=0;i<n;i++) {
	  v[0][i] = start[i];
        }
           
        for (i=1;i<=n;i++) {
	  for (j=0;j<n;j++) {
	    if (i-1 == j) {
	      v[i][j] = pn + start[j];
	    }
	    else {
	      v[i][j] = qn + start[j];
	    }
	  }
        }

	// New initialization
        for (i=0;i<n;i++) {
	  v[0][i] = 1.2*start[i];
        }

        for (i=1;i<=n;i++) {  // Which point
	  for (j=0;j<n;j++) {
	    if (i-1 == j) {
	      v[i][j] = 0.1*start[j];
	    }
	    else {
	      v[i][j] = 0.5*start[j];
	    }
	  }
        }

	
        /* find the initial function values */
        for (j=0;j<=n;j++) {
	  f[j] = func(v[j]);
        }
        
        k = n+1;

#ifdef debug       
        /* print out the initial values */
        printf("Initial Values\n");
        for (j=0;j<=n;j++) {
	  for (i=0;i<n;i++) {
	    printf("%f ",v[j][i]);
	  }
	  printf("%f\n",f[j]);
        }
#endif

        
/* begin the main loop of the minimization */
        for (itr=1;itr<=MAX_IT;itr++) {     
                /* find the index of the largest value */
                vg=0;
                for (j=0;j<=n;j++) {
                        if (f[j] > f[vg]) {
                                vg = j;
                                }
                }

                /* find the index of the smallest value */
                vs=0;
                for (j=0;j<=n;j++) {
                        if (f[j] < f[vs]) {
                                vs = j;
                                }
                }

                /* find the index of the second largest value */
                vh=vs;
                for (j=0;j<=n;j++) {
                        if (f[j] > f[vh] && f[j] < f[vg]) {
                                vh = j;
                        }
                }

                /* calculate the centroid */
                for (j=0;j<=n-1;j++) {
                        cent=0.0;
                        for (m=0;m<=n;m++) {
                                if (m!=vg) {
                                        cent += v[m][j];
                                }
                        }
                        vm[j] = cent/n;
                }

                /* reflect vg to new vertex vr */
                for (j=0;j<=n-1;j++) {
                        /*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
                        vr[j] = vm[j]+ALPHA*(vm[j]-v[vg][j]);
                }
		fr = func(vr);
                k++;
        
                if (fr < f[vh] && fr >= f[vs]) {
                        for (j=0;j<=n-1;j++) {
                                v[vg][j] = vr[j];
                        }
                        f[vg] = fr;
                }

                /* investigate a step further in this direction */
                if ( fr <  f[vs]) {
                        for (j=0;j<=n-1;j++) {
                                /*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
                                ve[j] = vm[j]+GAMMA*(vr[j]-vm[j]);
                        }
                        fe = func(ve);
                        k++;

                        /* by making fe < fr as opposed to fe < f[vs],                     
                           Rosenbrocks function takes 63 iterations as opposed 
                           to 64 when using double variables. */

                        if (fe < fr) {
                                for (j=0;j<=n-1;j++) {
                                        v[vg][j] = ve[j];
                                }
                                f[vg] = fe;
                        }
                        else {
                                for (j=0;j<=n-1;j++) {
                                        v[vg][j] = vr[j];
                                }
                                f[vg] = fr;
                        }
                }
        
                /* check to see if a contraction is necessary */
                if (fr >= f[vh]) {
                        if (fr < f[vg] && fr >= f[vh]) {
                        /* perform outside contraction */
                                for (j=0;j<=n-1;j++) {
                                        /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
                                        vc[j] = vm[j]+BETA*(vr[j]-vm[j]);
                                }
                                fc = func(vc);
                                k++;
                        }
                        else {
                        /* perform inside contraction */
                                for (j=0;j<=n-1;j++) {
                                        /*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
                                        vc[j] = vm[j]-BETA*(vm[j]-v[vg][j]);
                                }
                                fc = func(vc);
                                k++;
                        }


                        if (fc < f[vg]) {
                                for (j=0;j<=n-1;j++) {
                                        v[vg][j] = vc[j];
                                }
                                f[vg] = fc;
                        }
                        /* at this point the contraction is not successful,
                           we must halve the distance from vs to all the 
                           vertices of the simplex and then continue.
                           10/31/97 - modified to account for ALL vertices. 
                        */
                        else {
                                for (row=0;row<=n;row++) {
                                        if (row != vs) {
                                                for (j=0;j<=n-1;j++) {
                                                        v[row][j] = v[vs][j]+(v[row][j]-v[vs][j])/2.0;
                                                }
                                        }
                                }
                        f[vg] = func(v[vg]);
                        k++;
			f[vh] = func(v[vh]);
                        k++;

        
                        }
                }

#ifdef debug
                /* print out the value at each iteration */
                printf("Iteration %d\n",itr);
                for (j=0;j<=n;j++) {
                        printf("%f %f %f\n",v[j][0],v[j][1],f[j]);
                }
		fflush(0);
#endif

                /* test for convergence */
                fsum = 0.0;
                for (j=0;j<=n;j++) {
                        fsum += f[j];
                }
                favg = fsum/(n+1);
                s = 0.0;
                for (j=0;j<=n;j++) {
                        s += pow((f[j]-favg),2.0)/(n);
                }
                s = sqrt(s);
                if (s < EPSILON) break;
        }
/* end main loop of the minimization */

        /* find the index of the smallest value */
        vs=0;
        for (j=0;j<=n;j++) {
                if (f[j] < f[vs]) {
                        vs = j;
                        }
        }

#ifdef debug
        printf("The minimum was found at\n"); 
#endif

        for (j=0;j<n;j++) {
          start[j] = v[vs][j];
#ifdef debug
	  printf("%e\n",v[vs][j]);
#endif
        }


        min=func(v[vs]);
        k++;
#ifdef debug
        printf("%d Function Evaluations\n",k);
        printf("%d Iterations through program\n",itr);
#endif

        free(f);
        free(vr);
        free(ve);
        free(vc);
        free(vm);
        free(*v);
return min;
}

#ifdef debug
#undef debug
#endif

