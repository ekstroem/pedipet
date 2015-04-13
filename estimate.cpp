/*
 * PEDIPET
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


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "estimate.h"
#include "matrix.h"
#include "distr.h"

#define  strcmpl(s1,s2)        (strcasecmp(s1,s2))

#define TOL 0.00001
#define EDGE 0.0001
#define EPSILON 0.1
#define NEGINF  -10

#define multip
#define sibs

#define twostart

extern double allelefreq[MAXALLELES];
extern "C" int nChromosome;


//
// Calculates the kinship coeff matrix (\Phi)
//

MATRIX MMakeKinshipMatrix(individual *indlist)
{
  individual *ind, *ind2;
  int numpers, i, done, code;

  numpers = listlen(indlist);

  MATRIX res(numpers, numpers);

  // Sets tmpstr to 0 for all individual
  // 0 Means not set, 1 means set
  for(ind = indlist, i=1; ind; ind = ind->next, i++)
  {
    strcpy(ind->tmpstr1, "0");
    ind->localid = i;
  }

  // Should figure out the right order. A persons parents should
  // always appear first
  //
  // done is a boolean - are we finished?
  // i is the number to be placed
  // code is a booelan - did we place someone this round?
  done = 0;
  i = 1;
  while (!done)
  {
    // Noone has been placed
    code = 0;
    for(ind = indlist; ind; ind = ind->next)
    {
      // If founder and not previously placed
      if (!ind->father && !ind->mother)
      {
        if (atoi(ind->tmpstr1)==0)
        {
          sprintf(ind->tmpstr1, "%d",i);
          i++;
          code = 1;
	}
      }
      else if (atoip(ind->father->tmpstr1) && atoip(ind->mother->tmpstr1) && (atoip(ind->tmpstr1)==0))
      {
        // Non founder with both parents in data already
        sprintf(ind->tmpstr1, "%d",i);
        i++;
        code = 1;
      }
    }
    if (i==numpers+1)
    {
      done = 1;
      code = 1;
    }
    // If noone has been placed this turn
    if (!code)
    {
      printf("Error: Could not create kinship matrix.\nConsistency error in dataset\n");
      return res ;
    }
  }

  // Initialize
  for (i = 1; i<= numpers; i++) {
    for(ind = indlist; ind; ind = ind->next) {
      if (i==atoi(ind->tmpstr1))
        break;
    }
    // If founder
    if (!ind->father && !ind->mother)
      res(ind->localid, ind->localid) = 0.5;
    else
    {
      res(ind->localid, ind->localid) = 0.5 + 0.5*res(ind->father->localid, ind->mother->localid);
      // Should now fix the rest
      for (ind2 = indlist; ind2; ind2 = ind2->next)
      {
        if (atoi(ind2->tmpstr1)<i)
        {
          res(ind->localid, ind2->localid) = 0.5*res(ind2->localid, ind->father->localid)+0.5*res(ind2->localid, ind->mother->localid);
          res(ind2->localid, ind->localid) = res(ind->localid, ind2->localid);
        }
      }
    }
  }
  return res;
}




// Only works for non-inbred relatives
// Needs the kinship MATRIX as input
MATRIX MMakeDelta7Matrix(individual *indlist, const MATRIX kinship)
{

  individual *ind, *ind2;
  MATRIX Res(kinship.Rows(), kinship.Cols());
  int i, j, nObs;

  nObs = kinship.Cols();
  Res = 0; // Initializing the matrix

  // Use the local ID to fix matrix
  for (ind = indlist, i=1; ind; ind = ind->next, i++)
    ind->localid = i;

  for (ind = indlist, i=1; ind; ind = ind->next, i++)
  {
    Res(i,i) = 1;

    for (ind2 = ind->next, j=i+1; ind2; ind2 = ind2->next, j++)
    {
      if (!founder(ind) && !founder(ind2))
      {
	Res(i,j) = kinship(ind->father->localid, ind2->father->localid)*
	           kinship(ind->mother->localid, ind2->mother->localid)+
	           kinship(ind->father->localid, ind2->mother->localid)*
	           kinship(ind->mother->localid, ind2->father->localid);
	Res(j,i) = Res(i,j);
      }
    }
  }
  return Res;
}



MATRIX MakeIdentMatrix(int nObs)
{
  int i;

  MATRIX Res(nObs, nObs);
  Res = 0;
  for (i=1; i<=nObs; i++)
    Res(i,i) = 1;
  return Res;
}

MATRIX MakeRelationshipMatrix(individual *indlist)
{
  int i, j, code;
  individual *ind, *ind2;
  int nObs = listlen(indlist);
  char buf[60];

  MATRIX res(nObs, nObs);

  // Code checks if everything is OK
  code = 1;
  for (ind = indlist, i=1; ind; ind = ind->next, i++)
  {
    for (ind2 = ind, j=i; ind2; ind2 = ind2->next, j++)
    {
      res(i,j) = RelationType(ind,ind2);
      res(j,i) = res(i,j);

      if (res(i,j)<0)
      {
        sprintf(buf, "Relationship between %s and %s is not known. Please modify RelationType function\n", ind->id, ind2->id);
        printf("ERROR: %s",buf);
        WriteErrorMsg(buf);
	code = 0;
      }
    }
  }

  if (!code)
    printf("WARNING: Not all family relationships could be computed\n");
  
  return(res);
}


MATRIX MMakeXLinkedKinshipMatrix(individual *indlist) 
{
  individual *ind, *ind2;
  int numpers, i, done, code;

  numpers = listlen(indlist);

  MATRIX res(numpers, numpers);

  // Sets tmpstr to 0 for all individual
  // 0 Means not set, 1 means set
  for(ind = indlist, i=1; ind; ind = ind->next, i++)
  {
    strcpy(ind->tmpstr1, "0");
    ind->localid = i;
  }

  // Should figure out the right order. A persons parents should
  // always appear first
  //
  // done is a boolean - are we finished?
  // i is the number to be placed
  // code is a booelan - did we place someone this round?
  done = 0;
  i = 1;
  while (!done)
  {
    // Noone has been placed
    code = 0;
    for(ind = indlist; ind; ind = ind->next) {
      // If founder and not previously placed
      if (!ind->father && !ind->mother)
      {
        if (atoi(ind->tmpstr1)==0)
        {
          sprintf(ind->tmpstr1, "%d",i);
          i++;
          code = 1;
	}
      }
      else if (atoip(ind->father->tmpstr1) && atoip(ind->mother->tmpstr1) && (atoip(ind->tmpstr1)==0))
      {
        // Non founder with both parents in data already
        sprintf(ind->tmpstr1, "%d",i);
        i++;
        code = 1;
      }
    }
    if (i==numpers+1)
    {
      done = 1;
      code = 1;
    }
    // If noone has been placed this turn
    if (!code)
    {
      printf("Error: Could not create kinship matrix.\nConsistency error in dataset\n");
      return res ;
    }
  }

  // Initialize
  for (i = 1; i<= numpers; i++) {
    for(ind = indlist; ind; ind = ind->next) {
      if (i==atoi(ind->tmpstr1))
        break;
    }
    // If founder
    if (!ind->father && !ind->mother) {
      if (ind->sex == S_MALE)
	res(ind->localid, ind->localid) = 1;
      else
	res(ind->localid, ind->localid) = .5;
	}
    else {
      if (ind->sex == S_MALE) 
	res(ind->localid, ind->localid) = 1 + 0.5*res(ind->father->localid, ind->mother->localid);
      else {
	res(ind->localid, ind->localid) = 0.5 + 0.5*res(ind->father->localid, ind->mother->localid);
      }
      // Should now fix the rest
      for (ind2 = indlist; ind2; ind2 = ind2->next)
      {
        if (atoi(ind2->tmpstr1)<i) {
	  if (ind->sex == S_MALE) {
	    res(ind->localid, ind2->localid) = res(ind2->localid, ind->mother->localid);
	    res(ind2->localid, ind->localid) = res(ind->localid, ind2->localid);
	  }
	  else {
	    res(ind->localid, ind2->localid) = 0.5*res(ind2->localid, ind->father->localid)+0.5*res(ind2->localid, ind->mother->localid);
	    res(ind2->localid, ind->localid) = res(ind->localid, ind2->localid);
	  }
        }
      }
    }
  }
  return res;
}


MATRIX MMakeXLinkedKinshipMatrixBySex(individual *indlist, int sex1, int sex2) 
{
  individual *ind, *ind2;
  int i, j, nInd;
  MATRIX tempres = MMakeXLinkedKinshipMatrix(indlist);
  nInd = tempres.Rows();

  // Check for legal sex
  // Should never be true unless programming error occurs
  if (sex1>=MAXSEX || sex1 < 0 || sex2>=MAXSEX || sex2 < 0) {
    printf("ERROR: Trying to make X-linked kinship matrix for illegal sex");
    exit(1);
  }
    

  for (ind = indlist, i=1; ind; ind = ind->next, i++) {
    for (ind2 = ind, j=i; ind2; ind2 = ind2->next, j++) {
      // Delete the entry if not the right sex
      if ((ind->sex != sex1 || ind2->sex != sex2) && (ind->sex != sex2 || ind2->sex != sex1)) {
	tempres(i,j) = 0;
	tempres(j,i) = 0;
      }	
    }
  }
  return (tempres);
}



//
// SpecialFunction maximizes for 500 full sib pairs without parental info.
// ibd and k2 values must be specified
//
//
//

double SpecialFunction(individual *indlist, int traitnum, double position, double ibd[], double k2[])
{
  int nFamSize, nRealFamSize;
  int i, j, k;
  int nVC = 4;             // I, A, g, k
  int nCovar = 1;
  double full, no, add;

  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees;

  FILE *F;

  MATRIX *y, *x;
  MATRIX *VC;

  MATRIX Phi2;

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  //  printf("Made pedigree list (%d pedigrees)\n", nPedigrees);


  no = 0;

  // Prepares matrices for input to maximize
  y  = new MATRIX[nPedigrees];
  x  = new MATRIX[nPedigrees];
  VC = new MATRIX[nPedigrees*nVC]; // I + Phi2 + PI + K2
  int nSize[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++)
  {
    nFamSize = 0;

    ind2 = SelectPedigree(pedi->name, indlist);
    //    Phi2.resize(listlen(ind2), listlen(ind2));
    //    Phi2 = MMakeKinshipMatrix(ind2);
    nRealFamSize = listlen(ind2);
    FreeIndividualList(ind2);

    // Count the number of individuals
    for (ind = indlist; ind; ind = ind->next)
    {
      if (!strcmpl(ind->pedigree, pedi->name))
      {
	// From the correct pedigree
	
	// Check that no required covariates are missing
	if (!traitmiss(ind,traitnum) && !founder(ind))
	  nFamSize++;
      }
    }

    //    printf("Family %s has %d members\n", pedi->name, nFamSize);
    nSize[i] = nFamSize;
    y[i].Resize(nFamSize,1);
    x[i].Resize(nFamSize,nCovar);

    // Prepare the matrices
    for (j = 0; j< nVC; j++)
    {
      VC[i+j*nPedigrees].Resize(nFamSize, nFamSize);
      VC[i+j*nPedigrees] = 0.0;
    }

    j = 0;
    for (ind = indlist; ind; ind = ind->next)
    {
      if (!strcmpl(ind->pedigree, pedi->name))
      {
	// Check that no required covariates are missing
	if (!traitmiss(ind,traitnum) && !founder(ind))
	{

	  j++;
	  y[i](j,1) = trait(ind, traitnum);
	  x[i](j,1) = 1.0;
	  VC[i](j,j) = 1.0; // Residual error

	  for (k=1; k<=nFamSize; k++)
	  {
	    VC[i+nPedigrees](j,k) = .5;       // Residual additive effect
	    VC[i+2*nPedigrees](j,k) = ibd[i]; // Specific add. effect
	    VC[i+3*nPedigrees](j,k) = k2[i];  // Specific dom. effect
	  }
	  for (k = 0; k< nVC; k++)
	    VC[i+k*nPedigrees](j,j) = 1; // Diagonal
	}
      }
    }
  }

  MATRIX start(nVC,1);
  MATRIX beta(1, 1);

  beta = 1;

  start(1,1) = 1;
  start(2,1) = .1;
  start(3,1) = .1;
  start(4,1) = .1;

  //  start = 1;
  //  printf("Full model:\n");
  full = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC, VC, start, beta, METHOD_ML, 1, 0);


  if (position>0)
  {

    start(1,1) = 1;
    start(2,1) = .1;
    start(3,1) = .1;
    start(4,1) = .1;

    //  start = 1;
    printf("Full model:\n");
    full = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC, VC, start, beta, METHOD_ML, 1);

    start.Resize(3,1);
    start = .4;
    printf("Additive model (no k2):\n");
    add = MaximizeMixedModel(nPedigrees, nSize, y, x, start.Rows(), VC, start, beta, METHOD_ML, 1);

    start.Resize(2,1);
    start(1,1) = 1;
    start(2,1) = .2;
    printf("No effect- model (no pi and k2):\n");
    no = MaximizeMixedModel(nPedigrees, nSize, y, x, start.Rows(), VC, start, beta, METHOD_ML, 1);

    printf("No QTL LL is : %f \n", no);

  }


  delete[] y;
  delete[] x;
  delete[] VC;

  FreeNameList(Pedigrees);

  if (position>0)
  {
    F = fopen ("result.out","a");

    if (F == 0)
    {
      printf("Error in file result.out\n");
      exit(1);
    }

    fprintf(F, "%f %f %f %f %f %f %f\n", position, -full, -add, -1.0, -no, 33.0-position, pChi2(1,2*(full-no))*.5 + pChi2(2,2*(full-no))*.25);
    


    fclose(F);

  }

  return (full - no);
}


double SpecialFunction2(individual *indlist, int traitnum, double position, MATRIX ibd[], MATRIX k2[])
{
  int nFamSize, nRealFamSize;
  int i, j;
  int nVC = 5;             // I, A, D, g, k
  int nCovar = 1;
  double full, no, add, dom;

  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees;

  FILE *F;

  // General matrices for the families
  MATRIX GenI(9,9), GenA(9,9), GenD(9,9);
  
  GenI = MakeIdentMatrix(9);
  ind = SelectPedigree("F1", indlist);
  GenA = MMakeKinshipMatrix(ind)*2;
  FreeIndividualList(ind);
  GenD = GenI;
  GenD(5,6) = .25;
  GenD(6,5) = .25;
  GenD(7,8) = .25;
  GenD(8,7) = .25;


  MATRIX *y, *x;
  MATRIX *VC;

  MATRIX Phi2;

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  printf("Made pedigree list (%d pedigrees)\n", nPedigrees);

  // Prepares matrices for input to maximize
  y = new MATRIX[nPedigrees];
  x = new MATRIX[nPedigrees];
  VC = new MATRIX[nPedigrees*nVC]; // I + Phi2 + PI + K2
  int nSize[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++)
  {
    nFamSize = 0;

    ind2 = SelectPedigree(pedi->name, indlist);
    //    Phi2.resize(listlen(ind2), listlen(ind2));
    //    Phi2 = MMakeKinshipMatrix(ind2);
    nRealFamSize = listlen(ind2);
    FreeIndividualList(ind2);

    // Count the number of individuals
    for (ind = indlist; ind; ind = ind->next)
    {
      // From the correct pedigree
      if (!strcmpl(ind->pedigree, pedi->name))
      {	
	// Check that no required covariates are missing
	if (!traitmiss(ind,traitnum))
	  nFamSize++;
      }
    }

    //    printf("Family %s has %d members\n", pedi->name, nFamSize);
    nSize[i] = nFamSize;
    y[i].Resize(nFamSize,1);
    x[i].Resize(nFamSize,nCovar);

    // Prepare the matrices
    for (j = 0; j< nVC; j++)
    {
      VC[i+j*nPedigrees].Resize(nFamSize, nFamSize);
      VC[i+j*nPedigrees] = 0.0;
    }

    VC[i]              = GenI; // Residual error
    VC[i+1*nPedigrees] = GenA;
    VC[i+2*nPedigrees] = GenD;
    VC[i+3*nPedigrees] = ibd[i] + GenI;
    VC[i+4*nPedigrees] = k2[i]  + GenI;

    j = 0;
    for (ind = indlist; ind; ind = ind->next)
    {
      if (!strcmpl(ind->pedigree, pedi->name))
      {
	// Check that no required covariates are missing
	if (!traitmiss(ind,traitnum))
	{
	  j++;
	  y[i](j,1) = trait(ind, traitnum);
	  x[i](j,1) = 1.0;
	}
      }
    }
  }

    for (j = 0; j< nVC; j++)
    {
      //      VC[0+j*nPedigrees].Print();
      //      printf("------------------------------\n");
    }


  MATRIX start(nVC,1);
  MATRIX beta(1, 1);

  beta = 1;

  start(1,1) = 1;
  start(2,1) = .1;
  start(3,1) = .1;
  start(4,1) = .1;
  start(5,1) = .1;

  //  start = 1;
  printf("Full model:\n");
  full = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC, VC, start, beta, METHOD_ML, 1);
  
  if (position>0)
  {
    start.Resize(4,1);
    start(1,1) = 1;
    start = .1;
    printf("Additive model (no k2):\n");
    add = MaximizeMixedModel(nPedigrees, nSize, y, x, start.Rows(), VC, start, beta, METHOD_ML, 1);

    // Copies the dom effect into the add and estimates again

    for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++)
    {
      VC[i+3*nPedigrees] = VC[i+4*nPedigrees];
    }
    printf("Dominance model (no PI):\n");
    dom = MaximizeMixedModel(nPedigrees, nSize, y, x, start.Rows(), VC, start, beta, METHOD_ML, 1);

  
    start.Resize(3,1);
    start(1,1) = 1;
    start(2,1) = .001;
    start(3,1) = .001;

    printf("No effect- model (no pi and k2):\n");
    no = MaximizeMixedModel(nPedigrees, nSize, y, x, start.Rows(), VC, start, beta, METHOD_ML, 1);
  }

  delete[] y;
  delete[] x;
  delete[] VC;

  FreeNameList(Pedigrees);

  if (position>0)
  {
    F = fopen ("result.out","a");

    if (F == 0)
    {
      printf("Error in file result.out\n");
      exit(1);
    }

    fprintf(F, "%f %f %f %f %f %f %f\n", position, -full, -add, -dom, -no, 33.0-position, pChi2(1,2*(full-no))*.5 + pChi2(2,2*(full-no))*.25);
    


    fclose(F);

  }

  return (full);
}

//
//
// SingleMarker VC
//
//
double SingleMarkerVC(individual *indlist, int traitnum, int Dominance, MATRIX IBD[], MATRIX X[], int nCovar)
{
  int nFamSize;
  int i, j, k;
  int nVC; 
  double StartHyp, NullHyp, lldiff;

  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees;

  MATRIX *y, *x, *VC;
  MATRIX Phi2, Delta7, keep;

  // Figure out the number of variance components
  if (Dominance)
    nVC = 4;    // I 2*Phi Delta7 IBD
  else
    nVC = 3;    // I 2*PhiIBD

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  // Prepares matrices for input to maximize
  y = new MATRIX[nPedigrees];
  x = new MATRIX[nPedigrees];
  VC = new MATRIX[nPedigrees*nVC];

  int nSize[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    nFamSize = 0;

    // Select the current pedigree
    ind2 = SelectPedigree(pedi->name, indlist);
    int nPedSize = listlen(ind2);
    keep.Resize(nPedSize,1);
    
    // Figure out which individuals to keep
    // Only checks for response and NOT covariates
    keep = CompleteCases(ind2, traitnum);

    // The actual number of complet individuals
    int nFamSize = (int) keep.Sum();   


    Phi2.Resize(listlen(ind2), listlen(ind2));
    Phi2 = MMakeKinshipMatrix(ind2);
    if (Dominance)
      Delta7 = MMakeDelta7Matrix(indlist, Phi2);

    Phi2 = Phi2*2;

    FreeIndividualList(ind2);
    
    nSize[i] = nFamSize;
    y[i].Resize(nFamSize,1);
    x[i].Resize(nFamSize,nCovar);

    // Prepare the matrices
    for (j = 0; j< nVC; j++)
    {
      VC[i+j*nPedigrees].Resize(nFamSize, nFamSize);
      VC[i+j*nPedigrees] = 0.0;
    }

    VC[i+nPedigrees] = Phi2.Subset(keep,keep);
    if (Dominance) {
      VC[i+2*nPedigrees] = Delta7.Subset(keep,keep);
      VC[i+3*nPedigrees] = IBD[i].Subset(keep,keep);
    }
    else
      VC[i+2*nPedigrees] = IBD[i].Subset(keep,keep);

    x[i] = X[i].Subset(keep);

    j = 0;
    k = 0;
    for (ind = indlist; ind; ind = ind->next) {
      if (!strcmpl(ind->pedigree, pedi->name)) {
	k++;
	// Check that no required covariates are missing
	if (keep(k,1)) {
	  j++;
	  y[i](j,1) = trait(ind, traitnum);
	  VC[i](j,j) = 1.0; // Residual error
	}
      }
    }
  }

  freelist(Pedigrees);

  MATRIX start(nVC,1);
  MATRIX beta(nCovar, 1);

  start = 1;
  beta = 1;

  StartHyp = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC,  VC, start, beta, METHOD_ML, 1);
  start.Resize(nVC-1,1);
  start = 1;
  NullHyp = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC-1, VC, start, beta, METHOD_ML, 1);

  lldiff = -2*(NullHyp-StartHyp);
  printf("LR test statistic: %f  (p = %f)\n", lldiff, (lldiff<=0) ? 1 : 0.5*pChi2(1,lldiff));
  printf("LOD score        : %-5.2f\n", max(0,lldiff/(2*log(10))));

  delete[] y;
  delete[] x;
  delete[] VC;

  return 0;
}


/**********************************************

  Does multipoint analysis.
  Goes through the mibd files and 
			       
***********************************************/

double MultiMarkerVC(individual *indlist, int traitnum,  MATRIX IBD[], MATRIX X[], int nCovar, int Dominance, int Household)
{
  int nFamSize;
  int i, j, k, pos, bestpos;
  int nVC, currentVC=0, chromlength; 
  double StartHyp, NullHyp, lldiff, bestlod;
  FILE *fp;
  char buf[20];

  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees;

  MATRIX *y, *x, *VC, *keep;
  MATRIX Phi2, Delta7, HH;

  // Figure out the number of variance components
  nVC = 3;      // I 2*Phi IBD
  if (Dominance)
    nVC++;      // I 2*Phi Delta7 IBD
  if (Household)
    nVC++;      // I 2*Phi Delta7 HH IBD      - If not Dom then just remove Delta from this comment to get the order

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  // Prepares matrices for input to maximize
  y  = new MATRIX[nPedigrees];
  x  = new MATRIX[nPedigrees];
  keep = new MATRIX[nPedigrees];
  VC = new MATRIX[nPedigrees*nVC];

  int nSize[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    nFamSize = 0;

    // Select the current pedigree
    ind2 = SelectPedigree(pedi->name, indlist);

    int nPedSize = listlen(ind2);
    keep[i].Resize(nPedSize,1);

    
    // Figure out which individuals to keep
    // Only checks for response and NOT covariates
    keep[i] = CompleteCases(ind2, traitnum);

    // The actual number of complete individuals
    int nFamSize = (int) keep[i].Sum();

    Phi2.Resize(listlen(ind2), listlen(ind2));

    Phi2 = MMakeKinshipMatrix(ind2);

    if (Dominance)
      Delta7 = MMakeDelta7Matrix(indlist, Phi2);
    if (Household)
      HH = MakeIdentMatrix(nFamSize);

    Phi2 = Phi2*2;

    FreeIndividualList(ind2);
    
    nSize[i] = nFamSize;
    y[i].Resize(nFamSize,1);
    x[i].Resize(nFamSize,nCovar);

    // Prepare the matrices
    for (j = 0; j< nVC; j++) {
      VC[i+j*nPedigrees].Resize(nFamSize, nFamSize);
      VC[i+j*nPedigrees] = 0.0;
    }

    currentVC = 0;
    // Residual Variance, \sigma_I^2
    for (j = 0; j< VC[i].Rows(); j++) {
      VC[i](j+1,j+1) = 1;
    }
    currentVC++;

    // if (additive)
    VC[i+currentVC*nPedigrees] = Phi2.Subset(keep[i],keep[i]);
    currentVC++;

    if (Dominance) {
      VC[i+currentVC*nPedigrees] = Delta7.Subset(keep[i],keep[i]);
      currentVC++;
    }

    if (Household) {
      VC[i+currentVC*nPedigrees] = 1;
      currentVC++;
    }

    // if (IBD)
    //    VC[i+currentVC*nPedigrees] = IBD[i].Subset(keep[i],keep[i]);

    x[i] = X[i].Subset(keep[i]);

    j = 0;
    k = 0;
    for (ind = indlist; ind; ind = ind->next) {
      if (!strcmpl(ind->pedigree, pedi->name)) {
	k++;
	// Check that no required covariates are missing
	if (keep[i](k,1)) {
	  j++;
	  y[i](j,1) = trait(ind, traitnum);
	}
      }
    }
  }

  printf("Reading mibd\n");


  MATRIX start(nVC,1);
  MATRIX beta(nCovar, 1);



  // Should now go though the list of IBD files.
  // Add the VC's
  // Do analysis
  //


  // Calculate the null hypothesis
  start.Resize(nVC-1,1);
  start = 1;
  NullHyp = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC-1, VC, start, beta, METHOD_ML, 1, 0);
  start.Resize(nVC,1);
  start = 1;

  printf("Pos      LL         LOD\n");
  printf("--------------------------\n");

  // Go through all MIBD files in the right order
  chromlength = (int) (100*ChromosomeLength());

  bestpos = 0;
  bestlod = -1;
  // Go through all possible files at 1 cM 
  // WARNING: This should really be based on a sorted result of the directory listing
  for (pos = 0; pos <= chromlength; pos++) {
    sprintf(buf, "mibd.22.%d", pos);
    if ((fp = fopen(buf, "r"))) {
      fclose(fp);
      // Read the IBD file
      MATRIX *ibdres = ReadSolarIBDFile(indlist, buf, "pedindex.out");

      // Insert it in the dataset
      for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
	VC[i+currentVC*nPedigrees] = ibdres[i].Subset(keep[i],keep[i]);
      }
      delete[] ibdres;

      start = 1;
      beta = 1;
      StartHyp = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC,  VC, start, beta, METHOD_ML, 1, 0);
      lldiff = -2*(NullHyp-StartHyp);
      printf("%-3d   %5.3f     %6.4f\n", pos, StartHyp, max(0,lldiff/(2*log(10))));

      // Improved LOD score
      if (bestlod<max(0,lldiff/(2*log(10)))) {
	bestlod = max(0,lldiff/(2*log(10)));
	bestpos = pos;
      }
    }    
  }
  printf("\n    *** Highest LOD in pass ? was %6.4f at Chrom ?? Loc %d\n\n", bestlod, bestpos);

  freelist(Pedigrees);

  delete[] y;
  delete[] keep;
  delete[] x;
  delete[] VC;

  return 0;
}



//
// Calculate heritability for the individuals indexed by indlist
//

//
//
//         Dominance   : 0 exclude residual dominance, otherwise include
//

double Heritability(individual *indlist, int traitnum, int Dominance, MATRIX X[], int nCovar)
{
  int nFamSize, Household;
  int i, j, k, pos, bestpos;
  int nVC, currentVC, chromlength; 
  double StartHyp, NullHyp, lldiff, bestlod, NoDom;
  FILE *fp;
  char buf[20];

  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees;

  MATRIX *y, *x, *VC, *keep;
  MATRIX Phi2, Delta7, HH;

  Household = 0;

  //  printf("Here now\n");

  // Figure out the number of variance components
  nVC = 2;      // I 2*Phi 
  if (Dominance)
    nVC++;      // I 2*Phi Delta7 
  if (Household)
    nVC++;      // I 2*Phi Delta7 HH      - If not Dom then just remove Delta from this comment to get the order

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  // Prepares matrices for input to maximize
  y  = new MATRIX[nPedigrees];
  x  = new MATRIX[nPedigrees];
  keep = new MATRIX[nPedigrees];
  VC = new MATRIX[nPedigrees*nVC];

  int nSize[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    nFamSize = 0;

    // Select the current pedigree
    ind2 = SelectPedigree(pedi->name, indlist);

    int nPedSize = listlen(ind2);
    keep[i].Resize(nPedSize,1);
    
    // Figure out which individuals to keep
    // Only checks for response and NOT covariates
    keep[i] = CompleteCases(ind2, traitnum);

    // The actual number of complete individuals
    int nFamSize = (int) keep[i].Sum();

    Phi2.Resize(listlen(ind2), listlen(ind2));

    Phi2 = MMakeKinshipMatrix(ind2);

    if (Dominance) {
      Delta7.Resize(listlen(ind2), listlen(ind2));
      Delta7 = MMakeDelta7Matrix(ind2, Phi2);
    }
    if (Household)
      HH = MakeIdentMatrix(nFamSize);

    Phi2 = Phi2*2;

    FreeIndividualList(ind2);
    
    nSize[i] = nFamSize;
    y[i].Resize(nFamSize,1);
    x[i].Resize(nFamSize,nCovar);

    // Prepare the matrices
    for (j = 0; j< nVC; j++) {
      VC[i+j*nPedigrees].Resize(nFamSize, nFamSize);
      VC[i+j*nPedigrees] = 0.0;
    }

    currentVC = 0;
    // Residual Variance, \sigma_I^2
    for (j = 0; j< VC[i].Rows(); j++) {
      VC[i](j+1,j+1) = 1;
    }
    currentVC++;

    // if (additive)
    VC[i+currentVC*nPedigrees] = Phi2.Subset(keep[i],keep[i]);
    currentVC++;

    if (Dominance) {
      VC[i+currentVC*nPedigrees] = Delta7.Subset(keep[i],keep[i]);
      currentVC++;
    }

    if (Household) {
      VC[i+currentVC*nPedigrees] = 1;
      currentVC++;
    }

    // if (IBD)
    //    VC[i+currentVC*nPedigrees] = IBD[i].Subset(keep[i],keep[i]);

    x[i] = X[i].Subset(keep[i]);

    j = 0;
    k = 0;
    for (ind = indlist; ind; ind = ind->next) {
      if (!strcmpl(ind->pedigree, pedi->name)) {
	k++;
	// Check that no required covariates are missing
	if (keep[i](k,1)) {
	  j++;
	  y[i](j,1) = trait(ind, traitnum);
	}
      }
    }
  }


  MATRIX start(nVC,1);
  MATRIX beta(nCovar, 1);

  printf("Polygenic model:\n");
  start.Resize(nVC,1);
  start = .02;

  start(1,1) = 1;

  NullHyp = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC, VC, start, beta, METHOD_ML, 1, 1);
  
  /*
  double totalvar = start.Sum();
  printf("Narrow sense heritability: %5.3f\n", start(2,1)/totalvar);
  if (Dominance) {
    printf("Broad sense heritability : %5.3f\n", (start(2,1)+start(3,1))/totalvar);
  }
  */

  freelist(Pedigrees);

  delete[] y;
  delete[] keep;
  delete[] x;
  delete[] VC;

  return 0;
}

// New Genotype for full sibs
int NewGenoMin(int  gt1, int gt2)
{
  switch (gt1*10+gt2)
  {
  case 0: return 2;
  case 1: return 1;
  case 2: return 0;
  case 10: return 1;
  case 11: return 0;
  case 12: return 1;
  case 20: return 0;
  case 21: return 1;
  case 22: return 2;
  }
  return 0;  
}

// New Genotype for full sibs
int NewGenoMax(int  gt1, int gt2)
{
  switch (gt1*10+gt2)
  {
  case 0: return 2;
  case 1: return 1;
  case 2: return 0;
  case 10: return 1;
  case 11: return 2;
  case 12: return 1;
  case 20: return 0;
  case 21: return 1;
  case 22: return 2;
  }
  return 0;  
}

//
//  Model: 1 - full, 2 - add, 3 - dom, 4 - noeffect
//  1: I, A, D, s1, s2
//  2: I, A, D, s1 = s2/2
//  3: I, A, D, s1 = s2

double CalculateMixLL(individual *indlist, MATRIX RELAP[], MATRIX resp[], MATRIX param, int model)
{
  int nodom = 0;           // No residual dominance
  int nFamilyMembers = 8;  
  int nPedigrees;
  // Loop variables
  //  int i, i29, i39, i56, i57, i58, i59, i67, i68, i69, i78, i79;

  int i, i18, i28, i45, i46, i47, i48, i56, i57, i58, i67, i68;
  int min58, max58, min68, max68, min45, min56, min46, max45, max56, max46;

  double res, dLogDet, mean;
  double pprod, famres, meanvalue;

  MATRIX OrigOmega(nFamilyMembers,nFamilyMembers), Omega(nFamilyMembers,nFamilyMembers);
  MATRIX InvOmega(nFamilyMembers,nFamilyMembers);
  MATRIX Delta7(nFamilyMembers,nFamilyMembers), GenA(nFamilyMembers,nFamilyMembers);
  MATRIX GenD(nFamilyMembers,nFamilyMembers), GenI(nFamilyMembers,nFamilyMembers);
  MATRIX GenG(nFamilyMembers,nFamilyMembers);
  MATRIX y(nFamilyMembers,1), yT(1,nFamilyMembers), miniparam(3,1);
  namelist *Pedigrees;
  individual *ind;

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  GenI = MakeIdentMatrix(nFamilyMembers);
  ind = SelectPedigree("F1", indlist);
  GenA = MMakeKinshipMatrix(ind)*2;            // Res. Additive
  GenD = MMakeDelta7Matrix(ind, GenA*.5);      // Res. Dominance
  FreeIndividualList(ind);

  // Ok for the residual VC's

  GenG = 0;
  GenG(1, 4) = 1;
  GenG(4, 1) = 1;
  GenG(1, 5) = 1;
  GenG(5, 1) = 1;
  GenG(1, 6) = 1;
  GenG(6, 1) = 1;
  GenG(1, 7) = 1;
  GenG(7, 1) = 1;


  GenG(2, 4) = 1;
  GenG(4, 2) = 1;
  GenG(2, 5) = 1;
  GenG(5, 2) = 1;
  GenG(2, 6) = 1;
  GenG(6, 2) = 1;
  GenG(2, 7) = 1;
  GenG(7, 2) = 1;

  GenG(3, 8) = 1;
  GenG(8, 3) = 1;
  GenG(7, 8) = 1;
  GenG(8, 7) = 1;

  // 5 parametre:  I, A, D, 1, 2

  mean = param(1,1);
  miniparam = 0;

  if (!nodom) // No residual dominance
    {
      switch (model) {
      case 1: miniparam(2,1) = 0.5*param(5,1);
	miniparam(3,1) = param(5,1) + param(6,1); break;
      case 2: miniparam(2,1) = 0.5*param(5,1);
	miniparam(3,1) = param(5,1); break;
      case 4: miniparam = 0; break;
      }
      OrigOmega = GenI*(param(2,1) + miniparam(3,1)) + GenA*param(3,1) + 
	GenD*param(4,1) + GenG*miniparam(2,1); // I
    }
  else
    {
      switch (model) {
      case 1: miniparam(2,1) = 0.5*param(4,1);
	miniparam(3,1) = param(4,1) + param(5,1); break;
      case 2: miniparam(2,1) = 0.5*param(4,1);
	miniparam(3,1) = param(4,1); break;
      case 4: miniparam = 0; break;
      }

      OrigOmega = GenI*(param(2,1) + miniparam(3,1)) + GenA*param(3,1) + GenG*miniparam(2,1); // I
    }


  res = 0;
  meanvalue = 0;

  for (i=0; i<nPedigrees; i++)
  {
    y = resp[i]-mean;
    yT = Transpose(y);

    // Adds I and 2 in the diagonal

    res -= 0.5*log(2*PI)*nFamilyMembers;

    famres = 0;

  if (model < 4)
  {

    // Connects to individual 7
  for (i47 = 0; i47<=2; i47++)
  {
    if (RELAP[i](5, i47+1)==0) 
      continue;

  for (i57 = 0; i57<=2; i57++)
  {
    if (RELAP[i](8, i57+1)==0) 
      continue;

  for (i67 = 0; i67<=2; i67++)
  {
    if (RELAP[i](10, i67+1)==0) 
      continue;

  for (i18 = 0; i18<=1; i18++)
  {
    if (RELAP[i](1, i18+1)==0) 
      continue;

    i28 = 1 - i18;

    min56 = NewGenoMin(i57,i67);
    max56 = NewGenoMax(i57,i67);

    min45 = NewGenoMin(i47,i57);
    max45 = NewGenoMax(i47,i57);

  for (i56 = min56 ; i56<= max56; i56 +=2)
  {

  for (i45 = min45; i45<=max45; i45 +=2)
  {

    // Now prepares 46
    min46 = max(NewGenoMin(i45,i56), NewGenoMin(i47,i67));
    max46 = min(NewGenoMax(i45,i56), NewGenoMax(i47,i67));

    if (min46>max46)
      printf("OOps - error   %d  %d  [%d %d   %d %d]\n", min46, max46, i45, i56, i47, i67);

  for (i46 = min46; i46<=max46; i46 +=2)
  {

  // 2 => 1, 0 and 1=> 0
  for (i48 = max(i47-1,0); i48 <= min(i47,1); i48++)
  {

    // Figure out the min and max for i58
    switch(i57)
    {
    case 0: min58 = 0; max58 = 0; break;
    case 2: min58 = 1; max58 = 1; break;
    case 1: switch(i45) 
            {
              case 0: min58 = max58 = 1-i48; break;
              case 2: min58 = max58 = i48; break;
	      case 1: min58 = 0; max58 = 1; break;
	    }
    }

  for (i58 = min58 ; i58 <= max58; i58++)
  {
    //  {


    // Figure out the min and max for i68
    switch(i67)
    {
    case 0: min68 = 0; max68 = 0; break;
    case 2: min68 = 1; max68 = 1; break;
    case 1: switch(i46) 
            {
              case 0: min68 = max68 = 1-i48; break;
              case 2: min68 = max68 = i48; break;
	      case 1: switch(i56) 
		      {
		      case 0: min68 = max68 = 1-i58; break;
		      case 2: min68 = max68 = i58; break;
		      case 1: min68 = 0; max68 = 1; break;
		      }
	    }
    }

  for (i68 = min68 ; i68 <= max68; i68++)
  {

    pprod = 0;

    pprod += log(RELAP[i]( 5,i47+1));
    pprod += log(RELAP[i]( 8,i57+1));
    pprod += log(RELAP[i](10,i67+1));
    pprod += log(RELAP[i]( 1,i18+1));

    if (min45 < max45)
    {
      pprod += log(.5);
    }

    if (min46 < max46)
    {
      pprod += log(.5);
    }

    if (min56 < max56)
    {
      pprod += log(.5);
    }

    // i48
    if (max(i47-1,0) < min(i47,1))
      pprod += log(.5);

    // i58
    if (min58 < max58)
      pprod += log(.5);

    // i68
    if (min68 < max68)
      pprod += log(.5);

    Omega = OrigOmega;

    Omega(1,8) += miniparam(i18+1,1) ;
    Omega(8,1) = Omega(1,8);
    Omega(2,8) += miniparam(i28+1,1) ;
    Omega(8,2) = Omega(2,8);
    Omega(4,5) += miniparam(i45+1,1) ;
    Omega(5,4) = Omega(4,5);
    Omega(4,6) += miniparam(i46+1,1) ;
    Omega(6,4) = Omega(4,6);
    Omega(4,7) += miniparam(i47+1,1) ;
    Omega(7,4) = Omega(4,7);
    Omega(4,8) += miniparam(i48+1,1) ;
    Omega(8,4) = Omega(4,8);


    Omega(5,6) += miniparam(i56+1,1) ;
    Omega(6,5) = Omega(5,6);
    Omega(5,7) += miniparam(i57+1,1) ;
    Omega(7,5) = Omega(5,7);
    Omega(5,8) += miniparam(i58+1,1) ;
    Omega(8,5) = Omega(5,8);


    Omega(6,7) += miniparam(i67+1,1) ;
    Omega(7,6) = Omega(6,7);
    Omega(6,8) += miniparam(i68+1,1) ;
    Omega(8,6) = Omega(6,8);

    /*    printf("---------------------\n");
    printf("---------------------\n");
    printf("%d %d %d %d %d %d %d %d %d %d %d\n", i18, i28, i45, i46, i47, i48, i56, i57, i58, i67, i68);
    printf("---------------------\n");
    Omega.Print();
    printf("---------------------\n");
    */

    InvOmega = Inverse(Omega, &dLogDet);
    famres += exp(pprod - 0.5*dLogDet -0.5*(yT*InvOmega*y)(1,1));

    //    famres += exp(pprod);
    //    printf("Famres (%d): %f\n", i, famres);

  }
  } 
  } 
  }
  } 
  } 
  } 
  }
  } 
  } 
  }
  else  // Model == 4
    {
    Omega = OrigOmega;
    //  Omega.Print();
    InvOmega = Inverse(Omega, &dLogDet);
    famres += exp(- 0.5*dLogDet -0.5*(yT*InvOmega*y)(1,1));
    }
    
    res += log(famres);
  }

  return (-res);

}


double CalculateMixLLNucl(individual *indlist, MATRIX RELAP[], MATRIX resp[], MATRIX param, int model)
{
  int nodom = 1;           // No residual dominance
  int nFamilyMembers = 6;  
  int nPedigrees;
  // Loop variables
  //  int i, i29, i39, i56, i57, i58, i59, i67, i68, i69, i78, i79;

  int i, i46, i56, min45, max45, i45;
  int i34, i35, i36, min34, max34, min35,max35;

  double res, dLogDet, mean;
  double pprod, famres, meanvalue;

  MATRIX OrigOmega(nFamilyMembers,nFamilyMembers), Omega(nFamilyMembers,nFamilyMembers);
  MATRIX InvOmega(nFamilyMembers,nFamilyMembers);
  MATRIX Delta7(nFamilyMembers,nFamilyMembers), GenA(nFamilyMembers,nFamilyMembers);
  MATRIX GenD(nFamilyMembers,nFamilyMembers), GenI(nFamilyMembers,nFamilyMembers);
  MATRIX GenG(nFamilyMembers,nFamilyMembers);
  MATRIX y(nFamilyMembers,1), yT(1,nFamilyMembers), miniparam(3,1);
  namelist *Pedigrees;
  individual *ind;

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  GenI = MakeIdentMatrix(nFamilyMembers);
  ind = SelectPedigree("F1", indlist);
  GenA = MMakeKinshipMatrix(ind)*2;            // Res. Additive
  GenD = MMakeDelta7Matrix(ind, GenA*.5);      // Res. Dominance
  FreeIndividualList(ind);

  // Ok for the residual VC's

  GenG = 0;
  GenG(1, 3) = 1;
  GenG(3, 1) = 1;
  GenG(1, 4) = 1;
  GenG(4, 1) = 1;
  GenG(1, 5) = 1;
  GenG(5, 1) = 1;
  GenG(1, 6) = 1;
  GenG(6, 1) = 1;

  GenG(2, 3) = 1;
  GenG(3, 2) = 1;
  GenG(2, 4) = 1;
  GenG(4, 2) = 1;
  GenG(2, 5) = 1;
  GenG(5, 2) = 1;
  GenG(2, 6) = 1;
  GenG(6, 2) = 1;

  // 5 parametre:  I, A, D, 1, 2

  mean = param(1,1);
  miniparam = 0;

  if (!nodom) // No residual dominance
    {
      switch (model) {
      case 1: miniparam(2,1) = 0.5*param(5,1);
	miniparam(3,1) = param(5,1) + param(6,1); break;
      case 2: miniparam(2,1) = 0.5*param(5,1);
	miniparam(3,1) = param(5,1); break;
      case 4: miniparam = 0; break;
      }
      OrigOmega = GenI*(param(2,1) + miniparam(3,1)) + GenA*param(3,1) + 
	GenD*param(4,1) + GenG*miniparam(2,1); // I
    }
  else
    {
      switch (model) {
      case 1: miniparam(2,1) = 0.5*param(4,1);
	miniparam(3,1) = param(4,1) + param(5,1); break;
      case 2: miniparam(2,1) = 0.5*param(4,1);
	miniparam(3,1) = param(4,1); break;
      case 4: miniparam = 0; break;
      }

      OrigOmega = GenI*(param(2,1) + miniparam(3,1)) + GenA*param(3,1) + GenG*miniparam(2,1); // I
    }


  res = 0;
  meanvalue = 0;

  for (i=0; i<nPedigrees; i++)
  {
    y = resp[i]-mean;
    yT = Transpose(y);

    // Adds I and 2 in the diagonal

    res -= 0.5*log(2*PI)*nFamilyMembers;

    famres = 0;

  if (model < 4)
  {

    // Connects to individual 6
  for (i36 = 0; i36<=2; i36++)
  {
    if (RELAP[i](3, i36+1)==0) 
      continue;

  for (i46 = 0; i46<=2; i46++)
  {
    if (RELAP[i](5, i46+1)==0) 
      continue;

  for (i56 = 0; i56<=2; i56++)
  {
    if (RELAP[i](6, i56+1)==0) 
      continue;

    min34 = NewGenoMin(i36,i46);
    max34 = NewGenoMax(i36,i46);

    min45 = NewGenoMin(i46,i56);
    max45 = NewGenoMax(i46,i56);

  for (i34 = min34 ; i34<= max34; i34 +=2)
  {
    if (RELAP[i](1, i34+1)==0) 
      continue;

  for (i45 = min45; i45<=max45; i45 +=2)
  {
    if (RELAP[i](4, i45+1)==0) 
      continue;

    // Now prepares 35
    min35 = max(NewGenoMin(i36,i56), NewGenoMin(i34,i45));
    max35 = min(NewGenoMax(i36,i56), NewGenoMax(i34,i45));

    if (min35>max35)
      printf("OOps - error   %d  %d  [%d %d   %d %d]\n", min35, max35, i36, i56, i45, i34);

  for (i35 = min35; i35<=max35; i35 +=2)
  {
    //    if (RELAP[i](2, i35+1)==0) 
    //      continue;

    pprod = 0;

    pprod += log(RELAP[i]( 3,i36+1));
    pprod += log(RELAP[i]( 5,i46+1));
    pprod += log(RELAP[i]( 6,i56+1));

    if (min34 < max34)
    {
      pprod += log( RELAP[i](1,i34+1)  / (RELAP[i](1,1) + RELAP[i](1,3)) );
    }
    
    if (min45 < max45)
      pprod += log( RELAP[i](4,i45+1)  / (RELAP[i](4,1) + RELAP[i](4,3)) );

    //    if (min35 < max35) {
    //      pprod += log( RELAP[i](2,i35+1)  / (RELAP[i](2,1) + RELAP[i](2,3)) );
    //    }

    Omega = OrigOmega;

    Omega(3,4) += miniparam(i34+1,1) ;
    Omega(4,3) = Omega(3,4);
    Omega(3,5) += miniparam(i35+1,1) ;
    Omega(5,3) = Omega(3,5);
    Omega(3,6) += miniparam(i36+1,1) ;
    Omega(6,3) = Omega(3,6);
    Omega(4,5) += miniparam(i45+1,1) ;
    Omega(5,4) = Omega(4,5);
    Omega(4,6) += miniparam(i46+1,1) ;
    Omega(6,4) = Omega(4,6);
    Omega(5,6) += miniparam(i56+1,1) ;
    Omega(6,5) = Omega(5,6);

    InvOmega = Inverse(Omega, &dLogDet);
    famres += exp(pprod - 0.5*dLogDet -0.5*(yT*InvOmega*y)(1,1));

    //    famres += exp(pprod);
    //    printf("Famres (%d): %15f    %d %d %d %d %d %d\n", i, famres*100000, i34, i35, i36, i45, i46, i56);

  }
  } 
  } 
  }
  } 
  }
  } 
  else  // Model == 4
    {
    Omega = OrigOmega;
    //  Omega.Print();
    InvOmega = Inverse(Omega, &dLogDet);
    famres += exp(- 0.5*dLogDet -0.5*(yT*InvOmega*y)(1,1));
    }
    
    res += log(famres);
  }

  return (-res);

}


MATRIX *SpecialDeriv(MATRIX &param, individual *indlist, MATRIX RELAP[], MATRIX resp[], int model)
{
  int i, j, nParam;
  double h, savethis;

  //  printf("Parameters passed to Deriv function:\n");
  //  param.Print();
  //  printf("--- End of parameters\n");

  nParam = param.Rows();

  MATRIX newh(nParam, 1), dParam (nParam, 1), dSecOrd(nParam, nParam), xplus(nParam,1), yplus(nParam,1);
  MATRIX *res;

  res = new MATRIX[3](nParam, nParam);

  res[0].Resize(nParam,1);
  res[2].Resize(1,1);

  h = 0.0001;



  savethis = CalculateMixLLNucl(indlist, RELAP, resp, param, model);
  //  printf("Savethis : %f\n", savethis);

  for (i = 0; i < nParam; i++)
  {
    newh = 0;
    newh(i+1,1) = h;
    dParam(i+1,1) = CalculateMixLLNucl(indlist, RELAP, resp, param+newh, model) - savethis;

    //    printf("  Nr %d : %f\n", i+1, dParam(i+1,1)+savethis);

    for (j = 0; j<= i; j++)
    {
      if (i==j)
      {
	dSecOrd(i+1,j+1) =  (CalculateMixLLNucl(indlist, RELAP, resp, param+(newh*2), model) - 2*dParam(i+1,1) - savethis)/(h*h);
      }
      else
      {
	yplus = 0;
	yplus(j+1,1) = h;
	
	dSecOrd(i+1,j+1) =  (CalculateMixLLNucl(indlist, RELAP, resp, param+newh+yplus, model) - 
			     dParam(i+1,1) - dParam(j+1,1) - savethis)/(h*h);
	dSecOrd(j+1,i+1) =  dSecOrd(i+1,j+1);
      }
    } 

  }

  dParam *= 1.0/h;

  res[0] = dParam;
  res[1] = dSecOrd;
  res[2](1,1) = savethis;

  return res;
}

double SpecialFunction3(individual *indlist, int traitnum, int position, MATRIX RELAP[])
{
  //  int nVC = 5;             // I, A, D, g, k
  int i, j, iter;
  int model = position;
  int founders = 1; // INclude founders
  // Loop variables

  int nodom = 1;
  int nParam = 6-nodom;  // Mean, I, A, D; g, k
  int maxIter = 40;
  int code = 0;
  double downscale, step;
  double savethis;
  double oldfunctionvalue, functionvalue;
  double claus;

  int nPedSize = 6;
  int nPedigrees = listlen(indlist)/nPedSize;

  

  //  int nPedSize = 2;
  //  int nPedigrees = listlen(indlist)/4;

  individual *ind;

  MATRIX *resp, *deriv;
  MATRIX param(nParam, 1);

  resp = new MATRIX[nPedigrees](nPedSize,1);

  for (ind = indlist, i=0, j=1; ind; ind = ind->next)
  {
    if (founders || !founder(ind))
      {
    resp[i](j,1) = trait(ind,traitnum);
    j++;
    if (j>resp[i].Rows()) 
    {
      i++;
      j=1;
    }
      }
  }

  switch (model) {
  case 2:  // No dominance effect
    nParam -= 1;
    param.Resize(nParam,1);   
    break;
  case 4: // No locus effect
    nParam -= 2;
    param.Resize(nParam,1);   
    break;
  }

  MATRIX newh(nParam, 1);
  MATRIX dParam(nParam, 1),  SecOrd(nParam, nParam);
  MATRIX DeltaTheta(nParam, 1);
  MATRIX theta(nParam, 1);

  param = .2;
  param(1,1) = 1;

  param = .4;

  functionvalue = 900000;
  theta = param;

  //  CalculateNewMixLLNucl(indlist, RELAP, resp, theta, model);

  // Start value for parameters
  //  savethis = CalculateMixLLNucl(indlist, RELAP, resp, theta, model);

  // printf("Start likelihood: %f\n", savethis);

  for (iter = 1; iter<=maxIter; iter++)
  {
    if (iter==maxIter)
      code = 4;

    //    printf("Iteration : %d\n", iter);

    oldfunctionvalue = functionvalue;

    // Calculate the partial derivatives

    deriv = SpecialDeriv(theta, indlist, RELAP, resp, model);
    functionvalue = deriv[2](1,1);

    // Check
    if (fabs(oldfunctionvalue-functionvalue)<TOL)
    {
      code = 3;
      break;
    }

    // Modify diagonal
    for (i = 1; i<=nParam; i++)
      deriv[1](i,i) += 2.0/theta(i,1);
    DeltaTheta = Inverse(deriv[1])*deriv[0];

    downscale = 1;

    // Use 2 as 1 is the mean and it is unrestricted
    for (i = 2; i<=nParam; i++)
    {
      if (theta(i,1)<EDGE && DeltaTheta(i,1)>theta(i,1))
        DeltaTheta(i,1) = 0;

      if (DeltaTheta(i,1)>theta(i,1))
      {
        downscale = min(downscale, (1-EPSILON)*(fabs(theta(i,1)/DeltaTheta(i,1))));
      }  
    }

    DeltaTheta = DeltaTheta*downscale; 


    // Check stephalving
    step = 1;

    claus = CalculateMixLLNucl(indlist, RELAP, resp, theta-(DeltaTheta*step), model);

    //    printf("Suggested new theta:\n");
    //    (theta-(DeltaTheta*step)).Print();
    //    printf("Starting stephalving with above param resulting in %f\n", claus);
    
    while ( claus >= savethis)
    {
      step *= .5;

      printf("  Stephalving: %f %f\n", step, claus);

      if (step<1e-5)    
      {
        code = 2;
        break;
      }
      claus = CalculateMixLLNucl(indlist, RELAP, resp, theta-(DeltaTheta*step), model);
    }

    theta = theta - DeltaTheta*step;
    //    printf("-------------\nNew Theta\n-------------\n");
    //    (theta).Print();
    //    printf("-------------\n");

  
    if (code>0)
      break;

  }

  //  for (i = 1; i<=theta.Rows(); i++)
  //    {
  //      if (theta(i,1)<EDGE
  //    }

    claus = CalculateMixLLNucl(indlist, RELAP, resp, theta-(DeltaTheta*step), model);
    printf("Convergence code : %d   (Likelihood : %f)\n", code, claus);
    printf("-------------\n");
    theta.Print();
    printf("-------------\n");

    //    if (model==4)
    //      SpecialFunctionPiHat(indlist, traitnum);

    //    CalculateNewMixLLNucl(indlist, RELAP, resp, param, model);

  return claus;
}


double SpecialFunctionPiHat(individual *indlist, int traitnum, int markernum)
{
  int i, j, nPos, pos, nPed, bestpos;
  double maxdiff, difference, bestfull, bestadd, bestno;

  // N personer
  int nSize = 4;
  int noresdom = 1;

  int nVC = 5-noresdom;  // I, A, Pi, K2

  namelist *Pedigrees;
  char minibuf[50];
  FILE *F, *F2;

  MATRIX GenA(nSize,nSize), GenD(nSize,nSize);
  MATRIX *IBD, *K2;
  MATRIX *y, *x, *VC, start(nVC,1), beta(1,1);

  Pedigrees = MakePedigreeList(indlist);
  nPed = listlen(Pedigrees);
  freelist(Pedigrees);

  IBD   = new MATRIX[nPed](nSize,nSize);
  K2    = new MATRIX[nPed](nSize,nSize);
  y     = new MATRIX[nPed](nSize,1);
  x     = new MATRIX[nPed](nSize,1);

  VC    = new MATRIX[nPed*nVC](nSize,nSize);  // I, A, D, Pi, K2

  int pedsize[nPed];

  individual *ind;

  int rpos, rid1, rid2, rFID, oldrFID;
  double rp0, rp1, rp2;
  double full, nodom, noeff;

  ind = SelectPedigree("F1", indlist);
  //  GenA = MMakeKinshipMatrix(ind)*2;            // Res. Additive
  //  GenA.Resize(nSize, nSize);
  GenA = 0.5;
  GenA = GenA + MakeIdentMatrix(nSize)*0.5;
  //  GenD = GenA;      // Res. Dominance
  FreeIndividualList(ind);

  for (i=0; i<nPed ; i++)
  {
    pedsize[i] = nSize;
    x[i] = 1;
    VC[i]   = MakeIdentMatrix(nSize);
    VC[i + nPed] = GenA;
    VC[i + 2*nPed] = GenD;
  }

  for (ind = indlist, i=0, j=1; ind; ind = ind->next)
  {
    if (!founder(ind)) {
      y[i](j,1) = trait(ind,traitnum);
      j++;
      if (j>y[i].Rows()) 
	{
	  i++;
	  j=1;
	}
    }
  }

  bestpos = 0;
  maxdiff = -1000000.0;

  F2 = fopen ("pihat.out","a+");

  if (F2 == 0)
  {
    printf("Could not find file micropos.dat\n");
    exit(1);
  }

  #ifdef multip
  nPos = 34;
  {
  #else
  for (nPos = 6; nPos <= 60; nPos = nPos + 2)
  {
  
  #endif

    pos = (int) nPos;
    pos = (int) (markernum-1)*10;

    sprintf(minibuf, "awk '$1 == %d { print } ' micro.dat > micropos.dat\n", pos);
    system(minibuf);

 
    for (i=0; i<nPed ; i++)
    {
      IBD[i] = 0;
      K2[i] = 0;
      for (j=1; j<=nSize; j++)
      {
	IBD[i](j,j) = 1;
	K2[i](j,j) = 1;
      }
    }
    // Read the file
    F = fopen ("micropos.dat","r");

    if (F == 0)
    {
      printf("Could not find file micropos.dat\n");
      exit(1);
    }

    oldrFID = -1;

    while (fscanf(F, "%d %d %d %d %lf %lf %lf\n", &rpos, &rFID, &rid1, &rid2, &rp0, &rp1, &rp2 ) != EOF)
    {      
      if (rpos == pos)
      {
	if (oldrFID != rFID)
	  i = 0;

	oldrFID = rFID;

	i++;

	#ifdef sibs
	rid1 -=2;
	rid2 -=2;
	#endif

	// Reads the relationships and reduces precision to save a little time
	IBD[rFID-1](rid1, rid2) = 0.5*max(rp1,0) + max(0, rp2);
	IBD[rFID-1](rid2, rid1) = 0.5*max(rp1,0) + max(0, rp2);

	K2[rFID-1](rid1, rid2) = max(rp2,0);
	K2[rFID-1](rid2, rid1) = max(rp2,0);
      } 
    }
    fclose(F);

    // Fixes the rest.

    //    printf("------------\nPosition: %d\n", nPos);
    printf("------------\nMarker: %d\n", markernum);

    for (i=0; i<nPed ; i++)
    {
      if (noresdom) {
        VC[i+2*nPed] = IBD[i];
        VC[i+3*nPed] = K2[i];
      }
      else {
        VC[i+3*nPed] = IBD[i];
        VC[i+4*nPed] = K2[i];
      }
    }

    //  VC[3*nPed].Print();

    start = .4;
    full = MaximizeMixedModel(nPed, pedsize, y, x, nVC, VC, start, beta, METHOD_ML, 1);

    start.Resize(nVC-1,1);
    start = .4;

    nodom = MaximizeMixedModel(nPed, pedsize, y, x, nVC-1, VC, start, beta, METHOD_ML, 1);

    start.Resize(nVC-2,1);
    start = .4;

    noeff = MaximizeMixedModel(nPed, pedsize, y, x, nVC-2, VC, start, beta, METHOD_ML, 1);


    difference = full - noeff;

    printf("%d %f %f %f  : %f\n", nPos, full, nodom, noeff, difference);

    if (difference>maxdiff)
    {
      printf("------------------\n");
      bestpos  = pos;
      bestfull = full;
      bestadd  = nodom;
      bestno   = noeff;

      maxdiff = difference;
      
    }
    start.Resize(nVC,1);
  }

  fprintf (F2, "%d %f %f -1 %f %d \n", bestpos, bestfull, bestadd, bestno, bestpos - 33);
  fclose(F2);

  delete[] IBD;  
  delete[] K2;  
  delete[] VC;
  delete[] x;
  delete[] y;

  return 0;
}



// ----------------------------------------------------------
//
//  Functions for sibs
//
// ----------------------------------------------------------

double SpecialFunctionSibsPiHat(individual *indlist, int traitnum, int model)
{
  int i, j, nPed;

  // N personer
  int nSize = 2;

  namelist *Pedigrees;

  MATRIX GenA(nSize,nSize), GenI(nSize,nSize);
  MATRIX *RELAP, *IBD, *K2;
  MATRIX *y, *x, *VC, start(4,1), beta(1,1);

  Pedigrees = MakePedigreeList(indlist);
  nPed = listlen(Pedigrees);
  freelist(Pedigrees);

  IBD   = new MATRIX[nPed](nSize,nSize);
  K2    = new MATRIX[nPed](nSize,nSize);
  y     = new MATRIX[nPed](nSize,1);
  x     = new MATRIX[nPed](nSize,1);

  VC    = new MATRIX[nPed*4](nSize,nSize);  // I, A, D, Pi, K2

  int pedsize[nPed];

  individual *ind;

  double full;

  //  ind = SelectPedigree("F1", indlist);
  //  GenA = MMakeKinshipMatrix(ind)*2;            // Res. Additive
  //  GenD = MMakeDelta7Matrix(ind, GenA*.5);      // Res. Dominance
  //  FreeIndividualList(ind);

  GenI = MakeIdentMatrix(2);
  //  ind = SelectPedigree("F1", indlist);
  GenA = GenI;            // Res. Additive
  GenA(1,2) = .5;
  GenA(2,1) = .5;


  for (i=0; i<nPed ; i++)
  {
    pedsize[i] = nSize;
    x[i] = 1;
    VC[i]   = MakeIdentMatrix(nSize);
    VC[i + nPed] = GenA;

    // Reads the relationships and reduces precision to save a little time
    IBD[i+2*nPed](1, 2) = RELAP[i](1,3)+ 0.5*RELAP[i](1,2);
    IBD[i+2*nPed](2, 1) = RELAP[i](1,3)+ 0.5*RELAP[i](1,2);

    K2[i+3*nPed](1, 2) = RELAP[i](1,3);
    K2[i+3*nPed](2, 1) = RELAP[i](1,3);

  }

  for (ind = indlist, i=0, j=1; ind; ind = ind->next)
  {
    if (!founder(ind)) {
    y[i](j,1) = trait(ind,traitnum);
    j++;
    if (j>y[i].Rows()) 
    {
      i++;
      j=1;
    }
    }
  }

   //  VC[3*nPed].Print();

    switch (model) {
    case 1:  start = .4;
             full = MaximizeMixedModel(nPed, pedsize, y, x, 4, VC, start, beta, METHOD_ML, 1);
	     break;
    case 2:  start.Resize(3,1);
             start = .4;
	     full = MaximizeMixedModel(nPed, pedsize, y, x, 3, VC, start, beta, METHOD_ML, 1);
	     break;
    case 4:  start.Resize(2,1);
             start = .4;
	     full = MaximizeMixedModel(nPed, pedsize, y, x, 2, VC, start, beta, METHOD_ML, 1);
	     break;
    }

  delete[] IBD;  
  delete[] K2;  
  delete[] VC;
  delete[] x;
  delete[] y;

  return full;

}

double SpecialFunctionSibsMix(individual *indlist, int traitnum, int position, MATRIX RELAP[])
{
  //  int nVC = 5;             // I, A, D, g, k
  int i, j, iter;
  int model = position;
  // Loop variables

  int nodom = 1;
  int nParam = 6-nodom;  // Mean, I, A, D; g, k
  int maxIter = 40;
  int code = 0;
  double downscale, step=1;
  double savethis;
  double oldfunctionvalue, functionvalue;
  double claus;

  int nPedSize = 2;
  int nPedigrees = listlen(indlist)/(nPedSize+2);

  individual *ind;

  MATRIX *resp, *deriv;
  MATRIX param(nParam, 1);

  resp = new MATRIX[nPedigrees](nPedSize,1);

  for (ind = indlist, i=0, j=1; ind; ind = ind->next)
  {
    if (!founder(ind))
      {
    resp[i](j,1) = trait(ind,traitnum);
    j++;
    if (j>resp[i].Rows()) 
    {
      i++;
      j=1;
    }
      }
  }

  switch (model) {
  case 2:  // No dominance effect
    nParam -= 1;
    param.Resize(nParam,1);   
    break;
  case 4: // No locus effect
    nParam -= 2;
    param.Resize(nParam,1);   
    break;
  }

  MATRIX newh(nParam, 1);
  MATRIX dParam(nParam, 1),  SecOrd(nParam, nParam);
  MATRIX DeltaTheta(nParam, 1);
  MATRIX theta(nParam, 1);

  param = .2;
  param(1,1) = 1;

  param = .4;

  functionvalue = 900000;
  theta = param;

  CalculateMixLLSibs(indlist, RELAP, resp, theta, model);

  // Start value for parameters
  savethis = CalculateMixLLSibs(indlist, RELAP, resp, theta, model);

  // printf("Start likelihood: %f\n", savethis);

  for (iter = 1; iter<=maxIter; iter++)
  {
    if (iter==maxIter)
      code = 4;

    //    printf("Iteration : %d\n", iter);

    oldfunctionvalue = functionvalue;

    // Calculate the partial derivatives

    deriv = SpecialDerivSibs(theta, indlist, RELAP, resp, model);
    functionvalue = deriv[2](1,1);

    // Check
    if (fabs(oldfunctionvalue-functionvalue)<TOL)
    {
      code = 3;
      break;
    }

    // Modify diagonal
    for (i = 1; i<=nParam; i++)
      deriv[1](i,i) += 2.0/theta(i,1);
    DeltaTheta = Inverse(deriv[1])*deriv[0];

    downscale = 1;

    // Use 2 as 1 is the mean and it is unrestricted
    for (i = 2; i<=nParam; i++)
    {
      if (theta(i,1)<EDGE && DeltaTheta(i,1)>theta(i,1))
        DeltaTheta(i,1) = 0;

      if (DeltaTheta(i,1)>theta(i,1))
      {
        downscale = min(downscale, (1-EPSILON)*(fabs(theta(i,1)/DeltaTheta(i,1))));
      }  
    }

    DeltaTheta = DeltaTheta*downscale; 


    // Check stephalving
    step = 1;

    claus = CalculateMixLLSibs(indlist, RELAP, resp, theta-(DeltaTheta*step), model);

    //    printf("Suggested new theta:\n");
    //    (theta-(DeltaTheta*step)).Print();
    //    printf("Starting stephalving with above param resulting in %f\n", claus);
    
    while ( claus >= savethis)
    {
      step *= .5;

      printf("  Stephalving: %f %f\n", step, claus);

      if (step<1e-5)    
      {
        code = 2;
        break;
      }
      claus = CalculateMixLLSibs(indlist, RELAP, resp, theta-(DeltaTheta*step), model);
    }

    theta = theta - DeltaTheta*step;
    //    printf("-------------\nNew Theta\n-------------\n");
    //    (theta).Print();
    //    printf("-------------\n");

  
    if (code>0)
      break;

  }

  //  for (i = 1; i<=theta.Rows(); i++)
  //    {
  //      if (theta(i,1)<EDGE
  //    }

    claus = CalculateMixLLSibs(indlist, RELAP, resp, theta-(DeltaTheta*step), model);
    printf("Convergence code : %d   (Likelihood : %f)\n", code, claus);
    printf("-------------\n");
    theta.Print();
    printf("-------------\n");

    //    if (model==4)
    //      SpecialFunctionPiHat(indlist, traitnum);

    CalculateMixLLSibs(indlist, RELAP, resp, param, model);

  return claus;
}

MATRIX *SpecialDerivSibs(MATRIX &param, individual *indlist, MATRIX RELAP[], MATRIX resp[], int model)
{
  int i, j, nParam;
  double h, savethis;

  //  printf("Parameters passed to Deriv function:\n");
  //  param.Print();
  //  printf("--- End of parameters\n");

  nParam = param.Rows();

  MATRIX newh(nParam, 1), dParam (nParam, 1), dSecOrd(nParam, nParam), xplus(nParam,1), yplus(nParam,1);
  MATRIX *res;

  res = new MATRIX[3](nParam, nParam);

  res[0].Resize(nParam,1);
  res[2].Resize(1,1);

  h = 0.0001;

  savethis = CalculateMixLLSibs(indlist, RELAP, resp, param, model);
  //  printf("Savethis : %f\n", savethis);

  for (i = 0; i < nParam; i++)
  {
    newh = 0;
    newh(i+1,1) = h;
    dParam(i+1,1) = CalculateMixLLSibs(indlist, RELAP, resp, param+newh, model) - savethis;

    //    printf("  Nr %d : %f\n", i+1, dParam(i+1,1)+savethis);

    for (j = 0; j<= i; j++)
    {
      if (i==j)
      {
	dSecOrd(i+1,j+1) =  (CalculateMixLLSibs(indlist, RELAP, resp, param+(newh*2), model) - 2*dParam(i+1,1) - savethis)/(h*h);
      }
      else
      {
	yplus = 0;
	yplus(j+1,1) = h;
	
	dSecOrd(i+1,j+1) =  (CalculateMixLLSibs(indlist, RELAP, resp, param+newh+yplus, model) - 
			     dParam(i+1,1) - dParam(j+1,1) - savethis)/(h*h);
	dSecOrd(j+1,i+1) =  dSecOrd(i+1,j+1);
      }
    } 

  }

  dParam *= 1.0/h;

  res[0] = dParam;
  res[1] = dSecOrd;
  res[2](1,1) = savethis;

  return res;
}

double CalculateMixLLSibs(individual *indlist, MATRIX RELAP[], MATRIX resp[], MATRIX param, int model)
{
  int nFamilyMembers = 2;  
  int i18, i;
  int nPedigrees;
  // Loop variables
  //  int i, i29, i39, i56, i57, i58, i59, i67, i68, i69, i78, i79;

  double res, dLogDet, mean;
  double pprod, famres, meanvalue;

  MATRIX OrigOmega(nFamilyMembers,nFamilyMembers), Omega(nFamilyMembers,nFamilyMembers);
  MATRIX InvOmega(nFamilyMembers,nFamilyMembers);
  MATRIX Delta7(nFamilyMembers,nFamilyMembers), GenA(nFamilyMembers,nFamilyMembers);
  MATRIX GenD(nFamilyMembers,nFamilyMembers), GenI(nFamilyMembers,nFamilyMembers);
  MATRIX GenG(nFamilyMembers,nFamilyMembers);
  MATRIX y(nFamilyMembers,1), yT(1,nFamilyMembers), miniparam(3,1);
  namelist *Pedigrees;

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  GenI = MakeIdentMatrix(nFamilyMembers);
  //  ind = SelectPedigree("F1", indlist);
  GenA = GenI;            // Res. Additive
  GenA(1,2) = .5;
  GenA(2,1) = .5;
  //  GenD = MMakeDelta7Matrix(ind, GenA*.5);      // Res. Dominance
  //  FreeIndividualList(ind);

  mean = param(1,1);
  miniparam = 0;

  // Mean, I, A, pi, K2

  switch (model) {
    case 1: miniparam(2,1) = 0.5*param(4,1);
      miniparam(3,1) = param(4,1) + param(5,1); break;
      case 2: miniparam(2,1) = 0.5*param(4,1);
	miniparam(3,1) = param(4,1); break;
      case 4: miniparam = 0; break;
      }

  OrigOmega = GenI*(param(2,1) + miniparam(3,1)) + GenA*param(3,1); // I

  res = 0;
  meanvalue = 0;

  for (i18=0; i18<nPedigrees; i18++)
  {
    y = resp[i18]-mean;
    yT = Transpose(y);

    // Adds I and 2 in the diagonal

    res -= 0.5*log(2*PI)*nFamilyMembers;

    famres = 0;

  if (model < 4)
  {

  for (i = 0; i<=2; i++)
  {
    if (RELAP[i18](1, i+1)==0) 
      continue;

    pprod = 0;

    pprod += log(RELAP[i18]( 1,i+1));

    Omega = OrigOmega;

    Omega(1,2) += miniparam(i+1,1) ;
    Omega(2,1) = Omega(1,2);

    InvOmega = Inverse(Omega, &dLogDet);
    famres += exp(pprod - 0.5*dLogDet -0.5*(yT*InvOmega*y)(1,1));

    //    famres += exp(pprod);
    //    printf("Famres (%d): %15f    %d %d %d %d %d %d\n", i, famres*100000, i34, i35, i36, i45, i46, i56);

  }
  }   
  else  // Model == 4
    {
    Omega = OrigOmega;
    //  Omega.Print();
    InvOmega = Inverse(Omega, &dLogDet);
    famres += exp(- 0.5*dLogDet -0.5*(yT*InvOmega*y)(1,1));
    }
    
    res += log(famres);
  }

  return (res);

}



//----------------------------------------------------------
//
//  For 4 sibs, no parents
//
// ---------------------------------------------------------

int Code(long value) {
  int res = 0;
  int i = 0;
 
  while (value>0) {
    res += (value % 10) * (int) (exp(i*log(3))+.1);
    value /= 10;    
    i++;
  }
  return res+1;
}


long InvCode(int value) {
  int res = 0;
  int i = 0;
  
  value -=1;
 
  while (value>0) {
    res += (value % 3) * (long) (exp(i*log(10))+.01);
    value /= 3;    
    i++;
  }
  return res;

}


double CalculateFullSibsMixLL(individual *indlist, int model, int markernum, MATRIX LIKE[], MATRIX resp[], MATRIX param)
{
  long possible[8];
  int nFamilyMembers = 4, nFam;  
  int i, fam, combination, comb2;
  double like = 0, dLogDet, famres, mean;

  MATRIX OrigOmega(nFamilyMembers,nFamilyMembers), Omega(nFamilyMembers,nFamilyMembers);
  MATRIX InvOmega(nFamilyMembers,nFamilyMembers), GenA(nFamilyMembers,nFamilyMembers), GenI(nFamilyMembers,nFamilyMembers);
  MATRIX miniparam(3,1), y(4,1), yT(1,4);

  nFam = listlen(indlist)/6;

  possible[0] = 111111;
  possible[1] = 110100;
  possible[2] = 101010;
  possible[3] =  11001;
  possible[4] =    111;
  possible[5] = 100001;
  possible[6] =  10010;
  possible[7] =   1100;

  mean = param(1,1);
  miniparam = 0;
  switch (model) {
  case 1: miniparam(2,1) = 0.5*param(4,1);
          miniparam(3,1) = param(4,1) + param(5,1); break;
  case 2: miniparam(2,1) = 0.5*param(4,1);
          miniparam(3,1) = param(4,1); break;
  case 4: miniparam = 0; break;
  }

  like = 0;
  GenI = MakeIdentMatrix(4);
  GenA = 0.5;
  GenA = GenA + GenI*0.5;

  OrigOmega = GenI*(param(2,1) + miniparam(3,1)) + GenA*param(3,1); // I

  //  OrigOmega.Print();

  for (fam = 0; fam<nFam; fam++) {
    // For en familie
    y = resp[fam] - mean;
    yT = Transpose(y);

    like -= 0.5*log(2*PI)*nFamilyMembers;

    famres = 0;
    if (model<4) {
      for (i= 1; i<=LIKE[fam].Rows(); i++) {
	//	printf("   CHECK %d ", i);
	// Combination has positive likelihood
	if (LIKE[fam](i,1)>0) {

	  combination = InvCode(i);

	  //	  printf("Fam %d  (%d) --- %6d   %f    ooo  %f\n", fam, i, combination, famres, LIKE[fam](i,1));

	  comb2 = combination;
	  Omega = OrigOmega;

	  Omega(3,4) += miniparam((combination % 10)+1,1) ;
	  Omega(4,3) = Omega(3,4);
	  combination /= 10;

	  Omega(2,4) += miniparam((combination % 10)+1,1) ;
	  Omega(4,2) = Omega(2,4);
	  combination /= 10;

	  Omega(2,3) += miniparam((combination % 10)+1,1) ;
	  Omega(3,2) = Omega(2,3);
	  combination /= 10;

	  Omega(1,4) += miniparam((combination % 10)+1,1) ;
	  Omega(4,1) = Omega(1,4);
	  combination /= 10;

	  Omega(1,3) += miniparam((combination % 10)+1,1) ;
	  Omega(3,1) = Omega(1,3);
	  combination /= 10;

	  Omega(1,2) += miniparam((combination % 10) +1, 1) ;
	  Omega(2,1) = Omega(1,2);

	  //	  Omega.Print();

	  InvOmega = Inverse(Omega, &dLogDet);

	  famres += exp(log(LIKE[fam](i,1)) - 0.5*dLogDet -0.5*(yT*InvOmega*y)(1,1));
	  //	  famres += exp(log(LIKE[fam](Code(comb2),1)));
	  //  	  printf ("  **  %d - %f **  [[%f]]\n", Code(comb2), LIKE[fam](Code(comb2),1),famres);

	}     
      }
    }
    else  // Model == 4
      {
	Omega = OrigOmega;
	//  Omega.Print();
	InvOmega = Inverse(Omega, &dLogDet);
	famres += exp(- 0.5*dLogDet -0.5*(yT*InvOmega*y)(1,1));
      }
    like += log(famres);
  }
  return (-like);
}

MATRIX *SpecialDerivFullSibs(individual *indlist, int model, int markernum, MATRIX LIKE[], MATRIX RESP[], MATRIX &param)
{
  int i, j, nParam;
  double h, savethis;

  //  printf("Parameters passed to Deriv function:\n");
  //  param.Print();
  //  printf("--- End of parameters\n");

  nParam = param.Rows();

  MATRIX newh(nParam, 1), dParam (nParam, 1), dSecOrd(nParam, nParam), xplus(nParam,1), yplus(nParam,1);
  MATRIX *res;

  res = new MATRIX[3](nParam, nParam);

  res[0].Resize(nParam,1);
  res[2].Resize(1,1);

  h = 0.0001;

  savethis = CalculateFullSibsMixLL(indlist, model, markernum, LIKE, RESP, param);
  //  printf("Savethis : %f\n", savethis);

  for (i = 0; i < nParam; i++)
  {
    newh = 0;
    newh(i+1,1) = h;
    dParam(i+1,1) = CalculateFullSibsMixLL(indlist, model, markernum, LIKE, RESP, param+newh) - savethis;


    //    printf("  Nr %d : %f\n", i+1, dParam(i+1,1)+savethis);

    for (j = 0; j<= i; j++)
    {
      if (i==j)
      {
	dSecOrd(i+1,j+1) =  (CalculateFullSibsMixLL(indlist, model, markernum, LIKE, RESP, param+(newh*2)) - 2*dParam(i+1,1) - savethis)/(h*h);
      }
      else
      {
	yplus = 0;
	yplus(j+1,1) = h;
	
	dSecOrd(i+1,j+1) =  (CalculateFullSibsMixLL(indlist, model, markernum, LIKE, RESP, param+newh+yplus) - dParam(i+1,1) - dParam(j+1,1) - savethis)/(h*h);
	dSecOrd(j+1,i+1) =  dSecOrd(i+1,j+1);
      }
    } 

  }

  dParam *= 1.0/h;

  res[0] = dParam;
  res[1] = dSecOrd;
  res[2](1,1) = savethis;

  return res;
}

void PrintData(individual *indlist) {
  individual *ind;
  markerlist *marker;

  for (ind = indlist; ind; ind = ind->next) {
    printf("Individual: %s\n", ind->id);
    for (marker = ind->marker; marker; marker = marker->next) {
      printf("%d  %d\n", marker->allele1, marker->allele2);
    }   
  }  
}


double SpecialFunctionMixSibs(individual *indlist, int traitnum, int markernum, int model)
{
  long possible[8];
  int nFam;  
  int i, j, k, m, famno;
  long combination, combinations, comb;
  double like = 0, totallike;
  char minibuf[50];
  int curr[7], poss[7];
  double step=1, functionvalue, claus, downscale, savethis, oldfunctionvalue;
  individual *ind, *ind2, *gelist;
  markerlist *fmkr, *mmkr, *mkr, *mkr2;
  
  MATRIX *LIKE;
  MATRIX *RESP;

  if (listlen(indlist) == 0) {
    printf("No people in dataset\n");
    exit(1);
  }

  CopyAlleleFreq(markernum);

  nFam = listlen(indlist)/6;

  possible[0] = 111111;
  possible[1] = 110100;
  possible[2] = 101010;
  possible[3] =  11001;
  possible[4] =    111;
  possible[5] = 100001;
  possible[6] =  10010;
  possible[7] =   1100;

  LIKE = new MATRIX[nFam](730,1);
  RESP = new MATRIX[nFam](4,1);

  for (ind = indlist, i=0, j=1; ind; ind = ind->next) {
    if (!founder(ind)) {
      RESP[i](j,1) = trait(ind,traitnum);
      j++;
      if (j>RESP[i].Rows())    {
	i++;
	j=1;
      }
    }
  }


  for (famno = 0; famno<nFam; famno++) {

    // Start by doing genotype elimination
    sprintf(minibuf, "F%d", famno+1);

    gelist = GenotypeElimination(minibuf, markernum);

    //    printf("%s\n", minibuf);

    LIKE[famno] = 0;
    // Figure out the possible number of combinations
    combinations = 1;

    for (ind = gelist, j=1; ind; ind = ind->next, j++) {
      poss[j] = listlen(ind->marker);
      combinations *= poss[j];
      curr[j] = 1;
    }

    totallike = 0;

    for (comb = 1; comb<=combinations; comb++) {
      //      printf ("Combination %d (%ld)\n", comb, combinations);

      fmkr = markernumber(gelist, curr[1]);
      mmkr = markernumber(gelist->next, curr[2]);

      // Calculate the likelihood of this combination
      like = allelefreq[fmkr->allele1]*allelefreq[fmkr->allele2]*allelefreq[mmkr->allele1]*allelefreq[mmkr->allele2];

; // Multiply the allele frequencies

      if (heterozygous(gelist, curr[1])) {
	like *= exp(4*log(.5));
      }
      if (heterozygous(gelist->next, curr[2])) {
	like *= exp(4*log(.5));
      }

      totallike += like;

      //      printf("Likelihood : %f  %d \n", totallike, 2*homozygous(gelist, curr[1]) + homozygous(gelist->next, curr[2]));

      // Figure out the combination
      combination = 0;
      //      printf("%d  : %d %d   %d %d --- \n", famno, fmkr->allele1, fmkr->allele2, mmkr->allele1,mmkr->allele2 );
      switch((2*homozygous(gelist, curr[1]) + homozygous(gelist->next, curr[2]))) {
      case 0: // Both heterozygous
	      m = 5;  // Da vi har 6 relationer
	      for(ind = gelist->next->next, j = 3; ind->next; ind = ind->next, j++) {
		for(ind2 = ind->next, k = j+1; ind2; ind2 = ind2->next, k++) {
		  mkr = markernumber(ind, curr[j]);
		  mkr2= markernumber(ind2, curr[k]);

		  if (mkr->allele1 == mkr2->allele1) {
		    combination +=(int) (exp(m*log(10))+.1);
		  }
		  if (mkr->allele2 == mkr2->allele2) {
		    combination +=(int) (exp(m*log(10))+.1);
		  }
		  m--;
		}
	      }
	      //	      printf("IBD combo1: %6d (%d) for family %d\n",combination, Code(combination), famno);
	      LIKE[famno](Code(combination),1) += like;
	      
	      break;
      case 1: // Mother homozygous
	      m = 5;
	      for(ind = gelist->next->next, j = 3; ind->next; ind = ind->next, j++) {
		for(ind2 = ind->next, k = j+1; ind2; ind2 = ind2->next, k++) {
		  mkr = markernumber(ind, curr[j]);
		  mkr2= markernumber(ind2, curr[k]);

		  if (mkr->allele1 == mkr2->allele1) {
		    combination +=(int) (exp(m*log(10))+.1);
		  }
		  m--;
		}
	      }
	      for (m=0; m<8; m++) {
		LIKE[famno](Code(combination+possible[m]),1) += like / 8.0;
		//		printf("IBD combo2: %6d\n",combination+possible[m]);
	      }
	      break;
      case 2: // Father homozygous
	      m = 5;
	      for(ind = gelist->next->next, j = 3; ind->next; ind = ind->next, j++) {
		for(ind2 = ind->next, k = j+1; ind2; ind2 = ind2->next, k++) {
		  mkr = markernumber(ind, curr[j]);
		  mkr2= markernumber(ind2, curr[k]);

		  if (mkr->allele2 == mkr2->allele2) {
		    combination +=(int) (exp(m*log(10))+.1);
		  }
		  m--;
		}
	      }
	      for (m=0; m<8; m++) {
		//		printf("IBD combo3: %6d\n",combination+possible[m]);
		LIKE[famno](Code(combination+possible[m]),1) += like / 8.0;
	      }
	      break;
      case 3: // Both homozygous
	      for (m=0; m<8; m++) {
		for (j=0; j<8; j++) {
		  LIKE[famno](Code(possible[j]+possible[m]),1) += like / 64.0;
		}
	      }
	      break;
      }

      // Updates
      for (j=1; j<7; j++) {
	if (curr[j] < poss[j]) {
	  curr[j]++;
	  break;
	}
	else {
	  curr[j] = 1;
	}	  
      }
    }

    //    printf("End of combinations  : %f\n",totallike);

    LIKE[famno] = LIKE[famno] * (1.0/(totallike));

    //    PrintData(gelist);

    FreeIndividualList(gelist);
  }

  int iter, maxIter = 200, code=-1;
  int nParam = 5;  // mean, I,A,Q1,K2

  MATRIX param(nParam,1);

  switch (model) {
  case 2:  // No dominance effect
    nParam -= 1;
    param.Resize(nParam,1);   
    break;
  case 4: // No locus effect
    nParam -= 2;
    param.Resize(nParam,1);   
    break;
  }

  MATRIX newh(nParam, 1);
  MATRIX dParam(nParam, 1),  SecOrd(nParam, nParam);
  MATRIX DeltaTheta(nParam, 1);
  MATRIX theta(nParam, 1), *deriv;

  param = .4;
  
  param(1,1) = 0;    // Mean value
  param(2,1) = .5;   // Measurement error
  param(3,1) = .05;  // Residual additive
  param(param.Rows(),1) = .2; // Dominance effect
  //  param(4,1) = .3;   // Additive
  
  functionvalue = 900000;
  theta = param;

  // Start value for parameters
  savethis =   CalculateFullSibsMixLL(indlist, model, markernum, LIKE, RESP, theta);

  printf("Start likelihood: %f\n", savethis);

  for (iter = 1; iter<=maxIter; iter++)
  {
    if (iter==maxIter)
      code = 4;



    //    printf("Iteration : %d\n", iter);

    oldfunctionvalue = functionvalue;

    // Calculate the partial derivatives

    deriv = SpecialDerivFullSibs(indlist, model, markernum, LIKE, RESP, theta);

    functionvalue = deriv[2](1,1);

    // Check
    if (fabs(oldfunctionvalue-functionvalue)<TOL)
    {
      code = 3;
      break;
    }

    // Modify diagonal
    deriv[1](1,1) += 2.0/(theta(1,1) - NEGINF);
    for (i = 2; i<=nParam; i++)
      deriv[1](i,i) += 2.0/theta(i,1);

    DeltaTheta = Inverse(deriv[1])*deriv[0];

    delete[] deriv;

    downscale = 1;

    // Use modification for 1 is the mean and it is unrestricted
    for (i = 1; i<=nParam; i++)
    {
      if (i>1) {
	if (theta(i,1)<EDGE && DeltaTheta(i,1)>theta(i,1))
	  DeltaTheta(i,1) = 0;

	if (DeltaTheta(i,1)>theta(i,1))
	  {
	    downscale = min(downscale, (1-EPSILON)*(fabs(theta(i,1)/DeltaTheta(i,1))));
	  }  
      }
      else {  // Modifying the mean
	if (theta(i,1)<EDGE+NEGINF && DeltaTheta(i,1)>theta(i,1))
	  DeltaTheta(i,1) = 0;

	//	if (DeltaTheta(i,1)>theta(i,1))
	//	  {
	//	    downscale = min(downscale, (1-EPSILON)*(fabs(theta(i,1)/DeltaTheta(i,1))));
	//	  }  
      }
    }

    DeltaTheta = DeltaTheta*downscale; 

    // Check stephalving
    step = 1;

    claus = CalculateFullSibsMixLL(indlist, model, markernum, LIKE, RESP, theta-(DeltaTheta*step));
    //    printf("Suggested new theta:\n");
    //    (theta-(DeltaTheta*step)).Print();
    //    printf("Starting stephalving with above param resulting in %f\n", claus);
    
    while ( claus >= savethis)
    {
      step *= .5;

      printf("  Stephalving: %f %f\n", step, claus);

      if (step<1e-5)    
      {
        code = 2;
        break;
      }
    claus = CalculateFullSibsMixLL(indlist, model, markernum, LIKE, RESP, theta-(DeltaTheta*step));
    }

    theta = theta - DeltaTheta*step;
  
    if (code>0)
      break;

  }

  claus = CalculateFullSibsMixLL(indlist, model, markernum, LIKE, RESP, theta-(DeltaTheta*step));

  printf("Convergence code : %d   (Likelihood : -%f)\n", code, claus);
  printf("-------------\n");
  theta.Print();
  printf("-------------\n");

  delete[] LIKE;
  delete[] RESP;

  return claus;
}

void AnalyzeMatrix(individual *indlist, char *filename, int traitnum) {
  individual *ind, *ind2;
  namelist   *Pedigrees, *ped;
  int i, pedsize;
  double full, model;
 
  Pedigrees = MakePedigreeList(indlist);
 
  int nPers      = listlen(indlist);
  int nPed       = listlen(Pedigrees);
  int nSize[nPed];
 
  MATRIX *PiHat, *Phi2, *Delta7, *y, *VC, *x;
  MATRIX InMatrix(nPers,nPers);

  InMatrix.FileReadSymmetric(filename,nPers,nPers);

  PiHat = new MATRIX[nPed](1,1);
  Phi2  = new MATRIX[nPed](1,1);
  Delta7= new MATRIX[nPed](1,1);
  y     = new MATRIX[nPed](1,1);
  x     = new MATRIX[nPed](1,1);
  VC    = new MATRIX[nPed*3](1,1);

  for (ind = indlist, i=1; ind; ind=ind->next, i++) {
    ind->tmpint1 = i;
  }
  for (ped = Pedigrees, i=0; ped; ped = ped->next, i++) {
    pedsize = 0;
    for (ind = indlist; ind; ind=ind->next) {
      if (!strcmpl(ind->pedigree, ped->name)) {
        ind->tmpint2 = i;
        if (IsGenotyped(ind)) {
          pedsize++;
	  nSize[i] = pedsize;
          ind->tmpint3 = pedsize;
        }
      }
    }

    PiHat[i].Resize(pedsize, pedsize);
    Phi2[i].Resize(pedsize, pedsize);
    Delta7[i].Resize(pedsize, pedsize);
    y[i].Resize(pedsize,1);
    x[i].Resize(pedsize,1);
    VC[i].Resize(pedsize, pedsize);
    VC[i+nPed].Resize(pedsize, pedsize);
    VC[i+2*nPed].Resize(pedsize, pedsize);

    x[i] = 1;
    VC[i] = MakeIdentMatrix(pedsize);
  }

  for (ped = Pedigrees, i=0; ped; ped = ped->next, i++) {
    for (ind = indlist; ind; ind=ind->next) {
      if (!strcmpl(ind->pedigree, ped->name) && IsGenotyped(ind)) {
	y[i](ind->tmpint3,1) = trait(ind, traitnum);
	for (ind2 = ind; ind2; ind2=ind2->next) {
	  if (!strcmpl(ind->pedigree, ind2->pedigree) && IsGenotyped(ind2)) {
	    if (ind==ind2) {
	      PiHat[i](ind->tmpint3, ind2->tmpint3) = 1;
	      Phi2[i](ind->tmpint3, ind2->tmpint3) = 1;
	      Delta7[i](ind->tmpint3, ind2->tmpint3) = 1;
	    }
	    else {
	      PiHat[i](ind->tmpint3, ind2->tmpint3) = InMatrix(ind->tmpint1, ind2->tmpint1);
	      PiHat[i](ind2->tmpint3, ind->tmpint3) = InMatrix(ind->tmpint1, ind2->tmpint1);
	      Phi2[i](ind->tmpint3, ind2->tmpint3) = .5;
	      Phi2[i](ind2->tmpint3, ind->tmpint3) = .5;
	      Delta7[i](ind->tmpint3, ind2->tmpint3) = .25;
	      Delta7[i](ind2->tmpint3, ind->tmpint3) = .25;
	    }
	  }
	}
      }
    }
    VC[i+nPed] = Phi2[i];
    VC[i+2*nPed] = PiHat[i];
  }

  freelist(Pedigrees);

  MATRIX start(3,1), beta(1,1);

  beta = 1;
  start= 1;
  start(1,1) = 1;

  full = MaximizeMixedModel(nPed, nSize, y, x, 3, VC, start, beta, METHOD_ML, 1, 1);

  start.Resize(2,1);
  start = 1;
  start(1,1) = 1;

  model= MaximizeMixedModel(nPed, nSize, y, x, 2, VC, start, beta, METHOD_ML, 1, 1);


  delete[] y;
  delete[] x;
  delete[] VC;
  delete[] Phi2;
  delete[] Delta7;
  delete[] PiHat;
}


void AnalyzeAllMatrices(individual *indlist, char *filename, int traitnum) {
  individual *ind, *ind2;
  namelist   *Pedigrees, *ped;
  int i, pos, pedsize;
  double full, model;
  char buf[100];

  Pedigrees = MakePedigreeList(indlist);
 
  int nPers      = listlen(indlist);
  int nPed       = listlen(Pedigrees);
  int nSize[nPed];
 
  MATRIX *PiHat, *Phi2, *Delta7, *y, *VC, *x;
  MATRIX InMatrix(nPers,nPers);

  InMatrix.FileReadSymmetric(filename,nPers,nPers);

  PiHat = new MATRIX[nPed](1,1);
  Phi2  = new MATRIX[nPed](1,1);
  Delta7= new MATRIX[nPed](1,1);
  y     = new MATRIX[nPed](1,1);
  x     = new MATRIX[nPed](1,1);
  VC    = new MATRIX[nPed*3](1,1);

  for (ind = indlist, i=1; ind; ind=ind->next, i++) {
    ind->tmpint1 = i;
  }
  for (ped = Pedigrees, i=0; ped; ped = ped->next, i++) {
    pedsize = 0;
    for (ind = indlist; ind; ind=ind->next) {
      if (!strcmpl(ind->pedigree, ped->name)) {
        ind->tmpint2 = i;
        if (IsGenotyped(ind)) {
          pedsize++;
	  nSize[i] = pedsize;
          ind->tmpint3 = pedsize;
        }
      }
    }

    PiHat[i].Resize(pedsize, pedsize);
    Phi2[i].Resize(pedsize, pedsize);
    Delta7[i].Resize(pedsize, pedsize);
    y[i].Resize(pedsize,1);
    x[i].Resize(pedsize,1);
    VC[i].Resize(pedsize, pedsize);
    VC[i+nPed].Resize(pedsize, pedsize);
    VC[i+2*nPed].Resize(pedsize, pedsize);

    x[i] = 1;
    VC[i] = MakeIdentMatrix(pedsize);
  }

  for (ped = Pedigrees, i=0; ped; ped = ped->next, i++) {
    for (ind = indlist; ind; ind=ind->next) {
      if (!strcmpl(ind->pedigree, ped->name) && IsGenotyped(ind)) {
	y[i](ind->tmpint3,1) = trait(ind, traitnum);
	for (ind2 = ind; ind2; ind2=ind2->next) {
	  if (!strcmpl(ind->pedigree, ind2->pedigree) && IsGenotyped(ind2)) {
	    if (ind==ind2) {
	      PiHat[i](ind->tmpint3, ind2->tmpint3) = 1;
	      Phi2[i](ind->tmpint3, ind2->tmpint3) = 1;
	      Delta7[i](ind->tmpint3, ind2->tmpint3) = 1;
	    }
	    else {
	      PiHat[i](ind->tmpint3, ind2->tmpint3) = InMatrix(ind->tmpint1, ind2->tmpint1);
	      PiHat[i](ind2->tmpint3, ind->tmpint3) = InMatrix(ind->tmpint1, ind2->tmpint1);
	      Phi2[i](ind->tmpint3, ind2->tmpint3) = .5;
	      Phi2[i](ind2->tmpint3, ind->tmpint3) = .5;
	      Delta7[i](ind->tmpint3, ind2->tmpint3) = .25;
	      Delta7[i](ind2->tmpint3, ind->tmpint3) = .25;
	    }
	  }
	}
      }
    }
    VC[i+nPed] = Phi2[i];
    VC[i+2*nPed] = PiHat[i];
  }

  MATRIX start(3,1), beta(1,1);

  beta = 1;

  // Start by the null hypothesis
  start.Resize(2,1);
  start = .5;
  start(1,1) = 1;

  printf("ANALYZING THE NULL:\n");
  model= MaximizeMixedModel(nPed, nSize, y, x, 2, VC, start, beta, METHOD_ML, 1, 1);

  start.Resize(3,1);

  for (pos = 0; pos < (int) (80*ChromosomeLength()); pos++) {
    

    printf("Analyzing position %d\n", pos);
    // Reads the IBD-file
    sprintf(buf, "mibd-%d.pp",pos);
    InMatrix.FileReadSymmetric(buf,nPers,nPers);

    for (ped = Pedigrees, i=0; ped; ped = ped->next, i++) {
      for (ind = indlist; ind; ind=ind->next) {
        if (!strcmpl(ind->pedigree, ped->name) && IsGenotyped(ind)) {
	  for (ind2 = ind; ind2; ind2=ind2->next) {
	    if (!strcmpl(ind->pedigree, ind2->pedigree) && IsGenotyped(ind2)) {
  	      if (ind==ind2) {
	        PiHat[i](ind->tmpint3, ind2->tmpint3) = 1;
	      }
  	      else {
	        PiHat[i](ind->tmpint3, ind2->tmpint3) = InMatrix(ind->tmpint1, ind2->tmpint1);
	        PiHat[i](ind2->tmpint3, ind->tmpint3) = InMatrix(ind->tmpint1, ind2->tmpint1);
	      }
 	    }
	  }
        }
      }
      VC[i+2*nPed] = PiHat[i];
    }

    start= .5;
    start(1,1) = 1;

    full = MaximizeMixedModel(nPed, nSize, y, x, 3, VC, start, beta, METHOD_ML, 1, 1);
  }


  freelist(Pedigrees);


  delete[] y;
  delete[] x;
  delete[] VC;
  delete[] Phi2;
  delete[] Delta7;
  delete[] PiHat;
}

/*
 * Returns a vector of length (indlist), that holds a 1 for all the individuals
 * being complete cases (i.e. has a trait value for now)
 *
 */

MATRIX CompleteCases(individual *indlist, int traitnum)
{
  individual *ind;
  int i;

  int nPers = listlen(indlist);
  MATRIX res(nPers,1);

  res = 0;  
  for (ind = indlist, i=1; ind; ind = ind->next, i++) {
    if (!IsTraitMissing(ind, traitnum))
      res(i,1) = 1;
  }
  return res;
}


/******************************************

Rutine til at analysere heterogenitet

********************************************/


void AnalyzeMixtureModels(individual *indlist, char *filename, int traitnum) {
  individual *ind, *ind2, *ind3;
  namelist   *Pedigrees, *ped;
  int i, pedsize;
  double model;
  MATRIX *locusibd;
 
  Pedigrees = MakePedigreeList(indlist);
 
  int nPers      = listlen(indlist);
  int nPed       = listlen(Pedigrees);
  int nSize[nPed];
 
  MATRIX *PiHat, *Phi2, *Delta7, *y, *VC, *x;
  MATRIX Cases(nPers,1), Kinship, D7;

  locusibd = FileReadSymmetricKnownSize(filename,nPers/nPed);

  // Find complete cases
  Cases = CompleteCases(indlist, traitnum);

  PiHat = new MATRIX[nPed](1,1);
  Phi2  = new MATRIX[nPed](1,1);
  Delta7= new MATRIX[nPed](1,1);
  y     = new MATRIX[nPed](1,1);
  x     = new MATRIX[nPed](1,1);
  VC    = new MATRIX[nPed*5](1,1);

  // Sets local id for now for all possible pedigree members
  for (ind = indlist, i=1; ind; ind=ind->next, i++) {
    ind->localid = i; 
  }

  // Calculate the complete case pedigree sizes
  // Set tmpint1 to pedigree number
  // Set tmpint2 to internal number of complete cases
  for (ped = Pedigrees, i=0; ped; ped = ped->next, i++) {
    pedsize = 0;
    for (ind = indlist; ind; ind=ind->next) {
      if (!strcmpl(ind->pedigree, ped->name)) {
        ind->tmpint1 = i;
	if (Cases(ind->localid,1)) {
	  pedsize++;
	  nSize[i] = pedsize;
          ind->tmpint2 = pedsize;	  
	}
      }
    }

    PiHat[i].Resize(pedsize, pedsize);
    Phi2[i].Resize(pedsize, pedsize);
    Delta7[i].Resize(pedsize, pedsize);
    y[i].Resize(pedsize,1);
    x[i].Resize(pedsize,1);
    VC[i].Resize(pedsize, pedsize);
    VC[i+nPed].Resize(pedsize, pedsize);
    VC[i+2*nPed].Resize(pedsize, pedsize);
    VC[i+3*nPed].Resize(pedsize, pedsize);
    VC[i+4*nPed].Resize(pedsize, pedsize);

    x[i] = 1;
    VC[i] = MakeIdentMatrix(pedsize);
    VC[i+3*nPed] = MakeIdentMatrix(pedsize);
  }

  for (ped = Pedigrees, i=0; ped; ped = ped->next, i++) {

    ind3 = SelectPedigree(ped->name, indlist);
    Kinship.Resize(listlen(ind3), listlen(ind3));
    Kinship = MMakeKinshipMatrix(ind3);
    
    D7.Resize(listlen(ind3), listlen(ind3));
    D7 = MMakeDelta7Matrix(ind3, Kinship);


    FreeIndividualList(ind3);

    for (ind = indlist; ind; ind=ind->next) {
      // Check that correct pedigree
      if (!strcmpl(ind->pedigree, ped->name)) {
	y[i](ind->tmpint2,1) = trait(ind, traitnum);
	for (ind2 = ind; ind2; ind2=ind2->next) {
	  // Go through all other individuals
	  if (!strcmpl(ind->pedigree, ind2->pedigree)) {
	    PiHat[i](ind->tmpint2, ind2->tmpint2) = locusibd[i](ind->tmpint2, ind2->tmpint2);	    
	    Phi2[i](ind->tmpint2, ind2->tmpint2) = 2*Kinship(ind->tmpint2, ind2->tmpint2);
	    Delta7[i](ind->tmpint2, ind2->tmpint2) = D7(ind->tmpint2, ind2->tmpint2);

	    // Make symmetric
	    PiHat[i](ind2->tmpint2, ind->tmpint2)  = PiHat[i](ind->tmpint2, ind2->tmpint2);
	    Phi2[i](ind2->tmpint2, ind->tmpint2)   = Phi2[i](ind->tmpint2, ind2->tmpint2);  
	    Delta7[i](ind2->tmpint2, ind->tmpint2) = Delta7[i](ind->tmpint2, ind2->tmpint2);
	  }
	}
      }
    }
    VC[i+nPed] = Phi2[i];
    VC[i+4*nPed] = Phi2[i];
    VC[i+2*nPed] = PiHat[i];
  }

  delete[] locusibd;
  freelist(Pedigrees);


  MATRIX start(5,1), beta(2,1), include(5,1);

  beta = .1;
  beta(2,1) = 0;


  include = 0;
  include(4,1) = 1;   // I, non-carriers
  include(5,1) = 1;   // A, non-carriers

  printf("ANALYZING THE FULL MODEL:\n");

  // Start by the null hypothesis
  start = .1;
  start(1,1) = 1;     // I, carriers
  //  start(2,1) = max(RES2-RES3, 0.1);   // Additive
  start(2,1) = max(RES2-.5*RES3, 0.1);   // Additive
  start(3,1) = max(RES3,.1);   // QTL
  start(4,1) = 1;     // I, non-carriers
  start(5,1) = RES2;  // A, non-carriers
  // XXX Check dette .5 -> ,3 
  start(1,1) = 1;   // I
  start(2,1) = .1;  // Additive
  start(3,1) = 1;  // QTL
  start(4,1) = 1;   // 
  start(5,1) = .2;  // Additive


  // YYY
  model= MaximizeModifiedMixtureModel(nPed, nSize, y, x, 5, VC, start, beta, .25, include, METHOD_ML, 1, 1);

  delete[] y;
  delete[] x;
  delete[] VC;
  delete[] Phi2;
  delete[] Delta7;
  delete[] PiHat;
}




void AnalyzeMixtureModelsNULL(individual *indlist, char *filename, int traitnum) {
  individual *ind, *ind2, *ind3;
  namelist   *Pedigrees, *ped;
  int i, pedsize;
  double model;
 
  Pedigrees = MakePedigreeList(indlist);
 
  int nPers      = listlen(indlist);
  int nPed       = listlen(Pedigrees);
  int nSize[nPed];
 
  MATRIX *Phi2, *Delta7, *y, *VC, *x;
  MATRIX Cases(nPers,1), Kinship, D7;


  // Find complete cases
  Cases = CompleteCases(indlist, traitnum);

  Phi2  = new MATRIX[nPed](1,1);
  Delta7= new MATRIX[nPed](1,1);
  y     = new MATRIX[nPed](1,1);
  x     = new MATRIX[nPed](1,1);
  VC    = new MATRIX[nPed*4](1,1);

  // Sets local id for now for all possible pedigree members
  for (ind = indlist, i=1; ind; ind=ind->next, i++) {
    ind->localid = i; 
  }

  // Calculate the complete case pedigree sizes
  // Set tmpint1 to pedigree number
  // Set tmpint2 to internal number of complete cases
  for (ped = Pedigrees, i=0; ped; ped = ped->next, i++) {
    pedsize = 0;
    for (ind = indlist; ind; ind=ind->next) {
      if (!strcmpl(ind->pedigree, ped->name)) {
        ind->tmpint1 = i;
	if (Cases(ind->localid,1)) {
	  pedsize++;
	  nSize[i] = pedsize;
          ind->tmpint2 = pedsize;	  
	}
      }
    }

    Phi2[i].Resize(pedsize, pedsize);
    Delta7[i].Resize(pedsize, pedsize);
    y[i].Resize(pedsize,1);
    x[i].Resize(pedsize,1);
    VC[i].Resize(pedsize, pedsize);
    VC[i+nPed].Resize(pedsize, pedsize);
    VC[i+2*nPed].Resize(pedsize, pedsize);
    VC[i+3*nPed].Resize(pedsize, pedsize);

    x[i] = 1;
    VC[i] = MakeIdentMatrix(pedsize);
    VC[i+2*nPed] = MakeIdentMatrix(pedsize);
  }

  for (ped = Pedigrees, i=0; ped; ped = ped->next, i++) {

    ind3 = SelectPedigree(ped->name, indlist);
    Kinship.Resize(listlen(ind3), listlen(ind3));
    Kinship = MMakeKinshipMatrix(ind3);
    
    D7.Resize(listlen(ind3), listlen(ind3));
    D7 = MMakeDelta7Matrix(ind3, Kinship);


    FreeIndividualList(ind3);

    for (ind = indlist; ind; ind=ind->next) {
      // Check that correct pedigree
      if (!strcmpl(ind->pedigree, ped->name)) {
	y[i](ind->tmpint2,1) = trait(ind, traitnum);
	for (ind2 = ind; ind2; ind2=ind2->next) {
	  // Go through all other individuals
	  if (!strcmpl(ind->pedigree, ind2->pedigree)) {
	    Phi2[i](ind->tmpint2, ind2->tmpint2) = 2*Kinship(ind->tmpint2, ind2->tmpint2);
	    Delta7[i](ind->tmpint2, ind2->tmpint2) = D7(ind->tmpint2, ind2->tmpint2);

	    // Make symmetric
	    Phi2[i](ind2->tmpint2, ind->tmpint2)   = Phi2[i](ind->tmpint2, ind2->tmpint2);  
	    Delta7[i](ind2->tmpint2, ind->tmpint2) = Delta7[i](ind->tmpint2, ind2->tmpint2);
	  }
	}
      }
    }
    VC[i+nPed] = Phi2[i];
    VC[i+3*nPed] = Phi2[i];
  }

  freelist(Pedigrees);


  MATRIX start(4,1), beta(2,1), include(4,1);

  beta = .1;
  beta(2,1) = 0;


  include = 0;
  include(3,1) = 1;   // I, non-carriers
  include(4,1) = 1;   // A, non-carriers

  printf("ANALYZING THE NULL MODEL:\n");

  // Start by the null hypothesis
  start = .1;
  start(1,1) = 1;     // I, carriers
  //  start(2,1) = max(RES2-RES3, 0.1);   // Additive
  start(2,1) = max(RES2-.5*RES3, 0.1);   // Additive
  //  start(3,1) = max(RES3,.1);   // QTL
  start(3,1) = 1;     // I, non-carriers
  start(4,1) = RES2;  // A, non-carriers
  // XXX Check dette .5 -> ,3 
  start(1,1) = 1;   // I
  start(2,1) = .1;  // Additive
  start(3,1) = 1.2;   // 
  start(4,1) = .2;  // Additive


  // YYY
  model= MaximizeModifiedMixtureModelNULL(nPed, nSize, y, x, 4, VC, start, beta, .5, include, METHOD_ML, 1, 1);


  delete[] y;
  delete[] x;
  delete[] VC;
  delete[] Phi2;
  delete[] Delta7;
}






/*************


Returnerer en serie af matricer, der indeholde IBD-filerne


 ***********/

MATRIX* ReadSolarIBDFile(individual *indlist, char *fln, char *indexname) {
  MATRIX *res;
  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees, i, j, nFamSize, nCurrID, zipped;
  char name[NAMESIZE], buffer[256], oldfilename[256], filename[256];;
  float value;
  FILE *fptr;

  strcpy(filename, fln);
  strcpy(oldfilename, filename);

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  res = new MATRIX[nPedigrees];


  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    ind2 = SelectPedigree(pedi->name, indlist);
    nFamSize = listlen(ind2);
    res[i].Resize(nFamSize, nFamSize);

    for (ind = indlist, j=1; ind; ind = ind->next) {
      if (!strcmpl(ind->pedigree, pedi->name)) {
	ind->tmpint1 = i;  // pedigree
	ind->tmpint2 = j;  // localid
	j++;
      }
    }

    FreeIndividualList(ind2);
  }
  FreeNameList(Pedigrees);



  // Go though the index file and save the 
  // index numbers as localid 

  fptr = fopen(indexname,"r");
  if (fptr == NULL) {
    printf("WARNING: Could not open SOLAR indexfile: %s\n", indexname);
    return NULL;
  }

  // Run through the pedindex file and 
  // Assign the correct index numbers to the individuals
  while (fgets(buffer,256,fptr) != NULL) {
    sscanf(buffer, " %d %*d %*d %*d %*d %*d %*d %s", &nCurrID, &name);
    //    printf("%s   ---> %d\n", name, nCurrID);
    ind = FindIndListMember(indlist, name);
    if (ind==NULL) {
      printf("ERROR: Individual %s found in pedindex.out but not in data\n", name);
      return NULL;
    }
    ind->localid = nCurrID;
  }
  fclose(fptr);

  //
  // localid now holds the solar index
  //
  // Will now read the IBD file and make the pedigree matrices
  //
  zipped = 0;
  if (strstr(filename, ".gz")) {
    zipped = 1;
    sprintf(buffer, "gunzip %s", filename);
    system(buffer);
    filename[strlen(filename)-3] = '\0';
  }


  fptr = fopen(filename,"r");
  if (fptr == NULL) {
    printf("WARNING: Could not open SOLAR MIBD: %s\n", filename);
    return NULL;
  }

  while (fgets(buffer,256,fptr) != NULL) {
    sscanf(buffer, " %d %d %f", &i, &j, &value);
    if (i==j) {
      ind = FindIndListMemberLocalID(indlist, i);
      res[ind->tmpint1](ind->tmpint2, ind->tmpint2) = value;      
    }
    else {
      ind = FindIndListMemberLocalID(indlist, i);
      ind2 = FindIndListMemberLocalID(indlist, j);
      //      printf("Adding: %d %d   => %f\n", i, j, value);
      res[ind->tmpint1](ind->tmpint2, ind2->tmpint2) = value;      
      res[ind->tmpint1](ind2->tmpint2, ind->tmpint2) = value;      
      if (ind->tmpint1 != ind2->tmpint1) {
	printf("ERROR: Mess in pedigree for %s\n", ind->id);
      }
    }
  }
  fclose(fptr);

  //  res[1].Print();

  if (zipped) {
    sprintf(buffer, "gzip %s", filename);
    system(buffer);
    strcpy(filename, oldfilename);
  }

  return (res);

}


/*
 * Reads a pedipet IBD file and returns the result as a matrix
 *
 *
 *
 */

MATRIX* ReadPedipetIBDFile(individual *indlist, char *filename) {
  MATRIX *res;
  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees, i, j, nFamSize;
  char buffer[256];
  float value;
  FILE *fptr;

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  res = new MATRIX[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    ind2 = SelectPedigree(pedi->name, indlist);
    nFamSize = listlen(ind2);
    res[i].Resize(nFamSize, nFamSize);

    for (ind = indlist, j=1; ind; ind = ind->next) {
      if (!strcmpl(ind->pedigree, pedi->name)) {
	ind->tmpint1 = i;  // pedigree
	ind->tmpint2 = j;  // localid
	j++;
      }
    }

    FreeIndividualList(ind2);
  }
  FreeNameList(Pedigrees);

  for (ind = indlist, i=1; ind; ind = ind->next, i++) {
    ind->localid = i;
  }

  //
  // localid now holds the solar index
  //
  // Will now read the IBD file and make the pedigree matrices
  //
  fptr = fopen(filename,"r");
  if (fptr == NULL) {
    printf("WARNING: Could not open PEDIPET IBD file: %s\n", filename);
    return NULL;
  }

  while (fgets(buffer,256,fptr) != NULL) {
    sscanf(buffer, " %d %d %f", &i, &j, &value);
    if (i==j) {
      ind = FindIndListMemberLocalID(indlist, i);
      res[ind->tmpint1](ind->tmpint2, ind->tmpint2) = value;      
    }
    else {
      ind = FindIndListMemberLocalID(indlist, i);
      ind2 = FindIndListMemberLocalID(indlist, j);
      //      printf("Adding: %d %d   => %f\n", i, j, value);
      res[ind->tmpint1](ind->tmpint2, ind2->tmpint2) = value;      
      res[ind->tmpint1](ind2->tmpint2, ind->tmpint2) = value;      
      if (ind->tmpint1 != ind2->tmpint1) {
	printf("ERROR: Mess in pedigree for %s and %s  (%d, %d) when reading file %s\n", ind->id, ind2->id, ind->tmpint1, ind2->tmpint1, filename);
      }
    }
  }
  fclose(fptr);

  return (res);

}




/**********************************************

  Does multipoint analysis for X-linked loci
  Goes through the mibd files and 
			       
***********************************************/

double MultiXMarkerVC(individual *indlist, int traitnum,  MATRIX IBD[], MATRIX X[], int nCovar, int Dominance, int Household, int *covnum)
{
  int nFamSize;
  int i, j, k, pos, bestpos;
  int nVC, currentVC=0, chromlength; 
  double StartHyp, NullHyp, lldiff, bestlod, pval;
  char buf[20];
  FILE *F;

  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees;
  quanttrait *qttrait;
  MATRIX *y, *x, *VC, *keep;
  MATRIX Phi2, Delta7, HH, PhiMM, PhiFF, PhiMF;

  // Figure out the number of variance components
  nVC = 8;      // I 2*Phi + 3 times Phi_X + 3 times IBD_X
  if (Dominance)
    nVC++;      // I 2*Phi Delta7 IBD
  if (Household)
    nVC++;      // I 2*Phi Delta7 HH IBD      - If not Dom then just remove Delta from this comment to get the order

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  // Sets noprob of zero to missing
  for (ind=indlist; ind; ind=ind->next) {
    if (trait(ind,3) == 0) {
      qttrait = traitvalue(ind, 3);      
      qttrait->sysmiss = 1;      
    }
  }

  // Prepares matrices for input to maximize
  y  = new MATRIX[nPedigrees];
  x  = new MATRIX[nPedigrees];
  keep = new MATRIX[nPedigrees];
  VC = new MATRIX[nPedigrees*nVC];

  int nSize[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    nFamSize = 0;

    // Select the current pedigree
    ind2 = SelectPedigree(pedi->name, indlist);

    int nPedSize = listlen(ind2);
    keep[i].Resize(nPedSize,1);

    // Figure out which individuals to keep
    // Only checks for response and NOT covariates
    keep[i] = CompleteCases(ind2, traitnum);
    for (j = 2; j<nCovar; j++) {
      keep[i] = Hademard(keep[i], CompleteCases(ind2, covnum[j]));
    }

    // The actual number of complete individuals
    int nFamSize = (int) keep[i].Sum();

    Phi2.Resize(listlen(ind2), listlen(ind2));
    PhiMM.Resize(listlen(ind2), listlen(ind2));
    PhiFF.Resize(listlen(ind2), listlen(ind2));
    PhiMF.Resize(listlen(ind2), listlen(ind2));

    Phi2 = MMakeKinshipMatrix(ind2);
    PhiMM = MMakeXLinkedKinshipMatrixBySex(ind2, S_MALE, S_MALE);
    PhiFF = MMakeXLinkedKinshipMatrixBySex(ind2, S_FEMALE, S_FEMALE)*2;
    PhiMF = MMakeXLinkedKinshipMatrixBySex(ind2, S_MALE, S_FEMALE)*2;

    if (Dominance)
      Delta7 = MMakeDelta7Matrix(indlist, Phi2);
    if (Household)
      HH = MakeIdentMatrix(nFamSize);

    Phi2 = Phi2*2;

    FreeIndividualList(ind2);
    
    nSize[i] = nFamSize;
    y[i].Resize(nFamSize,1);
    x[i].Resize(nFamSize,nCovar);

    // Prepare the matrices
    for (j = 0; j< nVC; j++) {
      VC[i+j*nPedigrees].Resize(nFamSize, nFamSize);
      VC[i+j*nPedigrees] = 0.0;
    }

    currentVC = 0;
    // Residual Variance, \sigma_I^2
    for (j = 0; j< VC[i].Rows(); j++) {
      VC[i](j+1,j+1) = 1;
    }
    currentVC++;

    // if (additive)
    VC[i+currentVC*nPedigrees] = Phi2.Subset(keep[i],keep[i]);
    currentVC++;

    if (Dominance) {
      VC[i+currentVC*nPedigrees] = Delta7.Subset(keep[i],keep[i]);
      currentVC++;
    }

    if (Household) {
      VC[i+currentVC*nPedigrees] = 1;
      currentVC++;
    }

    VC[i+currentVC*nPedigrees] = PhiMM.Subset(keep[i],keep[i]);
    currentVC++;
    VC[i+currentVC*nPedigrees] = PhiFF.Subset(keep[i],keep[i]);
    currentVC++;
    VC[i+currentVC*nPedigrees] = PhiMF.Subset(keep[i],keep[i]);
    currentVC++;


    // if (IBD)
    //    VC[i+currentVC*nPedigrees] = IBD[i].Subset(keep[i],keep[i]);

    x[i] = X[i].Subset(keep[i]);

    j = 0;
    k = 0;
    for (ind = indlist; ind; ind = ind->next) {
      if (!strcmpl(ind->pedigree, pedi->name)) {
	k++;
	// Check that no required covariates are missing
	if (keep[i](k,1)) {
	  j++;
	  y[i](j,1) = trait(ind, traitnum);
	}
      }
    }
    //    x[i].Print();
    //    y[i].Print();
  }

  MATRIX start(nVC,1);
  MATRIX beta(nCovar, 1);
  MATRIX *ibdres;

  //  printf("------\n");
  //  x[0].Print();
  //  printf("------\n");

  // Should now go though the list of IBD files.
  // Add the VC's
  // Do analysis
  //

  // Calculate the null hypothesis
  start.Resize(nVC-3,1);
  beta = 1;
  start = 1;
  start(5,1) = 0;
  NullHyp = MaximizeMixedModelX(nPedigrees, nSize, y, x, nVC-3, VC, start, beta, METHOD_ML, 1, 0);

  //  beta.Print();

  start.Resize(nVC,1);
  start = 1;

  //  return 0;

  printf("Pos      LL         LOD\n");
  printf("--------------------------\n");

  // Go through all MIBD files in the right order
  chromlength = (int) (100*ChromosomeLength());

  bestpos = 0;
  bestlod = -1;

  // Go through all possible files at 1 cM 
  // WARNING: This should really be based on a sorted result of the directory listing
  for (pos = 0; pos <= chromlength; pos +=2) {

    // Verify that the MIBD matrices exist for the locus
    sprintf(buf, "mibd-mm.X.%d", pos);    
    if ((F = fopen (buf,"r"))) {
      fclose(F);
    }
    else 
      continue;
    sprintf(buf, "mibd-ff.X.%d", pos);    
    if ((F = fopen (buf,"r"))) {
      fclose(F);
    }
    else 
      continue;
    sprintf(buf, "mibd-mf.X.%d", pos);    
    if ((F = fopen (buf,"r"))) {
      fclose(F);
    }
    else 
      continue;

    sprintf(buf, "mibd-mm.X.%d", pos);
    // Read the three IBD files
    ibdres = ReadPedipetIBDFile(indlist, buf);
    // Insert it in the dataset
    for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
      VC[i+currentVC*nPedigrees] = ibdres[i].Subset(keep[i],keep[i]);
    }
    delete[] ibdres;
    currentVC++;  

    sprintf(buf, "mibd-ff.X.%d", pos);    
    ibdres = ReadPedipetIBDFile(indlist, buf);
    // Insert it in the dataset
    for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
      VC[i+currentVC*nPedigrees] = ibdres[i].Subset(keep[i],keep[i]);
    }
    delete[] ibdres;
    currentVC++;
    
    sprintf(buf, "mibd-mf.X.%d", pos);
    ibdres = ReadPedipetIBDFile(indlist, buf);
    // Insert it in the dataset
    for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
      VC[i+currentVC*nPedigrees] = ibdres[i].Subset(keep[i],keep[i]);      
    }
    delete[] ibdres;
    
    start = 1;
    start(2,1) = .2;
    start(3,1) = .1;
    start(4,1) = .1;
    start(5,1) = 0.05;
    start(6,1) = 1;
    start(7,1) = .4;
    start(8,1) = .2;
    beta = 0;
    beta(2,1) = .5;
    StartHyp = MaximizeMixedModelX(nPedigrees, nSize, y, x, nVC,  VC, start, beta, METHOD_ML, 1, 0);
    lldiff = max(0,-2*(NullHyp-StartHyp));
    pval = .25*pChi2(1, lldiff) + .5*pChi2(2, lldiff) + .25*pChi2(3, lldiff);

    // Check dette
    // Bør nok egentlig være pChi2Inv(1,2*pval); eller hvad??
    if (pval>0)
      lldiff = pChi2Inv(1,pval);
    else {
      lldiff = 62;
    }

    printf("%-3d   %5.3f     %6.4f\n", pos, StartHyp, max(0,lldiff/(2*log(10))));
    fflush(stdout);
    
    // Improved LOD score
    if (bestlod<max(0,lldiff/(2*log(10)))) {
      bestlod = max(0,lldiff/(2*log(10)));
      bestpos = pos;
    }

    currentVC -= 2;

    fflush(stdout);
    
  }
  printf("\n    *** Highest LOD in pass ? was %6.4f at Chrom X Loc %d\n\n", bestlod, bestpos);

  freelist(Pedigrees);
  fflush(stdout);

  delete[] y;
  delete[] keep;
  delete[] x;
  delete[] VC;

  return 0;
}




/**********************************************

  Does singlepoint analysis for X-linked loci
  Goes through the mibd files and 
			       
***********************************************/

double SingleXMarkerVC(individual *indlist, int traitnum,  MATRIX IBD[], MATRIX X[], int nCovar, int Dominance, int Household, int *covnum, int markernum)
{
  int nFamSize;
  int i, j, k;
  int nVC, currentVC=0; 
  double StartHyp, NullHyp, lldiff, pval;
  char buf[20];
  FILE *F;

  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees;
  quanttrait *qttrait;
  MATRIX *y, *x, *VC, *keep;
  MATRIX Phi2, Delta7, HH, PhiMM, PhiFF, PhiMF;

  // Figure out the number of variance components
  nVC = 8;      // I 2*Phi + 3 times Phi_X + 3 times IBD_X
  if (Dominance)
    nVC++;      // I 2*Phi Delta7 IBD
  if (Household)
    nVC++;      // I 2*Phi Delta7 HH IBD      - If not Dom then just remove Delta from this comment to get the order

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  // Sets noprob of zero to missing
  for (ind=indlist; ind; ind=ind->next) {
    if (trait(ind,3) == 0) {
      qttrait = traitvalue(ind, 3);      
      qttrait->sysmiss = 1;      
    }
  }

  // Prepares matrices for input to maximize
  y  = new MATRIX[nPedigrees];
  x  = new MATRIX[nPedigrees];
  keep = new MATRIX[nPedigrees];
  VC = new MATRIX[nPedigrees*nVC];

  int nSize[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    nFamSize = 0;

    // Select the current pedigree
    ind2 = SelectPedigree(pedi->name, indlist);

    int nPedSize = listlen(ind2);
    keep[i].Resize(nPedSize,1);

    // Figure out which individuals to keep
    // Only checks for response and NOT covariates
    keep[i] = CompleteCases(ind2, traitnum);
    for (j = 2; j<nCovar; j++) {
      keep[i] = Hademard(keep[i], CompleteCases(ind2, covnum[j]));
    }

    // The actual number of complete individuals
    int nFamSize = (int) keep[i].Sum();

    Phi2.Resize(listlen(ind2), listlen(ind2));
    PhiMM.Resize(listlen(ind2), listlen(ind2));
    PhiFF.Resize(listlen(ind2), listlen(ind2));
    PhiMF.Resize(listlen(ind2), listlen(ind2));

    Phi2 = MMakeKinshipMatrix(ind2);
    PhiMM = MMakeXLinkedKinshipMatrixBySex(ind2, S_MALE, S_MALE);
    PhiFF = MMakeXLinkedKinshipMatrixBySex(ind2, S_FEMALE, S_FEMALE)*2;
    PhiMF = MMakeXLinkedKinshipMatrixBySex(ind2, S_MALE, S_FEMALE)*2;

    if (Dominance)
      Delta7 = MMakeDelta7Matrix(indlist, Phi2);
    if (Household)
      HH = MakeIdentMatrix(nFamSize);

    Phi2 = Phi2*2;

    FreeIndividualList(ind2);
    
    nSize[i] = nFamSize;
    y[i].Resize(nFamSize,1);
    x[i].Resize(nFamSize,nCovar);

    // Prepare the matrices
    for (j = 0; j< nVC; j++) {
      VC[i+j*nPedigrees].Resize(nFamSize, nFamSize);
      VC[i+j*nPedigrees] = 0.0;
    }

    currentVC = 0;
    // Residual Variance, \sigma_I^2
    for (j = 0; j< VC[i].Rows(); j++) {
      VC[i](j+1,j+1) = 1;
    }
    currentVC++;

    // if (additive)
    VC[i+currentVC*nPedigrees] = Phi2.Subset(keep[i],keep[i]);
    currentVC++;

    if (Dominance) {
      VC[i+currentVC*nPedigrees] = Delta7.Subset(keep[i],keep[i]);
      currentVC++;
    }

    if (Household) {
      VC[i+currentVC*nPedigrees] = 1;
      currentVC++;
    }

    VC[i+currentVC*nPedigrees] = PhiMM.Subset(keep[i],keep[i]);
    currentVC++;
    VC[i+currentVC*nPedigrees] = PhiFF.Subset(keep[i],keep[i]);
    currentVC++;
    VC[i+currentVC*nPedigrees] = PhiMF.Subset(keep[i],keep[i]);
    currentVC++;


    x[i] = X[i].Subset(keep[i]);

    j = 0;
    k = 0;
    for (ind = indlist; ind; ind = ind->next) {
      if (!strcmpl(ind->pedigree, pedi->name)) {
	k++;
	// Check that no required covariates are missing
	if (keep[i](k,1)) {
	  j++;
	  y[i](j,1) = trait(ind, traitnum);
	}
      }
    }
    //    x[i].Print();
    //    y[i].Print();
  }

  MATRIX start(nVC,1);
  MATRIX beta(nCovar, 1);
  MATRIX *ibdres;

  //  printf("------\n");
  //  x[0].Print();
  //  printf("------\n");

  // Should now go though the list of IBD files.
  // Add the VC's
  // Do analysis
  //

  // Calculate the null hypothesis
  start.Resize(nVC-3,1);
  beta = 1;
  start = 1;
  start(5,1) = 0;
  NullHyp = MaximizeMixedModelX(nPedigrees, nSize, y, x, nVC-3, VC, start, beta, METHOD_ML, 1, 0);

  //  beta.Print();

  start.Resize(nVC,1);
  start = 1;

  // New stuff here

  // Verify that the MIBD matrices exist for the locus
  sprintf(buf, "marker%d.X.mm.ibd", markernum);
  if ((F = fopen (buf,"r"))) {
    fclose(F);
  }
  else {
    sprintf(buf, "Marker IBD file marker%d.X.mm.ibd not found.\n", markernum);
    printf("ERROR: %s",buf);
    WriteErrorMsg(buf);
  }

  sprintf(buf, "marker%d.X.ff.ibd", markernum);    
  if ((F = fopen (buf,"r"))) {
    fclose(F);
  }
  else {
    sprintf(buf, "Marker IBD file marker%d.X.ff.ibd not found.\n", markernum);
    printf("ERROR: %s",buf);
    WriteErrorMsg(buf);
  }

  sprintf(buf, "marker%d.X.ff.ibd", markernum);    
  if ((F = fopen (buf,"r"))) {
    fclose(F);
  }
  else {
    sprintf(buf, "Marker IBD file marker%d.X.mf.ibd not found.\n", markernum);
    printf("ERROR: %s",buf);
    WriteErrorMsg(buf);
  }

  sprintf(buf, "marker%d.X.mm.ibd", markernum);
  // Read the three IBD files
  ibdres = ReadPedipetIBDFile(indlist, buf);
  // Insert it in the dataset
  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    VC[i+currentVC*nPedigrees] = ibdres[i].Subset(keep[i],keep[i]);
  }
  delete[] ibdres;
  currentVC++;  
  
  sprintf(buf, "marker%d.X.ff.ibd", markernum);    
  ibdres = ReadPedipetIBDFile(indlist, buf);
  // Insert it in the dataset
  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    VC[i+currentVC*nPedigrees] = ibdres[i].Subset(keep[i],keep[i]);
  }
  delete[] ibdres;
  currentVC++;
  
  sprintf(buf, "marker%d.X.mf.ibd", markernum);
  ibdres = ReadPedipetIBDFile(indlist, buf);
  // Insert it in the dataset
  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    VC[i+currentVC*nPedigrees] = ibdres[i].Subset(keep[i],keep[i]);      
  }
  delete[] ibdres;
  currentVC++;  

  start = 1;
  start(2,1) = .2;
  start(3,1) = .1;
  start(4,1) = .1;
  start(5,1) = 0.05;
  start(6,1) = 1;
  start(7,1) = .4;
  start(8,1) = .2;

  beta = 1;
  //  beta(2,1) = .5;

  StartHyp = MaximizeMixedModelX(nPedigrees, nSize, y, x, nVC,  VC, start, beta, METHOD_ML, 1, 1);
  lldiff = max(0,-2*(NullHyp-StartHyp));
  pval = .25*pChi2(1, lldiff) + .5*pChi2(2, lldiff) + .25*pChi2(3, lldiff);

  // Check dette
  // Bør nok egentlig være pChi2Inv(1,2*pval); eller hvad??
  if (pval>0)
    lldiff = pChi2Inv(1,pval);
  else {
    lldiff = 62;
  }

  printf("Marker %-3d   %5.3f     %6.4f\n", markernum, StartHyp, max(0,lldiff/(2*log(10))));
    
  
  freelist(Pedigrees);
  fflush(stdout);

  delete[] y;
  delete[] keep;
  delete[] x;
  delete[] VC;

  return 0;
}



double MaximizePolygenicVC(individual *indlist, int traitnum,  MATRIX X[], int nCovar, int Dominance, int Household, int *covnum)
{
  int nFamSize;
  int i, j, k;
  int nVC, currentVC; 
  double NullHyp;


  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees;
  quanttrait *qttrait;
  MATRIX *y, *x, *VC, *keep;
  MATRIX Phi2, Delta7, HH, PhiMM, PhiFF, PhiMF;

  // Figure out the number of variance components
  nVC = 2;      // I 2*Phi 
  if (Dominance)
    nVC++;      // I 2*Phi Delta7 IBD
  if (Household)
    nVC++;      // I 2*Phi Delta7 HH IBD      - If not Dom then just remove Delta from this comment to get the order

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  // Sets noprob of zero to missing
  for (ind=indlist; ind; ind=ind->next) {
    if (trait(ind,3) == 0) {
      qttrait = traitvalue(ind, 3);      
      qttrait->sysmiss = 1;      
    }
  }

  // Prepares matrices for input to maximize
  y  = new MATRIX[nPedigrees];
  x  = new MATRIX[nPedigrees];
  keep = new MATRIX[nPedigrees];
  VC = new MATRIX[nPedigrees*nVC];

  int nSize[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    nFamSize = 0;

    // Select the current pedigree
    ind2 = SelectPedigree(pedi->name, indlist);

    int nPedSize = listlen(ind2);
    keep[i].Resize(nPedSize,1);

    // Figure out which individuals to keep
    // Only checks for response and NOT covariates
    keep[i] = CompleteCases(ind2, traitnum);
    for (j = 2; j<nCovar; j++) {
      keep[i] = Hademard(keep[i], CompleteCases(ind2, covnum[j]));
    }

    // The actual number of complete individuals
    int nFamSize = (int) keep[i].Sum();

    Phi2.Resize(listlen(ind2), listlen(ind2));
    Phi2 = MMakeKinshipMatrix(ind2);

    if (Dominance)
      Delta7 = MMakeDelta7Matrix(indlist, Phi2);
    if (Household)
      HH = MakeIdentMatrix(nFamSize);

    Phi2 = Phi2*2;

    FreeIndividualList(ind2);
    
    nSize[i] = nFamSize;
    y[i].Resize(nFamSize,1);
    x[i].Resize(nFamSize,nCovar);

    // Prepare the matrices
    for (j = 0; j< nVC; j++) {
      VC[i+j*nPedigrees].Resize(nFamSize, nFamSize);
      VC[i+j*nPedigrees] = 0.0;
    }

    currentVC = 0;
    // Residual Variance, \sigma_I^2
    for (j = 0; j< VC[i].Rows(); j++) {
      VC[i](j+1,j+1) = 1;
    }
    currentVC++;

    // if (additive)
    VC[i+currentVC*nPedigrees] = Phi2.Subset(keep[i],keep[i]);
    currentVC++;

    if (Dominance) {
      VC[i+currentVC*nPedigrees] = Delta7.Subset(keep[i],keep[i]);
      currentVC++;
    }

    if (Household) {
      VC[i+currentVC*nPedigrees] = 1;
      currentVC++;
    }

    x[i] = X[i].Subset(keep[i]);

    j = 0;
    k = 0;
    for (ind = indlist; ind; ind = ind->next) {
      if (!strcmpl(ind->pedigree, pedi->name)) {
	k++;
	// Check that no required covariates are missing
	if (keep[i](k,1)) {
	  j++;
	  y[i](j,1) = trait(ind, traitnum);
	}
      }
    }
    //    x[i].Print();
    //    y[i].Print();
  }

  MATRIX start(nVC,1);
  MATRIX beta(nCovar, 1);

  // Calculate the null hypothesis
  beta = 1;
  start = 1;
  NullHyp = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC, VC, start, beta, METHOD_ML, 1, 1);

  printf("Log-likelihood:   %6.4f\n", NullHyp);


  /*
  // Check dette
  // Bør nok egentlig være pChi2Inv(1,2*pval); eller hvad??
  if (pval>0)
    lldiff = pChi2Inv(1,pval);
  else {
    lldiff = 62;
  }

  printf("%-3d   %5.3f     %6.4f\n", pos, StartHyp, max(0,lldiff/(2*log(10))));
  fflush(stdout);
  
  // Improved LOD score
  if (bestlod<max(0,lldiff/(2*log(10)))) {
    bestlod = max(0,lldiff/(2*log(10)));
    bestpos = pos;
  }
  */
  
  fflush(stdout);
  

  freelist(Pedigrees);
  fflush(stdout);

  delete[] y;
  delete[] keep;
  delete[] x;
  delete[] VC;

  return 0;
}




double MaximizePolygenicXVC(individual *indlist, int traitnum,  MATRIX X[], int nCovar, int Dominance, int Household, int *covnum)
{
  int nFamSize;
  int i, j, k;
  int nVC, currentVC; 
  double NullHyp;

  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees;
  quanttrait *qttrait;
  MATRIX *y, *x, *VC, *keep;
  MATRIX Phi2, Delta7, HH, PhiMM, PhiFF, PhiMF;

  // Figure out the number of variance components
  nVC = 5;      // I 2*Phi + 3 times Phi_X + 3 times IBD_X
  if (Dominance)
    nVC++;      // I 2*Phi Delta7 IBD
  if (Household)
    nVC++;      // I 2*Phi Delta7 HH IBD      - If not Dom then just remove Delta from this comment to get the order

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  // Sets noprob of zero to missing
  for (ind=indlist; ind; ind=ind->next) {
    if (trait(ind,3) == 0) {
      qttrait = traitvalue(ind, 3);      
      qttrait->sysmiss = 1;      
    }
  }

  // Prepares matrices for input to maximize
  y  = new MATRIX[nPedigrees];
  x  = new MATRIX[nPedigrees];
  keep = new MATRIX[nPedigrees];
  VC = new MATRIX[nPedigrees*nVC];

  int nSize[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    nFamSize = 0;

    // Select the current pedigree
    ind2 = SelectPedigree(pedi->name, indlist);

    int nPedSize = listlen(ind2);
    keep[i].Resize(nPedSize,1);

    // Figure out which individuals to keep
    // Only checks for response and NOT covariates
    keep[i] = CompleteCases(ind2, traitnum);
    for (j = 2; j<nCovar; j++) {
      keep[i] = Hademard(keep[i], CompleteCases(ind2, covnum[j]));
    }

    // The actual number of complete individuals
    int nFamSize = (int) keep[i].Sum();

    Phi2.Resize(listlen(ind2), listlen(ind2));
    PhiMM.Resize(listlen(ind2), listlen(ind2));
    PhiFF.Resize(listlen(ind2), listlen(ind2));
    PhiMF.Resize(listlen(ind2), listlen(ind2));

    Phi2 = MMakeKinshipMatrix(ind2);
    PhiMM = MMakeXLinkedKinshipMatrixBySex(ind2, S_MALE, S_MALE);
    PhiFF = MMakeXLinkedKinshipMatrixBySex(ind2, S_FEMALE, S_FEMALE)*2;
    PhiMF = MMakeXLinkedKinshipMatrixBySex(ind2, S_MALE, S_FEMALE)*2;

    if (Dominance)
      Delta7 = MMakeDelta7Matrix(indlist, Phi2);
    if (Household)
      HH = MakeIdentMatrix(nFamSize);

    Phi2 = Phi2*2;

    FreeIndividualList(ind2);
    
    nSize[i] = nFamSize;
    y[i].Resize(nFamSize,1);
    x[i].Resize(nFamSize,nCovar);

    // Prepare the matrices
    for (j = 0; j< nVC; j++) {
      VC[i+j*nPedigrees].Resize(nFamSize, nFamSize);
      VC[i+j*nPedigrees] = 0.0;
    }

    currentVC = 0;
    // Residual Variance, \sigma_I^2
    for (j = 0; j< VC[i].Rows(); j++) {
      VC[i](j+1,j+1) = 1;
    }
    currentVC++;

    // if (additive)
    VC[i+currentVC*nPedigrees] = Phi2.Subset(keep[i],keep[i]);
    currentVC++;

    if (Dominance) {
      VC[i+currentVC*nPedigrees] = Delta7.Subset(keep[i],keep[i]);
      currentVC++;
    }

    if (Household) {
      VC[i+currentVC*nPedigrees] = 1;
      currentVC++;
    }

    VC[i+currentVC*nPedigrees] = PhiMM.Subset(keep[i],keep[i]);
    currentVC++;
    VC[i+currentVC*nPedigrees] = PhiFF.Subset(keep[i],keep[i]);
    currentVC++;
    VC[i+currentVC*nPedigrees] = PhiMF.Subset(keep[i],keep[i]);
    currentVC++;

    x[i] = X[i].Subset(keep[i]);

    j = 0;
    k = 0;
    for (ind = indlist; ind; ind = ind->next) {
      if (!strcmpl(ind->pedigree, pedi->name)) {
	k++;
	// Check that no required covariates are missing
	if (keep[i](k,1)) {
	  j++;
	  y[i](j,1) = trait(ind, traitnum);
	}
      }
    }
    //    x[i].Print();
    //    y[i].Print();
  }

  MATRIX start(nVC,1);
  MATRIX beta(nCovar, 1);

  // Calculate the null hypothesis
  beta = 1;
  start = 1;
  start(5,1) = 0;
  NullHyp = MaximizeMixedModelX(nPedigrees, nSize, y, x, nVC, VC, start, beta, METHOD_ML, 1, 1);

  printf("Log-likelihood:   %6.4f\n", NullHyp);
  //  printf("Mean parameters: ");
  //  beta.Print();
  //  printf("Variance component parameters: ");
  //  start.Print();

  
  fflush(stdout);
  

  freelist(Pedigrees);
  fflush(stdout);

  delete[] y;
  delete[] keep;
  delete[] x;
  delete[] VC;

  return 0;
}



double MultiGXE(individual *indlist, int traitnum,  MATRIX IBD[], MATRIX X[], MATRIX XX[], int nCovar, int Dominance, int Household, int *covnum)
{
  int nFamSize;
  int i, j, k, pos, bestpos;
  int nVC, currentVC = 0, chromlength; 
  double StartHyp, NullHyp, LinkHyp, lldiff, bestlod, pval;
  char buf[30];
  FILE *F;

  individual *ind, *ind2;
  namelist *Pedigrees, *pedi;
  int nPedigrees;
  quanttrait *qttrait;
  MATRIX *y, *x, *xx, *VC, *keep;
  MATRIX Phi2, Delta7, HH, PhiMM, PhiFF, PhiMF;

  // Figure out the number of variance components
  nVC = 3;      // I 2*Phi + 3 times Phi_X + 3 times IBD_X
  if (Dominance)
    nVC++;      // I 2*Phi Delta7 IBD
  if (Household)
    nVC++;      // I 2*Phi Delta7 HH IBD      - If not Dom then just remove Delta from this comment to get the order

  // Make a list of the pedigrees in the list of individuals
  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  // Sets noprob of zero to missing
  for (ind=indlist; ind; ind=ind->next) {
    if (trait(ind,3) == 0) {
      qttrait = traitvalue(ind, 3);      
      qttrait->sysmiss = 1;      
    }
  }

  // Prepares matrices for input to maximize
  y  = new MATRIX[nPedigrees];
  x  = new MATRIX[nPedigrees];
  xx  = new MATRIX[nPedigrees];
  keep = new MATRIX[nPedigrees];
  VC = new MATRIX[nPedigrees*nVC];

  int nSize[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    nFamSize = 0;

    // Select the current pedigree
    ind2 = SelectPedigree(pedi->name, indlist);

    int nPedSize = listlen(ind2);
    keep[i].Resize(nPedSize,1);

    // Figure out which individuals to keep
    // Only checks for response and NOT covariates
    keep[i] = CompleteCases(ind2, traitnum);
    for (j = 2; j<nCovar; j++) {
      keep[i] = Hademard(keep[i], CompleteCases(ind2, covnum[j]));
    }

    // The actual number of complete individuals
    int nFamSize = (int) keep[i].Sum();
    Phi2.Resize(listlen(ind2), listlen(ind2));
    Phi2 = MMakeKinshipMatrix(ind2);

    if (Dominance)
      Delta7 = MMakeDelta7Matrix(indlist, Phi2);
    if (Household)
      HH = MakeIdentMatrix(nFamSize);

    Phi2 = Phi2*2;

    FreeIndividualList(ind2);
    
    nSize[i] = nFamSize;
    y[i].Resize(nFamSize,1);
    x[i].Resize(nFamSize,nCovar);
    xx[i].Resize(nFamSize,1);

    // Prepare the matrices
    for (j = 0; j< nVC; j++) {
      VC[i+j*nPedigrees].Resize(nFamSize, nFamSize);
      VC[i+j*nPedigrees] = 0.0;
    }

    currentVC = 0;
    // Residual Variance, \sigma_I^2
    for (j = 0; j< VC[i].Rows(); j++) {
      VC[i](j+1,j+1) = 1;
    }
    currentVC++;

    // if (additive)
    VC[i+currentVC*nPedigrees] = Phi2.Subset(keep[i],keep[i]);
    currentVC++;

    if (Dominance) {
      VC[i+currentVC*nPedigrees] = Delta7.Subset(keep[i],keep[i]);
      currentVC++;
    }

    if (Household) {
      VC[i+currentVC*nPedigrees] = 1;
      currentVC++;
    }

    x[i] = X[i].Subset(keep[i]);
    xx[i] = XX[i].Subset(keep[i]);

    j = 0;
    k = 0;
    for (ind = indlist; ind; ind = ind->next) {
      if (!strcmpl(ind->pedigree, pedi->name)) {
	k++;
	// Check that no required covariates are missing
	if (keep[i](k,1)) {
	  j++;
	  y[i](j,1) = trait(ind, traitnum);
	}
      }
    }
    //    x[i].Print();
    //    y[i].Print();
  }

  MATRIX start(nVC,1);
  MATRIX beta(nCovar, 1);
  MATRIX *ibdres;

  // Should now go though the list of IBD files.
  // Add the VC's
  // Do analysis
  //

  // Calculate the null hypothesis
  start.Resize(nVC-1,1);
  beta = 1;
  start = 1;
  NullHyp = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC-1, VC, start, beta, METHOD_ML, 1, 0);
  start.Resize(nVC,1);
  start = 1;

  printf("Pos      LL         LOD(alm)   LOD(GxE)\n");
  printf("-----------------------------------------\n");

  // Go through all MIBD files in the right order
  chromlength = (int) (100*ChromosomeLength());

  bestpos = 0;
  bestlod = -1;

  chromlength = 400;
  // Go through all possible files at 1 cM 
  // WARNING: This should really be based on a sorted result of the directory listing
  for (pos = 0; pos <= chromlength; pos +=2) {
    // Verify that the MIBD matrices exist for the locus
    sprintf(buf, "ibd/mibd.%d.%d.gz", nChromosome, pos);    
    if ((F = fopen (buf,"r"))) {
      fclose(F);
    }
    else 
      continue;

    // Read the three IBD files
    //    ibdres = ReadSolarIBDFile(indlist, buf, "pedindex.out");
    ibdres = ReadSolarIBDFile(indlist, buf, "pedindex.out");
    // Insert it in the dataset

    for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
      VC[i+currentVC*nPedigrees] = ibdres[i].Subset(keep[i],keep[i]);
    }
    delete[] ibdres;
    //    currentVC++;  
    
    start = 1;
    beta = 1;

    StartHyp = MaximizeMixedGXEModel(nPedigrees, nSize, y, x, nVC, VC, 1, xx, start, beta, METHOD_ML, 1, 0, 1.0);
    LinkHyp = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC, VC, start, beta, METHOD_ML, 1, 0);
    //    StartHyp = MaximizeMixedModel(nPedigrees, nSize, y, x, nVC,  VC, start, beta, METHOD_ML, 1, 1);

    printf("%-3d   %5.3f     %6.4f     ", pos, StartHyp, max(0,max(0,-2*(NullHyp-LinkHyp))/(2*log(10))));

    lldiff = max(0,-2*(NullHyp-StartHyp));
    pval = .5*pChi2(1, lldiff) + .5*pChi2(2, lldiff);
    lldiff = pChi2Inv(1,pval);
    printf("%6.4f     %6.4f \n", max(0,lldiff/(2*log(10))), pChi2(1,max(0,-2*(LinkHyp-StartHyp))));


    
    // Improved LOD score
    if (bestlod<max(0,lldiff/(2*log(10)))) {
      bestlod = max(0,lldiff/(2*log(10)));
      bestpos = pos;
    }

  }
  printf("\n    *** Highest LOD in pass ? was %6.4f at Chrom ?? Loc %d\n\n", bestlod, bestpos);

  freelist(Pedigrees);

  delete[] y;
  delete[] keep;
  delete[] x;
  delete[] VC;

  return 0;
}


#undef multip
