/*
 * PEDIPET
 * 
 * bridge.cpp
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



#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

#include "matrix.h"
#include "bridge.h"
#include "estimate.h"

#define multip
#define BYPEDIGREE



//#define	 min(a,b)              ((a) < (b) ? (a) : (b))
//#define	 max(a,b)              ((a) > (b) ? (a) : (b))




individual *SimulateMultiQTLChromosome(int families, int nooff, MATRIX QTLloc, MATRIX QTLadd, MATRIX QTLdom, MATRIX QTLprev, double resadd, double resdom);


extern "C" double CalculateHeritability(individual *indlist, int traitnum, int Dominance)
{
  namelist *Pedigrees, *pedi;
  individual *ind;
  int nPedigrees, i;
  MATRIX *X;

  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  X = new MATRIX[nPedigrees];

  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    ind = SelectPedigree(pedi->name, indlist);
    
    X[i].Resize(listlen(ind),1);
    X[i] = 1;

    FreeIndividualList(ind);
  }

  freelist(Pedigrees);

  Heritability(indlist, traitnum, Dominance, X, 1);

  delete[] X;
  return 0;
}

extern "C" double SingleMarkerIBD(individual *ind, individual *ind2, int markernum);
extern "C" double InverseMapFunction(double dist);
extern "C" double InterMarkerRecomb(int marker1, int marker2);
extern "C" double MarkerDistance(double position, int markernum);


//
// Creates a selection matrix based on a dataset and covariates
//

//
// Only works for non-missing phenotype data
//
//
extern "C" double SinglePointVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov)
{
  MATRIX *IBD, *x;
  individual *ind;
  namelist *Pedigrees, *ped;
  quanttrait *qt;
  char buf[25];
  int nPed, nInd, i, k, nFamSize, fam;


  // Start by checking that the ibd matrix exists
  sprintf(buf, "marker%d.ibd", markernum);
  if (!FileExists(buf))
  {
    printf("Could not find IBD-file for marker %d.\n", markernum);
    return -1;
  }

  // Get the number of pedigrees and individuals
  Pedigrees = MakePedigreeList(indlist);
  nPed = listlen(Pedigrees);
  nInd = listlen(indlist);

  // Initialize IBD scores
  IBD = new MATRIX[nPed];
  x   = new MATRIX[nPed];

  for (i=0; i<nPed; i++) {
    IBD[i].Resize(1,1);
    x[i].Resize(1,1);
  }

  // Read the matrix file
#ifndef BYPEDIGREE
  MATRIX InputMatrix(nInd, nInd);
  InputMatrix.FileReadSymmetric(buf, nInd, nInd);
#endif
#ifdef BYPEDIGREE
  printf("WARNING: USING VC ESTIMATION FOR KNOWN PEDIGREE SIZES ONLY. ERRORS PRONE TO HAPPEN\n");
  MATRIX *InputMatrix;
  InputMatrix = FileReadSymmetricKnownSize(buf, 6);
#endif

  // Initializes local id
  for (ind = indlist, i=1; ind; ind = ind->next, i++)
    ind->localid = i;

  for (ped = Pedigrees, fam=0; ped; ped = ped->next)
  {
    i = 0;

    nFamSize = FamilySize(indlist, ped->name);
    if (nFamSize>0)
    {
      IBD[fam].Resize(nFamSize, nFamSize);
      x[fam].Resize(nFamSize,numcov+1);      
    }

#ifdef BYPEDIGREE
    IBD[fam] = InputMatrix[fam];
#endif

    for (ind = indlist; ind; ind = ind->next)
    {   
      // Found a family member
      if (!strcmpl(ind->pedigree, ped->name))
      {
	i++;

	x[fam](i,1) = 1; // Intercept
        for (k = 1 ; k <= numcov; k++)
	{
	  if (inclcov[k]>0)
	  {
	    qt = (quanttrait*)traitvalue(ind, inclcov[k]);
	    x[fam](i,k+1) = (double) qt->value;
	  }
	  else
	    x[fam](i,k+1) = ind->sex;
	}

#ifndef BYPEDIGREE
	for (ind2 = ind, j = i; ind2; ind2 = ind2->next)
	{
	  // Both from same family
	  if (!strcmpl(ind->pedigree, ind2->pedigree))
	  {

	    IBD[fam](i,j) = InputMatrix(ind->localid, ind2->localid);
	    IBD[fam](j,i) = InputMatrix(ind->localid, ind2->localid);
	    j++;
          }
	}
#endif
      }
    }  
    fam++;
  }

  freelist(Pedigrees);

#ifdef BYPEDIGREE
  delete[] InputMatrix;
#endif

  SingleMarkerVC(indlist, traitnum, 0, IBD, x, numcov+1);

  delete[] IBD;
  delete[] x;

  //  freelist(Pedigrees);

  return 0;
}



extern "C" double SinglePointXVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov) {
  namelist *Pedigrees, *pedi;
  individual *ind, *ind2;
  int nPedigrees, i, j, k;
  MATRIX *X;

  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  X = new MATRIX[nPedigrees];


  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    ind = SelectPedigree(pedi->name, indlist);    
    X[i].Resize(listlen(ind),numcov+1);
    for (ind2 = ind, k=1; ind2; ind2 = ind2->next, k++) {
      X[i](k,1) = 1; // Intercept
      for (j = 1; j <= numcov; j++) {
	if (inclcov[j] == 0) {
	  X[i](k,j+1) = ind2->sex-1;
	}
	else {
	  X[i](k,j+1) = trait(ind2,inclcov[j]);
	}
      }
    }

    FreeIndividualList(ind);
  }
  freelist(Pedigrees);

  SingleXMarkerVC(indlist, traitnum,  NULL, X, numcov+1, 0, 0, inclcov, markernum);

  delete[] X;

  return 0;
}



extern "C" double MultiPointVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov) {
  namelist *Pedigrees, *pedi;
  individual *ind, *ind2;
  int nPedigrees, i, j, k;
  MATRIX *X;


  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  X = new MATRIX[nPedigrees];


  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    ind = SelectPedigree(pedi->name, indlist);    
    X[i].Resize(listlen(ind),numcov+1);
    for (ind2 = ind, k=1; ind2; ind2 = ind2->next, k++) {
      X[i](k,1) = 1; // Intercept
      for (j = 1; j <= numcov; j++) {
	X[i](k,j+1) = trait(ind2,inclcov[j]);
      }
    }

    FreeIndividualList(ind);
  }

  freelist(Pedigrees);

  MultiMarkerVC(indlist, traitnum,  NULL, X, numcov, 0, 0);

  delete[] X;

  return 0;
}



extern "C" double MultiPointXVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov) {
  namelist *Pedigrees, *pedi;
  individual *ind, *ind2;
  int nPedigrees, i, j, k;
  MATRIX *X;

  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  X = new MATRIX[nPedigrees];


  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    ind = SelectPedigree(pedi->name, indlist);    
    X[i].Resize(listlen(ind),numcov+1);
    for (ind2 = ind, k=1; ind2; ind2 = ind2->next, k++) {
      X[i](k,1) = 1; // Intercept
      for (j = 1; j <= numcov; j++) {
	if (inclcov[j] == 0) {
	  X[i](k,j+1) = ind2->sex-1;
	}
	else {
	  X[i](k,j+1) = trait(ind2,inclcov[j]);
	}
      }
    }

    FreeIndividualList(ind);
  }
  freelist(Pedigrees);

  MultiXMarkerVC(indlist, traitnum,  NULL, X, numcov+1, 0, 0, inclcov);

  delete[] X;

  return 0;
}


extern "C" double PolygenicXVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov) {
  namelist *Pedigrees, *pedi;
  individual *ind, *ind2;
  int nPedigrees, i, j, k;
  MATRIX *X;

  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  X = new MATRIX[nPedigrees];


  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    ind = SelectPedigree(pedi->name, indlist);    
    X[i].Resize(listlen(ind),numcov+1);
    for (ind2 = ind, k=1; ind2; ind2 = ind2->next, k++) {
      X[i](k,1) = 1; // Intercept
      for (j = 1; j <= numcov; j++) {
	if (inclcov[j] == 0) {
	  X[i](k,j+1) = ind2->sex-1;
	}
	else {
	  X[i](k,j+1) = trait(ind2,inclcov[j]);
	}
      }
    }

    FreeIndividualList(ind);
  }
  freelist(Pedigrees);

  MaximizePolygenicXVC(indlist, traitnum, X, numcov+1, 0, 0, inclcov);

  delete[] X;

  return 0;
}



extern "C" double PolygenicVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov) {
  namelist *Pedigrees, *pedi;
  individual *ind, *ind2;
  int nPedigrees, i, j, k;
  MATRIX *X;

  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  X = new MATRIX[nPedigrees];


  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    ind = SelectPedigree(pedi->name, indlist);    
    X[i].Resize(listlen(ind),numcov+1);
    for (ind2 = ind, k=1; ind2; ind2 = ind2->next, k++) {
      X[i](k,1) = 1; // Intercept
      for (j = 1; j <= numcov; j++) {
	if (inclcov[j] == 0) {
	  X[i](k,j+1) = ind2->sex-1;
	}
	else {
	  X[i](k,j+1) = trait(ind2,inclcov[j]);
	}
      }
    }

    FreeIndividualList(ind);
  }
  freelist(Pedigrees);

  MaximizePolygenicVC(indlist, traitnum, X, numcov+1, 0, 0, inclcov);

  delete[] X;

  return 0;
}





extern "C" double Paper1Function(individual *indlist, int traitnum)
{
  double ibd[500];
  double K2[500];

  int nPos, nFam;
  double pos, maxdiff, bestpos, difference;

  FILE *F;
  int rpos, rid1, rid2, rFID;
  double rp0, rp1, rp2;
  char minibuf[50];

  nFam = listlen(indlist)/4;

  // Make data for GeneHunter
  ExportLinkagePedigree(indlist,1);
  MakeGHdatafile ();
 
  printf("Running GeneHunter ...");
  system("gh < rungh > /dev/null");
  system("reduce-ibd2");
  system("rm -f ibd_dist.out");
  printf(" done\n");



  bestpos = 0.0;
  maxdiff = -1000000.0;
  #ifndef multip
  nPos = 10;
  {

  #else
  for (nPos = 18; nPos <= 48; nPos = nPos + 2)
  {
    switch (nPos) {
    case 20:
    case 22:
    case 24:
    case 42:
    case 44:
    case 46: continue;
    }
  
  #endif

    pos = (int) nPos;

    sprintf(minibuf, "awk '$1 == %d { print } ' micro.dat > micropos.dat\n", nPos);
    system(minibuf);

    // Read the file
    F = fopen ("micropos.dat","r");

    if (F == 0)
    {
      printf("Could not find file micropos.dat\n");
      exit(1);
    }

    while (fscanf(F, "%d %d %d %d %lf %lf %lf\n", &rpos, &rFID, &rid1, &rid2, &rp0, &rp1, &rp2 ) != EOF)
    {      
      if (rpos == pos)
      {
	ibd[rFID-1] = rp2 + 0.5*rp1;
	K2[rFID-1] = rp2;

      } 
    }
    fclose(F);


    // Detect QTL

    //    printf("------------\nPosition: %d\n", nPos);

    difference = SpecialFunction(indlist, traitnum, 0, ibd, K2);
    //    printf("Difference in likelihoods: %f\n", difference);
    if (difference>maxdiff)
    {
      bestpos = pos;
      maxdiff = difference;
    }
  }

  printf("Best position: %5.3f cM from the left\n", bestpos);
  nPos = (int) bestpos;

  sprintf(minibuf, "awk '$1 == %d { print } ' micro.dat > micropos.dat\n", nPos);
  system(minibuf);

  // Read the file
  F = fopen ("micropos.dat","r");

  if (F == 0)
  {
    printf("Could not find file micropos.dat\n");
    exit(1);
  }

  while (fscanf(F, "%d %d %d %d %lf %lf %lf\n", &rpos, &rFID, &rid1, &rid2, &rp0, &rp1, &rp2 ) != EOF) {      
    if (rpos == nPos)
    {
      ibd[rFID-1] = rp2 + 0.5*rp1;
      K2[rFID-1] = rp2;
    } 
  }
  fclose(F);

  printf("------------\nPosition: %d\n", nPos);

    difference = SpecialFunction(indlist, traitnum, nPos, ibd, K2);

  printf("------------\nPosition: %d\n", nPos);

  return 0;
}


  //
  // Used for old special families
  //


extern "C" double Paper1Function2(individual *indlist, int traitnum)
{
  int i, nPos, pos, nPed;
  int nMarkers = numberofmarkers();
  double maxdiff, bestpos, difference;

  namelist *Pedigrees;
  char minibuf[50];
  FILE *F;

  MATRIX *IBD, *KK2;

  Pedigrees = MakePedigreeList(indlist);
  nPed = listlen(Pedigrees);
  freelist(Pedigrees);

  IBD = new MATRIX[nPed];
  KK2 = new MATRIX[nPed];

  for (i=0; i<nPed; i++) {
    IBD[i].Resize(9,9);
    KK2[i].Resize(9,9);
  }

  MATRIX VARI(nMarkers, nMarkers);              // Holds the V matrix
  MATRIX VARIK2(nMarkers, nMarkers);              // Holds the V matrix for K2
  MATRIX Beta(nMarkers,1), C(nMarkers,1), BetaK2(nMarkers,1), CK2(nMarkers,1);       // Holds weights for the markers

  int rpos, rid1, rid2, rFID;
  double rp0, rp1, rp2;


  bestpos = 0.0;
  maxdiff = -1000000.0;
  #ifndef multip
  nPos = 10;
  {

  #else
  for (nPos = 6; nPos <= 60; nPos = nPos + 2)
  {
  
  #endif

    pos = (int) nPos;

    sprintf(minibuf, "awk '$1 == %d { print } ' micro.dat > micropos.dat\n", pos);
    system(minibuf);

    // Read the file
    F = fopen ("micropos.dat","r");

    if (F == 0)
    {
      printf("Could not find file micropos.dat\n");
      exit(1);
    }

    while (fscanf(F, "%d %d %d %d %lf %lf %lf\n", &rpos, &rFID, &rid1, &rid2, &rp0, &rp1, &rp2 ) != EOF)
    {      
      if (rpos == pos)
      {
	//	printf("Family : %d  %d  %d    %f  %f  %f \n", rFID, rid1, rid2, rp0, rp1, rp2);
	IBD[rFID-1](rid1, rid2) = rp2 + 0.5*rp1;
	IBD[rFID-1](rid2, rid1) = IBD[rFID-1](rid1, rid2);

	KK2[rFID-1](rid1, rid2) = rp2;
	KK2[rFID-1](rid2, rid1) = KK2[rFID-1](rid1, rid2);
      } 
    }
    fclose(F);

    for (i = 0; i<nPed; i++)
    {
      IBD[i](1, 5) = 0.5;
      IBD[i](5, 1) = 0.5;
      IBD[i](1, 6) = 0.5;
      IBD[i](6, 1) = 0.5;

      IBD[i](2, 5) = 0.5;
      IBD[i](5, 2) = 0.5;
      IBD[i](2, 6) = 0.5;
      IBD[i](6, 2) = 0.5;
      IBD[i](2, 7) = 0.5;
      IBD[i](7, 2) = 0.5;
      IBD[i](2, 8) = 0.5;
      IBD[i](8, 2) = 0.5;
 
      IBD[i](3, 7) = 0.5;
      IBD[i](7, 3) = 0.5;
      IBD[i](3, 8) = 0.5;
      IBD[i](8, 3) = 0.5;

      IBD[i](8, 9) = 0.5;
      IBD[i](9, 8) = 0.5;

      IBD[i](4, 9) = 0.5;
      IBD[i](9, 4) = 0.5;
    }   

    // Detect QTL

    printf("------------\nPosition: %d\n", nPos);

    difference = SpecialFunction2(indlist, traitnum, 0, IBD, KK2);
    //    printf("Difference in likelihoods: %f\n", difference);
    if (difference>maxdiff)
    {
      bestpos = pos;
      maxdiff = difference;
    }
  }

  printf("Best position: %5.3f cM from the left\n", bestpos);
  pos = (int) bestpos;


    sprintf(minibuf, "awk '$1 == %d { print } ' micro.dat > micropos.dat\n", pos);
    system(minibuf);

    // Read the file
    F = fopen ("micropos.dat","r");

    if (F == 0)
    {
      printf("Could not find file micropos.dat\n");
      exit(1);
    }

    while (fscanf(F, "%d %d %d %d %lf %lf %lf\n", &rpos, &rFID, &rid1, &rid2, &rp0, &rp1, &rp2 ) != EOF)
    {      
      if (rpos == pos)
      {
	IBD[rFID-1](rid1, rid2) = rp2 + 0.5*rp1;
	IBD[rFID-1](rid2, rid1) = IBD[rFID-1](rid1, rid2);

	KK2[rFID-1](rid1, rid2) = rp2;
	KK2[rFID-1](rid2, rid1) = KK2[rFID-1](rid1, rid2);
      } 
    }
    fclose(F);
    
    for (i = 0; i<nPed; i++)
    {
      IBD[i](1, 5) = 0.5;
      IBD[i](5, 1) = 0.5;
      IBD[i](1, 6) = 0.5;
      IBD[i](6, 1) = 0.5;

      IBD[i](2, 5) = 0.5;
      IBD[i](5, 2) = 0.5;
      IBD[i](2, 6) = 0.5;
      IBD[i](6, 2) = 0.5;
      IBD[i](2, 7) = 0.5;
      IBD[i](7, 2) = 0.5;
      IBD[i](2, 8) = 0.5;
      IBD[i](8, 2) = 0.5;
 
      IBD[i](3, 7) = 0.5;
      IBD[i](7, 3) = 0.5;
      IBD[i](3, 8) = 0.5;
      IBD[i](8, 3) = 0.5;

      IBD[i](8, 9) = 0.5;
      IBD[i](9, 8) = 0.5;

      IBD[i](4, 9) = 0.5;
      IBD[i](9, 4) = 0.5;
    }   

  SpecialFunction2(indlist, traitnum, pos, IBD, KK2);

  delete[] IBD;
  delete[] KK2;

  return 0;
}


extern "C" double GHPiHat(individual *indlist, int traitnum, int markernum)
{
  return(SpecialFunctionPiHat(indlist, traitnum, markernum));
  //  SpecialFunctionSibsPiHat(indlist, traitnum);
}




extern "C" double ComplexMixed(individual *indlist, int traitnum)
{
  int i, nPos, pos, nPed;
  int nMarkers = numberofmarkers();
  double maxdiff, bestpos, difference;

  int nRela = 6;
  i = 0;
  difference = 0;
  //  int nRela = 1;

  namelist *Pedigrees;
  char minibuf[50];
  FILE *F;

  MATRIX *RELAP;

  Pedigrees = MakePedigreeList(indlist);
  nPed = listlen(Pedigrees);
  freelist(Pedigrees);

  RELAP = new MATRIX[nPed];

  for (i = 0 ; i<nPed; i++) {
    RELAP[i].Resize(nRela,3);
  }

  MATRIX VARI(nMarkers, nMarkers);              // Holds the V matrix
  MATRIX VARIK2(nMarkers, nMarkers);              // Holds the V matrix for K2
  MATRIX Beta(nMarkers,1), C(nMarkers,1), BetaK2(nMarkers,1), CK2(nMarkers,1);       // Holds weights for the markers
 
  int rpos, rid1, rid2, rFID, oldrFID;
  double rp0, rp1, rp2;
  double full, nodom, noeff;


  bestpos = 0.0;
  maxdiff = -1000000.0;
  #ifndef multip
  nPos = 34;
  {

  #else
  for (nPos = 34; nPos <= 34; nPos = nPos + 2)
  {
  
  #endif

    pos = (int) nPos;

    sprintf(minibuf, "awk '$1 == %d { print } ' micro.dat > micropos.dat\n", pos);
    system(minibuf);

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

	// Reads the relationships and reduces precision to save a little time

	RELAP[rFID-1](i,1) = rint(max(rp0,0)*100+.00001)/100.0;
	RELAP[rFID-1](i,2) = rint(max(rp1,0)*100+.00001)/100.0;
	RELAP[rFID-1](i,3) = max(1 - RELAP[rFID-1](i,1)-RELAP[rFID-1](i,2),0);
      } 
    }
    fclose(F);

    printf("------------\nPosition: %d\n", nPos);

    full  = SpecialFunction3(indlist, traitnum, 1, RELAP);
    nodom = SpecialFunction3(indlist, traitnum, 2, RELAP);
    noeff = SpecialFunction3(indlist, traitnum, 4, RELAP);
    printf("Mixed : %d %f %f %f\n", nPos, full, nodom, noeff);
    if (difference>maxdiff)
    {
      bestpos = pos;
      maxdiff = difference;
    }
  }

  delete[] RELAP;
  return 0;
}




extern "C" double SibsComplexMixed(individual *indlist, int traitnum, int markernum)
{
  double full, add;

  // trait marker
  printf("Staring special function for marker %d\n", markernum);

  full  = SpecialFunctionMixSibs(indlist, traitnum, markernum, 1);
  add   = SpecialFunctionMixSibs(indlist, traitnum, markernum, 2);

  printf("Full model    : %f\n", full);
  printf("Additive model: %f\n", add);

  printf ("Ending Complex model\n");

  return (full-add);
}



/*
 
  Calculates the multipoint IBD-score for full sibs
 
*/
 
void FullSibsIBD(individual *indlist, char *filename, int step) {
  individual *ind, *ind2;
  double pihat;
  int i, j, k, nPairs, location;
  char buf[100];
  int nMarkers = numberofmarkers();
  int nPers = listlen(indlist);
  MATRIX V(nMarkers,nMarkers), C(nMarkers,1), COk(nMarkers,1), InvV(nMarkers,nMarkers), Beta(nMarkers, 1);
  MATRIX VarMeans(nMarkers,1), K2Means(nMarkers, 1);
  MATRIX PiHatMean(nMarkers,1), PiHatSqr(nMarkers, 1);
  MATRIX K2Mean(nMarkers,1), K2Sqr(nMarkers, 1);


  MATRIX *IBD, *K2, result(nPers, nPers);

  // Calculates the IBD-matrices
  IBD = new MATRIX[numberofmarkers()];
  K2  = new MATRIX[numberofmarkers()];

  for (i=0; i < numberofmarkers(); i++) {
    IBD[i].Resize(nPers, nPers);
    K2[i].Resize(nPers, nPers);	
  }

  nPairs = 0;

  // Calculates the marker variances
  for (ind = indlist, i=1; ind->next; ind = ind->next, i++) {
    for (ind2 = ind->next, j=i+1; ind2; ind2 = ind2->next, j++) {
      if (fullsibs(ind, ind2)) {
	nPairs++;
	for (k = 1; k<= nMarkers; k++) {
	  pihat = SingleMarkerIBD(ind, ind2, k);

	  PiHatMean(k,1) += pihat;
	  PiHatSqr(k,1)  += SQR(pihat);

	  K2Mean(k,1)    += k2;
	  K2Mean(k,1)    += SQR(k2);

	  IBD[k-1](i,j)  = pihat;
	  K2[k-1](i,j)   = k2;

	}
      }
    }
  }

  PiHatMean *= (double) 1.0 / nPairs;

  //  PiHatMean.Print();

  // K2 is actually not used right now, so it will be deleted here
  delete[] K2;
  K2 = NULL;

  printf("HER - 3\n");

  // Have the means.
  for (i = 1; i<= nMarkers; i++) {
    V(i,i) = (PiHatSqr(i,1) - nPairs*SQR(PiHatMean(i,1))) /((double) nPairs-1.0);
    C(i,1) = V(i,i);
  }
  for (i = 1; i< nMarkers; i++) {
    for (j = i+1; j<= nMarkers; j++) {
      V(i,j) = 8*SQR(1-2*InterMarkerRecomb(i,j))*V(i,i)*V(j,j);
      V(j,i) = V(i,j);
    }
  }

  //  V.Print();

  InvV = Inverse(V);

  for (location = 0; location <=(int) (100*ChromosomeLength()) ; location += step) {
    result = 0.0;
    //  printf("Location : %d   %f\n", location, 100*ChromosomeLength());
    for (ind = indlist, i=1; ind->next; ind = ind->next, i++) {
      result(i,i) = 1.0;
      for (ind2 = ind->next, j=i+1; ind2; ind2 = ind2->next, j++) {
	result(j,j) = 1.0;
	// Checks for direct descendant
	if ( ind->father == ind2 || ind->mother == ind2 || ind2->father == ind || ind2->mother == ind ) {
	  result(i,j) = .5;
	  result(j,i) = .5;
	}
	else if (fullsibs(ind, ind2)) {
	  result(i,j) = .5;
	  for (k = 1; k<= nMarkers; k++) {
	    COk(k,1) = C(k,1)*SQR(1-2*InverseMapFunction(fabs(MarkerDistance(location,k)/100.0)));
	  }
	  Beta = InvV*COk;

	  for (k = 1; k<= nMarkers; k++) {
	    result(i,j) += Beta(k,1)*(IBD[k-1](i,j) - PiHatMean(k,1));
	    //	  printf("%d <-> %d : IBD: %f  [%f]   %f==> %f\n", i, j, IBD[k-1](i,j), Beta(k,1), PiHatMean(k,1),result(i,j));
	  }

	  result(i,j) = min(result(i,j),1);
	  result(i,j) = max(result(i,j),0);
	  result(j,i) = result(i,j);
	}
      }
    }
    sprintf(buf, "mibd-%d.pp", location);
    result.FileWriteSymmetric(buf);
  }

  delete[] IBD;
  delete[] K2;
}

int MMakeKinshipMatrix(individual *indlist, char *filename) {
  MATRIX res;
  int nInd;

  nInd = listlen(indlist);

  res.Resize(nInd, nInd);

  //
  // XXX   SHould change the stuff below to so not to multiply by 2
  // 
  if (options->Index[O_XLINKED]) {
    res = MMakeXLinkedKinshipMatrixBySex(indlist, S_MALE, S_MALE);
    res.FileWriteSymmetric("kinship.male");
    res = MMakeXLinkedKinshipMatrixBySex(indlist, S_FEMALE, S_FEMALE);
    res.FileWriteSymmetric("kinship.female");
    res = MMakeXLinkedKinshipMatrixBySex(indlist, S_MALE, S_FEMALE);
    res.FileWriteSymmetric("kinship.malefemale");
    res = MMakeXLinkedKinshipMatrix(indlist);
  }
  else 
    res = MMakeKinshipMatrix(indlist);

  res.FileWriteSymmetric(filename);
  return 0;
}



extern "C" void DataMining() {
  MATRIX QTLadd(3,1), QTLloc(3,1), QTLdom(3,1), QTLprev(3,1);

  QTLadd = 0;
  QTLdom = 0;
  QTLprev = .5;
  QTLloc(1,1) = 30;
  QTLloc(2,1) = 110;
  QTLloc(3,1) = 190;
  QTLadd(1,1) = 1.1155;
  QTLadd(2,1) = 0.8433;
  QTLadd(3,1) = 0.6992;

  individuals = SimulateMultiQTLChromosome(250, 4, QTLloc, QTLadd, QTLdom, QTLprev, 0, 0);

}


extern "C" double MultiPointGXE(individual *indlist, int traitnum, int envnum, int markernum, int numcov, int *inclcov) {
  namelist *Pedigrees, *pedi;
  individual *ind, *ind2;
  int nPedigrees, i, j, k;
  MATRIX *X, *XX;


  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);

  X = new MATRIX[nPedigrees];
  XX = new MATRIX[nPedigrees];


  for (pedi = Pedigrees, i=0; pedi; pedi = pedi->next, i++) {
    ind = SelectPedigree(pedi->name, indlist);    
    X[i].Resize(listlen(ind),numcov+1);
    XX[i].Resize(listlen(ind),1);
    for (ind2 = ind, k=1; ind2; ind2 = ind2->next, k++) {
      X[i](k,1) = 1; // Intercept
      XX[i](k,1) = trait(ind2, envnum); // environment
      for (j = 1; j <= numcov; j++) {
	if (inclcov[j] == 0) {
	  X[i](k,j+1) = ind2->sex-1;
	}
	else {
	  X[i](k,j+1) = trait(ind2,inclcov[j]);
	}
      }
    }

    FreeIndividualList(ind);
  }
  freelist(Pedigrees);

  MultiGXE(indlist, traitnum,  NULL, X, XX, numcov+1, 0, 0, inclcov);

  delete[] X;
  delete[] XX;

  return 0;
}
