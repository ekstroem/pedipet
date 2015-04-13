/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * This file is part of the statistical program Pedipet.
 * Functions to do ML and REML of mixed model
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

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "maximize.h"

// #define debug
#define constrained       // Do constrained maximization (i.e. sigma >= 0)
#define TOLERANCE         0.00001
#define CLOSETOEDGE       0.001       // When is a parameter close to the edge
#define PCLOSETOEDGE      0.001       // When is the mixing parameter close to the edge
#define SCALEEPSILON      0.1
#define MINIMUMSTEPSIZE   -30         // På log-skala
#ifndef PI                // Definition of PI when lacking in compiler
#define PI                3.141592654 
#endif

#ifndef SQR
static double sqrarg;
#define  SQR(a)                ((sqrarg=(a))==0.0 ? 0.0 : sqrarg*sqrarg)
#endif

#define REML

#ifndef min
#define min(a,b)          ((a) < (b) ? (a) : (b))
#define max(a,b)          ((a) > (b) ? (a) : (b))
#endif

#define constrainmu       .05
#define output

//#define newthings

//#define improveconv
#define fix_error
#define MIXTURECONVLAG    8
//#define REALCALC          // If a parameter is close to 0 than CLOSETOEDGE, then set it to 0 for likelihood calculations
#define GEM               // Use general EM
#define obsdatalikelihood
#define diffbeta
#define edgewalk          // If true then walk along edges without trying to get close to 0
#define newllh
#define showvariances

#define GXE
#define AITKEN

char *ConvergenceText[5] = {
  "No convergence",
  "No change in likelihoods",
  "No change in mean parameters",
  "Stepsize too small",
  "Max iterations reached",
};
 

/*

  MaximizeMixedModel maximizes a MM under the constraints, that the VC should be non-negative.
  No missing data should be present in the input data

  Input    : nPedigrees    number of pedigrees
             nPedSize      array describing the size of the nPedigrees pedigrees
             y             array of matrices (each of size nPedSize[i]*1) of responses
             x             array of matrices (each of size nPedSize[i]*p) of design matrices
	     nVC           number of variance components
             VC            array of variance components matrices. 
                           Indexed by families (the first nPedigrees matrices are alle VC_1, 
                           the next nPedigrees are VC_2 etc.)
	     start         a vector of starting values for the variance components
             beta          a vector of starting values for the mean parameters
	     Method        Used  maximization method (see below). 1 = ML
             Constrain     Should we use constrained optimization? 0 = no, otherwise yes
                           Isn't implemented. Allways constrained VC to be non-negative.
         


  Returns the negative log likelihood value at the minimum. 
  Start holds the estimated variance components paraemters

  REML doesn't work, as you can't work on families but have to work with the complete
  variance matrix (P isn't block-diagonal)
             
  (C) Claus Ekstrøm 1999--2000

*/


double MaximizeMixedModel(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, int Method, int Constrain, int PrintInfo)
{
  int i, j, nIter, nFam, nCode, nTotal=0;
  int nMeanParam, convergence;
  double dLogLike, dLogDet, dLogDet2, LogLike, downscale, change;
  double step, newll;

  convergence = 0;

  // Generel matrices
  MATRIX xT;
  MATRIX Omega, InvOmega, P, InvOmega2;
  MATRIX Fisher(nVC, nVC), InvFisher(nVC, nVC);
  MATRIX Deriv(nVC, 1);
  MATRIX theta(nVC, 1), DeltaTheta(nVC, 1), WorkingTheta(nVC, 1);;  
  MATRIX DeltaBeta = beta;
  MATRIX OldBeta = beta;
  MATRIX xOmegax(beta.Rows(), beta.Rows());

  // Matrices for REML estimation
  MATRIX *IOmega = NULL;
  MATRIX NewY, NewX, *NewVC = NULL;


  nMeanParam = x[0].Cols();

#ifdef REML
  if (Method != METHOD_ML && Method != METHOD_REML ) {
    printf("ERROR: Trying to maximize likelihood using non-(RE)ML methods.\nThis is not implemented yet.\n");
    exit(1);
  }
#else
  if (Method != METHOD_ML) {
    printf("ERROR: Trying to maximize likelihood using non-ML methods.\nThis is not implemented yet.\n");
    exit(1);
  }
#endif

  // Should combine all the matrices to one large set
  if (Method == METHOD_REML) {

#ifdef debug
  printf("Preparing REML code\n");
#endif

    NewVC = new MATRIX[nVC];

    nTotal=0;
    for (i = 0; i< nPedigrees; i++) {
      nTotal += nPedSize[i];
    }

    for (i = 0; i< nVC; i++) {
      NewVC[i].Resize(nTotal,nTotal);
      NewVC[i] = CombineMatrices(&VC[i*nPedigrees], nPedigrees);
    }

    NewY.Resize(nTotal,1);
    NewX.Resize(nTotal,beta.Rows());

    NewY = AppendMatrices(y, nPedigrees);
    NewX = AppendMatrices(x, nPedigrees);

    xT.Resize(beta.Rows(), nTotal);
    xT = Transpose(NewX);

    MATRIX tempmat(nTotal, nTotal);
    tempmat = 0;
    for (i = 1; i<=nTotal; i++) {
      tempmat(i,i) = 1;
    }
    NewY = (tempmat - NewX*Inverse(xT*NewX)*xT)*NewY;

    // Holds the pedigree variances below
    IOmega = new MATRIX[nPedigrees];
    for (nFam = 0; nFam < nPedigrees; nFam++){      
      IOmega[nFam].Resize(nPedSize[nFam], nPedSize[nFam]);
    }
    P.Resize(nTotal, nTotal);
    InvOmega.Resize(nTotal, nTotal);
    InvOmega2.Resize(nTotal, nTotal);


  }

#ifdef debug
  printf("Start maximizing\n");
#endif

  // Start checking that everythink looks ok

  theta = start;  

  // Starts iterating
  for (nIter = 1; nIter <= MIXEDMODELMAXITER; nIter++)  {

#ifdef debug
    printf("Iteration: %3d\n", nIter);
#endif
    Fisher = 0.0;
    Deriv  = 0.0;
    DeltaBeta = 0.0;

    dLogLike = 0.0;
    LogLike = 0.0;

    xOmegax = 0.0;
    

    if (Method == METHOD_ML) {
      // This part calculates the first derivative and the Fisher information
      // matrix using Maximum likelihood methods

#ifdef debug
    theta.Print();
#endif


      // Go through each family/pedigree
      // The block diagonal structure is kept when doing ML
      for (nFam = 0; nFam < nPedigrees; nFam++)	{
	// Calculate the current variance for the family
	Omega.Resize(nPedSize[nFam], nPedSize[nFam]);
	InvOmega.Resize(nPedSize[nFam], nPedSize[nFam]);
	P.Resize(nPedSize[nFam], nPedSize[nFam]);
	xT.Resize(x[nFam].Cols(), x[nFam].Rows());
	xT = Transpose(x[nFam]);
	
	Omega = 0.0;

	nCode = nFam*nPedigrees;

	for (i = 0; i < nVC; i++)
	  Omega += VC[nFam+nPedigrees*i]*theta(i+1,1);

	// Inverting Omega_i
	InvOmega = Inverse(Omega, &dLogDet);

	// Add to the log likelihood
	MATRIX mu = (x[nFam] * beta);
	
	xOmegax   += (xT*InvOmega)*x[nFam];
	DeltaBeta += xT*(InvOmega*(y[nFam]-mu));
	LogLike += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));
     
	InvOmega2.Resize(nPedSize[nFam], nPedSize[nFam]);
	MATRIX UsedMatrix(nPedSize[nFam], nPedSize[nFam]);
	UsedMatrix = InvOmega;

	InvOmega2 = UsedMatrix*UsedMatrix;
	// Calculating the derivatives and Fisher scoring matrix
	for (i = 0; i < nVC; i++) {
	  Deriv(i+1,1) += - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
	    + ((Transpose(y[nFam] - mu)*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-mu)))(1,1);

	  for (j = i; j < nVC; j++) {
	    Fisher(i+1,j+1) += Trace(InvOmega2*VC[nFam + nPedigrees*i]*VC[nFam + nPedigrees*j]);
	    Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	  }
	}
      }
    }
    else if (Method == METHOD_REML) {
      // This part calculates the first derivative and the Fisher information
      // matrix using Restricted maximum likelihood methods

      // Calculate the current inverse variance matrix
      dLogDet = 0;
      for (nFam = 0; nFam < nPedigrees; nFam++)	{      
	IOmega[nFam] = 0;
	for (i = 0; i < nVC; i++) {	
	  IOmega[nFam] += VC[nFam+nPedigrees*i]*theta(i+1,1);
	}
	IOmega[nFam] = Inverse(IOmega[nFam], &dLogDet2);
	dLogDet += dLogDet2;
      }
      InvOmega = CombineMatrices(IOmega, nPedigrees);

      // Calculate P
      xOmegax = Inverse((xT*InvOmega)*NewX, &dLogDet2);
      P = InvOmega - ((InvOmega*NewX)*xOmegax)*(xT*InvOmega);   
      InvOmega2 = P*P;

      MATRIX mu = (NewX * beta);

      //	DeltaBeta += xT*(InvOmega*(y[nFam]-mu));
      //	LogLike += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));

      for (i = 0; i < nVC; i++) {
	Deriv(i+1,1) += - Trace(P*NewVC[i]) 
	  + ((Transpose(NewY)*P)*NewVC[i]*(P*(NewY)))(1,1);
	  
	for (j = i; j < nVC; j++) {
	  Fisher(i+1,j+1) += Trace(InvOmega2*NewVC[i]*NewVC[j]);
	  Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	}
      }

      LogLike = - 0.5*(nTotal*log(2*PI) + dLogDet + dLogDet2 + (Transpose(NewY)*P*(NewY))(1,1));
      LogLike -= - 0.5*nMeanParam*log(2*PI);
    }


    // Remember to scale the 1st and 2nd derivatives by 1/2
    Deriv  *= .5;
    Fisher *= .5;

    OldBeta = beta;

    if (Method == METHOD_ML) {
      beta += Inverse(xOmegax)*DeltaBeta;
      DeltaBeta = Inverse(xOmegax)*DeltaBeta;
    }

#ifdef constrained
    // Modifies the diagonal of the Fisher Matrix
    for (i=1; i<=nVC; i++)
      Fisher(i,i) += constrainmu/theta(i,1);

#endif

    // Inverts the matrix
    InvFisher = Inverse(Fisher);

    // Calculates the change in delta
    DeltaTheta = InvFisher*Deriv;
    WorkingTheta = DeltaTheta;

    // This change should be constrained
#ifdef constrained

    downscale = 1.0;  // Maximum steplength to stay positive

    // Walk along edge for parameters close to the edge
    for (i=1; i<=nVC; i++) {
      if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	WorkingTheta(i,1)=0;
      // Now finds the maximum allowed steplength 
      if (WorkingTheta(i,1)<0)
	downscale = min(downscale, (1-SCALEEPSILON)*fabs(theta(i,1)/WorkingTheta(i,1)));
    }
    // Scales the step down
    DeltaTheta = DeltaTheta*downscale;
    WorkingTheta = WorkingTheta*downscale;
#endif

    step = 2.0;
    do  {
      step *= .5;
      WorkingTheta = DeltaTheta*step;

      for (i=1; i<=nVC; i++) {
	if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	  WorkingTheta(i,1)=0;
      }

      newll = 0.0;
      // Calculate the new Omega

      if (Method == METHOD_ML) {
	for (nFam = 0; nFam < nPedigrees; nFam++)	{
	  // Calculate the current variance for the family
	  MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX xT = Transpose(x[nFam]);
	  
	  Omega = 0.0;
	  
	  nCode = nFam*nPedigrees;
	  
	  for (i = 0; i < nVC; i++)
	    Omega += VC[nFam+nPedigrees*i]*(theta(i+1,1) + WorkingTheta(i+1,1));
	  
	  // Inverting Omega_i
	  InvOmega = Inverse(Omega, &dLogDet);
	  
	  // Add to the log likelihood
#ifdef fix_error
	  MATRIX mu = (x[nFam] * beta);
#else
	  MATRIX mu = (x[nFam] * OldBeta);
#endif
	  newll += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));
	}
	  
      } else if (Method == METHOD_REML) {
	// Calculate the current inverse variance matrix
	dLogDet = 0;
	for (nFam = 0; nFam < nPedigrees; nFam++)	{      
	  IOmega[nFam] = 0;
	  for (i = 0; i < nVC; i++) {	
	    IOmega[nFam] += VC[nFam+nPedigrees*i]*(theta(i+1,1) + WorkingTheta(i+1,1));
	  }
	  IOmega[nFam] = Inverse(IOmega[nFam], &dLogDet2);
	  dLogDet += dLogDet2;
	}
	InvOmega = CombineMatrices(IOmega, nPedigrees);
	xOmegax = Inverse((xT*InvOmega)*NewX, &dLogDet2);
	P = InvOmega - ((InvOmega*NewX)*xOmegax)*(xT*InvOmega);   

	newll = - 0.5*(nTotal*log(2*PI) + dLogDet + dLogDet2 + (Transpose(NewY)*P*(NewY))(1,1));
    	newll -= -0.5*nMeanParam*log(2*PI);	



      }    
    }
    while (newll < LogLike && fabs(newll-LogLike)>TOLERANCE*downscale*step);

    // Modifies Theta
    theta += WorkingTheta;

    // XXX
    // Should also check the change in beta for ML 
    // Exiting due to small increas in LL or theta
    change = 0.0;

    for (i=1; i<= nVC; i++)
      change += fabs(WorkingTheta(i,1));

    // if ML and no change in mean parameters
    if (Method == METHOD_ML) {
      for (i=1; i<= beta.Rows(); i++)
	change += fabs(DeltaBeta(i,1));
    }


    // Should perhaps multiply by step also?
    // If yes, then same criteria should be used above in the while loop
    if (fabs(newll-LogLike)<TOLERANCE*downscale*step || change<TOLERANCE*downscale*step) {
      if (fabs(newll-LogLike)<TOLERANCE*downscale*step) {
	convergence = 1;
      }
      else {
	convergence = 2;
      }

      // Calculating the mean parameter estimates based on the variances
      if (Method == METHOD_REML) {
	NewY = AppendMatrices(y, nPedigrees);
	beta = xOmegax*(xT*InvOmega)*NewY;

	// Fixes xOmegax
	// Above, when using REML xOmegax is the INVERSE of xomega x, but for
	// ML is is xomegax
	// Therefore, change it back sp the right variance of the
	// mean is printed

	xOmegax = Inverse(xOmegax);
	
      }

      RES1 = theta(1,1);
      if (theta.Rows()>=2)       
	RES2 = theta(2,1);
      if (theta.Rows()>=3) 
	RES3 = theta(3,1);
      break;
    }
  }
#ifdef output
  if (PrintInfo != 0) {
    printf("Convergence reached after %d iterations (%s)\n",nIter, ConvergenceText[convergence]);
    printf("Estimating using %s.\n", Method==METHOD_ML ? "ML" : "REML");
    printf("Estimated variance components:\n");
    theta.Print();
    printf("Estimated mean parameters:\n");
    beta.Print();
    printf("Log likelihood at maximum: %6.3f\n", newll);
#ifdef showvariances	
    printf("Inverse Fisher information matrix for VC:\n");
    InvFisher.Print();
    printf("Variances of the mean parameters:\n");
    Inverse(xOmegax).Print();
#endif
  }
    
#endif

  fflush(stdout);

  if (Method == METHOD_REML) {
    delete[] IOmega;
    delete[] NewVC;
  }


  return (newll);
}

/*********************************************************
 *
 * Special version for X-linked markers
 *
 *********************************************************/


double MaximizeMixedModelX(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, int Method, int Constrain, int PrintInfo)
{
  int i, j, nIter, nFam, nCode, nTotal=0;
  int nMeanParam, convergence;
  double dLogLike, dLogDet, dLogDet2, LogLike, downscale, change;
  double step, newll, edge;

  convergence = 0;

  // Generel matrices
  MATRIX xT;
  MATRIX Omega, InvOmega, P, InvOmega2;
  MATRIX Fisher(nVC, nVC), InvFisher(nVC, nVC);
  MATRIX Deriv(nVC, 1);
  MATRIX theta(nVC, 1), DeltaTheta(nVC, 1), WorkingTheta(nVC, 1);;  
  MATRIX DeltaBeta = beta;
  MATRIX OldBeta = beta;
  MATRIX xOmegax(beta.Rows(), beta.Rows());

  // Matrices for REML estimation
  MATRIX *IOmega = NULL;
  MATRIX NewY, NewX, *NewVC = NULL;

  nMeanParam = x[0].Cols();

#ifdef REML
  if (Method != METHOD_ML && Method != METHOD_REML ) {
    printf("ERROR: Trying to maximize likelihood using non-(RE)ML methods.\nThis is not implemented yet.\n");
    exit(1);
  }
#else
  if (Method != METHOD_ML) {
    printf("ERROR: Trying to maximize likelihood using non-ML methods.\nThis is not implemented yet.\n");
    exit(1);
  }
#endif

  // Should combine all the matrices to one large set
  if (Method == METHOD_REML) {

#ifdef debug
  printf("Preparing REML code\n");
#endif

    NewVC = new MATRIX[nVC];

    nTotal=0;
    for (i = 0; i< nPedigrees; i++) {
      nTotal += nPedSize[i];
    }

    for (i = 0; i< nVC; i++) {
      NewVC[i].Resize(nTotal,nTotal);
      NewVC[i] = CombineMatrices(&VC[i*nPedigrees], nPedigrees);
    }

    NewY.Resize(nTotal,1);
    NewX.Resize(nTotal,beta.Rows());

    NewY = AppendMatrices(y, nPedigrees);
    NewX = AppendMatrices(x, nPedigrees);

    xT.Resize(beta.Rows(), nTotal);
    xT = Transpose(NewX);

    MATRIX tempmat(nTotal, nTotal);
    tempmat = 0;
    for (i = 1; i<=nTotal; i++) {
      tempmat(i,i) = 1;
    }
    NewY = (tempmat - NewX*Inverse(xT*NewX)*xT)*NewY;

    // Holds the pedigree variances below
    IOmega = new MATRIX[nPedigrees];
    for (nFam = 0; nFam < nPedigrees; nFam++){      
      IOmega[nFam].Resize(nPedSize[nFam], nPedSize[nFam]);
    }
    P.Resize(nTotal, nTotal);
    InvOmega.Resize(nTotal, nTotal);
    InvOmega2.Resize(nTotal, nTotal);


  }

#ifdef debug
  printf("Start maximizing\n");
#endif

  // Start checking that everything looks ok

  theta = start;  

  // Starts iterating
  for (nIter = 1; nIter <= MIXEDMODELMAXITER; nIter++)  {

#ifdef debug
    printf("Iteration: %3d\n", nIter);
#endif
    Fisher = 0.0;
    Deriv  = 0.0;
    DeltaBeta = 0.0;

    dLogLike = 0.0;
    LogLike = 0.0;

    xOmegax = 0.0;
    

    if (Method == METHOD_ML) {
      // This part calculates the first derivative and the Fisher information
      // matrix using Maximum likelihood methods

      // Go through each family/pedigree
      // The block diagonal structure is kept when doing ML
      for (nFam = 0; nFam < nPedigrees; nFam++)	{
	// Calculate the current variance for the family
	Omega.Resize(nPedSize[nFam], nPedSize[nFam]);
	InvOmega.Resize(nPedSize[nFam], nPedSize[nFam]);
	P.Resize(nPedSize[nFam], nPedSize[nFam]);
	xT.Resize(x[nFam].Cols(), x[nFam].Rows());
	xT = Transpose(x[nFam]);
	
	Omega = 0.0;

	nCode = nFam*nPedigrees;

	for (i = 0; i < nVC; i++)
	  Omega += VC[nFam+nPedigrees*i]*theta(i+1,1);

	// Inverting Omega_i
	InvOmega = Inverse(Omega, &dLogDet);

	// Add to the log likelihood
	MATRIX mu = (x[nFam] * beta);
	
	xOmegax   += (xT*InvOmega)*x[nFam];
	DeltaBeta += xT*(InvOmega*(y[nFam]-mu));
	LogLike += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));
     
	InvOmega2.Resize(nPedSize[nFam], nPedSize[nFam]);
	MATRIX UsedMatrix(nPedSize[nFam], nPedSize[nFam]);
	UsedMatrix = InvOmega;

	InvOmega2 = UsedMatrix*UsedMatrix;
	// Calculating the derivatives and Fisher scoring matrix
	for (i = 0; i < nVC; i++) {
	  Deriv(i+1,1) += - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
	    + ((Transpose(y[nFam] - mu)*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-mu)))(1,1);

	  for (j = i; j < nVC; j++) {
	    Fisher(i+1,j+1) += Trace(InvOmega2*VC[nFam + nPedigrees*i]*VC[nFam + nPedigrees*j]);
	    Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	  }
	}
      }
    }
    else if (Method == METHOD_REML) {
      // This part calculates the first derivative and the Fisher information
      // matrix using Restricted maximum likelihood methods

      // Calculate the current inverse variance matrix
      dLogDet = 0;
      for (nFam = 0; nFam < nPedigrees; nFam++)	{      
	IOmega[nFam] = 0;
	for (i = 0; i < nVC; i++) {	
	  IOmega[nFam] += VC[nFam+nPedigrees*i]*theta(i+1,1);
	}
	IOmega[nFam] = Inverse(IOmega[nFam], &dLogDet2);
	dLogDet += dLogDet2;
      }
      InvOmega = CombineMatrices(IOmega, nPedigrees);

      // Calculate P
      xOmegax = Inverse((xT*InvOmega)*NewX, &dLogDet2);
      P = InvOmega - ((InvOmega*NewX)*xOmegax)*(xT*InvOmega);   
      InvOmega2 = P*P;

      MATRIX mu = (NewX * beta);

      //	DeltaBeta += xT*(InvOmega*(y[nFam]-mu));
      //	LogLike += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));

      for (i = 0; i < nVC; i++) {
	Deriv(i+1,1) += - Trace(P*NewVC[i]) 
	  + ((Transpose(NewY)*P)*NewVC[i]*(P*(NewY)))(1,1);
	  
	for (j = i; j < nVC; j++) {
	  Fisher(i+1,j+1) += Trace(InvOmega2*NewVC[i]*NewVC[j]);
	  Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	}
      }

      LogLike = - 0.5*(nTotal*log(2*PI) + dLogDet + dLogDet2 + (Transpose(NewY)*P*(NewY))(1,1));
      LogLike -= - 0.5*nMeanParam*log(2*PI);
    }


    // Remember to scale the 1st and 2nd derivatives by 1/2
    Deriv  *= .5;
    Fisher *= .5;

    OldBeta = beta;

    if (Method == METHOD_ML) {
      beta += Inverse(xOmegax)*DeltaBeta;
      DeltaBeta = Inverse(xOmegax)*DeltaBeta;
    }

#ifdef constrained
    // Modifies the diagonal of the Fisher Matrix
    for (i=1; i<=nVC; i++) {
      // Dont do it for the covariances
      if (i != 5 && i != 8)
	Fisher(i,i) += constrainmu/theta(i,1);
      else 
	Fisher(i,i) += constrainmu/(theta(i,1));
    }

    //    Deriv.Print();

#endif

    // Inverts the matrix
    InvFisher = Inverse(Fisher);

    // Calculates the change in delta
    DeltaTheta = InvFisher*Deriv;
    WorkingTheta = DeltaTheta;

    // This change should be constrained
#ifdef constrained

    downscale = 1.0;  // Maximum steplength to stay positive

    // Walk along edge for parameters close to the edge
    for (i=1; i<=nVC; i++) {
      if (i != 5 && i != 8) {      
	if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	  WorkingTheta(i,1)=0;
	// Now finds the maximum allowed steplength 
	if (WorkingTheta(i,1)<0)
	  downscale = min(downscale, (1-SCALEEPSILON)*fabs(theta(i,1)/WorkingTheta(i,1)));
      }
      else {
	// Tage højde for at de nærmer sig randen afhængig af de tidligere parametre
	// For de ægte kovarianser
	
	if (i == 5)
	  edge = sqrt((theta(3,1)+WorkingTheta(3,1))*(theta(4,1)+WorkingTheta(4,1)));
	else if (i == 8)
	  edge = sqrt((theta(6,1)+WorkingTheta(6,1))*(theta(7,1)+WorkingTheta(7,1)));
	
	// Check at den ikke bliver for negativ
	if (theta(i,1)<-(edge-CLOSETOEDGE) && theta(i,1) + WorkingTheta(i,1) <= -edge) {
	  WorkingTheta(i,1)=0;
	}
	else if (theta(i,1)>(edge-CLOSETOEDGE) && theta(i,1) + WorkingTheta(i,1) >= edge) {
	  WorkingTheta(i,1)=0;
	}

	//	if (WorkingTheta(i,1)<0)
	//	  downscale = min(downscale, (1-SCALEEPSILON)*fabs(theta(i,1)/WorkingTheta(i,1)));	  


	//	if (theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
      }
    }
    // Scales the step down
    DeltaTheta = DeltaTheta*downscale;
    WorkingTheta = WorkingTheta*downscale;

#endif

    step = 2.0;
    do  {
      step *= .5;
      WorkingTheta = DeltaTheta*step;

      //      printf("HEJ: step: %f\n",step);

      for (i=1; i<=nVC; i++) {
	if (i != 5 && i != 8) {      
	  if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	    WorkingTheta(i,1)=0;
	}
      }

      newll = 0.0;
      // Calculate the new Omega

      if (Method == METHOD_ML) {
	for (nFam = 0; nFam < nPedigrees; nFam++)	{
	  // Calculate the current variance for the family
	  MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX xT = Transpose(x[nFam]);
	  MATRIX ugge = theta+WorkingTheta;

	  if (nVC>=5) {
	    ugge(5,1) = min(ugge(5,1), sqrt(ugge(3,1)*ugge(4,1)));
	    ugge(5,1) = max(ugge(5,1), -sqrt(ugge(3,1)*ugge(4,1)));
	  }
	  if (nVC>=8) {
	    ugge(8,1) = min(ugge(8,1), sqrt(ugge(6,1)*ugge(7,1)));
	    ugge(8,1) = max(ugge(8,1), -sqrt(ugge(6,1)*ugge(7,1)));
	  }

	  
	  Omega = 0.0;
	  
	  nCode = nFam*nPedigrees;
	  
	  for (i = 0; i < nVC; i++)
	    Omega += VC[nFam+nPedigrees*i]*(ugge(i+1,1));
	  
	  // Inverting Omega_i
	  InvOmega = Inverse(Omega, &dLogDet);
	  
	  // Add to the log likelihood
#ifdef fix_error
	  MATRIX mu = (x[nFam] * beta);
#else
	  MATRIX mu = (x[nFam] * OldBeta);
#endif
	  newll += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));
	}
	  
      } else if (Method == METHOD_REML) {
	// Calculate the current inverse variance matrix
	dLogDet = 0;
	for (nFam = 0; nFam < nPedigrees; nFam++)	{      
	  IOmega[nFam] = 0;
	  for (i = 0; i < nVC; i++) {	
	    //	    IOmega[nFam] += VC[nFam+nPedigrees*i]*(ugge(i+1,1));
	  }
	  IOmega[nFam] = Inverse(IOmega[nFam], &dLogDet2);
	  dLogDet += dLogDet2;
	}
	InvOmega = CombineMatrices(IOmega, nPedigrees);
	xOmegax = Inverse((xT*InvOmega)*NewX, &dLogDet2);
	P = InvOmega - ((InvOmega*NewX)*xOmegax)*(xT*InvOmega);   

	newll = - 0.5*(nTotal*log(2*PI) + dLogDet + dLogDet2 + (Transpose(NewY)*P*(NewY))(1,1));
    	newll -= -0.5*nMeanParam*log(2*PI);	

      }    
      //      printf("--> %f\n", newll);
    }
    while (newll < LogLike && fabs(newll-LogLike)>fabs(TOLERANCE*downscale*step) && log(step) > -100);
    //    while (newll < LogLike && fabs(newll-LogLike)>fabs(TOLERANCE*downscale*step) && log(step) > -30);

    // Modifies Theta
    theta += WorkingTheta;

    if (nVC>=5) {
      theta(5,1) = min(theta(5,1), sqrt(theta(3,1)*theta(4,1)));
      theta(5,1) = max(theta(5,1), -sqrt(theta(3,1)*theta(4,1)));
    }
    if (nVC>=8) {
      theta(8,1) = min(theta(8,1), sqrt(theta(6,1)*theta(7,1)));
      theta(8,1) = max(theta(8,1), -sqrt(theta(6,1)*theta(7,1)));
    }

    //    theta.Print();

    // XXX
    // Should also check the change in beta for ML 
    // Exiting due to small increas in LL or theta
    change = 0.0;

    for (i=1; i<= nVC; i++)
      change += fabs(WorkingTheta(i,1));

    // if ML and no change in mean parameters
    if (Method == METHOD_ML) {
      for (i=1; i<= beta.Rows(); i++)
	change += fabs(DeltaBeta(i,1));
    }


    // Should perhaps multiply by step also?
    // If yes, then same criteria should be used above in the while loop
    // XXX if (fabs(newll-LogLike)<fabs(TOLERANCE*downscale*step) || change<TOLERANCE*downscale*step ||  step <= 0.000001) {
    if (fabs(newll-LogLike)<fabs(TOLERANCE*downscale*step) || change<TOLERANCE*downscale*step ||  log(step) <= -30) {
      if (fabs(newll-LogLike)<fabs(TOLERANCE*downscale*step)) {
	convergence = 1;
      }
      else if (change<TOLERANCE*downscale*step) {
	convergence = 2;
      }
      else 
	convergence = 3;

      // Calculating the mean parameter estimates based on the variances
      if (Method == METHOD_REML) {
	NewY = AppendMatrices(y, nPedigrees);
	beta = xOmegax*(xT*InvOmega)*NewY;

	// Fixes xOmegax
	// Above, when using REML xOmegax is the INVERSE of xomega x, but for
	// ML is is xomegax
	// Therefore, change it back sp the right variance of the
	// mean is printed

	xOmegax = Inverse(xOmegax);
	
      }

      RES1 = theta(1,1);
      if (theta.Rows()>=2)       
	RES2 = theta(2,1);
      if (theta.Rows()>=3) 
	RES3 = theta(3,1);
      break;
    }
    fflush(stdout);
  }

#ifdef output
  if (PrintInfo != 0) {
    printf("Convergence reached after %d iterations (%s)\n",nIter, ConvergenceText[convergence]);
    printf("Estimating using %s.\n", Method==METHOD_ML ? "ML" : "REML");
    printf("Estimated variance components:\n");
    theta.Print();
    printf("Estimated mean parameters:\n");
    beta.Print();
    printf("Log likelihood at maximum: %6.3f\n", newll);
#ifdef showvariances	
    printf("Inverse Fisher information matrix for VC:\n");
    InvFisher.Print();
    printf("Variances of the mean parameters:\n");
    Inverse(xOmegax).Print();
#endif
  }
    
#endif

  fflush(stdout);

  if (Method == METHOD_REML) {
    delete[] IOmega;
    delete[] NewVC;
  }
  return (newll);
}



/*

  Maximizes a mixture of two gaussians using the EM-algorithm.
  Uses the same criteria as for MaximizeModel.

  Assumes identical mean model for both gaussians, but different 
  variances. This gives a length of start that is larger than above.

  p       holds initial guess on mixture probability
  include is a nVC*1 vector, with 1 for the VC to include in both gaussians
          One gaussian will hold all, while the other only will hold those with a 1

*/

// double Likelihood

// Calculates the likelihood of a single pedigree
double Likelihood(MATRIX y, MATRIX x, MATRIX InvOmega, double logdet, MATRIX betaparam, double *loglike) {
  double res;

  MATRIX resid(y.Rows(), 1);
  
  resid = y - x*betaparam;
  res = -0.5*(y.Rows() * log(2*PI) + logdet + ((Transpose(resid)*InvOmega)*resid)(1,1));

  *loglike = res;

  return(exp(res));
}


double MaximizeMixtureModel(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, double p, MATRIX include, int Method, int Constrain, int PrintInfo)
{
  int i, j, k, nIter, nFam, nCode;
  int nMeanParam, convergence;
  double dLogLike, dLogDet, dLogDet2, LogLike, downscale, change;
  double dLogDetR, famp, tempres, dtemp1, dtemp2;
  int nEMIter;
  int nObs;
  double step, newll, sum, tempfamp;
  double EMll, newEMll, setzero;
  double *lag;
  double obsdatallh;
  double TotalLogDet;
  double TotalLogDetR;
  int gridpoints = 10; // Number of starting p's to examine
  double bestp, bestll;
  int nMixtures = 2;  // Only uses a mixture of 2 Gaussians, only works for that!
  int ConvergenceCode;
  char buf[50];

  convergence = 0;
  lag = new double[MIXTURECONVLAG];

  MATRIX Fisher(nVC, nVC), InvFisher(nVC, nVC);
  MATRIX Deriv(nVC, 1);
  MATRIX theta(nVC, 1), DeltaTheta(nVC, 1), WorkingTheta(nVC, 1);;
  MATRIX DeltaBeta = beta;
  MATRIX OldBeta = beta;
  MATRIX minibeta(1,1);
  MATRIX xOmegax(beta.Rows(), beta.Rows());
  MATRIX MixMatrix(nPedigrees,nMixtures-1);  // This matrix holds the probability of ped(i,1) 
  // having the disease gene

  // Keeping a set of copies of variance matrices for the reduced model
  int nVCReduced = (int) include.Sum();
  if (nVCReduced<1 || nVCReduced>=nVC) {
    printf("ERROR: %d variance components in the reduced mixture model. Illegal number\n", nVCReduced);
    exit(1);
  }

  // Match the length of include and start
  if (start.Rows() != include.Rows()) {
    printf("ERROR: Inconsistency in the length of vectors start and include\n");
    exit(1);
    
  }


  // Calculate the total number of individuals
  nObs = 0;
  for (i = 0; i<nPedigrees; i++) {
    nObs += nPedSize[i];    
  }

  Deriv.Resize(nVC,1);
  theta.Resize(nVC,1);
  DeltaTheta.Resize(nVC,1);
  WorkingTheta.Resize(nVC,1);
  Fisher.Resize(nVC,nVC);
  InvFisher.Resize(nVC,nVC);

  nMeanParam = x[0].Cols();

  if (Method != METHOD_ML) {
    printf("ERROR: Trying to maximize likelihood using non-ML methods.\nThis is not implemented yet!\n");
    exit(1);
  }

#ifdef debug
  printf("Start maximizing mixture model\n");
#endif

  // Start checking that everything looks ok
  theta = start;

  // Initial bet on mixing proportion
  if (p>.93)
    p = .93;
  if (p<.07)
    p = .07;

  for (nEMIter = 1; nEMIter <= EMMAXITER; nEMIter++) {

    //    printf("Iteration : %d\n", nEMIter);

#ifdef debug
    printf("E-step\n");
    printf("Initial p : %f\n",p);
#endif

    //    printf("Initial p : %f\n",p);

    obsdatallh = 0;
    TotalLogDet = 0;
    TotalLogDetR = 0;
    // Preparing the mixture matrix
    for (i=0; i<nPedigrees; i++) {
      MATRIX Omega(nPedSize[i], nPedSize[i]), ROmega(nPedSize[i], nPedSize[i]);
      MATRIX InvOmega(nPedSize[i], nPedSize[i]);
      MATRIX InvROmega(nPedSize[i], nPedSize[i]);

      // Generates covariance matrices for both distributions
      Omega = 0;
      ROmega = 0;
      for (k = 0; k < nVC; k++) {
	setzero = theta(k+1,1);
#ifdef REALCALC
	if (setzero<CLOSETOEDGE) {
	  setzero = 0;
	}
#endif
	if (include(k+1,1)) {
	  ROmega += VC[i+nPedigrees*k]*setzero;
	}
	else {
	  Omega += VC[i+nPedigrees*k]*setzero;
	}	
      }

      // Inverting Omega_i
      InvOmega = Inverse(Omega, &dLogDet);
      InvROmega = Inverse(ROmega, &dLogDetR);

      // POSSIBLE NUMERICAL PROBLEM  HERE due to not using logs XXX
      for (j=1; j< nMixtures; j++) {
#ifdef diffbeta
	minibeta(1,1) = beta(1,1);
	Likelihood(y[i], x[i], InvOmega, dLogDet, minibeta, &dtemp1);
	minibeta(1,1) = beta(2,1);
	Likelihood(y[i], x[i], InvROmega, dLogDetR, minibeta, &dtemp2);
#else
	Likelihood(y[i], x[i], InvOmega, dLogDet, beta, &dtemp1);
	Likelihood(y[i], x[i], InvROmega, dLogDetR, beta, &dtemp2);
#endif

#ifdef improveconv
	if (p<PCLOSETOEDGE) {
	  MixMatrix(i+1,j) = PCLOSETOEDGE;
	}
	else if (p>1-PCLOSETOEDGE) {
	  MixMatrix(i+1,j) = 1-PCLOSETOEDGE;
	}
	else {	  
	  MixMatrix(i+1,j) = 1.0 / (1.0 + exp(log(1.0-p) - log(p) + dtemp2 - dtemp1));
	}
#else 
	// Use log-transform to help stabilize computations
	MixMatrix(i+1,j) = 1.0 / (1.0 + exp(log(1.0-p) - log(p) + dtemp2 - dtemp1));
#endif
	obsdatallh += log(p*exp(dtemp1) + (1-p)*exp(dtemp2));
	TotalLogDet  += dLogDet;
	TotalLogDetR += dLogDetR;
      }
    }    

#ifdef debug
    // Print obs. data llh before going through the E and M step
    printf("Observed data log-likelihood : %f\n", obsdatallh);
    printf("Log Determinants             : %f   and %f\n", TotalLogDet, TotalLogDetR);
    printf("M-step\n");
    printf("Pedigree-probabilities of containing QTL\n");
    //    MixMatrix.Print();
#endif

    // Set lag[0] if first iteration
    if (nEMIter == 1) {
      for (i = 0; i < MIXTURECONVLAG; i++) {
	lag[i] = obsdatallh+1+i;
      }
    }

    // Start by recalculating p as the average
    p = MixMatrix.Sum()/MixMatrix.Rows();
#ifdef debug
    printf("New p  :%5.3f\n",p);
#endif
    //    printf("New p  :%5.3f\n",p);
    //    MixMatrix.Print();

    // Starts iterating
#ifdef GEM
    for (nIter = 1; nIter <= (MAXSCORINGITER); nIter++) {
#else
    for (nIter = 1; nIter <= MIXEDMODELMAXITER; nIter++) {
#endif

#ifdef debug
      printf("Numerical maximization step %d\n", nIter);
#endif

      // Going through the E-step
 
      Fisher = 0.0;
      Deriv  = 0.0;
      DeltaBeta = 0.0;

      dLogLike = 0.0;
      LogLike = 0.0;

      xOmegax = 0.0;

      //      theta.Print();

	// Go through each family/pedigree
      for (nFam = 0; nFam < nPedigrees; nFam++)	{
	  // Calculate the current variance for the family
	  MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvROmega(nPedSize[nFam], nPedSize[nFam]), ROmega(nPedSize[nFam], nPedSize[nFam]); // Covariance matrices for the reduced model
	  MATRIX xT = Transpose(x[nFam]);
	  MATRIX P(nPedSize[nFam], nPedSize[nFam]);

	  Omega = 0.0;
	  ROmega = 0.0;

	  //	  printf("Family %d\n", nFam);

	  nCode = nFam*nPedigrees;

	  // Generates covariance matrices for both distributions
	  for (i = 0; i < nVC; i++) {
	    setzero = theta(i+1,1);
#ifdef REALCALC
	    if (setzero<CLOSETOEDGE) {
	      setzero = 0;
	    }
#endif
	    if (include(i+1,1)) {
	      ROmega += VC[nFam+nPedigrees*i]*setzero;
	    }
	    else {
	      Omega += VC[nFam+nPedigrees*i]*setzero;
	    }
	  }

	  // Inverting Omega_i
	  InvOmega  = Inverse(Omega, &dLogDet);
	  InvROmega = Inverse(ROmega, &dLogDetR);

	  // Probability of family having disease gene
	  famp = MixMatrix(nFam+1,1);

	  // Add to the log likelihood
#ifdef diffbeta
	  MATRIX mu = (x[nFam] * beta(1,1));
	  MATRIX mu2 = (x[nFam] * beta(2,1));

	  if (Method == METHOD_ML) {
	    xOmegax(1,1)   += ((xT*(InvOmega))*x[nFam])(1,1)*famp;
	    xOmegax(2,2)   += (xT*(InvROmega)*x[nFam])(1,1)*(1-famp);

	    DeltaBeta(1,1) += (xT*((InvOmega)*(y[nFam]-mu)))(1,1)*famp;
	    DeltaBeta(2,1) += (xT*((InvROmega)*(y[nFam]-mu2)))(1,1)*(1.0-famp);

	    minibeta(1,1) = beta(1,1);
	    Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, minibeta, &dtemp1);
	    minibeta(1,1) = beta(2,1);
	    Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, minibeta, &dtemp2);
	    LogLike   += famp*dtemp1 + (1-famp)*dtemp2;    
	  }
#else
	  MATRIX mu = (x[nFam] * beta);

	  if (Method == METHOD_ML) {
	    xOmegax   += (xT*(InvOmega*famp + InvROmega*(1-famp)))*x[nFam];
	    DeltaBeta += xT*((InvOmega*famp + InvROmega*(1-famp))*(y[nFam]-mu));
	    Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, beta, &dtemp1);
	    Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, beta, &dtemp2);
	    LogLike   += famp*dtemp1 + (1-famp)*dtemp2;
	  }
#endif
          
	  MATRIX InvOmega2(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX UsedMatrix(nPedSize[nFam], nPedSize[nFam]);

	  UsedMatrix = InvOmega;     
	  InvOmega2 = UsedMatrix*UsedMatrix;

	  // Calculating the derivatives and Fisher scoring matrix
	  // YYY Possible to optimize this depending on the derivatives
	  for (i = 0; i < nVC; i++) {
	      if (Method == METHOD_ML) {
		// First (full) model
		if (i == (nVC-nVCReduced)) {
		  UsedMatrix = InvROmega;
		  InvOmega2 = UsedMatrix*UsedMatrix;
		}

		if (i < (nVC-nVCReduced)) {
		  tempres = - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
		    + ((Transpose(y[nFam] - (x[nFam] * beta(1,1)))*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-(x[nFam] * beta(1,1)))))(1,1);
		  Deriv(i+1,1) += tempres*famp;
		}
		else {
		  tempres = - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
		    + ((Transpose(y[nFam] - (x[nFam] * beta(2,1)))*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-(x[nFam] * beta(2,1)))))(1,1);
		  Deriv(i+1,1) += tempres*(1.0-famp);
		}
	      } 

	      // Udregner Fisher scoring matricen
	      for (j = i; j < nVC; j++) {
		//		Fisher(i+1,j+1) = 0;
		tempfamp = famp;
		
		if (i >= (nVC-nVCReduced))
		  tempfamp = 1 - famp;

		if (!(i<(nVC-nVCReduced) && j>=(nVC-nVCReduced))) {
		  Fisher(i+1,j+1) += tempfamp*Trace(InvOmega2*VC[nFam + nPedigrees*i]*VC[nFam + nPedigrees*j]);
		}
		// Make it symmetric
		Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	      }
	  }
      }
      //      printf("\n");

      // Remember to scale the 1st and 2nd derivatives by 1/2
      Deriv  *= .5;
      Fisher *= .5;

      //      printf("OldLOGLIKE: %f  (based on the following theta:)\n", LogLike);
      //      theta.Print();
      //      printf("+++++++++++++++++++++++++++++\n");

      //    Deriv.Print();
      //    Fisher.Print();

      OldBeta = beta;

#ifdef constrained
      // Modifies the diagonal of the Fisher Matrix
      for (i=1; i<=nVC; i++)
	Fisher(i,i) += constrainmu/theta(i,1);

#endif

      // Inverts the matrix
      InvFisher = Inverse(Fisher);

      // Calculates the change in delta
      DeltaTheta = InvFisher*Deriv;
      WorkingTheta = DeltaTheta;

	// This change should be constrained
#ifdef constrained

      downscale = 1.0;  // Maximum steplength to stay positive

	// Walk along edge for parameters close to the edge
      for (i=1; i<=nVC; i++) {
	if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	  WorkingTheta(i,1)=0;
	// Now finds the maximum allowed steplength 
	if (WorkingTheta(i,1)<0)
	  downscale = min(downscale, (1-SCALEEPSILON)*fabs(theta(i,1)/WorkingTheta(i,1)));
      }
      // Scales the step down
      DeltaTheta = DeltaTheta*downscale;
      WorkingTheta = WorkingTheta*downscale;
#endif

      newEMll = LogLike;
      if (nEMIter == 1) {
	EMll = newEMll - 10;
      }

      if (Method == METHOD_ML) {
	beta += Inverse(xOmegax)*DeltaBeta;
	DeltaBeta = Inverse(xOmegax)*DeltaBeta;
      }

      //      printf("Starting step-halving\n");

      step = 2.0;
      do {
	  step *= .5;
	  WorkingTheta = DeltaTheta*step;

	  for (i=1; i<=nVC; i++) {
	    if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	      WorkingTheta(i,1)=0;
	  }

	  //	  printf("Step : %f\n", step);
	  //	  printf("ooooo\n");
	  //	  WorkingTheta.Print();
	  //	  printf("ooooo\n");


	  newll = 0.0;

	  // Calculate the new Omega
	  for (nFam = 0; nFam < nPedigrees; nFam++) {
	      // Calculate the current variance for the family
	      MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	      MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	      MATRIX InvROmega(nPedSize[nFam], nPedSize[nFam]), ROmega(nPedSize[nFam], nPedSize[nFam]); // Covariance matrices for the reduced model
	      MATRIX xT = Transpose(x[nFam]);

	      Omega = 0.0;

		// YYY
	      nCode = nFam*nPedigrees;

	      // Generates covariance matrices for both distributions
	      for (i = 0; i < nVC; i++) {
		setzero = WorkingTheta(i+1,1);
#ifdef REALCALC
		if (setzero<CLOSETOEDGE) {
		  setzero = 0;
		}
#endif
		if (include(i+1,1)) {
		  ROmega += VC[nFam+nPedigrees*i]*(theta(i+1,1) + setzero);
		}
		else {
		  Omega += VC[nFam+nPedigrees*i]*(theta(i+1,1) + setzero);
		}
	      }

	      // Inverting Omega_i
	      InvOmega  = Inverse(Omega, &dLogDet);
	      InvROmega = Inverse(ROmega, &dLogDetR);

		// YYY Her skal rettes. Hvad??
		// Husk at checke matricerne bruge ovenfor i sporet og ...
	      if (Method == METHOD_ML) {
#ifdef diffbeta
		minibeta(1,1) = (OldBeta + DeltaBeta*step)(1,1);
		Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, minibeta, &dtemp1);
		minibeta(1,1) = (OldBeta + DeltaBeta*step)(2,1);
		Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, minibeta, &dtemp2);
#ifdef newllh
		newll += MixMatrix(nFam+1,1)*dtemp1 + (1-MixMatrix(nFam+1,1))*dtemp2;
#else
		newll += log(MixMatrix(nFam+1,1)*exp(dtemp1) + (1-MixMatrix(nFam+1,1))*exp(dtemp2));
#endif
#else
		newll += log(MixMatrix(nFam+1,1)*Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, OldBeta + DeltaBeta*step, &dtemp1) + (1-MixMatrix(nFam+1,1))*Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, OldBeta + DeltaBeta*step, &dtemp2));
#endif
	      }
	  //	  printf("*** %f ***     => %f\n", newll, LogLike);
	  //	  (theta+WorkingTheta).Print();
	  //	  printf("Above is current theta\n");
	  }
      }
      while (newll < LogLike);
      //      while (newll < LogLike && fabs(newll-LogLike)>fabs(TOLERANCE*downscale*LogLike));

      //      printf("LOGLIKE at Fisher scoring maximum: %f\n", newll);
      
      // Modifies Theta 
      theta += WorkingTheta;
      beta = OldBeta + DeltaBeta*step;

      //      printf("--------\n");
      //      theta.Print();

      // XXX
      // Should also check the change in beta for ML 
      // Exiting due to small increas in LL or theta
      change = fabs(WorkingTheta.Sum());

	// if ML and no change in mean parameters
      if (Method == METHOD_ML) {
	  for (i=1; i<= beta.Rows(); i++)
	    change += fabs(DeltaBeta(i,1));
	}

      // Should perhaps multiply by step also?
      // If yes, then same criteria should be used above in the while loop
      // Has the Newton-Raphson iterations converged? 
      // I.e., has Q been maximized?
      if (fabs(newll-LogLike)<fabs(step*TOLERANCE*downscale*LogLike*10) || change<step*TOLERANCE*downscale) {
	//	convergence = 1;
	break;
      }
    }  // End of maximize (NR of internal part)

    //    printf("Comparison:   %5.3f   %f\n", p, obsdatallh);
    //    theta.Print();

    // Fix the lag
    for (i = MIXTURECONVLAG-1; i; i--) {
      lag[i] = lag[i-1];      
    }
    lag[0] = obsdatallh;
  
    // Check convergence between EM-steps
    //    printf("UUUU: %f\n", fabs(TOLERANCE*downscale*10*lag[0]));
    double ak = (lag[0]-lag[1])/(lag[1]-lag[2]);
    double ak2 = (lag[1]-lag[2])/(lag[2]-lag[3]);
    double lak = lag[1] + (lag[0]-lag[1])/(1-ak);
    double lak2 = lag[2] + (lag[1]-lag[2])/(1-ak2);
    if (fabs(lak-lak2)<TOLERANCE) {
      convergence = 1;
      break;
    }
  } // End of EM-algorithm

  // Check convergence criterion. 
  // Did we reach the maximum number of iteration and no convergence 
  if (convergence == 0 && nEMIter == EMMAXITER) {
    convergence = 4;
  }

  if (PrintInfo) {
#ifdef debug
    printf("Pedigree probabilities:\n");
    MixMatrix.Print();
#endif
    printf("Convergence reached after %d EM-iterations (%s)\n",nEMIter, ConvergenceText[convergence]);

    // Check for small relative determinants
    // The 2.3 is for log(.1)
    if (fabs(TotalLogDet - TotalLogDetR) > 2.3*nObs) {
      printf("WARNING: Possible singular maximum.\n");
    }

    printf("Estimated disease-pedigree frequency: %f\n", p);
    printf("Estimated variance components:\n");
    theta.Print();
    printf("Estimated mean parameters:\n");
    beta.Print();
    printf("Log likelihood at maximum: %6.3f\n", obsdatallh);
#ifdef showvariances	
    printf("Inverse Fisher information matrix for VC:\n");
    InvFisher.Print();
    printf("Variances of the mean parameters:\n");
    Inverse(xOmegax).Print();
#endif

    printf("PERLGRAB: %f  %f  %f\n", p, obsdatallh, theta(3,1));
  }

  delete[] lag;


  fflush(stdout);

  return (newll);
  }





double MaximizeModifiedMixtureModel(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, double p, MATRIX include, int Method, int Constrain, int PrintInfo)
{
  double penalty = 2;
  int i, j, k, nIter, nFam, nCode;
  int nMeanParam, convergence;
  double dLogLike, dLogDet, dLogDet2, LogLike, downscale, change;
  double dLogDetR, famp, tempres, dtemp1, dtemp2;
  int nEMIter;
  int nObs;
  double step, newll, sum, tempfamp;
  double setzero;
  double *lag;
  double obsdatallh;
  double TotalLogDet;
  double TotalLogDetR;
  int gridpoints = 10; // Number of starting p's to examine
  double bestp, bestll;
  int nMixtures = 2;  // Only uses a mixture of 2 Gaussians, only works for that!
  int ConvergenceCode;
  char buf[50];

  convergence = 0;
  lag = new double[MIXTURECONVLAG];

  MATRIX Fisher(nVC, nVC), InvFisher(nVC, nVC);
  MATRIX Deriv(nVC, 1);
  MATRIX theta(nVC, 1), DeltaTheta(nVC, 1), WorkingTheta(nVC, 1);;
  MATRIX DeltaBeta = beta;
  MATRIX OldBeta = beta;
  MATRIX minibeta(1,1);
  MATRIX xOmegax(beta.Rows(), beta.Rows());
  MATRIX MixMatrix(nPedigrees,nMixtures-1);  // This matrix holds the probability of ped(i,1) 
  // having the disease gene

  // Keeping a set of copies of variance matrices for the reduced model
  int nVCReduced = (int) include.Sum();
  if (nVCReduced<1 || nVCReduced>=nVC) {
    printf("ERROR: %d variance components in the reduced mixture model. Illegal number\n", nVCReduced);
    exit(1);
  }

  // Match the length of include and start
  if (start.Rows() != include.Rows()) {
    printf("ERROR: Inconsistency in the length of vectors start and include\n");
    exit(1);
    
  }


  // Calculate the total number of individuals
  nObs = 0;
  for (i = 0; i<nPedigrees; i++) {
    nObs += nPedSize[i];    
  }

  Deriv.Resize(nVC,1);
  theta.Resize(nVC,1);
  DeltaTheta.Resize(nVC,1);
  WorkingTheta.Resize(nVC,1);
  Fisher.Resize(nVC,nVC);
  InvFisher.Resize(nVC,nVC);

  nMeanParam = x[0].Cols();

  if (Method != METHOD_ML) {
    printf("ERROR: Trying to maximize likelihood using non-ML methods.\nThis is not implemented yet!\n");
    exit(1);
  }

#ifdef debug
  printf("Start maximizing mixture model\n");
#endif

  // Start checking that everything looks ok
  theta = start;

  // Initial bet on mixing proportion
  if (p>.93)
    p = .93;
  if (p<.07)
    p = .07;

  for (nEMIter = 1; nEMIter <= EMMAXITER; nEMIter++) {

#ifdef debug
    printf("Iteration : %d\n", nEMIter);
    printf("E-step\n");
    printf("  Initial p : %f\n",p);
#endif

    obsdatallh = 0;
    TotalLogDet = 0;
    TotalLogDetR = 0;
    // Preparing the mixture matrix
    for (i=0; i<nPedigrees; i++) {
      MATRIX Omega(nPedSize[i], nPedSize[i]), ROmega(nPedSize[i], nPedSize[i]);
      MATRIX InvOmega(nPedSize[i], nPedSize[i]);
      MATRIX InvROmega(nPedSize[i], nPedSize[i]);

      // Generates covariance matrices for both distributions
      Omega = 0;
      ROmega = 0;
      for (k = 0; k < nVC; k++) {
	setzero = theta(k+1,1);
#ifdef REALCALC
	if (setzero<CLOSETOEDGE) {
	  setzero = 0;
	}
#endif
	if (include(k+1,1)) {
	  ROmega += VC[i+nPedigrees*k]*setzero;
	}
	else {
	  Omega += VC[i+nPedigrees*k]*setzero;
	}	
      }

      // Inverting Omega_i
      InvOmega = Inverse(Omega, &dLogDet);
      InvROmega = Inverse(ROmega, &dLogDetR);

      // POSSIBLE NUMERICAL PROBLEM  HERE due to not using logs XXX
      for (j=1; j< nMixtures; j++) {
#ifdef diffbeta
	minibeta(1,1) = beta(1,1);
	Likelihood(y[i], x[i], InvOmega, dLogDet, minibeta, &dtemp1);
	minibeta(1,1) = beta(2,1);
	Likelihood(y[i], x[i], InvROmega, dLogDetR, minibeta, &dtemp2);
#else
	Likelihood(y[i], x[i], InvOmega, dLogDet, beta, &dtemp1);
	Likelihood(y[i], x[i], InvROmega, dLogDetR, beta, &dtemp2);
#endif

#ifdef improveconv
	if (p<PCLOSETOEDGE) {
	  MixMatrix(i+1,j) = PCLOSETOEDGE;
	}
	else if (p>1-PCLOSETOEDGE) {
	  MixMatrix(i+1,j) = 1-PCLOSETOEDGE;
	}
	else {	  
	  MixMatrix(i+1,j) = 1.0 / (1.0 + exp(log(1.0-p) - log(p) + dtemp2 - dtemp1));
	}
#else 
	// Use log-transform to help stabilize computations
	MixMatrix(i+1,j) = 1.0 / (1.0 + exp(log(1.0-p) - log(p) + dtemp2 - dtemp1));
#endif
	obsdatallh += log(p*exp(dtemp1) + (1-p)*exp(dtemp2));
	TotalLogDet  += dLogDet;
	TotalLogDetR += dLogDetR;
      }
    }
    // YYY   - add penalty term to log likelihood
    obsdatallh += penalty*log(4*p*(1-p));


#ifdef debug
    // Print obs. data llh before going through the E and M step
    printf("Observed data log-likelihood : %f\n", obsdatallh);
    printf("Log Determinants             : %f   and %f\n", TotalLogDet, TotalLogDetR);
    printf("M-step\n");
    printf("Pedigree-probabilities of containing QTL\n");
    //    MixMatrix.Print();
#endif

    // Set lag[0] if first iteration
    if (nEMIter == 1) {
      for (i = 0; i < MIXTURECONVLAG; i++) {
	lag[i] = obsdatallh+10+3*i;
      }
    }

    // Start by recalculating p as the average

    // YYY
    // Update the estimate for the mixture proportion
    p = (MixMatrix.Sum() + penalty)/(MixMatrix.Rows()+2*penalty);

#ifdef debug
    printf("New update for p :%5.3f\n",p);
#endif

    // Starts iterating
#ifdef GEM
    for (nIter = 1; nIter <= (MAXSCORINGITER); nIter++) {
#else
    for (nIter = 1; nIter <= MIXEDMODELMAXITER; nIter++) {
#endif

#ifdef debug
      printf("Numerical maximization step %d\n", nIter);
#endif

      // Going through the E-step
 
      Fisher = 0.0;
      Deriv  = 0.0;
      DeltaBeta = 0.0;

      dLogLike = 0.0;
      LogLike = 0.0;

      xOmegax = 0.0;

      //      theta.Print();

	// Go through each family/pedigree
      for (nFam = 0; nFam < nPedigrees; nFam++)	{
	  // Calculate the current variance for the family
	  MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvROmega(nPedSize[nFam], nPedSize[nFam]), ROmega(nPedSize[nFam], nPedSize[nFam]); // Covariance matrices for the reduced model
	  MATRIX xT = Transpose(x[nFam]);
	  MATRIX P(nPedSize[nFam], nPedSize[nFam]);

	  Omega = 0.0;
	  ROmega = 0.0;

	  //	  printf("Family %d\n", nFam);

	  nCode = nFam*nPedigrees;

	  // Generates covariance matrices for both distributions
	  for (i = 0; i < nVC; i++) {
	    setzero = theta(i+1,1);
#ifdef REALCALC
	    if (setzero<CLOSETOEDGE) {
	      setzero = 0;
	    }
#endif
	    if (include(i+1,1)) {
	      ROmega += VC[nFam+nPedigrees*i]*setzero;
	    }
	    else {
	      Omega += VC[nFam+nPedigrees*i]*setzero;
	    }
	  }

	  // Inverting Omega_i
	  InvOmega  = Inverse(Omega, &dLogDet);
	  InvROmega = Inverse(ROmega, &dLogDetR);

	  // Probability of family having disease gene
	  famp = MixMatrix(nFam+1,1);

	  // Add to the log likelihood
#ifdef diffbeta
	  MATRIX mu = (x[nFam] * beta(1,1));
	  MATRIX mu2 = (x[nFam] * beta(2,1));

	  if (Method == METHOD_ML) {
	    xOmegax(1,1)   += ((xT*(InvOmega))*x[nFam])(1,1)*famp;
	    xOmegax(2,2)   += (xT*(InvROmega)*x[nFam])(1,1)*(1-famp);

	    DeltaBeta(1,1) += (xT*((InvOmega)*(y[nFam]-mu)))(1,1)*famp;
	    DeltaBeta(2,1) += (xT*((InvROmega)*(y[nFam]-mu2)))(1,1)*(1.0-famp);

	    minibeta(1,1) = beta(1,1);
	    Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, minibeta, &dtemp1);
	    minibeta(1,1) = beta(2,1);
	    Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, minibeta, &dtemp2);
	    LogLike   += famp*dtemp1 + (1-famp)*dtemp2;    
	  }
#else
	  MATRIX mu = (x[nFam] * beta);

	  if (Method == METHOD_ML) {
	    xOmegax   += (xT*(InvOmega*famp + InvROmega*(1-famp)))*x[nFam];
	    DeltaBeta += xT*((InvOmega*famp + InvROmega*(1-famp))*(y[nFam]-mu));
	    Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, beta, &dtemp1);
	    Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, beta, &dtemp2);
	    LogLike   += famp*dtemp1 + (1-famp)*dtemp2;
	  }
#endif
          
	  MATRIX InvOmega2(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX UsedMatrix(nPedSize[nFam], nPedSize[nFam]);

	  UsedMatrix = InvOmega;     
	  InvOmega2 = UsedMatrix*UsedMatrix;

	  // Calculating the derivatives and Fisher scoring matrix
	  // YYY Possible to optimize this depending on the derivatives
	  for (i = 0; i < nVC; i++) {
	      if (Method == METHOD_ML) {
		// First (full) model
		if (i == (nVC-nVCReduced)) {
		  UsedMatrix = InvROmega;
		  InvOmega2 = UsedMatrix*UsedMatrix;
		}

		if (i < (nVC-nVCReduced)) {
		  tempres = - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
		    + ((Transpose(y[nFam] - (x[nFam] * beta(1,1)))*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-(x[nFam] * beta(1,1)))))(1,1);
		  Deriv(i+1,1) += tempres*famp;
		}
		else {
		  tempres = - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
		    + ((Transpose(y[nFam] - (x[nFam] * beta(2,1)))*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-(x[nFam] * beta(2,1)))))(1,1);
		  Deriv(i+1,1) += tempres*(1.0-famp);
		}
	      } 

	      // Udregner Fisher scoring matricen
	      for (j = i; j < nVC; j++) {
		//		Fisher(i+1,j+1) = 0;
		tempfamp = famp;
		
		if (i >= (nVC-nVCReduced))
		  tempfamp = 1 - famp;

		if (!(i<(nVC-nVCReduced) && j>=(nVC-nVCReduced))) {
		  Fisher(i+1,j+1) += tempfamp*Trace(InvOmega2*VC[nFam + nPedigrees*i]*VC[nFam + nPedigrees*j]);
		}
		// Make it symmetric
		Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	      }
	  }
      }
      //      printf("\n");

      // Remember to scale the 1st and 2nd derivatives by 1/2
      Deriv  *= .5;
      Fisher *= .5;

      //      printf("OldLOGLIKE: %f  (based on the following theta:)\n", LogLike);
      //      theta.Print();
      //      printf("+++++++++++++++++++++++++++++\n");

      //    Deriv.Print();
      //    Fisher.Print();

      OldBeta = beta;

#ifdef constrained
      // Modifies the diagonal of the Fisher Matrix
      for (i=1; i<=nVC; i++)
	Fisher(i,i) += constrainmu/theta(i,1);

#endif

      // Inverts the matrix
      InvFisher = Inverse(Fisher);

      // Calculates the change in delta
      DeltaTheta = InvFisher*Deriv;
      WorkingTheta = DeltaTheta;

	// This change should be constrained
#ifdef constrained

      downscale = 1.0;  // Maximum steplength to stay positive

	// Walk along edge for parameters close to the edge
      for (i=1; i<=nVC; i++) {
	if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	  WorkingTheta(i,1)=0;
	// Now finds the maximum allowed steplength 
	if (WorkingTheta(i,1)<0)
	  downscale = min(downscale, (1-SCALEEPSILON)*fabs(theta(i,1)/WorkingTheta(i,1)));
      }
      // Scales the step down
      DeltaTheta = DeltaTheta*downscale;
      WorkingTheta = WorkingTheta*downscale;
#endif

      // YYY - penalty
      //      LogLike += penalty*log(4*p*(1-p));

      if (Method == METHOD_ML) {
	beta += Inverse(xOmegax)*DeltaBeta;
	DeltaBeta = Inverse(xOmegax)*DeltaBeta;
      }

      //      printf("Starting step-halving\n");

      step = 2.0;
      do {
	  step *= .5;
	  WorkingTheta = DeltaTheta*step;

	  for (i=1; i<=nVC; i++) {
	    if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	      WorkingTheta(i,1)=0;
	  }

	  //	  printf("Step : %f\n", step);
	  //	  printf("ooooo\n");
	  //	  WorkingTheta.Print();
	  //	  printf("ooooo\n");


	  newll = 0.0;

	  // Calculate the new Omega
	  for (nFam = 0; nFam < nPedigrees; nFam++) {
	      // Calculate the current variance for the family
	      MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	      MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	      MATRIX InvROmega(nPedSize[nFam], nPedSize[nFam]), ROmega(nPedSize[nFam], nPedSize[nFam]); // Covariance matrices for the reduced model
	      MATRIX xT = Transpose(x[nFam]);

	      Omega = 0.0;

		// YYY
	      nCode = nFam*nPedigrees;

	      // Generates covariance matrices for both distributions
	      for (i = 0; i < nVC; i++) {
		setzero = WorkingTheta(i+1,1);
#ifdef REALCALC
		if (setzero<CLOSETOEDGE) {
		  setzero = 0;
		}
#endif
		if (include(i+1,1)) {
		  ROmega += VC[nFam+nPedigrees*i]*(theta(i+1,1) + setzero);
		}
		else {
		  Omega += VC[nFam+nPedigrees*i]*(theta(i+1,1) + setzero);
		}
	      }

	      // Inverting Omega_i
	      InvOmega  = Inverse(Omega, &dLogDet);
	      InvROmega = Inverse(ROmega, &dLogDetR);

		// YYY Her skal rettes. Hvad??
		// Husk at checke matricerne bruge ovenfor i sporet og ...
	      if (Method == METHOD_ML) {
#ifdef diffbeta
		minibeta(1,1) = (OldBeta + DeltaBeta*step)(1,1);
		Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, minibeta, &dtemp1);
		minibeta(1,1) = (OldBeta + DeltaBeta*step)(2,1);
		Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, minibeta, &dtemp2);
#ifdef newllh
		newll += MixMatrix(nFam+1,1)*dtemp1 + (1-MixMatrix(nFam+1,1))*dtemp2;
#else
		newll += log(MixMatrix(nFam+1,1)*exp(dtemp1) + (1-MixMatrix(nFam+1,1))*exp(dtemp2));
#endif
#else
		newll += log(MixMatrix(nFam+1,1)*Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, OldBeta + DeltaBeta*step, &dtemp1) + (1-MixMatrix(nFam+1,1))*Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, OldBeta + DeltaBeta*step, &dtemp2));
#endif
	      }
	      //	  	  printf("*** %f ***     => %f\n", newll, LogLike);
	  //	  (theta+WorkingTheta).Print();
	  //	  printf("Above is current theta\n");
	  }
	  //	  newll += penalty*log(p*(1-p));
      }
      while (newll < LogLike);
      //      while (newll < LogLike && fabs(newll-LogLike)>fabs(TOLERANCE*downscale*LogLike));

      //      printf("LOGLIKE at Fisher scoring maximum: %f\n", newll);
      
      // Modifies Theta 
      theta += WorkingTheta;
      beta = OldBeta + DeltaBeta*step;

      //      printf("--------\n");
      //      theta.Print();

      // XXX
      // Should also check the change in beta for ML 
      // Exiting due to small increas in LL or theta
      change = fabs(WorkingTheta.Sum());

	// if ML and no change in mean parameters
      if (Method == METHOD_ML) {
	  for (i=1; i<= beta.Rows(); i++)
	    change += fabs(DeltaBeta(i,1));
	}

      // Should perhaps multiply by step also?
      // If yes, then same criteria should be used above in the while loop
      // Has the Newton-Raphson iterations converged? 
      // I.e., has Q been maximized?
      if (fabs(newll-LogLike)<fabs(step*TOLERANCE*downscale*LogLike*10) || change<step*TOLERANCE*downscale) {
	//	convergence = 1;
	break;
      }
    }  // End of maximize (NR of internal part)

    //    printf("Comparison:   %5.3f   %f\n", p, obsdatallh);
    //    theta.Print();

    // Fix the lag
    for (i = MIXTURECONVLAG-1; i; i--) {
      lag[i] = lag[i-1];      
    }
    lag[0] = obsdatallh;
  
    // Check convergence between EM-steps
    //    printf("UUUU: %f\n", fabs(TOLERANCE*downscale*10*lag[0]));
    double ak = (lag[0]-lag[1])/(lag[1]-lag[2]);
    double ak2 = (lag[1]-lag[2])/(lag[2]-lag[3]);
    double lak = lag[1] + (lag[0]-lag[1])/(1-ak);
    double lak2 = lag[2] + (lag[1]-lag[2])/(1-ak2);
    if (1000*fabs(lak-lak2)<TOLERANCE) {
      convergence = 1;
      break;
    }
  } // End of EM-algorithm

  // Check convergence criterion. 
  // Did we reach the maximum number of iteration and no convergence 
  if (convergence == 0 && nEMIter == EMMAXITER) {
    convergence = 4;
  }

  if (PrintInfo) {
#ifdef debug
    printf("Pedigree probabilities:\n");
    MixMatrix.Print();
#endif
    printf("Convergence reached after %d EM-iterations (%s)\n",nEMIter, ConvergenceText[convergence]);

    // Check for small relative determinants
    // The 2.3 is for log(.1)
    if (fabs(TotalLogDet - TotalLogDetR) > 2.3*nObs) {
      printf("WARNING: Possible singular maximum.\n");
    }

    printf("Estimated disease-pedigree frequency: %f\n", p);
    printf("Estimated variance components:\n");
    theta.Print();
    printf("Estimated mean parameters:\n");
    beta.Print();
    printf("Log likelihood at maximum: %6.3f\n", obsdatallh);
#ifdef showvariances	
    printf("Inverse Fisher information matrix for VC:\n");
    InvFisher.Print();
    printf("Variances of the mean parameters:\n");
    Inverse(xOmegax).Print();
#endif

    printf("PERLGRAB: %f  %f  %f\n", p, obsdatallh, theta(3,1));
  }

  delete[] lag;


  fflush(stdout);

  return (newll);
  }





double MaximizeModifiedMixtureModelNULL(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, double p, MATRIX include, int Method, int Constrain, int PrintInfo)
{
  double penalty = 2;
  int i, j, k, nIter, nFam, nCode;
  int nMeanParam, convergence;
  double dLogLike, dLogDet, dLogDet2, LogLike, downscale, change;
  double dLogDetR, famp, tempres, dtemp1, dtemp2;
  int nEMIter;
  int nObs;
  double step, newll, sum, tempfamp;
  double setzero;
  double *lag;
  double obsdatallh;
  double TotalLogDet;
  double TotalLogDetR;
  int gridpoints = 10; // Number of starting p's to examine
  double bestp, bestll;
  int nMixtures = 2;  // Only uses a mixture of 2 Gaussians, only works for that!
  int ConvergenceCode;
  char buf[50];

  convergence = 0;
  lag = new double[MIXTURECONVLAG];

  MATRIX Fisher(nVC, nVC), InvFisher(nVC, nVC);
  MATRIX Deriv(nVC, 1);
  MATRIX theta(nVC, 1), DeltaTheta(nVC, 1), WorkingTheta(nVC, 1);;
  MATRIX DeltaBeta = beta;
  MATRIX OldBeta = beta;
  MATRIX minibeta(1,1);
  MATRIX xOmegax(beta.Rows(), beta.Rows());
  MATRIX MixMatrix(nPedigrees,nMixtures-1);  // This matrix holds the probability of ped(i,1) 
  // having the disease gene

  // Keeping a set of copies of variance matrices for the reduced model
  int nVCReduced = (int) include.Sum();
  if (nVCReduced<1 || nVCReduced>=nVC) {
    printf("ERROR: %d variance components in the reduced mixture model. Illegal number\n", nVCReduced);
    exit(1);
  }

  // Match the length of include and start
  if (start.Rows() != include.Rows()) {
    printf("ERROR: Inconsistency in the length of vectors start and include\n");
    exit(1);
    
  }

  // Calculate the total number of individuals
  nObs = 0;
  for (i = 0; i<nPedigrees; i++) {
    nObs += nPedSize[i];    
  }


  Deriv.Resize(nVC,1);
  theta.Resize(nVC,1);
  DeltaTheta.Resize(nVC,1);
  WorkingTheta.Resize(nVC,1);
  Fisher.Resize(nVC,nVC);
  InvFisher.Resize(nVC,nVC);

  nMeanParam = x[0].Cols();

  if (Method != METHOD_ML) {
    printf("ERROR: Trying to maximize likelihood using non-ML methods.\nThis is not implemented yet!\n");
    exit(1);
  }

#ifdef debug
  printf("Start maximizing mixture model\n");
#endif

  // Start checking that everything looks ok
  theta = start;

  // Initial bet on mixing proportion
  p = .5;

  for (nEMIter = 1; nEMIter <= EMMAXITER; nEMIter++) {

    //    printf("Iteration : %d\n", nEMIter);

#ifdef debug
    printf("E-step\n");
    printf("Initial p : %f\n",p);
#endif

    //    printf("Initial p : %f\n",p);

    obsdatallh = 0;
    TotalLogDet = 0;
    TotalLogDetR = 0;
    // Preparing the mixture matrix
    for (i=0; i<nPedigrees; i++) {
      MATRIX Omega(nPedSize[i], nPedSize[i]), ROmega(nPedSize[i], nPedSize[i]);
      MATRIX InvOmega(nPedSize[i], nPedSize[i]);
      MATRIX InvROmega(nPedSize[i], nPedSize[i]);

      // Generates covariance matrices for both distributions
      Omega = 0;
      ROmega = 0;
      for (k = 0; k < nVC; k++) {
	setzero = theta(k+1,1);
#ifdef REALCALC
	if (setzero<CLOSETOEDGE) {
	  setzero = 0;
	}
#endif
	if (include(k+1,1)==1) {
	  ROmega += VC[i+nPedigrees*k]*setzero;
	}
	else if (include(k+1,1)==0){
	  Omega += VC[i+nPedigrees*k]*setzero;
	}	
      }

      // Inverting Omega_i
      InvOmega = Inverse(Omega, &dLogDet);
      InvROmega = Inverse(ROmega, &dLogDetR);

      // POSSIBLE NUMERICAL PROBLEM  HERE due to not using logs XXX
      for (j=1; j< nMixtures; j++) {
#ifdef diffbeta
	minibeta(1,1) = beta(1,1);
	Likelihood(y[i], x[i], InvOmega, dLogDet, minibeta, &dtemp1);
	minibeta(1,1) = beta(2,1);
	Likelihood(y[i], x[i], InvROmega, dLogDetR, minibeta, &dtemp2);
#else
	Likelihood(y[i], x[i], InvOmega, dLogDet, beta, &dtemp1);
	Likelihood(y[i], x[i], InvROmega, dLogDetR, beta, &dtemp2);
#endif

#ifdef improveconv
	if (p<PCLOSETOEDGE) {
	  MixMatrix(i+1,j) = PCLOSETOEDGE;
	}
	else if (p>1-PCLOSETOEDGE) {
	  MixMatrix(i+1,j) = 1-PCLOSETOEDGE;
	}
	else {	  
	  MixMatrix(i+1,j) = 1.0 / (1.0 + exp(log(1.0-p) - log(p) + dtemp2 - dtemp1));
	}
#else 
	// Use log-transform to help stabilize computations
	MixMatrix(i+1,j) = 1.0 / (1.0 + exp(log(1.0-p) - log(p) + dtemp2 - dtemp1));
#endif
	obsdatallh += log(p*exp(dtemp1) + (1-p)*exp(dtemp2));
	TotalLogDet  += dLogDet;
	TotalLogDetR += dLogDetR;
      }
    }
    // YYY   - add penalty term to log likelihood
    obsdatallh += penalty*log(4*p*(1-p));


#ifdef debug
    // Print obs. data llh before going through the E and M step
    printf("Observed data log-likelihood : %f\n", obsdatallh);
    printf("Log Determinants             : %f   and %f\n", TotalLogDet, TotalLogDetR);
    printf("M-step\n");
    printf("Pedigree-probabilities of containing QTL\n");
    //    MixMatrix.Print();
#endif

    // Set lag[0] if first iteration
    if (nEMIter == 1) {
      for (i = 0; i < MIXTURECONVLAG; i++) {
	lag[i] = obsdatallh+10+3*i;
      }
    }

    // Start by recalculating p as the average

    // YYY

#ifdef debug
    printf("New p  :%5.3f\n",p);
#endif
    //    printf("New p  :%5.3f\n",p);
    //    MixMatrix.Print();

    // Starts iterating
#ifdef GEM
    for (nIter = 1; nIter <= (MAXSCORINGITER); nIter++) {
#else
    for (nIter = 1; nIter <= MIXEDMODELMAXITER; nIter++) {
#endif

#ifdef debug
      printf("Numerical maximization step %d\n", nIter);
#endif

      // Going through the E-step
 
      Fisher = 0.0;
      Deriv  = 0.0;
      DeltaBeta = 0.0;

      dLogLike = 0.0;
      LogLike = 0.0;

      xOmegax = 0.0;

      //      theta.Print();

	// Go through each family/pedigree
      for (nFam = 0; nFam < nPedigrees; nFam++)	{
	  // Calculate the current variance for the family
	  MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvROmega(nPedSize[nFam], nPedSize[nFam]), ROmega(nPedSize[nFam], nPedSize[nFam]); // Covariance matrices for the reduced model
	  MATRIX xT = Transpose(x[nFam]);
	  MATRIX P(nPedSize[nFam], nPedSize[nFam]);

	  Omega = 0.0;
	  ROmega = 0.0;

	  //	  printf("Family %d\n", nFam);

	  nCode = nFam*nPedigrees;

	  // Generates covariance matrices for both distributions
	  for (i = 0; i < nVC; i++) {
	    setzero = theta(i+1,1);
#ifdef REALCALC
	    if (setzero<CLOSETOEDGE) {
	      setzero = 0;
	    }
#endif
	    if (include(i+1,1) == 1) {
	      ROmega += VC[nFam+nPedigrees*i]*setzero;
	    }
	    else if (include(i+1,1) == 0){
	      Omega += VC[nFam+nPedigrees*i]*setzero;
	    }
	  }

	  // Inverting Omega_i
	  InvOmega  = Inverse(Omega, &dLogDet);
	  InvROmega = Inverse(ROmega, &dLogDetR);

	  // Probability of family having disease gene
	  famp = MixMatrix(nFam+1,1);

	  // Add to the log likelihood
#ifdef diffbeta
	  MATRIX mu = (x[nFam] * beta(1,1));
	  MATRIX mu2 = (x[nFam] * beta(2,1));

	  if (Method == METHOD_ML) {
	    xOmegax(1,1)   += ((xT*(InvOmega))*x[nFam])(1,1)*famp;
	    xOmegax(2,2)   += (xT*(InvROmega)*x[nFam])(1,1)*(1-famp);

	    DeltaBeta(1,1) += (xT*((InvOmega)*(y[nFam]-mu)))(1,1)*famp;
	    DeltaBeta(2,1) += (xT*((InvROmega)*(y[nFam]-mu2)))(1,1)*(1.0-famp);

	    minibeta(1,1) = beta(1,1);
	    Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, minibeta, &dtemp1);
	    minibeta(1,1) = beta(2,1);
	    Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, minibeta, &dtemp2);
	    LogLike   += famp*dtemp1 + (1-famp)*dtemp2;    
	  }
#else
	  MATRIX mu = (x[nFam] * beta);

	  if (Method == METHOD_ML) {
	    xOmegax   += (xT*(InvOmega*famp + InvROmega*(1-famp)))*x[nFam];
	    DeltaBeta += xT*((InvOmega*famp + InvROmega*(1-famp))*(y[nFam]-mu));
	    Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, beta, &dtemp1);
	    Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, beta, &dtemp2);
	    LogLike   += famp*dtemp1 + (1-famp)*dtemp2;
	  }
#endif
          
	  MATRIX InvOmega2(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX UsedMatrix(nPedSize[nFam], nPedSize[nFam]);

	  UsedMatrix = InvOmega;     
	  InvOmega2 = UsedMatrix*UsedMatrix;

	  // Calculating the derivatives and Fisher scoring matrix
	  // YYY Possible to optimize this depending on the derivatives
	  for (i = 0; i < nVC; i++) {
	      if (Method == METHOD_ML) {
		// First (full) model
		if (i == (nVC-nVCReduced)) {
		  UsedMatrix = InvROmega;
		  InvOmega2 = UsedMatrix*UsedMatrix;
		}

		if (i < (nVC-nVCReduced)) {   // Carriers
		  tempres = - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
		    + ((Transpose(y[nFam] - (x[nFam] * beta(1,1)))*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-(x[nFam] * beta(1,1)))))(1,1);
		  Deriv(i+1,1) += tempres*famp;
		}
		else { // Non-carriers
		  tempres = - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
		    + ((Transpose(y[nFam] - (x[nFam] * beta(2,1)))*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-(x[nFam] * beta(2,1)))))(1,1);
		  Deriv(i+1,1) += tempres*(1.0-famp);
		}
	      } 

	      // Udregner Fisher scoring matricen
	      for (j = i; j < nVC; j++) {
		//		Fisher(i+1,j+1) = 0;
		tempfamp = famp;
		
		if (i >= (nVC-nVCReduced))
		  tempfamp = 1 - famp;

		if (!(i<(nVC-nVCReduced) && j>=(nVC-nVCReduced))) {
		  Fisher(i+1,j+1) += tempfamp*Trace(InvOmega2*VC[nFam + nPedigrees*i]*VC[nFam + nPedigrees*j]);
		}
		// Make it symmetric
		Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	      }
	  }
      }
      //      printf("\n");

      // Remember to scale the 1st and 2nd derivatives by 1/2
      Deriv  *= .5;
      Fisher *= .5;

      //      printf("OldLOGLIKE: %f  (based on the following theta:)\n", LogLike);
      //      theta.Print();
      //      printf("+++++++++++++++++++++++++++++\n");

      //    Deriv.Print();
      //    Fisher.Print();

      OldBeta = beta;

#ifdef constrained
      // Modifies the diagonal of the Fisher Matrix
      for (i=1; i<=nVC; i++) {
	if (theta(i,1)>0) {
	  Fisher(i,i) += constrainmu/theta(i,1);
	}
      }

#endif

      // Inverts the matrix
      InvFisher = Inverse(Fisher);

      // Calculates the change in delta
      DeltaTheta = InvFisher*Deriv;
      WorkingTheta = DeltaTheta;

	// This change should be constrained
#ifdef constrained

      downscale = 1.0;  // Maximum steplength to stay positive

	// Walk along edge for parameters close to the edge
      for (i=1; i<=nVC; i++) {
	if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	  WorkingTheta(i,1)=0;
	// Now finds the maximum allowed steplength 
	if (WorkingTheta(i,1)<0)
	  downscale = min(downscale, (1-SCALEEPSILON)*fabs(theta(i,1)/WorkingTheta(i,1)));
      }
      // Scales the step down
      DeltaTheta = DeltaTheta*downscale;
      WorkingTheta = WorkingTheta*downscale;
#endif

      if (Method == METHOD_ML) {
	beta += Inverse(xOmegax)*DeltaBeta;
	DeltaBeta = Inverse(xOmegax)*DeltaBeta;
      }

      //      printf("Starting step-halving\n");

      step = 2.0;
      do {
	  step *= .5;
	  WorkingTheta = DeltaTheta*step;

	  for (i=1; i<=nVC; i++) {
	    if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	      WorkingTheta(i,1)=0;
	  }

	  //	  printf("Step : %f\n", step);
	  //	  printf("ooooo\n");
	  //	  WorkingTheta.Print();
	  //	  printf("ooooo\n");


	  newll = 0.0;

	  // Calculate the new Omega
	  for (nFam = 0; nFam < nPedigrees; nFam++) {
	      // Calculate the current variance for the family
	      MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	      MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	      MATRIX InvROmega(nPedSize[nFam], nPedSize[nFam]), ROmega(nPedSize[nFam], nPedSize[nFam]); // Covariance matrices for the reduced model
	      MATRIX xT = Transpose(x[nFam]);

	      Omega = 0.0;

		// YYY
	      nCode = nFam*nPedigrees;

	      // Generates covariance matrices for both distributions
	      for (i = 0; i < nVC; i++) {
		setzero = WorkingTheta(i+1,1);
#ifdef REALCALC
		if (setzero<CLOSETOEDGE) {
		  setzero = 0;
		}
#endif
		if (include(i+1,1)==1) {
		  ROmega += VC[nFam+nPedigrees*i]*(theta(i+1,1) + setzero);
		}
		else if (include(i+1,1)==0) {
		  Omega += VC[nFam+nPedigrees*i]*(theta(i+1,1) + setzero);
		}
	      }

	      // Inverting Omega_i
	      InvOmega  = Inverse(Omega, &dLogDet);
	      InvROmega = Inverse(ROmega, &dLogDetR);

		// YYY Her skal rettes. Hvad??
		// Husk at checke matricerne bruge ovenfor i sporet og ...
	      if (Method == METHOD_ML) {
#ifdef diffbeta
		minibeta(1,1) = (OldBeta + DeltaBeta*step)(1,1);
		Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, minibeta, &dtemp1);
		minibeta(1,1) = (OldBeta + DeltaBeta*step)(2,1);
		Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, minibeta, &dtemp2);
#ifdef newllh
		newll += MixMatrix(nFam+1,1)*dtemp1 + (1-MixMatrix(nFam+1,1))*dtemp2;
#else
		newll += log(MixMatrix(nFam+1,1)*exp(dtemp1) + (1-MixMatrix(nFam+1,1))*exp(dtemp2));
#endif
#else
		newll += log(MixMatrix(nFam+1,1)*Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, OldBeta + DeltaBeta*step, &dtemp1) + (1-MixMatrix(nFam+1,1))*Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, OldBeta + DeltaBeta*step, &dtemp2));
#endif
	      }
	  //	  printf("*** %f ***     => %f\n", newll, LogLike);
	  //	  (theta+WorkingTheta).Print();
	  //	  printf("Above is current theta\n");
	  }
      }
      while (newll < LogLike);
      //      while (newll < LogLike && fabs(newll-LogLike)>fabs(TOLERANCE*downscale*LogLike));

      //      printf("LOGLIKE at Fisher scoring maximum: %f\n", newll);
      
      // Modifies Theta 
      theta += WorkingTheta;
      beta = OldBeta + DeltaBeta*step;

      //      printf("--------\n");
      //      theta.Print();

      // XXX
      // Should also check the change in beta for ML 
      // Exiting due to small increas in LL or theta
      change = fabs(WorkingTheta.Sum());

	// if ML and no change in mean parameters
      if (Method == METHOD_ML) {
	  for (i=1; i<= beta.Rows(); i++)
	    change += fabs(DeltaBeta(i,1));
	}

      // Should perhaps multiply by step also?
      // If yes, then same criteria should be used above in the while loop
      // Has the Newton-Raphson iterations converged? 
      // I.e., has Q been maximized?
      if (fabs(newll-LogLike)<fabs(step*TOLERANCE*downscale*LogLike*10) || change<step*TOLERANCE*downscale) {
	//	convergence = 1;
	break;
      }
    }  // End of maximize (NR of internal part)

    //    printf("Comparison:   %5.3f   %f\n", p, obsdatallh);
    //    theta.Print();

    // Fix the lag
    for (i = MIXTURECONVLAG-1; i; i--) {
      lag[i] = lag[i-1];      
    }
    lag[0] = obsdatallh;
  
    // Check convergence between EM-steps
    //    printf("UUUU: %f\n", fabs(TOLERANCE*downscale*10*lag[0]));
    double ak = (lag[0]-lag[1])/(lag[1]-lag[2]);
    double ak2 = (lag[1]-lag[2])/(lag[2]-lag[3]);
    double lak = lag[1] + (lag[0]-lag[1])/(1-ak);
    double lak2 = lag[2] + (lag[1]-lag[2])/(1-ak2);
    if (1000*fabs(lak-lak2)<TOLERANCE) {
      convergence = 1;
      break;
    }
  } // End of EM-algorithm

  // Check convergence criterion. 
  // Did we reach the maximum number of iteration and no convergence 
  if (convergence == 0 && nEMIter == EMMAXITER) {
    convergence = 4;
  }

  if (PrintInfo) {
#ifdef debug
    printf("Pedigree probabilities:\n");
    MixMatrix.Print();
#endif
    printf("Convergence reached after %d EM-iterations (%s)\n",nEMIter, ConvergenceText[convergence]);

    // Check for small relative determinants
    // The 2.3 is for log(.1)
    if (fabs(TotalLogDet - TotalLogDetR) > 2.3*nObs) {
      printf("WARNING: Possible singular maximum.\n");
    }

    printf("Estimated disease-pedigree frequency: %f\n", p);
    printf("Estimated variance components:\n");
    theta.Print();
    printf("Estimated mean parameters:\n");
    beta.Print();
    printf("Log likelihood at maximum: %6.3f\n", obsdatallh);
#ifdef showvariances	
    printf("Inverse Fisher information matrix for VC:\n");
    InvFisher.Print();
    printf("Variances of the mean parameters:\n");
    Inverse(xOmegax).Print();
#endif

    printf("PERLGRAB: %f  %f  %f\n", p, obsdatallh, theta(3,1));
  }

  delete[] lag;


  fflush(stdout);

  return (newll);
  }















































/*

  Maximizes a mixture of two gaussians using the EM-algorithm.
  Uses the same criteria as for MaximizeModel.

  Model is restricted to have all but grand mean and qtl variance in common.

  Please note, the the input for VC are different than from the
  general heterogeneity model. The VC's houls be I, A, and PI and 
  no include matrix. 


  WORKZONE


*/

#ifdef newthings

double MaximizeRestrictedMixtureModel(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, double p, int Method, int Constrain, int PrintInfo)
{
  int i, j, k, ii, nIter, nFam, nCode;
  int nMeanParam, convergence;
  double dLogLike, dLogDet, dLogDet2, LogLike, downscale, change;
  double dLogDetR, famp, tempres, dtemp1, dtemp2;
  int nEMIter;
  double step, newll, sum, tempfamp;
  double EMll, newEMll, setzero;
  double *lag;
  double obsdatallh;
  int gridpoints = 10; // Number of starting p's to examine
  double bestp, bestll;
  int nMixtures = 2;  // Only uses a mixture of 2 Gaussians, only works for that!
  int ConvergenceCode;
  char buf[50];

  convergence = 0;
  lag = new double[MIXTURECONVLAG];

  MATRIX Fisher(nVC, nVC), InvFisher(nVC, nVC);
  MATRIX Deriv(nVC, 1);
  MATRIX theta(nVC, 1), DeltaTheta(nVC, 1), WorkingTheta(nVC, 1);;
  MATRIX DeltaBeta = beta;
  MATRIX OldBeta = beta;
  MATRIX minibeta(1,1);
  MATRIX xOmegax(beta.Rows(), beta.Rows());
  MATRIX MixMatrix(nPedigrees,nMixtures-1);  // This matrix holds the probability of ped(i,1) 
  // having the disease gene


  Deriv.Resize(nVC,1);
  theta.Resize(nVC,1);
  DeltaTheta.Resize(nVC,1);
  WorkingTheta.Resize(nVC,1);
  Fisher.Resize(nVC,nVC);
  InvFisher.Resize(nVC,nVC);

  nMeanParam = x[0].Cols();

  if (Method != METHOD_ML)
    {
      printf("ERROR: Trying to maximize likelihood using non-ML methods.\nThis is not implemented yet!\n");
      exit(1);
    }

#ifdef debug
  printf("Start maximizing mixture model\n");
#endif

  // Start checking that everything looks ok
  theta = start;

  // Initial bet on mixing proportion
  if (p>.93)
    p = .93;
  if (p<.07)
    p = .07;

  for (nEMIter = 1; nEMIter <= EMMAXITER; nEMIter++) {

    //    printf("Iteration : %d\n", nEMIter);

#ifdef debug
    printf("E-step\n");
    printf("Initial p : %f\n",p);
#endif

    //    printf("Initial p : %f\n",p);

    obsdatallh = 0;
    // Preparing the mixture matrix
    for (i=0; i<nPedigrees; i++) {
      MATRIX Omega(nPedSize[i], nPedSize[i]), ROmega(nPedSize[i], nPedSize[i]);
      MATRIX InvOmega(nPedSize[i], nPedSize[i]);
      MATRIX InvROmega(nPedSize[i], nPedSize[i]);

      // Generates covariance matrices for both distributions
      // Omega er de IKKE-qtl-relaterede personer
      // ROmega er med qtl-delen
      Omega = 0;
      ROmega = 0;
      for (k = 0; k < nVC; k++) {
	setzero = theta(k+1,1);
#ifdef REALCALC
	if (setzero<CLOSETOEDGE) {
	  setzero = 0;
	}
#endif

	if (k == nVC-1) {
	  ROmega = Omega;
	  ROmega += VC[i+nPedigrees*k]*setzero;
	}
	else {
	  Omega += VC[i+nPedigrees*k]*setzero;
	}
      }

      // Inverting Omega_i
      InvOmega = Inverse(Omega, &dLogDet);
      InvROmega = Inverse(ROmega, &dLogDetR);

      // POSSIBLE NUMERICAL PROBLEM  HERE due to not using logs XXX
      for (j=1; j< nMixtures; j++) {
#ifdef diffbeta
	minibeta(1,1) = beta(1,1);
	Likelihood(y[i], x[i], InvOmega, dLogDet, minibeta, &dtemp1);
	minibeta(1,1) = beta(2,1);
	Likelihood(y[i], x[i], InvROmega, dLogDetR, minibeta, &dtemp2);
#else
	Likelihood(y[i], x[i], InvOmega, dLogDet, beta, &dtemp1);
	Likelihood(y[i], x[i], InvROmega, dLogDetR, beta, &dtemp2);
#endif
	// Use log-transform to help stabilize computations
	MixMatrix(i+1,j) = 1.0 / (1.0 + exp(log(1.0-p) - log(p) + dtemp2 - dtemp1));
	obsdatallh += log(p*exp(dtemp1) + (1-p)*exp(dtemp2));
      }
    }    

#ifdef debug
    // Print obs. data llh before going through the E and M step
    printf("Observed data log-likelihood : %f\n", obsdatallh);
    printf("M-step\n");
    printf("Pedigree-probabilities of containing QTL\n");
    //    MixMatrix.Print();
#endif

    // Set lag[0] if first iteration
    if (nEMIter == 1) {
      for (i = 0; i < MIXTURECONVLAG; i++) {
	lag[i] = obsdatallh+1;
      }
    }

    // Start by recalculating p as the average
    p = MixMatrix.Sum()/MixMatrix.Rows();
#ifdef debug
    printf("New p  :%5.3f\n",p);
#endif
    //    MixMatrix.Print();

    // Starts iterating
#ifdef GEM
    for (nIter = 1; nIter <= (MAXSCORINGITER); nIter++) {
#else
    for (nIter = 1; nIter <= MIXEDMODELMAXITER; nIter++) {
#endif

#ifdef debug
      printf("Numerical maximization step %d\n", nIter);
#endif

      // Going through the E-step
 
      Fisher = 0.0;
      Deriv  = 0.0;
      DeltaBeta = 0.0;

      dLogLike = 0.0;
      LogLike = 0.0;

      xOmegax = 0.0;

      //      theta.Print();

	// Go through each family/pedigree
      for (nFam = 0; nFam < nPedigrees; nFam++)	{
	  // Calculate the current variance for the family
	  MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvROmega(nPedSize[nFam], nPedSize[nFam]), ROmega(nPedSize[nFam], nPedSize[nFam]); // Covariance matrices for the reduced model
	  MATRIX xT = Transpose(x[nFam]);
	  MATRIX P(nPedSize[nFam], nPedSize[nFam]);

	  Omega = 0.0;
	  ROmega = 0.0;

	  //	  printf("Family %d\n", nFam);

	  nCode = nFam*nPedigrees;

	  // Generates covariance matrices for both distributions
	  for (i = 0; i < nVC; i++) {
	    setzero = theta(i+1,1);
#ifdef REALCALC
	    if (setzero<CLOSETOEDGE) {
	      setzero = 0;
	    }
#endif

	    if (i == nVC-1) {
	      ROmega = Omega;
	      ROmega += VC[nFam+nPedigrees*i]*setzero;
	    }
	    else {
	      Omega += VC[nFam+nPedigrees*i]*setzero;
	    }
	  }

	  // Inverting Omega_i
	  InvOmega  = Inverse(Omega, &dLogDet);
	  InvROmega = Inverse(ROmega, &dLogDetR);

	  // Probability of family having disease gene
	  famp = MixMatrix(nFam+1,1);

	  // Add to the log likelihood
#ifdef diffbeta
	  MATRIX mu = (x[nFam] * beta(1,1));
	  MATRIX mu2 = (x[nFam] * beta(2,1));

	  if (Method == METHOD_ML) {
	    xOmegax(1,1)   += ((xT*(InvOmega))*x[nFam])(1,1)*famp;
	    xOmegax(2,2)   += (xT*(InvROmega)*x[nFam])(1,1)*(1-famp);

	    DeltaBeta(1,1) += (xT*((InvOmega)*(y[nFam]-mu)))(1,1)*famp;
	    DeltaBeta(2,1) += (xT*((InvROmega)*(y[nFam]-mu2)))(1,1)*(1.0-famp);

	    minibeta(1,1) = beta(1,1);
	    Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, minibeta, &dtemp1);
	    minibeta(1,1) = beta(2,1);
	    Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, minibeta, &dtemp2);
	    LogLike   += famp*dtemp1 + (1-famp)*dtemp2;    
	  }
#else
	  MATRIX mu = (x[nFam] * beta);

	  if (Method == METHOD_ML) {
	    xOmegax   += (xT*(InvOmega*famp + InvROmega*(1-famp)))*x[nFam];
	    DeltaBeta += xT*((InvOmega*famp + InvROmega*(1-famp))*(y[nFam]-mu));
	    Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, beta, &dtemp1);
	    Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, beta, &dtemp2);
	    LogLike   += famp*dtemp1 + (1-famp)*dtemp2;
	  }
#endif
          
	  MATRIX InvOmega2(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX UsedMatrix(nPedSize[nFam], nPedSize[nFam]);

	  UsedMatrix = InvOmega;     
	  InvOmega2 = UsedMatrix*UsedMatrix;

	  // Calculating the derivatives and Fisher scoring matrix
	  // YYY Possible to optimize this depending on the derivatives
	  for (i = 0; i < nVC; i++) {
	      if (Method == METHOD_ML) {
		// First (full) model
		if (i == (nVC-1)) {
		  UsedMatrix = InvROmega;
		  InvOmega2 = UsedMatrix*UsedMatrix;
		}

		if (i < (nVC-nVCReduced)) {
		  tempres = - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
		    + ((Transpose(y[nFam] - (x[nFam] * beta(1,1)))*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-(x[nFam] * beta(1,1)))))(1,1);
		  Deriv(i+1,1) += tempres*famp;
		}
		else {
		  tempres = - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
		    + ((Transpose(y[nFam] - (x[nFam] * beta(2,1)))*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-(x[nFam] * beta(2,1)))))(1,1);
		  Deriv(i+1,1) += tempres*(1.0-famp);
		}
	      } 

	      // Udregner Fisher scoring matricen
	      for (j = i; j < nVC; j++) {
		//		Fisher(i+1,j+1) = 0;
		tempfamp = famp;
		
		if (i >= (nVC-nVCReduced))
		  tempfamp = 1 - famp;

		if (!(i<(nVC-nVCReduced) && j>=(nVC-nVCReduced))) {
		  Fisher(i+1,j+1) += tempfamp*Trace(InvOmega2*VC[nFam + nPedigrees*i]*VC[nFam + nPedigrees*j]);
		}
		// Make it symmetric
		Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	      }
	  }
      }
      //      printf("\n");

      // Remember to scale the 1st and 2nd derivatives by 1/2
      Deriv  *= .5;
      Fisher *= .5;

      //      printf("OldLOGLIKE: %f  (based on the following theta:)\n", LogLike);
      //      theta.Print();
      //      printf("+++++++++++++++++++++++++++++\n");

      //    Deriv.Print();
      //    Fisher.Print();

      OldBeta = beta;

#ifdef constrained
      // Modifies the diagonal of the Fisher Matrix
      for (i=1; i<=nVC; i++)
	Fisher(i,i) += constrainmu/theta(i,1);

#endif

      // Inverts the matrix
      InvFisher = Inverse(Fisher);

      // Calculates the change in delta
      DeltaTheta = InvFisher*Deriv;
      WorkingTheta = DeltaTheta;

	// This change should be constrained
#ifdef constrained

      downscale = 1.0;  // Maximum steplength to stay positive

	// Walk along edge for parameters close to the edge
      for (i=1; i<=nVC; i++)
	{
	  if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	    WorkingTheta(i,1)=0;
	  // Now finds the maximum allowed steplength 
	  if (WorkingTheta(i,1)<0)
	    downscale = min(downscale, (1-SCALEEPSILON)*fabs(theta(i,1)/WorkingTheta(i,1)));
	}
      // Scales the step down
      DeltaTheta = DeltaTheta*downscale;
      WorkingTheta = WorkingTheta*downscale;

#endif

      newEMll = LogLike;
      if (nEMIter == 1) {
	EMll = newEMll - 10;
      }

      if (Method == METHOD_ML) {
	beta += Inverse(xOmegax)*DeltaBeta;
	DeltaBeta = Inverse(xOmegax)*DeltaBeta;
      }

      //      printf("Starting step-halving\n");

      step = 2.0;
      do {
	  step *= .5;
	  WorkingTheta = DeltaTheta*step;

	  for (i=1; i<=nVC; i++) {
	    if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	      WorkingTheta(i,1)=0;
	  }

	  //	  printf("Step : %f\n", step);
	  //	  printf("ooooo\n");
	  //	  WorkingTheta.Print();
	  //	  printf("ooooo\n");


	  newll = 0.0;

	  // Calculate the new Omega
	  for (nFam = 0; nFam < nPedigrees; nFam++) {
	      // Calculate the current variance for the family
	      MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	      MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	      MATRIX InvROmega(nPedSize[nFam], nPedSize[nFam]), ROmega(nPedSize[nFam], nPedSize[nFam]); // Covariance matrices for the reduced model
	      MATRIX xT = Transpose(x[nFam]);

	      Omega = 0.0;

		// YYY
	      nCode = nFam*nPedigrees;

	      // Generates covariance matrices for both distributions
	      for (i = 0; i < nVC; i++) {
		setzero = WorkingTheta(i+1,1);
#ifdef REALCALC
		if (setzero<CLOSETOEDGE) {
		  setzero = 0;
		}
#endif
		if (include(i+1,1)) {
		  ROmega += VC[nFam+nPedigrees*i]*(theta(i+1,1) + setzero);
		}
		else {
		  Omega += VC[nFam+nPedigrees*i]*(theta(i+1,1) + setzero);
		}
	      }

	      // Inverting Omega_i
	      InvOmega  = Inverse(Omega, &dLogDet);
	      InvROmega = Inverse(ROmega, &dLogDetR);

		// YYY Her skal rettes. Hvad??
		// Husk at checke matricerne bruge ovenfor i sporet og ...
	      if (Method == METHOD_ML) {
#ifdef diffbeta
		minibeta(1,1) = (OldBeta + DeltaBeta*step)(1,1);
		Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, minibeta, &dtemp1);
		minibeta(1,1) = (OldBeta + DeltaBeta*step)(2,1);
		Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, minibeta, &dtemp2);
#ifdef newllh
		newll += MixMatrix(nFam+1,1)*dtemp1 + (1-MixMatrix(nFam+1,1))*dtemp2;
#else
		newll += log(MixMatrix(nFam+1,1)*exp(dtemp1) + (1-MixMatrix(nFam+1,1))*exp(dtemp2));
#endif
#else
		newll += log(MixMatrix(nFam+1,1)*Likelihood(y[nFam], x[nFam], InvOmega, dLogDet, OldBeta + DeltaBeta*step, &dtemp1) + (1-MixMatrix(nFam+1,1))*Likelihood(y[nFam], x[nFam], InvROmega, dLogDetR, OldBeta + DeltaBeta*step, &dtemp2));
#endif
	      }
	  //	  printf("*** %f ***     => %f\n", newll, LogLike);
	  //	  (theta+WorkingTheta).Print();
	  //	  printf("Above is current theta\n");
	  }
      }
      while (newll < LogLike);
      //      while (newll < LogLike && fabs(newll-LogLike)>fabs(TOLERANCE*downscale*LogLike));

      //      printf("LOGLIKE at Fisher scoring maximum: %f\n", newll);
      
      // Modifies Theta 
      theta += WorkingTheta;
      beta = OldBeta + DeltaBeta*step;

      //      theta.Print();

      // XXX
      // Should also check the change in beta for ML 
      // Exiting due to small increas in LL or theta
      change = fabs(WorkingTheta.Sum());

	// if ML and no change in mean parameters
      if (Method == METHOD_ML) {
	  for (i=1; i<= beta.Rows(); i++)
	    change += fabs(DeltaBeta(i,1));
	}

      // Should perhaps multiply by step also?
      // If yes, then same criteria should be used above in the while loop
      // Has the Newton-Raphson iterations converged? 
      // I.e. has Q been maximized?
      if (fabs(newll-LogLike)<fabs(TOLERANCE*downscale*LogLike*10) || change<TOLERANCE*downscale) {
	//	convergence = 1;
	break;
      }
    }  // End of maximize (NR of internal part)


    //    printf("Comparison: %f   %f\n", newll, LogLike);

    // Fix the lag
    for (i = MIXTURECONVLAG-1; i; i--) {
      lag[i] = lag[i-1];      
    }
    lag[0] = obsdatallh;
  
    // Check convergence between EM-steps
    //    printf("UUUU: %f\n", fabs(TOLERANCE*downscale*10*lag[0]));
    if (fabs(lag[MIXTURECONVLAG-1]-lag[0])<fabs(TOLERANCE*downscale*10) && nEMIter>MIXTURECONVLAG) {
      convergence = 1;
      break;
    }
  } // End of EM-algorithm

  // Check convergence criterion. 
  // Did we reach the maximum number of iteration and no convergence 
  if (convergence == 0 && nEMIter == EMMAXITER) {
    convergence = 3;
  }

  if (PrintInfo) {
#ifdef debug
    printf("Pedigree probabilities:\n");
    MixMatrix.Print();
#endif
    printf("Convergence reached after %d EM-iterations (%s)\n",nEMIter, ConvergenceText[convergence]);

    printf("Estimated disease-pedigree frequency: %f\n", p);
    printf("Estimated variance components:\n");
    theta.Print();
    printf("Estimated mean parameters:\n");
    beta.Print();
    printf("Log likelihood at maximum: %6.3f\n", obsdatallh);
#ifdef showvariances	
    printf("Inverse Fisher information matrix for VC:\n");
    InvFisher.Print();
    printf("Variances of the mean parameters:\n");
    Inverse(xOmegax).Print();
#endif

    printf("PERLGRAB: %f  %f  %f\n", p, obsdatallh, theta(3,1));
  }

  delete[] lag;


  fflush(stdout);

  return (newll);
  }
#endif




#ifdef GXE
  /****************************************
   *
   *
   *  Funktion til at udregne C matricen
   *
   *
   ****************************************/
MATRIX CalcC(MATRIX X, double mean, double parameter) 
{
  MATRIX res(X.Rows(),X.Rows());
  double calc1, calc2;
  int i, j;

  for (i=1; i<= X.Rows(); i++) {
    for (j=i; j<= X.Rows(); j++) {
      calc1 = exp((X(i,1)-mean)*parameter);
      calc1 = calc1/(calc1+1);
      calc2 = exp((X(j,1)-mean)*parameter);
      calc2 = calc2/(calc2+1);
      res(i,j) = calc1*calc2;
      res(j,i) = res(i,j);
    }
  }
  return(res);
}


MATRIX CalcCPrime(MATRIX X, double mean, double parameter) 
{
  MATRIX res(X.Rows(),X.Rows());
  double calc1, calc2;
  int i, j;

  for (i=1; i<= X.Rows(); i++) {
    for (j=i; j<= X.Rows(); j++) {
      calc2 = (exp(parameter*(X(i,1)-3*mean+2*X(j,1)))*X(i,1)+exp(parameter*(X(i,1)-2*mean+X(j,1)))*X(i,1) - exp(parameter*(2*X(i,1) - 3*mean + X(j,1)))*mean 
              - exp(parameter* (X(i,1) - 3* mean + 2* X(j,1)))* mean - 2* exp(parameter* (X(i,1) - 2* mean + X(j,1)))* mean
	       + exp(parameter* (2* X(i,1) - 3* mean + X(j,1)))* X(j,1) + exp(parameter* (X(i,1) - 2* mean + X(j,1)))* X(j,1));
      calc1 = (calc2)/(SQR(exp((X(i,1)-mean)*parameter)+1) * SQR(exp((X(j,1)-mean)*parameter)+1));
      res(i,j) = calc1;
      res(j,i) = res(i,j);
    }
  }
  return(res);
}



  /**************************************************
   *
   *  Funktion til at maximere GxE vekselvirkninger.
   *
   *  Tillader pt. kun en kovariat for env-effekten
   *
   *  Lader parametren i GXE-fordelingen være den sidste af 
   *  varianskomponent-parametrene
   *
   ***************************************************/

double MaximizeMixedGXEModel(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], int nGXEparam, MATRIX XX[], MATRIX start, MATRIX beta, int Method, int Constrain, int PrintInfo, double corfct)
{
  int i, j, nIter, nFam, nCode, nTotal, nInd;
  int nMeanParam, convergence;

  double dLogLike, dLogDet, dLogDet2, LogLike, downscale, change;
  double step, newll;
  double xi;
  double thingmean;

  convergence = 0;

  // Generel matrices
  MATRIX xT;
  MATRIX Omega, InvOmega, P, InvOmega2, F, Fprime, GxE;
  MATRIX Fisher(nVC+nGXEparam, nVC+nGXEparam), InvFisher(nVC+nGXEparam, nVC+nGXEparam);
  MATRIX Deriv(nVC+nGXEparam, 1);
  MATRIX theta(nVC+nGXEparam, 1), DeltaTheta(nVC+nGXEparam, 1), WorkingTheta(nVC+nGXEparam, 1);;  
  MATRIX DeltaBeta = beta;
  MATRIX OldBeta = beta;
  MATRIX xOmegax(beta.Rows(), beta.Rows());

  // Matrices for REML estimation
  MATRIX *IOmega;
  MATRIX *NewVC, NewY, NewX;


  nMeanParam = x[0].Cols();

#ifdef REML
  if (Method != METHOD_ML && Method != METHOD_REML ) {
    printf("ERROR: Trying to maximize likelihood using non-(RE)ML methods.\nThis is not implemented yet.\n");
    exit(1);
  }
#else
  if (Method != METHOD_ML) {
    printf("ERROR: Trying to maximize likelihood using non-ML methods.\nThis is not implemented yet.\n");
    exit(1);
  }
#endif

  // Should combine all the matrices to one large set
  if (Method == METHOD_REML) {

#ifdef debug
  printf("Preparing REML code\n");
#endif


    nTotal = 0;
    NewVC = new MATRIX[nVC];

    nTotal=0;
    for (i = 0; i< nPedigrees; i++) {
      nTotal += nPedSize[i];
    }

    for (i = 0; i< nVC; i++) {
      NewVC[i].Resize(nTotal,nTotal);
      NewVC[i] = CombineMatrices(&VC[i*nPedigrees], nPedigrees);
    }

    NewY.Resize(nTotal,1);
    NewX.Resize(nTotal,beta.Rows());

    NewY = AppendMatrices(y, nPedigrees);
    NewX = AppendMatrices(x, nPedigrees);

    xT.Resize(beta.Rows(), nTotal);
    xT = Transpose(NewX);

    MATRIX tempmat(nTotal, nTotal);
    tempmat = 0;
    for (i = 1; i<=nTotal; i++) {
      tempmat(i,i) = 1;
    }
    NewY = (tempmat - NewX*Inverse(xT*NewX)*xT)*NewY;

    // Holds the pedigree variances below
    IOmega = new MATRIX[nPedigrees];
    for (nFam = 0; nFam < nPedigrees; nFam++){      
      IOmega[nFam].Resize(nPedSize[nFam], nPedSize[nFam]);
    }
    P.Resize(nTotal, nTotal);
    InvOmega.Resize(nTotal, nTotal);
    InvOmega2.Resize(nTotal, nTotal);


  }

#ifdef debug
  printf("Start maximizing\n");
#endif

  // Start checking that everythink looks ok

  for (i = 1; i<=start.Rows(); i++)
    theta(i,1) = start(i,1);  
  theta(i,1) = 0;
  xi = 0;

  thingmean = 0;
  nInd = 0;
  for (nFam = 0; nFam < nPedigrees; nFam++)	{
    thingmean +=XX[nFam].Sum();
    nInd += nPedSize[nFam];    
  }
  thingmean /= (double)nInd;

  // Starts iterating
  for (nIter = 1; nIter <= MIXEDMODELMAXITER; nIter++)  {

#ifdef debug
    printf("Iteration: %3d\n", nIter);
#endif
    Fisher = 0.0;
    Deriv  = 0.0;
    DeltaBeta = 0.0;

    dLogLike = 0.0;
    LogLike = 0.0;

    xOmegax = 0.0;    
    

    if (Method == METHOD_ML) {
      // This part calculates the first derivative and the Fisher information
      // matrix using Maximum likelihood methods

      // Go through each family/pedigree
      // The block diagonal structure is kept when doing ML
      for (nFam = 0; nFam < nPedigrees; nFam++)	{
	// Calculate the current variance for the family
	Omega.Resize(nPedSize[nFam], nPedSize[nFam]);
	InvOmega.Resize(nPedSize[nFam], nPedSize[nFam]);
	F.Resize(nPedSize[nFam], nPedSize[nFam]);
	Fprime.Resize(nPedSize[nFam], nPedSize[nFam]);
	GxE.Resize(nPedSize[nFam], nPedSize[nFam]);
	P.Resize(nPedSize[nFam], nPedSize[nFam]);
	xT.Resize(x[nFam].Cols(), x[nFam].Rows());
	xT = Transpose(x[nFam]);
	
	Omega = 0.0;

	nCode = nFam*nPedigrees;

	// CALC F
	// CAlc Fprime
	F = CalcC(XX[nFam], thingmean, theta(nVC+1,1));
	//	theta.Print();
	//	F.Print();
	Fprime = CalcCPrime(XX[nFam], thingmean, theta(nVC+1,1));
	//	Fprime.Print();
	GxE = Hademard(F, VC[nFam+nPedigrees*(nVC-1)]);

	// Det er QTL'en der ganges på
	for (i = 0; i < nVC-1; i++) {
	  Omega += VC[nFam+nPedigrees*i]*theta(i+1,1);
	}
	Omega += GxE*theta(nVC,1);

	// Inverting Omega_i
	InvOmega = Inverse(Omega, &dLogDet);

	// Add to the log likelihood
	MATRIX mu = (x[nFam] * beta);
	
	xOmegax   += (xT*InvOmega)*x[nFam];
	DeltaBeta += xT*(InvOmega*(y[nFam]-mu));
	LogLike += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));
     
	InvOmega2.Resize(nPedSize[nFam], nPedSize[nFam]);
	MATRIX UsedMatrix(nPedSize[nFam], nPedSize[nFam]);
	UsedMatrix = InvOmega;

	InvOmega2 = UsedMatrix*UsedMatrix;
	// Calculating the derivatives and Fisher scoring matrix
	for (i = 0; i <= nVC; i++) {
	  // sigma_q
	  if (i==nVC) {
	    Deriv(i+1,1) += - Trace(Fprime*UsedMatrix) 
	      + ((Transpose(y[nFam] - mu)*UsedMatrix)*Fprime*(UsedMatrix*(y[nFam]-mu)))(1,1);
	  }
	  else {
	    Deriv(i+1,1) += - Trace(VC[nFam + nPedigrees*i]*UsedMatrix) 
	      + ((Transpose(y[nFam] - mu)*UsedMatrix)*VC[nFam + nPedigrees*i]*(UsedMatrix*(y[nFam]-mu)))(1,1);
	  }

	  for (j = i; j <= nVC; j++) {
	    if (j==nVC) {
	      if (i==j) {
		Fisher(i+1,j+1) += Trace(InvOmega2*Fprime*Fprime);
	      }
	      else 
		Fisher(i+1,j+1) += Trace(InvOmega2*VC[nFam + nPedigrees*i]*Fprime);
	      Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	    }
	    else {
	      Fisher(i+1,j+1) += Trace(InvOmega2*VC[nFam + nPedigrees*i]*VC[nFam + nPedigrees*j]);
	      Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	    }
	  }
	}
      }
    }
    else if (Method == METHOD_REML) {
      // This part calculates the first derivative and the Fisher information
      // matrix using Restricted maximum likelihood methods

      // Calculate the current inverse variance matrix
      dLogDet = 0;
      for (nFam = 0; nFam < nPedigrees; nFam++)	{      
	IOmega[nFam] = 0;
	for (i = 0; i < nVC; i++) {	
	  IOmega[nFam] += VC[nFam+nPedigrees*i]*theta(i+1,1);
	}
	IOmega[nFam] = Inverse(IOmega[nFam], &dLogDet2);
	dLogDet += dLogDet2;
      }
      InvOmega = CombineMatrices(IOmega, nPedigrees);

      // Calculate P
      xOmegax = Inverse((xT*InvOmega)*NewX, &dLogDet2);
      P = InvOmega - ((InvOmega*NewX)*xOmegax)*(xT*InvOmega);   
      InvOmega2 = P*P;

      MATRIX mu = (NewX * beta);

      //	DeltaBeta += xT*(InvOmega*(y[nFam]-mu));
      //	LogLike += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));

      for (i = 0; i < nVC; i++) {
	Deriv(i+1,1) += - Trace(P*NewVC[i]) 
	  + ((Transpose(NewY)*P)*NewVC[i]*(P*(NewY)))(1,1);
	  
	for (j = i; j < nVC; j++) {
	  Fisher(i+1,j+1) += Trace(InvOmega2*NewVC[i]*NewVC[j]);
	  Fisher(j+1,i+1)  = Fisher(i+1,j+1);
	}
      }

      LogLike = - 0.5*(nTotal*log(2*PI) + dLogDet + dLogDet2 + (Transpose(NewY)*P*(NewY))(1,1));
      LogLike -= - 0.5*nMeanParam*log(2*PI);
    }
    // End of REML

    // Remember to scale the 1st and 2nd derivatives by 1/2
    Deriv  *= .5;
    Fisher *= .5;

    OldBeta = beta;

    if (Method == METHOD_ML) {
      beta += Inverse(xOmegax)*DeltaBeta;
      DeltaBeta = Inverse(xOmegax)*DeltaBeta;
    }

#ifdef constrained
    // Modifies the diagonal of the Fisher Matrix
    for (i=1; i<=nVC; i++)
      Fisher(i,i) += constrainmu/theta(i,1);

#endif


    // Inverts the matrix
    InvFisher = Inverse(Fisher);

    // Calculates the change in delta
    DeltaTheta = InvFisher*Deriv;
    WorkingTheta = DeltaTheta;

    // This change should be constrained
#ifdef constrained

    downscale = 1.0;  // Maximum steplength to stay positive

    // Walk along edge for parameters close to the edge
    for (i=1; i<=nVC; i++) {
      if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	WorkingTheta(i,1)=0;
      // Now finds the maximum allowed steplength 
      if (WorkingTheta(i,1)<0)
	downscale = min(downscale, (1-SCALEEPSILON)*fabs(theta(i,1)/WorkingTheta(i,1)));
    }
    // Scales the step down
    DeltaTheta = DeltaTheta*downscale;
    WorkingTheta = WorkingTheta*downscale;
#endif

    step = 2.0;
    do  {
      step *= .5;
      WorkingTheta = DeltaTheta*step;

      for (i=1; i<=nVC; i++) {
	if (theta(i,1)>0 && theta(i,1)<CLOSETOEDGE && theta(i,1)<=-WorkingTheta(i,1))
	  WorkingTheta(i,1)=0;
      }

      newll = 0.0;
      // Calculate the new Omega
      if (Method == METHOD_ML) {
	for (nFam = 0; nFam < nPedigrees; nFam++)	{
	  // Calculate the current variance for the family
	  MATRIX Omega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX InvOmega(nPedSize[nFam], nPedSize[nFam]);
	  MATRIX xT = Transpose(x[nFam]);

	  F.Resize(nPedSize[nFam], nPedSize[nFam]);
	  Fprime.Resize(nPedSize[nFam], nPedSize[nFam]);
	  GxE.Resize(nPedSize[nFam], nPedSize[nFam]);
	  
	  Omega = 0.0;
	  
	  nCode = nFam*nPedigrees;

	  // CALC F
	  // CAlc Fprime
	  F = CalcC(XX[nFam], thingmean, theta(nVC+1,1) + WorkingTheta(nVC+1,1));
	  Fprime = CalcCPrime(XX[nFam], thingmean, theta(nVC+1,1) + WorkingTheta(nVC+1,1));
	  GxE = Hademard(F, VC[nFam+nPedigrees*(nVC-1)]);

	  // Det er QTL'en der ganges på
	  for (i = 0; i < nVC-1; i++) {
	    Omega += VC[nFam+nPedigrees*i]*(theta(i+1,1)+WorkingTheta(i+1,1));
	  }
	  Omega += GxE*(theta(nVC,1) + WorkingTheta(nVC,1));

	  // Inverting Omega_i
	  InvOmega = Inverse(Omega, &dLogDet);
	  
	  // Add to the log likelihood
#ifdef fix_error
	  MATRIX mu = (x[nFam] * beta);
#else
	  MATRIX mu = (x[nFam] * OldBeta);
#endif
	  newll += - 0.5*(nPedSize[nFam]*log(2*PI) + dLogDet + (Transpose(y[nFam]-mu)*InvOmega*(y[nFam]-mu))(1,1));
	}
	  
      } else if (Method == METHOD_REML) {
	// Calculate the current inverse variance matrix
	dLogDet = 0;
	for (nFam = 0; nFam < nPedigrees; nFam++)	{      
	  IOmega[nFam] = 0;
	  for (i = 0; i < nVC; i++) {	
	    IOmega[nFam] += VC[nFam+nPedigrees*i]*(theta(i+1,1) + WorkingTheta(i+1,1));
	  }
	  IOmega[nFam] = Inverse(IOmega[nFam], &dLogDet2);
	  dLogDet += dLogDet2;
	}
	InvOmega = CombineMatrices(IOmega, nPedigrees);
	xOmegax = Inverse((xT*InvOmega)*NewX, &dLogDet2);
	P = InvOmega - ((InvOmega*NewX)*xOmegax)*(xT*InvOmega);   

	newll = - 0.5*(nTotal*log(2*PI) + dLogDet + dLogDet2 + (Transpose(NewY)*P*(NewY))(1,1));
    	newll -= -0.5*nMeanParam*log(2*PI);	



      }    
    }
    while (newll < LogLike && fabs(newll-LogLike)>TOLERANCE*downscale*step && step > 0.000001);

    // Modifies Theta
    theta += WorkingTheta;

    // XXX
    // Should also check the change in beta for ML 
    // Exiting due to small increas in LL or theta
    change = 0.0;

    for (i=1; i<= nVC; i++)
      change += fabs(WorkingTheta(i,1));

    // if ML and no change in mean parameters
    if (Method == METHOD_ML) {
      for (i=1; i<= beta.Rows(); i++)
	change += fabs(DeltaBeta(i,1));
    }


    // Should perhaps multiply by step also?
    // If yes, then same criteria should be used above in the while loop
    if (fabs(newll-LogLike)<TOLERANCE*downscale*step || change<TOLERANCE*downscale*step || step < 0.000001) {
      if (fabs(newll-LogLike)<TOLERANCE*downscale*step) {
	convergence = 1;
      }
      else {
	convergence = 2;
      }

      // Calculating the mean parameter estimates based on the variances
      if (Method == METHOD_REML) {
	NewY = AppendMatrices(y, nPedigrees);
	beta = xOmegax*(xT*InvOmega)*NewY;

	// Fixes xOmegax
	// Above, when using REML xOmegax is the INVERSE of xomega x, but for
	// ML is is xomegax
	// Therefore, change it back sp the right variance of the
	// mean is printed

	xOmegax = Inverse(xOmegax);
	
      }

      RES1 = theta(1,1);
      if (theta.Rows()>=2)       
	RES2 = theta(2,1);
      if (theta.Rows()>=3) 
	RES3 = theta(3,1);
      break;
    }
  }
#ifdef output
  if (PrintInfo != 0) {
    printf("Convergence reached after %d iterations (%s)\n",nIter, ConvergenceText[convergence]);
    printf("Estimating using %s.\n", Method==METHOD_ML ? "ML" : "REML");
    printf("Estimated variance components:\n");
    theta.Print();
    printf("Estimated mean parameters:\n");
    beta.Print();
    printf("Log likelihood at maximum: %6.3f\n", newll);
#ifdef showvariances	
    printf("Inverse Fisher information matrix for VC:\n");
    InvFisher.Print();
    printf("Variances of the mean parameters:\n");
    Inverse(xOmegax).Print();
#endif
  }
    
#endif

  fflush(stdout);

  if (Method == METHOD_REML) {
    delete[] IOmega;
    delete[] NewVC;
  }


  return (newll);
}


#endif




#undef debug
#undef output
