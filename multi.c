/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "nrutil.h"
#include "pedipet.h"
#include "ibd.c"

#undef debug


#ifndef MULTIPOINTVC

#define MULTIPOINTVC
#define SHOWINFO

// #define debug
/*
#define AMOSMAXITER        40
#define MINIMUMMEASURE     0.00001  // Minimum values for constrained REML for measurement error
#define constrained
#define VCTOL             0.0001
#define CLOSETOEDGE       0.0001
#define SCALEEPSILON      0.1

#define stephalving
#define newmatrices
#define constrainmu       1
#define scalemult         // If the length of the new direction was scaled down, then allow for a smaller increase in 
                          // likelihood before breaking
 
*/

//
//
//

double VarianceOfMarker(OLDMATRIX *ibdmatrix, int relation, OLDMATRIX *relationshipmtx)
{
  double mean, var;
  int i, j, n;

  mean = 0.0;
  n = 0;
  for (i=1; i<=ibdmatrix->rows; i++)
  {
    for (j=i; j<=ibdmatrix->cols; j++)
    {
      if (relationshipmtx->element[i][j]==relation)
      {
        mean += ibdmatrix->element[i][j];
        n++;
      }
    }
  }
  if (n>0)
    mean /= n ;

  var = 0.0;
  for (i=1; i<=ibdmatrix->rows; i++)
  {
    for (j=i; j<=ibdmatrix->cols; j++)
    {
      if (relationshipmtx->element[i][j]==relation)
      {
        var += SQR(ibdmatrix->element[i][j]-mean);
      }
    }
  }

  if (n==0)
    return n;
  else
    return(var/n);  
}


//
// For now, change all info to families instead of one large matrix
//

//
// YYY Could save a LOT OF TIME by doing it by pedigree instead of total
//

// XXX Some mess here as kinship and delta are based on individuals  - should modify them
// ot only use indlist or to reduce them appropriately

void PrepareMultipoint(individual *indlist) {
  OLDMATRIX *dominance, *kinship, *markerIBD, *relationship, *tmpmat, *mat, *mat2;
  OLDMATRIX *betas, *VC, *C, *BaseC;
  OLDMATRIX *PIHAT;
  OLDMATRIX *InterMarkerVariance, *MarkerVariance;
  int i, j, k, relation;
  int NumberOfRelations, RelationshipPresent[REL_MAXREL];
  double location, var;

  // Calculates and prepares the required matrices
  #ifdef SHOWINFO 
    printf("Creating kinship, delta and relationships matrix\n"); 
  #endif
  kinship   = MakeKinshipMatrix();
  dominance = MakeDelta7Matrix(kinship);
  relationship = CalcDataRelationship(individuals);
  markerIBD = 0;

  // Go through all markers and calculate the IBD matrix for each
  for (i = 1; i<= numberofmarkers(); i++)
  {
    #ifdef SHOWINFO 
      printf("Calculating IBD-matrix for marker %d\n",i); 
    #endif
    // Use MC estimation
    tmpmat = CalcMCPIHat(indlist, i);
    addlist(&markerIBD, tmpmat);
  }

  // Figure out which relationships, 
  // that are present in the current dataset
  NumberOfRelations = 0;
  memset(RelationshipPresent, 0, sizeof(RelationshipPresent));
  for (i = 1; i< relationship->rows; i++)
  {
    for (j = i+1; j<= relationship->cols; j++)
    {
      RelationshipPresent[(int)relationship->element[i][j]] = 1;
    }
  }
  // Counting the numbers, that are present. 
  // This EXCLUDES the REL_NONE, REL_ID and REL_PARENT relationships
  // as there is no info from these
  for (i=0; i<REL_MAXREL; i++)
    if (RelationshipPresent[i] && (i != REL_NONE && i != REL_ID && i != REL_PARENT))
      NumberOfRelations++;

  // Preparing the covariance matrices between markers for different relationtypes
  InterMarkerVariance = 0;
  MarkerVariance = 0;
  
  for (k = 1; k<= NumberOfRelations; k++)
  {

    // UDREGNE HVILKEN RELÆATION VI ER I GANG MED!
    relation = REL_FULL_SIBS;
    #ifdef SHOWINFO 
      printf("Going through relationship %s  [%d]\n",RelationshipNames[relation], k); 
    #endif

    VC = MtxNew(numberofmarkers(), numberofmarkers());
    BaseC  = MtxNew(numberofmarkers(), 1); 

    for (i=1, mat = markerIBD; i<= numberofmarkers(); i++, mat = mat->next)
      BaseC->element[i][1] = VarianceOfMarker(mat, relation, relationship);
  

    for (i=1, mat = markerIBD; i<= numberofmarkers(); i++, mat = mat->next)
    {
      var = VarianceOfMarker(mat, relation, relationship);
      printf("The variance: %f\n", var);
      for (j=i, mat2 = mat; j<= numberofmarkers(); j++, mat2 = mat2->next)
      { 
        if (i==j)
          VC->element[i][j] = BaseC->element[i][1];
        else
        {
          VC->element[i][j] = CorrCoeffIBD(InterMarkerRecomb(i,j), relation)*
	                      BaseC->element[i][1]*BaseC->element[j][1]/varpi(relation);
          VC->element[j][i] = VC->element[i][j];
        }
	//  C->element[i][1] = CorrCoeffIBD(InverseMapFunction(fabs(MarkerDistance(location, i))/100), relation)*var;
      }
    }
    MtxInver(VC);

    addlist(&InterMarkerVariance, VC);
    addlist(&MarkerVariance, BaseC);
  }
  // Have now calculated the stuff needed for the regressions
  // Formula 16 in Almasy and Blangero
  // Need to multiply the location<->marker recomb frequency on the BaseC before using


  location = 40;
  //for (location)
  {

    // Going through each relationship
    for (k = 1, VC = InterMarkerVariance, BaseC = MarkerVariance; k<= NumberOfRelations; k++, VC = VC->next, BaseC = BaseC->next)
    {
      // First calculate the betas
      betas = MtxNew(numberofmarkers(), NumberOfRelations);
      C     = MtxNew(numberofmarkers(), NumberOfRelations);
      // PIHAT holds the new Pi-matrix at this specific location
      PIHAT = MtxNew(listlen(indlist), listlen(indlist));

      // Should this be baseC and not just C?

      // Calculate the distance to each marker for this relationship

      // Multiply this result to baseC

      // Calculate the betas


      MtxMulti(VC, BaseC, betas);
      for (i=1; i<=listlen(indlist); i++)
      {
	for (j=i; j<=listlen(indlist); j++)
	{
	  if (i==j)
	  {
	    PIHAT->element[i][j] = 1.0;	      
	  }
	  else
	  {
	    // HMmm. This i s wrong
	    if (relationship->element[i][j]==k)
	    {
   	      PIHAT->element[i][j] = 0.0;
	      //              for (ii=1; ii <= numberofmarkers(); ii++)
		//   	        PIHAT->element[i][j] += betas->element[ii][k]*();
	    }	    
	  }
	}
      }


      MtxDel(betas);
      MtxDel(C);
      MtxDel(PIHAT);
    }
  }


  /*
  // betas will hold the regression coefficients for each marker
  // and for each relation type

  vector = MtxNew(numberofmarkers(),1);
  VC = MtxNew(numberofmarkers(), numberofmarkers());
  C  = MtxNew(numberofmarkers(), 1); 
  // Holds all comp except multiplying with the marker-location recomb
  BaseC  = MtxNew(numberofmarkers(), 1); 


  location = 40;

  // Should only go through the relationships that are in the dataset
  for (relation = REL_FULL_SIBS; relation<REL_HALF_SIBS; relation++)  
  {
    // YYY
    // Could save time here by calculating all the variance for each 
    // relationtype first

    for (i=1, mat = markerIBD; i<= numberofmarkers(); i++, mat = mat->next)
    {
      var = VarianceOfMarker(mat, relation, relationship);
      printf("The variance: %f\n", var);
      for (j=i, mat2 = mat; j<= numberofmarkers(); j++, mat2 = mat2->next)
      { 
        if (i==j)
          VC->element[i][j] = var;
        else
        {
          VC->element[i][j] = CorrCoeffIBD(InterMarkerRecomb(i,j), relation)*var*VarianceOfMarker(mat2, relation, relationship)/varpi(relation);
          VC->element[j][i] = VC->element[i][j];
        }
        BaseC->element[i][1] = var;
        C->element[i][1] = CorrCoeffIBD(InverseMapFunction(fabs(MarkerDistance(location, i))/100), relation)*var;
      }
    }
    MtxInver(VC);
    MtxMulti(VC, C, vector);

    // Moving the result over to the betas
    for (i=1; i<=numberofmarkers(); i++)
      betas->element[i][relation] = vector->element[i][1];

  }

  MtxPrint(betas);

  MtxDel(VC);
  MtxDel(C);
  MtxDel(betas);
  MtxDel(vector);
  */

  MtxFreeList(MarkerVariance);
  MtxFreeList(InterMarkerVariance);

}


#endif


#undef debug
