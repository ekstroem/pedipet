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

#ifndef IBD_HEADER_FILE
#define IBD_HEADER_FILE

#include "structure.h"
#include "mathfct.h"

// ---------------------
//  Relationsships
// ---------------------

int RelationType(individual *ind, individual *ind2);
           // Returns the relationship between ind and ind2
           // NOT GENERAL


// ---------------------
//  Genotype Elimination
// ---------------------

individual *GenotypeElimination(char *pedigree, int markernum);      
           // Returns a copy of pedigree with genotype elimination
           // for marker markernum
           // The marker list holds possible ordered genotypes for each
           // individual

individual *XGenotypeElimination(char *pedigree, int markernum);      
           // Returns a copy of pedigree with genotype elimination
           // for marker markernum
           // The marker list holds possible ordered genotypes for each
           // individual

double ReduceGenotypes(individual *indlist);
           // Reduces a list of ordered genotypes by genotype elimination
           // indlist is modified !
           // Returns the number of combinations in the resulting dataset


// ----------------
//  IBD Calculation
// ----------------

OLDMATRIX *PedigreeIBD (individual *pedigree);
           // Inputs a genotype-eliminated list of individuals
           // from a single family
           // Returns an IBD-matrix of exact calculations

OLDMATRIX *CalcPIHat(individual *indlist, int markernum);
           // Calculates the complete PI Hat matrix by the smartest method

double CalcLikelihoodOfPedigree(individual *indlist);
           // Calculates the likelihood of a pedigree with
           // ordered, non-missing genotypes
           // Requires correct allele frequencies!


#endif
