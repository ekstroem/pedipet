/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * This file sets the appropriate global
 * structures for pedipet.
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


#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "defs.h"
#include "enums.h"

struct individual;
typedef struct individual individual;

typedef struct quanttrait
{
  struct quanttrait *next;
  double value;
  char sysmiss;  // System missing value
} quanttrait;

typedef struct markerlist
{
  struct markerlist *next;
  int allele1;
  int allele2;

} markerlist;

typedef struct ParentalGenotypesList
{
  struct ParentalGenotypesList *next;
  // Father 1, father 2, mother 1, mother 2
  int allele[4];  
} ParentalGenotypesList;

typedef struct MarkerTree
{
  struct MarkerTree *next;
  struct MarkerTree *parent;
  struct MarkerTree *child;
  
  int allele[2];
  int homozygote;
  int fromparent;

} MarkerTree;


typedef struct IDlist
{
  struct IDlist *next;
  struct individual *ind;
} IDlist;

// Reverse ID list
typedef struct RIDlist
{
  struct RIDlist *next;
  struct RIDlist *prev;
  struct individual *ind;
} RIDlist;


struct individual
{
  struct individual *next;

  individual *father;
  individual *mother;
  IDlist *sib;
  IDlist *mate;
  IDlist *offspring;

  markerlist *marker;
  quanttrait *qttrait;

  char sex;                  // Individuals sex 
  char pedigree[NAMESIZE];   // Family code 
  char name[NAMESIZE];       // Individuals name 
  char id[NAMESIZE];         // Individuals id (unique!) 
  char idinped[NAMESIZE];    // Internal id in family (default = id) 

  unsigned int globalid;     // Global, internal ID number
  unsigned int localid;      // Local internal ID number

  // Temporary variables for use in different rutines
  char tmpstr1[NAMESIZE];
  char tmpstr2[NAMESIZE];
  int  tmpint1;
  int  tmpint2;
  int  tmpint3;
  int  tmpint4;
  double tmpdouble1;
  //  ParentalGenotypesList *pParentGenotypes;
};

typedef struct namelist
{
  struct namelist *next;
  char name[NAMESIZE];
} namelist;

typedef struct pedigreelist
{
  struct pedigreelist *next;
  char name[NAMESIZE];
  int combinations;
} pedigreelist;

typedef struct FreqList
{
  struct FreqList *next;
  int num_alleles;
  double frequency[MAXALLELES];
} FreqList;

typedef struct optionlist
{
  char Index[MAXOPTIONS];  // Array of pointers to options

  /*  
  char IDcasesensitive;  // Are IDs case sensitive?        
  char ConvertLinkageID; // Convert linkage IDs to fam-id  
  char ReduceAlleles;    // Automatically reduce alleles   
  char NewErrorFile;     // Delete the errorfile when importing new data  
  char NoBound;          // No boundary restriction when maximizing  
  char WriteErrorFile;   // Should the errorfile be used at all  
  char IterationInfo;    // Should the errorfile be used at all 
  char colour;           // Use ANSI colour 
  char xlinked;          // Chromosome is x-linked
  */
} optionlist;

typedef struct IntegerList
{
  int length;
  int *values;
} IntegerList;

typedef struct FloatList
{
  struct FloatList *next;
  double *value;
} FloatList;

// Holds info about each datatype
typedef struct DataDefinitions
{
  struct DataDefinitions *next;
  int datatype;
  char name[25];
  FreqList  *allelefrequency;
  double distance;
} DataDefinitions;


typedef struct StructModelDef
{
  int nTraits;           // Number of traits in model
  int nCovariates;       // NUmber of covariates

  int traits[MAXTRAITS];
  int covariates[MAXTRAITS];
  int fixed[MAXTRAITS];
} ModelDef;


#ifdef __cplusplus
}
#endif


#endif
