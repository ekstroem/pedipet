/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * Contains all data structures af basic functions
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


#ifndef _PEDIPET_H_
#define _PEDIPET_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "defs.h"
#include "nrutil.h"

static double sqrarg;
#define	 addlist2(l,p)	       (*l = p, l = &p->next)
#define	 min(a,b)              ((a) < (b) ? (a) : (b))
#define	 max(a,b)              ((a) > (b) ? (a) : (b))
#define  absolut(a)            ((a) < 0.0 ? (-a) : (a))
#define	 xisdigit(c)	       ((c) == '-' || ((c) >= '0' && (c) <= '9'))
#define	 addptr(p,i)	       ((void *)(((char *)p) + i))
#define  SQR(a)                ((sqrarg=(a))==0.0 ? 0.0 : sqrarg*sqrarg)
#define  strcmpl(s1,s2)        (strcasecmp(s1,s2))
#define  forind                for (ind = individuals; ind; ind = ind->next)

#define  MISSING               -1
#define  OPTIONSFILE           "options.pp"
#define  SEEDFILE              "seed.pp"
#define  ERRORFILE             "errors.pp"
#define  DESCFILE              "desc.pp"
#define  HEFILE                "haselst.s"
#define  SYSMISS               "."       // Symbol for the system missing value
#define  MATCHKEYWORD          3  // Number of letters used for command mathcing

#include "enums.h"
#include "structure.h"


/*

  Global variables

*/

#ifndef __cplusplus
FILE *F;
individual *individuals;
#endif

char buf[BUFFERSIZE];
char buf2[BUFFERSIZE];
int rndno;
int usedmap;
int UseColour;

// Holds the populations allele frequency
double allelefreq[MAXALLELES];

int nChromosome;
int outi;
char outbuf[256];

optionlist *options;
namelist *pedigrees;
DataDefinitions *datainfo;
IntegerList *IntList;
FreqList *allfreq;      // List of allele frequencies for each marker

namelist *markernames;  // List of marker names
namelist *traitnames;   // List of trait names
int      *order;        // Vector of marker order
double   *distance;     // Vector of marker distances
int      *invorder;     // Inverse vector of marker order

double RES1, RES2, RES3;

// Model variables
int ModelTrait;
int nModelCovariates;
int ModelCovariates[MAXTRAITS];

//-----------------------
// Info about individuals
//-----------------------

/*
int founder(individual *ind);                     // Is person a founder
int fullsibs(individual *ind1, individual *ind2);  // Ind and ind2 full sibs?
int halfsibs(individual *ind1, individual *ind2);  // Ind and ind2 half sibs?
int SamePedigree(individual *ind1, individual *ind2);  // Same pedigree?
int homozygous (individual *ind, int n);
*/


//-----------------------


// Rutines to write error messages


int WriteErrorMsg(char *s);

// Returns 1 if no data on father and mother exists
int founder(individual *ind);

// Returns 1 if ind1 and ind2 have the same parents, 0 otherwise
int fullsibs(individual *ind1, individual *ind2);

// Returns 1 if ind1 and ind2 are half sibs, 0 otherwise
int halfsibs(individual *ind1, individual *ind2);

int SamePedigree(individual *ind1, individual *ind2);

int numberofpedigrees (void);

// Gets marker n from individual ind
markerlist *markernumber (individual *ind, int n);

// Gets info from marker n 
FreqList *FrequencyNumber (int n);

// converts a string to an integer and negative values to 0
int atoip (char *s);


// Gets trait n from individual ind
double trait(individual *ind, int n);
quanttrait *traitvalue(individual *ind, int traitnumber);
int IsTraitMissing(individual *ind, int n);

void *cmalloc (int n);

// Returns 1 if added else 0
// Add p1 to list
int adduniquelist (void *l1,void *p1, int structsize);

// Only adds if it does not already exist in list
void adduniqueidlist (IDlist **l1, individual *p1);


/* ---------------------------- */

/* Checking if the individual is homoz at marker n */
int homozygous (individual *ind, int n);

/* Checking if the individual is heteroz at marker n */
int heterozygous (individual *ind, int n);

individual *findperson (char *s);
individual *FindIndListMember(individual *l, char *id);
individual *FindIndListMemberFromPedigree(individual *l, char *id, char *ped);

void RemoveMarker(individual *indlist, int mkrnum);

/*

Returns the largest allelenumber of a typed person - 0 otherwise

*/

int istypedformarker(individual *ind, int markernum);


int maxallelenumber (int number);


/*

------ INFO ABOUT THE DATASET ------

*/

/*

  Returns the number of different alleles for marker 'markernum'

*/

int numberofmarkers (void);

int NumberOfAlleles(int markernum);


//  Makes a list of the pedigrees
namelist *MakePedigreeList(individual *indlist);

// Makes pedigree list of total dataset
void MakeCompletePedigreeList();


// Creates a copy of indlist
individual *CopyIndividualList(individual *indlist);

// Creates a copy of indlist of all people in pedigree pedname
individual *SelectPedigree(char *pedname, individual *indlist);

void cfopen (char *filename,char *mode);

char *igetstr (char *s1);

char *getstr (void);

void getbuf (void);

char *itostr (int i);

char *cstrdup (char *s1);

int geti (void);

double getf (void);

int strpcmp (const void *s1, const void *s2);

// Gets the n'th name in a namelist
char *GetName(namelist *checklist, int number);

void CopyAlleleFreq(int marker);
IDlist *SortIndividuals(individual *indlist);


/*

 Rutines for handling quantitative traits

*/

// Return 1 if the trait is missing, 0 otherwise, -1 if trait does not exists
int traitmiss(individual *ind, int traitnumber);

// Returns a pointer to a quantitative trait
quanttrait *traitvalue(individual *ind, int traitnumber);
int numberoftraits (void);

//
// Returns the number of pedigrees in dataset
//


double MapFunction(double rf);

double InverseMapFunction(double dist);

//  Returns the recombination fraction between marker 1 and marker 2
double InterMarkerRecomb(int marker1, int marker2);

// Calculates the distance (in cM) from position to marker
double MarkerDistance(double position, int markernum);

double ChromosomeLength(void);

void InitializeFrequencies(int nmarkers);


void FreeIndividualList(individual *indlist);
void FreeNameList(namelist *list);

// Deletes the pedigree data
void DeletePedigreeData (void);

// Deletes the parameter data
void DeleteParameterData (void);

// Orders the markers according to string and fills out the rest
void OrderMarkers(char *string);

void TextColour(int colour);

void BackgroundColour(int colour);

char *InputLine(char *buffer, int bufsize);

int FamilySize(individual *indlist, char *pedid);


#undef debug


#ifdef __cplusplus
}
#endif


#endif
