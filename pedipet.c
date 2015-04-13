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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stddef.h>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#include <ctype.h>


#include "lists.h"
#include "distr.h"
#include "fileio.h"
#include "pedipet.h"
#include "simulate.h"

#include "bridge.h"
//#include "allfreq.h"
#include "linkage.h"

#include "mathfct.c"
#include "haselst.c"
#include "ibd.c"
#include "multi.c"

//#define norandom

// Possible datafile keywords
enum
{
  K_NAME,
  K_MARKER,
  K_QUANT,
  K_QUALI,
  MAXDATATYPES,
};


char *optionnames[MAXOPTIONS] =
{
  "colour",
  "convertlinkageid",
  "dominance", 
  "idcase",
  "iterationinfo",
  "newerrorfile",
  "nobound",
  "noseed",
  "reducealleles",
  "savedir",
  "seed",
  "useerrorfile",
  "xlinked",
};

// Command keywords
enum
{
  C_ARTIKEL1,
  C_CEPH,
  C_CONSISTENCY,
  C_CHROM,
  C_CLAUS,
  C_CLONE,
  C_COUNT,
  C_COVARIATES,
  C_CRIMAP,
  C_DATA,
  C_DELETE,
  C_DISTANCE,
  C_ECHO,
  C_ELLEN,
  C_EQUAL,
  C_ERROR,
  C_EXPAND,
  C_FATHER,
  C_FIX,
  C_FREQ,
  C_GENOTYPES,
  C_GH,
  C_GXE,
  C_HASEMAN,
  C_HELP,
  C_HERITABILITY,
  C_IBD,
  C_IBS,
  C_IMPORT,
  C_INDIVIDUAL,
  C_INPUT,
  C_KILL,
  C_LICENSE,
  C_LOAD,
  C_MAKE,
  C_MAP,
  C_MARKER,
  C_MENDEL,
  C_MERLIN,
  C_MIBD,
  C_MISSING,
  C_MIXTURE,
  C_ML,
  C_MOTHER,
  C_MP,
  C_MULTI,
  C_NAME,
  C_OPTION,
  C_ORDER,
  C_PARAM,
  C_PEDIGREE,
  C_PHI,
  C_POLYGENIC,
  C_QUIT,
  C_RECODE,
  C_RELPAIR, 
  C_REMOVE,
  C_SCAN,
  C_SEED,
  C_SELECT,
  C_SHOW,
  C_SIBS,
  C_SIMULATE,
  C_SOLAR,
  C_SPECIAL,
  C_STENO,
  C_TEST,
  C_TRAIT,
  C_VC,
  C_XLINK,
  MAXKEYWORDS
};

char *keywords[MAXKEYWORDS] =
{
  "artikel",
  "ceph",
  "check",
  "chromosome",
  "claus",
  "clone",
  "count",
  "covariates",
  "crimap",
  "data",
  "delete",
  "distance",
  "echo",
  "ellen",
  "equal",
  "error",
  "expand",
  "father",
  "fix",
  "freq",
  "genotypes", 
  "gh",
  "gxe",
  "haseman",
  "help",
  "heritability", 
  "ibd",
  "ibs",
  "import",
  "individual",
  "input",
  "kill",
  "license",
  "load",
  "make",
  "map",
  "marker",
  "mendel",
  "merlin",
  "mibd",
  "missing",
  "mixture",
  "ml",
  "mother",
  "mp",
  "multi",
  "name",
  "option",
  "order",
  "param",
  "pedigree",
  "phi",
  "polygenic",
  "quit",
  "recode",
  "relpair", 
  "remove", 
  "scan",
  "seed",
  "select",
  "show",
  "sibs",
  "simulate",
  "solar",
  "special",
  "steno",
  "test",
  "trait",
  "vc",
  "xlinked",
};

char *sexnames[MAXSEX] =
{
  "Unknown",
  "Male",
  "Female",
};

char *mapfctnames[MAXMAP] =
{
  "Haldane",
  "Kosambi",
};


char *leading (char *s)
{
  s[0]=toupper(s[0]);
  return(s);
}

char *noleading (char *s)
{
  s[0]=tolower(s[0]);
  return(s);
}

int keywordcmp (const void *s1, const void *s2)
{
	return strncasecmp (*(char **)s1,*(char **)s2, MATCHKEYWORD);
}


int findoptionskeyword (char *s)
{
   char **sp;

   sp = (char **) bsearch (&s, optionnames, MAXOPTIONS, sizeof s, keywordcmp);
   if (sp == 0)
   {
      return -1;
   }
   return sp - optionnames;
}


int findstr (char **v,char *s,int n)
{
	int i;

	for (i = 0; i != n; i++)
		if (!strcmpl (v[i],s))
			return i;

	return -1;
}

// -----------------------------------------------------
// Functions from pedipet.h
// -----------------------------------------------------

// Rutines to write error messages


int WriteErrorMsg(char *s)
{
  FILE *f;

  time_t now;
  time(&now);

  if (options->Index[O_USEERRORFILE])
  {
    strcpy(buf,asctime(localtime(&now)));
    buf[24] = ':';

    if ((f = fopen (ERRORFILE,"a")))
    {
      fprintf(f, "%s %s",buf,s);
      fclose(f);
    }
  }
  return 0;
}


// Returns 1 if no data on father and mother exists
int founder(individual *ind)
{
  if (!(ind->father) && !(ind->mother))
    return 1;
  return 0;
}


// Returns 1 if ind1 and ind2 have the same parents, 0 otherwise
// AND are two different individuals
int fullsibs(individual *ind1, individual *ind2)
{
  if ((ind1->father==ind2->father) && (ind1->mother == ind2->mother) &&
      ind1->father && ind2->mother && ind1 != ind2)
    return 1;

  return 0;
}

// Returns 1 if ind1 and ind2 are half sibs, 0 otherwise
int halfsibs(individual *ind1, individual *ind2)
{
  if (((ind1->father==ind2->father) && (ind1->mother != ind2->mother)) ||
      ((ind1->father!=ind2->father) && (ind1->mother == ind2->mother)))
    return 1;
  return 0;
}

int SamePedigree(individual *ind1, individual *ind2)
{
  if (!strcmpl(ind1->pedigree, ind2->pedigree))
    return 1;
  return 0;
}


// Gets marker n from individual ind
markerlist *markernumber (individual *ind, int n)
{
  int count;
  markerlist *marker;

  count = 1;
  for (marker = ind->marker; marker; marker = marker->next)
  {
    if (count==n)
      return marker;
    count++;
  }
  printf("WARNING: Truing to access a non-existing marker\n");
  return 0;
}


// Gets info from marker n 
FreqList *FrequencyNumber (int n)
{
  FreqList *fl;
  int count;

  count = 1;
  for (fl = allfreq; fl; fl = fl->next)
  {
    if (count==n)
      return fl;
    count++;
  }
  printf("WARNING: Internal error: referencing illegal marker (%d)\n",n);
  return 0;
}


// converts a string to an integer and negative values to 0
int atoip (char *s)
{
        int n;
 
        n = atoi (s);
 
        if (n < 0)
          n = 0;
 
        return n;
}



// Gets trait n from individual ind
double trait(individual *ind, int n)
{
  int count;
  quanttrait *qt;

  count = 1;
  for (qt = ind->qttrait; qt; qt = qt->next)
  {
    if (count==n)
      return qt->value;
    count++;
  }
  return MISSING;
}


// Return 1 if missing, 0 if ok and -1 if not found
int IsTraitMissing(individual *ind, int n)
{
  int count;
  quanttrait *qt;

  for (qt = ind->qttrait, count = 1; qt; qt = qt->next, count++)  {
    if (count==n) {
      if (qt->sysmiss)
	return 1;
      return 0;
    }
  }
  return -1;
}


void *cmalloc (int n)
{
	void *p;

	if (n == 0)
		n = 1;

	p = malloc (n);

	if (p == 0)
	{
		puts ("Out of memory.");
		exit (1);
	}

	return p;
}

// Returns 1 if added else 0
// Add p1 to list
int adduniquelist (void *l1,void *p1, int structsize)
{
	list **l;
	list *p,*q;
        int code;

	l = (list **) l1;
	p = (list *) p1;

        // Does the thing exist?
        code = 0;

	p->next = 0;

	if (*l)	{
	  for (q = *l; q->next; q = q->next) {
	    // If element already exists in list
	    if (!memcmp(q, p, structsize)) {
	      code = 1;
	    }
	  }
	  // If item does not exist in list then add to list
	  q->next = p;
	  return 0;
	}
	else
	  *l = p;	
   return 1;
}



// Only adds if it does not already exist in list
void adduniqueidlist (IDlist **l1, individual *p1)
{
  IDlist *p, *idlist;
  int addid;

  addid = 1;

  for (p = *l1; p; p = p->next)
    if (!strcmpl(p->ind->id, p1->id))
      addid=0;

  if (addid)
  {
    idlist = (IDlist *) cmalloc (sizeof (IDlist));
    memset (idlist, 0, sizeof (IDlist));
    idlist->ind = p1;
    addlist(l1, idlist);
  }
}

/* ---------------------------- */

/* Checking if the individual is homoz at marker n */
int homozygous (individual *ind, int n)
{
  markerlist *marker;
  
  marker = markernumber(ind, n);

  if (marker->allele1 == marker->allele2 && marker->allele1>0)
    return 1;
  else if (marker->allele1 == 0)
      return -1;

  return 0;
}

/* Checking if the individual is heteroz at marker n */
int heterozygous (individual *ind, int n)
{
  int res;
  
  res = homozygous(ind,n);
  if (res<0)
    return res;
  else
    return (1 - res);
}

individual *findperson (char *s)
{
  individual *indi;

  for (indi = individuals; indi; indi = indi->next)
    if (!strcmpl(indi->id,s))
      return indi;
  return 0;
}

//
// Finds an individual in a list l
//
individual *FindIndListMember(individual *l, char *id)
{
  individual *ind;
  
  for (ind = l; ind; ind = ind->next)
  {
    if (!strcmpl(ind->id, id))
      return ind;
  }
  return 0;
}


//
// Finds an individual in a list l from pedigree ped
//
individual *FindIndListMemberFromPedigree(individual *l, char *id, char *ped)
{
  individual *ind;
  
  for (ind = l; ind; ind = ind->next)
  {
    if (!strcmpl(ind->id, id) && !strcmpl(ind->pedigree, ped))
      return ind;
  }
  return 0;
}




//
// Finds an individual in a list l
//
individual *FindIndListMemberLocalID(individual *l, int localid)
{
  individual *ind;
  
  for (ind = l; ind; ind = ind->next) {
    if (ind->localid==localid)
      return ind;
  }
  return 0;
}



/*

Returns the largest allelenumber of a typed person - 0 otherwise

*/

int istypedformarker(individual *ind, int markernum)
{
  markerlist *marker;

  marker = markernumber(ind,markernum);
  if (!marker->allele1 || !marker->allele2)
    return 0;
  else
    return max(marker->allele1,marker->allele2);
}


/*

------ INFO ABOUT THE DATASET ------

*/

/*

  Returns the number of different alleles for marker 'markernum'

*/

// NumberOfAllelesPresent returns the number of different alleles 
// (not counting the '0' allele) present in the dataset for the marker

int NumberOfAllelesPresent(int markernum)
{
  individual *ind;
  markerlist *marker;
  int allarray[MAXALLELES], i, j;

  memset(allarray, 0, sizeof(allarray));

  forind
  {
    marker = markernumber(ind, markernum);
    allarray[marker->allele1] = 1;
    allarray[marker->allele2] = 1;
  }

  /* Now counting */
  j = 0;
  for (i=1; i<MAXALLELES; i++)
  {
    if (allarray[i])
      j++;
  }

  return(j);
}

// NumberOfAlleles returns the number of different alleles 
// (not counting the '0' allele) possible for a marker in the dataset
// based on the frequency info

int NumberOfAlleles(int markernum)
{
  FreqList *fl;
  int res;
  //  int j;

  fl = FrequencyNumber(markernum);
  res = fl->num_alleles;

  // Remove the alleles that have frequency 0
  //  for (j=1; j<=fl->num_alleles; j++)
  //    if (fl->frequency[j] ==0)
  //      res--;
  
  return (res);
}



/*
  Makes a list of the pedigrees
*/

namelist *MakePedigreeList(individual *indlist)
{
  individual *ind;
  namelist *pedinfo, *pedigree, *retlist;
  int newpedigree;

  retlist = NULL;

  for (ind = indlist; ind; ind = ind->next)
  {
    newpedigree=1;
    for (pedigree = retlist; pedigree; pedigree=pedigree->next)
    {
      if (!strcmpl(ind->pedigree,pedigree->name))
        newpedigree=0;
    }
    // Newpedigree - add to list
    if (newpedigree)
    {
      pedinfo = cmalloc (sizeof (namelist));
      memset (pedinfo, 0, sizeof (namelist));

      strcpy(pedinfo->name, ind->pedigree);
      addlist(&retlist, pedinfo);
    }
  }
  return retlist;
}

// Makes pedigree list of total dataset
void MakeCompletePedigreeList()
{
  if (pedigrees) {
    freelist(pedigrees);
    pedigrees = NULL;
  }

  pedigrees =  MakePedigreeList(individuals);
}


// Opens a file
void cfopen (char *filename, char *mode) {

  F = fopen (filename,mode);
  
  if (F == 0) {
    printf ("Can't open file %s in mode %s.\n",filename,mode);
    exit (1);
  }
}

char *igetstr (char *s1)
{
	int i;
	static char *s;
	static char buf[1256];

	if (s1)
		s = s1;
	while (*s == ' ')
		s++;
	i = 0;

	while (*s && *s != ' ')
	{
		buf[i] = *s;
		//		if (*s == '_')
           	//			buf[i] = ' ';
		s++;
		i++;
	}

	buf[i] = 0;
	return buf;
}

char *getstr (void)
{
	return igetstr (0);
}

void getbuf (void)
{
  int i;
  int c;
  
  i = 0;
  
  for (;;)
    {
      c = fgetc (F);

      if (c == EOF)
	{
	  buf[0] = EOF;
	  return;
	}

      if (c == '\n')
	{
	  buf[i] = 0;
	  return;
	}

      if (i == sizeof buf - 1)
	{
	  buf[i] = 0;
	  while (c != EOF && c != '\n')
	    c = fgetc (F);
	  if (c == EOF)
	    buf[0] = EOF;
	  return;
	}

      buf[i++] = c;
    }
}

char *itostr (int i)
{
	static char buf[20];

	sprintf (buf,"%d",i);
	return buf;
}

char *cstrdup (char *s1)
{
	char *s;

	if (!s1[0])
		return 0;

	s = (char *) cmalloc (strlen (s1) + 1);
	strcpy (s,s1);
	return s;
}

int geti (void)
{
	return atoip (getstr ());
}

double getf (void)
{
	return atof (getstr ());
}

int strpcmp (const void *s1, const void *s2)
{
	return strcmpl (*(char **)s1,*(char **)s2);
}

// Gets the n'th name in a namelist
char *GetName(namelist *checklist, int number)
{
  int i;
  namelist *nl;

  i = 1;
  for (nl = checklist; nl; nl = nl->next)
  {
    if (i == number)
      return(nl->name);
    i++;
  }
  return(NULL);
}


// Gets the n'th name in a namelist
void SetName(namelist *checklist, int number, char *string)
{
  int i;
  namelist *nl;

  i = 1;
  for (nl = checklist; nl; nl = nl->next)
  {
    if (i == number)
      break;
    i++;
  }
  strcpy(nl->name, string);
}


void CopyAlleleFreq(int marker)
{
  int j;
  FreqList *fl;

  fl = FrequencyNumber(marker);
  for (j=1; j<=fl->num_alleles; j++)
    allelefreq[j] = fl->frequency[j];
}



/*

 Rutines for handling quantitative traits

*/

// Return 1 if the trait is missing, 0 otherwise, -1 if trait does not exists
int traitmiss(individual *ind, int traitnumber)
{
  quanttrait *qt;
  int count = 1;

  for (qt = ind->qttrait; qt; qt = qt->next)
  {
    if (count == traitnumber)
    {
      if (qt->sysmiss)
        return 1;
      else
        return 0;
    }
    count++;
  }
  return -1;
}

// Returns a pointer to a quantitative trait
quanttrait *traitvalue(individual *ind, int traitnumber)
{
  quanttrait *qt;
  int count;

  for (qt = ind->qttrait, count = 1; qt; qt = qt->next, count++) {
    if (count == traitnumber)
      return qt;
  }
  return 0;
}

void ClearString(char *string) {
  memset (string, 0, sizeof (char));
}

//
// Returns the number of pedigrees in dataset
//


double MapFunction(double rf)
{
  if (rf<0 || rf>0.5)
  {
    printf("Error: recombination fraction not in [0,.5] ==> %f\n",rf);
    return -1;
  }
  if (usedmap == M_HALDANE)
    return (-log(1-2*rf)/2);
  else
    return(log((1+2*rf)/(1-2*rf))/4);
}

//
// Returns the recombination fraction from a distance (in Morgans)
//
double InverseMapFunction(double dist)
{
  if (dist<0)
  {
    printf("Error: distance can not be negative\n");
    return -1;
  }
  if (usedmap == M_HALDANE)
    return ((1-exp(-2*dist))/2);
  else
    return((exp(4*dist)-1)/(2*(1+exp(4*dist))));
}

//  Returns the recombination fraction between marker 1 and marker 2
double InterMarkerRecomb(int marker1, int marker2)
{
  int minmar, maxmar, i;
  double dist;

  minmar = min(marker1, marker2);
  maxmar = max(marker1, marker2);

  if (minmar<1 || maxmar>numberofmarkers())
  {
    printf("ERROR: Internal error in InterMarkerRecomb\n");
    return 0;
  }

  if (minmar==maxmar)
    return 0.0;

  dist = 0.0;
  for (i=minmar; i<maxmar; i++)
  {
    dist += MapFunction(distance[i]);    
  }
  return(InverseMapFunction(dist));
}

// Calculates the distance (in cM) from position to marker
double MarkerDistance(double position, int markernum)
{
  int i;
  double dist;

  if (markernum<1 || markernum>numberofmarkers())
  {
    printf("Internal error: Marker Distance\n");
    return 0;
  }

  dist  = 0.0;
  for (i=1; i< markernum; i++)
    dist += 100*MapFunction(distance[i]);
  return (dist - position);
}

//
// Returns the length of the chromosome in Morgans!
//
double ChromosomeLength() {
  double dist = 0.0;
  int i;

  for (i = 1; i < numberofmarkers(); i++) {
    dist += MapFunction(distance[i]);    
  }

  return(dist);
}

void InitializeFrequencies(int nmarkers)
{
  FreqList *fl;
  int i;

  freelist(allfreq);
  allfreq = NULL;

  for (i=1; i<=nmarkers; i++)  {
    fl = (FreqList *) cmalloc (sizeof (FreqList));
    memset (fl, 0, sizeof (FreqList));

    addlist(&allfreq, fl);
  }
}

void RemoveIndividual(individual *indlist, individual *ind) 
{
  IDlist *person, *person2;
  int code;

  freelist(ind->marker);
  freelist(ind->qttrait);

  ind->marker = NULL;
  ind->qttrait = NULL;

  // Fixes mates
  for (person = ind->mate; person; person = person->next) {
    for (person2 = person->ind->mate; person2; person2 = person2->next) {
      if (person2->ind == ind) {
	removelist(&person->ind->mate, person2);
	break;
      }
    }
  }
  freelist(ind->mate);
  ind->mate = NULL;


  // Remove siblings
  for (person = ind->sib; person; person = person->next) {
    for (person2 = person->ind->sib; person2; person2 = person2->next) {
      if (person2->ind == ind) {
	removelist(&person->ind->sib, person2);
	break;
      }
    }
  }
  freelist(ind->sib);
  ind->sib = NULL;

  // Fixes offspring
  // Removes the links for the offspring
  for (person = ind->offspring; person; person = person->next) {
    if (person->ind->father == ind)
      person->ind->father = NULL;
    else 
      person->ind->mother = NULL;
  }
  freelist(ind->offspring);
  ind->offspring = NULL;

  // Fixes parents
  // Removing an offspring, so remove from parents lists
  if (ind->father) {
    for (person = ind->father->offspring; person; person = person->next) {
      if (person->ind == ind) {
	removelist(&ind->father->offspring, person);
	break;
      }
    }    
  }
  if (ind->mother) {
    for (person = ind->mother->offspring; person; person = person->next) {
      if (person->ind == ind) {
	removelist(&ind->mother->offspring, person);
	break;
      }
    }    
  }

  // Check that parents still are mates  
  code = 0;
  if (ind->father) {
    for (person = ind->father->offspring; person; person = person->next) {
      if (person->ind->mother == ind->mother) {
	code = 1;
	break;
      }
    }
    // No common offspring
    if (code == 0) {
      for (person = ind->father->mate; person; person = person->next) {
	if (person->ind == ind->mother) {
	  removelist(&ind->father->mate, person);
	  break;
	}
      }
    }
  }
  code = 0;
  if (ind->mother) {
    for (person = ind->mother->offspring; person; person = person->next) {
      if (person->ind->father == ind->father) {
	code = 1;
	break;
      }
    }
    // No common offspring
    if (code == 0) {
      for (person = ind->mother->mate; person; person = person->next) {
	if (person->ind == ind->father) {
	  removelist(&ind->mother->mate, person);
	  break;
	}
      }
    }
  }

#ifdef debug
  printf("Removing %s from pedigree %s\n", ind->id, ind->pedigree);
#endif
  removelist(&indlist, ind);
}


void FreeIndividualList(individual *indlist)
{
  individual *ind;

  for (ind = indlist; ind; ind = ind->next)
  {
    freelist(ind->marker);
    freelist(ind->qttrait);
    freelist(ind->sib);
    freelist(ind->mate);
    freelist(ind->offspring);
  }  
  freelist(indlist);
  indlist = NULL;
}

void FreeNameList(namelist *list)
{
  freelist(list);
  list = NULL;
}


/******************
 *
 * Deletes an individual list
 *
 ******************/
// Deletes the pedigree data
void DeleteIndividualList (individual *indlist)
{
  individual *ind;

  // Frees individual info
  for (ind=indlist; ind; ind=ind->next)
  {
    freelist(ind->marker);
    freelist(ind->qttrait);
    freelist(ind->sib);
    freelist(ind->mate);
    freelist(ind->offspring);
  }

  freelist(indlist);
  indlist = NULL;
}




// Deletes the pedigree data
void DeletePedigreeData (void)
{
  individual *ind;

  // Frees individual info
  forind
  {
    freelist(ind->marker);
    freelist(ind->qttrait);
    freelist(ind->sib);
    freelist(ind->mate);
    freelist(ind->offspring);
  }

  freelist(individuals);
  individuals = 0;
}

// Deletes the parameter data
void DeleteParameterData (void)
{
  int nmarkers;

  nmarkers = listlen(markernames);

  if (nmarkers > 0)
    free_ivector(order, 1, nmarkers);
  if (nmarkers > 1)
  free_vector(distance, 1, nmarkers-1);

  freelist(allfreq);
  allfreq = NULL;

  freelist(markernames);
  markernames = NULL;
  freelist(traitnames);
  traitnames = NULL;
}

// Orders the markers according to string and fills out the rest
void OrderMarkers(char *string)
{
  int i, k, code;
  int usedorder[MAXMARKERS];

  // Sets the new order
  memset(usedorder, 0, sizeof(usedorder));

  i = 1;
  code = atoip(igetstr(string));

  while (code)
  {
    if (code <= numberofmarkers())
    {
      // Has this marker been placed?
      if (!usedorder[code])
      {
        order[i] = code;
        invorder[code] = i;
        i++;
      }

      usedorder[code] = 1;
    }

    code = atoip(getstr());
  }

  // Should now enter the remaining orders
  for (k=1; k<=numberofmarkers(); k++)
  {
    if (!usedorder[k])
    {
      usedorder[k] = 1;
      order[i] = k;
      invorder[k] = i;
      i++;
    }
  }
}


void TextColour(int colour)
{
  if (UseColour)
    printf("\x1b[%dm",colour+30);
}

void BackgroundColour(int colour)
{
  if (UseColour)
    printf("\x1b[%dm",colour+40);
}


// -----------------------------------------------------

int numberoftraits (void)
{
  quanttrait *qt;
  individual *ind;
  int count = 0;

  ind = individuals;

  if (ind)
  {
    for (qt = ind->qttrait; qt; qt = qt->next)
      count++;
  }
  return count;
}


int numberofmarkers (void)
{
//  markerlist *marker;
//  individual *ind;
//  int count = 0;

//  ind = individuals;


//  if (ind)
//  {
//    for (marker = ind->marker; marker; marker = marker->next)
//      count++;
//  }
//  return count;
  return listlen(allfreq);
}


//
// Is the person for any of the markers genotyped? 1 = yes, 0 = no
// 
int IsGenotyped(individual *ind)
{
  markerlist *marker;

  for (marker = ind->marker; marker; marker = marker->next)
  {
    if (marker->allele1 > 0 || marker->allele2>0) {
      return 1;
    }
  }
  return 0;
}

int IsPhenotyped(individual *ind)
{
  quanttrait *qt;

  for (qt = ind->qttrait; qt; qt = qt->next)
  {
    if (!qt->sysmiss) {
      return 1;
    }
  }
  return 0;
}


// Returns the number of missing alleles for the marker
int missingalleles (markerlist *marker)
{
  return ((marker->allele1 ? 0 : 1) + (marker->allele2 ? 0 : 1));
}


// Returns the number of pedigrees in the dataset
int numberofpedigrees (void)
{
  int count, code;
  individual *ind;
  namelist *pedi, *ped;

  count = 0;
  pedi = 0;
  for (ind = individuals; ind; ind = ind->next)
  {
    code = 1;
    for (ped = pedi; ped; ped = ped->next)
    {
      if (!strcmp(ped->name, ind->pedigree))
        code = 0;
    }
    // New pedigree
    if (code)
    {
      ped = (namelist *) cmalloc (sizeof (namelist));
      memset (ped,0,sizeof (namelist));
      strcpy(ped->name,ind->pedigree);
      addlist(&pedi, ped);
      count++;
    }
  }

  freelist(pedi);
  pedi = NULL;

  return count;
}

// Returns the value of the highest allele for the marker
int maxallelenumber (int number)
{
  markerlist *marker;
  individual *ind;
  int i;

  i = 0;

  forind
  {
    marker = markernumber(ind, number);

    i = max(marker->allele1, i);
    i = max(marker->allele2, i);
  }

  return i;
}

int checkmarker (individual *ind, int n)
{
  markerlist *marker, *fmarker, *mmarker;
  int okallele, notyped;

  marker = markernumber (ind, n);
  fmarker = markernumber (ind->father, n);
  mmarker = markernumber (ind->mother, n);

  notyped = (marker->allele1 ? 1 : 0)+(marker->allele2 ? 1 : 0);

  switch (notyped)
  {
    case 0 : // Nothing typed - matches all
             return 1;
             break;
    case 1 : if (marker->allele1)
               okallele = marker->allele1;
             else
               okallele = marker->allele2;
             // Does the allele fit one of the parents alleles?
             if ((okallele==fmarker->allele1) || (okallele == fmarker->allele2) ||
                (okallele==mmarker->allele1) || (okallele == mmarker->allele2))
                return 1;
             // At least one of the parents has a missing marker
             if (missingalleles(fmarker)+missingalleles(mmarker))
             {
               // If the allele matches on of the
               if (missingalleles(fmarker))
               {
                 if (okallele == atoi (ind->father->tmpstr1) ||
                     okallele == atoi (ind->father->tmpstr2) ||
                     okallele == atoi (ind->mother->tmpstr1) ||
                     okallele == atoi (ind->mother->tmpstr2) )
                   return 1;
               }
               else
               {
               }
             }

             break;
  }


  // If just 1 typed then match
  if ((marker->allele1==0) || (marker->allele2==0))
  {
  }
  // Both markers typed

  // Both parents typed
  if (((marker->allele1 == fmarker->allele1) || (marker->allele1 == fmarker->allele2)) &&
      ((marker->allele2 == mmarker->allele1) || (marker->allele2 == mmarker->allele2)))
    return 1;
  if (((marker->allele2 == fmarker->allele1) || (marker->allele2 == fmarker->allele2)) &&
      ((marker->allele1 == mmarker->allele1) || (marker->allele1 == mmarker->allele2)))
    return 1;

  return 0;
}

// ------------------------------------
//  Allele frequency estimation
// ------------------------------------


void SetEqualFreq(int markernum)
{
  int number, j;
  FreqList *fl;

  fl = FrequencyNumber(markernum);

  number = NumberOfAllelesPresent(markernum);
   
  fl->num_alleles = number;
  for (j=1; j<=number; j++)
    fl->frequency[j] = (double) 1 / number;
 
}

void SetCountFreq(int markernum)
{
  int number, i, count, tmp[MAXALLELES];
  individual *ind;
  markerlist *mkr;
  FreqList *fl;

  for (i=1; i<MAXALLELES; i++)
    tmp[i] = 0;
  count = 0;

  forind
  {
    mkr = markernumber(ind,markernum);
    if (mkr->allele1 && mkr->allele2)
    {
      count +=2;
      tmp[mkr->allele1]++;
      tmp[mkr->allele2]++;
    }
  }

  fl = FrequencyNumber(markernum);
  number = NumberOfAlleles(markernum);
   
  fl->num_alleles = number;
  for (i=1; i<=number; i++)
    fl->frequency[i] = (double) tmp[i] / count;
 
}

void SetInputFreq(int markernum)
{
  int number, i, count, tmp[MAXALLELES];
  double res;
  char minibuf[30];
  FreqList *fl;

  if (1>markernum || markernum>numberofmarkers())
  {
    printf("Illegal marker number\n");
    return;
  }

  for (i=1; i<MAXALLELES; i++)
    tmp[i] = 0;
  count = 0;

  number = NumberOfAlleles(markernum);

  printf("Input frequencies for first %d alleles for marker %d (%s)\n",number-1,markernum,
         GetName(markernames,markernum));

  res = 0.0;
  fl = FrequencyNumber(markernum);
  for (i=1; i<number; i++) {
    printf("  Allele %d frequency: ",i);
    fl->frequency[i] = atof(InputLine(minibuf, 30));
    res += fl->frequency[i];
  }

  if (res>1) {
    printf("WARNING: Sum of first %d alleles is greater than 1.\nSetting last allele frequency to 0 and makes sum to 1\n", number);
    fl->frequency[number] = 0;
    for (i=1; i<number; i++)
      fl->frequency[i] /= res;
  }
  else
    fl->frequency[number] = 1-res;

}

// ------------------------------------

int int_compare(const void *i, const void *j) 
{
  return((*(int *)i)-(*(int *)j));
}

void ReduceAlleles (int number)
{
  individual *ind;
  markerlist *marker;
  int allarray[MAXALLELES], i, count, found;

  memset(allarray,0,sizeof(allarray));

  // First figure out the number of different alleles in the dataset
  count = 0;
  forind {
    marker = markernumber(ind,number);
    // Is the person genotyped
    if ((marker->allele1) || (marker->allele2)) {
      // Check if the alleles already exists
      found = 0;
      for (i = 0; i<MAXALLELES && allarray[i]>0; i++) {
	if (marker->allele1==allarray[i]) {
	  found = 1;
	  break;
	}
      }
      if (!found) {
	allarray[count] = marker->allele1;
	count++;
      }
      
      if (count>MAXALLELES) {
	printf("ERROR: Found more than %d different alleles for marker %d\n", MAXALLELES, number);
	exit(1);
      }

      // Do the same thing for allele2 if it is different from allele1
      if (marker->allele2 != marker->allele1) {
	found = 0;
	for (i = 0; i<MAXALLELES && allarray[i]>0; i++) {
	  if (marker->allele2==allarray[i]) {
	    found = 1;
	    break;
	  }
	}
	if (!found) {
	  allarray[count] = marker->allele2;
	  count++;
	}
      }
      if (count>MAXALLELES) {
	printf("ERROR: Found more than %d different alleles for marker %d\n", MAXALLELES, number);
	exit(1);
      }    
    }
  }

  // Have now set up the table from old to new allele numbers
  // Then convert the original allele values to the new
  // This could have been done in the loop above, but I want to 
  // sort the values in ascending order

  qsort(allarray, count, sizeof(int), int_compare);

  #ifdef debug
  for (i = 0; i< count; i++) {
    printf("%2d %2d\n", allarray[i], i+1);
  }
  #endif

  // Now do the actual replacement
  forind {
    marker = markernumber(ind,number);
    for (i = 0; i<MAXALLELES && allarray[i]>0; i++) {
      if (marker->allele1 == allarray[i]) {
	marker->allele1 = i+1;
	break;
      }
    }
    for (i = 0; i<MAXALLELES && allarray[i]>0; i++) {
      if (marker->allele2 == allarray[i]) {
	marker->allele2 = i+1;
	break;
      }
    }    
  }
}

void inputredall (void) 
{
  int i;

  for (i = 1; i<=numberofmarkers(); i++) {
    ReduceAlleles(i);
  }
  if (numberofmarkers())
    printf("Recoded alleles for all markers\n");
  else
    printf("No marker data in memory\n");
}

void ReadDataset (void)
{
  int i;
  importlinkagepedigree();

  if (options->Index[O_REDUCEALLELES]) {
    inputredall();
  }

  printf("WARNING: Setting allele frequencies to equal values for all markers\n");
  for (i=1; i<=numberofmarkers(); i++)
    SetEqualFreq(i);                    
}

/*
  Calculates the marker heterozygosity based on the typed founders
  in the dataset

*/

double MarkerHeterozygosity(int markernum, double *variance)
{
  individual *ind;
  markerlist *marker;
  double number, different;

  different = 0.0;
  number = 0.0;
  forind
  {
    if (founder(ind))
    {
      marker = markernumber(ind, markernum);
      // If the person has been typed
      if (marker->allele1 && marker->allele2)
      {
        number++;
        if (marker->allele1-marker->allele2)
        {
          different++;
        }
      }
    }
  }
  *variance = (double) different/number * (1.0-((double) different/number)) / number;
  return different/number;
}


//
// Calculates the theoretical Marker polymorphism information content 
//
double MarkerPIC(int markernum)
{
  int i, j;
  double pic;

  FreqList *fl;

  pic = 1;
  fl = FrequencyNumber(markernum);

  for (i=1; i<= fl->num_alleles; i++)
    pic -= fl->frequency[i]*fl->frequency[i];

  for (i=1; i< fl->num_alleles; i++) {
    for (j=i+1; j<= fl->num_alleles; j++) {
      pic -= 2*fl->frequency[i]*fl->frequency[i]*fl->frequency[j]*fl->frequency[j];
    }
  }

  return pic;
}

char *printyesno(int ok)
{
  char *s;


  if (ok)
    s= "yes";
  else
    s = "no";
  return s;
}

void ShowOptions(void)
{
  int i;

  printf("Compiled-in constants:\n");
  printf("  LINESIZE    %4d  []\n", LINESIZE);
  printf("  MAXALLELES  %4d  []\n", MAXALLELES);
  printf("  MAXMARKERS  %4d  [Maximum number of markers]\n", MAXMARKERS);
  printf("  MAXPERSONS  %4d  [Maximum number of alleles pr. marker]\n", MAXPERSONS);
  printf("  MAXTRAITS   %4d  [Maximum number of traits]\n", MAXTRAITS);
  printf("  NAMESIZE    %4d  [Maximum length of names]\n", NAMESIZE);

  printf("Showing options:\n");

  for (i = 0; i<MAXOPTIONS; i++) {
    printf("%-20s : %s\n",optionnames[i], options->Index[i] ? "yes" : "no");
  }
}

void ShowMarkerInfo (int markernum)
{
  individual *ind;
  FreqList *fl;
  double variance, theo;
  int i, nofound;

  if (numberofmarkers()==0)
  {

    printf("No data in memory\n");
    return;
  }

  printf("Marker order:\n");
  for (i = 1; i<=numberofmarkers(); i++)
    printf("%2d    ", order[i]);
  printf("\n");

  //  for (i = 1; i<=numberofmarkers(); i++)
  //    printf("%2d    ", invorder[i]);
  //  printf("\n");


  for (i = 1; i<numberofmarkers(); i++)
    printf("  %4.1f", 100*MapFunction(distance[i]));
  printf("   cM\n[Map Function : %s]\n\n", mapfctnames[usedmap]);

  if (markernum<1 || markernum>numberofmarkers())
  {
    puts("Marker not found");
    return;
  }

  printf("Observed    heterozygosity for marker %-2d (%-12s): %1.3f ", markernum, GetName(markernames,markernum),
  MarkerHeterozygosity(markernum, &variance));
  printf("(%1.3f)\n",sqrt(variance));

  theo = 0.0;
  fl = FrequencyNumber(markernum);
  for (i=1; i<= fl->num_alleles; i++)
    theo += fl->frequency[i]*fl->frequency[i];
  nofound = 0;
  forind
  {
    if (founder(ind))
      nofound++;
  }

  printf("Theoretical heterozygosity for marker %-2d (%-12s): %1.3f (%1.3f)\n", markernum, GetName(markernames,markernum), 1.0-theo, sqrt(theo*(1.-0-theo)/nofound));
  printf("Polymorphism information content (PIC)                 : %1.3f\n", MarkerPIC(markernum));

  printf("Number of alleles: %d\n",fl->num_alleles);
  printf("Alleles   : ");
  for (i=1; i<= fl->num_alleles; i++)
    printf("%5d ",i);
  printf("\n");

  printf("Frequency : ");
  for (i=1; i<= fl->num_alleles; i++)
    printf("%5.3f ",fl->frequency[i]);
  printf("\n");

}

void showindividual (char *string)
{
  individual *ind;
  markerlist *marker;
  quanttrait *qt;
  int i;

  ind = findperson(string);
  if (ind)
  {
    printf("ID          : %s (%s)\n"
	   "Pedigree    : %s\n"
           "Sex         : %s\n"
           "Founder     : %s\n"
           "# offspring : %d\n"
           "# mates     : %d\n"
           "# sibs      : %d\n", 
    ind->id, ind->name, ind->pedigree,
    sexnames[(int)ind->sex], founder(ind) ? "Yes" : "No",
    listlen(ind->offspring),
    listlen(ind->mate),
    listlen(ind->sib));

    printf("Markers:\n");
    for (i = 1; i<= numberofmarkers(); i++)
    {
      marker = markernumber(ind,order[i]);
      printf("   %-12s  %2d %2d\n", GetName(markernames, order[i]), marker->allele1, marker->allele2);
    }

    printf("Traits:\n");
    for (qt = ind->qttrait; qt; qt = qt ->next)
    {
      if (qt->sysmiss)
        printf("   NaN\n");
      else
        printf("   %f\n", qt->value);
    }

  }
  else
    printf("Individual not found\n");
}

//
// Shows info about the dataset
//
void ShowDataset(individual *indlist)
{
  individual *ind, *ind2;
  namelist *Pedigrees;
  int i, j, nObs, nGenotyped, nPhenotyped, nGenoPheno, nFounder, nNonFounder, nPedigrees;
  int nGFounder, nGNonFounder, nPFounder; 
  OLDMATRIX *Rela;
  int nRela[REL_MAXREL], nPhenoPair[REL_MAXREL];

  // Check that data is present
  if ((nObs = listlen(indlist)) == 0)
    return;

  nGenotyped = 0;
  nPhenotyped= 0;
  nGenoPheno = 0;
  nFounder   = 0;
  nPFounder  = 0;
  nGFounder  = 0;
  nGNonFounder=0;
  for (ind = indlist; ind; ind = ind->next)
  {
    if (founder(ind))
      nFounder++;
    if (IsGenotyped(ind))
      nGenotyped++;
    if (IsPhenotyped(ind))
      nPhenotyped++;

    if (IsGenotyped(ind) && IsPhenotyped(ind))
      nGenoPheno++;

    if (IsGenotyped(ind) && founder(ind))
      nGFounder++;

    if (IsPhenotyped(ind) && founder(ind))
      nPFounder++;

  }
  nNonFounder = nObs - nFounder;
  nGNonFounder = nGenotyped - nGFounder;

  Pedigrees = MakePedigreeList(indlist);
  nPedigrees = listlen(Pedigrees);
  FreeNameList(Pedigrees);


  printf("\nNumber of pedigrees in dataset     : %d\n", nPedigrees);
  printf("Number of individuals in dataset   : %d\n", nObs);
  printf("  Founders                         : %d\n", nFounder);
  printf("  Non-founders                     : %d\n", nNonFounder);
  printf("  Average pedigree size            : %-6.2f\n\n", (double)nObs/nPedigrees);

  printf("Genotyped individuals              : %d\n", nGenotyped);
  printf("  Founders                         : %d\n", nGFounder);
  printf("  Non-founders                     : %d\n\n", nGNonFounder);

  printf("Phenotyped individuals             : %d\n", nPhenotyped);
  printf("  Founders                         : %d\n", nPFounder);
  printf("  Non-founders                     : %d\n\n", nPhenotyped - nPFounder);

  printf("Both geno- and phenotyped          : %d\n", nGenoPheno);
  printf("  Average # gtyped ind pr pedigree : %-6.2f\n", (double)nGenotyped/nPedigrees);
  printf("  Average # ptyped ind pr pedigree : %-6.2f\n", (double)nPhenotyped/nPedigrees);
  printf("  Average # both t ind pr pedigree : %-6.2f\n", (double)nGenoPheno/nPedigrees);

  // Calculate and print all relationships in the dataset
  memset(nRela, 0, REL_MAXREL*sizeof(int));
  memset(nPhenoPair, 0, REL_MAXREL*sizeof(int));

  Rela = MtxNew(nObs, nObs);
  Rela = CalcDataRelationship(indlist);
  for (i=1, ind = indlist; i<nObs; i++, ind = ind->next) {
    for (j=i+1, ind2 = ind->next; j<=nObs; j++, ind2 = ind2->next) {
      nRela[(int)Rela->element[i][j]]++;
      if (IsPhenotyped(ind) && IsPhenotyped(ind2))
	nPhenoPair[(int)Rela->element[i][j]]++;
    }
  }


  printf("\nRelationships in dataset (present and phenotyped pairs):\n");
  for (i=0; i<REL_MAXREL; i++) {
    if (nRela[i]>0) {
      printf("  %-32s : %5d  %5d\n", RelationshipNames[i], nRela[i], nPhenoPair[i]);
    }
  }
  MtxDel(Rela);


 
}


void printok(int ok)
{
  if (ok)
  {
    if (UseColour)
      TextColour(COL_GREEN);
    printf("ok\n");
  }
  else
  {
    if (UseColour)
      TextColour(COL_RED);
    printf("error!\n");
  }
  if (UseColour==2)
    TextColour(COL_CYAN);
  else
    TextColour(COL_WHITE);
}


//
// Locate Mendelian culprit
//

double DetectMendelianCulprit (char *pedname, int markernum, individual *indlist) {
  individual *originalFamily, *genolist, *ind, *ind2;
  markerlist *mkr;
  int foundone = 0;
  int nFamSize = 0;
  int alleleA, alleleB;

  // Start by creating a copy of the list
  originalFamily = SelectPedigree(pedname, indlist);
  nFamSize = listlen(originalFamily);

  for (ind = originalFamily; ind; ind = ind->next) {
    if (istypedformarker(ind, markernum)) {
      // Make a copy and set this person as missing
      mkr = markernumber(ind, markernum);
      alleleA = mkr->allele1;
      alleleB = mkr->allele2;

      // Now set the alleles as missing
      mkr->allele1 = 0;
      mkr->allele2 = 0;

      foundone = 1;
      genolist = GenotypeEliminationFromList(ind->pedigree, markernum, originalFamily);
      for (ind2 = genolist; ind2; ind2 = ind2->next) { 

        // Genotypes not consistent
        if (!listlen(ind2->marker)) {
	  foundone = 0;
	  break;
	}
      }
      FreeIndividualList(genolist);

      if (foundone) {
          sprintf(buf2, "Possible culprit that makes the pedigree (%s) non-Mendelian: %s\n", ind->pedigree, ind->id);
          WriteErrorMsg(buf2);
      }

      mkr->allele1 = alleleA;
      mkr->allele2 = alleleB;


    }
  }


  FreeIndividualList(originalFamily);
  return 0;
}

//
//
// Check consistency of dataset
//
//
// This function goes through the pedigree and checks for possible errors
// using a variety of (simple) checks.
//
void CheckConsistency (individual *indlist)
{
  IDlist *mate, *sibs, *fattree, *mottree, *id, *newid;
  individual *ind, *ind2, *ind3, *genolist;
  int tmpok, code, nofounders, notyped, notypefound, i, j, match, okcheck;
  markerlist *marker;
  int sibmarkers[5], alnum, allelecount, newall;

  if (!indlist)  {
    printf("No data in memory\n");
    return;
  }

  printf("A total of %d individuals from %d pedigree%s to be checked\n\n",
         listlen(individuals), numberofpedigrees(), numberofpedigrees()>1 ? "s":"");

  code = 0;
  nofounders = 0;

  printf("Simple checks:\n");

  // Checking that no individuals with unknown sex exist in dataset
  printf("  Checking for no unknown sexes        : ");
  okcheck = 1;
  for (ind = indlist; ind; ind=ind->next)  {
    if (ind->sex == S_UNKNOWN)  {
      sprintf(buf2, "Individual %s has an unknown sex\n", ind->id);
      WriteErrorMsg(buf2);
      okcheck = 0;
    }
  }
  printok(okcheck);
  
  // Checking that fathers and mothers are the right sex
  printf("  Checking sexes of parents are correct: ");
  okcheck = 1;
  for (ind = indlist; ind; ind=ind->next) {
    // Checking sex of mother and father
    if (ind->father && ind->father->sex != S_MALE) {
      sprintf(buf2, "Individual %s has a father that is not male\n", ind->id);
      WriteErrorMsg(buf2);
      okcheck = 0;
    }
    if (ind->mother && ind->mother->sex != S_FEMALE) {
      sprintf(buf2, "Individual %s has a mother that is not female\n", ind->id);
      WriteErrorMsg(buf2);
      okcheck = 0;
    }
  }
  printok(okcheck);

  // Checking for no mess in IDnumbers
  printf("  Checking for no identical ID's       : ");
  okcheck = 1;

  for (ind = indlist; ind; ind=ind->next) {
    match = 1;
    for (ind2 = ind->next; ind2; ind2 = ind2->next) {
      if (!strcmp(ind->id, ind2->id))
        match++;
    }
    if (match != 1) {
      okcheck = 0;
      sprintf(buf2, "The ID: %s occurs %d times\n", ind->id, match);
      WriteErrorMsg(buf2);
    }
  }
  printok(okcheck);

  printf("  Checking mates for different sexes   : ");
  okcheck = 1;
  // Checking sex of mates
  for (ind = indlist; ind; ind=ind->next) {    
    for (mate = ind->mate; mate; mate = mate->next) {
      if (!(ind->sex-mate->ind->sex)) {
        sprintf(buf2, "%s has a mate of the same sex\n",ind->id);
        WriteErrorMsg(buf2);
        okcheck = 0;
      }
    }
  }
  printok(okcheck);

  // Checking that all persons have both or none of their parents in the dataset
  printf("  Checking both or no parents exist    : ");
  okcheck = 1;

  for (ind = indlist; ind; ind=ind->next) {    
    if (!((!ind->father && !ind->mother) || (ind->father && ind->mother) )) {
      sprintf(buf2, "%s has only one of the two parents in the dataset\n",ind->id);
      WriteErrorMsg(buf2);
      okcheck = 0;
    }
  }
  printok(okcheck);


  // Checking that all persons have both or none of their parents in the dataset
  printf("  Checking markers are not half-typed  : ");

  okcheck = 1;

  for (ind = indlist; ind; ind=ind->next) {  
    i = 1;  
    for (marker = ind->marker; marker ; marker = marker->next) {
      if ((marker->allele1 == 0 && marker->allele2 != 0) || (marker->allele2 == 0 && marker->allele1 != 0)) {
        sprintf(buf2, "%s is not completely typed for marker %d\n",ind->id, i);
        WriteErrorMsg(buf2);
        okcheck = 0;
      }
      i++;
    }
  }
  printok(okcheck);

  // Checks for inbreeding
  // This is done by taking each non-founder individual with at least one parent being a non-founder
  // then expand the parental trees and look for individuals in both fathers and mothers tree
  printf("  Checking for inbreeding              : ");  
  okcheck = 1;
  for (ind = indlist; ind; ind=ind->next)
  {
    if (!founder(ind) && (!founder(ind->father) || !founder(ind->mother)))
    {
      // Make the two trees
      fattree = 0;
      mottree = 0;

      id = cmalloc(sizeof(IDlist));
      memset(id, 0, sizeof(IDlist));

      // Start the paternal tree
      id->ind = ind->father;
      addlist(&fattree, id);
      
      // Expand the tree
      for (id = fattree; id; id = id->next)
      {
	if (!founder(id->ind))
	{
	  newid = cmalloc(sizeof(IDlist));
	  memset(newid, 0, sizeof(IDlist));

	  newid->ind = id->ind->father;
	  addlist(&fattree, newid);

	  newid = cmalloc(sizeof(IDlist));
	  memset(newid, 0, sizeof(IDlist));

	  newid->ind = id->ind->mother;
	  addlist(&fattree, newid);	  
	}
      }


      // Now do the same for the maternal tree

      id = cmalloc(sizeof(IDlist));
      memset(id, 0, sizeof(IDlist));

      // Start the paternal tree
      id->ind = ind->mother;
      addlist(&mottree, id);
      
      // Expand the tree
      for (id = mottree; id; id = id->next)
      {
	if (!founder(id->ind))
	{
	  newid = cmalloc(sizeof(IDlist));
	  memset(newid, 0, sizeof(IDlist));

	  newid->ind = id->ind->father;
	  addlist(&mottree, newid);

	  newid = cmalloc(sizeof(IDlist));
	  memset(newid, 0, sizeof(IDlist));

	  newid->ind = id->ind->mother;
	  addlist(&mottree, newid);	  
	}
      }

      // Have now created the trees
      // Should then check for matches in both

      for (id = fattree; id; id = id->next)
      {
	for (newid = mottree; newid; newid = newid->next)
	{
	  if (id->ind == newid->ind)
	  {
	    // Inbreeding
	    okcheck = 0;
	    sprintf(buf2, "%s is inbred\n", ind->id);
	    WriteErrorMsg(buf2);
	  }
	}
      }
      freelist(fattree);
      fattree = NULL;
      freelist(mottree);
      mottree = NULL;
    }   
  }   
  printok(okcheck);

  // End of inbreeding check


  // Checks for loops
  // This is done by taking each non-founder individual with at least one parent being a non-founder
  // then expand the parental trees and look for individuals in both fathers and mothers tree
  printf("  Checking for marriage loops          : ");  
  /*
  okcheck = 1;
  for (ind = indlist; ind; ind=ind->next)
  {
    if (!founder(ind) && !ind->offspring && (!founder(ind->father) || !founder(ind->mother)))
    {
      // Make the two trees
      fattree = 0;
      mottree = 0;
      spouselist = 0;

      id = cmalloc(sizeof(IDlist));
      memset(id, 0, sizeof(IDlist));

      // Start the paternal tree
      id->ind = ind->father;
      addlist(&fattree, id);
      
      // Expand the tree
      for (id = fattree; id; id = id->next)
      {
	if (!founder(id->ind))
	{
	  newid = cmalloc(sizeof(IDlist));
	  memset(newid, 0, sizeof(IDlist));

	  newid->ind = id->ind->father;
	  addlist(&fattree, newid);

	  newid = cmalloc(sizeof(IDlist));
	  memset(newid, 0, sizeof(IDlist));

	  newid->ind = id->ind->mother;
	  addlist(&fattree, newid);	  
	}
      }


      // Now do the same for the maternal tree

      id = cmalloc(sizeof(IDlist));
      memset(id, 0, sizeof(IDlist));

      // Start the paternal tree
      id->ind = ind->mother;
      addlist(&mottree, id);
      
      // Expand the tree
      for (id = mottree; id; id = id->next)
      {
	if (!founder(id->ind))
	{
	  newid = cmalloc(sizeof(IDlist));
	  memset(newid, 0, sizeof(IDlist));

	  newid->ind = id->ind->father;
	  addlist(&mottree, newid);

	  newid = cmalloc(sizeof(IDlist));
	  memset(newid, 0, sizeof(IDlist));

	  newid->ind = id->ind->mother;
	  addlist(&mottree, newid);	  
	}
      }

      // Have now created the trees
      // Should then check for matches in both

      for (id = fattree; id; id = id->next)
      {
	for (newid = mottree; newid; newid = newid->next)
	{
	  if (id->ind == newid->ind)
	  {
	    // Inbreeding
	    okcheck = 0;
	    sprintf(buf2, "%s is inbred\n", ind->id);
	    WriteErrorMsg(buf2);
	  }
	}
      }
      freelist(fattree);
      freelist(mottree);
    }   
  }   
  printok(okcheck);

  */
  printf("Not implemented yet!\n");

  // End of loop checking


  // Check that males for xlinked chromosomes are homozygous

  if (options->Index[O_XLINKED]) {
    okcheck = 1;
    printf("  Checking for homoz. males (X-linked) : ");  
    for (ind = indlist; ind; ind=ind->next) {
      if (ind->sex==S_MALE) {
	for (marker = ind->marker, i=1; marker ; marker = marker->next, i++) {
	  if (marker->allele1 != marker->allele2) {
	    sprintf(buf2, "%s is not homozygous for marker %d (male for X-linked marker)\n",ind->id, i);
	    WriteErrorMsg(buf2);
	    okcheck = 0;	    
	  }
	}      
      }
    }
    printok(okcheck);
  }





  // Counting founders
  nofounders  = 0;
  notyped     = 0;
  notypefound = 0;
  for (ind = indlist; ind; ind=ind->next)
  {
    if (founder(ind))
      nofounders++;
    if (IsGenotyped(ind))
      notyped++;
    if (founder(ind) && IsGenotyped(ind))
      notypefound++;
  }
  printf ("\nNumber of founders in dataset          : %d\n",nofounders);
  printf ("Number of persons typed in dataset     : %d\n",notyped);
  printf ("Number of genotyped founders in dataset: %d\n\n",notypefound);

  printf("Checking markers:\n");

  // Checking that each full sibship has at most 4 different markers
  printf("  Max 4 alleles in each full sibship   : ");

  // Setting the temp1 string to 0 for all individuals
  // A 1 will then indicate, that the family has more than 4 alleles
  // This is to reduce output (only write once for each pedigree instead of
  // each person)
  for (ind = indlist; ind; ind=ind->next)
    strcpy(ind->tmpstr1, "\0");

  okcheck = 1;
  for (j=1; j<=numberofmarkers(); j++)
  {
    for (ind = indlist; ind; ind=ind->next)
    {
      alnum = 0; 
      marker = markernumber(ind, j);

      if (marker->allele1)
      {
        sibmarkers[alnum] = marker->allele1;
        alnum++;
      }
      if ((marker->allele2) && (marker->allele2 != marker->allele1))
      { 
        sibmarkers[alnum] = marker->allele2;
        alnum++;
      }
      allelecount = alnum;


      for (sibs = ind->sib; sibs; sibs = sibs->next)
      {
        // Only checking full sibs
        if (fullsibs(sibs->ind, ind))
        {
          marker = markernumber(sibs->ind, j);
          newall = 1;
          for (i=0; i<alnum; i++)
          {
            if (marker->allele1 == sibmarkers[i])
              newall = 0;
          }
          if (newall && marker->allele1)
          {
            allelecount++;
            alnum = min(3, alnum);
            sibmarkers[alnum] = marker->allele1;
            alnum++;
          }

          newall = 1;
          for (i=0; i<alnum; i++)
          {
            if (marker->allele2 == sibmarkers[i])
              newall = 0;
          }
          if (newall && marker->allele2)
          {
            allelecount++;
            alnum = min(3, alnum);
            sibmarkers[alnum] = marker->allele2;
            alnum++;
          }
        }        
      } 

    if (allelecount>4)
    {
      strcpy(ind->father->tmpstr1, "1");
      strcpy(ind->mother->tmpstr1, "1");
    }

    if (allelecount>4)
    {
      okcheck = 0;
      sprintf(buf2, "More than 4 alleles among full sibs individual %s (father: %s, mother: %s) for family %s at marker %d\n",
              ind->id, ind->father->id, ind->mother->id, ind->pedigree, j);
      WriteErrorMsg(buf2);
    }
    }  // end individual

    // Now writing to the errorfile if more than 4 alleles occur
    for (ind = indlist; ind; ind=ind->next) {
      if (ind->father) {
        if (strlen(ind->father->tmpstr1)) {

        }
        strcpy(ind->father->tmpstr1, "\0");
      }      
    }

  }
  printok(okcheck);
  printf("\n");


  // Will now check consistency of markers by Genotype Elimination
  for (i = 1; i<=numberofmarkers(); i++)
  {
  printf("  Checking marker %2d consistency       : ", i);  
  okcheck = 1;

  for (ind = indlist; ind; ind=ind->next)
    ind->tmpint1 = 1;

  for (ind = indlist; ind; ind=ind->next)
  {
    if (ind->tmpint1)
    {
      tmpok = 1;
      genolist = GenotypeElimination(ind->pedigree, i);
      for (ind2 = genolist; ind2; ind2 = ind2->next)
      { 

//	printf("Just checking the walls: %s:\n", ind2->id);
//	for (marker = ind2->marker; marker; marker=marker->next) {
//	  printf("   %d  %d\n", marker->allele1, marker->allele2);
//	}


        // Genotypes not consistent
        if (!listlen(ind2->marker))
        {
	  tmpok = 0;
          okcheck = 0;
          sprintf(buf2, "%s is not consistent with the rest of the pedigree for marker %d %s\n", ind2->id, i, options->Index[O_XLINKED] ? "(X-linked)" : "");
          WriteErrorMsg(buf2);
        }
        ind3 = findperson(ind2->id);
        ind3->tmpint1 = 0; 
      }
      FreeIndividualList(genolist);
      if (!tmpok)
	DetectMendelianCulprit (ind->pedigree, i, individuals);
    }
  }
  printok(okcheck);
  }


}


void deletedata (void)
{
  DeletePedigreeData();
  DeleteParameterData();

  freelist(datainfo);  
  datainfo = 0;

}

//
// Sorts the individuals so a persons parents are allways 
// appearing before a person
//
IDlist *SortIndividuals(individual *indlist)
{
  int i, done, code, numpers, famsize;
  individual *ind;
  IDlist *res, *idlist;
  namelist *nl, *pedlist;

  // Start by getting all the pedigree names
  pedlist = 0;
  pedlist = MakePedigreeList(indlist);

  res = 0;

  numpers = listlen(indlist);

  for(ind = indlist; ind; ind = ind->next)
  {
    strcpy(ind->tmpstr1, "0");
  }


  // Sets tmpstr to 0 for all individual
  // 0 Means not set, 1 means set
  for (nl = pedlist; nl; nl = nl->next)
  {
    i = 0;

    // done is a boolean - are we finished for this family?
    done = 0;

    while (!done)
    {
      // code is a boolean - did we place someone this round?
      code = 0;
      famsize = 0;

      // Take a run through the list and order the persons
      for (ind = indlist; ind; ind = ind->next)
      {
	// Check that it is the right pedigree
	if (strcmpl(ind->pedigree, nl->name))
	  continue;

	famsize++;

        // If founder and not previously placed
        if (!ind->father && !ind->mother)
        {
          if (atoi(ind->tmpstr1)==0)
          {
	    idlist = (IDlist *) cmalloc (sizeof (IDlist));
	    memset (idlist, 0, sizeof (IDlist));
	    idlist->ind = ind;
	    addlist(&res, idlist);
            i++;
            code = 1;
	    strcpy(ind->tmpstr1, "1");	    
	  }
        }
        else if (atoip(ind->father->tmpstr1) && atoip(ind->mother->tmpstr1) && (atoip(ind->tmpstr1)==0))
        {
          // Non founder with both parents in data already
	  // Should not be necessary to check for both parents present, but doing it just
	  // in case
	  idlist = (IDlist *) cmalloc (sizeof (IDlist));
	  memset (idlist, 0, sizeof (IDlist));
	  idlist->ind = ind;
	  addlist(&res, idlist);
          i++;
          code = 1;
	  strcpy(ind->tmpstr1, "1");
        }
      }

      if (i==famsize)
      {
        done = 1;
        code = 1;
      }

      // If noone has been placed this turn
      if (!code)
      {
        printf("Error: Could not sort dataset so parents appear before children in pedigree %s\n", nl->name);

	freelist(pedlist);
	pedlist = NULL;
        return res ;
      }
    }
  }
  freelist(pedlist);
  pedlist = NULL;
  return res;
}

// ---------------------------------------------------------

//  OPTION rutines

// ---------------------------------------------------------

void createoptions (void)
{
  // write new optionsfile
  printf ("Writing new options file\n");
  cfopen (OPTIONSFILE, "w");

  fprintf(F, "%-20s ; ID\'s are case sensitive \n",
          optionnames[O_IDCASESENSITIVE]);
  fprintf(F, "%-20s ; Convert ID\'s to FAM-ID when importing linkage files\n",
          optionnames[O_CONVERTLINKAGEID]);
  fprintf(F, "%-20s ; Delete error file when importing new data\n",
          optionnames[O_NEWERRORFILE]);
  fprintf(F, "%-20s ; Reduce number of alleles automatically when importing data\n",
          optionnames[O_REDUCEALLELES]);
  fprintf(F, "%-20s ; Write errors to the error file\n",
          optionnames[O_USEERRORFILE]);
  fprintf(F, "%-20s ; Write iteration info\n",
          optionnames[O_ITERATIONINFO]);
  fprintf(F, "%-20s ; Use ANSI colour\n",
          optionnames[O_COLOUR]);
  fclose (F);
}


//
// Reads the options list
//

void PrepareOptions(optionlist *option)
{
  int i;
  for (i=0; i<MAXOPTIONS; i++) {
    option->Index[i] = 1;
  }  
}

void ReadOptions (void)
{
  int i;
  options = (optionlist *) cmalloc (sizeof (optionlist));
  memset (options,0,sizeof (optionlist));

  if (!(F = fopen (OPTIONSFILE,"r"))) {
    createoptions();
  }
  else {
    getbuf (); 

    while (buf[0] != EOF) {
      // Forget blank lines
      if (!buf[0])
      {
        getbuf ();
        continue;
      }

      i = findoptionskeyword(igetstr(buf));
      if (i>=0 && i<MAXOPTIONS)
	options->Index[i] = 1;

      getbuf ();
    }
    fclose (F);
  }  
  if (options->Index[O_COLOUR] == 1)
    UseColour = 1;
}

void setoptions (void)
{
  int i;

  printf("Options for PediPet v%s\n\n", VERSION);

  if (!(F = fopen (OPTIONSFILE,"r")))
  {
    printf("No options file found - creating new\n");
    createoptions();
    cfopen(OPTIONSFILE, "r");
  }

  for (i = 0; i<MAXOPTIONS; i++)
  {
    printf("%-20s : ",optionnames[i]);
    if (options->Index[O_IDCASESENSITIVE])
       printf("yes\n");
    else
       printf("no\n");
  }
  

  fclose(F);
}

void FounderAlleleFreq (int markernum)
{
  int allarray[MAXALLELES], i,count;
  individual *ind;
  markerlist *marker;

  printf("Calculating allele frequencies based on founders\n");
  printf("Remember to recode alleles before using\n");

  memset(allarray,0, sizeof(allarray));
  count = 0;
  forind
  {
    if (founder(ind))
    {
      marker = markernumber(ind, markernum);
      allarray[marker->allele1]++;
      allarray[marker->allele2]++;
      if (marker->allele1)
        count++;
      if (marker->allele2)
        count++;
    }
  }

  for (i=1; i<=maxallelenumber(markernum); i++)
  {
      printf("%0.4f  ",(double) allarray[i]/count);
  }
  printf("\n");
}

void RemoveAllData(void)
{
  individual *ind;

  forind
  {
    freelist(ind->sib);
    freelist(ind->mate);
    freelist(ind->offspring);
    freelist(ind->marker);
    freelist(ind->qttrait);
  }

  freelist(individuals);
}

//
// Initializes various data
// and start the program
//
void StartUp()
{
  char phrase[81]="";
  unsigned long seed1, seed2;

  // Clear all infor about alleles
  memset(allelefreq, 0, sizeof(allelefreq));

  // Reseting data definitions
  datainfo = (DataDefinitions *) cmalloc (sizeof (DataDefinitions));
  memset (datainfo,0,sizeof (DataDefinitions));

  individuals = 0;
  pedigrees = 0;

  // Reset the errorfile
  if (options->Index[O_NEWERRORFILE])
    remove(ERRORFILE);    

  // Using a combination of time, used id and process id to 
  // generate the seed. Should probably use a high precision 
  // timer instead of the standard time function
  // 
  F = fopen (SEEDFILE, "r");
  if (F>0) {
    getbuf ();
    seed2 = (long) getf();
    seed1 = (long) getf();
    fclose(F);
  }
  else {
    // Use a random starting point if no seed exist
    sprintf(phrase, "%d%ld%d%ld%ld%x", getpid(), random(), (int)getppid(), (long) time(NULL), random(), (int)getuid());
  }

  // If we have the noseed option then start the RNG from the seed file
  if (!options->Index[O_NOSEED]) {
    sprintf(phrase, "%d%ld%d%ld%ld%x", getpid(), seed1, (int)getppid(), (long) time(NULL), seed2, (int)getuid());
    phrtsd(phrase, &seed1, &seed2);
  }
  setall(seed1,seed2);

  allfreq = 0;

  allelefreq[1]=.3;
  allelefreq[2]=.23;
  allelefreq[3]=.23;
  allelefreq[4]=.19;
  allelefreq[5]=.05;

  IntList = 0;

  usedmap = M_KOSAMBI;

  // Model stuff
  ModelTrait = 0;
  nModelCovariates = 0;
}

int ShowFile(char *filename)
{
  char buffer[128];
	
  F = fopen (filename, "r");

  if (F == 0) {
    printf ("Can't show the %s file. Please check your distribution\n",filename);
    return 1;
  }
  else  {
    while (fgets(&buffer[0], 128, F)) {
      printf("%s",buffer);
    }
    
    fclose(F);
  } 
  return 0;
}

void SaveSeedFile()
{
  long seed1, seed2;

  // Save the seed file unless we've asked not to
  if (!options->Index[O_NOSEED]) {
    cfopen (SEEDFILE, "w");
    getsd(&seed1,&seed2);
    fprintf(F, "%ld %ld\n", seed1, seed2);
    fclose (F);
  }
}

void RemoveMarker(individual *indlist, int mkrnum) {
  individual *ind;
  markerlist *mkr;

  for (ind = indlist; ind ; ind = ind->next) {
    mkr = markernumber(ind, mkrnum);
    removelist(&ind->marker, mkr);
  }
}

void KillFounders(int which)
{
  individual *ind, *ind2;
  markerlist *mkr;
  quanttrait *trait;

  printf("Killing off ");
  switch (which) {    
    case S_MALE: printf("male"); break;
    case S_FEMALE: printf("female"); break;
    case -1: printf("all"); break;
  default: printf("ERROR IN REMOVING"); break;
  }
  printf(" founders\n");

  for(ind = individuals; ind; ) {
    ind2 = ind->next;
    if (founder(ind) && (which == -1 || ind->sex==which)) {
      // Sets their alleles to missing
      for (mkr = ind->marker; mkr; mkr = mkr->next)
      {
        mkr->allele1 = 0;
        mkr->allele2 = 0;
      }
      // Sets their trait values to missing
      for (trait = ind->qttrait; trait; trait = trait->next)
      {
	trait->sysmiss = 1;
      }
    }
    ind = ind2;
  }
}


individual *SelectLargestNuclearFamiliesFromPedigree(individual *indlist, char *pedname)
{
  IDlist *sibs;
  individual *ind, *ind2, *ind3, *indkeep;
  int count;
  int maxsibsize = 0;

  indkeep = NULL;

  // Start by finding the largest typed nuclear family
  ind2 = SelectPedigree(pedname, indlist);

  for (ind = ind2; ind; ind=ind->next) {
    // Find the non-founders
    if (!founder(ind)) {
      // Found a larger sibship
      // XXX Potential problem here???
      // If there are large untyped families then what?
      if (listlen(ind->sib)>maxsibsize) {
	// Is the individual genotyped for any of the markers?
	count = IsGenotyped(ind);
	for (sibs = ind->sib; sibs; sibs = sibs->next) {
	  // Now count how many of them are typed for at least one marker
	  count += IsGenotyped(sibs->ind);	    
	}
	if (count>maxsibsize) {
	  maxsibsize = count; 
	  indkeep = ind;
	  // Will keep the pedigre for which indkeep is part of
	}
      }
    }
  }

  
  // No genotyped nuclear families were found
  // return
  if (!indkeep)
    return(NULL);

  // Should now remove everyone but the ones related to indkeep

  for (ind = ind2; ind; ) {
    ind3 = ind->next;

    // Check that not father, mother or individual self
    if (ind != indkeep && ind != indkeep->father && ind != indkeep->mother) {
      // Check that individual is not a sibling
      if (ind->father != indkeep->father || ind->mother != indkeep->mother) {
	RemoveIndividual(ind2, ind);
	// If we remove the first individual then change this
	if (ind2 == ind)
	  ind2 = ind3;
      }
    }    
    ind = ind3;
  }  

  return(ind2);
}

individual *SelectLargestNuclearFamilies(individual *indlist)
{
  namelist *pedlist, *ped;
  individual *ind, *result;

  result = NULL;
  
  pedlist =  MakePedigreeList(indlist);

  for (ped = pedlist; ped; ped=ped->next) {
    ind = SelectLargestNuclearFamiliesFromPedigree(indlist, ped->name);
    if (!ind)
      continue;
    addtolist(&result, ind);
  } 
  freelist(pedlist);

  FreeIndividualList(individuals);
  
  individuals = result; 

  return(result);
}


void deductmean(int traitnum)
{
  quanttrait *qt;
  individual *ind;
  double mean;
  int i, n;

  i = 0;
  n = 0;
  mean = 0.0;
  forind
  {
    for (qt = ind->qttrait; qt && i<traitnum; qt = qt ->next) {
      i++;
      if (i==traitnum && !qt->sysmiss) {
        n++;
        mean += qt->value;
      }
    }
  }
  mean /= n;

  forind
  {
    for (qt = ind->qttrait; qt && i<traitnum; qt = qt ->next)
    {
      i++;
      if (i==traitnum && !qt->sysmiss)
      {        
        qt->value -= mean;
      }
    }
  }
}


void takelog(int traitnum)
{
  quanttrait *qt;
  individual *ind;
  int i, j;

  j=0;
  forind
  {
    i=0;
    for (qt = ind->qttrait; qt && i<traitnum; qt = qt ->next)
    {
      i++;
      if (i==traitnum && !qt->sysmiss)
      {
        if (qt->value>0)
          qt->value =  log(qt->value);
        else
        {
          qt->value = 0;
          qt->sysmiss = 1;
          j++;
        }
      }
    }
  }
  printf("Logarithm taken on trait %d.\n", traitnum);
  if (j>0)
    printf("WARNING: %d values were negative and set to missing.\n", j);

  //  deductmean(traitnum);
  //  printf("Mean deducted\n");
}

void SetMissing(int traitnum)
{
  quanttrait *qt, *qt2;
  individual *ind;
  int j;

  j=0;
  forind
  {
    qt = traitvalue(ind,traitnum);
    // If value 1 is diabetic
    qt2 = traitvalue(ind,1);

    if (qt2->value > 2)
    {
      qt->value = 0;
      qt->sysmiss = 1;
    }
  }
  printf("All non-1 values sat to missing trait %d.\n", traitnum);
}


void GetNumberList(char *streng)
{
  int number, i, code;
  char *string;

  number = 0;

  IntList = (IntegerList *) cmalloc (sizeof (IntegerList));
  memset (IntList, 0, sizeof (IntegerList));

  IntList->length = 0;

  // Count the numbers on in the string
  if (strlen(string = igetstr(streng)))
  {
    if (atoi(string)>=0)
    {
      number = 1;

      while ( strlen(string = getstr()) )
      {
        if (atoi(string)>=0)
          number++;
      }
    }
  }

  // Sets the model
  nModelCovariates = number;

  // Fill up a vector
  if (number>0)
  {
    IntList->length = number;   
    IntList->values = ivector(1,number);
    i = 1;

    if (strlen(string = igetstr(streng)))
    {
      if ((code = atoi(string))>=0)
      {
        IntList->values[i] = code; 
        i++;
      }

      while ( strlen(string = getstr()) )
      {
        if ((code = atoi(string))>=0)
        {
          IntList->values[i] = code; 
          i++;
        }
      } 
    }
  }

  // Sets the model parameters
  for (i=0; i<number; i++)
    ModelCovariates[i] = IntList->values[i+1];
}

int findkeyword (char *s)
{
  char **sp;

  // Special shortcuts
  if (!strncasecmp (s,"q",1))
    return C_QUIT;
  if (!strncasecmp (s,"?",1))
    return C_HELP;

  sp = (char **) bsearch (&s, keywords, MAXKEYWORDS, sizeof s, keywordcmp);
  if (sp == 0)
  {
    return -1;
  }

  return sp - keywords;
}


// Make a copy of an individual list
/**********************
 *
 * Creates an identical copy of an individual list.
 *
 *
 **********************/
individual *CopyIndividualList(individual *indlist)
{
  individual *ind, *ind2, *minilist;
  markerlist *mkr, *marker;
  quanttrait *qt, *newqt;

  // First make a copy of the examined pedigree
  minilist = 0;

  for (ind = indlist; ind; ind = ind->next)
  {
    // Is the individual from the correct pedigree?
    // if (!strcmpl(pedigree, ind->pedigree))
    {
      ind2 = cmalloc (sizeof (individual));
      memset(ind2, 0, sizeof(individual));

      strcpy(ind2->id, ind->id);
      strcpy(ind2->pedigree, ind->pedigree);
      ind2->sex = ind->sex;
      ind2->tmpint1 = ind->tmpint1;
      ind2->tmpint2 = ind->tmpint2;
      ind2->tmpint3 = ind->tmpint3;
      ind2->tmpint4 = ind->tmpint4;

      if (ind->father)
        strcpy(ind2->tmpstr1, ind->father->id);
      else
        strcpy(ind2->tmpstr1,"\0");
      if (ind->mother)
        strcpy(ind2->tmpstr2, ind->mother->id);
      else
        strcpy(ind2->tmpstr2,"\0");

      addlist(&minilist, ind2);

      // Copies the marker data
      for (mkr = ind->marker; mkr; mkr  = mkr->next) {
        marker = cmalloc (sizeof (markerlist));
        memset (marker,0,sizeof (markerlist));

        marker->allele1 = mkr->allele1;
        marker->allele2 = mkr->allele2;

        addlist(&ind2->marker, marker);
      }

      // Copies the phenotypes
      for (qt = ind->qttrait; qt; qt = qt->next) {
	newqt = cmalloc (sizeof (quanttrait));
	memset(newqt, 0, sizeof(quanttrait));

	newqt->value = qt->value;
    	newqt->sysmiss = qt->sysmiss;

	addlist(&ind2->qttrait, newqt);
      }
    }
  }

  // Should now fix fathers, mothers and sibs
  for (ind = minilist; ind; ind = ind->next)
  {
    for (ind2 = minilist; ind2; ind2 = ind2->next)
    {
      if (!strcmpl(ind2->id, ind->tmpstr1))
        ind->father = ind2;
      if (!strcmpl(ind2->id, ind->tmpstr2))
        ind->mother = ind2;
    }
    if (ind->father)
      adduniqueidlist (&ind->father->offspring,ind);
    if (ind->mother)
      adduniqueidlist (&ind->mother->offspring,ind);
  }

  for (ind = minilist; ind; ind = ind->next)
  {

    if (!founder(ind))
    {

      // Adds spouses for parents
      adduniqueidlist(&ind->father->mate, ind->mother);
      adduniqueidlist(&ind->mother->mate, ind->father);


      // Fixing sibs
      for (ind2 = ind->next; ind2; ind2 = ind2->next)
      {
        if (((ind->father == ind2->father) && (ind->father)) ||
            ((ind->mother == ind2->mother) && (ind->mother)))
        {
          adduniqueidlist(&ind->sib, ind2);
          adduniqueidlist(&ind2->sib, ind);
	}
      }
    }
  }
  // The pedigree structures should be correct by now. 
  return minilist;
}



// Make a copy of an individual list for a specific pedigree
individual *CopyIndividualListPedigree(individual *indlist, char *pedigree)
{
  individual *ind, *ind2, *minilist;
  markerlist *mkr, *marker;
  quanttrait *qt, *newqt;

  // First make a copy of the examined pedigree
  minilist = 0;

  for (ind = indlist; ind; ind = ind->next)  {
    // Is the individual from the correct pedigree?
    if (!strcmpl(pedigree, ind->pedigree))
    {
      ind2 = cmalloc (sizeof (individual));
      memset(ind2, 0, sizeof(individual));

      strcpy(ind2->id, ind->id);
      strcpy(ind2->pedigree, ind->pedigree);
      ind2->sex = ind->sex;
      ind2->tmpint1 = ind->tmpint1;
      ind2->tmpint2 = ind->tmpint2;
      ind2->tmpint3 = ind->tmpint3;
      ind2->tmpint4 = ind->tmpint4;

      if (ind->father)
        strcpy(ind2->tmpstr1, ind->father->id);
      else
        strcpy(ind2->tmpstr1,"\0");
      if (ind->mother)
        strcpy(ind2->tmpstr2, ind->mother->id);
      else
        strcpy(ind2->tmpstr2,"\0");

      addlist(&minilist, ind2);

      // Copies the marker data
      for (mkr = ind->marker; mkr; mkr  = mkr->next) {
        marker = cmalloc (sizeof (markerlist));
        memset (marker,0,sizeof (markerlist));

        marker->allele1 = mkr->allele1;
        marker->allele2 = mkr->allele2;

        addlist(&ind2->marker, marker);
      }

      // Copies the phenotypes
      for (qt = ind->qttrait; qt; qt = qt->next) {
	newqt = cmalloc (sizeof (quanttrait));
	memset(newqt, 0, sizeof(quanttrait));

	newqt->value = qt->value;
    	newqt->sysmiss = qt->sysmiss;

	addlist(&ind2->qttrait, newqt);
      }
    }
  }

  // Should now fix fathers, mothers and sibs
  for (ind = minilist; ind; ind = ind->next)
  {
    for (ind2 = minilist; ind2; ind2 = ind2->next)
    {
      if (!strcmpl(ind2->id, ind->tmpstr1))
        ind->father = ind2;
      if (!strcmpl(ind2->id, ind->tmpstr2))
        ind->mother = ind2;
    }
    if (ind->father)
      adduniqueidlist (&ind->father->offspring,ind);
    if (ind->mother)
      adduniqueidlist (&ind->mother->offspring,ind);
  }

  for (ind = minilist; ind; ind = ind->next)
  {

    if (!founder(ind))
    {

      // Adds spouses for parents
      adduniqueidlist(&ind->father->mate, ind->mother);
      adduniqueidlist(&ind->mother->mate, ind->father);


      // Fixing sibs
      for (ind2 = ind->next; ind2; ind2 = ind2->next)
      {
        if (((ind->father == ind2->father) && (ind->father)) ||
            ((ind->mother == ind2->mother) && (ind->mother)))
        {
          adduniqueidlist(&ind->sib, ind2);
          adduniqueidlist(&ind2->sib, ind);
	}
      }
    }
  }
  // The pedigree structures should be correct by now. 
  return minilist;
}


// Selects a specific pedigree from indlist, and makes a copy of the list
individual *SelectPedigree(char *pedname, individual *indlist)
{
  individual *res = NULL;

  res = CopyIndividualListPedigree(indlist, pedname);

  return res;
}

//
// Returns the number of individuals in indlist belonging to pedigree pedid
//
int FamilySize(individual *indlist, char *pedid)
{
  individual *ind;
  int res;

  res = 0;
  for (ind = indlist; ind; ind = ind->next)
  {
    if (!strcmpl(pedid, ind->pedigree))
      res++;
  }
  return res;
}

void GXEanalyse() 
{
  int code2, code1;
  char buf2[100];

  //  ExportLinkagePedigree(individuals,1); 
  printf("Which trait : ");
  code1 = atoi(InputLine(buf2, BUFFERSIZE));
  printf("Which env   : ");
  code2 = atoi(InputLine(buf2, BUFFERSIZE));

  printf("Input covariates : ");
  InputLine(buf2, 99);
  GetNumberList(buf2);
  
  MultiPointGXE(individuals, code1, code2, 0, IntList->length, IntList->values);

}


void ClausStuff()
{
  /*
  individual *ind, *ind2;
  markerlist *mkr;
  double  res;
  int i, code2;
  char buf2[100];
  */

  ThreepointLinkageAnalysis(individuals, order[1], order[2], order[3], .1, .1, .1) ;

  /*
  //  ExportLinkagePedigree(individuals,1); 
  printf("Which trait : ");
  code2 = atoi(InputLine(buf2, BUFFERSIZE));

  printf("Input covariates : ");
  InputLine(buf2, 99);
  GetNumberList(buf2);
  
  MultiPointVC(individuals, code2, 0, IntList->length, IntList->values);
  */


  /*


  //  ExportLinkagePedigree(individuals,1); 
  printf("Which trait : ");
  code2 = atoi(InputLine(buf2, BUFFERSIZE));

  printf("Input covariates : ");
  InputLine(buf2, 99);
  GetNumberList(buf2);
  
  MultiPointXVC(individuals, code2, 0, IntList->length, IntList->values);

  // make data for genehunter  
  ExportLinkagePedigree(individuals,1); 
  MakeGHdatafile ();

  system("gh < rungh > /dev/null");
  system("reduce-ibd");

  Paper1Function2(individuals, 1);
  
  //  TestMatrix(individuals);
  
  //  ReadSolarIBDFile(individuals, "mibd.22.4", "pedindex.out");


  for (ind = individuals; ind; ind = ind->next) {
    for (ind2 = ind; ind2; ind2 = ind2->next) {
      printf("%s <-> %s   : %5.3f\n", ind->id, ind2->id, xlinkedphi2(ind, ind2));
    }
  }

  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .15));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .20));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .25));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .30));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .35));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .40));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .45));

  */

  /*
  res = TwopointLinkageAnalysis(individuals, 1, 2,.45) ;
  printf("Theta = %f\n", res);

  //  printf("LLH = %f\n", ThreepointLinkage(individuals, 1, 2, 3, .09, .05, .01));
  //  printf("LLH = %f\n", ThreepointLinkage(individuals, 1, 2, 3, .10, .05, .01));
  //  printf("LLH = %f\n", ThreepointLinkage(individuals, 1, 2, 3, .11, .05, .01));
  //  printf("LLH = %f\n", ThreepointLinkage(individuals, 1, 2, 3, .12, .05, .01));

  ThreepointLinkageAnalysis(individuals, 1, 2, 3, .1, .1, .1) ;


  printf("Start of ClausStuff\n");
  i = 1;

  printf("Keeping only genotypes from individuals with trait3==%d\n",i);
  forind {
    if (trait(ind,2) == i || trait(ind,2) > 3) {
      for (mkr = ind->marker; mkr; mkr = mkr->next) {
	mkr->allele1 = 0;
	mkr->allele2 = 0;
      }      
    }
  }
  */

  //
  //
  //
  //

  /*
  for (i = 0; i< 1024; i++)
    if (IBDforInheritanceVector(individuals->pedigree, i)) 
      printf("DIVINE: %d\n", i);

  */

  printf("End of ClausStuff\n");

}

void TestMapExpander() 
{
  //  double res;

  /*
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .15));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .20));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .25));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .30));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .35));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .40));
  printf("LINKAGE: %f\n", TwopointLinkage(individuals, 1, 2, .45));

  res = TwopointLinkageAnalysis(individuals, 1, 2,.45) ;
  printf("Theta = %f\n", res);

  */

  //  printf("LLH = %f\n", ThreepointLinkage(individuals, 1, 2, 3, .09, .05, .01));
  //  printf("LLH = %f\n", ThreepointLinkage(individuals, 1, 2, 3, .10, .05, .01));
  //  printf("LLH = %f\n", ThreepointLinkage(individuals, 1, 2, 3, .11, .05, .01));
  //  printf("LLH = %f\n", ThreepointLinkage(individuals, 1, 2, 3, .12, .05, .01));


  ThreepointLinkageAnalysis(individuals, order[3], order[4], order[5], .1, .1, .1) ;
//  ThreepointLinkageAnalysis(individuals, order[1], order[2], order[3], .1, .1, .1) ;

}

void TestStuff()
{
  //  int markernum = 1;
  

  CalcPIHat(individuals, 4);

  /*
  // make data for genehunter  
  ExportLinkagePedigree(individuals,1); 
  MakeGHdatafile ();

  printf("Running GeneHunter ...");
  system("rm -f ibd_dist.out");
  system("gh < rungh > /dev/null");
  system("convert-ibd");
  printf(" done\n");

  // list, trait, marker
  SibsComplexMixed(individuals, 1, markernum);
  printf("HEJ\n");

  GHPiHat(individuals, 1, markernum);
  printf("End of sibs\n");

  //  ComplexMixed(individuals, 1);
  */

  printf("End of TestStuff\n");
}

void MakeSolarIBDFile(int marker)
{
  individual *ind, *ind2;

  sprintf(buf, "ibd.Marker%d", marker);
  cfopen (buf, "w");

  for (ind = individuals; ind; ind = ind->next) {
    for (ind2 = individuals; ind2 != ind->next && ind2->next; ind2 = ind2->next) {
      if (SamePedigree(ind,ind2)) {
	if (ind == ind2)
	  fprintf(F, "%5d %5d  %9.7f  %9.7f\n", 
                ind->globalid, ind->globalid,1.0, 1.0); 
        else
	  {
	    if (!founder(ind) && founder(ind2))
	      fprintf(F, "%5d %5d  %9.7f  %9.7f\n", 
                ind->globalid, ind2->globalid, 0.5, 0.0);       
	    if ((!founder(ind) && !founder(ind2)))
          fprintf(F, "%5d %5d  %9.7f  %9.7f\n", 
              ind->globalid, ind2->globalid,SingleMarkerIBD(ind,ind2,marker), k2);
	  }
      }
    }
  }

  fclose(F);
  sprintf(buf, "gzip ibd.Marker%d", marker);
  system(buf);
}





void ReadCommandLineOptions(int argc, char *argv[])
{
  int i;

  for (i=1; i<argc; i++)
  {
    if (!strcmpl(argv[i], "--colour"))
      UseColour = 1;
    if (!strcmpl(argv[i], "--color"))
      UseColour = 1;

    if (!strcmpl(argv[i], "--heavy"))
      UseColour = 2;

    if (!strcmpl(argv[i], "--noseed"))
      options->Index[O_NOSEED] = 1;

    printf("LKNLKJ\n");

    if (strspn(argv[i], "<>|"))
      break;
  }
}

//
// Used for command input
//

char *InputLine(char *buffer, int bufsize)
{
  int i;

  strcpy(buffer, "\0");

  fgets(buffer, bufsize, stdin);

  // Remove trailing spaces and CR / LF 
  if (buffer && ((i = strlen(buffer)) != 0))
  {
    while (--i>=0)
    {     
      if (!isspace(buffer[i]))
	break;
    }    
    buffer[++i] = '\0';
  }

  return buffer;
}


int main (int argc, char *argv[])
{
  int code, mincode, maxcode, code2, i, j;
  //  int *covariates;
  char *string;
  double res1;
  double a, b, c, d, e;
  long seed1, seed2;
  individual *ind, *ind2;
  markerlist *marker;
  OLDMATRIX *ibd, *ibd2;

  int nModelTraits, *ModelTraits;

  code2 = 0;

  // Flush Output Buffer (force unbuffered output)
  setbuf(stdout, (char *) 0);

  // Initializing
  nModelTraits = 0;    // No traits selected
  ModelTraits  = NULL;  // No traits selected

  UseColour = 0;

  // Read the options file
  ReadOptions();

  if (argc>1)
    ReadCommandLineOptions(argc, argv);

  StartUp();


  BackgroundColour(COL_BLACK);
  TextColour(COL_BLUE);
  puts ("PediPet v" VERSION "  (" __DATE__ ")\n"
        "Copyright (C) 1997-2006 by Claus Ekstroem.\n\n"
        "PediPet is free software and comes with ABSOLUTELY NO WARRANTY.\n"
	"You are welcome to redistribute it under certain conditions;\n"
	"type `show license' for details.\n" 

        "Type ? for list of commands.\n");
  TextColour(COL_WHITE);

  for (;;)
  {
    if (UseColour==2) TextColour(COL_CYAN);

    if (!feof(stdin)) {
      printf ("\nInput command > ");
      InputLine(buf, BUFFERSIZE);
    }
    else  // Set the order to quit if EOF
      strcpy(buf, "quit\0");
    if (UseColour==2) TextColour(COL_WHITE);

    switch (findkeyword(igetstr(buf)))
    {
      case C_QUIT   : fflush(0);
                      RemoveAllData();  // Remove data and
                      TextColour(COL_WHITE);
                      BackgroundColour(COL_BLACK);
		      SaveSeedFile();
		      sync();           // Sync buffers/files
                      exit(0);          // quit program
                      break;
      case C_RECODE : inputredall ();   // Recode alleles for all markers
                      break;
      case C_DELETE : deletedata ();    // Deletes all data
                      break;
      case C_CONSISTENCY: CheckConsistency (individuals); // Checks dataset for consistencies
                      break;
      case C_KILL   : switch (findkeyword(getstr()))
	              {
		      case C_FATHER: KillFounders(S_MALE);
			break;
		      case C_MOTHER: KillFounders(S_FEMALE);
			break;
		      default:
			KillFounders(-1); break;
		      }
                      break;

      case C_LOAD   : ReadLinkageParameterFile(); break;
      case C_IMPORT : ReadDataset ();
                      break;
      case C_CLONE  : code = atoi(getstr());
	if (code>0) {
	  ind = CloneDataset(individuals, code);
	  DeleteIndividualList(individuals);
	  individuals = ind;
	}
	else {
	  printf("Must specify a (positive) number of cloning to perform.\n");
	}
	break;
      case C_ELLEN :  //for (i = 0; i<(int) (100*ChromosomeLength()); i++) {
                      //   sprintf(minibuf, "mibd-%d.pp", i);
                      //   printf("Analyzing position %d\n", i);
	                AnalyzeAllMatrices(individuals, "mibd-0.pp", 1);
		      // }
	              break;
      case C_MIXTURE : AnalyzeMixtureModels(individuals, "marker4.ibd", 1);
                       AnalyzeMixtureModelsNULL(individuals, "marker4.ibd", 1);
	              break;
      case C_SIMULATE:switch (findkeyword(getstr()))
                      {
                        case C_CEPH : printf("CEPH families not implemented\n"); break;
		      case C_CLAUS: GenotypeError(individuals->next->next, 1, 1); break;
                        case C_SIBS : // Gets number
			              code = atoi(getstr());
				      if (code) {
					a = sqrt(2.0/1.0);
//					a = sqrt(16.0/(9*9.0));
//					a = sqrt(4.0/(3*9.0));
					individuals = SimulateChromosome(code, 2, 30, a, 0*a, .5, 0, 0); 
				      }
			              else
                                        printf("Must specify at least one family\n");
				      break;
                      case C_SPECIAL : code = atoi(getstr());
                                       if (code)
                                           individuals = SimulateChromosome(code, 4, 30, 0*sqrt(2.0), 0*sqrt(.5), .5, 0, 0);
			                 else
                                           printf("Must specify at least one family\n");
                                         break;
		      case C_ECHO : code = atoi(getstr());
                                    if (code)  // Kosambi mapping function
                                           individuals = SimulateKosambiChromosome(code, 2, 80, sqrt(2.0), 0*sqrt(.5), .5, 0, 0);
			                 else
                                           printf("Must specify at least one family\n");
                                         break;
		      case C_TEST : code = atoi(getstr());
			            if (code) {
			              a = atof(getstr());
			              b = atof(getstr());
			              individuals = SimulatePedigreeStructure(code, 0, 2, 1);
				      SimulateDataset(individuals, 30, a, 0, .5, b, 0);
		                    }
 		                    else
				      printf("Must specify at least one family\n");
		                    break;
		      case C_ERROR : string = getstr();
                                     if (!strcmpl(string, "all")) { 
                                       for (i=1; i<=numberofmarkers(); i++)
					 RandomGenotypeError(individuals, i);
				     }
                                     else
				     {
				       code=atoip(string);
				       if (code)
					 RandomGenotypeError(individuals, code);
				       else
					 puts("Illegal marker number");
				     }
			             break;
       		      case C_MISSING : string = getstr();
                                     if (!strcmpl(string, "all"))
				     { 
				       res1 = atof(getstr());

                                       for (i=1; i<=numberofmarkers(); i++)
					 SimulateRandomMissingGenotypes(individuals, i, res1);
				     }
                                     else
				     {
				       code=atoip(string);
				       if (code)
					 SimulateRandomMissingGenotypes(individuals, code, atof(getstr()));
				       else
					 puts("Illegal marker number");
				     }
			             break;
		      case C_DATA:   a = 0; b = 0, c=30, d=.5, e=0;
			             a = atof(getstr());
			             b = atof(getstr());
			             c = atof(getstr());
			             d = atof(getstr());
			             e = atof(getstr());
				     if (c == 0)
				       c = 30;
				     if (d == 0)
				       d = .5;
				     // Call the correct function if data are X-linked
				     if (options->Index[O_XLINKED]) {
				       SimulateXDataset(individuals, c, a, 0, d, b, e);
				     }
				     else {
				       SimulateDataset(individuals, c, a, 0, d, b, e);
				     }

			             break;
		      case C_STENO:    a = 0; b = 0; b = atof(getstr());
			                 a = atof(getstr());
					 SimulateStenoDataset(individuals, 30, b, 0, .2, a, 0, 1);
			             break;
		      case C_MARKER: if ((code = atoi(getstr())) && (code2 = atoi(getstr())) && (a = atof(getstr())) && (b = atof(getstr()))) {
                                         SimulateMarkerSet(individuals, code, code2, a, b);
		                     } 
			             break;
		        case C_GENOTYPES: SimulateGenotypes(individuals);
			             break;
		      default: printf("Simulate what??\n"); break;
				       
                      }
                      break;

      case C_ARTIKEL1: Paper1Function(individuals, 1);
                      break;
      case C_REMOVE:  printf("Removes founders genotypes\n");
                      code = numberofmarkers();
                      forind
		      {
                        if (founder(ind))
			{
                          for (marker = ind->marker; marker; marker = marker->next)
			  {
			    marker->allele1 = 0;
			    marker->allele2 = 0;
			  }
			}
		      }
                      break;

    case C_OPTION:  i = findoptionskeyword(getstr());
                    if (i>=0 && i < MAXOPTIONS)
		      options->Index[i] = 1 - options->Index[i];

                      break;
    case C_EXPAND : TestMapExpander();
	              break;

      case C_FIX    : SetMissing(2);
	              break;
      case C_FREQ   : switch (findkeyword(getstr()))
                      {
                        case C_EQUAL : string = getstr();
                                       if (!strcmpl(string, "all")) { 
                                         for (i=1; i<=numberofmarkers(); i++)
 			                   SetEqualFreq(i);                    
				       }
                                       else {
					 if (atoip(string)>0 && atoip(string)<=numberofmarkers())
					   SetEqualFreq(atoip(string)); 
				       }
                                       break;
                        case C_COUNT : string = getstr();
                                       if (!strcmpl(string, "all"))
				       { 
                                         for (i=1; i<=numberofmarkers(); i++)
 			                   SetCountFreq(i);
				       }
                                       else
				       {
                                         if (atoip(string))
   			                   SetCountFreq(atoip(string)); 
                                         else
                                         puts("Syntax error: Allele frequency for which marker?");
                                       }
                                       break;
		        case C_ML :    string = getstr();
			               if (atoip(string))
				       {
					 MLAlleleFrequencyEstimation(individuals, atoip(string));
				       }
                                       else
                                         puts("Syntax error: Estimate allele frequencies for which marker?");
				       break;
			               
                        case C_INPUT : SetInputFreq(atoip(getstr())); 
                                       break;
                        default : puts ("Syntax error: Set allele frequencies as what?"); break;
                      }
                      break;
      case C_DISTANCE: if (numberofmarkers()<2)
     	               {
                         printf("Can not input distances with %d marker in dataset\n", numberofmarkers());
                         break;
		       }
                       printf("Input the %d distances in cM (e.g. 10 5.6 3.4 16.9):\n",numberofmarkers()-1);
		       InputLine(buf, BUFFERSIZE);
                       res1 = atof(igetstr(buf));
                       for (i=1; i<numberofmarkers(); i++)
		       {
                         if (res1 == 0.0)
			 {
                           printf("  Distance %d set to 10 cM\n",i);
                           distance[i] = InverseMapFunction(.1);
			 }
			 else
                           distance[i] = InverseMapFunction(res1/100);
                         res1 = getf(); 
		       }
                      break;
      case C_MAKE   : switch (findkeyword(getstr()))
                      {
                        case C_PARAM : ExportLinkageParameter(individuals);
                                       break;
                        case C_SOLAR : ExportSolarPedigree(individuals);
			               ExportSolarMarker(individuals);
				       ExportSolarMap();
				       ExportSolarPhenotype(individuals);
			               break;
                        case C_CRIMAP: ExportCRIMAPGen (individuals, options->Index[O_XLINKED]);
			               ExportCRIMAPPar (individuals);
                                       break;
                        case C_PEDIGREE: ExportLinkagePedigree(individuals,0); 
                                       break;
                        case C_MENDEL: ExportMendelPedigree(individuals);
			               ExportMendelLocus(individuals);
				       ExportMendelMap();
                                       break;
                        case C_RELPAIR: ExportRelpairControl();
			               printf("Which chromosome is this: ");
			               ExportRelpairLocus(atoi(InputLine(buf2, BUFFERSIZE)));
			               ExportRelpairPedigree(individuals);
                                       break;
                        case C_GH: ExportLinkagePedigree(individuals,1); MakeGHdatafile (); break;
                        case C_MERLIN: ExportMerlinPedigree(individuals); 
			               ExportMerlinDataFile(); 
				       ExportMerlinMap();
				       ExportMerlinFreq();
				       break;
                        default : puts ("Syntax error: Make what?"); break;
                      }
                      break;
      case C_HASEMAN: code  = atoi(getstr());
                      printf("Which trait : ");
                      code2 = atoi(InputLine(buf2, BUFFERSIZE));

                      printf("Input covariates : ");
		      InputLine(buf2, BUFFERSIZE);
                      GetNumberList(buf2);

                      HasemanElston (1, code, code2,
                      IntList->values, IntList->length);

                      // Frees data used for covariates
                      if (IntList->length>0)
                        free_ivector(IntList->values, 1, IntList->length);
                      free(IntList);
                      break;
      case C_IBD    : code = atoip(getstr());
	              if (code>0 && code <= numberofmarkers()) {
			mincode = code;
			maxcode = code;
		      }
		      else {
			mincode = 1;
			maxcode = numberofmarkers();			
		      }

		      for (code=mincode; code <= maxcode; code++) {
                        ibd = CalcPIHat(individuals, code);		      
                        sprintf(buf, "marker%d.ibd",code);
#ifndef BYPEDIGREE			
 		        MtxFSymmetricPrint(buf, ibd);
			
			if (options->Index[O_XLINKED]) {
			  // Make a version of the IBD matrix that depends 
			  // on the sex of the individuals
			  // Assuming below, that ibd is calculated 
			  // for all individuals

			  // Male-male
			  ibd2 = MtxNew(ibd->rows,ibd->cols);
			  MtxCopy(ibd, ibd2);

			  for (ind=individuals, i=1; ind; ind=ind->next, i++) {
			    for (ind2=ind, j=i; ind2; ind2=ind2->next, j++) {
			      if (ind->sex!=S_MALE || ind2->sex!=S_MALE) {
				ibd2->element[i][j] = 0;
				ibd2->element[j][i] = ibd2->element[i][j];
			      }
			    }			    
			  }
			  sprintf(buf, "marker%d-mm.ibd",code);
			  MtxFSymmetricPrint(buf, ibd2);


			  // Female-female
			  MtxCopy(ibd, ibd2);

			  for (ind=individuals, i=1; ind; ind=ind->next, i++) {
			    for (ind2=ind, j=i; ind2; ind2=ind2->next, j++) {
			      if (ind->sex!=S_FEMALE || ind2->sex!=S_FEMALE) {
				ibd2->element[i][j] = 0;
				ibd2->element[j][i] = ibd2->element[i][j];
			      }
			    }			    
			  }
			  sprintf(buf, "marker%d-ff.ibd",code);
			  MtxFSymmetricPrint(buf, ibd2);

			  // Male-female
			  MtxCopy(ibd, ibd2);
			  for (ind=individuals, i=1; ind; ind=ind->next, i++) {
			    for (ind2=ind, j=i; ind2; ind2=ind2->next, j++) {
			      if (!((ind->sex==S_MALE && ind2->sex==S_FEMALE) || (ind->sex==S_FEMALE && ind2->sex==S_MALE))) {
				ibd2->element[i][j] = 0;
				ibd2->element[j][i] = ibd2->element[i][j];
			      }
			    }			    
			  }
			  sprintf(buf, "marker%d-mf.ibd",code);
			  MtxFSymmetricPrint(buf, ibd2);


			  MtxDel(ibd2);
			}

                        MtxDel(ibd);
#endif
		      }
                      break;
      case C_MIBD   : // Starts by checking that all ibd-filer are there
		      for (code=1; code <= numberofmarkers(); code++) {

		      }
                      // Calculate and save the relationsship-matrix
		      ibd = CalcDataRelationship(individuals);
		      MtxFSymmetricPrint("relation.pp", ibd);
		      MtxDel(ibd);
		      FullSibsIBD(individuals, "fulkerpi.out", 1);

	              break;
      case C_VC     : code = atoip(getstr());
                      printf("Which trait : ");
                      code2 = atoi(InputLine(buf2,BUFFERSIZE));

                      printf("Input covariates : ");
                      InputLine(buf2, BUFFERSIZE);
                      GetNumberList(buf2);

		      if (options->Index[O_XLINKED]) {  
			SinglePointXVC(individuals, code2, code, IntList->length, IntList->values);
		      }
		      else {
			SinglePointVC(individuals, code2, code, IntList->length, IntList->values);
		      }
                      // Frees data used for covariates
                      if (IntList->length>0)
                        free_ivector(IntList->values, 1, IntList->length);
                      free(IntList);
                      break;
      case C_SCAN   : // Do a complete scan using a specified method
                      switch (findkeyword(getstr()))
                      {
                        case C_HASEMAN : printf("Which trait : ");
                                         code = atoi(InputLine(buf2,BUFFERSIZE));

                                         printf("Input covariates : ");
                                         InputLine(buf2, BUFFERSIZE);
                                         GetNumberList(buf2);

                                         for (i=1; i<=numberofmarkers(); i++)
                                         {
  			                   HasemanElston (1, order[i], code,
                                           IntList->values, IntList->length);
                                           printf("\n");
                                         }
                                         // Frees data used for covariates
                                         if (IntList->length>0)
                                           free_ivector(IntList->values, 1, IntList->length);
                                         free(IntList);

                                         break;
                        case C_VC : printf("Which trait : ");
                                    code = atoi(InputLine(buf2,BUFFERSIZE));

                                    printf("Input covariates : ");
                                    InputLine(buf2, BUFFERSIZE);
                                    GetNumberList(buf2);

                                    for (i=1; i<=numberofmarkers(); i++)
                                    {                    
				      SinglePointVC(individuals, code2, i, IntList->length, IntList->values);
                                    }
                                    // Frees data used for covariates
                                    if (IntList->length>0)
                                      free_ivector(IntList->values, 1, IntList->length);
                                    free(IntList);

                                    break;
                        default : puts("Syntax Error: Use what scanning method?"); break;
                      }
                      // End of scan
                      break;
      case C_ECHO    : string = strstr(buf,igetstr(buf));
                       string +=4; // Skip the "ech" part of the string      
                       if (*string == ' ')
                         string++;
                       printf("%s",string);
                       break;
    case C_SELECT : SelectLargestNuclearFamilies(individuals);
                       break;
      case C_IBS     : if (!strcmpl(getstr(), "all")) {
                         VerifyIndependence(individuals, 1);
                       } 
                       else {
                         VerifyIndependence(individuals, 0);
                       }
                       break;
      case C_MAP     : string = getstr();                          
                       if (!strncasecmp (string,"kosambi",3))
                         usedmap = M_KOSAMBI;
                       if (!strncasecmp (string,"haldane",3))
                         usedmap = M_HALDANE;
                       printf("Using %s map function\n",mapfctnames[usedmap]);
                       break;
      case C_NAME    : string = getstr();                          
                       if (!strncasecmp (string,"trait",3))
		       { 
                         for (i=1; i<=numberoftraits(); i++)
                           SetName(traitnames, i, InputLine(buf, BUFFERSIZE));
		       }
                       else if (!strncasecmp (string,"marker",3))
		       {
                         for (i=1; i<=numberofmarkers(); i++)
                           SetName(markernames, i, InputLine(buf, BUFFERSIZE));
		       }
	               break;
    case C_HERITABILITY: printf("Which trait : ");
                         code = atoi(InputLine(buf2, BUFFERSIZE));

                         printf("Input covariates : ");
                         InputLine(buf2, BUFFERSIZE);
                         GetNumberList(buf2);

			 if (code > 0 ) {
			   if (options->Index[O_DOMINANCE]) {
			     CalculateHeritability(individuals, code, 1);
			   }
			   else 
			     CalculateHeritability(individuals, code, 0);
			 }
			 break;
      case C_MP      :
      case C_MULTI   : // Do a complete scan using a specified method
                      switch (findkeyword(getstr()))
                      {
                        case C_HASEMAN : printf("Which trait : ");
                                         code = atoi(InputLine(buf2, BUFFERSIZE));

                                         printf("Input covariates : ");
                                         InputLine(buf2, BUFFERSIZE);
                                         GetNumberList(buf2);
 
					 MultiPointHasemanElston(1, code,
                                           IntList->values, IntList->length);
                                         printf("\n");
                                         
                                         // Frees data used for covariates
                                         if (IntList->length>0)
                                           free_ivector(IntList->values, 1, IntList->length);
                                         free(IntList);

                                         break;
                        case C_VC : PrepareMultipoint(individuals);
			                 break;
 		      case C_XLINK: printf("Which trait : ");
			  code2 = atoi(InputLine(buf2, BUFFERSIZE));

			  printf("Input covariates : ");
			  InputLine(buf2, 99);
			  GetNumberList(buf2);
			  if (options->Index[O_XLINKED]) {  
			    MultiPointXVC(individuals, code2, 0, IntList->length, IntList->values);
			  }
			else {
			  MultiPointVC(individuals, code2, 0, IntList->length, IntList->values);
			}
			  break;
					 
                        default : puts("Syntax Error: Use what scanning method?"); break;
                      }
                      // End of multipoint
                      break;
    case C_SEED : if (!strcmpl(getstr(), "reset")) {
                    // The following two lines are a hack and should be removed
                    phrtsd("Bruges eksemplet", &seed1, &seed2);
		    setall(seed1,seed2);
                  }
      // Gets and prints the seeds
      getsd(&seed1,&seed2);
      printf("Current seeds: %ld  %ld\n", seed1, seed2);
      break;
    case C_PHI : MMakeKinshipMatrix(individuals, "kinship.dat");
      break;
    case C_POLYGENIC : 
      printf("Which trait : ");
      code2 = atoi(InputLine(buf2, BUFFERSIZE));
      printf("Input covariates : ");
      InputLine(buf2, 99);
      GetNumberList(buf2);
      if (options->Index[O_XLINKED]) {  
        printf("KKJLKJLK\n");
	PolygenicXVC(individuals, code2, 0, IntList->length, IntList->values);
      }
      else {
	PolygenicVC(individuals, code2, 0, IntList->length, IntList->values);
      }
      break;
    case C_GXE : GXEanalyse(); break;
      break;
      case C_SHOW   : // Show either individual or marker info
                      switch (findkeyword(getstr()))
                      {
                        case C_MARKER     : ShowMarkerInfo(atoi(getstr()));
                                            break;
                        case C_INDIVIDUAL : showindividual (getstr());
                                            break;
                        case C_DATA       : ShowDataset (individuals);
                                            break;
                        case C_OPTION     : ShowOptions ();
                                            break;
		        case C_LICENSE    : // Show the COPYING file
			                    ShowFile("COPYING");
			                    break;
                        default : puts("Syntax error: Show what?"); break;
                      }
                      break;
      case C_CLAUS  : ClausStuff(); break;
      case C_CHROM  : code = atoip(getstr());
                      if (code>0 && code<23)
			nChromosome = code;
		      break;
      case C_ORDER  : printf("Input new order (eg.: 2 3 1 4) : ");
                      OrderMarkers(InputLine(buf, BUFFERSIZE)); break;
      case C_TEST   : TestStuff(); break;
      case C_TRAIT  : string = getstr();
	              if (!strcmpl(string, "")) {
			// Print model trait if it exists. Otherwise
			if (ModelTrait) {
			  printf("Trait: %d (%s)\n", ModelTrait, GetName(traitnames, ModelTrait));
			}
			else {
			  printf("No trait selected\n");
			}	
			break;
		      }
		      code = atoip(string);
		      if (code == 0 || code>numberoftraits()) {
			printf("Illegal trait selected.\n");
			ModelTrait = 0;
		      }
		      else
			ModelTrait = code;
	              break;
      case C_COVARIATES  : string = getstr();
	                   if (!strcmpl(string, "")) {
			   // Print model trait if it exists. Otherwise
   			   if (nModelCovariates>0) {
			     printf("Covariates: ");
			     for (i = 0; i < nModelCovariates; i++) 
			       if (ModelCovariates[i]==0)
				 printf("  %d (%s)", ModelCovariates[i], "Sex");
			     else
			       printf("  %d (%s)", ModelCovariates[i], GetName(traitnames, ModelCovariates[i]));
			     printf("\n");
			   }
			   else {
			     printf("No covaraites selected\n");
			   }
			   break;
		          }
			   string = strstr(buf,igetstr(buf));
			   while (!isspace(*string))
			     string++;
			   
			   GetNumberList(string);
			   // Frees data used for covariates
			   // This is really sloppy programming. Should make a new function instead of
			   // GetNUmberList
			   if (IntList->length>0)
			     free_ivector(IntList->values, 1, IntList->length);
			   free(IntList);

	              break;
      case C_HELP   : TextColour(COL_BLUE);
                      puts ("PediPet v" VERSION "  " __DATE__ "\n");
                      TextColour(COL_WHITE);
                      puts ("CHECK         Checks the consistency of the dataset");
                      puts ("DELETE        Deletes all data in memory");
                      puts ("DISTANCE      Input the distance between markers");
                      puts ("ECHO          Echos a line of text");
                      puts ("FREQ method n Sets allele frequency for marker n using method");
                      puts ("HASEMAN n     Do Haseman-Elston regression on marker n");
                      puts ("IMPORT        Read a LINKAGE pedigree file");
                      puts ("KILL          Removes all traits and genotypes from founders");
                      puts ("MAKE filetype Make an output file of specified type");
                      puts ("MAP function  Use function as map function");
                      puts ("ORDER         Input new order for the markers");
                      puts ("RECODE        Recodes all markers");
                      puts ("SCAN method   Scan entire dataset using method");
                      puts ("SHOW MARKER n Shows marker info");
                      puts ("SHOW INDIV n  Shows person info");
                      puts ("SIMULATE fa n Simulate n families of type fam");
                      puts ("VC n          Do variance component analysis on marker n");
                      puts ("QUIT          Deletes data and quits pedipet");
                      break;
//-------------------------------------------------------------------------------------
      default       : if (strlen(buf)) {
	                puts ("PediPet v" VERSION "  " __DATE__);
			puts ("Type 'help' for help\n");
                      }
                      break;

    }
  }
 

  TextColour(COL_WHITE);
  BackgroundColour(COL_BLACK);
 
  return 0;
}
