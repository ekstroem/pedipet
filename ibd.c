/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * Functions to calculate IBD sharing between relatives
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


#ifndef IBD
#include <math.h>
#include <string.h>

#include "ibd.h"
#include "pedipet.h"
#include "haselst.h"

#define IBD
//#define debug

#define MAXGEN   10    // USED FOR RELATIONSHIP MATCHING. Maximum number of generations
#define MCSIM    1500  // Number of MC simulations
#define MAXEXACT 500000
#define FASTEXACT
// #define VERYFASTEXACT // Doesn't work either
#define IMPSAMPLE  // Doesn't work!
// #define SLOWLY

#define XLINKED        // Include the xlink hack
#define NEWXLINKED     // Include the xlink hack

// Allows IBD calculation to be performed by pedigree. Is NOT completed, and should
// only be included for the paper.
// #define BYPEDIGREE     

// Relation pair
enum {
  REL_NONE, 
  REL_ID,
  REL_FULL_SIBS, // Should be the first relationship after none and id
  REL_HALF_SIBS,
  REL_PARENT,
  REL_GRAND_PARENT,
  REL_GREAT_GRAND_PARENT,
  REL_GREAT_GREAT_GRAND_PARENT,
  REL_GREAT_GREAT_GREAT_GRAND_PARENT,
  REL_AVUNCULAR,
  REL_HALF_AVUNCULAR,
  REL_GRAND_AVUNCULAR,
  REL_HALF_GRAND_AVUNCULAR,
  REL_GREAT_GRAND_AVUNCULAR,
  REL_HALF_GREAT_GRAND_AVUNCULAR,
  REL_FIRST_COUSIN,
  REL_HALF_FIRST_COUSIN,
  REL_FIRST_COUSIN_ONCE_REMOVED,
  REL_HALF_FIRST_COUSIN_ONCE_REMOVED,
  REL_FIRST_COUSIN_TWICE_REMOVED,
  REL_HALF_FIRST_COUSIN_TWICE_REMOVED,
  REL_SECOND_COUSIN,
  REL_HALF_SECOND_COUSIN,
  REL_SECOND_COUSIN_ONCE_REMOVED,
  REL_HALF_SECOND_COUSIN_ONCE_REMOVED,
  REL_THIRD_COUSIN,
  REL_HALF_THIRD_COUSIN,
  REL_MAXREL,
};

char *RelationshipNames[REL_MAXREL] = 
{
  "None",
  "Same",
  "Full sibs",
  "Half-sibs",
  "Parent-offspring",
  "Grandparent-grandchild",
  "Great-grandparent / offspring",
  "Great-great-grandparent / offspring",
  "Great-great-great-grandparent / offspring",
  "Avuncular",
  "Half-avuncular",
  "Grand-avuncular",
  "Half-grand-avuncular",
  "Great-grand-avuncular",
  "Half-great-grand-avuncular",
  "First cousins",
  "Half-first cousins",
  "First cousin once removed",
  "Half-first cousin once removed",
  "First cousin twice removed",
  "Half-first cousin twice removed",
  "Second cousin",
  "Half-second cousin",
  "Second cousin once removed",
  "Half-second cousin once removed",
  "Third cousin",
  "Half-third cousin",
};


//
// Finds the common ancestors of ind and ind2
// Ancestral father is placed in father and likewise for mother
// Function returns the number of ancestors found
//
int FindCommonAncestor(individual *ind1, individual *ind2, individual **father, individual **mother)
{
  IDlist *anc1, *anc2, *id, *id2, *newid;
  int i, stop, newanc, maxgen, code;
  individual *ind;

  // Noone found initially
  *father = 0;
  *mother = 0;
  anc1 = 0;
  anc2 = 0;

  if (!SamePedigree(ind1, ind2))       // not from same pedigree
    return 0;

  // Adding the individuals to the list
  id = cmalloc (sizeof (IDlist));
  memset (id,0,sizeof (IDlist));
  id->ind = ind1;
  id->ind->tmpint3 = 0;
  addlist(&anc1, id);
  sprintf(ind1->tmpstr1, "1");
  sprintf(ind1->tmpstr2, "0");

  id = cmalloc (sizeof (IDlist));
  memset (id,0,sizeof (IDlist));
  id->ind = ind2;
  id->ind->tmpint3 = 0;
  addlist(&anc2, id);
  sprintf(ind2->tmpstr1, "0");
  sprintf(ind2->tmpstr2, "1");

  forind
    ind->tmpint3 = 0;

  i = 1;
  newanc = 1;
  stop = 0;
  while (!stop && newanc)
  {
    newanc = 0;

    // Adding new ancestors
    for (id = anc1; id; id = id->next)
    {
      if (atoi(id->ind->tmpstr1)==i)
      {
        if (id->ind->father)
        {
          newid = cmalloc (sizeof (IDlist));
          memset (newid,0,sizeof (IDlist));
          newid->ind = id->ind->father;
          newid->ind->tmpint3 = 0;
	  // XXX Skal muligvis checke at tmpstr == 0
          sprintf(newid->ind->tmpstr1, "%d",i+1);
          addlist(&anc1, newid);    
          newanc = 1;
        }
        if (id->ind->mother)
        {
          newid = cmalloc (sizeof (IDlist));
          memset (newid,0,sizeof (IDlist));
          newid->ind = id->ind->mother;
          newid->ind->tmpint3 = 0;
          sprintf(newid->ind->tmpstr1, "%d",i+1);
          addlist(&anc1, newid);    
          newanc = 1;
        }
      }
    }

    for (id = anc2; id; id = id->next)
    {
      if (atoi(id->ind->tmpstr2)==i)
      {
        if (id->ind->father)
        {
          newid = cmalloc (sizeof (IDlist));
          memset (newid,0,sizeof (IDlist));
          newid->ind = id->ind->father;
          sprintf(newid->ind->tmpstr2, "%d",i+1);
          addlist(&anc2, newid);    
          newanc = 1;
        }
        if (id->ind->mother)
        {
          newid = cmalloc (sizeof (IDlist));
          memset (newid,0,sizeof (IDlist));
          newid->ind = id->ind->mother;
          sprintf(newid->ind->tmpstr2, "%d",i+1);
          addlist(&anc2, newid);    
          newanc = 1;
        }

      }
    }

    i++;

    // Checking that the two ancestral lists have nothing in common
    for (id = anc1; id; id = id->next)
    {
      for (id2 = anc2; id2; id2 = id2->next)
      {
        if (id->ind==id2->ind)
        {
          stop = 1;  // Found common ancestor
          id->ind->tmpint3 = 1;
        }
      }        
    }
  }

  // Makes common set ie. keeps thoses inds with tmpint = 1
  for (id = anc1; id; )
  {
    id2 = id->next;
    if (!id->ind->tmpint3)
      removelist(&anc1, id);
    id = id2;
  }
  for (id = anc2; id; )
  {
    id2 = id->next;
    if (!id->ind->tmpint3)
      removelist(&anc2, id);
    id = id2;
  }

  // Should now only look at the youngest ancestors from the common set
  maxgen = 200000;
  for (id = anc1; id; id = id->next)
    maxgen = min(maxgen, atoi(id->ind->tmpstr1));

  // Found the min generation.
  // Removes all other generations
  for (id = anc1; id; )
  {
    id2 = id->next;
    if (atoi(id->ind->tmpstr1)>maxgen)
      removelist(&anc1, id);
    id = id2;
  }


  maxgen = 200000;
  for (id = anc2; id; id = id->next)
    maxgen = min(maxgen, atoi(id->ind->tmpstr2));

  // Found the min generation.
  // Removes all other generations
  for (id = anc2; id; )
  {
    id2 = id->next;
    if (atoi(id->ind->tmpstr2)>maxgen)
      removelist(&anc2, id);
    id = id2;
  }

  code = 0;
  for (id = anc1; id; id = id->next)
  {
    if (id->ind->sex==S_MALE)
    {
      *father = id->ind;
      code++;
    }
    if (id->ind->sex==S_FEMALE)
    {
      *mother = id->ind;
      code++;
    }
  }


  freelist(anc1);
  freelist(anc2);

#ifdef debug
  if (code>3)
    printf("ERROR: Something is wrong in rutine FindCommonAncestor\n");
#endif

  // Returns the number of found ancestors
  return code;
}

// Creates a linked list of individuals, that connects
// an individual to his/her ancestor
// Input: (ind) individual and his/her ancestor (anc)
IDlist *LinkToAncestor(individual *ind, individual *anc)
{
  RIDlist *tree, *id, *newid, *final;
  IDlist *res, *resid;
  int stop, newanc;

  tree = 0;
  res = 0;
  final = NULL;

  // Adding the individuals to the list
  id = cmalloc (sizeof (RIDlist));
  memset (id, 0, sizeof (RIDlist));
  id->ind = ind;
  addlist(&tree, id);

  if (ind==anc)
  {
    resid = cmalloc (sizeof (IDlist));
    memset (resid,0,sizeof (IDlist));
    resid->ind = id->ind;
    addlist(&res, resid);   

    freelist(tree);
    return res;
  }

  newanc = 1;
  stop = 1;

  while (newanc && stop)
  {
    newanc = 0;

    // Adding new ancestors
    for (id = tree; id; id = id->next)
    {
      if (id->ind->father)
      {
        newid = cmalloc (sizeof (RIDlist));
        memset (newid,0,sizeof (RIDlist));
        newid->ind = id->ind->father;
        addlist(&tree, newid);   
        newid->prev = id;
        newanc = 1;
        if (newid->ind==anc)
        {
          stop = 0;
          final = newid;
        }
      }
      if (id->ind->mother)
      {
        newid = cmalloc (sizeof (RIDlist));
        memset (newid,0,sizeof (RIDlist));
        newid->ind = id->ind->mother;
        addlist(&tree, newid);   
        newid->prev = id;
        newanc = 1;
        if (newid->ind==anc)
        {
          stop = 0;
          final = newid;
        }
      }
    }
  }

  // Now has a tree with both directions complete
  // Sets the persons to keep
  // Assume all should be removed
  for (id = tree; id; id = id->next)
    id->ind->tmpint3 = 0;

  // Keep those in line from ind to anc
  for (id = final; id; id = id->prev)
    id->ind->tmpint3 = 1;

  // Removes the rest
  for (id = tree; id; )
  {
    newid = id->next;
    if (!id->ind->tmpint3)
      removelist(&tree, id);
    id = newid;
  }

  res = 0;
  for (id = tree; id; id = id->next)
  {
    resid = cmalloc (sizeof (IDlist));
    memset (resid,0,sizeof (IDlist));
    resid->ind = id->ind;
    // Adds the info on the relationship to the next person in the linked list
    // 1 for fathers and 2 for mothers
    resid->ind->tmpint2 = 0;

    if (id->next)
    {
      if (id->next->ind == id->ind->father)
        resid->ind->tmpint2 = 1;
      else
        resid->ind->tmpint2 = 2;
    }
    addlist(&res, resid);   
  }

  freelist(tree);
  return(res);
}




// Return the generation between ind and ind2. -1 means no descendant
// ind is assumed to be the youngest ie. is ind a direct desc of ind2

// XXXX Right now only return 1 for DD and -1 for not DD
int DirectDescendant(individual *ind, individual *ind2)
{
  individual *fat, *mot;
  IDlist *id;
  int res, n;

  fat = 0;
  mot = 0;
  res = FindCommonAncestor(ind, ind2, &fat, &mot);

  if (res != 1)
    return -1;
  else
  {
    if (ind==fat || ind==mot)
    {
      if (fat)
        id = LinkToAncestor(ind2, fat);
      else
        id = LinkToAncestor(ind2, mot);
      n = listlen(id)-1;
      freelist(id);
      return n;
    }
    if (ind2==fat || ind2==mot)
    {
      if (fat)
        id = LinkToAncestor(ind, fat);
      else
        id = LinkToAncestor(ind, mot);
      n = listlen(id)-1;
      freelist(id);
      return (n);
    }
    return -1;
  }
}



// Creates a copy of a pedigree with the possible ordered genotypes for
// a specific marker, and then does genotype elimination
// Input: pedigree: Name of pedigree to create
//        markernum: Number of marker
individual *GenotypeEliminationFromList(char *pedigree, int markernum, individual *indlist)
{
  markerlist *mkr, *marker;
  individual *ind, *ind2, *minilist, *lind;
  int i, j, k, numberofalleles;
  FreqList *fl;

  numberofalleles = NumberOfAlleles(markernum);

  fl = FrequencyNumber(markernum);

//  printf("Starting genotype elimination for pedigree %s (%d)\n", pedigree, numberofalleles);

  // First make a copy of the examined pedigree
  minilist = 0;

  for (ind=indlist; ind; ind=ind->next) {
    // Is the individual from the correct pedigree?
    if (!strcmpl(pedigree, ind->pedigree)) {
      ind2 = cmalloc (sizeof (individual));
      memset(ind2, 0, sizeof(individual));

      strcpy(ind2->id, ind->id);
      strcpy(ind2->pedigree, ind->pedigree);
      ind2->sex = ind->sex;
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

      mkr = markernumber(ind, markernum);
      // If the person is typed
      // Then let the possible alleles be the observed alleles
      if ((mkr->allele1>0) && (mkr->allele2>0)) {
        // If homozygous
        if (homozygous(ind,markernum)) {
          marker = cmalloc (sizeof (markerlist));
          memset (marker,0,sizeof (markerlist));

          marker->allele1 = mkr->allele1;
          marker->allele2 = mkr->allele2;

#ifdef XLINKED
	  // Setting one of the alleles for males as a pseudo marker
	  if (options->Index[O_XLINKED] && ind2->sex == S_MALE) {
	    marker->allele1 = numberofalleles+1;
	  }
#endif
          addlist(&ind2->marker, marker);
        }
	else // Heterozygote
        {
          marker = cmalloc (sizeof (markerlist));
          memset (marker,0,sizeof (markerlist));

          marker->allele1 = mkr->allele1;
          marker->allele2 = mkr->allele2;
          addlist(&ind2->marker, marker);
 
          marker = cmalloc (sizeof (markerlist));
          memset (marker,0,sizeof (markerlist));

          marker->allele1 = mkr->allele2;
          marker->allele2 = mkr->allele1;

#ifdef XLINKED
	  // Found a heterozygous male for an xlinked marker. Ouch!
	  if (options->Index[O_XLINKED] && ind2->sex == S_MALE) {
	    printf("WARNING: Heterozygous male (%s) for marker %d\n", ind2->id, markernum);
	  }
#endif
          addlist(&ind2->marker, marker);
        }
      }
      else { // The person is not typed and has all possible genotypes
        for (i=1; i<=numberofalleles; i++) {
	  // Goes through all possible combinations of alleles
	  // That includes numberofalleles for allele1 and
	  // numberofalleles for allele2
	  // For X-linked markers in males, the options for allele1 (the paternal allele)
	  // is automatically set to the pseudomarker
	  k = numberofalleles; // default

#ifdef XLINKED
	  if (options->Index[O_XLINKED] && ind2->sex == S_MALE) {
	    k = 1;
	  }
#endif

	  if (fl->frequency[i] ==0)
	    continue;

          for (j=1; j<=k; j++) {

	    if (fl->frequency[j] == 0 && !(options->Index[O_XLINKED] && ind2->sex == S_MALE))
	      continue;


            marker = cmalloc (sizeof (markerlist));
            memset (marker,0,sizeof (markerlist));
            marker->allele1 = i;
            marker->allele2 = j;

#ifdef XLINKED
	    if (options->Index[O_XLINKED] && ind2->sex == S_MALE) {
	      marker->allele2 = marker->allele1;
	      marker->allele1 = numberofalleles+1;
	    }
#endif

            addlist(&ind2->marker, marker);
          }
        }
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
    // Fixing sibs
    
    if (!founder(ind))
    {

      // Adds spouses for parents
      adduniqueidlist(&ind->father->mate, ind->mother);
      adduniqueidlist(&ind->mother->mate, ind->father);

      for (lind = ind->next; lind; lind = lind->next)
      {
        if (((ind->father == lind->father) && (ind->father)) ||
            ((ind->mother == lind->mother) && (ind->mother)))
        {
          adduniqueidlist(&ind->sib, lind);
          adduniqueidlist(&lind->sib, ind);
        }
      }
    }
  }

  // Have now prepared the dataset and should start eliminating genotypes
  ReduceGenotypes(minilist);

  return (minilist);
}

individual *GenotypeElimination(char *pedigree, int markernum)
{
  return(GenotypeEliminationFromList(pedigree, markernum, individuals));
}



//
// Reduce an ordered genotype list by genotype elimination
//
double ReduceGenotypes(individual *indlist)
{
  individual *ind, *lind;
  markerlist *fmkr, *mmkr, *marker, *mkrptr, *mkr2;
  IDlist *sib;
  int endelim, keep, childok, allchildok;
  double combinations;

//  markerlist *mkr, *mkr2;

  combinations = 0;
  endelim = 0;
  allchildok = 1;

  do
  {
    endelim = 1;

    for (ind = indlist; ind; ind = ind->next) {
      // Looking at each tripel (indexed by the children [non-founders])
      // This is extremely redundant, as we only need to look at 
      // nuclear families.
      // However - this is not being fixed right now

      if (!founder(ind))
      {

        // Starts by examining fathers genotypes
        for (fmkr = ind->father->marker; fmkr; )
        {

          mkrptr = fmkr->next;

	  // Wants to see if fmkr should be kept or not.
	  // If combined with one of mothers genotypes and
	  // all children then are legal, then keep fathers genotype
          for (mmkr = ind->mother->marker; mmkr; mmkr = mmkr->next)
          {
            allchildok = 1; // Assume all children ok
            childok    = 0; // Assume this child is not ok

	    // Could the child come from parental genotype combination
            for (marker = ind->marker; marker; marker = marker->next)
            {
              if ((marker->allele1==fmkr->allele1 || marker->allele1==fmkr->allele2) &&
                  (marker->allele2==mmkr->allele1 || marker->allele2==mmkr->allele2))
              {
                childok = 1;
                break;
              }
            }

            allchildok *= childok;

	    // The child could not be from than combination, so check if
	    // fathers genotype is better together with next maternal marker
            if (!allchildok)
              continue; 

            // Shold now check the other sibs of the nuclear family
            for (sib = ind->sib; sib; )
            {
              lind = sib->ind;
              if (fullsibs(ind, lind))
              {

                childok = 0;
                for (mkr2 = lind->marker; mkr2; mkr2 = mkr2->next)
                {
                  if ((mkr2->allele1==fmkr->allele1 || mkr2->allele1==fmkr->allele2) &&
                      (mkr2->allele2==mmkr->allele1 || mkr2->allele2==mmkr->allele2))
                  {
                     childok=1;
                     break;
                  }
                }
                allchildok *= childok;
		// allchildok is 1 if all children could come from parents and
		// 0 otherwise
              }

              if (!allchildok)
                sib = NULL;
              else
                sib = sib->next;
            }

            if (!allchildok)
              continue; // Go to next maternal marker

            // If we are here, then each child in the nuclear family 
            // are ok, and the fathers genotype should be saved

            if (allchildok)
              break;


          } // end-of-for mmkr

          if (!allchildok)
	  {
            removelist(&ind->father->marker, fmkr);	    
	    endelim = 0;
	  }

          fmkr = mkrptr;
        }


        // Do the same for the mothers genotypes
        for (mmkr = ind->mother->marker; mmkr; )
        {
          mkrptr = mmkr->next;

          for (fmkr = ind->father->marker; fmkr; fmkr = fmkr->next)
          {
            allchildok = 1; // Assume all children ok
            childok    = 0; // Assume this child is not ok


            for (marker = ind->marker; marker; marker = marker->next)
            {
              if ((marker->allele1==fmkr->allele1 || marker->allele1==fmkr->allele2) &&
                  (marker->allele2==mmkr->allele1 || marker->allele2==mmkr->allele2))
              {
                childok = 1;
                break;
              }
            }

            allchildok *= childok;
            
            if (!allchildok)
              continue; // Go to next maternal marker

            // Should now check the other sibs
            for (sib = ind->sib; sib; )
            {
              lind = sib->ind;
              if (fullsibs(ind, lind))
              {

                childok = 0;
                for (mkr2 = lind->marker; mkr2; mkr2 = mkr2->next)
                {
                  if ((mkr2->allele1==fmkr->allele1 || mkr2->allele1==fmkr->allele2) &&
                      (mkr2->allele2==mmkr->allele1 || mkr2->allele2==mmkr->allele2))
                  {
                     childok=1;
                     break;
                  }
                }
                allchildok *= childok;
              }

              if (!allchildok)
                sib = NULL;
              else
                sib = sib->next;
            }

            if (!allchildok)
            {
              continue; // Go to next maternal marker
            }

            // If we are here, then each child in the nuclear family 
            // are ok, and the fathers genotype should be saved

            if (allchildok)
              break;


          } // end-of-for fmkr

          if (!allchildok)
	  {
            removelist(&ind->mother->marker, mmkr);
	    endelim = 0;
	  }

          mmkr = mkrptr;
        }


        // Removes ordered genotypes for non-founders if they are
        // incompatible with the parents ordered genotypes

	for (marker = ind->marker; marker;)
        {
          keep = 0;
          mkr2 = marker->next;

          for (fmkr = ind->father->marker; fmkr; fmkr = fmkr->next)
          {
            if (marker->allele1==fmkr->allele1 || marker->allele1==fmkr->allele2)
               keep=1;
          }
          if (!keep)
          {
	    removelist(&ind->marker,marker);
            endelim = 0;
          }
          marker = mkr2;
        }

	for (marker = ind->marker; marker; )
        {
          keep = 0;
          mkr2 = marker->next; 

          for (mmkr = ind->mother->marker; mmkr; mmkr = mmkr->next)
          {
            if (marker->allele2==mmkr->allele1 || marker->allele2==mmkr->allele2)
               keep=1;
          }
          if (!keep)
          {
            fmkr = marker->next;
	    removelist(&ind->marker,marker);
            endelim = 0;
          }
          marker = mkr2;
        }
      }
    }
  } while (!endelim);

  combinations = 1;
  for (ind = indlist; ind; ind = ind->next)
  {
    combinations *= listlen(ind->marker);
  }

  return (combinations);
}

/*

  All this tree notation is not too smart since a parental tree (e.g. the root)
  corresponds to the youngest person (ie. the person from the newest generation), while
  a tree's child is the persons parent :-/

  Shame on you Claus

 */

double RecursiveCalc(MarkerTree *parent, MarkerTree *child)
{
  MarkerTree *mt;
  int res;

  res = 0;

  if (parent->allele[parent->fromparent-1]==child->allele[child->fromparent-1])
  {
    // If the end
    if (!child->child)
      return 1;
    for (mt = child->child; mt; mt = mt->next)
      res += RecursiveCalc(child, child->child);
    return res;
  }
  else 
    return 0;
}

void AddChildToTree(MarkerTree *parent, MarkerTree *child)
{
  child->parent = parent;  
  addlist(&parent->child, child);
}

double CalcIBDSharing(individual *ind, individual *ind2, int markernum)
{
  MarkerTree *mt, *subtree, *generation, *root;
  markerlist *mkr;
  IDlist *chain1, *chain2, *id;
  double res; // The probability of sharing 1 allele
  int i;
  int ind1komb, ind2komb;  // # of combinations for each individual
  int ch1komb, ch2komb;    // # of combinations for eachchain, without the ind and the ancestor
  int nanc, chain1len, chain2len, homozygote;
  long totalcomb;
  individual *genolist, *fat, *mot, *anc;


  // Gets the allowed genotypes
  genolist = GenotypeElimination(ind->pedigree, markernum);

  // Finds the common ancestors
  nanc = FindCommonAncestor(ind, ind2, &fat, &mot);
  printf("Common ancestors: %d\n", nanc);

  // If no common ancestors
  if (nanc == 0)
    return 0.0;

  if (fat)
   anc = fat;
  else
   anc = mot;

  // Figures out the linkage chains
  chain1 = LinkToAncestor(ind, anc);
  chain2 = LinkToAncestor(ind2, anc);

  ind1komb = listlen(FindIndListMember(genolist, ind->id)->marker);
  ind2komb = listlen(FindIndListMember(genolist, ind2->id)->marker);

  ch1komb = 1;
  ch2komb = 1;

  for (id = chain1->next; id && id->next; id = id->next)
    ch1komb *= listlen(id->ind->marker);
  for (id = chain2->next; id && id->next; id = id->next)
    ch2komb *= listlen(id->ind->marker);

  chain1len = listlen(chain1);
  chain2len = listlen(chain2);

  printf("Chainlen: %d  %d   -   \n", chain1len, chain2len);

  root = cmalloc (sizeof (MarkerTree));
  memset (root, 0, sizeof (MarkerTree));

  generation = root;

  // Should now fill up the tree
  totalcomb = 1;
  homozygote = 1;
  for (id = chain1; id && id->next; id = id->next)
  {
    // Makes a subtree
    subtree = 0;
    i = 1;
    for (mkr = FindIndListMember(genolist, id->ind->id)->marker; mkr; mkr = mkr->next)
    {
      mt = cmalloc (sizeof (MarkerTree));
      memset (mt,0,sizeof (MarkerTree));

      mt->allele[0] = mkr->allele1;
      mt->allele[1] = mkr->allele2;

      if (!strcmpl(id->next->ind->id, id->ind->father->id))
        mt->fromparent = 1;
      else
        mt->fromparent = 2;
      if (mkr->allele1 == mkr->allele2)
      {
        mt->homozygote = 1;
        homozygote *=2;
      }
       

      mt->next = subtree;
      subtree = mt;
      i++;
    }
    totalcomb *= (i-1);

    //    printf("%d\n",listlen(subtree));
    for (mt = subtree; mt; mt = mt->next)
      AddChildToTree(generation, mt);
    generation = subtree;
    //    printf("End of this part\n");
  }

  res = 1.0;
  for (mt = root->child; mt; mt = mt->next)
    res += RecursiveCalc(root, mt);

  

  /*
    // Have now made a subtree with linear links to "sibs"

    for (mt2 = tree; mt2; mt2 = mt2->sib)    
      mt2->child = subtree;
    tree = subtree;
  }
  

  */

  freelist(genolist);

  return ((float)res/(totalcomb*homozygote));
}



OLDMATRIX *PedigreeIBD (individual *pedigree)
{
  int i, j, update, nanc;
  unsigned long elements, longi;
  int smartarray[MAXPERSONS], maxarray[MAXPERSONS];
  int savetime[MAXPERSONS][5];
  double tempres, tempres2;
  double founderprob, extracheck;
  IDlist *chain1, *chain2, *id, *chain3, *chain4;
  markerlist *mkr, *fmkr, *mmkr, *mkr2;
  individual *fat, *mot, *anc;
  individual *ind, *ind1, *ind2;
  OLDMATRIX *resmat, *tmpmat, *share1, *share2;

#ifdef VERYFASTEXACT
  int VFsmartarray[MAXPERSONS], VFmaxarray[MAXPERSONS];
#endif

  // Calculating the number of combinations required for this
  // pedigree
  elements   = 1;
  extracheck = 1;
  j = 0;

  i = 0;

#ifdef FASTEXACT
  // Reduce the number of required combinations by checking if founders
  // have both copies of an ordered pair present.

  for (ind = pedigree; ind; ind = ind->next) {
    if (founder(ind) && listlen(ind->marker)>1) {
      for (mkr = ind->marker; mkr; ) {
        mkr2 = mkr->next;
	for (fmkr = mkr->next; fmkr; ) {
	  mmkr = fmkr->next;
          if (mkr->allele1 == fmkr->allele2 && mkr->allele2 == fmkr->allele1) {
            if (fmkr==mkr2)
              mkr2 = fmkr->next;
	    removelist(&ind->marker, fmkr);
	  }
	  fmkr = mmkr;
	}
        mkr = mkr2;
      }
    }
  }
#endif

  // Preparing the individuals
  for (ind = pedigree; ind; ind = ind->next)
  {
    maxarray[j] = listlen(ind->marker);
    elements   *= maxarray[j];
    extracheck *= maxarray[j];
    j++;
    // tmpint holds the internal number: 1 for first person, 2 for second etc.
    ind->tmpint1 = j;
  }

#ifdef debug
  printf("Calculated %ld required computations for family %s\n", elements, pedigree->name);
#endif

  if (extracheck != (double)elements)
    printf("WARNING: Too many untyped individuals in pedigree %s to do exact IBD computations\n", pedigree->pedigree);

  // Preparing the matrices to hold the results
  // tmpmat holds the total mass possible, while share1 and share2
  // holds the mass for sharing 1 and 2 alleles IBD
  resmat = MtxNew(listlen(pedigree), listlen(pedigree));
  tmpmat = MtxNew(listlen(pedigree), listlen(pedigree));
  share1 = MtxNew(listlen(pedigree), listlen(pedigree));
  share2 = MtxNew(listlen(pedigree), listlen(pedigree));

  for (ind1 = pedigree; ind1; ind1 = ind1->next)
  {
    for (ind2 = ind1; ind2; ind2 = ind2->next)
    {

	// If the individuals are the same then IBD sharing is exactly 1
        if (ind1==ind2)
	{
          share2->element[ind1->tmpint1][ind2->tmpint1] = 1;
          tmpmat->element[ind1->tmpint1][ind2->tmpint1] = 1;
          continue;
	}

	// If direct descendent then exactly 1/2
        if (ind1->father==ind2 || ind1->mother==ind2 || ind1==ind2->father || ind1==ind2->mother)
	{
          share1->element[ind1->tmpint1][ind2->tmpint1] = 1;
          tmpmat->element[ind1->tmpint1][ind2->tmpint1] = 1;
          continue;
	}

	// Two founders are always (assumed to be) unrelated
        if (founder(ind1) && founder(ind2))
	{
          share1->element[ind1->tmpint1][ind2->tmpint1] = 0;
          tmpmat->element[ind1->tmpint1][ind2->tmpint1] = 1;
          continue;
	}

	// Do the persons have any common ancestors
        nanc = FindCommonAncestor(ind1, ind2, &fat, &mot);

        switch (nanc)
	{
	  case 0 : // No common ancestors, so IBD = 0
                   tmpmat->element[ind1->tmpint1][ind2->tmpint1] = 1;
		   break;
	  case 1 : // One common ancestor

	           // Figuring out if the father of the mother is the common ancestor
                   if (fat)
                     anc = fat;
                   else
                     anc = mot;

		   // Calculate the chain from each individual to the common ancestor
		   chain1 = 0;
		   chain2 = 0;
                   chain1 = LinkToAncestor(ind1, anc);
                   chain2 = LinkToAncestor(ind2, anc);

		   // Readying the smart array
		   memset(smartarray,0, sizeof(smartarray));

#ifdef VERYFASTEXACT

		   elements = 1;
		   i = 0;
                   for (id = chain1; id; id = id->next)
		   {
		     maxarray[i] = listlen(id->ind->marker);
		     elements *= VFmaxarray[i];
		     i++;
		     id->ind->tmpint1 = i;
		   }
                   for (id = chain2; id && id->next; id = id->next)
		   {
		     maxarray[i] = listlen(id->ind->marker);
		     elements *= VFmaxarray[i];
		     i++;
		     id->ind->tmpint1 = i;
		   }

#endif

	           // Goes through all possible genotype combinations for the pedigree
                   for (longi = 0; longi< elements; longi++) {
		     tempres = 1.0;

                     // Calculates the probability of the pedigree given the ordered pairs of genotypes
		     for (ind = pedigree; ind; ind = ind->next)
		     {
		       mkr = markernumber(ind, smartarray[ind->tmpint1-1]+1);
		       if (founder(ind)) {
#ifdef NEWXLINKED
			 if (options->Index[O_XLINKED] && ind->sex == S_MALE) {
			   tempres *= allelefreq[mkr->allele2];
			 }
			 else {
			   tempres *= allelefreq[mkr->allele1]*allelefreq[mkr->allele2];
			 }
#else
			 tempres *= allelefreq[mkr->allele1]*allelefreq[mkr->allele2];
#endif
#ifdef FASTEXACT
#ifdef NEWXLINKED
			 if (options->Index[O_XLINKED] && ind->sex == S_FEMALE && mkr->allele1 != mkr->allele2) {
#else
 		         if (mkr->allele1 != mkr->allele2) {
#endif
                           tempres *=2.0;
			 }
#endif
		       }
		       else  // Check that the offsprings genotype is legal
		       {
			 fmkr = markernumber(ind->father, smartarray[ind->father->tmpint1-1]+1);
			 mmkr = markernumber(ind->mother, smartarray[ind->mother->tmpint1-1]+1);           
			 // if not consistent with the parents 
			 if (!((mkr->allele1 == fmkr->allele1 || mkr->allele1 == fmkr->allele2) &&
			    (mkr->allele2 == mmkr->allele1 || mkr->allele2 == mmkr->allele2)))
			   tempres = 0.0;
                         else
			 {
                           if (fmkr->allele1 != fmkr->allele2)
                             tempres /= 2;
                           if (mmkr->allele1 != mmkr->allele2)
                             tempres /= 2;
			 }
		       }
		     }
		     // tmpmat holds the prob of the likelihood
		     tmpmat->element[ind1->tmpint1][ind2->tmpint1] += tempres;
		     founderprob = tempres;

		     // Get the ordered genotypes from ind1 and in2 for this combination
                     mkr  = markernumber(chain1->ind, smartarray[chain1->ind->tmpint1-1]+1);
                     mkr2 = markernumber(chain2->ind, smartarray[chain2->ind->tmpint1-1]+1);

                     // Verifying that the two alleles examined are IBS
                     j = 10*(chain1->ind->tmpint2)+(chain2->ind->tmpint2);
                     switch(j)
		     {
		       case  0: // These first cases should only occur, when
		       case  1: // one of the two individuals is the common ancestor
		       case  2: 
		       case 10:
		       case 20: break;
		       case 11: // Both links through fathers
			        if (mkr->allele1 != mkr2->allele1) tempres = 0.0; break;
		       case 12: // Fath Moth
			        if (mkr->allele1 != mkr2->allele2) tempres = 0.0; break;
		       case 21: // Moth Father
			        if (mkr->allele2 != mkr2->allele1) tempres = 0.0; break;
		       case 22: // Both from mothers
			        if (mkr->allele2 != mkr2->allele2) tempres = 0.0; break;
		     }

		     // End the two persons are not even IBS then they share 0
                     if (tempres > 0.0)
		     {

                       for (id = chain1; id && id->next; id = id->next)
                       {  
                         mkr = markernumber(id->ind, smartarray[id->ind->tmpint1-1]+1);
	   	         // fmkr is not necessarily the father.
		         // Just the next person in the list
                         fmkr = markernumber(id->next->ind, smartarray[id->next->ind->tmpint1-1]+1);

	 	         // If that person is homozygote the divide the result by 2
		         // For all non-start and non-end persons in the chain
                         if (id->next->next && fmkr->allele1 == fmkr->allele2)
                           tempres /= 2;

                         if (id->next->ind==anc)
		         {
			   // Only add after correcting for the other chain
                           // share1->element[ind1->tmpint1][ind2->tmpint1] += tempres;
                           break;
  		         }

		         // printf("Ordered genotypes : %d %d    <->    %d %d\n", mkr->allele1, 
		         //        mkr->allele2, fmkr->allele1, fmkr->allele2);

                         if (id->ind->tmpint2==1)
  		         {
			   // If from fathers father and not alike
                           if ((id->next->ind->tmpint2==1 && mkr->allele1 != fmkr->allele1) ||
                              (id->next->ind->tmpint2==2 && mkr->allele1 != fmkr->allele2))
                             tempres = 0.0;
	  		 }
			 else
		         {
                           if ((id->next->ind->tmpint2==1 && mkr->allele2 != fmkr->allele1) ||
                              (id->next->ind->tmpint2==2 && mkr->allele2 != fmkr->allele2))
                             tempres = 0.0;
			 }
		       }
		     

                       for (id = chain2; id && id->next; id = id->next)
                       {  
                         mkr = markernumber(id->ind, smartarray[id->ind->tmpint1-1]+1);
		         // fmkr is not necessarily the father. 
		         // Just the next person in the list
                         fmkr = markernumber(id->next->ind, smartarray[id->next->ind->tmpint1-1]+1);

		         // If that person is homozygote the divide the result by 2
                         if (id->next->next && fmkr->allele1 == fmkr->allele2)
                           tempres /= 2;

                         if (id->next->ind==anc)
	 	         {
                           share1->element[ind1->tmpint1][ind2->tmpint1] += tempres;
                           break;
		         }

                         if (id->ind->tmpint2==1)
    		         {
			   // If from fathers father and not alike
                           if ((id->next->ind->tmpint2==1 && mkr->allele1 != fmkr->allele1) ||
                              (id->next->ind->tmpint2==2 && mkr->allele1 != fmkr->allele2))
                             tempres = 0.0;
		         }
		         else
	  	         {
                           if ((id->next->ind->tmpint2==1 && mkr->allele2 != fmkr->allele1) ||
                              (id->next->ind->tmpint2==2 && mkr->allele2 != fmkr->allele2))
                             tempres = 0.0;
		         }
                       }		   
		     }
		     
		     // If none of the two persons are the common ancestor
		     // and the ancestor is homozygous then divide the prob by 2
                     if (ind1 != anc && ind2 != anc && (listlen(anc->marker)==1))
		     {
                       mkr  = markernumber(anc, 1);
                       if (mkr->allele1 == mkr->allele2)
                         tempres /= 2;
		     }

                     // Should now update the smartarray
		     update=0;
		     j=0; // First element
		     while (!update)
		     {
		       smartarray[j]++;
		       if (smartarray[j]==maxarray[j])
		       {
			 smartarray[j]=0;
			 if (i<(elements-1))
			   j++;
			 else
			   update=1;
		       }
		       else
			 update=1;
		     }
		   }

	           freelist(chain1);
	           freelist(chain2);
		   chain1 = NULL;
		   chain2 = NULL;

		   break;
	  case 2:  // Readying the smart array

		   // Calculate the chain from each individual to the common ancestor
		   chain1 = 0;
		   chain2 = 0;
		   chain3 = 0;
		   chain4 = 0;

                   chain1 = LinkToAncestor(ind1, fat);
                   chain2 = LinkToAncestor(ind2, fat);

                   for (id = chain1; id; id = id->next)
                     savetime[id->ind->tmpint1][1] = id->ind->tmpint2;
                   for (id = chain2; id; id = id->next)
                     savetime[id->ind->tmpint1][2] = id->ind->tmpint2;

		   chain3 = LinkToAncestor(ind1, mot);
		   chain4 = LinkToAncestor(ind2, mot);

                   for (id = chain3; id; id = id->next)
                     savetime[id->ind->tmpint1][3] = id->ind->tmpint2;	     
                   for (id = chain4; id; id = id->next)
                     savetime[id->ind->tmpint1][4] = id->ind->tmpint2;

		   memset(smartarray,0, sizeof(smartarray));

	           // Goes through all possible genotype combinations for the pedigree
                   for (longi = 0; longi< elements; longi++) {
		     tempres  = 1.0;
                     tempres2 = 1.0;

                     // Calculates the probability of the pedigree given the ordered pairs of genotypes
		     for (ind = pedigree; ind; ind = ind->next)
		     {
		       mkr = markernumber(ind, smartarray[ind->tmpint1-1]+1);
		       if (founder(ind)) {
#ifdef NEWXLINKED
			 if (options->Index[O_XLINKED] && ind->sex == S_MALE) {
			   tempres *= allelefreq[mkr->allele2];
			 }
			 else {
			   tempres *= allelefreq[mkr->allele1]*allelefreq[mkr->allele2];
			 }
#else
			 tempres *= allelefreq[mkr->allele1]*allelefreq[mkr->allele2];
#endif

#ifdef FASTEXACT
#ifdef NEWXLINKED
			 if (options->Index[O_XLINKED] && ind->sex == S_FEMALE && mkr->allele1 != mkr->allele2) {
#else
 		         if (mkr->allele1 != mkr->allele2) {
#endif
                           tempres *=2.0;
			 }
#endif
		       }
		       else  // Check that the offsprings genotype is legal
		       {
			 fmkr = markernumber(ind->father, smartarray[ind->father->tmpint1-1]+1);
			 mmkr = markernumber(ind->mother, smartarray[ind->mother->tmpint1-1]+1);           
			 if (!((mkr->allele1 == fmkr->allele1 || mkr->allele1 == fmkr->allele2) &&
			    (mkr->allele2 == mmkr->allele1 || mkr->allele2 == mmkr->allele2)))
			   {
			     tempres = 0.0;         
			     tempres2 = 0.0;
			   }
                         else
			 {
			   // If the parent is heterozygous
                           if (fmkr->allele1 != fmkr->allele2)
                             tempres /= 2;
                           if (mmkr->allele1 != mmkr->allele2)
                             tempres /= 2;
			 }
		       }
		     }

		     // tmpmat holds the prob of the observed likelihood
		     tmpmat->element[ind1->tmpint1][ind2->tmpint1] += tempres;
		     founderprob = tempres;

                     mkr  = markernumber(chain1->ind, smartarray[chain1->ind->tmpint1-1]+1);
                     mkr2 = markernumber(chain2->ind, smartarray[chain2->ind->tmpint1-1]+1);

                     // Verifying that the two alleles examined are IBS
                     j = 10*(savetime[chain1->ind->tmpint1][1])+(savetime[chain2->ind->tmpint1][2]);

		     switch(j)
	             {
		       case  0: // These first cases should only occur, when
		       case  1: // one of the two individuals is the common ancestor
		       case  2: 
		       case 10:
		       case 20: break;
		       case 11: // Both links from fathers
			       if (mkr->allele1 != mkr2->allele1) tempres = 0.0; break;
		       case 12: // Moth Fath
			       if (mkr->allele1 != mkr2->allele2) tempres = 0.0; break;
		       case 21: // Fath Moth
			       if (mkr->allele2 != mkr2->allele1) tempres = 0.0; break;
		       case 22: // Both from mothers
			       if (mkr->allele2 != mkr2->allele2) tempres = 0.0; break;
		     }

		     // If they're IBS then continue, otherwise they're IBD=0
                     if (tempres > 0.0)
		     {
                       for (id = chain1; id && id->next; id = id->next)
                       {  
                         mkr = markernumber(id->ind, smartarray[id->ind->tmpint1-1]+1);
			 // fmkr is not necessarily the father. 
			 // Just the next person in the list
                         fmkr = markernumber(id->next->ind, smartarray[id->next->ind->tmpint1-1]+1);

			 // If that person is homozygote the divide the result by 2
			 // For all non-start and non-end persons in the chain
                         if (id->next->next && fmkr->allele1 == fmkr->allele2)
                           tempres /= 2;

                         if (id->next->ind==fat)
			 {
                           break;
			 }

                         if (savetime[id->ind->tmpint1][1]==1)
  		         {
			   // If from fathers father and not alike
                           if ((savetime[id->next->ind->tmpint1][1] ==1 && mkr->allele1 != fmkr->allele1) ||
                               (savetime[id->next->ind->tmpint1][1] ==2 && mkr->allele1 != fmkr->allele2))
                             tempres = 0.0;
			 }
			 else
		         {
                           if ((savetime[id->next->ind->tmpint1][1] ==1 && mkr->allele2 != fmkr->allele1) ||
                               (savetime[id->next->ind->tmpint1][1]==2 && mkr->allele2 != fmkr->allele2))
                             tempres = 0.0;
			 }
		       }

                       for (id = chain2; id && id->next; id = id->next)
                       {  
                         mkr = markernumber(id->ind, smartarray[id->ind->tmpint1-1]+1);
			 // fmkr is not necessarily the father. 
			 // Just the next person in the list
                         fmkr = markernumber(id->next->ind, smartarray[id->next->ind->tmpint1-1]+1);

			 // If that person is homozygote the divide the result by 2
                         if (id->next->next && fmkr->allele1 == fmkr->allele2)
                           tempres /= 2;

                         if (id->next->ind==fat)
			 {
                           // share1->element[ind1->tmpint1][ind2->tmpint1] += tempres;
                           break;
			 }

                         if (savetime[id->ind->tmpint1][2]==1)
  		         {
			   // If from fathers father and not alike
                           if ((savetime[id->next->ind->tmpint1][2] ==1 && mkr->allele1 != fmkr->allele1) ||
                               (savetime[id->next->ind->tmpint1][2] ==2 && mkr->allele1 != fmkr->allele2))
                             tempres = 0.0;
			 }
			 else
		         {
                           if ((savetime[id->next->ind->tmpint1][2] ==1 && mkr->allele2 != fmkr->allele1) ||
                               (savetime[id->next->ind->tmpint1][2] ==2 && mkr->allele2 != fmkr->allele2))
                             tempres = 0.0;
			 }
		       }
		       // If none of the two persons are the common ancestor
		       // and the ancestor is homozygous then divide the prob by 2
		       // OLD                       if (ind1 != fat && ind2 != fat && (listlen(fat->marker)==1))
                       if (ind1 != fat && ind2 != fat)
		       {
			 // OLD                         mkr  = markernumber(fat, 1);

			 mkr = markernumber(fat, smartarray[fat->tmpint1-1]+1);

                         if (mkr->allele1 == mkr->allele2)
                           tempres /= 2;
		       }
		     }



		     // Have now calculated the path to the father
		     // Starting on doing the same for the mother

                     mkr  = markernumber(chain3->ind, smartarray[chain3->ind->tmpint1-1]+1);
                     mkr2 = markernumber(chain4->ind, smartarray[chain4->ind->tmpint1-1]+1);

                     // Verifying that the two alleles examined are IBS
                     j = 10*(savetime[chain3->ind->tmpint1][3])+(savetime[chain4->ind->tmpint1][4]);

		     switch(j)
	             {
		       case  0: // These first cases should only occur, when
		       case  1: // one of the two individuals is the common ancestor
		       case  2: 
		       case 10:
		       case 20: break;
		       case 11: // Both from mothers 
			       if (mkr->allele1 != mkr2->allele1) tempres2 = 0.0; break;
		       case 12: // Moth Fath
			       if (mkr->allele1 != mkr2->allele2) tempres2 = 0.0; break;
		       case 21: // Fath Moth
			       if (mkr->allele2 != mkr2->allele1) tempres2 = 0.0; break;
		       case 22: // Both from fathers
			       if (mkr->allele2 != mkr2->allele2) tempres2 = 0.0; break;
		     }

 		     // If they're IBS then continue, otherwise they're IBD=0
                     if (tempres2 > 0.0)
		     {

                       for (id = chain3; id && id->next; id = id->next)
                       {  
                         mkr = markernumber(id->ind, smartarray[id->ind->tmpint1-1]+1);
			 // fmkr is not necessarily the father. 
			 // Just the next person in the list
                         fmkr = markernumber(id->next->ind, smartarray[id->next->ind->tmpint1-1]+1);

			 // If that person is homozygote the divide the result by 2
			 // For all non-start and non-end persons in the chain
                         if (id->next->next && fmkr->allele1 == fmkr->allele2)
                           tempres2 /= 2;

                         if (id->next->ind==mot)
			 {
                           break;
			 }

                         if (savetime[id->ind->tmpint1][3]==1)
  		         {
			   // If from fathers father and not alike
                           if ((savetime[id->next->ind->tmpint1][3] ==1 && mkr->allele1 != fmkr->allele1) ||
                               (savetime[id->next->ind->tmpint1][3] ==2 && mkr->allele1 != fmkr->allele2))
                             tempres2 = 0.0;
			 }
			 else
		         {
                           if ((savetime[id->next->ind->tmpint1][3] ==1 && mkr->allele2 != fmkr->allele1) ||
                               (savetime[id->next->ind->tmpint1][3] ==2 && mkr->allele2 != fmkr->allele2))
                             tempres2 = 0.0;
			 }
		       }

                       for (id = chain4; id && id->next; id = id->next)
                       {  
                         mkr = markernumber(id->ind, smartarray[id->ind->tmpint1-1]+1);
			 // fmkr is not necessarily the father. 
			 // Just the next person in the list
                         fmkr = markernumber(id->next->ind, smartarray[id->next->ind->tmpint1-1]+1);

			 // If that person is homozygote the divide the result by 2
                         if (id->next->next && fmkr->allele1 == fmkr->allele2)
                           tempres2 /= 2;

                         if (id->next->ind==mot)
			 {
                           // share1->element[ind1->tmpint1][ind2->tmpint1] += tempres;
                           break;
			 }

                         if (savetime[id->ind->tmpint1][4]==1)
  		         {
			   // If from fathers father and not alike
                           if ((savetime[id->next->ind->tmpint1][4] ==1 && mkr->allele1 != fmkr->allele1) ||
                               (savetime[id->next->ind->tmpint1][4] ==2 && mkr->allele1 != fmkr->allele2))
                             tempres2 = 0.0;
			 }
			 else
		         {
                           if ((savetime[id->next->ind->tmpint1][4] ==1 && mkr->allele2 != fmkr->allele1) ||
                               (savetime[id->next->ind->tmpint1][4] ==2 && mkr->allele2 != fmkr->allele2))
                             tempres2 = 0.0;
			 }
		       }
     

		       // If none of the two persons are the common ancestor
		       // and the ancestor is homozygous then divide the prob by 2
		       //  OLD                     if (ind1 != mot && ind2 != mot && (listlen(mot->marker)==1))
                       if (ind1 != mot && ind2 != mot)
		       {
			 //   OLD                      mkr  = markernumber(mot, 1);

			 mkr = markernumber(mot, smartarray[mot->tmpint1-1]+1);

                         if (mkr->allele1 == mkr->allele2)
                           tempres2 /= 2;
		       }
		     }

		     if (tempres>0 && tempres2>0)
                       share2->element[ind1->tmpint1][ind2->tmpint1] += (tempres+founderprob*tempres2)/2;
		     else
		     {
		       if (tempres>0)
                         share1->element[ind1->tmpint1][ind2->tmpint1] += tempres;
		       else if (tempres2>0)
                         share1->element[ind1->tmpint1][ind2->tmpint1] += (founderprob*tempres2);
		     }

		     // The following IF is just for debugging
		     /*		     if (tempres>0 || tempres2 > 0)
		     {
		       mkr = markernumber(fat, smartarray[fat->tmpint1-1]+1);
		       mkr2 = markernumber(mot, smartarray[mot->tmpint1-1]+1);
		       mmkr = markernumber(ind1, smartarray[ind1->tmpint1-1]+1);
		       fmkr = markernumber(ind2, smartarray[ind2->tmpint1-1]+1);
		       printf("Test:  Share2: %f (total: %f)  Individuals: (%d <-> %d) | Fat: %d %d | Mot: %d %d | %d %d   %d %d   |  %f\n", 
                     share2->element[ind1->tmpint1][ind2->tmpint1], 
		     tmpmat->element[ind1->tmpint1][ind2->tmpint1],  ind1->tmpint1, ind2->tmpint1, mkr->allele1, mkr->allele2,
mkr2->allele1, mkr2->allele2,
mmkr->allele1, mmkr->allele2,
fmkr->allele1, fmkr->allele2, share2->element[ind1->tmpint1][ind2->tmpint1]/tmpmat->element[ind1->tmpint1][ind2->tmpint1]

);
		     }
		     */

                     // Should now update the smartarray
		     update=0;
		     j=0; // First element
		     while (!update)
		     {
		       smartarray[j]++;
		       if (smartarray[j]==maxarray[j])
		       {
			 smartarray[j]=0;
			 if (i<(elements-1))
			   j++;
			 else
			   update=1;
		       }
		       else
			 update=1;
		     }
		   }

	           freelist(chain1);
	           freelist(chain2);
	           freelist(chain3);
	           freelist(chain4);
		   
		   chain1 = NULL;
		   chain2 = NULL;
		   chain3 = NULL;
		   chain4 = NULL;

	           break;
	}
    }
  }

  // Returns the IBD sharing matrix

  for (i = 1; i<=tmpmat->rows; i++)
    for (j = i; j<=tmpmat->cols; j++)
      {
        if (tmpmat->element[i][j]>0)
	{
          resmat->element[i][j] = (share2->element[i][j] + 
               0.5*share1->element[i][j]) / tmpmat->element[i][j];
          resmat->element[j][i] = resmat->element[i][j];
	}
      }

  MtxDel(tmpmat);
  MtxDel(share1);
  MtxDel(share2);

  return (resmat);
}


// Calculates the complete PI Hat matrix
OLDMATRIX *CalcExactPIHat(individual *indlist, int markernum)
{
  int i, j, nobs;
  individual *ind, *ind2, *ind3, *genolist;
  OLDMATRIX *res, *tempres;

  nobs = listlen(indlist);
  CopyAlleleFreq(markernum);

  for (ind = indlist, i=1; ind; ind = ind->next, i++)
    ind->tmpint4 = -i;

  res = MtxNew(nobs, nobs);


  for (ind = indlist; ind; ind = ind->next)
  {
    // If not used before
    if (ind->tmpint4<0)
    {
      // Gets the allowed genotypes
      genolist = GenotypeElimination(ind->pedigree, markernum);

      // Calculates the IBDmatrix
      tempres = PedigreeIBD(genolist); 

      // Fixes the persons
      for (ind2 = genolist; ind2; ind2 = ind2->next)
      {
        ind3 = FindIndListMember(indlist,ind2->id);
        ind3->tmpint4 = -ind3->tmpint4;
      }

      for (ind2 = genolist, i=1; ind2; ind2 = ind2->next, i++)
      {
        for (ind3 = ind2, j=i; ind3; ind3 = ind3->next, j++)
	{
          res->element[abs(ind2->tmpint4)][abs(ind3->tmpint4)] = tempres->element[i][j];
          res->element[abs(ind3->tmpint4)][abs(ind2->tmpint4)] = tempres->element[j][i];
	}
      }
    
      MtxDel(tempres);
      FreeIndividualList(genolist);

    }
  }
  
  return(res);
}


//
// Calculates the likelihood of a pedigree with 
// ordered, non-missing genotyped.
// This function does no checking at all !
// 
double CalcLikelihoodOfPedigree(individual *indlist)
{
  double res;
  individual *ind;
  markerlist *mkr, *fmkr, *mmkr;

  res = 1.0;

  printf("GOGO\n");

  for (ind = indlist; ind; ind = ind->next)
  {
    printf("  Ind: %s\n", ind->id);
    mkr = markernumber(ind,1);
    printf("ok\n");
    if (founder(ind))
    {
      res *= allelefreq[mkr->allele1]*allelefreq[mkr->allele2];
    }
    else
    {
      fmkr = markernumber(ind->father,1);
      mmkr = markernumber(ind->mother,1);
  
      if (fmkr->allele1 != fmkr->allele2)
        res /= 2;
      if (mmkr->allele1 != mmkr->allele2)
        res /= 2;

    }
  }
  printf("GOGOGO\n");
  return res;
}



//
// Do MC estimation of pi-hat on the genotype eliminated pedigree given by
// indlist for the relevant marker
// Should only be called when marker freq is ok
//
OLDMATRIX *CalcMCPIHatForPedigree(individual *indlist, int markernum)
{
  individual *ind, *workingdata;
  markerlist *mkr, *mkr3;
  int i, ii, nOffspring;
  double totalprob, combinations, u;
  double weight, totalweight, parweight, modweight, totalmass;
  double impweight;
  OLDMATRIX *res, *tempres;

  if (options->Index[O_XLINKED]) {
    printf("ERROR: The importance sampling/MC IBD estimation does not work for x-linked loci\n");
    return 0;
  }

  // Used for checking
  int checked[MAXPERSONS];
  memset(checked, 0, MAXPERSONS*sizeof(int));

  res = MtxNew(listlen(indlist),listlen(indlist));

  CopyAlleleFreq(markernum);

  totalweight = 0;
  combinations = 0;

  printf("NOW HERE\n");

  for (i = 1; i <= MCSIM; i++)
  {
    printf("MC iteration %d\n", i);
    weight = 1.0;
    modweight = 1.0;
    workingdata = CopyIndividualList(indlist);

    for (ind = workingdata, ii=1; ind; ind = ind->next, ii++)
    {
      // Checked is set to 1
      // This is strictly neccesary 
      checked[ii] = 1;
    }

    printf("Part 1\n");
    // Pick a random possible ordered genotype for each founder
    // and remove the rest
    for (ind = workingdata, ii=1; ind; ind = ind->next, ii++)
    {
      if (founder(ind))
      {
	printf("Founder %s\n", ind->id);
	// Pick a random possibility if more than one exists
	if (listlen(ind->marker)>1)
	{      
	  // Calculate total probability
	  nOffspring = listlen(ind->offspring);
	  if (options->Index[O_XLINKED] && ind->sex == S_MALE) 
	    impweight = 1;
	  else
	    impweight = pow(2, nOffspring);

	  // Calculate the probability of the nuclear pedigree
	  // Xlinked: There is no randomness in male founders segregation of alleles. Female are like for autosomes

	  totalprob = 0.0;
	  for (mkr = ind->marker; mkr; mkr = mkr->next)
	  {
	    printf("  Marker %2d %2d\n", mkr->allele1, mkr->allele2);
#ifdef XLINKED
	  if (options->Index[O_XLINKED] && ind->sex == S_MALE) {
	      totalprob += allelefreq[mkr->allele2];
	  }
	  else {
	    // Homozygous
	    if (mkr->allele1 == mkr->allele2)
	      totalprob += allelefreq[mkr->allele1]*allelefreq[mkr->allele2];
	    else
	      totalprob += allelefreq[mkr->allele1]*allelefreq[mkr->allele2]/impweight;
	  }
#else
	    if (mkr->allele1 == mkr->allele2)
	      totalprob += allelefreq[mkr->allele1]*allelefreq[mkr->allele2];
	    else
	      totalprob += allelefreq[mkr->allele1]*allelefreq[mkr->allele2]/impweight;
#endif
	  }
	  totalmass = totalprob;

	  checked[ii] = listlen(ind->marker);

	  // Draw random variable
	  u = genunf(0, totalprob);

	  for (mkr = ind->marker; mkr; )
	  {
            mkr3 = mkr->next;

	    if (mkr->allele1 == mkr->allele2)
	      totalprob -= allelefreq[mkr->allele1]*allelefreq[mkr->allele2];
	    else
	      totalprob -= allelefreq[mkr->allele1]*allelefreq[mkr->allele2]/impweight;

	    if (u > totalprob)
	    {
	      u = -totalprob;
	      if (mkr->allele1 == mkr->allele2)
	        modweight *= allelefreq[mkr->allele1]*allelefreq[mkr->allele2]/totalmass;
	      else
	        modweight *= allelefreq[mkr->allele1]*allelefreq[mkr->allele2]/(impweight*totalmass);
	    }
	    else
              removelist(&ind->marker, mkr);

	    mkr = mkr3;
	  }

	  // Have now picked a genotype for a founder
	  // Reduce possible genotypes based on this result

	  combinations = ReduceGenotypes(workingdata);

	}
      }
    }
    printf("Part 2\n");
    // Do the same thing for each non-founder with both parents having
    // only a single ordered genotype
    // Keep doing it until there is exactly one genotype for each person
    while (combinations>1)
    {
      for (ind = workingdata, ii=1; ind; ind = ind->next, ii++)
      {        
        if (!founder(ind) && (listlen(ind->marker)>1) && listlen(ind->father->marker)==1 && listlen(ind->mother->marker)==1)
        {
          checked[ii] = listlen(ind->marker);

          if (listlen(ind->marker)>1)        
          {
	    #ifdef IMPSAMPLE

  	    totalprob = 0.0;
	    nOffspring = listlen(ind->offspring);
	    impweight = pow(2, nOffspring);

	    for (mkr = ind->marker; mkr; mkr = mkr->next)
	    {
	      if (mkr->allele1 == mkr->allele2)
		totalprob += impweight;
	      else
		totalprob += 1;
	    }
	    totalmass = totalprob;

	    // XXX Should probably typecast something else than an int
	    u = genunf(0, totalprob);

  	    for (mkr = ind->marker; mkr; )
	    {
              mkr3 = mkr->next;
	      if (mkr->allele1 == mkr->allele2)
		totalprob -= impweight;
	      else
		totalprob -= 1.0;

	      if (u > totalprob)
	      {
	        u = -totalprob;

 	        if (mkr->allele1 == mkr->allele2)
		  modweight *= impweight / totalmass;
	        else
		  modweight *= 1.0 / totalmass;		
	      }
	      else
                removelist(&ind->marker, mkr);

	      mkr = mkr3;
	    }


	    #else

	    // The genotypes for a non-founder are equally probable
	    // if both parents have only one genotype
	    k = (int)ignuin(1,listlen(ind->marker));

	    // modweight *= 1.0/ (double) listlen(ind->marker);

	    for (mkr = ind->marker, j=1; mkr; j++)
	    {
	      mkr3 = mkr->next;

	      if (k != j)
		removelist(&ind->marker, mkr);

	      mkr = mkr3;
	    }

	    #endif
	    combinations = ReduceGenotypes(workingdata);
	  }
	}
      }
    }
    printf("Part 3\n");
    //    weight = modweight*CalcLikelihoodOfPedigree(workingdata);
    weight = CalcLikelihoodOfPedigree(workingdata);
    printf("Part 4\n");
    for (ind = workingdata, ii=1; ind; ind = ind->next, ii++)
    {
      if (founder(ind))
	{
	  //          modweight /= (double) checked[ii];
	}
      else  // It is a non-founder
      {
        parweight = 1.0;

	//	if (heterozygous(ind->father,markernum))
	//          parweight *=2.0;
	//        if (heterozygous(ind->mother,markernum))
	//          parweight *=2.0;

	if (heterozygous(ind->father,1))
          parweight *=2.0;
        if (heterozygous(ind->mother,1))
          parweight *=2.0;

	//   modweight *= (double)checked[ii]/parweight;
	//	modweight /=(double) parweight;
	//       	 weight *=1.0 / parweight;
	//	modweight *= 1.0 / (double) checked[ii];
      }
    }

    printf("Part 5\n");
    weight /= modweight;

    totalweight += weight;

    // Should now calculate the exact IBD distribution
    // From the reduced dataset
    tempres = PedigreeIBD(workingdata);

    MtxScale(tempres, weight);
  
    // Add to result
    MtxAdd(res,tempres,res);   

    MtxDel(tempres);  

    FreeIndividualList(workingdata);
  }
  printf("AND NOW HERE\n");

  MtxScale(res, (double) 1.0/res->element[1][1] );

  return (res);
}



// Calculates the complete PI Hat matrix by MC
OLDMATRIX *CalcMCPIHat(individual *indlist, int markernum)
{
  int i, j, nobs;
  individual *ind, *ind2, *ind3, *genolist;
  OLDMATRIX *res, *tempres;

  nobs = listlen(indlist);
  CopyAlleleFreq(markernum);

  for (ind = indlist, i=1; ind; ind = ind->next, i++)
    ind->tmpint4 = -i;

  res = MtxNew(nobs, nobs);


  for (ind = indlist; ind; ind = ind->next)
  {
    // If not used before
    if (ind->tmpint4<0)
    {
      // Gets the allowed genotypes
      genolist = GenotypeElimination(ind->pedigree, markernum);

      // Calculates the IBDmatrix
      tempres = CalcMCPIHatForPedigree(genolist, markernum); 

      // Fixes the persons
      for (ind2 = genolist; ind2; ind2 = ind2->next)
      {
        ind3 = FindIndListMember(indlist,ind2->id);
        ind3->tmpint4 = -ind3->tmpint4;
      }

      for (ind2 = genolist, i=1; ind2; ind2 = ind2->next, i++)
      {
        for (ind3 = ind2, j=i; ind3; ind3 = ind3->next, j++)
	{
          res->element[abs(ind2->tmpint4)][abs(ind3->tmpint4)] = tempres->element[i][j];
          res->element[abs(ind3->tmpint4)][abs(ind2->tmpint4)] = tempres->element[j][i];
	}
      }
    
      MtxDel(tempres);
      FreeIndividualList(genolist);

    }
  }
  
  return(res);
}

// Calculates the complete PI Hat matrix by the smartest method
OLDMATRIX *CalcPIHat(individual *indlist, int markernum)
{
  int i, j, nobs, nAncestors, code;
  individual *ind, *ind2, *ind3, *genolist, *father, *mother;
  OLDMATRIX *res, *tempres;
  double combinations;
  IDlist *chain, *id;
#ifndef NEWXLINKED
  FILE *fp;
  double sumfreq;
#endif

  nobs = listlen(indlist);
  CopyAlleleFreq(markernum);

#ifdef BYPEDIGREE
  printf("WARNING: USES AN UNTESTED VERSION OF IBD CALCULATION. USE WITH CARE\n");
  sprintf(buf, "marker%d.ibd", markernum);

  if (!(fp = fopen (buf,"w"))) {
    printf("ERROR: Could not open marker4.ibd for writing. Exiting.\n");    
    exit(1);
  }
#endif



#ifndef NEWXLINKED
  // Only include if I haven't done the fix correctly
#ifdef XLINKED
  // Need to fix the allele frequencies as the pseudo-allele for the males
  // does not really exist. This is done (tsk ... tsk) by adjusting the
  // remaining allele frequencies and given the pseudomarker a frequency 
  // of 0.001

  if (options->Index[O_XLINKED]) {
    sumfreq = 0;
    for (i = 1; i<=NumberOfAlleles(markernum); i++) {
      sumfreq += allelefreq[i];
    }
    sumfreq /= .999;
    for (i = 1; i<=NumberOfAlleles(markernum); i++) {
      allelefreq[i] *= sumfreq;
    }
    allelefreq[NumberOfAlleles(markernum)+1] = .001;    
  }
#endif
#endif

  for (ind = indlist, i=1; ind; ind = ind->next, i++)
    ind->tmpint4 = -i;

#ifdef BYPEDIGREE
  res = MtxNew(1,1);
#else
  res = MtxNew(nobs, nobs);
#endif

  for (ind = indlist; ind; ind = ind->next) {
    // If not used before
    if (ind->tmpint4<0) {
      // Gets the allowed genotypes
      genolist = GenotypeElimination(ind->pedigree, markernum);

      combinations = 0;
      for (ind2 = genolist; ind2; ind2 = ind2->next)
	combinations += log(listlen(ind2->marker));

#ifdef debug
      printf("Calculating IBD for pedigree %s ", ind->pedigree);
#endif

      // Calculates the IBDmatrix
      if (combinations < log(MAXEXACT))  {
#ifdef debug
	printf("(exact)\n");
#endif
	tempres = PedigreeIBD(genolist); 
      }
      else  {
#ifdef debug
	printf("(MC)\n");
#endif
	tempres = CalcMCPIHatForPedigree(genolist, markernum); 
      }

      // Fixes the persons
      // I.e. go through the pedigree and find (in indlist) the id's
      // of all the individuals we've just handled in the pedigree.
      // Change the sign for each of them
      for (ind2 = genolist; ind2; ind2 = ind2->next) {
        ind3 = FindIndListMember(indlist,ind2->id);
        ind3->tmpint4 = -ind3->tmpint4;
      }

      for (ind2 = genolist, i=1; ind2; ind2 = ind2->next, i++) {
        for (ind3 = ind2, j=i; ind3; ind3 = ind3->next, j++) {
#ifdef BYPEDIGREE 
	  if (tempres->element[i][j]>0)
	    fprintf(fp, "%5d %5d  %9.7f\n", abs(ind2->tmpint4), abs(ind3->tmpint4), tempres->element[i][j]);
#else 
          res->element[abs(ind2->tmpint4)][abs(ind3->tmpint4)] = tempres->element[i][j];
          res->element[abs(ind3->tmpint4)][abs(ind2->tmpint4)] = tempres->element[j][i];
#endif
	}
      }
    
      MtxDel(tempres);
      FreeIndividualList(genolist);
    }
  }

#ifdef BYPEDIGREE
  fclose(fp);
  MtxDel(res);
  res = NULL;
#endif

#ifdef XLINKED
  if (options->Index[O_XLINKED]) {
    // This part fixed the calculations for the different combined pairs
    for (ind = indlist, i=1; ind; ind = ind->next, i++) {
      for (ind2 = ind, j=i; ind2; ind2 = ind2->next, j++) {
	if (ind==ind2)
	  continue;
	// Male-male
	// Remove .5 and multiply by two if the two individuals
	// have an ancestral father in common.
	if (ind->sex == S_MALE && ind2->sex == S_MALE) {	  
	  nAncestors = FindCommonAncestor(ind, ind2, &father, &mother);
	  if (father != NULL) {   	    
	    // Have a paternal ancestor in common.
	    // Go through the linked chain from each individual to his ancestor.
	    // If two males occur in a row the discard the effect of the paternal allele
	    // This is easier computed as follows:
	    // If a single female exists on the chain between the two males then
	    // the IBD-score will be correct as she makes sure that the pseudo-allele
	    // is not segregated on. Conversely, if the chain is purely male, each male 
	    // will have the pseudo allele and the IBD-score will be .5 too high
	    code = 0;

	    chain = LinkToAncestor(ind, father);	    
	    for (id = chain; id; id = id->next) {
	      if (id->ind->sex==S_FEMALE) {
		code = 1;
		break;
	      }
	    }
	    // Check the other chain
	    if (code == 0) {
	      chain = LinkToAncestor(ind2, father);
	      for (id = chain; id; id = id->next) {
		if (id->ind->sex==S_FEMALE) {
		  code = 1;
		  break;
		}
	      }
	    }
	    if (code == 0) {
	      // Both chains are only male, so the IBD score is off by .5
	      res->element[i][j] -= .5;
	      if (res->element[i][j]<0) 
		res->element[i][j] = 0;
	      res->element[j][i] = res->element[i][j];	     
	    }
	  }
	  res->element[i][j] *= 2;	    
	  res->element[j][i] *= 2;	    
	}
	else if (ind->sex == S_FEMALE && ind2->sex == S_FEMALE) {
	  // Female-female pair
	  // Do nothing
	}
	else {
	  // Male-female pair
	  // Multiply by two
	  res->element[i][j] *= 2;	    
	  res->element[j][i] *= 2;	    
	}
      }
    }
  }

#endif
  return(res);
}





int RelationType(individual *ind, individual *ind2)
{
  int nancest, l1, l2, code;
  individual *fat, *mot;
  IDlist *id;

  nancest = FindCommonAncestor(ind, ind2, &fat, &mot);
  l1 = 0;
  l2 = 0;

  if (!nancest)
    return REL_NONE;  

  // If the two persons have any common ancestors
  if (fat)
  {
    id = LinkToAncestor(ind, fat);    
    l1 = listlen(id)-1;
    freelist(id);
    id = NULL;
    id = LinkToAncestor(ind2, fat);    
    l2 = listlen(id)-1;
    freelist(id);
    id = NULL;
  }
  else if (mot)
  {
    id = LinkToAncestor(ind, mot);    
    l1 = listlen(id)-1;
    freelist(id);
    id = NULL;
    id = LinkToAncestor(ind2, mot);
    l2 = listlen(id)-1;
    freelist(id);
    id = NULL;
  }

  code = min(l1,l2)*MAXGEN+max(l1,l2);

  switch (code)
  {
    case    0 : return REL_ID; break;
    case    1 : return REL_PARENT; break;
    case    2 : return REL_GRAND_PARENT; break;
    case    3 : return REL_GREAT_GRAND_PARENT; break;
    case    4 : return REL_GREAT_GREAT_GRAND_PARENT; break;
    case    5 : return REL_GREAT_GREAT_GREAT_GRAND_PARENT; break;
    case   11 : if (nancest==2) return REL_FULL_SIBS; else return REL_HALF_SIBS;  break;
    case   12 : if (nancest==2) return REL_AVUNCULAR; else return REL_HALF_AVUNCULAR;  break;
    case   13 : if (nancest==2) return REL_GRAND_AVUNCULAR; else return REL_HALF_GRAND_AVUNCULAR;  break;
    case   14 : if (nancest==2) return REL_GREAT_GRAND_AVUNCULAR; else return REL_HALF_GREAT_GRAND_AVUNCULAR;  break;
    case   22 : if (nancest==2) return REL_FIRST_COUSIN; else return REL_HALF_FIRST_COUSIN;  break;
    case   23 : if (nancest==2) return REL_FIRST_COUSIN_ONCE_REMOVED; else return REL_HALF_FIRST_COUSIN_ONCE_REMOVED;  break;    
    case   24 : if (nancest==2) return REL_FIRST_COUSIN_TWICE_REMOVED; else return REL_HALF_FIRST_COUSIN_TWICE_REMOVED;  break;
    case   33 : if (nancest==2) return REL_SECOND_COUSIN; else return REL_HALF_SECOND_COUSIN;  break;
    case   34 : if (nancest==2) return REL_SECOND_COUSIN_ONCE_REMOVED; else return REL_HALF_SECOND_COUSIN_ONCE_REMOVED;  break;
    case   44 : if (nancest==2) return REL_THIRD_COUSIN; else return REL_HALF_THIRD_COUSIN;  break;
    default : return -1;
  }
}


// Creates a matrix with the relationships between all persons in the dataset
OLDMATRIX *CalcDataRelationship(individual *indlist)
{
  int i, code;
  individual *ind, *ind2;
  OLDMATRIX *res;

  for (ind = indlist, i=1; ind; ind = ind->next, i++)
    ind->tmpint4 = i;

  res = MtxNew(listlen(indlist), listlen(indlist));

  code = 1;
  for (ind = indlist; ind; ind = ind->next)
  {
    for (ind2 = ind; ind2; ind2 = ind2->next)
    {
      res->element[ind->tmpint4][ind2->tmpint4] = RelationType(ind,ind2);
      res->element[ind2->tmpint4][ind->tmpint4] = res->element[ind->tmpint4][ind2->tmpint4];
      if (res->element[ind->tmpint4][ind2->tmpint4]==-1)
      {
        printf("ERROR: Relationship between %s and %s is not known. Please modify RelationType function\n", ind->id, ind2->id);
        sprintf(buf, "Relationship between %s and %s is not known. Please modify RelationType function\n", ind->id, ind2->id);
        WriteErrorMsg(buf);
        code = 0;
      }
    }
  }

  if (!code)
    printf("WARNING: Not all family relationships could be computed\n");
  
  return(res);
}

// Returns the 
// Table 3 in Almasy and Blangero
double CorrCoeffIBD(double theta, int relationtype)
{
  switch (relationtype)
  {
    case REL_NONE: return (0); break;
    case REL_ID: return (0); break; // XXX Hmm check these simple cases
    case REL_PARENT: return (0); break;
    case REL_GRAND_PARENT: return (1-2*theta); break;
    case REL_GREAT_GRAND_PARENT: return (1-8.0/3*theta+4.0/3*theta*theta); break;
    case REL_GREAT_GREAT_GRAND_PARENT: return (1-24.0/7*(theta-theta*theta) - 8.0/7*pow(theta,3)); break;
    case REL_FULL_SIBS: 
    case REL_HALF_SIBS: return (1-4*theta+4*theta*theta); break;
    case REL_AVUNCULAR: return (1-5*theta+8*theta*theta-4*theta*theta*theta); break;
    case REL_HALF_AVUNCULAR: return (1-4*theta+16.0/3*theta*theta-8.0/3*theta*theta*theta); break;
    case REL_FIRST_COUSIN: return (1-16.0/3*theta+10*pow(theta,2)-8*pow(theta,3)+8.0/3*pow(theta,4)); break;
    case REL_HALF_FIRST_COUSIN: return (1-32.0/7*theta+8*pow(theta,2)-48.0/7*pow(theta,3)+16.0/7*pow(theta,4)); break;
  }
  printf("WARNING: Relationship (%s) not found in CorrCoeffIBD\n", RelationshipNames[relationtype]);
  return -1;
}

// The average estimate of pi for a given relationship
// Taken from almasy and Blangero table 3
double meanpi(int relationtype)
{
  switch (relationtype)
  {
    case REL_NONE: return (0); break;
    case REL_ID: return (1); break; 
    case REL_PARENT: return (.5); break;
    case REL_GRAND_PARENT: return (.25); break;
    case REL_GREAT_GRAND_PARENT: return (.125); break;
    case REL_GREAT_GREAT_GRAND_PARENT: return (0.0625); break;
    case REL_GREAT_GREAT_GREAT_GRAND_PARENT: return (0.03125); break;
    case REL_FULL_SIBS: return .5; break;
    case REL_HALF_SIBS: return (.25); break;
    case REL_AVUNCULAR: return (.25); break;
    case REL_HALF_AVUNCULAR: return (.125); break;
    case REL_FIRST_COUSIN: return (.125); break;
    case REL_HALF_FIRST_COUSIN: return (.0625); break;
    case REL_GRAND_AVUNCULAR: return (.125); break;
    case REL_HALF_GRAND_AVUNCULAR: return (.0625); break;
    case REL_GREAT_GRAND_AVUNCULAR: return (.0625); break;
    case REL_HALF_GREAT_GRAND_AVUNCULAR: return (.03125); break;
    case REL_FIRST_COUSIN_ONCE_REMOVED: return (.0625); break;
    case REL_HALF_FIRST_COUSIN_ONCE_REMOVED: return (.03125); break;
    case REL_FIRST_COUSIN_TWICE_REMOVED: return (.03125); break;
    case REL_HALF_FIRST_COUSIN_TWICE_REMOVED: return (.015625); break;
    case REL_SECOND_COUSIN: return (.03125); break;
    case REL_HALF_SECOND_COUSIN: return (.015625); break;
    case REL_SECOND_COUSIN_ONCE_REMOVED: return (1.0/64.0); break;
    case REL_HALF_SECOND_COUSIN_ONCE_REMOVED: return (1.0/128.0); break;
    case REL_THIRD_COUSIN: return (1.0/128.0); break;
    case REL_HALF_THIRD_COUSIN: return (1.0/256.0); break;
  }
  printf("WARNING: Relationship (%s) not found in meanpi\n", RelationshipNames[relationtype]);
  return -1;
}


// The average variance of pi for a given relationship
// Taken from almasy and Blangero table 3
double varpi(int relationtype)
{
  switch (relationtype)
  {
    case REL_NONE:
    case REL_ID: 
    case REL_PARENT: return (0); break;
    case REL_FULL_SIBS: return (1.0/8); break;

    case REL_HALF_SIBS: 
    case REL_AVUNCULAR: 
    case REL_GRAND_PARENT: return (1.0/16); break;

    case REL_FIRST_COUSIN: 
    case REL_HALF_AVUNCULAR: 
    case REL_GREAT_GRAND_PARENT: return (3.0/64); break;

    case REL_GREAT_GREAT_GRAND_PARENT: return (7.0/256); break;
    case REL_HALF_FIRST_COUSIN: return (7.0/256); break;
  }
  printf("WARNING: Relationship (%s) not found in varpi\n", RelationshipNames[relationtype]);
  return -1;
}


// Undefines debug for this file onwards
#undef debug

/*******************************
 *
 * Returns E(\pi) for x-linked markers for two individuals
 *
 *******************************/

double xlinkedphi2(individual *ind, individual *ind2) 
{
  // Return 0 if not same pedigree
  if (strcmpl(ind->pedigree, ind2->pedigree))
    return 0;

  // Return 1 if same individula
  if (!strcmpl(ind->id, ind2->id))
    return 1;

  // Check sexes of individuals
  if (ind->sex == S_MALE && ind2->sex==S_MALE) {
    // Male-male
    switch (RelationType(ind,ind2)) {
    case REL_NONE: return 0;
    case REL_FULL_SIBS: return .5;
    case REL_PARENT: return 0;
    case REL_HALF_SIBS: if (ind->father == ind2->father) {
                          // Paternal HS
                          return 0;
                        }
                        else
			  return .5;
    default: printf("ERROR: Unavailable relationship (%s) for male-male pairs\n", RelationshipNames[RelationType(ind,ind2)]);
      return -1;
    }
  }
  else if (ind->sex == S_FEMALE && ind2->sex==S_FEMALE) {
    switch (RelationType(ind,ind2)) {
    case REL_NONE: return 0;
    case REL_PARENT: return .5;
    case REL_FULL_SIBS: return .75;
    case REL_HALF_SIBS: if (ind->father == ind2->father) {
                          // Paternal HS
                          return .5;
                        }
                        else
			  return .25;
    default: printf("ERROR: Unavailable relationship (%s) for female-female pairs\n", RelationshipNames[RelationType(ind,ind2)]);
      return -1;
    }
  }
  else if ((ind->sex == S_MALE && ind2->sex==S_FEMALE) || (ind->sex == S_FEMALE && ind2->sex==S_MALE)) {
    switch (RelationType(ind,ind2)) {
    case REL_NONE: return 0;
    case REL_PARENT: return 1;
    case REL_FULL_SIBS: return .5;
    case REL_HALF_SIBS: if (ind->father == ind2->father) {
                          // Paternal HS
                          return 0;
                        }
                        else
			  return .5;
    default: printf("ERROR: Unavailable relationship (%s) for male-female pairs\n", RelationshipNames[RelationType(ind,ind2)]);
      return -1;
    }
  }
  else {
    printf("ERROR: Illegal sexes when computing X-linked kinship matrix for individuals %s and %s\n", ind->id, ind2->id);
  }

  return -1;
}

#undef XLINKED

/*
double rho(int rela, double theta, int sex1, int sex2) {

  if (sex1==S_MALE && sex2==S_MALE) {
    switch (relationtype) {
    case REL_NONE: return (0); break;
    case REL_ID: return (1); break;       
    }
    printf("WARNING: Relationship (%s) not found in meanpi\n", RelationshipNames[relationtype]);
  } 
  else if (sex1==S_FEMALE && sex2==S_FEMALE) {
    switch (relationtype) {
    case REL_NONE: return (0); break;
    case REL_ID: return (1); break;       
    }
    printf("WARNING: Relationship (%s) not found in meanpi\n", RelationshipNames[relationtype]);
  }
  printf("WARNING: Relationship (%s) not found in meanpi\n", RelationshipNames[relationtype]);
  return -1;
}


int CalculateMultipointIBD(individual *indlist, double distance)
{
    OLDMATRIX *relationship;

    relationship = CalcDataRelationship(individuals);
  


    MtxDel(relationship);
}
*/

int IBDforInheritanceVector(char *pedigree, unsigned long int vector)
{
  individual *ind, *indlist;
  markerlist *marker;
  int nBits = 0;
  int all = 1;
  int i, j, k;
  FreqList *fl;
  OLDMATRIX *tempres;

  for (ind = individuals; ind; ind = ind->next) {
    if (!founder(ind)) {
      if (options->Index[O_XLINKED]) {
	nBits++;
      }
      else 
	nBits += 2;
    }
  }
  if (log(vector+.1)>nBits*log(2)) 
    return 0;

  indlist = SelectPedigree(pedigree, individuals);

  for (ind = indlist; ind; ind = ind->next) {
    //    printf("Individual: %s\n", ind->id);
    
    freelist(ind->marker);
    ind->marker = NULL;

    marker = cmalloc (sizeof (markerlist));
    memset(marker, 0, sizeof(markerlist));

    if (founder(ind)) {
      if (options->Index[O_XLINKED] && ind->sex == S_MALE) {
	marker->allele1 = 65;
	marker->allele2 = all;
	all++;
      }
      else {
	marker->allele1 = all;
	all++;
	marker->allele2 = all;
	all++;
      }
    }
    else {
      if (vector & 1) 
	marker->allele2 = ind->mother->marker->allele1;
      else 
	marker->allele2 = ind->mother->marker->allele2;
      vector >>= 1;


      if (options->Index[O_XLINKED]) {
	if (ind->sex == S_FEMALE) {
	  marker->allele1 = ind->father->marker->allele2;
	  //	  vector >>= 1;
	}
	else {
	  marker->allele1 = 65;
	}
      }
      else {
	if (vector & 1) 
	  marker->allele1 = ind->father->marker->allele1;
	else 
	  marker->allele1 = ind->father->marker->allele2;
	vector >>= 1;
      }      
    }
    addlist(&ind->marker, marker);
  }  

  fl = FrequencyNumber(1);
  for (j=0; j<MAXALLELES; j++) {
    fl->frequency[j] = (double) 1.0 / (float)MAXALLELES;
    allelefreq[j] = fl->frequency[j];
  }


    for (ind = indlist, k=1; ind ; ind = ind->next, k++) 
      printf("%s: %d %d \n", ind->id, ind->marker->allele1, ind->marker->allele2);

  tempres = PedigreeIBD(indlist); 
  //  MtxPrint(tempres);
  i = 0; 
  j = 0;
  for (ind = individuals, k=1; ind && !j; ind = ind->next, k++) {    
    if (ind->marker->allele1 == 2) {
      if (i == 0)
	i = k;
      else 
	j = k;
    }
  }
  //  printf("Bits required     : %d \n", nBits);
  //  printf("Inheritance vector: %ld \n", oldvector);
  printf("Relationship %2d-%2d: %f \n", i, j, tempres->element[i][j]);

  MtxDel(tempres);

  FreeIndividualList(indlist);

  if (tempres->element[i][j]>0)
    return 1;
  return 0;
}





/*

int MIBDestimation(individual *indlist) {
  int rela, i, j, code;
  double loca, chromlength;
  OLDMATRIX *relationship;
  
  relationship = CalcDataRelationship(individuals);

  // Go through each possible relationship
  for (rela = REL_NONE+1; rela<REL_MAXREL; rela++) {
    code = 0;

    // Does this relationship exist?
    for (i = 1; i< relationship->rows && !code; i++) {
      for (j = i+1; j<= relationship->cols && !code; j++) {
	if (relationship->element[i][j]==rela)
	  code = 1;
      }
    }
    if (code) {
      printf("Skipping relationship %d\n", rela);
    }


    VC = MtxNew(numberofmarkers(), numberofmarkers());

    for (i=1; i<= numberofmarkers(); i++) {
      for (j=i; j<= numberofmarkers(); j++) { 
        if (i==j)
          VC->element[i][j] = 1;
        else
        {
          VC->element[i][j] = rho(rela, InterMarkerRecomb(i,j), relation)*
	                      BaseC->element[i][1]*BaseC->element[j][1]/varpi(relation);
          VC->element[j][i] = VC->element[i][j];
        }
	//  C->element[i][1] = CorrCoeffIBD(InverseMapFunction(fabs(MarkerDistance(location, i))/100), relation)*var;
      }
    }
    MtxInver(VC);



    for (loca = 0; ChromosomeLength()*100; loca +=2) {
      
    }





  }
  
  return 0;
}
*/
#endif
