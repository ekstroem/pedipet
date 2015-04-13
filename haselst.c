/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * Functions to do an analysis using Haseman-Elston method
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

#include <math.h>
#include "haselst.h"

// Corresponds to the g function in the SAGE manual
int HEpenetrance(int pheno1, int pheno2, int gen1, int gen2)
{
  if ((gen1==pheno1 && gen2==pheno2) || (gen1==pheno2 && gen2==pheno1))
    return 1;
  return 0;
}


double HEFSLikelihoodFamily(individual *ind, int markernum)
{
  int numberofalleles, i, j;
  individual *father, *mother;
  markerlist *marker, *mmkr, *fmkr;
  phenoset *fphenoset, *mphenoset, *pheno, *pheno2;
  double res, tempres;
  IDlist *sibs;

  // Calculating the number of alleles for this marker
  numberofalleles = NumberOfAlleles(markernum);

  father = ind->father;
  mother = ind->mother;
  fmkr = markernumber(father, markernum);
  mmkr = markernumber(mother, markernum);
  fphenoset = 0;
  mphenoset = 0;

  if (istypedformarker(father,markernum))
  {
    pheno = cmalloc (sizeof (phenoset));
    memset (pheno,0,sizeof (phenoset));
    pheno->allele1 = min(fmkr->allele1, fmkr->allele2);
    pheno->allele2 = max(fmkr->allele1, fmkr->allele2);

    addlist(&fphenoset,pheno);
  }
  else
  {
    for (i =1; i<=numberofalleles; i++)
    {
      for (j =i; j<=numberofalleles; j++)
      {
        pheno = cmalloc (sizeof (phenoset));
        memset (pheno,0,sizeof (phenoset));
        pheno->allele1 = i;
        pheno->allele2 = j;
        addlist(&fphenoset,pheno);
      }
    }
  }

  // Creating the mothers possible phenotype
  if (istypedformarker(mother,markernum))
  {
    pheno = cmalloc (sizeof (phenoset));
    memset (pheno,0,sizeof (phenoset));
    pheno->allele1 = min(mmkr->allele1, mmkr->allele2);
    pheno->allele2 = max(mmkr->allele1, mmkr->allele2);

    addlist(&mphenoset,pheno);
  }
  else
  {
    for (i =1; i<=numberofalleles; i++)
    {
      for (j =i; j<=numberofalleles; j++)
      {
        pheno = cmalloc (sizeof (phenoset));
        memset (pheno,0,sizeof (phenoset));
        pheno->allele1 = i;
        pheno->allele2 = j;
        addlist(&mphenoset,pheno);
      }
    }
  }

  // Should now calculate the family likelihood (1) in SAGE manual

  res = 0.0;
  for (pheno = fphenoset; pheno; pheno = pheno->next)
  {
    for (pheno2 = mphenoset; pheno2; pheno2 = pheno2->next)
    {
      tempres = allelefreq[pheno->allele1]*allelefreq[pheno->allele2]*
                allelefreq[pheno2->allele1]*allelefreq[pheno2->allele2];

      // If phenotype ab then multiply by two (ab and ba are equi. prop under HW)
      if (pheno->allele1 != pheno->allele2)
        tempres *= 2;
      if (pheno2->allele1 != pheno2->allele2)
        tempres *= 2;

      for (sibs = ind->sib; sibs; sibs = sibs->next)
      {
        // Checking that they are full sibs and that the sibling is typed
        if (fullsibs(ind,sibs->ind) && istypedformarker(sibs->ind,markernum))
        {
          marker = markernumber(sibs->ind,markernum);
          tempres = tempres * 0.25*(
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele1) +
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele2) +
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele1) +
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele2));
        }
      }
      marker = markernumber(ind,markernum);
      tempres = tempres * 0.25*(
      HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele1) +
      HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele2) +
      HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele1) +
      HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele2));
      res = res + tempres;
    }
  }
  freelist(fphenoset);
  freelist(mphenoset);
  return(res);
}

double HELikelihoodOneAllele(individual *ind, individual *ind2, int markernum)
{
  int numberofalleles, i, j;
  individual *father, *mother;
  markerlist *marker, *mmkr, *fmkr, *marker1, *marker2;
  phenoset *fphenoset, *mphenoset, *pheno, *pheno2;
  double tempres, res;
  IDlist *sibs;

  // Calculating the number of alleles for this marker
  numberofalleles = NumberOfAlleles(markernum);

  father = ind->father;
  mother = ind->mother;
  fmkr = markernumber(father, markernum);
  mmkr = markernumber(mother, markernum);

  fphenoset = 0;
  mphenoset = 0;

  if (istypedformarker(father,markernum))
  {
    pheno = cmalloc (sizeof (phenoset));
    memset (pheno,0,sizeof (phenoset));
    pheno->allele1 = min(fmkr->allele1, fmkr->allele2);
    pheno->allele2 = max(fmkr->allele1, fmkr->allele2);

    addlist(&fphenoset,pheno);
  }
  else
  {
    for (i =1; i<=numberofalleles; i++)
    {
      for (j =i; j<=numberofalleles; j++)
      {
        pheno = cmalloc (sizeof (phenoset));
        memset (pheno,0,sizeof (phenoset));
        pheno->allele1 = i;
        pheno->allele2 = j;
        addlist(&fphenoset,pheno);
      }
    }
  }

  // Creating the mothers possible phenotype
  if (istypedformarker(mother,markernum))
  {
    pheno = cmalloc (sizeof (phenoset));
    memset (pheno,0,sizeof (phenoset));
    pheno->allele1 = min(mmkr->allele1, mmkr->allele2);
    pheno->allele2 = max(mmkr->allele1, mmkr->allele2);

    addlist(&mphenoset,pheno);
  }
  else
  {
    for (i =1; i<=numberofalleles; i++)
    {
      for (j =i; j<=numberofalleles; j++)
      {
        pheno = cmalloc (sizeof (phenoset));
        memset (pheno,0,sizeof (phenoset));
        pheno->allele1 = i;
        pheno->allele2 = j;
        addlist(&mphenoset,pheno);
      }
    }
  }

  // Should now calculate the family likelihood (1) in SAGE manual

  marker1 = markernumber(ind,markernum);
  marker2 = markernumber(ind2,markernum);

  res = 0.0;
  for (pheno = fphenoset; pheno; pheno = pheno->next)
  {
    for (pheno2 = mphenoset; pheno2; pheno2 = pheno2->next)
    {
      tempres = allelefreq[pheno->allele1]*allelefreq[pheno->allele2] *
                allelefreq[pheno2->allele1]*allelefreq[pheno2->allele2];

      // If phenotype ab then multiply by two (ab and ba are equi. prop under HW)
      if (pheno->allele1 != pheno->allele2)
        tempres *= 2;
      if (pheno2->allele1 != pheno2->allele2)
        tempres *= 2;

      for (sibs = ind->sib; sibs; sibs = sibs->next)
      {
        // Should be full sibs and not the sibling
        if (fullsibs(ind,sibs->ind) && (sibs->ind != ind2) && istypedformarker(sibs->ind,markernum))
        {
          marker = markernumber(sibs->ind,markernum);
          tempres = tempres * 0.25*(
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele1) +
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele2) +
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele1) +
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele2));
        }
      }
          tempres = tempres *(
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele1, pheno2->allele1) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele1, pheno2->allele2) +
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele1, pheno2->allele2) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele1, pheno2->allele1) +
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele2, pheno2->allele2) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele2, pheno2->allele1) +
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele2, pheno2->allele1) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele2, pheno2->allele2) +
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele1, pheno2->allele1) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele2, pheno2->allele1) +
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele2, pheno2->allele1) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele1, pheno2->allele1) +
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele1, pheno2->allele2) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele2, pheno2->allele2) +
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele2, pheno2->allele2) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele1, pheno2->allele2))/16.0;
      res = res + tempres;
    }
  }
  freelist(fphenoset);
  freelist(mphenoset);

  return(res);
}

double HELikelihoodTwoAlleles(individual *ind, individual *ind2, int markernum)
{
  int numberofalleles, i, j;
  individual *father, *mother;
  markerlist *marker, *mmkr, *fmkr, *marker1, *marker2;
  phenoset *fphenoset, *mphenoset, *pheno, *pheno2;
  double res, tempres;
  IDlist *sibs;

  // Calculating the number of alleles for this marker
  numberofalleles = NumberOfAlleles(markernum);

  father = ind->father;
  mother = ind->mother;
  fmkr = markernumber(father, markernum);
  mmkr = markernumber(mother, markernum);

  fphenoset = 0;
  mphenoset = 0;

  if (istypedformarker(father,markernum))
  {
    pheno = cmalloc (sizeof (phenoset));
    memset (pheno,0,sizeof (phenoset));
    pheno->allele1 = min(fmkr->allele1, fmkr->allele2);
    pheno->allele2 = max(fmkr->allele1, fmkr->allele2);

    addlist(&fphenoset,pheno);
  }
  else
  {
    for (i =1; i<=numberofalleles; i++)
    {
      for (j =i; j<=numberofalleles; j++)
      {
        pheno = cmalloc (sizeof (phenoset));
        memset (pheno,0,sizeof (phenoset));
        pheno->allele1 = i;
        pheno->allele2 = j;
        addlist(&fphenoset,pheno);
      }
    }
  }

  // Creating the mothers possible phenotype
  if (istypedformarker(mother,markernum))
  {
    pheno = cmalloc (sizeof (phenoset));
    memset (pheno,0,sizeof (phenoset));
    pheno->allele1 = min(mmkr->allele1, mmkr->allele2);
    pheno->allele2 = max(mmkr->allele1, mmkr->allele2);

    addlist(&mphenoset,pheno);
  }
  else
  {
    for (i =1; i<=numberofalleles; i++)
    {
      for (j =i; j<=numberofalleles; j++)
      {
        pheno = cmalloc (sizeof (phenoset));
        memset (pheno,0,sizeof (phenoset));
        pheno->allele1 = i;
        pheno->allele2 = j;
        addlist(&mphenoset,pheno);
      }
    }
  }

  // Should now calculate the fanily likelihood (1) in SAGE manual

  marker1 = markernumber(ind,markernum);
  marker2 = markernumber(ind2,markernum);

  res = 0.0;
  for (pheno = fphenoset; pheno; pheno = pheno->next)
  {
    for (pheno2 = mphenoset; pheno2; pheno2 = pheno2->next)
    {
      tempres = allelefreq[pheno->allele1]*allelefreq[pheno->allele2]*
                allelefreq[pheno2->allele1]*allelefreq[pheno2->allele2];

      // If phenotype ab then multiply by two (ab and ba are equi. prop under HW)
      if (pheno->allele1 != pheno->allele2)
        tempres *= 2;
      if (pheno2->allele1 != pheno2->allele2)
        tempres *= 2;

      for (sibs = ind->sib; sibs; sibs = sibs->next)
      {
        // Should be full sibs and not the sibling
        if (fullsibs(ind,sibs->ind) && (sibs->ind != ind2) && istypedformarker(sibs->ind,markernum))
        {
          marker = markernumber(sibs->ind,markernum);
          tempres = tempres * 0.25*(
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele1) +
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele2) +
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele1) +
          HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele2));
        }
      }    
          tempres = tempres *(
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele1, pheno2->allele1) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele1, pheno2->allele1) +
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele2, pheno2->allele2) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele2, pheno2->allele2) +
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele1, pheno2->allele2) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele1, pheno2->allele2) +
          HEpenetrance(marker1->allele1, marker1->allele2, pheno->allele2, pheno2->allele1) *
          HEpenetrance(marker2->allele1, marker2->allele2, pheno->allele2, pheno2->allele1))/16.0;
      res = res + tempres;
    }
  }
  freelist(fphenoset);
  freelist(mphenoset);

  return(res);
}


double SingleMarkerIBD(individual *ind, individual *ind2, int markernum)
{
  double f1, f2, fam;
  fam = 0.0;

  fam = HEFSLikelihoodFamily(ind, markernum);
  f1  = HELikelihoodOneAllele(ind, ind2, markernum) / fam;
  f2  = HELikelihoodTwoAlleles(ind, ind2, markernum) / fam;

  k2  = f2;

  return (f2+0.5*f1);
}


//
// Rutines for exact calculation of IBD scores for half sibs
// From SAGE manual

double HEHS(individual *ind, individual *ind2, int markernum)
{
  int numberofalleles, i, j;
  individual *commonpar, *halfp1, *halfp2;
  markerlist *marker, *marker2, *commkr, *h1mkr, *h2mkr;
  phenoset *comps, *h1ps, *h2ps, *pheno, *pheno2, *pheno3;
  IDlist *sibs;
  double nom, denom, tempres;

  // Calculating the number of alleles for this marker
  numberofalleles = NumberOfAlleles(markernum);

  // Uses the correct allele frequency
  CopyAlleleFreq(markernum);

  // Fixing the parents
  if (ind->father==ind2->father)
  {
    commonpar = ind->father;
    halfp1    = ind->mother;
    halfp2    = ind2->mother;
  }
  else
  {
    commonpar = ind->mother;
    halfp1    = ind->father;
    halfp2    = ind2->father;
  }

  commkr = markernumber(commonpar, markernum);
  h1mkr  = markernumber(halfp1, markernum);
  h2mkr  = markernumber(halfp2, markernum);

  comps = 0;
  h1ps  = 0;
  h2ps  = 0;

  // Creating the parents possible genotypes

  if (istypedformarker(commonpar,markernum))
  {
    pheno = cmalloc (sizeof (phenoset));
    memset (pheno,0,sizeof (phenoset));
    pheno->allele1 = min(commkr->allele1, commkr->allele2);
    pheno->allele2 = max(commkr->allele1, commkr->allele2);

    addlist(&comps,pheno);
  }
  else
  {
    for (i =1; i<=numberofalleles; i++)
    {
      for (j =i; j<=numberofalleles; j++)
      {
        pheno = cmalloc (sizeof (phenoset));
        memset (pheno,0,sizeof (phenoset));
        pheno->allele1 = i;
        pheno->allele2 = j;
        addlist(&comps,pheno);
      }
    }
  }
  if (istypedformarker(halfp1,markernum))
  {
    pheno = cmalloc (sizeof (phenoset));
    memset (pheno,0,sizeof (phenoset));
    pheno->allele1 = min(h1mkr->allele1, h1mkr->allele2);
    pheno->allele2 = max(h1mkr->allele1, h1mkr->allele2);

    addlist(&h1ps,pheno);
  }
  else
  {
    for (i =1; i<=numberofalleles; i++)
    {
      for (j =i; j<=numberofalleles; j++)
      {
        pheno = cmalloc (sizeof (phenoset));
        memset (pheno,0,sizeof (phenoset));
        pheno->allele1 = i;
        pheno->allele2 = j;
        addlist(&h1ps,pheno);
      }
    }
  }
  if (istypedformarker(halfp2,markernum))
  {
    pheno = cmalloc (sizeof (phenoset));
    memset (pheno,0,sizeof (phenoset));
    pheno->allele1 = min(h2mkr->allele1, h2mkr->allele2);
    pheno->allele2 = max(h2mkr->allele1, h2mkr->allele2);

    addlist(&h2ps,pheno);
  }
  else
  {
    for (i =1; i<=numberofalleles; i++)
    {
      for (j =i; j<=numberofalleles; j++)
      {
        pheno = cmalloc (sizeof (phenoset));
        memset (pheno,0,sizeof (phenoset));
        pheno->allele1 = i;
        pheno->allele2 = j;
        addlist(&h2ps,pheno);
      }
    }
  }

  // The pheno-sets should no be in order
  nom = 0.0;
  denom = 0.0;


  for (pheno = h1ps; pheno; pheno = pheno->next)
  {
    for (pheno2 = comps; pheno2; pheno2 = pheno2->next)
    {
      for (pheno3 = h2ps; pheno3; pheno3 = pheno3->next)
      {
        tempres = allelefreq[pheno->allele1]*allelefreq[pheno->allele2]*
                  allelefreq[pheno2->allele1]*allelefreq[pheno2->allele2]*
                  allelefreq[pheno3->allele1]*allelefreq[pheno3->allele2];

        // If phenotype ab then multiply by two (ab and ba are equi. prop under HW)
        if (pheno->allele1 != pheno->allele2)
          tempres *= 2;
        if (pheno2->allele1 != pheno2->allele2)
          tempres *= 2;
        if (pheno3->allele1 != pheno3->allele2)
          tempres *= 2;

        for (sibs = ind->sib; sibs; sibs = sibs->next)
        {
          // Checking that they are full sibs and that the sibling is typed
          if (fullsibs(ind,sibs->ind) && istypedformarker(sibs->ind,markernum))
          {
            marker = markernumber(sibs->ind,markernum);
            tempres = tempres * 0.25*(
            HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele1) +
            HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele2) +
            HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele1) +
            HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele2));
	  }
        }
        for (sibs = ind2->sib; sibs; sibs = sibs->next)
        {
          // Checking that they are full sibs and that the sibling is typed
          if (fullsibs(ind2,sibs->ind) && istypedformarker(sibs->ind,markernum))
          {
            marker = markernumber(sibs->ind,markernum);
            tempres = tempres * 0.25*(
            HEpenetrance(marker->allele1, marker->allele2, pheno2->allele1, pheno3->allele1) +
            HEpenetrance(marker->allele1, marker->allele2, pheno2->allele1, pheno3->allele2) +
            HEpenetrance(marker->allele1, marker->allele2, pheno2->allele2, pheno3->allele1) +
            HEpenetrance(marker->allele1, marker->allele2, pheno2->allele2, pheno3->allele2));
	  }
        }      
        marker = markernumber(ind,markernum);
        marker2 = markernumber(ind2,markernum);

        denom += tempres*0.25*(      
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele1) +
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele2) +
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele1) +
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele2))*
               0.25*(      
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele1, pheno3->allele1) +
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele1, pheno3->allele2) +
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele2, pheno3->allele1) +
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele2, pheno3->allele2));
        nom += tempres*0.0625*(
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele1)*
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele1, pheno3->allele1) +
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele1)*
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele1, pheno3->allele2) +
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele2)*
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele2, pheno3->allele1) +
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele1, pheno2->allele2)*
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele2, pheno3->allele2) +
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele1)*
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele1, pheno3->allele1) +
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele1)*
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele1, pheno3->allele2) +
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele2)*
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele2, pheno3->allele1) +
             HEpenetrance(marker->allele1, marker->allele2, pheno->allele2, pheno2->allele2)*
             HEpenetrance(marker2->allele1, marker2->allele2, pheno2->allele2, pheno3->allele2));
      }
    }
  }
  freelist(h2ps);
  freelist(h1ps);
  freelist(comps);

  return(0.5*nom/denom);
}




/*
  Input:

    datatype : 1 = Full sibs
               2 = Half sibs
    markernum: Number of marker used for analysis
    traitnum : Number of trait used as response for analysis
    inclcov  : Vector representing the covariates used in regression
    numcov   : Number of  covariates used

*/

// datatype 1: Full sibs, 2: Half sibs and 3: Both. Only 1 works
float HasemanElston (int datatype, int markernum, int traitnum, int *inclcov, int numcov)
{
  HasemanElstonData *dataset, *obs;
  individual *ind, *ind2;
  quanttrait *qt, *qt2, *qt3;
  OLDMATRIX *x,*transx, *xTx, *xy, *y, *modmean;
  float b, pimean, stdpi;
  float multresvar;
  double se, ttest, testssh;
  int df, sibsize, addsize;
  int i, j, code, datasetsize, effsamplesize, newperson, hasmissing;

  dataset = 0;
  pimean = 0;
  effsamplesize = 0;

  // Checking that data exists and not half/both
  if (!individuals)
  {
    printf("WARNING: No data in memory for Haseman-Elston\n");
    return -1;
  } 
  if (datatype != 1)
  {
    printf("ERROR: Pedipet can only do Haseman-Elston on full sibs\n");
    return -1;
  }

  // Checks that response is ok and that marker is ok
  if (markernum<1 || markernum>numberofmarkers())
  {
    printf("Marker (marker no.: %d) used for Haseman Elston not found\n", markernum);
    return -1;
  }

  // Checks that response is ok and that marker is ok
  if (traitnum<1 || traitnum>numberoftraits())
  {
    printf("Response trait used for Haseman Elston not found\n");
    return -1;
  }


  // Uses the correct allele frequency
  CopyAlleleFreq(markernum);

  forind
  {
    strcpy(ind->tmpstr1, "1");
  }

  df = 0;

  forind
  {
    // newperson is used for calculating the df
    newperson = 0;

    sibsize = 0;

    if (strchr(ind->tmpstr1,'1'))
    {
      strcpy(ind->tmpstr1, "0");
      sibsize = 1;
    }

    for (ind2 = ind->next; ind2; ind2 = ind2->next)
    {
      if (fullsibs(ind, ind2))
      {
        code = 0;
        addsize = 0;

        if (!istypedformarker(ind,markernum) || !istypedformarker(ind2,markernum))
          continue; 
  
	// New observation
        obs = cmalloc (sizeof (HasemanElstonData));
        memset (obs, 0, sizeof (HasemanElstonData));

        obs->ibd   = SingleMarkerIBD(ind, ind2, markernum);

        // Calculating the trait difference
        qt = traitvalue(ind, traitnum);
        qt2 = traitvalue(ind2, traitnum);

        // If either value is missing
        if (qt->sysmiss || qt2->sysmiss)
          code = 1;
        else
          obs->traitdiff = SQR(qt->value-qt2->value);

        // Should maybe check that the trait isn't missing
	// For now: Let the statistical program do that 

        if (!code) // No missing trait difference so add observation to dataset
	{
          #ifdef debug
            printf("Adding pair: %s and %s: (IBD: %1.4f)\n",ind->id, ind2->id, obs->ibd);
          #endif

          if (strchr(ind2->tmpstr1,'1'))
          {  
            sibsize++;
            addsize = 1;
            strcpy(ind2->tmpstr1, "0");          
          }

          addlist(&dataset, obs);
          newperson = 1;

          hasmissing = 0;

      	  // Should now add covariates ()
	  // For now covariates are all quant traits except the one used as response
          for (i = 1 ; i <= numcov; i++)
          {
            // If it is a trait
	    if (inclcov[i] != traitnum && inclcov[i]>0)
	    {
              qt  = traitvalue(ind, inclcov[i]);
              qt2 = traitvalue(ind2, inclcov[i]);

              qt3 = cmalloc (sizeof (quanttrait));
              memset (qt3, 0, sizeof (quanttrait));

              // Add value if none are missing, otherwise sysmiss
              if (qt->sysmiss || qt2->sysmiss)
              {
                qt3->sysmiss = 1;
                hasmissing = 1;
              }
              else
                qt3->value = fabs(qt->value-qt2->value);

              addlist(&obs->covariates, qt3);  
            }
            else if (inclcov[i] == 0) // Put the sex difference in where the response was
	    {
   	      // Now adding sex from the pedigree file info
	      qt3 = cmalloc (sizeof (quanttrait));
	      memset (qt3, 0, sizeof (quanttrait));

      	      qt3->value = fabs(ind->sex-ind2->sex);
	      addlist(&obs->covariates, qt3);
	    }
          }

          if (hasmissing)
          {
            freelist(obs->covariates);
            removelist(&dataset,obs);
            if (addsize)
              sibsize--; 
            newperson = 0;
          }
          else
	  {
            pimean += obs->ibd;
	  }
        }
        else
        {
          free(obs);
	}
      }
    }
    if (sibsize)
    {
      df = df + sibsize - 1;
    }
    
    if (newperson==1)
      effsamplesize++;
  }

  // Should now perform a multiple regression using the relevant covariates
  // Remember to just look at complete cases
  // Assumes that all covariates are relevant by now

  // Creating matrices
  x      = MtxNew(listlen(dataset), numcov+2);
  transx = MtxNew(numcov + 2, listlen(dataset));
  y      = MtxNew(listlen(dataset), 1);
  xTx    = MtxNew(numcov+2, numcov+2);
  xy     = MtxNew(numcov+2, 1);
  modmean= MtxNew(listlen(dataset), 1);

  // Filling out matrices
  i = 1;
  for (obs = dataset; obs; obs = obs ->next)
  {
    y->element[i][1] = obs->traitdiff;        // Squared difference
    x->element[i][numcov+1] = 1;              // Intercept
    x->element[i][numcov+2] = obs->ibd;       // IBD effect
    
    j = 1; // How many markers ?
    for (qt = obs->covariates; qt; qt = qt->next)
    {
      x->element[i][j] = (float)qt->value;    // Abs(difference)
      j++;
    }
    i++;
  }

  MtxTrans(x,transx);
  // Calculating xTy and xTx
  MtxMulti(transx,y,xy);
  MtxMulti(transx, x, xTx);

  //  MtxPrint(x);

  // Should now calculate INV(xTx)
  MtxInver(xTx);
  // Calculating (xTx)^{-1} xTy
  MtxMulti(xTx,xy,xy);  
  // xy now holds the parameter estimates

  #ifdef debug
  MtxPrint(xy);
  #endif

  MtxMulti(x, xy, modmean);

  multresvar = 0.0;
  for (i=1; i<=modmean->rows; i++)
  {
    multresvar += SQR(y->element[i][1]-modmean->element[i][1]);
  }
  multresvar /= (listlen(dataset)-1);  

  // The effective df is number of full sibs minus number of pedigrees
  //  effectivesamplesize -= numberofpedigrees();
  effsamplesize = df - (numcov + 2);


  // Check that there is actually data in the dataset
  if (listlen(dataset) == 0)
  {
    printf("WARNING: No sibpairs to do Haseman Elston analysis on\n");
    return 0;
  }
  else
  {
    printf("Haseman-Elston regression analysis of marker %d\n", markernum);
    printf("%d sibpairs included in the dataset for Haseman-Elston regression\n", listlen(dataset));
  }

  // Takes the mean of all the \hat(\pi) values in the dataset
  pimean = pimean/listlen(dataset);

  // Calculating SD of estimated pi's
  stdpi  = 0;
  for (obs = dataset; obs; obs = obs ->next)
  {
    stdpi += SQR(obs->ibd-pimean);
  }
  stdpi = sqrt(stdpi / (listlen(dataset)-1));

  // Testing the overall mean proportion of alleles ibd
  // To compare with a t-distribution
  if (datatype == 1)       // Full sibs
    ttest = sqrt(listlen(dataset))*((pimean -.5)/stdpi);
  else if (datatype == 2)  // Half sibs
    ttest = sqrt(listlen(dataset))*((pimean -.25)/stdpi);
  else {
    // Should never happen
    ttest = 0;
    exit(255);
  }

  printf("Test for overall mean proportion of alleles ibd: %f (p = %f) \n", ttest, 2*pTdistr(listlen(dataset),fabs(ttest)));

  se = sqrt(xTx->element[numcov+2][numcov+2]*multresvar);
  b  = xy->element[numcov+2][1];

  datasetsize=listlen(dataset);

  // Checks if negative squared variance effect can occur, ie. if b+a < 0
  // If that is true, reestimate the lot without intercept and with 1-ibd as effect
  // XXX Should be based on the residuals !!!  
  /*
  if (xy->element[numcov+2][1]+xy->element[numcov+1][1] < 0)
  {
    printf("Reestimating to have non-negative ibd");
    nyx      = MtxNew(listlen(dataset), numcov+1);
    nytransx = MtxNew(numcov + 1, listlen(dataset));
    nyy      = MtxNew(listlen(dataset), 1);
    nyxTx    = MtxNew(numcov+1, numcov+1);
    nyxy     = MtxNew(numcov+1, 1);

    for (i=1; i<=listlen(dataset); i++)
    {
      for (j=1; j<=numcov+2; j++)
      {
        if (j = numcov+2)
          nyx->element[i][j-1] = 1 - x->element[i][j];
        else
          nyx->element[i][j] = x->element[i][j];
      }
    }

    MtxTrans(nyx,nytransx);
    // Calculating xTy and xTx
    MtxMulti(nytransx,nyy,nyxy);
    MtxMulti(nytransx, nyx, nyxTx);

    // Should now calculate INV(xTx)
    MtxInver(nyxTx);
    // Calculating (xTx)^{-1} xTy
    MtxMulti(nyxTx,nyxy,nyxy);  
    // xy now holds the parameter estimates

    MtxMulti(x, xy, modmean);
    multresvar = 0.0;
    for (i=1; i<=modmean->rows; i++)
    {
      multresvar += SQR(y->element[i][1]-modmean->element[i][1]);
    }
    multresvar /= (listlen(dataset)-1);  

    // Changes sign on the estimated slope
    xy->element[numcov+1][1] = -xy->element[numcov+1][1];

    se = sqrt(xTx->element[numcov+1][numcov+1]*multresvar);
    b  = xy->element[numcov+1][1];

    MtxDel(nyx);
    MtxDel(nyxTx);  
    MtxDel(nytransx);
    MtxDel(nyxy);
    MtxDel(nyy);  
  }
  */
  printf("Estimated slope  : %f (%f)\n",b, se);

  // If the slope is negative
  if (b<0)
    testssh = 1 - pTdistr(effsamplesize,b/se);
  else
    testssh = 1.0;
  printf("T-test statistic : %f (p = %f) [one-sided on %d df]\n", b/se, testssh, effsamplesize);

  // Calculates the sample size correctly now
  //  printf("WARNING: df not correctly modified\n");

  // Cleans up after use
  // Removes the data from memory and NULLs the pointer
  obs = dataset;

  for (obs = dataset; obs; obs = obs ->next)
    freelist(obs->covariates);
  freelist(dataset);
  dataset = 0;

  MtxDel(x);
  MtxDel(xTx);  
  MtxDel(transx);
  MtxDel(xy);
  MtxDel(y);  
  MtxDel(modmean);

  fflush(F);

  return testssh;
}


float MultiPointHasemanElston(int datatype, int traitnum, int *inclcov, int numcov)
{
  HasemanElstonData *dataset, *obs;
  individual *ind, *ind2;
  quanttrait *qt, *qt2, *qt3;
  OLDMATRIX *x,*transx, *xTx, *xy, *y, *modmean;
  OLDMATRIX *pimatrix, *V, *beta, *C, *scaledC, *avpis;
  float b, pimean, stdpi;
  float multresvar, dist, dist2, position;
  double se;
  int df, sibsize, addsize, code2;
  int count, i, j, k, m, code, effsamplesize, newperson, hasmissing;

  pimatrix = NULL;
  effsamplesize = 0;

  dataset = 0;

  avpis = MtxNew(numberofmarkers(), 1);

  // Checking that data exists and not half/both
  if (!individuals)
  {
    printf("WARNING: No data in memory for Haseman-Elston\n");
    return -1;
  } 
  if (datatype != 1)
  {
    printf("ERROR: Pedipet can only do Haseman-Elston on full sibs\n");
    return -1;
  }

  // Checks that response is ok and that marker is ok
  if (traitnum<1 || traitnum>numberoftraits())
  {
    printf("Response trait used for Haseman Elston not found\n");
    return -1;
  }

  // Estimating IBD for all pairs for all markers now
  for (k=1; k<=numberofmarkers(); k++) {
#ifdef debug
    printf("Doing marker %d\n",k);
#endif
    
  pimean = 0;
  effsamplesize = 0;
  code2 = 0;

  // Uses the correct allele frequency
  CopyAlleleFreq(k);

  forind
  {
    strcpy(ind->tmpstr1, "1");
  }

  df = 0;

  forind
  {
    // newperson is used for calculating the df
    newperson = 0;

    sibsize = 0;

    if (strchr(ind->tmpstr1,'1'))
    {
      strcpy(ind->tmpstr1, "0");
      sibsize = 1;
    }

    for (ind2 = ind->next; ind2; ind2 = ind2->next)
    {
      if (fullsibs(ind, ind2))
      {
        code = 0;
        addsize = 0;

        j = 0;
        m = 0;
	// Both persons should be typed for at least two markers
        for (i=1; i<=numberofmarkers(); i++)
        {
          if (istypedformarker(ind,i))
            j++;
          if (istypedformarker(ind2,i))
            m++;
	}
        if (j<2 || m<2)
          continue; 

  
	// New observation
        obs = cmalloc (sizeof (HasemanElstonData));
        memset (obs, 0, sizeof (HasemanElstonData));

        if (istypedformarker(ind,k) && istypedformarker(ind2,k))
          obs->ibd   = SingleMarkerIBD(ind, ind2, k);
        else
          obs->ibd   = -1 ; // The "not avaliable" code

        // Calculating the trait difference
        qt = traitvalue(ind, traitnum);
        qt2 = traitvalue(ind2, traitnum);

        // If either value is missing
        if (qt->sysmiss || qt2->sysmiss)
          code = 1;
        else
          obs->traitdiff = SQR(qt->value-qt2->value);

        // Should maybe check that the trait isn't missing
	// For now: Let the statistical program do that 

        if (!code) // No missing trait difference so add observation to dataset
	{
          #ifdef debug
            printf("Adding pair: %s and %s: (IBD: %1.4f)\n",ind->id, ind2->id, obs->ibd);
          #endif

          if (strchr(ind2->tmpstr1,'1'))
          {  
            sibsize++;
            addsize = 1;
            strcpy(ind2->tmpstr1, "0");          
          }

          addlist(&dataset, obs);
          newperson = 1;

          hasmissing = 0;
      	  // Should now add covariates ()
	  // For now covariates are all quant traits except the one used as response
          for (i = 1 ; i <= numcov; i++)
          {
            // If it is a trait
	    if (inclcov[i] != traitnum && inclcov[i]>0)
	    {
              qt  = traitvalue(ind, inclcov[i]);
              qt2 = traitvalue(ind2, inclcov[i]);

              qt3 = cmalloc (sizeof (quanttrait));
              memset (qt3, 0, sizeof (quanttrait));

              // Add value if none are missing, otherwise sysmiss
              if (qt->sysmiss || qt2->sysmiss)
              {
                qt3->sysmiss = 1;
                hasmissing = 1;
              }
              else
                qt3->value = fabs(qt->value-qt2->value);

              addlist(&obs->covariates, qt3);  
            }
            else if (inclcov[i] == 0) // Put the sex difference in where the response was
	    {
   	      // Now adding sex from the pedigree file info
	      qt3 = cmalloc (sizeof (quanttrait));
	      memset (qt3, 0, sizeof (quanttrait));

      	      qt3->value = fabs(ind->sex-ind2->sex);
	      addlist(&obs->covariates, qt3);
	    }
          }
          if (hasmissing)
          {
            freelist(obs->covariates);
            removelist(&dataset,obs);
            if (addsize)
              sibsize--; 
            newperson = 0;
          }
          else
	  {
            if (obs->ibd>=0)
            {
              pimean += obs->ibd;
              code2++;
            }
	  }
        }
        else
        {
          free(obs);
	}
      }
    }
    if (sibsize)
    {
      df = df + sibsize - 1;
    }
    
    if (newperson==1)
      effsamplesize++;
  }

  // If the first marker then
  if (k==1) {
    pimatrix = MtxNew(listlen(dataset),numberofmarkers());
  }
  i = 1;
  for (obs = dataset; obs; obs = obs->next)
  {
    pimatrix->element[i][invorder[k]] = obs->ibd;
    i++;
  }
  avpis->element[invorder[k]][1] = pimean / code2;

  // Delete dataset before last marker
  if (k<numberofmarkers()) {
    for (obs = dataset; obs; obs = obs ->next)
      freelist(obs->covariates);
    freelist(dataset);
    dataset = 0;
  }
}
 
  x      = MtxNew(listlen(dataset), numcov+2);
  transx = MtxNew(numcov + 2, listlen(dataset));
  y      = MtxNew(listlen(dataset), 1);
  xTx    = MtxNew(numcov+2, numcov+2);
  xy     = MtxNew(numcov+2, 1);
  modmean= MtxNew(listlen(dataset), 1);


  // Filling out the y and x matrices
  i = 1;
  for (obs = dataset; obs; obs = obs ->next)
  {
    y->element[i][1] = obs->traitdiff;        // Squared difference
    x->element[i][numcov+1] = 1;              // Intercept
    x->element[i][numcov+2] = obs->ibd;       // IBD effect
    
    j = 1; // How many markers ?
    for (qt = obs->covariates; qt; qt = qt->next)
    {
      x->element[i][j] = (float)qt->value;    // Abs(difference)
      j++;
    }
    i++;
  }


  // Preparing V, C, beta and pimatrix

  V = MtxNew(numberofmarkers(),numberofmarkers());
  C = MtxNew(numberofmarkers(), 1);
  scaledC = MtxNew(numberofmarkers(), 1);
  beta = MtxNew(numberofmarkers(), 1);

  // Calculating the observed marker variances
  // Using beta to keep for now
  for (i=1; i<= numberofmarkers(); i++)
  {
    pimean = 0;
    count = 0;
    for (j=1; j<=pimatrix->rows; j++)
    {    
      if (pimatrix->element[j][i]>-1)
      {
        pimean += pimatrix->element[j][i];
        count ++;
      }
    }
    pimean /= count;
    stdpi = 0;
    count = 0;
    for (j=1; j<=pimatrix->rows; j++)
    {
      if (pimatrix->element[j][i]>-1)
      {
        stdpi += SQR(pimatrix->element[j][i]-pimean);
        count ++;
      }
    }
    stdpi /= (count-1);
    beta->element[i][1] = stdpi;
  }

  for (i=1; i<= numberofmarkers(); i++)
  {
    for (j=i; j<= numberofmarkers(); j++)
    {
      if (i==j)
      {
        V->element[i][j] = beta->element[i][1];
      }
      else
      {
        dist = 0.0;
        for (k=i; k<j; k++)
        {
          dist += MapFunction(distance[k]);
        }
        dist = InverseMapFunction(dist);
        V->element[i][j] = 8*SQR(1-2*dist)*beta->element[j][1]*beta->element[i][1];
        V->element[j][i] = V->element[i][j];
      }
    }
  }
  // V and beta are set up in correct order
  MtxCopy(beta, C);
//  MtxPrint(V);
  MtxInver(V);

  // Have now calculated the inverse of V

  // Calculates the total length of the chromosome
  dist2 = 0.0;
  for (i=1; i< numberofmarkers(); i++)
    dist2 += 100*MapFunction(distance[i]);

  #ifdef EXPORT
  cfopen("he-mp.r", "w");
  fprintf(F,"hedata <- matrix(c(");
  #endif

  printf("Pos.(cM)  Neg.t-test  p-value   LOD\n------------------------------------\n");

  for (m = -6; m<= dist2+6; m = m + 2.0) {

  // Position is in cM
  position = m;

  // Estimates betas
  MtxCopy (C, scaledC);

  for (i=1; i<= numberofmarkers(); i++)
  {
    dist = InverseMapFunction(fabs(MarkerDistance(position,i))/100);
    scaledC->element[i][1] *= SQR(1-2*dist);
  }
  MtxMulti(V, scaledC, beta);

  for (i = 1; i<=pimatrix->rows; i++)
  {
    pimean = 0.0;
    stdpi = 0.0;

    for (j=1; j<=pimatrix->cols ; j++)
      if (pimatrix->element[i][j]>-1)
        stdpi  += beta->element[j][1]; 
    // stdpi now holds the total weight used

    count = 0;
    for (j=1; j<=pimatrix->cols ; j++)
    {     
      if (pimatrix->element[i][j]>-1)
      {
        pimean += beta->element[j][1]*(pimatrix->element[i][j]-avpis->element[j][1]);
        count++;
      }
    }

    // Constraining to [0,1]
    switch(datatype)
    {
      case 1: x->element[i][numcov+2] = 0.5 + pimean;
              x->element[i][numcov+2] = min(x->element[i][numcov+2],1);
              x->element[i][numcov+2] = max(x->element[i][numcov+2],0);
              break;
      case 2: x->element[i][numcov+2] = 0.25 + pimean;
              x->element[i][numcov+2] = min(x->element[i][numcov+2],.5);
              x->element[i][numcov+2] = max(x->element[i][numcov+2],0);
              break;
    }

  }
  //  MtxPrint(x);
//  MtxPrint(beta);

  MtxTrans(x,transx);
  // Calculating xTy and xTx
  MtxMulti(transx,y,xy);
  MtxMulti(transx, x, xTx);

  // Should now calculate INV(xTx)
  MtxInver(xTx);

  // Calculating (xTx)^{-1} xTy
  MtxMulti(xTx,xy,xy);  

  // xy now holds the parameter estimates
  MtxMulti(x, xy, modmean);

  multresvar = 0.0;
  for (i=1; i<=modmean->rows; i++)
  {
    multresvar += SQR(y->element[i][1]-modmean->element[i][1]);
  }
  multresvar /= (listlen(dataset)-1);  

  se = sqrt(xTx->element[numcov+2][numcov+2]*multresvar);
  b  = xy->element[numcov+2][1];

//  printf("%5.1f  %5.3f\n",position,b/se);

  printf("%4d     %8.3f     %6.4f   %5.2f\n", m, -b/se, b<0 ? 1 - pTdistr(effsamplesize,b/se) : 1.0, SQR(max(-b/se,0))/4.61);

  #ifdef EXPORT
  if (m > -6)
    fprintf(F,",\n");
  #ifdef chi2
  fprintf(F,"%f, %f", position, SQR(max(-b/se,0)));
  #else
  fprintf(F,"%f, %f", position, b/se);
  #endif
  #endif

}  

  MtxDel(V);
  MtxDel(C);
  MtxDel(scaledC);
  MtxDel(beta);
  MtxDel(pimatrix);

  MtxDel(x);
  MtxDel(xTx);  
  MtxDel(transx);
  MtxDel(xy);
  MtxDel(y);  
  MtxDel(modmean);
  MtxDel(avpis);

  for (obs = dataset; obs; obs = obs ->next)
    freelist(obs->covariates);
  freelist(dataset);
  dataset = 0;

  fflush(F);

  #ifdef EXPORT
  fprintf(F, "),byrow=T, ncol=2)\n\n");
  fprintf(F,"plot(hedata[,1], -hedata[,2], xlim=c(min(hedata[,1]),max(hedata[,1])),\n");
  fprintf(F,"ylim=c(-1,max(-hedata[,2],4)),xlab=\"Position (cM)\", ylab=\"-t\", type=\"l\")\n");
  fprintf(F,"title(main=\"Multipoint Haseman-Elston\")\n");
  fprintf(F,"mtext(\"Pedipet v%s\", side=1, line=3, cex=.3, adj=1)\n", VERSION);

// If p-lines
  fprintf(F, "lines(c(min(hedata[,1]),max(hedata[,1])), rep(qt(0.95,%d),2), lty=3)\n", effsamplesize);
  fprintf(F, "lines(c(min(hedata[,1]),max(hedata[,1])), rep(qt(0.99,%d),2), lty=3)\n", effsamplesize);
  fprintf(F, "lines(c(min(hedata[,1]),max(hedata[,1])), rep(qt(0.999,%d),2), lty=3)\n", effsamplesize);


// Plot marker names on plot
  dist = 0.0;
  for (i=1; i<=numberofmarkers(); i++)
  {
    fprintf(F, "text(%6.2f, -1.1, \"%s\", srt=90, crt=90, cex=.5, adj=0)\n",dist,GetName(markernames,order[i]));
    if (i<numberofmarkers())
      dist += 100*MapFunction(distance[i]);
  }

  fclose(F);
#endif
  return 0;
}

