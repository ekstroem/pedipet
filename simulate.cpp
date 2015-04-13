/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * This file is part of the statistical program Pedipet.  
 *
 * Functions to simulate genetic marker data and families for QTL
 * linkage analysis
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



/*



   Claus Ekstrøm 1997-2003

 */

#include <cstdlib>
#include <cstring>
#include <ctime>

#include <sys/types.h>
#include <unistd.h>


#include "simulate.h"
#include "estimate.h"
#include "matrix.h"

extern "C" {
#include "linpack.c"
#include "com.c"
}


//#define DOS  // If compiled under DOS
//#define fullyinformative

// Change seed during calculations
// #define seed         

// Allow for X-linked datasets
#define xlinked

// Generates 4-offspring nuclear families with exactly two males and two females
#define twoofeach

// Generates 4-offspring nuclear families with exactly two males and two females
// #define foundershomoz


// Corrections for Kosambi definition
#define KOSAMBI

#define DOSAGECOMPENSATION

#ifndef DOSAGECOMPENSATION
#define DOMINANCE
#endif

 // Simulates a dataset of 'families' nuclear families each having 'nooff'
 // offspring (full sibs). It uses the global allele frequencies
 // and creates a single QTL dist (rec.frac) away from the marker, 
 // with an additive effect of effect SD. 
 // The QTL is found in prev % of the population
 // Mean of unaffected individuals is 0, SD is 1

individual *simulateQTL(int families, int nooff, double dist, double addeffect, double domeffect, double prev)
{
  int i, j, nogenes, count, newmarkers, persons;
  individual *ind, *ind2, *father, *mother, *finaldata;
  markerlist *marker;
  quanttrait *qt;
  namelist *nl;
  double randomnumber, res;
  float tmp, heritab;

  finaldata = 0;

  tmp = addeffect + domeffect*(1-2*prev);
  heritab = 2*prev*(1-prev)*SQR(tmp);
  tmp = 2*prev*(1-prev)*domeffect;
  heritab += SQR(tmp);

  persons = 0;

  for (i=1; i<=families; i++)
    {
      // Generates the founders

      // First the father
      father = (individual *) cmalloc (sizeof (individual));
      memset (father,0,sizeof (individual));

      sprintf(father->pedigree, "F%d",i);
      sprintf(father->id, "F%d-1",i);
      father->sex = S_MALE;

      // Gives the person one copy of the QTL gene with probability prev
      // and two copies with probability prev*prev
      if (genunf(0,1)<prev)
	strcpy(father->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(father->tmpstr2,"1");

      marker = (markerlist *) cmalloc (sizeof (markerlist));
      memset (marker,0,sizeof (markerlist));

      count = 1;
      res = allelefreq[count];

     // Gives the father the
      randomnumber = genunf(0,1);
      while(randomnumber>res)
	{
	  count++;
	  res += allelefreq[count];
	}
      marker->allele1 = count;

      count = 1;
      res = allelefreq[count];
      randomnumber = genunf(0,1);
      while(randomnumber>res)
	{
	  count++;
	  res += allelefreq[count];
	}
      marker->allele2 = count;

      addlist(&father->marker, marker);


      // This qt represents the observed trait
      qt = (quanttrait *) cmalloc (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));

      nogenes = strlen(father->tmpstr1)+strlen(father->tmpstr2);

      qt->sysmiss = 0;
      qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
      if (nogenes==1)
	qt->value += domeffect;


      addlist(&father->qttrait,qt);

      // Makes an additional irrelevant covariate
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 2);
      addlist(&father->qttrait,qt);

      addlist (&finaldata,father);
      father->globalid = ++persons;
      father->localid = persons;


     // Then the mother
      mother = (individual *) cmalloc (sizeof (individual));
      memset (mother,0,sizeof (individual));

      sprintf(mother->pedigree, "F%d",i);
      sprintf(mother->id, "F%d-2",i);
      mother->sex = S_FEMALE;

      // Gives the person one copy of the QTL gene with probability prev
      // and two copies with probability prev*prev
      if (genunf(0,1)<prev)
	strcpy(mother->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(mother->tmpstr2,"1");


      marker = (markerlist *) cmalloc (sizeof (markerlist));
      memset (marker,0,sizeof (markerlist));

      count = 1;
      res = allelefreq[count];
      randomnumber = genunf(0,1);
      while(randomnumber>res)
	{
	  count++;
	  res += allelefreq[count];
	}
      marker->allele1 = count;

      count = 1;
      res = allelefreq[count];
      randomnumber = genunf(0,1);
      while(randomnumber>res)
	{
	  count++;
	  res += allelefreq[count];
	}
      marker->allele2 = count;

      addlist(&mother->marker, marker);

      // This qt represents the observed trait
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));

      nogenes = strlen(mother->tmpstr1)+strlen(mother->tmpstr2);

      qt->sysmiss = 0;
      qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
      if (nogenes==1)
	qt->value += domeffect;


      addlist(&mother->qttrait,qt);

      // Makes an additional irrelevant covariate
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 2);

      addlist(&mother->qttrait,qt);

      addlist (&finaldata,mother);
      mother->globalid = ++persons;
      mother->localid = persons;

     // Generates the offspring
      for (j=0; j<nooff; j++)
	{
	  ind = (individual *) cmalloc (sizeof (individual));
	  memset (ind,0,sizeof (individual));

	  sprintf(ind->pedigree, "F%d",i);
	  sprintf(ind->id, "F%d-%d",i,j+3);

	  if (genunf(0,1)<.5)
	    ind->sex = S_MALE;
	  else
	    ind->sex = S_FEMALE;

	  ind->father = father;
	  ind->mother = mother;

	  adduniqueidlist (&ind->father->offspring, ind);
	  adduniqueidlist (&ind->father->mate, ind->mother);
	  adduniqueidlist (&ind->mother->offspring, ind);
	  adduniqueidlist (&ind->mother->mate,ind->father); 

	  marker = (markerlist *) cmalloc (sizeof (markerlist));
	  memset (marker,0,sizeof (markerlist));
	  if (ignbin(1,.5))
	    {
	      marker->allele1 = father->marker->allele1;
	      // If recombination has occured
	      if (genunf(0,1)<dist)
		strcpy(ind->tmpstr1,father->tmpstr2);
	      else
		strcpy(ind->tmpstr1,father->tmpstr1);
	    }
	  else
	    {
	      marker->allele1 = father->marker->allele2;      
	      // If recombination has occured
	      if (genunf(0,1)<dist)
		strcpy(ind->tmpstr1,father->tmpstr1);
	      else
		strcpy(ind->tmpstr1,father->tmpstr2);
	    }

	  if (ignbin(1,.5))
	    {
	      marker->allele2 = mother->marker->allele1;
	      // If recombination has occured
	      if (genunf(0,1)<dist)
		strcpy(ind->tmpstr2,mother->tmpstr2);
	      else
		strcpy(ind->tmpstr2,mother->tmpstr1);

	    }
	  else
	    {
	      marker->allele2 = mother->marker->allele2;      
	      // If recombination has occured
	      if (genunf(0,1)<dist)
		strcpy(ind->tmpstr2,mother->tmpstr1);
	      else
		strcpy(ind->tmpstr2,mother->tmpstr2);

	    }
	  addlist(&ind->marker, marker);

	  nogenes = strlen(ind->tmpstr1)+strlen(ind->tmpstr2);
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0 + addeffect*(nogenes-1), 1);
	  if (nogenes==1)
	    qt->value += domeffect;


	  addlist(&ind->qttrait,qt);

	  // Makes an additional irrelevant covariate
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0; 
	  qt->value = gennor(0, 2);
	  addlist(&ind->qttrait,qt);

	  addlist (&finaldata,ind);
	  ind->globalid = ++persons;
	  ind->localid = persons;
	}
    }

  // Fixing sibs
  for (ind = finaldata; ind; ind = ind->next)
    {
      for (ind2 = ind->next; ind2; ind2 = ind2->next)
	{
	  // If same father or mother then add to list of both
	  if (((ind->father == ind2->father) && (ind->father)) || 
	      ((ind->mother == ind2->mother) && (ind->mother)))
	    {
	      adduniqueidlist(&ind->sib, ind2);
	      adduniqueidlist(&ind2->sib, ind);
	    }
	}
    }


  printf("Simulated %d families with %d offspring each.\n", families, nooff);
  // Fixes the order etc.

  InitializeFrequencies(1);

  newmarkers = 1;

  // First set the marker names:
  for (i=1; i<=newmarkers; i++)
    {
      nl = (namelist *) cmalloc (sizeof (namelist));
      memset (nl,0,sizeof (namelist));

      sprintf(nl->name, "Marker %d",i);
      addlist(&markernames, nl);
    }

  // Then add trait names
  for (i=1; i<=2; i++)
    {
      nl = (namelist *) cmalloc (sizeof (namelist));
      memset (nl,0,sizeof (namelist));

      sprintf(nl->name, "Trait %d",i);
      addlist(&traitnames, nl);
    }

  // XXXXX
  // Should fix markernumbers
  // Initialize order and distances
  order = ivector(1,newmarkers);
  for (i=1; i<=1; i++)
    order[i] = i;

  if (newmarkers>1)
    {
      distance = vector(1, newmarkers-1);
      for (i=1; i<newmarkers; i++)
	distance[i] = 0.1;
    }
  else
    distance = 0;

  printf("The heritability of the simulated QTL is: %5.3f (Add. effect: %5.3f, dom. dev. effect: %5.3f)\n", heritab / (1 + heritab), addeffect, domeffect);

  return finaldata; 
}



int roll_1Dx(int n)
{
  static int      old_n = 0;
  static long     max;
  long            rnd;
  /*
   * Go
   */
  if (old_n != n) {
    old_n = n;
    max = n;
    while ((max & 0xFF000000) == 0)
      max <<= 1;
  }
  while (1) {
#ifdef DOS
    rnd = rand();
#else
    rnd = lrand48() >> 6;
#endif
    if (rnd < max)
      break;
  }
  return rnd % n;
}


// Assumes order is ok. Should be fixed
// families each with nooff offspring
// A single QTL a position location with an effect of effect
//
// The QTL must be inside the chromosome or else it won't work
individual *SimulateChromosome(int families, int nooff, double location, double addeffect, double domeffect, double prev, double resadd, double resdom)
{
  int i, j, nogenes, count, mkrnum, NMarkers;
  individual *ind, *ind2, *father, *mother, *finaldata;
  markerlist *marker, *fmkr, *mmkr, *mkr;
  quanttrait *qt;
  namelist *nl;
  double randomnumber, res, recomb;
  markerlist *fixdata[MAXMARKERS+1];
  float dist, mindist, dgenefreq, heritab, tmp;
  int closestmarker, diseasecopies, persons;

  if (!allfreq)
    {
      printf("WARNING: No marker allele frequency information avaliable to simulate data.\n         Remember to load a parameter file first\n");
      return 0;
    }

  persons = 0;
  finaldata = 0;

  tmp = addeffect + domeffect*(1-2*prev);
  heritab = 2*prev*(1-prev)*SQR(tmp);
  tmp = 2*prev*(1-prev)*domeffect;
  heritab += SQR(tmp);

  // finds the marker closest to the left of the position
  mindist = 100000;
  dist = 0.0;
  closestmarker = 1;
  for (i=1; i<numberofmarkers(); i++)  {
    dist += MapFunction(distance[i])*100;
    if (location > dist)
      closestmarker = i+1;
  }
  // dgenefreq holds the recomb fraction between left marker and QTL
  dgenefreq = InverseMapFunction(fabs(MarkerDistance(location, closestmarker))/100);

  for (i=1; i<=families; i++)  {
      // Generates the founders

      // First the father
      father = (individual *) cmalloc (sizeof (individual));
      memset (father,0,sizeof (individual));

      sprintf(father->pedigree, "F%d",i);
      sprintf(father->id, "%d-1",i);
      father->sex = S_MALE;

      // Gives the person one copy of the QTL gene with probability 2*prev(1-prev)
      // and two copies with probability prev*prev
      if (genunf(0,1)<prev)
	strcpy(father->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(father->tmpstr2,"1");

#ifdef foundershomoz
      strcpy(father->tmpstr1,"1");
      strcpy(father->tmpstr2,"1");
#endif      

      NMarkers = listlen(allfreq);
      for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {
	marker = (markerlist *) cmalloc (sizeof (markerlist));
	memset (marker,0,sizeof (markerlist));
	
	// Gets the relevant allele frequency for the marker
	CopyAlleleFreq(mkrnum);
	
	count = 1;
	res = allelefreq[count];
	// Gives the father allele 1
	randomnumber = genunf(0,1);
	while(randomnumber>res)  {
	  count++;
	  res += allelefreq[count];
	}
	marker->allele1 = count;

	  // and allele 2
	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res)
	    {
	      count++;
	      res += allelefreq[count];
	    }
	  marker->allele2 = count;

#ifdef fullyinformative
	  marker->allele1 = 1;
	  marker->allele2 = 2;      
#endif

	  addlist(&father->marker, marker);
	}

      // This qt represents the observed trait
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));

      nogenes = strlen(father->tmpstr1)+strlen(father->tmpstr2);

      qt->sysmiss = 0;
      qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
      if (nogenes==1)
	qt->value += domeffect;

      addlist(&father->qttrait,qt);

      // Makes an additional irrelevant covariate
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 2);
      addlist(&father->qttrait,qt);

      addlist (&finaldata,father);
      father->globalid = ++persons;
      father->localid = persons;


      // Then the mother
      mother = (individual *) cmalloc (sizeof (individual));
      memset (mother,0,sizeof (individual));

      sprintf(mother->pedigree, "F%d",i);
      sprintf(mother->id, "%d-2",i);
      mother->sex = S_FEMALE;

      // Gives the person one copy of the QTL gene with probability prev
      // and two copies with probability prev*prev
      if (genunf(0,1)<prev)
	strcpy(mother->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(mother->tmpstr2,"1");

#ifdef foundershomoz
      strcpy(mother->tmpstr1,"");
      strcpy(mother->tmpstr2,"");
#endif      



      for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++)
	{

	  // Gets the relevant allele frequency for the marker
	  CopyAlleleFreq(mkrnum);

	  marker = (markerlist *) cmalloc (sizeof (markerlist));
	  memset (marker,0,sizeof (markerlist));

	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res)
	    {
	      count++;
	      res += allelefreq[count];
	    }
	  marker->allele1 = count;

	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res)
	    {
	      count++;
	      res += allelefreq[count];
	    }
	  marker->allele2 = count;

#ifdef fullyinformative
	  marker->allele1 = 3;
	  marker->allele2 = 4;      
#endif


	  addlist(&mother->marker, marker);
	}

      // This qt represents the observed trait
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));

      nogenes = strlen(mother->tmpstr1)+strlen(mother->tmpstr2);

      qt->sysmiss = 0;
      qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
      if (nogenes==1)
	qt->value += domeffect;


      addlist(&mother->qttrait,qt);

      // Makes an additional irrelevant covariate
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 2);

      addlist(&mother->qttrait,qt);

      addlist (&finaldata,mother);
      mother->globalid = ++persons;
      mother->localid = persons;


      // Should now generate the children. 

      // Generates the offspring
      for (j=0; j<nooff; j++)
	{
	  ind = (individual *) cmalloc (sizeof (individual));
	  memset (ind,0,sizeof (individual));

	  diseasecopies = 0;

	  sprintf(ind->pedigree, "F%d",i);
	  sprintf(ind->id, "%d-%d",i,j+3);

	  if (genunf(0,1)<.5)
	    ind->sex = S_MALE;
	  else
	    ind->sex = S_FEMALE;

#ifdef twoofeach
	  if (j<2)
	    ind->sex = S_MALE;
	  else
	    ind->sex = S_FEMALE;
#endif

	  ind->father = father;
	  ind->mother = mother;

	  adduniqueidlist (&ind->father->offspring, ind);
	  adduniqueidlist (&ind->father->mate, ind->mother);
	  adduniqueidlist (&ind->mother->offspring, ind);
	  adduniqueidlist (&ind->mother->mate,ind->father); 

	  // tmpstr 1 and 2 hold information about which haplotype
	  // the person received last time. Starts random
	  sprintf(ind->tmpstr1, "%d", (int) ignbin(1,.5));
	  sprintf(ind->tmpstr2, "%d", (int) ignbin(1,.5));

	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++)
	    {

	      marker = (markerlist *) cmalloc (sizeof (markerlist));
	      memset (marker,0,sizeof (markerlist));

	      fmkr = markernumber(father, order[mkrnum]);
	      mmkr = markernumber(mother, order[mkrnum]);

	      if (mkrnum>1)
		{

		  if (mkrnum==closestmarker+1)
		    recomb = InverseMapFunction(fabs(MarkerDistance(location, mkrnum))/100);
		  else
		    recomb = distance[mkrnum-1];

		  // Checks for recombination for father
		  if (genunf(0,1)<recomb)
		    sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
		  // and mother
		  if (genunf(0,1)<recomb)
		    sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));
		}
	      // Give the person the allele
	      if (!strcmpl(ind->tmpstr1,"0"))
		marker->allele1 = fmkr->allele1;
	      else
		marker->allele1 = fmkr->allele2;

	      if (!strcmpl(ind->tmpstr2,"0"))
		marker->allele2 = mmkr->allele1;
	      else
		marker->allele2 = mmkr->allele2;

	      addlist(&ind->marker, marker);


	      // Adds disease genes
	      if (mkrnum == closestmarker)
		{
		  // If no recombination between marker and disease QTL
		  if (genunf(0,1)>=dgenefreq)
		    {
		      if (!strcmpl(ind->tmpstr1,"0") && strlen(ind->father->tmpstr1))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr1,"1") && strlen(ind->father->tmpstr2))
			diseasecopies++;
		    }
		  else
		    {
		      if (!strcmpl(ind->tmpstr1,"0") && strlen(ind->father->tmpstr2))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr1,"1") && strlen(ind->father->tmpstr1))
			diseasecopies++;
		      // Makes a recombination
		      sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
		    }


		  // If no recombination between marker and disease QTL
		  if (genunf(0,1)>=dgenefreq)
		    {
		      if (!strcmpl(ind->tmpstr2,"0") && strlen(ind->mother->tmpstr1))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr2,"1") && strlen(ind->mother->tmpstr2))
			diseasecopies++;
		    }
		  else
		    {
		      if (!strcmpl(ind->tmpstr2,"0") && strlen(ind->mother->tmpstr2))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr2,"1") && strlen(ind->mother->tmpstr1))
			diseasecopies++;
		      // Makes a recombination
		      sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));
		    }

		}

	    }

	  nogenes = diseasecopies;
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0 + addeffect*(nogenes-1), 1);
	  if (nogenes==1)
	    qt->value += domeffect;

	  addlist(&ind->qttrait,qt);

	  // Makes an additional irrelevant covariate
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0, 2);
	  addlist(&ind->qttrait,qt);

	  addlist (&finaldata,ind);
	  ind->globalid = ++persons;
	  ind->localid = persons;

	}

    }


  // Fixing sibs
  for (ind = finaldata; ind; ind = ind->next)
    {
      for (ind2 = ind->next; ind2; ind2 = ind2->next)
	{
	  // If same father or mother then add to list of both
	  if (((ind->father == ind2->father) && (ind->father)) || 
	      ((ind->mother == ind2->mother) && (ind->mother)))
	    {
	      adduniqueidlist(&ind->sib, ind2);
	      adduniqueidlist(&ind2->sib, ind);
	    }
	}
    }

  // Needs to reorder the markers back for non-founders
  // so everything else works fine

  for (ind = finaldata; ind; ind = ind->next)
    {
      if (!founder(ind))
	{
	  i = 2;
	  memset(fixdata, 0, sizeof(fixdata));
	  fixdata[1] = ind->marker;
	  for (mkr = ind->marker; mkr; mkr = mkr->next)
	    {
	      fixdata[i] = mkr->next;
	      i++;
	    }

	  ind->marker = fixdata[invorder[1]];
	  for (i = 1; i<= numberofmarkers(); i++)
	    {
	      mkr = fixdata[i];
	      if (order[i] < numberofmarkers())
		mkr->next = fixdata[invorder[order[i]+1]];
	      else
		mkr->next = 0;
	    }
	}
    }

  // Then add trait names
  // Not necessary to add marker names as they've been read from the
  // parameter file already
  for (i=1; i<=2; i++)
    {
      nl = (namelist *) cmalloc (sizeof (namelist));
      memset (nl,0,sizeof (namelist));

      sprintf(nl->name, "Trait %d",i);
      addlist(&traitnames, nl);
    }


#ifdef twoofeach
  puts("WARNING: Simulating two of each sex !!! \n");
#endif
#ifdef foundershomoz
  puts("WARNING: Simulating doubly homozygous parents !!! \n");
#endif      



#ifdef fullyinformative
  printf("Simulated %d fully informative families with %d offspring each\n", families, nooff);
#else
  printf("Simulated %d families with %d offspring each\n", families, nooff);
#endif
  printf("The heritability of the simulated QTL is: %5.3f (Add. effect: %5.3f, dom. dev. effect: %5.3f) positioned at %5.3f\n", heritab / (1 + heritab), addeffect, domeffect, location);

  // Fixes the order etc.

  return finaldata; 
}












// Assumes order is ok. Should be fixed
// families each with nooff offspring
// A single QTL a position location with an effect of effect
//
// The QTL must be inside the chromosome or else it won't work
individual *SimulateKosambiChromosome(int families, int nooff, double location, double addeffect, double domeffect, double prev, double resadd, double resdom)
{
  int i, j, nogenes, count, mkrnum, NMarkers;
  individual *ind, *ind2, *father, *mother, *finaldata;
  markerlist *marker, *fmkr, *mmkr, *mkr;
  quanttrait *qt;
  namelist *nl;
  double randomnumber, res, gamma11, gamma10;
  // Set the initial distance to the first distance
  double recomb = distance[1];
  markerlist *fixdata[MAXMARKERS+1];
  float dist, mindist, dgenefreq, heritab, tmp, tal;
  int closestmarker, diseasecopies, persons, gamma, mamma;

  if (!allfreq) {
    printf("WARNING: No marker allele frequency information avaliable to simulate data.\n         Remember to load a parameter file first\n");
    return 0;
  }

  printf("KOSAMBI data\n");

  persons = 0;
  finaldata = 0;

  tmp = addeffect + domeffect*(1-2*prev);
  heritab = 2*prev*(1-prev)*SQR(tmp);
  tmp = 2*prev*(1-prev)*domeffect;
  heritab += SQR(tmp);

  // Finds the marker closest to the left of the position
  mindist = 100000;
  dist = 0.0;
  closestmarker = 1;
  for (i=1; i<numberofmarkers(); i++) {
    dist += MapFunction(distance[i])*100;
    if (location > dist)
      closestmarker = i+1;
  }

  // dgenefreq holds the recomb fraction between left marker and QTL
  dgenefreq = InverseMapFunction(fabs(MarkerDistance(location, closestmarker))/100);

  for (i=1; i<=families; i++) {
    // Generates the founders

    // First the father
    father = (individual *) cmalloc (sizeof (individual));
    memset (father,0,sizeof (individual));

    sprintf(father->pedigree, "F%d",i);
    sprintf(father->id, "%d-1",i);
    father->sex = S_MALE;

    // Gives the person one copy of the QTL gene with probability 2*prev(1-prev)
    // and two copies with probability prev*prev
    if (genunf(0,1)<prev)
      strcpy(father->tmpstr1,"1");
    if (genunf(0,1)<prev)
      strcpy(father->tmpstr2,"1");
    
    // Assign markers to the father
    NMarkers = listlen(allfreq);
    for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {
      marker = (markerlist *) cmalloc (sizeof (markerlist));
      memset (marker,0,sizeof (markerlist));
      
	  // Gets the relevant allele frequency for the marker
	  CopyAlleleFreq(mkrnum);

	  count = 1;
	  res = allelefreq[count];
	  // Gives the father allele 1
	  randomnumber = genunf(0,1);
	  while(randomnumber>res)
	    {
	      count++;
	      res += allelefreq[count];
	    }
	  marker->allele1 = count;

	  // and allele 2
	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res)
	    {
	      count++;
	      res += allelefreq[count];
	    }
	  marker->allele2 = count;

#ifdef fullyinformative
	  marker->allele1 = 1;
	  marker->allele2 = 2;      
#endif

	  addlist(&father->marker, marker);
	}

      // This qt represents the observed trait
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));

      nogenes = strlen(father->tmpstr1)+strlen(father->tmpstr2);

      qt->sysmiss = 0;
      qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
      if (nogenes==1)
	qt->value += domeffect;

      addlist(&father->qttrait,qt);

      // Makes an additional irrelevant covariate
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 2);
      addlist(&father->qttrait,qt);

      addlist (&finaldata,father);
      father->globalid = ++persons;
      father->localid = persons;


      // Then the mother
      mother = (individual *) cmalloc (sizeof (individual));
      memset (mother,0,sizeof (individual));

      sprintf(mother->pedigree, "F%d",i);
      sprintf(mother->id, "%d-2",i);
      mother->sex = S_FEMALE;

      // Gives the person one copy of the QTL gene with probability prev
      // and two copies with probability prev*prev
      if (genunf(0,1)<prev)
	strcpy(mother->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(mother->tmpstr2,"1");


      for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++)
	{

	  // Gets the relevant allele frequency for the marker
	  CopyAlleleFreq(mkrnum);

	  marker = (markerlist *) cmalloc (sizeof (markerlist));
	  memset (marker,0,sizeof (markerlist));

	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res)
	    {
	      count++;
	      res += allelefreq[count];
	    }
	  marker->allele1 = count;

	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res)
	    {
	      count++;
	      res += allelefreq[count];
	    }
	  marker->allele2 = count;

#ifdef fullyinformative
	  marker->allele1 = 3;
	  marker->allele2 = 4;      
#endif


	  addlist(&mother->marker, marker);
	}

      // This qt represents the observed trait
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));

      nogenes = strlen(mother->tmpstr1)+strlen(mother->tmpstr2);

      qt->sysmiss = 0;
      qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
      if (nogenes==1)
	qt->value += domeffect;


      addlist(&mother->qttrait,qt);

      // Makes an additional irrelevant covariate
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 2);

      addlist(&mother->qttrait,qt);

      addlist (&finaldata,mother);
      mother->globalid = ++persons;
      mother->localid = persons;


      // Should now generate the children. 

      // Generates the offspring
      for (j=0; j<nooff; j++) {
	ind = (individual *) cmalloc (sizeof (individual));
	memset (ind,0,sizeof (individual));

	diseasecopies = 0;

	sprintf(ind->pedigree, "F%d",i);
	sprintf(ind->id, "%d-%d",i,j+3);

	if (genunf(0,1)<.5)
	  ind->sex = S_MALE;
	else
	  ind->sex = S_FEMALE;

	ind->father = father;
	ind->mother = mother;

	adduniqueidlist (&ind->father->offspring, ind);
	adduniqueidlist (&ind->father->mate, ind->mother);
	adduniqueidlist (&ind->mother->offspring, ind);
	adduniqueidlist (&ind->mother->mate,ind->father); 

	// tmpstr 1 and 2 hold information about which haplotype
	// the person received last time. Starts random
	sprintf(ind->tmpstr1, "%d", (int) ignbin(1,.5));
	sprintf(ind->tmpstr2, "%d", (int) ignbin(1,.5));

	//
	// Figure out the group
	//
	// Assumes that there is the same distance between markers
	  gamma = 0;
	  tal = genunf(0,1);
	  // Formula (3.2) in dissertation
	  gamma11 = 0.5*(2*recomb - InverseMapFunction(2*MapFunction(recomb)));
	  gamma10 = recomb - gamma11;
	  if (tal < gamma11) {
	    gamma = 11;
	  }
	  else if (tal < gamma11+gamma10 ) {
	    gamma = 10;
	  }	  
	  else if (tal < gamma11+2*gamma10 ) {
	    // Have here assumed that gamma01 = gamma10 
	    gamma = 1;
	  }

	  mamma = 0;
	  tal = genunf(0,1);
	  if (tal < gamma11) {
	    mamma = 11;
	  }
	  else if (tal < gamma11+gamma10 ) {
	    mamma = 10;
	  }
	  else if (tal < gamma11+2*gamma10 ) {
	    mamma = 1;
	  }

	  //	  printf("%2d   %2d\n", gamma, mamma);

	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {

	      marker = (markerlist *) cmalloc (sizeof (markerlist));
	      memset (marker,0,sizeof (markerlist));

	      fmkr = markernumber(father, order[mkrnum]);
	      mmkr = markernumber(mother, order[mkrnum]);

	      if (mkrnum>1) {
		
		if (mkrnum==closestmarker+1)
		  recomb = InverseMapFunction(fabs(MarkerDistance(location, mkrnum))/100);
		else
		  recomb = distance[mkrnum-1];
		
		// Checks for recombination for father
		if (mkrnum == 2) {
		  if (gamma==11 || gamma==10)
		    sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
		  // and mother
		  if (mamma==11 || mamma==10)
		    sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));
		}
		if (mkrnum == 3) {
		  if (gamma==11 || gamma==1)
		    sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
		  // and mother
		  if (mamma==11 || mamma==1)
		    sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));
		}
	      }

	      // Give the person the allele
	      if (!strcmpl(ind->tmpstr1,"0"))
		marker->allele1 = fmkr->allele1;
	      else
		marker->allele1 = fmkr->allele2;

	      if (!strcmpl(ind->tmpstr2,"0"))
		marker->allele2 = mmkr->allele1;
	      else
		marker->allele2 = mmkr->allele2;

	      addlist(&ind->marker, marker);


	      // Adds disease genes
	      if (mkrnum == closestmarker)
		{
		  // If no recombination between marker and disease QTL
		  if (genunf(0,1)>=dgenefreq)
		    {
		      if (!strcmpl(ind->tmpstr1,"0") && strlen(ind->father->tmpstr1))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr1,"1") && strlen(ind->father->tmpstr2))
			diseasecopies++;
		    }
		  else
		    {
		      if (!strcmpl(ind->tmpstr1,"0") && strlen(ind->father->tmpstr2))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr1,"1") && strlen(ind->father->tmpstr1))
			diseasecopies++;
		      // Makes a recombination
		      sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
		    }


		  // If no recombination between marker and disease QTL
		  if (genunf(0,1)>=dgenefreq)
		    {
		      if (!strcmpl(ind->tmpstr2,"0") && strlen(ind->mother->tmpstr1))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr2,"1") && strlen(ind->mother->tmpstr2))
			diseasecopies++;
		    }
		  else
		    {
		      if (!strcmpl(ind->tmpstr2,"0") && strlen(ind->mother->tmpstr2))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr2,"1") && strlen(ind->mother->tmpstr1))
			diseasecopies++;
		      // Makes a recombination
		      sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));
		    }

		}

	    }

	  nogenes = diseasecopies;
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0 + addeffect*(nogenes-1), 1);
	  if (nogenes==1)
	    qt->value += domeffect;

	  addlist(&ind->qttrait,qt);

	  // Makes an additional irrelevant covariate
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0, 2);
	  addlist(&ind->qttrait,qt);

	  addlist (&finaldata,ind);
	  ind->globalid = ++persons;
	  ind->localid = persons;

	}

    }


  // Fixing sibs
  for (ind = finaldata; ind; ind = ind->next)
    {
      for (ind2 = ind->next; ind2; ind2 = ind2->next)
	{
	  // If same father or mother then add to list of both
	  if (((ind->father == ind2->father) && (ind->father)) || 
	      ((ind->mother == ind2->mother) && (ind->mother)))
	    {
	      adduniqueidlist(&ind->sib, ind2);
	      adduniqueidlist(&ind2->sib, ind);
	    }
	}
    }

  // Needs to reorder the markers back for non-founders
  // so everything else works fine

  for (ind = finaldata; ind; ind = ind->next)
    {
      if (!founder(ind))
	{
	  i = 2;
	  memset(fixdata, 0, sizeof(fixdata));
	  fixdata[1] = ind->marker;
	  for (mkr = ind->marker; mkr; mkr = mkr->next)
	    {
	      fixdata[i] = mkr->next;
	      i++;
	    }

	  ind->marker = fixdata[invorder[1]];
	  for (i = 1; i<= numberofmarkers(); i++)
	    {
	      mkr = fixdata[i];
	      if (order[i] < numberofmarkers())
		mkr->next = fixdata[invorder[order[i]+1]];
	      else
		mkr->next = 0;
	    }
	}
    }

  // Then add trait names
  // Not necessary to add marker names as they've been read from the
  // parameter file already
  for (i=1; i<=2; i++)
    {
      nl = (namelist *) cmalloc (sizeof (namelist));
      memset (nl,0,sizeof (namelist));

      sprintf(nl->name, "Trait %d",i);
      addlist(&traitnames, nl);
    }



#ifdef fullyinformative
  printf("Simulated %d fully informative families with %d offspring each\n", families, nooff);
#else
  printf("Simulated %d families with %d offspring each\n", families, nooff);
#endif
  printf("The heritability of the simulated QTL is: %5.3f (Add. effect: %5.3f, dom. dev. effect: %5.3f) positioned at %5.3f\n", heritab / (1 + heritab), addeffect, domeffect, location);

  // Fixes the order etc.

  return finaldata; 
}

/*

// Assumes order is ok. Should be fixed
// families each with nooff offspring
// A single QTL a position location with an effect of effect
//
// The QTL must be inside the chromosome or else it won't work
individual *SimulateCEPHChromosome(int families, int nooff, double location, double addeffect, double domeffect, double prev, double resadd, double resdom)
{
  int i, j, nogenes, count, mkrnum, NMarkers;
  individual *ind, *ind2, *father1, *mother1, *father2, *mother2, *finaldata;
  markerlist *marker, *fmkr, *mmkr, *mkr;
  quanttrait *qt;
  namelist *nl;
  double randomnumber, res, recomb;
  markerlist *fixdata[MAXMARKERS+1];
  float dist, mindist, dgenefreq, heritab, tmp;
  int closestmarker, diseasecopies, persons;

  if (!allfreq) {
    printf("WARNING: No marker allele frequency information avaliable to simulate data.\n         Remember to load a parameter file first\n");
    return 0;
  }

  persons = 0;
  finaldata = 0;

  tmp = addeffect + domeffect*(1-2*prev);
  heritab = 2*prev*(1-prev)*SQR(tmp);
  tmp = 2*prev*(1-prev)*domeffect;
  heritab += SQR(tmp);

  // finds the marker closest to the left of the position
  mindist = 100000;
  dist = 0.0;
  closestmarker = 1;
  for (i=1; i<numberofmarkers(); i++)  {
    dist += MapFunction(distance[i])*100;
    if (location > dist)
      closestmarker = i+1;
  }
  // dgenefreq holds the recomb fraction between left marker and QTL
  dgenefreq = InverseMapFunction(fabs(MarkerDistance(location, closestmarker))/100);

  for (i=1; i<=families; i++)  {
      // Generates the founders (i.e., first generation)

      // First the first father
      father1 = (individual *) cmalloc (sizeof (individual));
      memset (father1,0,sizeof (individual));

      sprintf(father1->pedigree, "F%d",i);
      sprintf(father1->id, "%d-1",i);
      father1->sex = S_MALE;

      // Gives the person one copy of the QTL gene with probability 2*prev(1-prev)
      // and two copies with probability prev*prev
      if (genunf(0,1)<prev)
	strcpy(father1->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(father1->tmpstr2,"1");

      NMarkers = listlen(allfreq);
      for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {
	marker = (markerlist *) cmalloc (sizeof (markerlist));
	memset (marker,0,sizeof (markerlist));
	
	// Gets the relevant allele frequency for the marker
	CopyAlleleFreq(mkrnum);
	
	count = 1;
	res = allelefreq[count];
	// Gives the father1 allele 1
	randomnumber = genunf(0,1);
	while(randomnumber>res)  {
	  count++;
	  res += allelefreq[count];
	}
	marker->allele1 = count;

	  // and allele 2
	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res)
	    {
	      count++;
	      res += allelefreq[count];
	    }
	  marker->allele2 = count;

#ifdef fullyinformative
	  marker->allele1 = 1;
	  marker->allele2 = 2;      
#endif

	  addlist(&father1->marker, marker);
	}

      // This qt represents the observed trait
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));

      nogenes = strlen(father1->tmpstr1)+strlen(father1->tmpstr2);

      qt->sysmiss = 0;
      qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
      if (nogenes==1)
	qt->value += domeffect;

      addlist(&father1->qttrait,qt);

      // Makes an additional irrelevant covariate
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 2);
      addlist(&father1->qttrait,qt);

      addlist (&finaldata,father1);
      father1->globalid = ++persons;
      father1->localid = persons;


      // Then the mother
      mother1 = (individual *) cmalloc (sizeof (individual));
      memset (mother1,0,sizeof (individual));

      sprintf(mother1->pedigree, "F%d",i);
      sprintf(mother1->id, "%d-2",i);
      mother1->sex = S_FEMALE;

      // Gives the person one copy of the QTL gene with probability prev
      // and two copies with probability prev*prev
      if (genunf(0,1)<prev)
	strcpy(mother1->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(mother1->tmpstr2,"1");


      for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {

	  // Gets the relevant allele frequency for the marker
	  CopyAlleleFreq(mkrnum);

	  marker = (markerlist *) cmalloc (sizeof (markerlist));
	  memset (marker,0,sizeof (markerlist));

	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res)
	    {
	      count++;
	      res += allelefreq[count];
	    }
	  marker->allele1 = count;

	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res)
	    {
	      count++;
	      res += allelefreq[count];
	    }
	  marker->allele2 = count;

#ifdef fullyinformative
	  marker->allele1 = 3;
	  marker->allele2 = 4;      
#endif


	  addlist(&mother1->marker, marker);
	}

      // This qt represents the observed trait
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));

      nogenes = strlen(mother1->tmpstr1)+strlen(mother1->tmpstr2);

      qt->sysmiss = 0;
      qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
      if (nogenes==1)
	qt->value += domeffect;


      addlist(&mother1->qttrait,qt);

      // Makes an additional irrelevant covariate
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 2);

      addlist(&mother1->qttrait,qt);

      addlist (&finaldata,mother1);
      mother1->globalid = ++persons;
      mother1->localid = persons;




      // Then the second set of founders
      father2 = (individual *) cmalloc (sizeof (individual));
      memset (father2,0,sizeof (individual));

      sprintf(father2->pedigree, "F%d",i);
      sprintf(father2->id, "%d-3",i);
      father2->sex = S_MALE;

      // Gives the person one copy of the QTL gene with probability 2*prev(1-prev)
      // and two copies with probability prev*prev
      if (genunf(0,1)<prev)
	strcpy(father2->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(father2->tmpstr2,"1");

      NMarkers = listlen(allfreq);
      for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {
	marker = (markerlist *) cmalloc (sizeof (markerlist));
	memset (marker,0,sizeof (markerlist));
	
	// Gets the relevant allele frequency for the marker
	CopyAlleleFreq(mkrnum);
	
	count = 1;
	res = allelefreq[count];
	// Gives the father2 allele 1
	randomnumber = genunf(0,1);
	while(randomnumber>res)  {
	  count++;
	  res += allelefreq[count];
	}
	marker->allele1 = count;

	// and allele 2
	count = 1;
	res = allelefreq[count];
	randomnumber = genunf(0,1);
	while(randomnumber>res) {
	  count++;
	  res += allelefreq[count];
	}
	marker->allele2 = count;

#ifdef fullyinformative
	marker->allele1 = 1;
	marker->allele2 = 2;      
#endif
	
	addlist(&father2->marker, marker);
      }

      // This qt represents the observed trait
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));

      nogenes = strlen(father2->tmpstr1)+strlen(father2->tmpstr2);

      qt->sysmiss = 0;
      qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
      if (nogenes==1)
	qt->value += domeffect;

      addlist(&father2->qttrait,qt);

      // Makes an additional irrelevant covariate
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 2);
      addlist(&father2->qttrait,qt);

      addlist (&finaldata,father2);
      father2->globalid = ++persons;
      father2->localid = persons;


      // Then the second mother
      mother2 = (individual *) cmalloc (sizeof (individual));
      memset (mother2,0,sizeof (individual));

      sprintf(mother2->pedigree, "F%d",i);
      sprintf(mother2->id, "%d-4",i);
      mother2->sex = S_FEMALE;

      // Gives the person one copy of the QTL gene with probability prev
      // and two copies with probability prev*prev
      if (genunf(0,1)<prev)
	strcpy(mother2->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(mother2->tmpstr2,"1");


      for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {

	  // Gets the relevant allele frequency for the marker
	  CopyAlleleFreq(mkrnum);

	  marker = (markerlist *) cmalloc (sizeof (markerlist));
	  memset (marker,0,sizeof (markerlist));

	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res) {
	    count++;
	    res += allelefreq[count];
	  }
	  marker->allele1 = count;
	  
	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res) {
	    count++;
	    res += allelefreq[count];
	  }
	  marker->allele2 = count;

#ifdef fullyinformative
	  marker->allele1 = 3;
	  marker->allele2 = 4;      
#endif


	  addlist(&mother2->marker, marker);
      }

      // This qt represents the observed trait
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));

      nogenes = strlen(mother2->tmpstr1)+strlen(mother2->tmpstr2);

      qt->sysmiss = 0;
      qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
      if (nogenes==1)
	qt->value += domeffect;


      addlist(&mother2->qttrait,qt);

      // Makes an additional irrelevant covariate
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 2);

      addlist(&mother2->qttrait,qt);

      addlist (&finaldata,mother2);


      mother2->globalid = ++persons;
      mother2->localid = persons;



      // Now the second generation is created











      // Should now generate the children. 

      // Generates the offspring
      for (j=0; j<nooff; j++)
	{
	  ind = (individual *) cmalloc (sizeof (individual));
	  memset (ind,0,sizeof (individual));

	  diseasecopies = 0;

	  sprintf(ind->pedigree, "F%d",i);
	  sprintf(ind->id, "%d-%d",i,j+5);

	  if (genunf(0,1)<.5)
	    ind->sex = S_MALE;
	  else
	    ind->sex = S_FEMALE;

	  ind->father = father;
	  ind->mother = mother;

	  adduniqueidlist (&ind->father->offspring, ind);
	  adduniqueidlist (&ind->father->mate, ind->mother);
	  adduniqueidlist (&ind->mother->offspring, ind);
	  adduniqueidlist (&ind->mother->mate,ind->father); 

	  // tmpstr 1 and 2 hold information about which haplotype
	  // the person received last time. Starts random
	  sprintf(ind->tmpstr1, "%d", (int) ignbin(1,.5));
	  sprintf(ind->tmpstr2, "%d", (int) ignbin(1,.5));

	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++)
	    {

	      marker = (markerlist *) cmalloc (sizeof (markerlist));
	      memset (marker,0,sizeof (markerlist));

	      fmkr = markernumber(father, order[mkrnum]);
	      mmkr = markernumber(mother, order[mkrnum]);

	      if (mkrnum>1)
		{

		  if (mkrnum==closestmarker+1)
		    recomb = InverseMapFunction(fabs(MarkerDistance(location, mkrnum))/100);
		  else
		    recomb = distance[mkrnum-1];

		  // Checks for recombination for father
		  if (genunf(0,1)<recomb)
		    sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
		  // and mother
		  if (genunf(0,1)<recomb)
		    sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));
		}
	      // Give the person the allele
	      if (!strcmpl(ind->tmpstr1,"0"))
		marker->allele1 = fmkr->allele1;
	      else
		marker->allele1 = fmkr->allele2;

	      if (!strcmpl(ind->tmpstr2,"0"))
		marker->allele2 = mmkr->allele1;
	      else
		marker->allele2 = mmkr->allele2;

	      addlist(&ind->marker, marker);


	      // Adds disease genes
	      if (mkrnum == closestmarker)
		{
		  // If no recombination between marker and disease QTL
		  if (genunf(0,1)>=dgenefreq)
		    {
		      if (!strcmpl(ind->tmpstr1,"0") && strlen(ind->father->tmpstr1))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr1,"1") && strlen(ind->father->tmpstr2))
			diseasecopies++;
		    }
		  else
		    {
		      if (!strcmpl(ind->tmpstr1,"0") && strlen(ind->father->tmpstr2))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr1,"1") && strlen(ind->father->tmpstr1))
			diseasecopies++;
		      // Makes a recombination
		      sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
		    }


		  // If no recombination between marker and disease QTL
		  if (genunf(0,1)>=dgenefreq)
		    {
		      if (!strcmpl(ind->tmpstr2,"0") && strlen(ind->mother->tmpstr1))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr2,"1") && strlen(ind->mother->tmpstr2))
			diseasecopies++;
		    }
		  else
		    {
		      if (!strcmpl(ind->tmpstr2,"0") && strlen(ind->mother->tmpstr2))
			diseasecopies++;
		      if (!strcmpl(ind->tmpstr2,"1") && strlen(ind->mother->tmpstr1))
			diseasecopies++;
		      // Makes a recombination
		      sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));
		    }

		}

	    }

	  nogenes = diseasecopies;
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0 + addeffect*(nogenes-1), 1);
	  if (nogenes==1)
	    qt->value += domeffect;

	  addlist(&ind->qttrait,qt);

	  // Makes an additional irrelevant covariate
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0, 2);
	  addlist(&ind->qttrait,qt);

	  addlist (&finaldata,ind);
	  ind->globalid = ++persons;
	  ind->localid = persons;

	}

    }


  // Fixing sibs
  for (ind = finaldata; ind; ind = ind->next)
    {
      for (ind2 = ind->next; ind2; ind2 = ind2->next)
	{
	  // If same father or mother then add to list of both
	  if (((ind->father == ind2->father) && (ind->father)) || 
	      ((ind->mother == ind2->mother) && (ind->mother)))
	    {
	      adduniqueidlist(&ind->sib, ind2);
	      adduniqueidlist(&ind2->sib, ind);
	    }
	}
    }

  // Needs to reorder the markers back for non-founders
  // so everything else works fine

  for (ind = finaldata; ind; ind = ind->next)
    {
      if (!founder(ind))
	{
	  i = 2;
	  memset(fixdata, 0, sizeof(fixdata));
	  fixdata[1] = ind->marker;
	  for (mkr = ind->marker; mkr; mkr = mkr->next)
	    {
	      fixdata[i] = mkr->next;
	      i++;
	    }

	  ind->marker = fixdata[invorder[1]];
	  for (i = 1; i<= numberofmarkers(); i++)
	    {
	      mkr = fixdata[i];
	      if (order[i] < numberofmarkers())
		mkr->next = fixdata[invorder[order[i]+1]];
	      else
		mkr->next = 0;
	    }
	}
    }

  // Then add trait names
  // Not necessary to add marker names as they've been read from the
  // parameter file already
  for (i=1; i<=2; i++)
    {
      nl = (namelist *) cmalloc (sizeof (namelist));
      memset (nl,0,sizeof (namelist));

      sprintf(nl->name, "Trait %d",i);
      addlist(&traitnames, nl);
    }



#ifdef fullyinformative
  printf("Simulated %d fully informative families with %d offspring each\n", families, nooff);
#else
  printf("Simulated %d families with %d offspring each\n", families, nooff);
#endif
  printf("The heritability of the simulated QTL is: %5.3f (Add. effect: %5.3f, dom. dev. effect: %5.3f) positioned at %5.3f\n", heritab / (1 + heritab), addeffect, domeffect, location);

  // Fixes the order etc.

  return finaldata; 
}



*/








/*

 Simulates the family I need for my analysis

 */


individual *SimulateSpecialChromosome(int families, int nooff, double location, double addeffect, double domeffect, double prev)
{
  int i, ii, j, nogenes, count, mkrnum, NMarkers;
  individual *ind, *ind2, *finaldata;
  individual *father1, *father2, *mother1, *mother2, *spfather = NULL;
  markerlist *marker, *fmkr, *mmkr, *mkr;
  quanttrait *qt;
  namelist *nl;
  double randomnumber, res, recomb;
  markerlist *fixdata[MAXMARKERS+1];
  float dist, mindist, dgenefreq, heritab, tmp;
  int closestmarker, diseasecopies, persons;
  int discopy1, discopy2;


  if (!allfreq)
    {
      printf("WARNING: No marker allele frequency information avaliable to simulate data.\n         Remember to load a parameter file first\n");
      return 0;
    }

  persons = 0;
  finaldata = 0;


  tmp = addeffect + domeffect*(1-2*prev);
  heritab = 2*prev*(1-prev)*SQR(tmp);
  tmp = 2*prev*(1-prev)*domeffect;
  heritab += SQR(tmp);

  // finds the marker closest to the left of the position
  mindist = 100000;
  dist = 0.0;
  closestmarker = 1;
  for (i=1; i<numberofmarkers(); i++)
    {
      dist += MapFunction(distance[i])*100;
      if (location > dist)
	closestmarker = i+1;
    }
  // dgenefreq holds the recomb fraction between left marker and QTL
  dgenefreq = InverseMapFunction(fabs(MarkerDistance(location, closestmarker))/100);

  for (i=1; i<=families; i++)
    {
      // Generate the 4 founders
      father1 = (individual *) cmalloc (sizeof (individual));
      memset (father1,0,sizeof (individual));

      sprintf(father1->pedigree, "F%d",i);
      sprintf(father1->id, "%d-1",i);
      father1->sex = S_MALE;

      // Gives the person one copy of the QTL gene with probability 2*prev(1-prev)
      // and two copies with probability prev*prev
      if (genunf(0,1)<prev)
	strcpy(father1->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(father1->tmpstr2,"1");

      mother1 = (individual *) cmalloc (sizeof (individual));
      memset (mother1,0,sizeof (individual));

      sprintf(mother1->pedigree, "F%d",i);
      sprintf(mother1->id, "%d-2",i);
      mother1->sex = S_FEMALE;

      if (genunf(0,1)<prev)
	strcpy(mother1->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(mother1->tmpstr2,"1");

      father2 = (individual *) cmalloc (sizeof (individual));
      memset (father2,0,sizeof (individual));

      sprintf(father2->pedigree, "F%d",i);
      sprintf(father2->id, "%d-3",i);
      father2->sex = S_MALE;

      if (genunf(0,1)<prev)
	strcpy(father2->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(father2->tmpstr2,"1");

      mother2 = (individual *) cmalloc (sizeof (individual));
      memset (mother2,0,sizeof (individual));

      sprintf(mother2->pedigree, "F%d",i);
      sprintf(mother2->id, "%d-4",i);
      mother2->sex = S_FEMALE;

      if (genunf(0,1)<prev)
	strcpy(mother2->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(mother2->tmpstr2,"1");

      // Have now made the 4 founders


      // Assigns markers and alleles to the founders

      for (ii = 1; ii<=4; ii++)
	{
	  NMarkers = listlen(allfreq);
	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++)
	    {
	      marker = (markerlist *) cmalloc (sizeof (markerlist));
	      memset (marker,0,sizeof (markerlist));

	      // Gets the relevant allele frequency for the marker
	      CopyAlleleFreq(mkrnum);

	      count = 1;
	      res = allelefreq[count];
	      // Gives the father allele 1
	      randomnumber = genunf(0,1);
	      while(randomnumber>res)
		{
		  count++;
		  res += allelefreq[count];
		}
	      marker->allele1 = count;

	      // and allele 2
	      count = 1;
	      res = allelefreq[count];
	      randomnumber = genunf(0,1);
	      while(randomnumber>res)
		{
		  count++;
		  res += allelefreq[count];
		}
	      marker->allele2 = count;

	      switch(ii)
		{
		case 1: addlist(&father1->marker, marker); break;
		case 2: addlist(&mother1->marker, marker); break;
		case 3: addlist(&father2->marker, marker); break;
		case 4: addlist(&mother2->marker, marker); break;
		}      
	    }

	  // This qt represents the observed trait
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));

	  switch(ii)
	    {
	    case 1: spfather = father1; break;
	    case 2: spfather = mother1; break;
	    case 3: spfather = father2; break;
	    case 4: spfather = mother2; break;
	    }

	  nogenes = strlen(spfather->tmpstr1)+strlen(spfather->tmpstr2);

	  qt->sysmiss = 0;
	  qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
	  if (nogenes==1)
	    qt->value += domeffect;

	  switch(ii)
	    {
	    case 1: addlist(&father1->qttrait, qt); break;
	    case 2: addlist(&mother1->qttrait, qt); break;
	    case 3: addlist(&father2->qttrait, qt); break;
	    case 4: addlist(&mother2->qttrait, qt); break;
	    }      

	  // Makes an additional irrelevant covariate
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0, 2);


	  switch(ii)
	    {
	    case 1: addlist(&father1->qttrait, qt); break;
	    case 2: addlist(&mother1->qttrait, qt); break;
	    case 3: addlist(&father2->qttrait, qt); break;
	    case 4: addlist(&mother2->qttrait, qt); break;
	    }      
	}

      addlist (&finaldata,father1);
      father1->globalid = ++persons;
      father1->localid = persons;

      addlist (&finaldata,mother1);
      mother1->globalid = ++persons;
      mother1->localid = persons;

      addlist (&finaldata,father2);
      father2->globalid = ++persons;
      father2->localid = persons;

      addlist (&finaldata,mother2);
      mother2->globalid = ++persons;
      mother2->localid = persons;

      // Should now generate the children. 

      // Generates the offspring
      for (j=0; j<5; j++)
	{

	  diseasecopies = 0;
	  discopy1 = 0;
	  discopy2 = 0;

	  ind = (individual *) cmalloc (sizeof (individual));
	  memset (ind,0,sizeof (individual));

	  if (j==3)
	    spfather = ind;

	  sprintf(ind->pedigree, "F%d",i);
	  sprintf(ind->id, "%d-%d",i,j+5);

	  ind->sex = S_MALE;

	  ind->mother = mother1;
	  if (j < 2)
	    ind->father = father1;
	  else
	    ind->father = father2;
	  if (j==4)
	    {
	      ind->father = spfather;
	      ind->mother = mother2;
	    }

	  adduniqueidlist (&ind->father->offspring, ind);
	  adduniqueidlist (&ind->father->mate, ind->mother);
	  adduniqueidlist (&ind->mother->offspring, ind);
	  adduniqueidlist (&ind->mother->mate,ind->father); 

	  // Randomly select which parental haplotype to start with
	  sprintf(ind->tmpstr1, "%d", (int) ignbin(1,.5));
	  sprintf(ind->tmpstr2, "%d", (int) ignbin(1,.5));

	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++)
	    {
	      marker = (markerlist *) cmalloc (sizeof (markerlist));
	      memset (marker,0,sizeof (markerlist));

	      fmkr = markernumber(ind->father, order[mkrnum]);
	      mmkr = markernumber(ind->mother, order[mkrnum]);

	      if (mkrnum>1)
		{

		  if (mkrnum==closestmarker+1)
		    recomb = InverseMapFunction(fabs(MarkerDistance(location, mkrnum))/100);
		  else 
		    recomb = distance[mkrnum-1];

		  // Checks for recombination for father
		  if (genunf(0,1)<=recomb)
		    sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
		  // and mother
		  if (genunf(0,1)<=recomb)
		    sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));

		}

	      // Give the person the correct allele for the marker
	      if (!strcmpl(ind->tmpstr1,"0"))
		marker->allele1 = fmkr->allele1;
	      else
		marker->allele1 = fmkr->allele2;

	      if (!strcmpl(ind->tmpstr2,"0"))
		marker->allele2 = mmkr->allele1;
	      else
		marker->allele2 = mmkr->allele2;

	      addlist(&ind->marker, marker);

	      // Adds disease genes
	      if (mkrnum == closestmarker)
		{
		  // If a recombination between marker and disease QTL for paternal and maternal haplotype
		  if (genunf(0,1)<=dgenefreq)
		    sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
		  if (genunf(0,1)<=dgenefreq)
		    sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));


		  if ( (!strcmpl(ind->tmpstr1,"0") && strlen(ind->father->tmpstr1)) || 
		       (!strcmpl(ind->tmpstr1,"1") && strlen(ind->father->tmpstr2)))
		    {
		      diseasecopies++;
		      discopy1 = 1;
		    }

		  if ( (!strcmpl(ind->tmpstr2,"0") && strlen(ind->mother->tmpstr1)) ||
		       (!strcmpl(ind->tmpstr2,"1") && strlen(ind->mother->tmpstr2)))
		    {
		      diseasecopies++;
		      discopy2 = 1;
		    }          
		}
	    }

	  nogenes = diseasecopies;
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0 + addeffect*(nogenes-1), 1);
	  if (nogenes==1)
	    qt->value += domeffect;

	  addlist(&ind->qttrait,qt);

	  // Makes an additional irrelevant covariate
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0, 2);
	  addlist(&ind->qttrait,qt);

	  addlist (&finaldata,ind);
	  ind->globalid = ++persons;
	  ind->localid = persons;

	  // Fixes the disease copies
	  if (discopy1>0)
	    sprintf(ind->tmpstr1, "1");
	  else
	    sprintf(ind->tmpstr1, "");

	  if (discopy2>0)
	    sprintf(ind->tmpstr2, "1");
	  else
	    sprintf(ind->tmpstr2, "");

	}
    }

  // Fixing sibs
  for (ind = finaldata; ind; ind = ind->next)
    {
      for (ind2 = ind->next; ind2; ind2 = ind2->next)
	{
	  // If same father or mother then add to list of both
	  if (((ind->father == ind2->father) && (ind->father)) || 
	      ((ind->mother == ind2->mother) && (ind->mother)))
	    {
	      adduniqueidlist(&ind->sib, ind2);
	      adduniqueidlist(&ind2->sib, ind);
	    }
	}
    }

  // Needs to reorder the markers back for non-founders
  // so everything else works fine

  for (ind = finaldata; ind; ind = ind->next)
    {
      if (!founder(ind))
	{
	  i = 2;
	  memset(fixdata, 0, sizeof(fixdata));
	  fixdata[1] = ind->marker;
	  for (mkr = ind->marker; mkr; mkr = mkr->next)
	    {
	      fixdata[i] = mkr->next;
	      i++;
	    }

	  ind->marker = fixdata[invorder[1]];
	  for (i = 1; i<= numberofmarkers(); i++)
	    {
	      mkr = fixdata[i];
	      if (order[i] < numberofmarkers())
		mkr->next = fixdata[invorder[order[i]+1]];
	      else
		mkr->next = 0;
	    }
	}
    }

  // Then add trait names
  // Not necessary to add marker names as they've been read from the
  // parameter file already
  for (i=1; i<=2; i++)
    {
      nl = (namelist *) cmalloc (sizeof (namelist));
      memset (nl,0,sizeof (namelist));

      sprintf(nl->name, "Trait %d",i);
      addlist(&traitnames, nl);
    }

  //  #ifdef fullyinformative
  //  printf("Simulated %d fully informative special families\n", families);
  //  #else
  printf("Simulated %d special families\n", families);
  //  #endif

  printf("The heritability of the simulated QTL is: %5.3f (Add. effect: %5.3f, dom. dev. effect: %5.3f)\n", heritab / (1 + heritab), addeffect, domeffect);

  // Fixes the order etc.

  return finaldata; 
}


individual *SimulateSpecialChromosome2(int families, int nooff, double location, double addeffect, double domeffect, double prev)
{
  int i, ii, j, nogenes, count, mkrnum, NMarkers;
  individual *ind, *ind2, *finaldata;
  individual *father1, *mother1, *mother2, *spfather = NULL;
  markerlist *marker, *fmkr, *mmkr, *mkr;
  quanttrait *qt;
  namelist *nl;
  double randomnumber, res, recomb;
  markerlist *fixdata[MAXMARKERS+1];
  float dist, mindist, dgenefreq, heritab, tmp;
  int closestmarker, diseasecopies, persons;
  int discopy1, discopy2;


  if (!allfreq)
    {
      printf("WARNING: No marker allele frequency information avaliable to simulate data.\n         Remember to load a parameter file first\n");
      return 0;
    }

  persons = 0;
  finaldata = 0;


  tmp = addeffect + domeffect*(1-2*prev);
  heritab = 2*prev*(1-prev)*SQR(tmp);
  tmp = 2*prev*(1-prev)*domeffect;
  heritab += SQR(tmp);

  // finds the marker closest to the left of the position
  mindist = 100000;
  dist = 0.0;
  closestmarker = 1;
  for (i=1; i<numberofmarkers(); i++)
    {
      dist += MapFunction(distance[i])*100;
      if (location > dist)
	closestmarker = i+1;
    }
  // dgenefreq holds the recomb fraction between left marker and QTL
  dgenefreq = InverseMapFunction(fabs(MarkerDistance(location, closestmarker))/100);

  for (i=1; i<=families; i++)
    {
      // Generate the 3 founders
      father1 = (individual *) cmalloc (sizeof (individual));
      memset (father1,0,sizeof (individual));

      sprintf(father1->pedigree, "F%d",i);
      sprintf(father1->id, "%d-1",i);
      father1->sex = S_MALE;

      // Gives the person one copy of the QTL gene with probability 2*prev(1-prev)
      // and two copies with probability prev*prev
      if (genunf(0,1)<prev)
	strcpy(father1->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(father1->tmpstr2,"1");

      mother1 = (individual *) cmalloc (sizeof (individual));
      memset (mother1,0,sizeof (individual));

      sprintf(mother1->pedigree, "F%d",i);
      sprintf(mother1->id, "%d-2",i);
      mother1->sex = S_FEMALE;

      if (genunf(0,1)<prev)
	strcpy(mother1->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(mother1->tmpstr2,"1");

      mother2 = (individual *) cmalloc (sizeof (individual));
      memset (mother2,0,sizeof (individual));

      sprintf(mother2->pedigree, "F%d",i);
      sprintf(mother2->id, "%d-3",i);
      mother2->sex = S_FEMALE;

      if (genunf(0,1)<prev)
	strcpy(mother2->tmpstr1,"1");
      if (genunf(0,1)<prev)
	strcpy(mother2->tmpstr2,"1");

      // Have now made the 3 founders


      // Assigns markers and alleles to the founders

      for (ii = 1; ii<=3; ii++)
	{
	  NMarkers = listlen(allfreq);
	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++)
	    {
	      marker = (markerlist *) cmalloc (sizeof (markerlist));
	      memset (marker,0,sizeof (markerlist));

	      // Gets the relevant allele frequency for the marker
	      CopyAlleleFreq(mkrnum);

	      count = 1;
	      res = allelefreq[count];
	      // Gives the father allele 1
	      randomnumber = genunf(0,1);
	      while(randomnumber>res)
		{
		  count++;
		  res += allelefreq[count];
		}
	      marker->allele1 = count;

	      // and allele 2
	      count = 1;
	      res = allelefreq[count];
	      randomnumber = genunf(0,1);
	      while(randomnumber>res)
		{
		  count++;
		  res += allelefreq[count];
		}
	      marker->allele2 = count;

#ifdef fullyinformative
	      switch(ii)
		{
		case 1: marker->allele1 = 1; marker->allele2 = 2; break;
		case 2: marker->allele1 = 3; marker->allele2 = 4; break;
		case 3: marker->allele1 = 5; marker->allele2 = 2; break;
		}      
#endif

	      switch(ii)
		{
		case 1: addlist(&father1->marker, marker); break;
		case 2: addlist(&mother1->marker, marker); break;
		case 3: addlist(&mother2->marker, marker); break;
		}      
	    }

	  // This qt represents the observed trait
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));

	  switch(ii) {
	  case 1: spfather = father1; break;
	  case 2: spfather = mother1; break;
	  case 3: spfather = mother2; break;
	  }

	  nogenes = strlen(spfather->tmpstr1)+strlen(spfather->tmpstr2);

	  qt->sysmiss = 0;
	  qt->value = gennor(0 + (nogenes-1)*addeffect, 1);
	  if (nogenes==1)
	    qt->value += domeffect;

	  switch(ii)
	    {
	    case 1: addlist(&father1->qttrait, qt); break;
	    case 2: addlist(&mother1->qttrait, qt); break;
	    case 3: addlist(&mother2->qttrait, qt); break;
	    }      

	  // Makes an additional irrelevant covariate
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0, 2);


	  switch(ii)
	    {
	    case 1: addlist(&father1->qttrait, qt); break;
	    case 2: addlist(&mother1->qttrait, qt); break;
	    case 3: addlist(&mother2->qttrait, qt); break;
	    }      
	}

      addlist (&finaldata,father1);
      father1->globalid = ++persons;
      father1->localid = persons;

      addlist (&finaldata,mother1);
      mother1->globalid = ++persons;
      mother1->localid = persons;

      addlist (&finaldata,mother2);
      mother2->globalid = ++persons;
      mother2->localid = persons;

      // Should now generate the children. 

      // Generates the offspring
      for (j=0; j<5; j++)
	{

	  diseasecopies = 0;
	  discopy1 = 0;
	  discopy2 = 0;

	  ind = (individual *) cmalloc (sizeof (individual));
	  memset (ind,0,sizeof (individual));

	  if (j==3)
	    spfather = ind;

	  sprintf(ind->pedigree, "F%d",i);
	  sprintf(ind->id, "%d-%d",i,j+4);

	  ind->sex = S_MALE;

	  ind->mother = mother1;
	  ind->father = father1;

	  if (j==4)
	    {
	      ind->father = spfather;
	      ind->mother = mother2;
	    }

	  adduniqueidlist (&ind->father->offspring, ind);
	  adduniqueidlist (&ind->father->mate, ind->mother);
	  adduniqueidlist (&ind->mother->offspring, ind);
	  adduniqueidlist (&ind->mother->mate,ind->father); 

	  // Randomly select which parental haplotype to start with
	  sprintf(ind->tmpstr1, "%d", (int) ignbin(1,.5));
	  sprintf(ind->tmpstr2, "%d", (int) ignbin(1,.5));

	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++)
	    {
	      marker = (markerlist *) cmalloc (sizeof (markerlist));
	      memset (marker,0,sizeof (markerlist));

	      fmkr = markernumber(ind->father, order[mkrnum]);
	      mmkr = markernumber(ind->mother, order[mkrnum]);

	      if (mkrnum>1)
		{

		  if (mkrnum==closestmarker+1)
		    recomb = InverseMapFunction(fabs(MarkerDistance(location, mkrnum))/100);
		  else 
		    recomb = distance[mkrnum-1];

		  // Checks for recombination for father
		  if (genunf(0,1)<=recomb)
		    sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
		  // and mother
		  if (genunf(0,1)<=recomb)
		    sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));

		}

	      // Give the person the correct allele for the marker
	      if (!strcmpl(ind->tmpstr1,"0"))
		marker->allele1 = fmkr->allele1;
	      else
		marker->allele1 = fmkr->allele2;

	      if (!strcmpl(ind->tmpstr2,"0"))
		marker->allele2 = mmkr->allele1;
	      else
		marker->allele2 = mmkr->allele2;

	      addlist(&ind->marker, marker);

	      // Adds disease genes
	      if (mkrnum == closestmarker)
		{
		  // If a recombination between marker and disease QTL for paternal and maternal haplotype
		  if (genunf(0,1)<=dgenefreq)
		    sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
		  if (genunf(0,1)<=dgenefreq)
		    sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));


		  if ( (!strcmpl(ind->tmpstr1,"0") && strlen(ind->father->tmpstr1)) || 
		       (!strcmpl(ind->tmpstr1,"1") && strlen(ind->father->tmpstr2)))
		    {
		      diseasecopies++;
		      discopy1 = 1;
		    }

		  if ( (!strcmpl(ind->tmpstr2,"0") && strlen(ind->mother->tmpstr1)) ||
		       (!strcmpl(ind->tmpstr2,"1") && strlen(ind->mother->tmpstr2)))
		    {
		      diseasecopies++;
		      discopy2 = 1;
		    }          
		}
	    }

	  nogenes = diseasecopies;
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0 + addeffect*(nogenes-1), 1);
	  if (nogenes==1)
	    qt->value += domeffect;

	  addlist(&ind->qttrait,qt);

	  // Makes an additional irrelevant covariate
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = gennor(0, 2);
	  addlist(&ind->qttrait,qt);

	  addlist (&finaldata,ind);
	  ind->globalid = ++persons;
	  ind->localid = persons;

	  // Fixes the disease copies
	  if (discopy1>0)
	    sprintf(ind->tmpstr1, "1");
	  else
	    sprintf(ind->tmpstr1, "");

	  if (discopy2>0)
	    sprintf(ind->tmpstr2, "1");
	  else
	    sprintf(ind->tmpstr2, "");

	}
    }

  // Fixing sibs
  for (ind = finaldata; ind; ind = ind->next)
    {
      for (ind2 = ind->next; ind2; ind2 = ind2->next)
	{
	  // If same father or mother then add to list of both
	  if (((ind->father == ind2->father) && (ind->father)) || 
	      ((ind->mother == ind2->mother) && (ind->mother)))
	    {
	      adduniqueidlist(&ind->sib, ind2);
	      adduniqueidlist(&ind2->sib, ind);
	    }
	}
    }

  // Needs to reorder the markers back for non-founders
  // so everything else works fine

  for (ind = finaldata; ind; ind = ind->next)
    {
      if (!founder(ind))
	{
	  i = 2;
	  memset(fixdata, 0, sizeof(fixdata));
	  fixdata[1] = ind->marker;
	  for (mkr = ind->marker; mkr; mkr = mkr->next)
	    {
	      fixdata[i] = mkr->next;
	      i++;
	    }

	  ind->marker = fixdata[invorder[1]];
	  for (i = 1; i<= numberofmarkers(); i++)
	    {
	      mkr = fixdata[i];
	      if (order[i] < numberofmarkers())
		mkr->next = fixdata[invorder[order[i]+1]];
	      else
		mkr->next = 0;
	    }
	}
    }

  // Then add trait names
  // Not necessary to add marker names as they've been read from the
  // parameter file already
  for (i=1; i<=2; i++)
    {
      nl = (namelist *) cmalloc (sizeof (namelist));
      memset (nl,0,sizeof (namelist));

      sprintf(nl->name, "Trait %d",i);
      addlist(&traitnames, nl);
    }

#ifdef fullyinformative
  printf("Simulated %d fully informative special families\n", families);
#else
  printf("Simulated %d special families\n", families);
#endif

  printf("The heritability of the simulated QTL is: %5.3f (Add. effect: %5.3f, dom. dev. effect: %5.3f)\n", heritab / (1 + heritab), addeffect, domeffect);

  // Fixes the order etc.

  return finaldata; 
}


// ---------------------------------------------------------------
// Simulates genotypes and phenotypes for a dataset in memory
// Removes old genotyping data, uses allele frequencies etc.
//
// Sets the persons with missing data to missings
// ---------------------------------------------------------------

individual *SimulateDataset(individual *indlist, double location, double addeffect, double domeffect, double prev, double resadd, double resdom)
{
  int i, j, nogenes, count, mkrnum, NMarkers, nTraits;
  individual *ind, *minilist;
  namelist *ped, *Pedigrees;
  IDlist *idlist, *sortlist;
  markerlist *marker, *fmkr, *mmkr, *mkr;
  quanttrait *qt;
  double randomnumber, res, recomb;
  markerlist *fixdata[MAXMARKERS+1];
  float dist, mindist, dgenefreq, heritab, tmp;
  int closestmarker, diseasecopies, families;
  int discopy1, discopy2;

#ifdef seed
  long seed1, seed2;
  char phrase[81];
#endif

  if (!allfreq) {
    printf("WARNING: No marker allele frequency information avaliable to simulate data.\n");
    printf("Remember to load a parameter file first\n");
    return 0;
  }
  
  NMarkers = listlen(allfreq);
  nTraits  = max(numberoftraits(),2); // Generate at least two traits
  sortlist = SortIndividuals(indlist);

  tmp = addeffect + domeffect*(1-2*prev);
  heritab = 2*prev*(1-prev)*SQR(tmp);
  tmp = 2*prev*(1-prev)*domeffect;
  heritab += SQR(tmp) + resadd + resdom;

  // finds the marker closest to the left of the position
  mindist = 100000;
  dist = 0.0;
  closestmarker = 1;
  for (i=1; i<numberofmarkers(); i++) {
    dist += MapFunction(distance[i])*100;
    if (location > dist)
      closestmarker = i+1;
  }
  // dgenefreq holds the recomb fraction between left marker and QTL
  dgenefreq = InverseMapFunction(fabs(MarkerDistance(location, closestmarker))/100);


  families = 1;
  // Go through the sorted list
  for (idlist = sortlist; idlist; idlist = idlist->next)
    {
      // tmpstr holds the markers typed
      strcpy(idlist->ind->tmpstr1, "");

      for (i=1; i<=NMarkers; i++)
	{
	  if (istypedformarker(idlist->ind, i))
	    strcat(idlist->ind->tmpstr1, "1");
	  else
	    strcat(idlist->ind->tmpstr1, "0");
	}

      // tmpstr2 holds the phenotypes present
      strcpy(idlist->ind->tmpstr2, "");
      for (i=1; i<=nTraits; i++)
	{
	  if (!traitmiss(idlist->ind, i))
	    strcat(idlist->ind->tmpstr2, "1");
	  else
	    strcat(idlist->ind->tmpstr2, "0");
	}


      // Start by removing the old genotypes and traits ...
      if (idlist->ind->marker)
	{
	  freelist(idlist->ind->marker);
	  idlist->ind->marker=0;
	}

      if (idlist->ind->qttrait)
	{
	  freelist(idlist->ind->qttrait);
	  idlist->ind->qttrait = 0;
	}

      if (idlist->next && strcmpl(idlist->ind->pedigree, idlist->next->ind->pedigree))
	families++;

      // ... and clearing the tmp variables
      // tmpint holds phase-info or QTL info
      idlist->ind->tmpint1 = 0;
      idlist->ind->tmpint2 = 0;

#ifdef seed
      //  sprintf(phrase,"%ux%lx%f%x%u", getpid(), (long) time(NULL), cos((double) time(NULL)), getpid(), getuid());
      //  phrtsd(phrase, &seed1, &seed2);
      //  setall(seed1,seed2);  
#endif

      // if the person is a founder
      if (founder(idlist->ind))
	{
	  // Give the person one copy of the QTL gene with probability 2*prev(1-prev)
	  // and two copies with probability prev*prev
	  if (genunf(0,1)<prev)
	    idlist->ind->tmpint1 = 1;
	  if (genunf(0,1)<prev)
	    idlist->ind->tmpint2 = 1;

#ifdef foundershomoz
	  if (idlist->ind->sex==S_MALE) {
	    idlist->ind->tmpint1 = 1;
	    idlist->ind->tmpint2 = 1;
	  }
	  else if (idlist->ind->sex==S_FEMALE) {
	    idlist->ind->tmpint1 = 0;
	    idlist->ind->tmpint2 = 0;
	  }
#endif      

	  // Assigns markers
	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++)
	    {
	      marker = (markerlist *) cmalloc (sizeof (markerlist));
	      memset (marker,0,sizeof (markerlist));

	      // Gets the relevant allele frequency for the marker
	      CopyAlleleFreq(mkrnum);

	      count = 1;
	      res = allelefreq[count];
	      // Gives the father allele 1
	      randomnumber = genunf(0,1);
	      while(randomnumber>res)
		{
		  count++;
		  res += allelefreq[count];
		}
	      marker->allele1 = count;

	      // and allele 2
	      count = 1;
	      res = allelefreq[count];
	      randomnumber = genunf(0,1);

	      while(randomnumber>res) {
		  count++;
		  res += allelefreq[count];
		}
	      marker->allele2 = count;

#ifdef fullyinformative
	      if (idlist->ind->sex==S_MALE) {
		marker->allele1 = 1;
		marker->allele2 = 2;
	      }
	      else {
		marker->allele1 = 3;
		marker->allele2 = 4;
	      }
#endif

	      // Add the marker
	      addlist(&idlist->ind->marker, marker);
	    }

	  // Make a QTL-influenced trait
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  nogenes = idlist->ind->tmpint1+idlist->ind->tmpint2;

	  qt->sysmiss = 0;
	  qt->value = (nogenes-1)*addeffect;
	  if (nogenes==1)
	    qt->value += domeffect;

	  addlist(&idlist->ind->qttrait, qt);

	  // Makes an additional irrelevant covariate
	  for (i=2; i<=nTraits; i++) {
	    qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	    memset (qt,0,sizeof (quanttrait));
	    qt->sysmiss = 0;
	    qt->value = 0; 
	    addlist(&idlist->ind->qttrait, qt);
	  }
	}
      else  // If not a founder
	{
	  // Since parents appear first, they have been assigned new genotypes

	  diseasecopies = 0;
	  discopy1 = 0;
	  discopy2 = 0;

	  // Randomly select which parental haplotype to start with
	  idlist->ind->tmpint1 = ignbin(1,.5);
	  idlist->ind->tmpint2 = ignbin(1,.5);

	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++)
	    {
	      marker = (markerlist *) cmalloc (sizeof (markerlist));
	      memset (marker,0,sizeof (markerlist));

	      // If the parent is a founder then use the order to find
	      // the correct marker, otherwise use the same
	      if (founder(idlist->ind->father))
		fmkr = markernumber(idlist->ind->father, order[mkrnum]);
	      else
		fmkr = markernumber(idlist->ind->father, mkrnum);

	      if (founder(idlist->ind->mother))
		mmkr = markernumber(idlist->ind->mother, order[mkrnum]);
	      else
		mmkr = markernumber(idlist->ind->mother, mkrnum);

	      if (mkrnum>1)
		{

		  if (mkrnum==closestmarker+1)
		    recomb = InverseMapFunction(fabs(MarkerDistance(location, mkrnum))/100);
		  else 
		    recomb = distance[mkrnum-1];

		  // Checks for recombination for father
		  if (genunf(0,1)<=recomb)
		    idlist->ind->tmpint1 = 1 - idlist->ind->tmpint1;
		  // and mother
		  if (genunf(0,1)<=recomb)
		    idlist->ind->tmpint2 = 1 - idlist->ind->tmpint2;
		}

	      // Give the person the correct allele for the marker
	      if (idlist->ind->tmpint1 == 0)
		marker->allele1 = fmkr->allele1;
	      else
		marker->allele1 = fmkr->allele2;

	      if (idlist->ind->tmpint2 == 0)
		marker->allele2 = mmkr->allele1;
	      else
		marker->allele2 = mmkr->allele2;

	      addlist(&idlist->ind->marker, marker);

	      // Adds disease genes
	      if (mkrnum == closestmarker)
		{
		  // If a recombination between marker and disease QTL for paternal and maternal haplotype
		  if (genunf(0,1)<=dgenefreq)
		    idlist->ind->tmpint1 = 1 - idlist->ind->tmpint1;
		  if (genunf(0,1)<=dgenefreq)
		    idlist->ind->tmpint2 = 1 - idlist->ind->tmpint2;


		  if ( (idlist->ind->tmpint1 == 0 && idlist->ind->father->tmpint1 == 1) || 
		       (idlist->ind->tmpint1 == 1 && idlist->ind->father->tmpint2 == 1))
		    {
		      diseasecopies++;
		      discopy1 = 1;
		    }

		  if ( (idlist->ind->tmpint2 == 0 && idlist->ind->mother->tmpint1 == 1) || 
		       (idlist->ind->tmpint2 == 1 && idlist->ind->mother->tmpint2 == 1))
		    {
		      diseasecopies++;
		      discopy2 = 1;
		    }          
		}
	    }

	  nogenes = diseasecopies;
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = addeffect*(nogenes-1);
	  if (nogenes==1)
	    qt->value += domeffect;

	  addlist(&idlist->ind->qttrait,qt);

	  // Makes an additional irrelevant covariate
	  for (i=2; i<=nTraits; i++) {
	    qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	    memset (qt,0,sizeof (quanttrait));
	    qt->sysmiss = 0;
	    qt->value = 0; 
	    addlist(&idlist->ind->qttrait,qt);
	  }

	  // Fixes the disease copies
	  if (discopy1>0)
	    idlist->ind->tmpint1 = 1;
	  else
	    idlist->ind->tmpint1 = 0;

	  if (discopy2>0)
	    idlist->ind->tmpint2 = 1;
	  else
	    idlist->ind->tmpint2 = 0;
	}
    }


  // Needs to reorder the markers back for non-founders
  // so everything else works fine

  for (ind = indlist; ind; ind = ind->next)
    {
      if (!founder(ind))
	{
	  i = 2;
	  memset(fixdata, 0, sizeof(fixdata));
	  fixdata[1] = ind->marker;
	  for (mkr = ind->marker; mkr; mkr = mkr->next)
	    {
	      fixdata[i] = mkr->next;
	      i++;
	    }

	  ind->marker = fixdata[invorder[1]];
	  for (i = 1; i<= numberofmarkers(); i++)
	    {
	      mkr = fixdata[i];
	      if (order[i] < numberofmarkers())
		mkr->next = fixdata[invorder[order[i]+1]];
	      else
		mkr->next = 0;
	    }
	}
    }

  // Now clears marker and phenotype info as for the original dataset
  for (ind = indlist; ind; ind = ind->next)  {
    for (mkr = ind->marker, i=0; mkr; mkr = mkr->next, i++)  {
      if (ind->tmpstr1[i] == '0')
	{
	  mkr->allele1 = 0;
	  mkr->allele2 = 0;
	}	
    }
    for (qt = ind->qttrait, i=0; qt; qt = qt->next, i++)
      {      
	if (ind->tmpstr2[i] == '0')
	  {
	    qt->sysmiss = 1;
	    qt->value = 0;
	  }	
      }
  }

  // Also clears people with missing phenotypes
#ifdef seed
  sprintf(phrase,"%ux%lx%f%x%u", getpid(), (long) time(NULL), cos((double) time(NULL)), getpid(), getuid());
  phrtsd(phrase, &seed1, &seed2);
  setall(seed1,seed2);  
#endif

  // Finally adds variation by measurement error, resadd and resdom
  MATRIX Res;
  Pedigrees = MakePedigreeList(indlist);

  for (ped = Pedigrees; ped; ped = ped->next) {
    minilist = SelectPedigree(ped->name, indlist);

    Res.Resize(listlen(minilist), 1);

    MATRIX Error = MakeIdentMatrix(listlen(minilist));
    if (resadd>0) {
      MATRIX ResA = MMakeKinshipMatrix(minilist);
      Error += ResA*(resadd*2);
      
      // Check to see if residual dominance is used
      if (resdom>0) {
	MATRIX ResD = MMakeDelta7Matrix(minilist, ResA);
	Error += ResD*resdom; 
      }     
    }
    Error = Cholesky(Error);

    for (j = 1; j<=nTraits; j++) {
      // Generates standard normals
      for (i=1; i<=Res.Rows(); i++) {
	Res(i,1) = gennor(0,1);
      }

      Res = Error*Res;

      for (ind = minilist, i=1; ind; ind = ind->next, i++) {
	qt = traitvalue(FindIndListMemberFromPedigree(indlist,ind->id,ped->name), j);
	qt->value += Res(i,1);
      }
    }

    FreeIndividualList(minilist);
  }
  FreeNameList(Pedigrees);

#ifdef foundershomoz
  puts("WARNING: Simulating doubly homozygous parents !!! \n");
#endif      

#ifdef fullyinformative
  printf("Simulated %d individuals from %d pedigrees with fully informative parents\n", listlen(indlist), families);
#else
  printf("Simulated %d individuals from %d pedigrees\n", listlen(indlist), families);
#endif

  printf("Keeping untyped individuals untyped\n");

  freelist(sortlist);
  

  printf("  Disease prev.  : %5.3f\n", prev);
  printf("  QTL position   : %4.1f cM\n", location);
  printf("  QTL add. eff.  : %5.3f\n", addeffect);
  printf("  QTL dom. eff.  : %5.3f\n", domeffect);
  printf("  Res. add. eff. : %5.3f\n", resadd);
  printf("  Res. dom. eff. : %5.3f\n", resdom);

  printf("The (broad sense) heritability of the simulated QTL is: %5.3f\n", heritab / (1 + heritab));
  printf("The (narrow sense) heritability of the simulated QTL is: %5.3f\n", (heritab-resdom) / (1 + heritab));

  // Does actually not return anything as the original dataset has been modified
  return 0;

}



// ---------------------------------------------------------------
// Simulates genotypes and phenotypes for a dataset in memory
// Removes old genotyping data, uses allele frequencies etc.
//
// Sets the persons with missing data to missings
// ---------------------------------------------------------------

individual *SimulateXDataset(individual *indlist, double location, double addeffect, double domeffect, double prev, double resadd, double resdom)
{
  int i, j, nogenes, count, mkrnum, NMarkers, nTraits;
  individual *ind, *minilist;
  namelist *ped, *Pedigrees;
  IDlist *idlist, *sortlist;
  markerlist *marker, *fmkr, *mmkr, *mkr;
  quanttrait *qt;
  double randomnumber, res, recomb;
  markerlist *fixdata[MAXMARKERS+1];
  float dist, mindist, dgenefreq, heritab, tmp;
  int closestmarker, diseasecopies, families;
  int discopy1, discopy2;

#ifdef seed
  long seed1, seed2;
  char phrase[81];
#endif

  if (!allfreq) {
    printf("WARNING: No marker allele frequency information avaliable to simulate data.\n");
    printf("Remember to load a parameter file first\n");
    return 0;
  }

  // Check that we indeed have set the correct options
  if (!options->Index[O_XLINKED]) {
    printf("Trying to simulate X-linked data for an autosomal locus. I am confused.\n");
    return 0;    
  }
  
  NMarkers = listlen(allfreq);
  nTraits  = max(numberoftraits(),2); // Generate at least two traits
  sortlist = SortIndividuals(indlist);

  tmp = addeffect + domeffect*(1-2*prev);
  heritab = 2*prev*(1-prev)*SQR(tmp);
  tmp = 2*prev*(1-prev)*domeffect;
  heritab += SQR(tmp) + resadd + resdom;

  // finds the marker closest to the left of the position
  mindist = 100000;
  dist = 0.0;
  closestmarker = 1;
  for (i=1; i<numberofmarkers(); i++) {
    dist += MapFunction(distance[i])*100;
    if (location > dist)
      closestmarker = i+1;
  }
  // dgenefreq holds the recomb fraction between left marker and QTL
  dgenefreq = InverseMapFunction(fabs(MarkerDistance(location, closestmarker))/100);


  families = 1;
  // Go through the sorted list
  for (idlist = sortlist; idlist; idlist = idlist->next) {
      // tmpstr holds the markers typed
      strcpy(idlist->ind->tmpstr1, "");

      // Save a list of missing information. A 1 means typed and a 0 means missing
      for (i=1; i<=NMarkers; i++) {
	if (istypedformarker(idlist->ind, i))
	  strcat(idlist->ind->tmpstr1, "1");
	else
	  strcat(idlist->ind->tmpstr1, "0");
      }

      // Save the list of missing phenotypes
      // tmpstr2 holds the phenotypes present
      strcpy(idlist->ind->tmpstr2, "");
      for (i=1; i<=nTraits; i++) {
	if (!traitmiss(idlist->ind, i))
	  strcat(idlist->ind->tmpstr2, "1");
	else
	  strcat(idlist->ind->tmpstr2, "0");
      }


      // Start by removing the old genotypes and traits ...
      if (idlist->ind->marker) {
	freelist(idlist->ind->marker);
	idlist->ind->marker=0;
      }

      if (idlist->ind->qttrait)	{
	freelist(idlist->ind->qttrait);
	idlist->ind->qttrait = 0;
      }
      
      // Count the number of families
      if (idlist->next && strcmpl(idlist->ind->pedigree, idlist->next->ind->pedigree))
	families++;

      // ... and clearing the tmp variables
      // tmpint holds phase-info or QTL info
      idlist->ind->tmpint1 = 0;
      idlist->ind->tmpint2 = 0;

#ifdef seed
      //  sprintf(phrase,"%ux%lx%f%x%u", getpid(), (long) time(NULL), cos((double) time(NULL)), getpid(), getuid());
      //  phrtsd(phrase, &seed1, &seed2);
      //  setall(seed1,seed2);  
#endif

      // if the person is a founder
      if (founder(idlist->ind))	{
	  // Give males one copy of the QTL gene with probability prev
	  if (options->Index[O_XLINKED] && idlist->ind->sex == S_MALE) {
	    if (genunf(0,1)<prev) {
	      idlist->ind->tmpint1 = 1;
	      idlist->ind->tmpint2 = 1;
	    }
	  }
	  else {
	    // Give females one copy of the QTL gene with probability 2*prev(1-prev)
	    // and two copies with probability prev*prev
	    if (genunf(0,1)<prev)
	      idlist->ind->tmpint1 = 1;
	    if (genunf(0,1)<prev)
	      idlist->ind->tmpint2 = 1;
	  }

	  // Assigns markers
	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {
	    marker = (markerlist *) cmalloc (sizeof (markerlist));
	    memset (marker,0,sizeof (markerlist));

	    // Gets the relevant allele frequency for the marker
	    CopyAlleleFreq(mkrnum);


	    count = 1;
	    res = allelefreq[count];
	    // Gives the father allele 1
	    randomnumber = genunf(0,1);
	    while(randomnumber>res) {
	      count++;
	      res += allelefreq[count];
	    }
	    marker->allele1 = count;

	    // If a male then copy the allele otherwise generate a new for females
	    if (idlist->ind->sex == S_MALE) {
	      marker->allele2 = marker->allele1;
	    }
	    else {
	      // and allele 2
	      count = 1;
	      res = allelefreq[count];
	      randomnumber = genunf(0,1);

	      while(randomnumber>res) {
		count++;
		res += allelefreq[count];
	      }
	      marker->allele2 = count;
	    }

#ifdef fullyinformative
	      if (idlist->ind->sex==S_MALE) {
		marker->allele1 = 1;
		marker->allele2 = 1;
	      }
	      else {
		marker->allele1 = 3;
		marker->allele2 = 4;
	      }
#endif

	      // Add the marker
	      addlist(&idlist->ind->marker, marker);
	    }

	  // Make a QTL-influenced trait
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  nogenes = idlist->ind->tmpint1+idlist->ind->tmpint2;

	  qt->sysmiss = 0;

#ifdef DOSAGECOMPENSATION
	if (idlist->ind->sex == S_MALE) {
	  qt->value = addeffect*(idlist->ind->tmpint2);
	}
	else {
	  qt->value = 0.5*addeffect*(idlist->ind->tmpint1 + idlist->ind->tmpint2);
	}
#endif
#ifdef DOMINANCE
	if (idlist->ind->sex == S_MALE) {
	  qt->value = addeffect*(idlist->ind->tmpint2);
	}
	else {
	  qt->value = 0;
	  if (idlist->ind->tmpint1 + idlist->ind->tmpint2 >= 2) {
	    qt->value = addeffect;
	  }
	}
#endif
	addlist(&idlist->ind->qttrait, qt);

	// Makes an additional irrelevant covariate
	for (i=2; i<=nTraits; i++) {
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = 0; 
	  addlist(&idlist->ind->qttrait, qt);
	}
      }
      else { 
	// If not a founder
	// Since parents appear first, they have already been assigned new genotypes

	diseasecopies = 0;
	discopy1 = 0;
	discopy2 = 0;
	
	// Randomly select which parental haplotype to start with
	idlist->ind->tmpint1 = ignbin(1,.5);
	idlist->ind->tmpint2 = ignbin(1,.5);

	for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {
	  marker = (markerlist *) cmalloc (sizeof (markerlist));
	  memset (marker,0,sizeof (markerlist));
	  
	  // If the parent is a founder then use the order to find
	  // the correct marker, otherwise use the same
	  if (founder(idlist->ind->father))
	    fmkr = markernumber(idlist->ind->father, order[mkrnum]);
	  else
	    fmkr = markernumber(idlist->ind->father, mkrnum);
	  
	  if (founder(idlist->ind->mother))
	    mmkr = markernumber(idlist->ind->mother, order[mkrnum]);
	  else
	    mmkr = markernumber(idlist->ind->mother, mkrnum);
	  
	  if (mkrnum>1) {

	    if (mkrnum==closestmarker+1)
	      recomb = InverseMapFunction(fabs(MarkerDistance(location, mkrnum))/100);
	    else 
	      recomb = distance[mkrnum-1];

	    // Checks for recombination for father
	    if (genunf(0,1)<=recomb)
	      idlist->ind->tmpint1 = 1 - idlist->ind->tmpint1;
	    // and mother
	    if (genunf(0,1)<=recomb)
	      idlist->ind->tmpint2 = 1 - idlist->ind->tmpint2;
	  }

	  // Give the person the correct allele for the marker
	  if (idlist->ind->tmpint1 == 0)
	    marker->allele1 = fmkr->allele1;
	  else
	    marker->allele1 = fmkr->allele2;
	  
	  if (idlist->ind->tmpint2 == 0)
	    marker->allele2 = mmkr->allele1;
	  else
	    marker->allele2 = mmkr->allele2;

	  // Now fix males by setting all their paternal-inherited alleles equal to their mater
	  if (idlist->ind->sex == S_MALE) {
	    marker->allele1 = marker->allele2;
	  }
	  
	  addlist(&idlist->ind->marker, marker);

	  // Adds disease genes
	  if (mkrnum == closestmarker) {
	    // If a recombination between marker and disease QTL for paternal and maternal haplotype
	    if (genunf(0,1)<=dgenefreq)
	      idlist->ind->tmpint1 = 1 - idlist->ind->tmpint1;
	    if (genunf(0,1)<=dgenefreq)
	      idlist->ind->tmpint2 = 1 - idlist->ind->tmpint2;


	    // Only females inherit disease alleles from their father

	    if (idlist->ind->sex == S_FEMALE) {
	      // Always inherits the disease allele from the fathers mother
	      if (idlist->ind->father->tmpint2 == 1) {
		diseasecopies++;
		discopy1 = 1;
	      }
	    }
	    
	    // Disease allele from the mother
	    if ( (idlist->ind->tmpint2 == 0 && idlist->ind->mother->tmpint1 == 1) || 
		 (idlist->ind->tmpint2 == 1 && idlist->ind->mother->tmpint2 == 1)) {
		diseasecopies++;
		discopy2 = 1;
	      }          
	  }
	}


	//
	// Generates a QTL-trait
	//

	nogenes = diseasecopies;
	qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	memset (qt,0,sizeof (quanttrait));
	qt->sysmiss = 0;
	// To modify the X-linked model change the lines below.
	qt->value = addeffect*(nogenes-1);
	if (nogenes==1)
	  qt->value += domeffect;

#ifdef DOSAGECOMPENSATION
	if (idlist->ind->sex == S_MALE) {
	  qt->value = addeffect*(discopy2);
	}
	else {
	  qt->value = 0.5*addeffect*(discopy1 + discopy2);
	}
#endif
#ifdef DOMINANCE
	if (idlist->ind->sex == S_MALE) {
	  qt->value = addeffect*(discopy2);
	}
	else {
	  qt->value = 0;
	  if (diseasecopies >= 2) {
	    qt->value = addeffect;
	  }
	}
#endif

	
	addlist(&idlist->ind->qttrait,qt);
	
	// Makes an additional irrelevant covariate
	for (i=2; i<=nTraits; i++) {
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  // XXXX  qt->value = discopy2*10*addeffect; 
	  qt->value = 0; 
	  addlist(&idlist->ind->qttrait,qt);
	}
	
	// Fixes the disease copies
	if (discopy1>0)
	  idlist->ind->tmpint1 = 1;
	else
	  idlist->ind->tmpint1 = 0;
	
	if (discopy2>0)
	  idlist->ind->tmpint2 = 1;
	else
	  idlist->ind->tmpint2 = 0;
      }
  }


  // Needs to reorder the markers back for non-founders
  // so everything else works fine

  for (ind = indlist; ind; ind = ind->next) {
    if (!founder(ind)) {
      i = 2;
      memset(fixdata, 0, sizeof(fixdata));
      fixdata[1] = ind->marker;
      for (mkr = ind->marker; mkr; mkr = mkr->next) {
	fixdata[i] = mkr->next;
	i++;
      }
      
      ind->marker = fixdata[invorder[1]];
      for (i = 1; i<= numberofmarkers(); i++) {
	mkr = fixdata[i];
	if (order[i] < numberofmarkers())
	  mkr->next = fixdata[invorder[order[i]+1]];
	else
	  mkr->next = 0;
      }
    }
  }

  // Now clears marker and phenotype info as for the original dataset
  for (ind = indlist; ind; ind = ind->next)  {
    for (mkr = ind->marker, i=0; mkr; mkr = mkr->next, i++)  {
      if (ind->tmpstr1[i] == '0') {
	mkr->allele1 = 0;
	mkr->allele2 = 0;
      }	
    }
    for (qt = ind->qttrait, i=0; qt; qt = qt->next, i++) {      
      if (ind->tmpstr2[i] == '0') {
	qt->sysmiss = 1;
	qt->value = 0;
      }	
    }
  }

  // Also clears people with missing phenotypes
#ifdef seed
  sprintf(phrase,"%ux%lx%f%x%u", getpid(), (long) time(NULL), cos((double) time(NULL)), getpid(), getuid());
  phrtsd(phrase, &seed1, &seed2);
  setall(seed1,seed2);  
#endif

  // Finally adds variation by measurement error, resadd and resdom
  MATRIX Res;
  Pedigrees = MakePedigreeList(indlist);
  for (ped = Pedigrees; ped; ped = ped->next) {
    minilist = SelectPedigree(ped->name, indlist);

    Res.Resize(listlen(minilist), 1);

    MATRIX Error = MakeIdentMatrix(listlen(minilist));
    if (resadd>0) {
      MATRIX ResA = MMakeKinshipMatrix(minilist);
      Error += ResA*(resadd*2);

      // Check to see if residual dominance is used
      if (resdom>0) {
	MATRIX ResD = MMakeDelta7Matrix(minilist, ResA*2);
	Error += ResD*resdom; 
      }     
    }
    Error = Cholesky(Error);

    for (j = 1; j<=nTraits; j++) {
      // Generates standard normals
      for (i=1; i<=Res.Rows(); i++) {
	Res(i,1) = gennor(0,1);
      }

      Res = Error*Res;

      for (ind = minilist, i=1; ind; ind = ind->next, i++) {
	qt = traitvalue(FindIndListMember(indlist,ind->id), j);
	qt->value += Res(i,1);
      }
    }

    FreeIndividualList(minilist);
  }
  FreeNameList(Pedigrees);

#ifdef fullyinformative
  printf("Simulated %d individuals from %d pedigrees with fully informative parents\n", listlen(indlist), families);
#else
  printf("Simulated %d individuals from %d pedigrees\n", listlen(indlist), families);
#endif

  printf("Keeping untyped individuals untyped\n");
  printf("NOTE: simulated X-linked data\n");  
#ifdef DOSAGECOMPENSATION
  printf("NOTE: Dosage compensation model\n");
#endif
#ifdef DOMINANCE
  printf("NOTE: Dominance model\n");
#endif

  freelist(sortlist);
  

  printf("  Disease prev.  : %5.3f\n", prev);
  printf("  QTL position   : %4.1f cM\n", location);
  printf("  QTL add. eff.  : %5.3f\n", addeffect);
  printf("  QTL dom. eff.  : %5.3f\n", domeffect);
  printf("  Res. add. eff. : %5.3f\n", resadd);
  printf("  Res. dom. eff. : %5.3f\n", resdom);

  printf("The (broad sense) heritability of the simulated QTL is: %5.3f\n", heritab / (1 + heritab));
  printf("The (narrow sense) heritability of the simulated QTL is: %5.3f\n", (heritab-resdom) / (1 + heritab));

  // Does actually not return anything as the original dataset has been modified
  return 0;

}


// ---------------------------------------------------------------
// Simulates genotypes and phenotypes for a dataset in memory
// Meant for STENO power study
// Removes old genotyping data, uses allele frequencies etc.
//
// Sets the persons with missing data to missings
// ---------------------------------------------------------------

individual *SimulateStenoDataset(individual *indlist, double location,
				 double addeffect, double domeffect, double prev, double resadd, double resdom, int probandtrait) 
{
  int i, nogenes, count, mkrnum, NMarkers, nTraits, TraitsToMake;
  individual *ind, *minilist;
  namelist *ped, *Pedigrees, *oldped;
  IDlist *idlist, *sortlist;
  markerlist *marker, *fmkr, *mmkr, *mkr;
  quanttrait *qt;
  namelist *nl;
  double randomnumber, res, recomb, a = 0;
  markerlist *fixdata[MAXMARKERS+1];
  float dist, mindist, dgenefreq, heritab, tmp;
  int closestmarker, diseasecopies;
  int discopy1, discopy2;
  int codeok;
  double cutoff = 1.0;



  // Number of traits to make
  TraitsToMake = 1;

  if (!allfreq) {
    printf("WARNING: No marker allele frequency information avaliable to simulate data.\n         Remember to load a parameter file first\n");
    return 0;
  }

  NMarkers = listlen(allfreq);
  nTraits  = max(numberoftraits(),2); // Generate at least two traits
  sortlist = SortIndividuals(indlist);

  tmp = addeffect + domeffect*(1-2*prev);
  heritab = 2*prev*(1-prev)*SQR(tmp);
  tmp = 2*prev*(1-prev)*domeffect;
  heritab += SQR(tmp) + resadd + resdom;

  // finds the marker closest to the left of the position
  mindist = 100000;
  dist = 0.0;
  closestmarker = 1;
  for (i=1; i<numberofmarkers(); i++)
    {
      dist += MapFunction(distance[i])*100;
      if (location > dist)
	closestmarker = i+1;
    }
  // dgenefreq holds the recomb fraction between left marker and QTL
  dgenefreq = InverseMapFunction(fabs(MarkerDistance(location, closestmarker))/100);


  // Fixes traitnames
  freelist(traitnames);
  traitnames = 0;

  nl = (namelist *) cmalloc (sizeof (namelist));
  memset (nl,0,sizeof (namelist));
  sprintf(nl->name, "proband");
  addlist(&traitnames, nl);
  nl = (namelist *) cmalloc (sizeof (namelist));
  memset (nl,0,sizeof (namelist));
  sprintf(nl->name, "status");
  addlist(&traitnames, nl);
  nl = (namelist *) cmalloc (sizeof (namelist));
  memset (nl,0,sizeof (namelist));
  sprintf(nl->name, "noprob");
  addlist(&traitnames, nl);

  for (i = 1; i<=TraitsToMake ; i++) {
    nl = (namelist *) cmalloc (sizeof (namelist));
    memset (nl,0,sizeof (namelist));
    sprintf(nl->name, "trait%d",i);
    addlist(&traitnames, nl);
  }



  // Saves the proband info in tmpint3
  for (ind=indlist; ind; ind=ind->next) {
    ind->tmpint3 = (int) trait(ind, probandtrait);

    // tmpstr holds a 01-string of markers typed (i.e. missing genotypes)
    strcpy(ind->tmpstr1, "");
    for (i=1; i<=NMarkers; i++)  {
      if (istypedformarker(ind, i))
	strcat(ind->tmpstr1, "1");
      else
	strcat(ind->tmpstr1, "0");
    }

    // tmpstr2 holds the phenotypes present (i.e. missing phenotypes)
    strcpy(ind->tmpstr2, "");
    for (i=1; i<=nTraits; i++) {
      if (!traitmiss(ind, i))
	strcat(ind->tmpstr2, "1");
      else
	strcat(ind->tmpstr2, "0");
    }
  }

  // Go through each pedigree and make the data pedigree by pedigree
  Pedigrees = MakePedigreeList(indlist);

  for (ped = Pedigrees; ped;) {
    oldped = ped->next;
    // Go through the sorted list
    for (idlist = sortlist; idlist; idlist = idlist->next) {      
      if (strcmpl(idlist->ind->pedigree, ped->name))
	continue;

      // The individual belongs to the correct pedigree, so proceed


      // Start by removing the old genotypes and traits ...
      if (idlist->ind->marker) {
	freelist(idlist->ind->marker);
	idlist->ind->marker=0;
      }

      if (idlist->ind->qttrait) {
	// Free all but he first three traits. 
	// This is hardcoded and corresponds to STENO info
	freelist(idlist->ind->qttrait->next->next->next);
	idlist->ind->qttrait->next->next->next = 0;
      }

      // ... and clearing the tmp variables
      // tmpint holds phase-info or QTL info
      idlist->ind->tmpint1 = 0;
      idlist->ind->tmpint2 = 0;


      // Will now assign a QTL and markers


      // if the person is a founder
      if (founder(idlist->ind)) {
	// Give the person one copy of the QTL gene with probability 2*prev(1-prev)
	// and two copies with probability prev*prev
	if (genunf(0,1)<prev)
	  idlist->ind->tmpint1 = 1;
	if (genunf(0,1)<prev)
	  idlist->ind->tmpint2 = 1;

	// If it is a proband, then at least one allele
	// This is a dirty hack
	if (idlist->ind->tmpint3) {
	  if (idlist->ind->tmpint1)
	    idlist->ind->tmpint2 = 1;
	  else
	    idlist->ind->tmpint1 = 1;
	}

	// Assigns markers
	for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {
	  marker = (markerlist *) cmalloc (sizeof (markerlist));
	  memset (marker,0,sizeof (markerlist));

	  // Gets the relevant allele frequency for the marker
	  CopyAlleleFreq(mkrnum);

	  count = 1;
	  res = allelefreq[count];
	  // Gives the father allele 1
	  randomnumber = genunf(0,1);
	  while(randomnumber>res) {
	    count++;
	    res += allelefreq[count];
	  }
	  marker->allele1 = count;

	  // and allele 2
	  count = 1;
	  res = allelefreq[count];
	  randomnumber = genunf(0,1);
	  while(randomnumber>res) {
	    count++;
	    res += allelefreq[count];
	  }
	  marker->allele2 = count;

	  // Add the marker
	  addlist(&idlist->ind->marker, marker);
	}

	nogenes = idlist->ind->tmpint1 + idlist->ind->tmpint2;


	// Makes 66 different covariates
	for (i=0; i<TraitsToMake; i++) {
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = nogenes-1;

	  if (nogenes==1)
	    qt->value += domeffect;

	  addlist(&idlist->ind->qttrait,qt);
	}
      }
      else { // If not a founder
	// Since parents appear first, thet have been assigned new genotypes

	diseasecopies = 0;
	discopy1 = 0;
	discopy2 = 0;

	// Randomly select which parental haplotype to start with
	idlist->ind->tmpint1 = ignbin(1,.5);
	idlist->ind->tmpint2 = ignbin(1,.5);

	for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {
	  marker = (markerlist *) cmalloc (sizeof (markerlist));
	  memset (marker,0,sizeof (markerlist));

	  // If the parent is a founder then use the order to find
	  // the correct marker, otherwise use the same
	  if (founder(idlist->ind->father))
	    fmkr = markernumber(idlist->ind->father, order[mkrnum]);
	  else
	    fmkr = markernumber(idlist->ind->father, mkrnum);

	  if (founder(idlist->ind->mother))
	    mmkr = markernumber(idlist->ind->mother, order[mkrnum]);
	  else
	    mmkr = markernumber(idlist->ind->mother, mkrnum);

	  if (mkrnum>1) {

	    if (mkrnum==closestmarker+1)
	      recomb = InverseMapFunction(fabs(MarkerDistance(location, mkrnum))/100);
	    else 
	      recomb = distance[mkrnum-1];

	    // Checks for recombination for father
	    if (genunf(0,1)<=recomb)
	      idlist->ind->tmpint1 = 1 - idlist->ind->tmpint1;
	    // and mother
	    if (genunf(0,1)<=recomb)
	      idlist->ind->tmpint2 = 1 - idlist->ind->tmpint2;
	  }

	  // Give the person the correct allele for the marker
	  if (idlist->ind->tmpint1 == 0)
	    marker->allele1 = fmkr->allele1;
	  else
	    marker->allele1 = fmkr->allele2;

	  if (idlist->ind->tmpint2 == 0)
	    marker->allele2 = mmkr->allele1;
	  else
	    marker->allele2 = mmkr->allele2;

	  addlist(&idlist->ind->marker, marker);

	  // Adds disease genes
	  if (mkrnum == closestmarker){
	    // If a recombination between marker and disease QTL for paternal and maternal haplotype
	    if (genunf(0,1)<=dgenefreq)
	      idlist->ind->tmpint1 = 1 - idlist->ind->tmpint1;
	    if (genunf(0,1)<=dgenefreq)
	      idlist->ind->tmpint2 = 1 - idlist->ind->tmpint2;


	    if ( (idlist->ind->tmpint1 == 0 && idlist->ind->father->tmpint1 == 1) || 
		 (idlist->ind->tmpint1 == 1 && idlist->ind->father->tmpint2 == 1)) {
	      diseasecopies++;
	      discopy1 = 1;
	    }

	    if ( (idlist->ind->tmpint2 == 0 && idlist->ind->mother->tmpint1 == 1) || 
		 (idlist->ind->tmpint2 == 1 && idlist->ind->mother->tmpint2 == 1)) {
	      diseasecopies++;
	      discopy2 = 1;
	    }          
	  }
	}

	nogenes = diseasecopies;

	// Makes 66 different covariates
	for (i=0; i<TraitsToMake; i++) {
	  qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
	  memset (qt,0,sizeof (quanttrait));
	  qt->sysmiss = 0;
	  qt->value = nogenes-1;

	  if (nogenes==1)
	    qt->value += domeffect;

	  addlist(&idlist->ind->qttrait,qt);
	}

	// Fixes the disease copies
	if (discopy1>0)
	  idlist->ind->tmpint1 = 1;
	else
	  idlist->ind->tmpint1 = 0;

	if (discopy2>0)
	  idlist->ind->tmpint2 = 1;
	else
	  idlist->ind->tmpint2 = 0;
      }
    }  // End of sortlist


    MATRIX Res;

    minilist = SelectPedigree(ped->name, indlist);
    Res.Resize(listlen(minilist), 1);

    MATRIX ResA  = MMakeKinshipMatrix(minilist);

    // Calculates sigma_q**2
    a = 2.0*prev*(1-prev)*addeffect*addeffect;


    MATRIX Error = MakeIdentMatrix(listlen(minilist));
    Error        = Error + (ResA*(2 * resadd)); 

    Error = Cholesky(Error);


    // This part of the code runs until we have extreme probands
    // Could potentially run forever with bad choices of 
    // cutoff and 
    codeok = 1;
    for (i=1; i<=Res.Rows(); i++) {
      Res(i,1) = gennor(0,1);	  
    }
    Res = Error*Res;

    for (ind = minilist, i=1; ind; ind = ind->next, i++) {
      if (ind->tmpint3) { 
	// If it is a proband
	qt = traitvalue(FindIndListMember(indlist,ind->id), 1+3);
	if (qt->value*addeffect+Res(i,1) <= (1-2*prev)*addeffect) {
	  codeok = 0;	    	    
	}
      }
    }
    if (codeok)
      ped = oldped;

    for (ind = minilist, i=1; ind; ind = ind->next, i++) {
      qt = traitvalue(FindIndListMember(indlist,ind->id), 1+3);
      qt->value = qt->value*addeffect + Res(i,1);
    }
    FreeIndividualList(minilist);
  }

  // WORKZONE !!!


  // Needs to reorder the markers back for non-founders
  // so everything else works fine

  for (ind = indlist; ind; ind = ind->next)
    {
      if (!founder(ind))
	{
	  i = 2;
	  memset(fixdata, 0, sizeof(fixdata));
	  fixdata[1] = ind->marker;
	  for (mkr = ind->marker; mkr; mkr = mkr->next)
	    {
	      fixdata[i] = mkr->next;
	      i++;
	    }

	  ind->marker = fixdata[invorder[1]];
	  for (i = 1; i<= numberofmarkers(); i++)
	    {
	      mkr = fixdata[i];
	      if (order[i] < numberofmarkers())
		mkr->next = fixdata[invorder[order[i]+1]];
	      else
		mkr->next = 0;
	    }
	}
    }

  // Now clears marker and phenotype info as for the original dataset
  for (ind = indlist; ind; ind = ind->next) {
    for (mkr = ind->marker, i=0; mkr; mkr = mkr->next, i++) {
      if (ind->tmpstr1[i] == '0') {
	mkr->allele1 = 0;
	mkr->allele2 = 0;
      }	
    }
    for (qt = ind->qttrait, i=0; qt; qt = qt->next, i++) {      
      // Do not remove the first 4
      // Base removal on trait number 13 --- supposed to be gluc120
      if (ind->tmpstr2[12] == '0' && i>2) {
	qt->sysmiss = 1;
	qt->value = 0;
      }	
    }
  }

#ifdef fullyinformative
  printf("Simulated %d individuals from %d pedigrees with fully informative parents\n", listlen(indlist), listlen(Pedigrees));
#else
  printf("Simulated %d individuals from %d pedigrees\n", listlen(indlist), listlen(Pedigrees));
#endif

  printf("Keeping untyped individuals untyped\n");

  freelist(sortlist);
  FreeNameList(Pedigrees);

  printf("The (broad sense) heritability of the simulated QTL is: %5.3f (Add. effect: %5.3f, dom. dev. effect: %5.3f)\n", heritab / (1 + heritab), addeffect, domeffect);
  printf("  Cutoff used    : %5.3f\n", cutoff);
  printf("  Disease all frq: %5.3f\n", prev);
  printf("  Res. add. eff. : %5.3f\n", resadd);
  printf("  Res. dom. eff. : %5.3f\n", resdom);
  printf("  QTL %%          : %5.3f\n", a / (a+resadd+.000000001));

  // Does actually not return anything as the original dataset has been modified
  return 0;

}


// ---------------------------------------------------------------
// Simulates n markers of a specific type
// All markers assumed to be on the same chromosome
// All markers assumed to be in LE and HWE
// Using predetermined allele frequencies based on type
// Uses function SimulateGenotypes (and a dirty hack)
// If chromosomelength == 0 then same chromosome else length (in cM)
// ---------------------------------------------------------------

individual *SimulateMarkerSet(individual *indlist, int markertype, int nmarkers, double markerdist, double chromosomelength) {
  int i, oldn;
  FreqList *fq;
  individual *ind;
  namelist *nl;
  markerlist *marker;
  double currentlength;

  // Check that the #define constants are sufficiently large
  if (nmarkers>=NAMESIZE) {
    printf("WARNING: NAMESIZE constant too small to generate markers\n");
    return 0;
  }

  // Removes the old frequency information
  oldn = numberofmarkers()-1;
  InitializeFrequencies(nmarkers);

  switch (markertype) {
  case 1: printf("Simulating %d STR, 5 equifrequent alleles, %5.3f cM distance, chromosome is %5.3f cM long.\n", nmarkers, 100*markerdist, chromosomelength);
    break;
  case 2: printf("Simulating %d SNP, 2 equifrequent alleles, %5.3f cM distance, chromosome is %5.3f cM long\n", nmarkers, 100*markerdist, chromosomelength);
    break;
  default: printf("Illegal marker type (%d). Nothing simulated\n", markertype);
    return (NULL);
  }


  // Starts by generating the marker and allele frequency info
  for (fq = allfreq; fq; fq = fq->next) {
    switch (markertype) {
    case 1 : // Microsat. 5 alleles, 80% heterozyg
      fq->num_alleles = 5;
      for (i = 1; i<=5; i++) {
	fq->frequency[i] = .2;
      }
      break;
    case 2:  // SNP ... .5
      fq->num_alleles = 2;
      for (i = 1; i<=2; i++) {
	fq->frequency[i] = .5;
      }
      break;
    }
  }

  for (ind = indlist; ind; ind = ind->next) {
    freelist(ind->marker);
    ind->marker = NULL;

    for (i = 0; i< nmarkers; i++) {
      marker = (markerlist *) cmalloc (sizeof (markerlist));
      memset (marker,0,sizeof (markerlist));

      marker->allele1 = 1;
      marker->allele2 = 2;

      addlist(&ind->marker, marker);
    }
  }

  if (oldn>0) {
    free_ivector(order, 1, oldn);
    free_ivector(invorder, 1, oldn);
    if (oldn>1) {
      free_vector(distance,1,oldn-1);
    }
  }

  order = ivector(1,nmarkers);
  invorder = ivector(1,nmarkers);
  distance = vector(1, nmarkers-1);
  for (i=1; i<=nmarkers; i++)
    {
      order[i] = i;
      invorder[i] = i;
    }

  currentlength = 0 ;
  for (i=1; i<nmarkers; i++) {
    currentlength += 100*InverseMapFunction(markerdist);
    if (currentlength>chromosomelength && chromosomelength>0) {
      currentlength = 0;
      distance[i] = InverseMapFunction(.499);      
    }
    else {
      distance[i] = InverseMapFunction(markerdist);
    }
  }

  freelist(markernames);
  markernames = NULL;

  // First set the marker names:
  for (i=1; i<=nmarkers; i++) {
    nl = (namelist *) cmalloc (sizeof (namelist));
    memset (nl,0,sizeof (namelist));
    
    sprintf(nl->name, "Marker%d",i);
    addlist(&markernames, nl);
  }
  return(SimulateGenotypes(indlist));
}


// ---------------------------------------------------------------
// Simulates genotypes for a dataset in memory
// Removes old genotyping data, uses allele frequencies etc.
// Sets a person with missing genotypes to missing
// ---------------------------------------------------------------

individual *SimulateGenotypes(individual *indlist)
{
  int i, count, mkrnum, NMarkers;
  individual *ind;
  IDlist *idlist, *sortlist;
  markerlist *marker, *fmkr, *mmkr, *mkr;
  double randomnumber, res, recomb;
  markerlist *fixdata[MAXMARKERS+1];
  int families;


  if (!allfreq)
    {
      printf("WARNING: No marker allele frequency information avaliable to simulate data.\n         Remember to load a parameter file first\n");
      return 0;
    }

  NMarkers = listlen(allfreq);
  sortlist = SortIndividuals(indlist);

  families = 1;
  // Go through the sorted list
  for (idlist = sortlist; idlist; idlist = idlist->next)
    {
      // tmpstr holds the markers typed
      strcpy(idlist->ind->tmpstr1, "");

      for (i=1; i<=NMarkers; i++)
	{
	  if (istypedformarker(idlist->ind, i))
	    strcat(idlist->ind->tmpstr1, "1");
	  else
	    strcat(idlist->ind->tmpstr1, "0");
	}

      // Start by removing the old genotypes
      if (idlist->ind->marker) {
	  freelist(idlist->ind->marker);
	  idlist->ind->marker=0;
	}

      // Hmmm. I wonder what this is doing? Counting families? For what?
      if (idlist->next && strcmpl(idlist->ind->pedigree, idlist->next->ind->pedigree))
	families++;

      // if the person is a founder
      if (founder(idlist->ind))	{

	  // Assigns markers
	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {
	      marker = (markerlist *) cmalloc (sizeof (markerlist));
	      memset (marker,0,sizeof (markerlist));

	      // Gets the relevant allele frequency for the marker
	      CopyAlleleFreq(mkrnum);

	      count = 1;
	      res = allelefreq[count];
	      // Gives the father allele 1
	      randomnumber = genunf(0,1);
	      while(randomnumber>res) {
		  count++;
		  res += allelefreq[count];
		}
	      marker->allele1 = count;

	      // and allele 2
	      count = 1;
	      res = allelefreq[count];
	      randomnumber = genunf(0,1);
	      while(randomnumber>res) {
		count++;
		res += allelefreq[count];
	      }
	      marker->allele2 = count;

	      // Add the marker
	      addlist(&idlist->ind->marker, marker);
	    }
	}
      else  // If not a founder
	{
	  // Since parents appear first, thet have been assigned new genotypes

	  // Randomly select which parental haplotype to start with
	  idlist->ind->tmpint1 = ignbin(1,.5);
	  idlist->ind->tmpint2 = ignbin(1,.5);
      
	  for (mkrnum = 1; mkrnum <= NMarkers; mkrnum++) {
	      marker = (markerlist *) cmalloc (sizeof (markerlist));
	      memset (marker,0,sizeof (markerlist));

	      // If the parent is a founder then use the order to find
	      // the correct marker, otherwise use the same
	      if (founder(idlist->ind->father))
		fmkr = markernumber(idlist->ind->father, order[mkrnum]);
	      else
		fmkr = markernumber(idlist->ind->father, mkrnum);

	      if (founder(idlist->ind->mother))
		mmkr = markernumber(idlist->ind->mother, order[mkrnum]);
	      else
		mmkr = markernumber(idlist->ind->mother, mkrnum);

	      if (mkrnum>1) {
		  recomb = distance[mkrnum-1];

		  // Checks for recombination for father
		  if (genunf(0,1)<=recomb)
		    idlist->ind->tmpint1 = 1 - idlist->ind->tmpint1;
		  // and mother
		  if (genunf(0,1)<=recomb)
		    idlist->ind->tmpint2 = 1 - idlist->ind->tmpint2;
		}

	      // Give the person the correct allele for the marker
	      if (idlist->ind->tmpint1 == 0)
		marker->allele1 = fmkr->allele1;
	      else
		marker->allele1 = fmkr->allele2;

	      if (idlist->ind->tmpint2 == 0)
		marker->allele2 = mmkr->allele1;
	      else
		marker->allele2 = mmkr->allele2;

	      addlist(&idlist->ind->marker, marker);
	    }
	}
    }


  // Needs to reorder the markers back for non-founders
  // so everything else works fine

  for (ind = indlist; ind; ind = ind->next)
    {
      if (!founder(ind))
	{
	  i = 2;
	  memset(fixdata, 0, sizeof(fixdata));
	  fixdata[1] = ind->marker;
	  for (mkr = ind->marker; mkr; mkr = mkr->next)
	    {
	      fixdata[i] = mkr->next;
	      i++;
	    }

	  ind->marker = fixdata[invorder[1]];
	  for (i = 1; i<= numberofmarkers(); i++)
	    {
	      mkr = fixdata[i];
	      if (order[i] < numberofmarkers())
		mkr->next = fixdata[invorder[order[i]+1]];
	      else
		mkr->next = 0;
	    }
	}
    }

  // Now clears marker info as for the original dataset
  for (ind = indlist; ind; ind = ind->next)
    {
      for (mkr = ind->marker, i=0; mkr; mkr = mkr->next, i++)
	{
	  if (ind->tmpstr1[i] == '0')
	    {
	      mkr->allele1 = 0;
	      mkr->allele2 = 0;
	    }	
	}
    }

  printf("Simulated new genotypes for %d individuals from %d pedigrees\n", listlen(indlist), families);
  printf("Keeping untyped individuals untyped\n");

  freelist(sortlist);

  // Does actually not return anything as the original dataset has been modified
  return 0;

}



// -------------------------
//
// Genotyping error rutines
//
// -------------------------


//
// Introduce genotyping error for ind that are consistent with Mendelian inheritance
//
void GenotypeErrorMendelian(individual *ind, int markernum, int *changed)
{
  int allele1, allele2, rannum;
  markerlist *origmarker, *mkr;
  individual *elim, *ind2, *ind3, *ind4;


  *changed = 0;
  // Check that person is genotyped at all
  // If not then move on

  if (!istypedformarker(ind,markernum)) {
    return;
  }


  // Save the two alleles
  origmarker = markernumber(ind,markernum);
  allele1 = origmarker->allele1;
  allele2 = origmarker->allele2;

  // Set the original genotype to missing
  origmarker->allele1 = 0;
  origmarker->allele2 = 0;

  CopyAlleleFreq(markernum);

  // Do genotype elimination
  elim = (individual *) GenotypeElimination(ind->pedigree, markernum);

  // Find the individual 
  for (ind2 = elim; ind2; ind2=ind2->next) {
    if (!strcmpl(ind->id, ind2->id))
      break;
  }

  while (listlen(ind2->marker)) {
      rannum = ignuin(1,listlen(ind2->marker));
      mkr = markernumber(ind2, rannum);

      // Check to see if it is different from the original allele
      if ( (min(mkr->allele1, mkr->allele2) != min(allele1, allele2)) ||
	   (max(mkr->allele1, mkr->allele2) != max(allele1, allele2))) {

	origmarker->allele1 = mkr->allele1;
	origmarker->allele2 = mkr->allele2;

	// Assume we changed
	*changed = 1;

	// Check to see that the new genotype is legal (i.e. still Mendelian)
	ind4 = (individual *) GenotypeElimination(ind->pedigree, markernum);
	for (ind3 = ind4; ind3; ind3 = ind3->next) {
	  if (listlen(ind3->marker) == 0) {
	    // NOT MENDELIAN GENOTYPE
	    // Set the original genotype to missing
	    origmarker->allele1 = 0;
	    origmarker->allele2 = 0;

	    removelist(&ind2->marker, mkr);
	    *changed = 0;
	    break;
	  }
	}
	FreeIndividualList(ind4);
      }
      else {
	removelist(&ind2->marker, mkr);
      }
      if (*changed == 1)
	break;
    }
  
  // Check to see if any new marker was found
  if (listlen(ind2->marker)==0) {
    origmarker->allele1 = allele1;
    origmarker->allele2 = allele2;
    *changed = 0;
  }

  FreeIndividualList(elim);
}

//
// Make a random genotyping error
//
void RandomGenotypeError(individual *indlist, int markernum)
{
  individual *ind;
  int nObs, *pnPerm, i, nRandom, temp, nChange, nTotalChange, nError;
  float freq;

  printf("Input error frequency: ");
  InputLine(buf, BUFFERSIZE);
  freq = atof(igetstr(buf));
  if (freq<= 0 || freq>=1.0) {
      printf("Error frequency must be between 0 and 1\n");
      return;
    }
  printf("Simulates a random genotyping error with frequency %5.3f and Mendelian consistent\n", freq);

  nTotalChange = 0;
  nError = 0;
  nObs = listlen(indlist);
  pnPerm = (int *) cmalloc(nObs*sizeof(int));

  for (i=0; i<nObs; i++)
    pnPerm[i] = i+1;

  // Permutate the way people are selected
  for (i=0; i<nObs; i++) {
    nRandom = ignuin(0, nObs-1);
    temp = pnPerm[i];
    pnPerm[i] = pnPerm[nRandom];
    pnPerm[nRandom] = temp;    
  }  

  // Create errors
  for (i=0; i<nObs; i++) {
    // if an error occurs
    if (genunf(0,1)<= freq) {
      nError++;
      for (ind = indlist, temp=1; ind; ind = ind->next, temp++) {
	if (temp==pnPerm[i])
	  break;
      }
      GenotypeErrorMendelian(ind, markernum, &nChange);
      nTotalChange += nChange;      
    }
  }

  printf("%4.2f%% (%d individuals) of the individuals were sampled\n", 100.0*(double) nError/ (double)nObs, nError);
  printf("%4.2f%% (%d individuals) of samples could not get a different genotype\n", nError>0 ? (double) 100.0*(nError-nTotalChange)/(double) nError : 0, (nError-nTotalChange));
  printf("Real error-rate: %4.2f%%\n", (double)100.0*nTotalChange/(double) nObs);

  free(pnPerm);
}

// ---------------------------------------------------------
// Creates a completely random genotype for individual ind
// Genotype Error
//
// ---------------------------------------------------------
void GenotypeError(individual *ind, int markernum, int errorfct) {
  markerlist *marker;
  int allele1, allele2;
  int newallele;

  marker = markernumber(ind, markernum);
  allele1 = marker->allele1;
  allele2 = marker->allele2;

  switch (errorfct) {
  case 0 : // Jumps up and down. Equal frequencies
    marker->allele1 = 1;
  default: 
    newallele = ignuin(1,4);
    if (marker->allele1<=newallele) {
      newallele++;
    }
    marker->allele1 = newallele;
    
  }

}

//
//
// Randomly 
//
//

void RemovePhenotypeData(individual *ind, int traitnum)
{
  quanttrait *trait;

  trait = traitvalue(ind,traitnum);
  trait->sysmiss = 1;
}

void RemoveGenotypeData(individual *ind, int markernum)
{
  markerlist *marker;

  marker = markernumber(ind, markernum);

  marker->allele1 = 0;
  marker->allele2 = 0;
}


void SimulateRandomMissingGenotypes(individual *indlist, int markernum, double freq)
{
  individual *ind;

  for (ind = indlist; ind; ind = ind->next)
    {
      if (genunf(0,1)<=freq)
	{
	  RemoveGenotypeData(ind, markernum);
	}
    }

  printf("Simulated missing genotypes in %4.1f%% of typed markers for marker %d\n", freq*100, markernum);
}



/*

 NEW STUFF


 */

individual *PrepareNewIndividual(char* indid, char* pedid, int sex, char *father, char* mother) 
{
  individual *ind;
  char minibuf[25];

  ind = (individual *) cmalloc (sizeof (individual));
  memset (ind,0,sizeof (individual));
  strcpy(ind->pedigree, pedid);  
  strcpy(minibuf, pedid);
  strcat(minibuf, "-");
  strcat(minibuf, indid);
  sprintf(ind->id, minibuf);
  if (sex == -1) {
    if (ignbin(1,.5))
      ind->sex = S_MALE;
    else
      ind->sex = S_FEMALE;
  }
  else
    ind->sex = sex;

  strcpy(ind->tmpstr1, pedid);
  strcat(ind->tmpstr1, "-");
  strcat(ind->tmpstr1, father);

  strcpy(ind->tmpstr2, pedid);
  strcat(ind->tmpstr2, "-");
  strcat(ind->tmpstr2, mother);

  return (ind);
}


/*****************
 *
 * SimulatePedigreeStructure simulates a specific pedigree structure. 0 = Nuclear families, 1 = CEPH
 *
 *
 *****************/
individual *SimulatePedigreeStructure(int nPedigrees, int nPedType, int nTraits, int nQTL)
{
  int nFam, mkrnum, count, nMarkers, i;
  individual *ind, *ind2, *finaldata;
  markerlist *marker;
  quanttrait *qt;
  namelist *nl;
  char minibuf[24];


  // Check that we have available marker/allele info
  if (!allfreq)
    {
      printf("WARNING: No marker allele frequency information avaliable to simulate data.\n");
      printf("Remember to load a parameter file first\n");
      return 0;
    }

  // Check that known pedigree structure
  if (nPedType<0 || nPedType>2) {
    printf("WARNING: Trying to simulate unknown pedigree structure\n");
    return 0;
  }

  finaldata = 0;

  for (nFam=1; nFam<=nPedigrees; nFam++) {
    switch (nPedType) {
    case 0: // Nuclear families
      sprintf(minibuf, "F%d", nFam);
      ind = PrepareNewIndividual("1", minibuf, S_MALE, "", "");
      addlist (&finaldata,ind);
      ind = PrepareNewIndividual("2", minibuf, S_FEMALE, "", "");
      addlist (&finaldata,ind);
      ind = PrepareNewIndividual("3", minibuf, -1, "1", "2");
      addlist (&finaldata,ind);
      ind = PrepareNewIndividual("4", minibuf, -1, "1", "2");
      addlist (&finaldata,ind);      
      ind = PrepareNewIndividual("5", minibuf, -1, "1", "2");
      addlist (&finaldata,ind);
      ind = PrepareNewIndividual("6", minibuf, -1, "1", "2");
      addlist (&finaldata,ind);      
      break;
    case 1: // CEPH families
      sprintf(minibuf, "F%d", nFam);
      ind = PrepareNewIndividual("1", minibuf, S_MALE, "", "");
      addlist (&finaldata,ind);
      ind = PrepareNewIndividual("2", minibuf, S_FEMALE, "", "");
      addlist (&finaldata,ind);
      ind = PrepareNewIndividual("3", minibuf, S_MALE, "", "");
      addlist (&finaldata,ind);
      ind = PrepareNewIndividual("4", minibuf, S_FEMALE, "", "");
      addlist (&finaldata,ind);
      ind = PrepareNewIndividual("5", minibuf, S_MALE, "1", "2");
      addlist (&finaldata,ind);
      ind = PrepareNewIndividual("6", minibuf, S_FEMALE, "3", "4");
      addlist (&finaldata,ind);      
      break;
    }
  }

  count = 0;
  // Fix the individual according to their
  // Should now fix fathers and mothers
  for (ind = finaldata; ind; ind = ind->next) {
    ind->father = FindIndListMember(finaldata, ind->tmpstr1);
    ind->mother = FindIndListMember(finaldata, ind->tmpstr2);

    // if a person has a father and a mother, then add the
    //  person to the parents lists
    if (ind->father && ind->mother)
      {
	adduniqueidlist (&ind->father->offspring, ind);
	adduniqueidlist (&ind->father->mate, ind->mother);
	adduniqueidlist (&ind->mother->offspring,ind);
	adduniqueidlist (&ind->mother->mate,ind->father);
      }
    else if (ind->father || ind->mother) { // The person has only one parent
      count = 1;
      sprintf(buf2, "%s has only one parent in pedigree file\n", ind->id);
      WriteErrorMsg(buf2);
    }    
  }

  if (count)
    printf("ERROR: A least one parent not found in pedigree file\n");

  // Fixing sibs
  for (ind = finaldata; ind; ind = ind->next) {
    for (ind2 = ind->next; ind2; ind2 = ind2->next) {
      // If same father or mother then add to list of both
      if (((ind->father == ind2->father) && (ind->father)) || 
          ((ind->mother == ind2->mother) && (ind->mother))) {
        adduniqueidlist(&ind->sib, ind2);
        adduniqueidlist(&ind2->sib, ind);
      }
    }
  }


  // Add markers and phenotypes

  for (i = 1; i<=nTraits ; i++) {
    nl = (namelist *) cmalloc (sizeof (namelist));
    memset (nl,0,sizeof (namelist));

    sprintf(nl->name, "Trait %d",i);
    addlist(&traitnames, nl);
  }

  for (ind = finaldata; ind; ind = ind->next) {

    nMarkers = listlen(allfreq);
    for (mkrnum = 1; mkrnum <= nMarkers; mkrnum++) {
      // Gets the relevant allele frequency for the marker

      marker = (markerlist *) cmalloc (sizeof (markerlist));
      memset (marker,0,sizeof (markerlist));

      marker->allele1 = 1;
      marker->allele2 = 1;
      
      addlist(&ind->marker, marker);
    }

    // This qt represents the observed trait
    for (i = 0; i< nTraits; i++) {
      qt = (quanttrait *) cmalloc (sizeof (quanttrait));
      memset (qt, 0, sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 1);
      addlist(&ind->qttrait,qt);
    }
  }

  return finaldata;
}




// Assumes order is ok. Should be fixed
// families each with nooff offspring
// A single QTL a position location with an effect of effect

//
// The QTL must be inside the chromosome or else it won't work
/*******************************************************
 *
 *   Tager ikke højde for Marker order
 *
 *  Assumes QTLs are independent.
 *
 ********************************************************/

individual *SimulateMultiQTLChromosome(int families, int nooff, MATRIX QTLloc, MATRIX QTLadd, MATRIX QTLdom, MATRIX QTLprev, double resadd, double resdom)
{
  int i, j, k, nogenes, count, mkrnum, NMarkers;
  individual *ind, *ind2, *father, *mother, *finaldata, *minilist;
  markerlist *marker, *fmkr, *mmkr, *mkr;
  quanttrait *qt;
  namelist *nl;
  namelist *ped, *Pedigrees;
  double randomnumber, res, recomb;
  markerlist *fixdata[MAXMARKERS+1];
  float dist, mindist, dgenefreq, heritab;
  int diseasecopies, persons, truemarker, trueqtl, nTraits = 1;
  int nQTL = QTLloc.Rows();
  int nPseudo = numberofmarkers()+QTLloc.Rows();

  if (!allfreq) {
      printf("WARNING: No marker allele frequency information avaliable to simulate data.\n         Remember to load a parameter file first\n");
      return 0;
    }

  persons = 0;
  finaldata = 0;

  //
  // Start error checking
  //
  if ((QTLloc.Rows() != QTLadd.Rows()) && (QTLloc.Rows() != QTLdom.Rows()) && (QTLloc.Rows() != QTLprev.Rows())) {
    printf("WARNING: Inconsistency in QTL information\n");
    return 0;
  }


  // Start by sorting the QTL positions so they are in ascending order
  MATRIX InfoPos(nPseudo,2);
  InfoPos = 0;
  k=1;
  for (i=1; i<numberofmarkers(); i++)  {      
    // Include marker
    InfoPos(k,1) = MarkerDistance(0, i);
    k++;
    for (j=1; j<= QTLloc.Rows(); j++) {
      if (QTLloc(j,1)>=MarkerDistance(0, i) && QTLloc(j,1)<MarkerDistance(0, i+1)) {
	InfoPos(k,1) = QTLloc(j,1);
	InfoPos(k,2) = 1;
	k++;
      }      
    }  
  }      
  // Include last marker
  InfoPos(k,1) = MarkerDistance(0, numberofmarkers());

  //  InfoPos.Print();

  for (i=1; i<=families; i++) {
    // Generates the founders

    // First the father
    father = (individual *) cmalloc (sizeof (individual));
    memset (father,0,sizeof (individual));
    
    sprintf(father->pedigree, "F%d",i);
    sprintf(father->id, "%d-1",i);
    father->sex = S_MALE;
    
    for (mkrnum = 1, truemarker = 0, trueqtl=0; mkrnum <= nPseudo; mkrnum++) {
      marker = (markerlist *) cmalloc (sizeof (markerlist));
      memset (marker,0,sizeof (markerlist));

      // This is a regular marker
      if (InfoPos(mkrnum,2)==0) {
	truemarker++;
	CopyAlleleFreq(truemarker);
	
	count = 1;
	res = allelefreq[count];
	// Gives the father allele 1
	randomnumber = genunf(0,1);
	while(randomnumber>res) {
	  count++;
	  res += allelefreq[count];
	}
	marker->allele1 = count;
	
	// and allele 2
	count = 1;
	res = allelefreq[count];
	randomnumber = genunf(0,1);
	while(randomnumber>res)  {
	  count++;
	  res += allelefreq[count];
	}
	marker->allele2 = count;

#ifdef fullyinformative
	marker->allele1 = 1;
	marker->allele2 = 2;      
#endif

      }
      else {
	// This a QTL
	trueqtl++;
	marker->allele1 = 0;
	marker->allele2 = 0;

	if (genunf(0,1)<QTLprev(trueqtl,1))
	  marker->allele1 = 1;
	if (genunf(0,1)<QTLprev(trueqtl,1))
	  marker->allele2 = 1;
      }

      addlist(&father->marker, marker);
    }

    // This qt represents the observed trait
    qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
    memset (qt,0,sizeof (quanttrait));
    
    qt->sysmiss = 0;
    qt->value = 0; 
    addlist(&father->qttrait,qt);

    // Makes an additional irrelevant covariate
    qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
    memset (qt,0,sizeof (quanttrait));
    qt->sysmiss = 0;
    qt->value = gennor(0, 2);
    addlist(&father->qttrait,qt);

    addlist (&finaldata,father);
    father->globalid = ++persons;
    father->localid = persons;

    // Then the mother
    mother = (individual *) cmalloc (sizeof (individual));
    memset (mother,0,sizeof (individual));
    
    sprintf(mother->pedigree, "F%d",i);
    sprintf(mother->id, "%d-2",i);
    mother->sex = S_FEMALE;

    for (mkrnum = 1, truemarker = 0, trueqtl=0; mkrnum <= nPseudo; mkrnum++) {
      marker = (markerlist *) cmalloc (sizeof (markerlist));
      memset (marker,0,sizeof (markerlist));

      // This is a regular marker
      if (InfoPos(mkrnum,2)==0) {
	truemarker++;
	CopyAlleleFreq(truemarker);
	
	count = 1;
	res = allelefreq[count];
	// Gives the father allele 1
	randomnumber = genunf(0,1);
	while(randomnumber>res) {
	  count++;
	  res += allelefreq[count];
	}
	marker->allele1 = count;
	
	// and allele 2
	count = 1;
	res = allelefreq[count];
	randomnumber = genunf(0,1);
	while(randomnumber>res)  {
	  count++;
	  res += allelefreq[count];
	}
	marker->allele2 = count;

#ifdef fullyinformative
	marker->allele1 = 1;
	marker->allele2 = 2;      
#endif

      }
      else {
	// This a QTL
	trueqtl++;
	marker->allele1 = 0;
	marker->allele2 = 0;

	if (genunf(0,1)<QTLprev(trueqtl,1))
	  marker->allele1 = 1;
	if (genunf(0,1)<QTLprev(trueqtl,1))
	  marker->allele2 = 1;
      }

      addlist(&mother->marker, marker);
    }

    // This qt represents the observed trait
    qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
    memset (qt,0,sizeof (quanttrait));
    
    qt->sysmiss = 0;
    qt->value = 0;
    addlist(&mother->qttrait,qt);

    // Makes an additional irrelevant covariate
    qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
    memset (qt,0,sizeof (quanttrait));
    qt->sysmiss = 0;
    qt->value = gennor(0, 2);
    addlist(&mother->qttrait,qt);

    addlist (&finaldata,mother);
    mother->globalid = ++persons;
    mother->localid = persons;

    // Should now generate the children. 
    
    // Generates the offspring
    for (j=0; j<nooff; j++) {
      ind = (individual *) cmalloc (sizeof (individual));
      memset (ind,0,sizeof (individual));

      sprintf(ind->pedigree, "F%d",i);
      sprintf(ind->id, "%d-%d",i,j+3);

      if (genunf(0,1)<.5)
	ind->sex = S_MALE;
      else
	ind->sex = S_FEMALE;

      ind->father = father;
      ind->mother = mother;

      adduniqueidlist (&ind->father->offspring, ind);
      adduniqueidlist (&ind->father->mate, ind->mother);
      adduniqueidlist (&ind->mother->offspring, ind);
      adduniqueidlist (&ind->mother->mate,ind->father); 

      // tmpstr 1 and 2 hold information about which haplotype
      // the person received last time. Starts random
      sprintf(ind->tmpstr1, "%d", (int) ignbin(1,.5));
      sprintf(ind->tmpstr2, "%d", (int) ignbin(1,.5));
      
      for (mkrnum = 1; mkrnum <= nPseudo; mkrnum++) {

	marker = (markerlist *) cmalloc (sizeof (markerlist));
	memset (marker,0,sizeof (markerlist));
	
	fmkr = markernumber(father, mkrnum);
	mmkr = markernumber(mother, mkrnum);

	// Check to see if a recombination has occured for all but
	// the first marker
	if (mkrnum>1) {
     	  recomb = InverseMapFunction((InfoPos(mkrnum,1)-InfoPos(mkrnum-1,1))/100);

	  //	  printf("Recomb frac between mkr %d and %d : %f\n", mkrnum-1, mkrnum, recomb);

	  // Checks for recombination for father
	  if (genunf(0,1)<recomb)
	    sprintf(ind->tmpstr1, "%d", 1 - atoi(ind->tmpstr1));
	  // and mother
	  if (genunf(0,1)<recomb)
	    sprintf(ind->tmpstr2, "%d", 1 - atoi(ind->tmpstr2));
	}
	// Give the person the allele
	if (!strcmpl(ind->tmpstr1,"0"))
	  marker->allele1 = fmkr->allele1;
	else
	  marker->allele1 = fmkr->allele2;
	
	if (!strcmpl(ind->tmpstr2,"0"))
	  marker->allele2 = mmkr->allele1;
	else
	  marker->allele2 = mmkr->allele2;
	
	addlist(&ind->marker, marker);	
      }

      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = 0;
      addlist(&ind->qttrait,qt);
      
      // Makes an additional irrelevant covariate
      qt = (quanttrait *) cmalloc  (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      qt->sysmiss = 0;
      qt->value = gennor(0, 2);
      addlist(&ind->qttrait,qt);
      
      addlist (&finaldata,ind);
      ind->globalid = ++persons;
      ind->localid = persons;
      
    }
    
  }

  // Fixing sibs
  for (ind = finaldata; ind; ind = ind->next)
    {
      for (ind2 = ind->next; ind2; ind2 = ind2->next)
	{
	  // If same father or mother then add to list of both
	  if (((ind->father == ind2->father) && (ind->father)) || 
	      ((ind->mother == ind2->mother) && (ind->mother)))
	    {
	      adduniqueidlist(&ind->sib, ind2);
	      adduniqueidlist(&ind2->sib, ind);
	    }
	}
    }

  // Needs to reorder the markers back for non-founders
  // so everything else works fine
  /*
  for (ind = finaldata; ind; ind = ind->next)
    {
      if (!founder(ind))
	{
	  i = 2;
	  memset(fixdata, 0, sizeof(fixdata));
	  fixdata[1] = ind->marker;
	  for (mkr = ind->marker; mkr; mkr = mkr->next)
	    {
	      fixdata[i] = mkr->next;
	      i++;
	    }

	  ind->marker = fixdata[invorder[1]];
	  for (i = 1; i<= numberofmarkers(); i++)
	    {
	      mkr = fixdata[i];
	      if (order[i] < numberofmarkers())
		mkr->next = fixdata[invorder[order[i]+1]];
	      else
		mkr->next = 0;
	    }
	}
    }

  */

  // Finally adds variation by measurement error, resadd and resdom
  MATRIX Res;
  Pedigrees = MakePedigreeList(finaldata);
  for (ped = Pedigrees; ped; ped = ped->next) {
    minilist = SelectPedigree(ped->name, finaldata);

    Res.Resize(listlen(minilist), 1);

    MATRIX Error = MakeIdentMatrix(listlen(minilist));
    if (resadd>0) {
      MATRIX ResA = MMakeKinshipMatrix(minilist);
      MATRIX ResD = MMakeDelta7Matrix(minilist, ResA*2);
      Error = Error + ResA*(resadd*2)+ResD*resdom; 
    }
    Error = Cholesky(Error);

    for (j = 1; j<=nTraits; j++) {
      // Generates standard normals
      for (i=1; i<=Res.Rows(); i++) {
	Res(i,1) = gennor(0,1);
      }

      Res = Error*Res;

      for (ind = minilist, i=1; ind; ind = ind->next, i++) {
	qt = traitvalue(FindIndListMember(finaldata,ind->id), j);
	qt->value += Res(i,1);
	trueqtl = 0;
	for (k = 1; k<=nPseudo; k++) {
	  if (InfoPos(k,2) == 1) {
	    trueqtl++;
	    mkr = markernumber(ind, k);
	    qt->value += QTLadd(trueqtl,1)*(mkr->allele1+mkr->allele2-1);
	    // Add dominance effects
	    if (mkr->allele1+mkr->allele2 ==1) {
	      qt->value += QTLdom(trueqtl,1)*(mkr->allele1+mkr->allele2-1);
	    }
	  }
	}
      }
    }

    FreeIndividualList(minilist);
  }
  FreeNameList(Pedigrees);


  // Then add trait names
  // Not necessary to add marker names as they've been read from the
  // parameter file already
  for (i=1; i<=2; i++)
    {
      nl = (namelist *) cmalloc (sizeof (namelist));
      memset (nl,0,sizeof (namelist));

      sprintf(nl->name, "Trait %d",i);
      addlist(&traitnames, nl);
    }


  //  tmp = addeffect + domeffect*(1-2*prev);
  //  heritab = 2*prev*(1-prev)*SQR(tmp);
  //  tmp = 2*prev*(1-prev)*domeffect;
  //  heritab += SQR(tmp);

  // Removes markers corresponding to QTLs
  for (i = nPseudo; i>0; i--) {
    if (InfoPos(i,2)==1) {
      RemoveMarker(finaldata, i);
    }
  }

#ifdef fullyinformative
  printf("Simulated %d fully informative families with %d offspring each\n", families, nooff);
#else
  printf("Simulated %d families with %d offspring each\n", families, nooff);
#endif
  //  printf("The heritability of the simulated QTL is: %5.3f (Add. effect: %5.3f, dom. dev. effect: %5.3f) positioned at %5.3f\n", heritab / (1 + heritab), addeffect, domeffect, location);

  // Fixes the order etc.

  return finaldata; 
}




/*******************************************************************
 *
 * CloneDataset clones the current existing dataset in nclone copies
 *
 *******************************************************************/

individual *CloneDataset(individual *indlist, int nclones)
{
  individual *finaldata, *tmp, *ind;
  char name[NAMESIZE];
  int clone;

  finaldata = NULL;
  for (clone=1; clone <= nclones; clone++) {
    tmp = CopyIndividualList(indlist);
    // Start by fixing the pedigree id's
    // and individual id's
    for (ind = tmp; ind; ind=ind->next) {      
      strncpy(name, ind->pedigree, NAMESIZE-1);
      sprintf(ind->pedigree, "C%d-%s", clone, name);
      strncpy(name, ind->id, NAMESIZE-1);
      sprintf(ind->id, "C%d-%s", clone, name);
    }

    // Now add the cloned data to the current data set
    // Find last element
    if (clone==1) {
      finaldata=tmp;
    }
    else {
      for (ind = finaldata; ind->next; ind=ind->next) {    
      }
      // Add link to new list
      ind->next=tmp;  
    }
  }
  printf("Made %d clones of current dataset\n", nclones);
  return finaldata;
}



#undef seed
#undef xlinked
