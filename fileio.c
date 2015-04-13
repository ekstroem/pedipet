/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 *  Various functions for importing and exporting data in a variety of formats
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
#include "lists.h"
#include "fileio.h"

//#define debug

#ifndef FILEIO
#define FILEIO


// I am not sure the XGH think works and can't remember why I need it
// #define XGH


//
//  Export data to other formats
//


/*
     SOLAR
*/

int ExportSolarPedigree (individual *indlist)
{
  individual *ind;
  

  // Checks that data exists
  if (!listlen(indlist))
  {
    // XXX Should create a proper output rutine for warnings
    //    WriteWarning("No pedigree data avaliable for creating SOLAR file\n");
    return -1;
  }

  // Makes pedigree file
  cfopen ("solar.ped", "w");
  fprintf(F, "id,fa,mo,sex\n");  // The first line
  for (ind = indlist; ind; ind = ind->next)
  {
    if (founder(ind))
      fprintf(F, "%s,0,0,%d\n",ind->id, ind->sex);
    else
      fprintf(F, "%s,%s,%s,%d\n",ind->id, ind->father->id, ind->mother->id,ind->sex);
  }
  fclose(F);
  return 0;
}

int ExportSolarMarker (individual *indlist)
{
  individual *ind;
  markerlist *marker;
  char mkrname[NAMESIZE], tempstr[NAMESIZE];
  int i;
  

  // Makes the marker file
  cfopen ("solar.mkr", "w");
  fprintf(F, "id");
  for (i = 1; i<=numberofmarkers(); i++)
  {
    // Remove white spaces in marker name
    RemoveWhiteSpace(tempstr, GetName(markernames, order[i]));
    strncpy(mkrname, tempstr, NAMESIZE);
    fprintf(F, ",%s", mkrname);
  }
  fprintf(F, "\n");

  for (ind = indlist; ind; ind = ind->next)
  {
    fprintf(F, "%s", ind->id);
    for (i = 1; i<=numberofmarkers(); i++)
    {
      marker = markernumber(ind, order[i]);
      fprintf(F, ",%d/%d", marker->allele1, marker->allele2);
    }
    fprintf(F, "\n");
  }

  fclose(F);
  puts("Made Solar marker file");
  return 0;
}


int ExportSolarMap (void)
{
  char mkrname[NAMESIZE], tempstr[NAMESIZE];
  int i;
  double currentpos;
  

  printf("Which chromosome is this: ");
  InputLine(buf, BUFFERSIZE);

  currentpos = 0;
  // Makes the marker file
  cfopen ("solar.map", "w");
  fprintf(F, "%s\n",buf);
  for (i = 1; i<=numberofmarkers(); i++)
  {
    // Remove white spaces in marker name
    RemoveWhiteSpace(tempstr, GetName(markernames, order[i]));
    strncpy(mkrname, tempstr, NAMESIZE);

    fprintf(F, "%-18s  %f\n", mkrname, MarkerDistance(0,i));
  }
  fclose(F);
  puts("Made Solar map file");
  return 0;
}


int ExportSolarFrequency (void)
{
  char mkrname[NAMESIZE], tempstr[NAMESIZE];
  int i, j;
  

  // Makes the frequency file
  cfopen ("solar.frq", "w");

  for (i = 1; i<=numberofmarkers(); i++)
  {
    // Remove white spaces in marker name
    RemoveWhiteSpace(tempstr, GetName(markernames, order[i]));
    strncpy(mkrname, tempstr, NAMESIZE);

    CopyAlleleFreq(order[i]);

    fprintf(F, "%-18s  ", mkrname);
    for (j=1; j<= NumberOfAlleles(order[i]); j++)
      fprintf(F, "%2d %5.3f  ", j, allelefreq[j]);

    fprintf(F, "\n");
  }

  fclose(F);

  puts("Made Solar frequency file");
  return 0;
}

int ExportSolarPhenotype (individual *indlist)
{
  individual *ind;
  char mkrname[NAMESIZE], tempstr[NAMESIZE];
  int i;
  double currentpos, value;
  

  currentpos = 0;
  // Makes the marker file
  cfopen ("solar.phe", "w");

  // Prints the trait names
  fprintf(F,"id");
  for (i = 1; i<=numberoftraits(); i++)
  {
    // Remove white spaces in marker name
    RemoveWhiteSpace(tempstr, GetName(traitnames, i));
    strncpy(mkrname, tempstr, 18);

    fprintf(F, ",%s", mkrname);
  }
  fprintf(F,"\n");

  for (ind = indlist; ind; ind = ind->next)
  {
    fprintf(F, "%s", ind->id);
    for (i = 1; i<=numberoftraits(); i++)
    {
      fprintf(F, ",");
      if (!traitmiss(ind, i)) {
	value = trait(ind, i);
	// Print only integer part if it is an integer
	if ( (int) value == value) {
	  fprintf(F, "%d", (int) value);
	}
	else {
	  fprintf(F, "%f", value);
	}
      }
    }
    fprintf(F,"\n");
  }

  fclose(F);
  puts("Made Solar phenotype file");
  return 0;
}


/*

    MENDEL and friends.
    This is created for RELPAIR

 */

int ExportRelpairControl()
{
  

  cfopen ("relpair.ctl", "w");
  fprintf(F, "relpair.loc\n"); 
  fprintf(F, "relpair.ped\n"); 
  fprintf(F, "relpair.out\n");
  fprintf(F, "all\n"); 
  fprintf(F, "y\n");                 // Map information echoed
  fprintf(F, "n\n"); 
  fprintf(F, "F\n"); 
  fprintf(F, "M\n"); 
  fprintf(F, "1\n");  
  fclose(F);

  printf("Created relpair control file\n");

  return 1;
}


int ExportRelpairLocus(int chromosome)
{
  int i, j;
  char mkrname[NAMESIZE], tempstr[NAMESIZE];
  

  cfopen ("relpair.loc", "w");

  for (i = 1; i<=numberofmarkers(); i++)
  {
    // Remove white spaces in marker name
    RemoveWhiteSpace(tempstr, GetName(markernames, order[i]));
    strncpy(mkrname, tempstr, NAMESIZE);
    // Prints the locus information
    fprintf(F, "%-8s%8s%2d%2d%4d%8.5f\n", mkrname, "AUTOSOME", 
	    NumberOfAlleles(order[i]), 0, chromosome, MarkerDistance(0,i)/100);
    // Should now print the frequency information

    CopyAlleleFreq(order[i]);

    for (j=1; j<= NumberOfAlleles(order[i]); j++)
      fprintf(F, "%c       %8.5f\n", 64+j, allelefreq[j]);
  }
    
  fclose(F);

  printf("Created relpair locus file\n");

  return 1;
}


int ExportRelpairPedigree(individual *indlist)
{
  int i, famsize;
  namelist *pedlist, *pedigr;
  individual *ind, *ind2;
  markerlist *marker;
  

  cfopen ("relpair.ped", "w");

  fprintf(F, "(I4,1X,A8)\n");
  fprintf(F, "(3A8,2A1");

  if (numberofmarkers())
      fprintf(F, ",A3");

  if (numberofmarkers()>1) {
      fprintf(F, ",%d(1X,A3))\n", numberofmarkers()-1);    
  }

  pedigr = MakePedigreeList(indlist);

  for (pedlist = pedigr; pedlist; pedlist = pedlist ->next) {
    // Start by calculating the family size
    famsize = 0;
    for (ind2 = indlist; ind2; ind2 = ind2->next) {
      if (!strcmpl(pedlist->name, ind2->pedigree)) {
	famsize++;
      }
    }
    fprintf(F, "%4d %s \n",famsize, pedlist->name);

    for (ind = indlist; ind; ind = ind->next) {
      // Pick correct pedigree
      if (!strcmpl(ind->pedigree,pedlist->name)) {

	fprintf(F, "%-8s", ind->id);
	if (ind->father)
	  fprintf(F, "%-8s", ind->father->id);
	else
	  fprintf(F, "%-8s", "0");
	if (ind->mother)
	  fprintf(F, "%-8s", ind->mother->id);
	else
	  fprintf(F, "%-8s", "0");

	if (ind->sex == S_MALE)
	  fprintf(F, "M");
	else
	  fprintf(F, "F");

	fprintf(F, " ");  // MZ status

	// Prints the genotypes
	for (i = 1; i<=numberofmarkers(); i++) {
	  marker = markernumber(ind, order[i]);
	  if (marker->allele1*marker->allele2)
	    fprintf(F, "%c/%c ", marker->allele1 + 64, marker->allele2 + 64);
	  else
	    fprintf(F, "    ");
	}
	fprintf(F, "\n");
      }
    }
  }
    
  fclose(F);


  freelist(pedigr);

  printf("Created relpair pedigree file\n");

  return 1;
}




/*
   CRIMAP
*/

int ExportCRIMAPGen(individual *indlist, int xlinked)
{
  individual *ind;
  
  namelist *pedigree, *ped;
  markerlist *mkr;
  char mkrname[NAMESIZE], tempstr[NAMESIZE];
  int nFamilies, nLoci, i, nInternalID, nMembers, al1, al2, pseudoal;

  nLoci = listlen(indlist->marker);
  pedigree = MakePedigreeList(indlist);
  nFamilies = listlen(pedigree);

  cfopen ("crimap.gen", "w");
  fprintf(F, "%d %d", nFamilies, nLoci);  // The first line

  // Prints the markers in the correct order
  for (i = 1; i<= nLoci; i++)
  {
    // Remove white spaces in marker name
    RemoveWhiteSpace(tempstr, GetName(markernames, order[i]));
    strncpy(mkrname, tempstr, 15);
    fprintf(F, " %s", mkrname);
  }
  fprintf(F, "\n");

  // Fixes internal id's for use below
  for (ped = pedigree; ped; ped = ped->next)
  {
    nInternalID = 1;
    for (ind = indlist; ind; ind = ind->next)
    {
      // From the right pedigree
      if (!strcmpl(ind->pedigree, ped->name))
      {
        ind->localid = nInternalID;
        nInternalID++;
      }
    }    
  }  

  // Writes each family 
  for (ped = pedigree; ped; ped = ped->next)
  {
    // Family name
    fprintf(F, "%s\n", ped->name); 

    nMembers = 0;
    // Count the number of family members
    for (ind = indlist; ind; ind = ind->next)
    {
      // From the right pedigree
      if (!strcmpl(ind->pedigree, ped->name))
	nMembers++;
    }
    fprintf(F, "%d\n", nMembers);

    for (ind = indlist; ind; ind = ind->next)
    {
      // From the right pedigree
      if (!strcmpl(ind->pedigree, ped->name))
      {
        if (founder(ind))
          fprintf(F, "%d 0 0 ", ind->localid);
        else
          fprintf(F, "%d %d %d ", ind->localid, ind->father->localid, ind->mother->localid);

        switch (ind->sex)
	{
	  case S_MALE:   fprintf(F, "1"); break;
	  case S_FEMALE: fprintf(F, "0"); break;
 	  default:       fprintf(F, "3"); break;
	}

	// Should now print the alleles for each locus

	// If xlinked, then create a pseudo allele for each marker
	// to make the males heterozygous
        for (i = 1; i<= nLoci; i++)
	{
          mkr = markernumber(ind, order[i]);
	  al1 = mkr->allele1;
	  al2 = mkr->allele2;
	  
	  if (xlinked && ind->sex == S_MALE && al1>0 && al2>0) {
	    pseudoal = NumberOfAlleles(order[i])+1;
	    al2 = pseudoal;
	  }
          fprintf(F, " %d %d", al1, al2);
	}
        fprintf(F, "\n");
      }
    }
  }


  fclose(F);

  printf("Created crimap gen file\n");
  if (xlinked) {
    printf("Pseudo-allele generated for all males\n");
  }

  return 0;  
}

int ExportCRIMAPPar(individual *indlist)
{
  int nLoci, i;
  

  nLoci = listlen(indlist->marker);

  cfopen ("crimap.par", "w");
  fprintf(F, "dat_file crimap.dat *\n");
  fprintf(F, "gen_file crimap.gen *\n");
  fprintf(F, "ord_file crimap.ord *\n");
  fprintf(F, "nb_our_alloc 3000000 *\n");
  fprintf(F, "SEX_EQ  1 *\n");
  fprintf(F, "TOL .010000 *\n");
  fprintf(F, "PUK_NUM_ORDERS_TOL  6 *\n");
  fprintf(F, "PK_NUM_ORDERS_TOL  8 *\n");
  fprintf(F, "PUK_LIKE_TOL 3.000 *\n");
  fprintf(F, "PK_LIKE_TOL 3.000 *\n");
  fprintf(F, "use_ord_file  0 *\n");
  fprintf(F, "write_ord_file  1 *\n");
  fprintf(F, "use_haps  1 *\n");
  fprintf(F, "ordered_loci 0 1  *\n");

  fprintf(F, "inserted_loci ");
  for (i=2; i<nLoci; i++)
    fprintf(F, "%d ",i);
  fprintf(F, "  *\n");

  fprintf(F, "END\n");
  fclose(F);

  printf("Created crimap parameter file\n");
  return 0;

}



/*
    SimWalk2 / MENDEL
*/


void ExportMendelPedigree (individual *indlist)
{
  int i;
  namelist *pedigree, *pedlist;
  markerlist *marker;

  individual *ind;
  pedlist = MakePedigreeList(indlist);

  cfopen ("mendel.ped", "w");
  fprintf(F, "(I4,A15,A8)\n");
  fprintf(F, "(3A11,3A2");

  // Add marker info
  if (numberofmarkers()>0)
    fprintf(F, ",%d(1X,A7)", numberofmarkers());
  if (numberoftraits()>0)
    fprintf(F, ",%d(1X,A8)", numberoftraits());

  fprintf(F, ")\n");

  for (pedigree = pedlist; pedigree; pedigree = pedigree ->next) {
    fprintf(F, "%-4d%-15s\n",FamilySize(indlist, pedigree->name), pedigree->name);

    for (ind = indlist; ind; ind = ind->next)
    {
      if (!strcmpl(pedigree->name, ind->pedigree)) {
	if (founder(ind)) {
	  fprintf(F, "%-11s%-11s%-11s", ind->id, "","");
	}
	else {
 	  fprintf(F, "%-11s%-11s%-11s", ind->id, ind->father->id, ind->mother->id);
	}
	// Spaces here for non-existing disease status and monozygotic twin status?
	fprintf(F, "%2s    ", ind->sex==S_MALE ? "M" : "F");

        // Now write the marker information
     
	for (i=1; i<=numberofmarkers(); i++) {
	  marker = markernumber(ind, order[i]);
	  if (marker->allele1 == 0 || marker->allele2 == 0) {
	    fprintf(F, " %7s", "");
	  }
	  else {
	    //	    fprintf(F, " %3d/%-3d", marker->allele1, marker->allele2);
	    fprintf(F, " %3c/%-3c", 64+marker->allele1, 64+marker->allele2);
	  }
	}

        // Then write the trait information
	for (i=1; i<=numberoftraits(); i++) {
	  if (!traitmiss(ind, i)) 
	    fprintf(F, " %8f", trait(ind, i));
	  else
	    fprintf(F, " %8s", "");

	}
        fprintf(F, "\n");
      }
    }
  }

  freelist(pedlist);  
  fclose(F);
}

void ExportMendelLocus(individual *indlist) {
  int i, j, nAlleles;
  char minibuf[10], tempstr[NAMESIZE];  

  cfopen ("mendel.loc", "w");

  // Make an artificial DISEASE variable
  fprintf(F, "DISEASE AUTOSOME 2 0\n+       .99\n-       .01\n");

  for (i = 0; i < numberofmarkers(); i++) {
    nAlleles = NumberOfAlleles(order[i+1]);

    RemoveWhiteSpace(tempstr, GetName(markernames, order[i+1]));
    strncpy(minibuf, tempstr, 8);
    minibuf[8] = '\0';

    // Prints the locus information
    //    fprintf(F, "%-8s%-8s%2d%2d\n", minibuf,"AUTOSOME", nAlleles, nAlleles*(nAlleles+1)/2);
    fprintf(F, "%-8s%-8s%2d%2d\n", minibuf,"AUTOSOME", nAlleles, 0);

    // Should now print the frequency information

    CopyAlleleFreq(order[i+1]);
    
    for (j=1; j<= NumberOfAlleles(order[i+1]); j++)
      fprintf(F, "%c       %8.5f\n", 64+j, allelefreq[j]);

  }

  fclose(F);

}


void ExportMendelMap() {
  int i;
  char minibuf[10], tempstr[NAMESIZE];

  cfopen ("mendel.map", "w");

  for (i = 1; i<=numberofmarkers(); i++) {
    RemoveWhiteSpace(tempstr, GetName(markernames, order[i]));
    strncpy(minibuf, tempstr, 8);
    minibuf[8] = '\0';
    if (i<numberofmarkers()) {
      fprintf(F, "%-8s\n        %8.5f\n", minibuf, distance[i]);
    }
    else {
      fprintf(F, "%-8s\n", minibuf);
    }
  }

  fclose(F);
  
}




/*

  Linkage

*/


void ExportLinkagePedigree (individual *indlist, int aff)
{
  individual *ind;
  markerlist *marker;
  namelist *pedlist, *pedigr;
  quanttrait *qt;
  int fndr, i;
  

  if (aff==1) 
    cfopen ("gh.ped","w");
  else
    cfopen ("linkage.ped","w");

  pedigr = MakePedigreeList(indlist);

  for (pedlist = pedigr; pedlist; pedlist = pedlist ->next)
  {
    for (fndr = 1; fndr >= 0; fndr--) {
    for (ind = indlist; ind; ind = ind->next)
    {
      // First prints the founders from right pedigree
      if ((founder(ind)==fndr) && (!strcmpl(ind->pedigree, pedlist->name)))
      {
	  // Converts GeneHunter ID number to simple numbers
        if (aff == 1) {
            fprintf(F, "%-7s %-12d", ind->pedigree, ind->globalid);
	    if (ind->father)
	      fprintf(F, "%-12d", ind->father->globalid);
	    else
	      fprintf(F, "%-12d", 0);
	    if (ind->mother)
	      fprintf(F, "%-12d", ind->mother->globalid);
	    else
	      fprintf(F, "%-12d", 0);
	}
	else { 
	  fprintf(F, "%-7s %-12s", ind->pedigree, ind->id);
	  if (ind->father)
	    fprintf(F, "%-12s", ind->father->id);
	  else
	    fprintf(F, "%-12d", 0);
	  if (ind->mother)
	    fprintf(F, "%-12s", ind->mother->id);
	  else
	    fprintf(F, "%-12d", 0);
	}
        fprintf(F, "%-2d", ind->sex);

        if (aff)
          fprintf(F, "1 ");

        for (marker = ind->marker, i=1; marker; marker = marker->next, i++)
        {
#ifdef XGH	 
	  if (options->Index[O_XLINKED] && ind->sex == S_MALE ) {
	    fprintf(F, " %2d %2d ", NumberOfAlleles(i)+1, marker->allele2);
	  }
	  else 
	    fprintf(F, " %2d %2d ", marker->allele1, marker->allele2);
#else
          fprintf(F, " %2d %2d ", marker->allele1, marker->allele2);
#endif
        }
        for (qt = ind->qttrait; qt; qt = qt->next)
        {
          if (qt->sysmiss)
            fprintf(F, " %s ", SYSMISS);
          else
            fprintf(F, " %f ", qt->value);
        }

        fprintf(F, "\n");    
      }
    }
    }
  }

  fclose(F);
  if (aff)
    printf("Created GeneHunter pedigree file\n");
  else
    printf("Created LINKAGE pedigree file\n");

  freelist(pedigrees);
  pedigrees = 0;
}

void ExportLinkageParameter (individual *indlist)
{
  individual *ind;
  FreqList *fl;
  int i, j, k;
  

  // Check that data exists
  if (listlen(indlist)==0)
  {
    printf("WARNING: No data avaliable for creating linkage data file\n");
    return;
  }

  ind = indlist;

  cfopen ("linkage.par", "w");

  fprintf(F, "%d 0 0 3\n", listlen(ind->marker) + numberoftraits());
  fprintf(F, "0 0.0 0.0 0\n");

  // Now prints the marker order
  for (i=1; i<= numberofmarkers(); i++)
  {
    fprintf(F, "%d ",order[i]);
  }
  fprintf(F, "\n");

  for (i = 1; i<= numberofmarkers(); i++)
  {
    // This is a genetic marker with so many alleles (requires reducing)
    j = maxallelenumber(i);
    fprintf(F, "3 %d # %s\n", j,  GetName(markernames, i));

    fl = FrequencyNumber(i);
    for (k=1; k<= fl->num_alleles; k++)
      fprintf(F, "%-6.4f ",fl->frequency[k]);
    fprintf(F, "\n");
  }
  for (i=1; i<=numberoftraits(); i++)
  {
    fprintf(F, "4 0 # %s \n", GetName(traitnames, i));
  }
  fprintf(F, "0 0\n");

  // Printing recombination fractions
  if (numberofmarkers()>1)
  {
    for (i = 1; i< numberofmarkers(); i++)
      fprintf(F, "%-5.3f ",distance[i]);
    fprintf(F, "\n"); 
  }
  else
    fprintf(F, "0\n ");

  fprintf(F, "0\n1\n");

  fflush(F);  
  fclose(F);

  puts("Made LINKAGE parameter file");
}


/*
     Merlin
*/

int ExportMerlinPedigree (individual *indlist)
{
  individual *ind;
  markerlist *marker;
  namelist *pedlist, *pedigr;
  quanttrait *qt;
  int fndr, i, pednr;


  // Checks that data exists
  if (!listlen(indlist))
  {
    return -1;
  }  

  cfopen ("merlin.ped","w");

  pedigr = MakePedigreeList(indlist);

  // Run through each pedigree
  for (pedlist = pedigr, pednr=1; pedlist; pedlist = pedlist ->next, pednr++) {
    // First select all founders then all offspring
    for (fndr = 1; fndr >= 0; fndr--) {
      for (ind = indlist; ind; ind = ind->next)	{
	// First prints the founders from right pedigree
	if ((founder(ind)==fndr) && (!strcmpl(ind->pedigree, pedlist->name))) {
	  // Converts ID's to simple numbers
	  fprintf(F, "%-5d %-12d", pednr, ind->globalid);
	  if (ind->father)
	    fprintf(F, "%-12d", ind->father->globalid);
	  else
	    fprintf(F, "%-12d", 0);
	  if (ind->mother)
	    fprintf(F, "%-12d", ind->mother->globalid);
	  else
	    fprintf(F, "%-12d", 0);
	  fprintf(F, "%-2d", ind->sex);

	  // Print marker info
	  for (marker = ind->marker, i=1; marker; marker = marker->next, i++) {
#ifdef XGH	 
	    if (options->Index[O_XLINKED] && ind->sex == S_MALE ) {
	      fprintf(F, " %2d %2d ", NumberOfAlleles(i)+1, marker->allele2);
	    }
	    else 
	      fprintf(F, " %2d %2d ", marker->allele1, marker->allele2);
#else
	    fprintf(F, " %2d %2d ", marker->allele1, marker->allele2);
#endif
	  }	  


	  // Print traits
	  for (qt = ind->qttrait; qt; qt = qt->next) {
	    if (qt->sysmiss)
	      fprintf(F, " %-6s ", "x");
	    else
	      fprintf(F, " %-6f ", qt->value);
	  }

	  fprintf(F, "\n");    
	}
      }
    }
  }
  fclose(F);
  printf("Created MERLIN pedigree file\n");

  freelist(pedigrees);
  pedigrees = 0;

  return(0);
}


int ExportMerlinDataFile ()
{
  int i;
  char mkrname[NAMESIZE], tempstr[NAMESIZE];

  // Makes pedigree file
  cfopen ("merlin.dat", "w");


  for (i = 1; i<=numberofmarkers(); i++)
  {
    // Remove white spaces in marker name
    RemoveWhiteSpace(tempstr, GetName(markernames, order[i]));
    strncpy(mkrname, tempstr, NAMESIZE);
    fprintf(F, "M %s\n", mkrname);
  }

  for (i = 1; i<=numberoftraits(); i++)
  {
    // Remove white spaces in marker name
    RemoveWhiteSpace(tempstr, GetName(traitnames, i));
    strncpy(mkrname, tempstr, 18);

    fprintf(F, "T %s\n", mkrname);
  }


  fclose(F);
  return 0;
}


int ExportMerlinMap (void)
{
  char mkrname[NAMESIZE], tempstr[NAMESIZE];
  int i;
  double currentpos;

  printf("Which chromosome is this: ");
  InputLine(buf, BUFFERSIZE);

  currentpos = 0;
  // Makes the marker file
  cfopen ("merlin.map", "w");
  fprintf(F, "CHROMOSOME   MARKER          LOCATION\n");
  for (i = 1; i<=numberofmarkers(); i++)
  {
    // Remove white spaces in marker name
    RemoveWhiteSpace(tempstr, GetName(markernames, order[i]));
    strncpy(mkrname, tempstr, NAMESIZE);

    fprintf(F, "%-5s  %-18s  %f\n", buf, mkrname, MarkerDistance(0,i));
  }
  fclose(F);
  puts("Made Merlin map file");
  return 0;
}

int ExportMerlinFreq (void) {
  int i, k;
  FreqList *fl;
  char mkrname[NAMESIZE], tempstr[NAMESIZE];

  cfopen ("merlin.freq", "w");

  for (i = 1; i<=numberofmarkers(); i++) {
    // Remove white spaces in marker name
    RemoveWhiteSpace(tempstr, GetName(markernames, order[i]));
    strncpy(mkrname, tempstr, NAMESIZE);
    fprintf(F, "M %s\nF ", mkrname);

    fl = FrequencyNumber(i);
    for (k=1; k<= fl->num_alleles; k++)
      fprintf(F, "%-6.4f ",fl->frequency[k]);
    fprintf(F, "\n");
  }

  fclose(F);
  puts("Made Merlin frequency file");
  

  return 0;
}




/*-----------------------*/
/* End of Export rutines */
/*-----------------------*/



individual *newperson (void)
{
  individual *ind;

  ind = cmalloc (sizeof (individual));
  memset (ind,0,sizeof (individual));

  ind->sex = S_UNKNOWN;

  if (findperson(ind->id))
  {
    printf("Not possible to create a new person as ID already exists\n");
    return findperson(ind->id);
  }
  else
    addlist (&individuals,ind);
  return ind;
}


/*-----------------------*/

void ReadATrait(quanttrait *traitp)
{
  char *s;

  s = getstr();

  // Checking if it has the sysmiss symbol
  if (!strcmpl(SYSMISS,s))
  {
    traitp->sysmiss = 1;
    traitp->value   = 0;
  }
  else
  {
    traitp->sysmiss = 0;
    traitp->value   = atof(s);
  }

  return;
}


void readinputinfo(char *filename)
{
  
  cfopen (DESCFILE,"r");

  getbuf ();

  while (buf[0] != EOF)
  {

  }
  fclose(F);
}

void importlinkagepedigree (void)
{
  char minibuf[NAMESIZE];
  individual *ind, *ind2;
  

  //  IDlist *idlist;
  markerlist *marker;
  namelist *nl;
  quanttrait *qt;
  int i, families, persons, nmarkers, ntraits, code;

  printf("Number of markers in file: ");
  InputLine(buf, BUFFERSIZE);
  nmarkers = atoip(buf);

  printf("Number of traits in file : ");
  InputLine(buf, BUFFERSIZE);
  ntraits = atoip(buf);

  printf ("Name of linkage pedigree file to input: ");
  InputLine(buf, BUFFERSIZE);
  cfopen (buf,"r");

  getbuf ();

  persons = 0;
  families = 0;

  if (options->Index[O_CONVERTLINKAGEID])
    printf("WARNING: ID's are converted from number to fam-number\n");

  while (buf[0] != EOF)
  {

    // If empty line then read next line
    if (!buf[0] || !strpbrk(buf, "0123456789"))
    {
      getbuf ();
      continue;
    }

    ind = newperson ();
    persons++;
    ind->globalid = persons;
    ind->localid = persons;

    strcpy(ind->pedigree, igetstr (buf));
    strcpy(ind->id, getstr ());

    strcpy(ind->tmpstr1, getstr ()); // Father ID
    strcpy(ind->tmpstr2, getstr ()); // Mother ID

    // If conversion of IDs
    if (options->Index[O_CONVERTLINKAGEID]) {
      sprintf(minibuf,"%s-%s", ind->pedigree,ind->id);
      strcpy(ind->id,minibuf);
      sprintf(minibuf, "%s-%s", ind->pedigree,ind->tmpstr1);
      strcpy(ind->tmpstr1,minibuf);
      sprintf(minibuf, "%s-%s", ind->pedigree,ind->tmpstr2);
      strcpy(ind->tmpstr2,minibuf);
    }
    switch (geti ())
    {
      case  1 : ind->sex = S_MALE; break;
      case  2 : ind->sex = S_FEMALE; break;
      default : ind->sex = S_UNKNOWN;
    }
//    atof(getstr());

    for (i = 0 ; i<nmarkers; i++)
    {
      marker = cmalloc (sizeof (markerlist));
      memset (marker,0,sizeof (markerlist));
      marker->allele1 = geti();
      marker->allele2 = geti ();

      addlist(&ind->marker, marker);
    }

    // Reading traits
    for (i = 0 ; i<ntraits; i++)
    {
      qt = cmalloc (sizeof (quanttrait));
      memset (qt,0,sizeof (quanttrait));
      ReadATrait(qt);

      addlist(&ind->qttrait,qt);
    }

//    printf("\n");
    getbuf ();
  }

  fclose (F);
  printf ("%d persons from %d pedigree(s) read\n", 
          persons, numberofpedigrees());

  code = 0;
  // Should now fix fathers and mothers
  forind
  {
    ind->father = findperson(ind->tmpstr1);
    ind->mother = findperson(ind->tmpstr2);

    // if a person has a father and a mother, then add the
    //  person to the parents lists
    if (ind->father && ind->mother)
    {
      adduniqueidlist (&ind->father->offspring, ind);
      adduniqueidlist (&ind->father->mate, ind->mother);
      adduniqueidlist (&ind->mother->offspring, ind);
      adduniqueidlist (&ind->mother->mate,ind->father);
    }
    else if (ind->father || ind->mother) // The person has only one parent
      {
        code = 1;
        sprintf(buf2, "%s has only one parent in pedigree file\n", ind->id);
        WriteErrorMsg(buf2);
      }
     
  }

  if (code)
    printf("ERROR: A least one parent not found in pedigree file\n");

  // Fixing sibs
  forind
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

  // Have now read the file and made the dataset
  // Should fix the missing parts
  // This part definately needs some error checking

  // First set the marker names:
  for (i=1; i<=nmarkers; i++)
  {
    nl = cmalloc (sizeof (namelist));
    memset (nl,0,sizeof (namelist));

    sprintf(nl->name, "Marker %d",i);
    addlist(&markernames, nl);
  }

  // Then add trait names
  for (i=1; i<=ntraits; i++)
  {
    nl = cmalloc (sizeof (namelist));
    memset (nl,0,sizeof (namelist));

    sprintf(nl->name, "Trait %d",i);
    addlist(&traitnames, nl);
  }

  InitializeFrequencies(nmarkers);

  // Initialize order and distances
  order = ivector(1,nmarkers);
  invorder = ivector(1,nmarkers);
  distance = vector(1, nmarkers-1);
  for (i=1; i<=nmarkers; i++)
  {
    order[i] = i;
    invorder[i] = i;
  }

  for (i=1; i<nmarkers; i++)
    distance[i] = InverseMapFunction(0.1);
  
}

void ReadNextNonEmptyLine(void)
{
  getbuf ();

  while (!buf[0] && !feof(F))
  {
    getbuf ();
  }
}


/*

  Reads a standard LINKAGE parameter file

*/

void FreeDataInfo(void)
{
  freelist(datainfo->allelefrequency);

  memset(datainfo, 0, sizeof(DataDefinitions));
}

char *GetRestOfLine(void)
{
  char *s2, *retstr;
  static char tmpbuf[100];
  int newword;

  tmpbuf[0] = 0;

  retstr = &tmpbuf[0];

  newword = 0;
  s2 = getstr();

  while (*s2)
  {
    if (newword)
      strcat(retstr, " ");
    strcat(retstr, s2);
    s2 = getstr();
    newword = 1;
  }
  //  printf("Resulting : .....|%s|.....\n",retstr);
  return retstr;
}

void ReadLinkageParameterFile(void)
{
  int i;
  int linenr, locinumber;
  int numloci, risk, xlink, program;
  int locus_type, num_alleles, nmarkers;
  char *traitname;
  char tmpbuf[256];
  FreqList *fl;
  namelist *nl;
  

  printf("Name of parameter file: ");
  InputLine(buf, BUFFERSIZE);

  // Deletes old info
  FreeDataInfo();

  DeletePedigreeData();
  DeleteParameterData();

  cfopen (buf,"r");

  linenr = 0;
  locinumber = 1;
  nmarkers = 0;


  // Read line 1
  ReadNextNonEmptyLine();
  numloci = atoip (igetstr (buf));
  risk = geti();
  xlink = geti();
  program = geti();

  if (numloci==0)
  {
    printf("WARNING: No loci in parameter file\n");
    return;
  }

  InitializeFrequencies(numloci);
  freelist(markernames);
  markernames = 0;

  // Read line 2 and 3. Not used yet
  ReadNextNonEmptyLine();

  ReadNextNonEmptyLine();
  // Saves the order for later
  strncpy(tmpbuf, buf, sizeof(tmpbuf));

  // Read all locus/traits
  while (locinumber <= numloci)
  {
    ReadNextNonEmptyLine();

    locus_type = atoip (igetstr (buf));
    num_alleles = geti();
    traitname = GetRestOfLine();

    switch (locus_type)
    {
      case 1 : break;
      case 3 : // Read a marker locus
               fl = FrequencyNumber(locinumber);

               fl->num_alleles = num_alleles;

               // Read in the allele frequencies
	       ReadNextNonEmptyLine();
               fl->frequency[1] = atof(igetstr(buf));
               for (i=2; i<= num_alleles; i++)
               {
                 fl->frequency[i] = getf();
               }

               // Increase the number of markers in the dataset
               nmarkers++;

               // Saves the marker name
               nl = cmalloc(sizeof(namelist));
               memset(nl, 0, sizeof(namelist));
               
               if (strlen(traitname)==0)
                 sprintf(nl->name, "Marker %d",locinumber);
               else
               {
                 // remove initial # if any
                 traitname +=strspn(traitname, "# ") ;
                 strncpy(nl->name, traitname, sizeof(nl->name));
               }
               addlist(&markernames, nl);
               break;
      case 4 : // Read a quantitative trait
               break;
      default: printf("Unknown locus type (%d)\nExiting\n", locus_type);
    }
    locinumber++;

    #ifdef debug
    printf("Found locus (%s) type %d with %d alleles\n", nl->name, locus_type, num_alleles);
    #endif
  }

  ReadNextNonEmptyLine();
  // Now reads the recombination fractions for markers
  ReadNextNonEmptyLine();
  if (nmarkers>1)
  {
    distance = vector(1,nmarkers-1);
    distance[1] = atof(igetstr(buf));
    for (i = 2; i<nmarkers; i++)
      distance[i] = atof(getstr());
  }

  ReadNextNonEmptyLine();
  ReadNextNonEmptyLine();

  fclose(F);

  // Fixes order
  order = ivector(1,nmarkers);
  invorder = ivector(1,nmarkers);
  OrderMarkers(tmpbuf);

}

void MakeGHdatafile (void)
{
  individual *ind;
  FreqList *fl;
  int i, j, k;
  

  // Check that data exists
  if (!individuals)
  {
    printf("WARNING: No data avaliable for creating linkage data file\n");
    return;
  }
  

  ind = individuals;

  cfopen ("gh.data", "w");

  fprintf(F, "%d 0 0 5\n", listlen(ind->marker)+1+numberoftraits());
  fprintf(F, "0 0.0 0.0 0\n");

  // Now prints the marker order
  fprintf(F, "1 ");
  for (i=1; i<= numberofmarkers(); i++)
  {
    fprintf(F, "%d ",order[i]+1);
  }
  fprintf(F, "\n");

  fprintf(F, "1  2  << AFFECTION, NO OF ALLELES\n");
  fprintf(F, "0.99 0.01  << GENE FREQ\n");
  fprintf(F, "1\n");

  fprintf(F, "0.001 0.999 0.999  << GENE FREQ\n");
  //  if (options->Index[O_XLINKED]) {
  //    fprintf(F, "0.001 0.999   << GENE FREQ\n");
  //  }



  for (i = 1; i<= numberofmarkers(); i++)
  {
    // This is a genetic marker with so many alleles (requires reducing)
    j = NumberOfAlleles(i);
#ifdef XGH
    if (options->Index[O_XLINKED]) {
      fprintf(F, "3 %d\n", j+1);
    }
    else {
      fprintf(F, "3 %d\n", j);
    }
#else
    fprintf(F, "3 %d\n", j);
#endif

    fl = FrequencyNumber(i);
    for (k=1; k<= fl->num_alleles; k++)
      fprintf(F, "%-6.4f ",fl->frequency[k]);
#ifdef XGH
    if (options->Index[O_XLINKED]) {
      fprintf(F, "0.001 ");
    }
#endif    
    fprintf(F, "\n");
  }
  for (i=1; i<=numberoftraits(); i++)
  {
    fprintf(F, "4 0 # %s \n", GetName(traitnames, i));
  }

  fprintf(F, "0 0\n");

  // Printing recombination fractions
  if (numberofmarkers()>1)
  {
    fprintf(F, ".1 ");
    for (i = 1; i< numberofmarkers(); i++)
      fprintf(F, "%-7.5f ",distance[i]);
    fprintf(F, "\n"); 
  }
  else
    fprintf(F, "0\n ");

  fprintf(F, "0\n1\n");

  fflush(F);  
  fclose(F);

  puts("Made GeneHunter 2 parameter file");
}


/*

 General File operations

*/

// 1 if file exists, 0 otherwise
int FileExists(char *name)
{
  FILE *F;
  
  if ( (F = fopen (name,"r")) == 0)
    return 0;

  fclose(F);
  return 1;
}

#ifdef debug
#undef debug
#endif

#endif
