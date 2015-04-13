/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * Functions for genotype and allele frequency estimation
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


void FastAlleleFrequencyEstimation(individual *indlist, int markernum) 
{
  namelist *pedigreelist, *ped;
  individual *genolist, *ind;
  markerlist *mkr, *mkr2, *mkr3;
  int nAlleles, all_founders_known, nTemp, i;
  double *frequencies, *tempfrequencies, res, checkcombinations;
  unsigned long combinations;


  nAlleles = NUmberOfAlleles(markernum);

  frequencies = new double[nAlleles];    // Used for keeping the results
  tempfrequencies = new double[nAlleles];    // Used for keeping the results

  for (i = 0; i< nAlleles; i++) {
    frequencies[i] = 0;
  }


  pedigreelist = MakePedigreeList(indlist);
  for (ped = pedigreelist; ped; ped = ped -> next) {

    genolist = GenotypeElimination(ind->pedigree, markernum);
    // Assume only one possible genotype exists for each founder
    all_founder_known = 1;

    // Start by speeding it up by a factor 2 for founders
    // Must remember to multiply by two later on
    for (ind = genolist; ind; ind = ind->next) {
      if (founder(ind)) {
	// Remove one copy of ordered heterozygous markers for founders
	if (listlen(ind->marker>1)) {
	  for (mkr = ind->marker; mkr->next; mkr = mkr->next) {
	    for (mkr2 = mkr->next; mkr2; ) {
	      mkr3 = mkr2->next;
	      // If they are identical then remove
	      if ((mkr2->allele1 == mkr->allele1 && mkr2->allele2 == mkr->allele2) || (mkr2->allele1 == mkr->allele2 && mkr2->allele2 == mkr->allele1)) {
		// Remove one copy
		choplist(ind->marker, mkr2);
	      }	      
	      mkr2 = mkr3;
	    }	  
	  }
	}
	if (listlen(ind->marker)>1) {
	  all_founders_known = 0;
	}
      }
    }   
    // Could do something like the above for pedigree ends (nf without offspring)


    // Start calculating the likelihood of the pedigree
    res = 0;

    if (all_founders_known) {
      // Easy case.
    }

    // Sets up the smart array
    int *smartarray, *maxarray;
    smartarray = new int[listlen(genolist)];
    maxarray = new int[listlen(genolist)];

    combinations = 1;
    checkcombinations = 0;
    for(ind = genolist, i=0; ind; ind = ind->next, i++) {
      nTemp = listlen(ind);
      combinations *= nTemp;
      checkcombinations += log(nTemp);      
      maxarray[i] = nTemp;
      smartarray[i] = 1;
      ind->tmpint4 = i+1;
    }
    if (log(combinations) != checkcombinations) {
      cout << "Overflow error in calculation of combinations. Too complicated" << endl;
    }

    for (i = 0; i< nAlleles; i++) {
      tempfrequencies[i] = 0;
    }

    res = 0;
    tempres = 0;
    weight = 0;
    for (i = 0 ; i < combinations; i++) {
      for (ind = genolist; ind; ind = ind->next) {
	if (founder(ind)) {
	  mkr = markernumber(ind, smartarray[ind->tmpint4]);
	  tempres = tempres*(double)param[mkr->allele1]*(double)param[mkr->allele2];
	}
      }

      
      
    }


    
    

  }

  freelist(genolist);
  freelist(pedigreelist);

  delete[] frequencies;

}
