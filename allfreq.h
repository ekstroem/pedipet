/*
 * PEDIPET
 * 
 * allfreq.h
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


#include "pedigree.h"
#include "matrix.h"

extern "C" optionlist *options;

extern "C" {
#include "lists.h"
  int founder(individual *ind);
  int heterozygous (individual *ind, int n);
  void CopyAlleleFreq(int marker);
  int NumberOfAlleles(int markernum);
  int istypedformarker(individual *ind, int markernum);
  void RemoveIndividual(individual *indlist, individual *ind) ;
  void FreeNameList(namelist *list);
  void FreeIndividualList(individual *indlist);

}

extern "C" individual *GenotypeElimination(char *pedigree, int markernum);
extern "C" markerlist *markernumber (individual *ind, int n);
extern "C" int listlen (void *l);
extern "C" individual *CopyIndividualList(individual *indlist);
extern "C" namelist *MakePedigreeList(individual *indlist);
extern "C" double allelefreq[MAXALLELES];
extern "C" FreqList *FrequencyNumber (int n);
extern "C" void ShowMarkerInfo (int markernum);
extern "C" double ReduceGenotypes(individual *indlist);
extern "C" void *cmalloc (int n);

#ifdef GCC

void MLAlleleFrequencyEstimation(individual *listofpersons, int markernum);

#else

extern "C" void MLAlleleFrequencyEstimation(individual *listofpersons, int markernum);


#endif
