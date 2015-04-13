/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * Functions for two and threepoint linkage analysis
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


#ifndef _LINKAGE_H_
#define _LINKAGE_H_

#define MAX_IT      1000      /* maximum number of iterations */
#define ALPHA       1.0       /* reflection coefficient */
#define BETA        0.5       /* contraction coefficient */
#define GAMMA       2.0       /* expansion coefficient */


#define LINKAGE_PRECISION 0.001
#define LINKAGE_SIGNIFICANCE 0.05

#ifdef GCC
double LogTwopointLikelihood(individual *indlist1, individual *indlist2, int markernum, int markernum2, double theta);
double TwopointLinkage(individual *indlist, int marker1, int marker2, double theta);
double TwopointLinkageAnalysis(individual *indlist, int marker1, int marker2, double theta);

double LogThreepointLikelihood(individual *indlist1, individual *indlist2, individual *indlist3, int markernum, int markernum2, int markernum3, double gamma01, double gamma10, double gamma11);
double ThreepointLinkage(individual *indlist, int marker1, int marker2, int marker3, double gamma01, double gamma10, double gamma11);
double ThreepointLinkageAnalysis(individual *indlist, int marker1, int marker2, int marker3, double gamma01, double gamma10, double gamma11);

double simplex(double (*func)(double[]), double start[],int n, double EPSILON,
	       double scale);
double rosen(double x[]);
#else

extern "C" {
#include "pedigree.h"
#include "lists.h"
#include "distr.h"

  extern "C" int usedmap;

  extern "C" double MapFunction(double rf);

  extern "C" int NumberOfAlleles(int markernum);
  extern "C" FreqList *FrequencyNumber (int n);

  extern double allelefreq[MAXALLELES];
  individual *FindIndListMember(individual *l, char *id);
  individual *GenotypeElimination(char *pedigree, int markernum); 
  markerlist *markernumber (individual *ind, int n); 
  namelist *MakePedigreeList(individual *indlist);
  void FreeIndividualList(individual *indlist);
  void FreeNameList(namelist *list);
  void CopyAlleleFreq(int marker);

  double TwopointLinkageAnalysis(individual *indlist, int marker1, int marker2, double theta);
  double LogTwopointLikelihood(individual *indlist1, individual *indlist2, int markernum, int markernum2, double theta);
  double TwopointLinkage(individual *indlist, int marker1, int marker2, double theta);


  double LogThreepointLikelihood(individual *indlist1, individual *indlist2, individual *indlist3, int markernum, int markernum2, int markernum3, double gamma01, double gamma10, double gamma11);
  double ThreepointLinkage(individual *indlist, int marker1, int marker2, int marker3, double gamma01, double gamma10, double gamma11);
  double ThreepointLinkageAnalysis(individual *indlist, int marker1, int marker2, int marker3, double gamma01, double gamma10, double gamma11);

double simplex(double (*func)(double[]), double start[],int n, double EPSILON,
	       double scale);
  double rosen(double x[]);

}
#endif


#endif

