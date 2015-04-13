/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * Functions to check for errors
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

#ifndef _ERRORCHK_H_
#define _ERRORCHK_H_

#include "matrix.h"

// #include "pedipet.h"
#include "pedigree.h"
#include "distr.h"

#ifdef __cplusplus
extern "C" {
#endif


  
  // Fra andre filer
  int numberofmarkers (void);
  void CopyAlleleFreq(int marker);
  int NumberOfAlleles(int markernum);
  int founder(individual *ind);
  int IsGenotyped(individual *ind);
  int istypedformarker(individual *ind, int markernum);
  markerlist *markernumber (individual *ind, int n);
  void TextColour(int colour);
  
  
  extern double allelefreq[MAXALLELES];
  
  void VerifyIndependence(individual *indlist, int all);
  MATRIX Paper10Function(int a, int b, int c, int d);

#ifdef __cplusplus
}
#endif

#endif
