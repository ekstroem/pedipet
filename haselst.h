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

#ifndef HASEMANELSTON

#include "pedigree.h"
#define EXPORT
// #define CHI2                // Print Chi^2 instead of -t. Doesn't work
#define HASEMANELSTON
// #define NOMULTI
// #define debug

typedef struct HasemanElstonData
{
  struct HasemanElstonData *next;
  float traitdiff;
  float ibd;
  quanttrait *covariates;

} HasemanElstonData;

typedef struct phenoset
{
  struct phenoset *next;
  int allele1;
  int allele2;
} phenoset;


double k2;


double SingleMarkerIBD(individual *ind, individual *ind2, int markernum);

#endif
