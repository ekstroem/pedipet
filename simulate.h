/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * This file is part of the statistical program Pedipet.
 * Functions to simulate genetic marker data
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


#ifndef _SIMULATE_H_
#define _SIMULATE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "ranlib.h"
#include "ibd.h"
#include "pedipet.h"
#include "lists.h"
  // #include "pedigree.h"

//
// Global variables
//


individual *simulateQTL(int families, int nooff, double dist, double addeffect, double domeffect, double prev);

int roll_1Dx(int n);

individual *SimulateChromosome(int families, int nooff, double location, double addeffect, double domeffect, double prev, double resadd, double resdom);
individual *SimulateKosambiChromosome(int families, int nooff, double location, double addeffect, double domeffect, double prev, double resadd, double resdom);
individual *SimulateSpecialChromosome(int families, int nooff, double location, double addeffect, double domeffect, double prev);
individual *SimulateSpecialChromosome2(int families, int nooff, double location, double addeffect, double domeffect, double prev);
individual *SimulateDataset(individual *indlist, double location, double addeffect, double domeffect, double prev, double resadd, double resdom);
individual *SimulateXDataset(individual *indlist, double location, double addeffect, double domeffect, double prev, double resadd, double resdom);
individual *SimulateStenoDataset(individual *indlist, double location, double addeffect, double domeffect, double prev, double resadd, double resdom, int probandtrait);
individual *SimulateMarkerSet(individual *indlist, int markertype, int nmarkers, double markerdist, double chromosomelength);
individual *SimulateGenotypes(individual *indlist);
individual *SimulatePedigreeStructure(int nPedigrees, int nPedType, int nTraits, int nQTL);
individual *CloneDataset(individual *indlist, int nclones);

void GenotypeError(individual *ind, int markernum, int errorfct);
void GenotypeErrorMendelian(individual *ind, int markernum, int *changed);
void RandomGenotypeError(individual *indlist, int markernum);

void RemovePhenotypeData(individual *ind, int traitnum);
void RemoveGenotypeData(individual *ind, int markernum);
void SimulateRandomMissingGenotypes(individual *indlist, int markernum, double freq);

#ifdef __cplusplus
}
#endif


#endif
