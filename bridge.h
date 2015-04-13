/*
 * PEDIPET
 * 
 * bridge.h
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

extern double k2;
extern optionlist *options;

#ifdef __cplusplus
extern "C" {
individual *individuals;
int fullsibs(individual *ind1, individual *ind2);
void MakeGHdatafile (void);
int numberofmarkers();
void ExportLinkagePedigree (individual *indlist, int aff);
int FileExists(char *name);
quanttrait *traitvalue(individual *ind, int traitnumber);
void freelist (void *p1);
int FamilySize(individual *indlist, char *pedid);
double ChromosomeLength(void);
double CalculateHeritability(individual *indlist, int traitnum, int Dominance);
double Paper1Function(individual *indlist, int traitnum);
double Paper1Function2(individual *indlist, int traitnum);
double ComplexMixed(individual *indlist, int traitnum);
double SibsComplexMixed(individual *indlist, int traitnum, int markernum);
void FullSibsIBD(individual *indlist, char *filename, int step);
int MMakeKinshipMatrix(individual *indlist, char *filename);

void CalcPIHatByPedigree(individual *indlist, int markernum, char *filename);
}

#else

void FullSibsIBD(individual *indlist, char *filename, int step);
int MMakeKinshipMatrix(individual *indlist, char *filename);
double Paper1Function(individual *indlist, int traitnum);
double CalculateHeritability(individual *indlist, int traitnum, int Dominance);
double Paper1Function2(individual *indlist, int traitnum);
double SinglePointVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov);
double SinglePointXVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov);
double MultiPointXVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov);
double PolygenicXVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov);
double PolygenicVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov);
double MultiPointVC(individual *indlist, int traitnum, int markernum, int numcov, int *inclcov);
void DataMining();

double MultiPointGXE(individual *indlist, int traitnum, int envnum, int markernum, int numcov, int *inclcov);

#endif



