/*
 * PEDIPET
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
#include "maximize.h"

extern "C" {
#include "lists.h"
}

//extern "C" optionlist *options;

// Rutines from pedipet used in estimate
extern "C" {
int founder(individual *ind);                     // Is person a founder
double trait(individual *ind, int n);
quanttrait *traitvalue(individual *ind, int traitnumber);
int traitmiss(individual *ind, int traitnumber);
namelist *MakePedigreeList(individual *indlist);
individual *SelectPedigree(char *pedname, individual *indlist);
void FreeIndividualList(individual *indlist);
individual *GenotypeElimination(char *pedigree, int markernum);
markerlist *markernumber (individual *ind, int n);
double ChromosomeLength ();
  int RelationType(individual *ind, individual *ind2);
  int homozygous (individual *ind, int n);
  int heterozygous (individual *ind, int n);
  int IsTraitMissing(individual *ind, int n);
  int IsGenotyped(individual *ind);
  void CopyAlleleFreq(int marker);
  void FreeNameList(namelist *list);
  int numberofmarkers (void);
  int WriteErrorMsg(char *s);
  int atoip (char *s);
  individual *FindIndListMember(individual *l, char *id);
  individual *FindIndListMemberLocalID(individual *l, int localid);
}

int RelationType(individual *ind, individual *ind2);

double Heritability(individual *indlist, int traitnum, int Dominance, MATRIX X[], int nCovar);
double SingleMarkerVC(individual *indlist, int traitnum, int Dominance, MATRIX IBD[], MATRIX X[], int nCovar);
double SingleXMarkerVC(individual *indlist, int traitnum,  MATRIX IBD[], MATRIX X[], int nCovar, int Dominance, int Household, int *covnum, int markernum);
double MultiMarkerVC(individual *indlist, int traitnum,  MATRIX IBD[], MATRIX X[], int nCovar, int Dominance, int Household);
double MultiXMarkerVC(individual *indlist, int traitnum,  MATRIX IBD[], MATRIX X[], int nCovar, int Dominance, int Household, int *covnum);
double MaximizePolygenicXVC(individual *indlist, int traitnum, MATRIX X[], int nCovar, int Dominance, int Household, int *covnum);
double MaximizePolygenicVC(individual *indlist, int traitnum,  MATRIX X[], int nCovar, int Dominance, int Household, int *covnum);
double MultiGXE(individual *indlist, int traitnum,  MATRIX IBD[], MATRIX X[], MATRIX XX[], int nCovar, int Dominance, int Household, int *covnum);

double SpecialFunction(individual *indlist, int traitnum, double position, double ibd[], double k2[]);
double SpecialFunction2(individual *indlist, int traitnum, double position, MATRIX ibd[], MATRIX k2[]);
double SpecialFunction3(individual *indlist, int traitnum, int position, MATRIX RELAP[]);
double SpecialFunctionPiHat(individual *indlist, int traitnum, int markernum);

double CalculateMixLLSibs(individual *indlist, MATRIX RELAP[], MATRIX resp[], MATRIX param, int model);
double SpecialFunctionSibsMix(individual *indlist, int traitnum, int position, MATRIX RELAP[]);
double SpecialFunctionSibsPiHat(individual *indlist, int traitnum, int model);
MATRIX *SpecialDerivSibs(MATRIX &param, individual *indlist, MATRIX RELAP[], MATRIX resp[], int model);

MATRIX MMakeKinshipMatrix(individual *indlist);
MATRIX MMakeXLinkedKinshipMatrix(individual *indlist);
MATRIX MMakeXLinkedKinshipMatrixBySex(individual *indlist, int sex1, int sex2);
MATRIX MMakeDelta7Matrix(individual *indlist, const MATRIX kinship);
MATRIX MakeIdentMatrix(int nObs);
MATRIX CompleteCases(individual *indlist, int traitnum);

double CalculateNewMixLLNucl(individual *indlist, MATRIX RELAP[], MATRIX resp[], MATRIX param, int model);

double SpecialFunctionMixSibs(individual *indlist, int traitnum, int markernum, int model);

extern "C" void AnalyzeMatrix(individual *indlist, char *filename, int traitnum);
extern "C" void AnalyzeAllMatrices(individual *indlist, char *filename, int traitnum);
extern "C" void AnalyzeMixtureModels(individual *indlist, char *filename, int traitnum);
extern "C" void AnalyzeMixtureModelsNULL(individual *indlist, char *filename, int traitnum);

extern "C" MATRIX* ReadSolarIBDFile(individual *indlist, char *filename, char *indexname);

