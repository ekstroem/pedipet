/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * This file is part of the statistical program Pedipet.
 * The functions in maximize are for maximizing mixed models
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

#ifndef MAXIMIZE_H
#define MAXIMIZE_H

#include "matrix.h"
//#include "pedigree.h"

#define MIXEDMODELMAXITER   50        // Max number of iterations
#define MAXSCORINGITER      1         // Max number of iterations
#define EMMAXITER           1300       // Max number of EM-iterations

enum
{
  METHOD_ML,
  METHOD_REML,
  MAX_METHOD
};


/*

  MaximizeMixedModel maximizes a MM under the constraints, that the VC should be non-negative.
  No missing data should be present in the input data

  Input    : response  - an n*1 vector of responses
             design    - an n*p design matrix 
             nVarParam - # of variance parameters (k)
             VC        - a linked list of the (k) variances

             start

  Returns the negative log likelihood value in the minimum
             
  (C) Claus Ekstrøm 1999

*/

double MaximizeMixedModel(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, int Method = METHOD_ML, int Constrain = 1, int PrintInfo = 1);

double MaximizeMixedModelX(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, int Method = METHOD_ML, int Constrain = 1, int PrintInfo = 1);

double MaximizeMixtureModel(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, double p, MATRIX include, int Method, int Constrain, int PrintInfo);

double MaximizeModifiedMixtureModel(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, double p, MATRIX include, int Method, int Constrain, int PrintInfo);
double MaximizeModifiedMixtureModelNULL(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, double p, MATRIX include, int Method, int Constrain, int PrintInfo);

double MaximizeRestrictedMixtureModel(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], MATRIX start, MATRIX beta, double p, int Method, int Constrain, int PrintInfo);

double MaximizeMixedGXEModel(int nPedigrees, int nPedSize[], MATRIX y[], MATRIX x[], int nVC, MATRIX VC[], int nGXEparam, MATRIX XX[], MATRIX start, MATRIX beta, int Method, int Constrain, int PrintInfo, double corfct);

extern "C" double RES1, RES2, RES3;

#endif


