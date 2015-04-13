/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
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



#include <cstdlib>
#include "model.h"

// Constructor
CModel::CModel()
{
  m_nTraits = 0;
  m_nCovariates = 0;

  m_traitnumbers = NULL;
  m_covariatenumbers = NULL;
}


CModel::CModel(int ntraits, int traits[], int ncov, int cov[])
{
  int i;

  m_intercept = 1;
  m_nTraits = ntraits;
  m_nCovariates = ncov;

  m_traitnumbers = new int[ntraits];
  m_covariatenumbers = new int[ncov];

  for (i = 0; i< ntraits; i++) {
    m_traitnumbers[i] = traits[i];
  }

  for (i = 0; i< ncov; i++) {
    m_covariatenumbers[i] = cov[i];
  }
}

CModel::~CModel()
{
  if (m_traitnumbers != NULL) {
    delete[] m_traitnumbers;
    m_traitnumbers = NULL;
  }
  if (m_covariatenumbers != NULL) {
    delete[] m_covariatenumbers;
    m_covariatenumbers = NULL;
  }

}


void CModel::Clear() {
  m_intercept = 1;
  if (m_traitnumbers != NULL) {
    delete[] m_traitnumbers;
    m_traitnumbers = NULL;
  }
  if (m_covariatenumbers != NULL) {
    delete[] m_covariatenumbers;
    m_covariatenumbers = NULL;
  }
}


int CModel::nTraits() {
  return m_nTraits;

}


int AddTrait(CModel model) 
{
  int *newdata = new int[model.nTraits()+1];
  //  int i;

  delete[] newdata;
  return 0;
}
