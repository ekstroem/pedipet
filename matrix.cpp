/*
 * PEDIPET
 * 
 * Copyright (C) 1998--2004 Claus Ekstrøm
 *
 * (pedipet@ekstroem.dk)
 *
 * matrix functions
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



#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cassert>

#include "matrix.h"

#define CHECK_BOUNDS

/*

  This file contains standard mathematical functions.

  (C) Claus Ekstr°m 1996-1999

*/


// Constructor
MATRIX::MATRIX()
{
  m_nRows = 0;
  m_nCols = 0;
  m_nMatrixType = MAT_GENERAL;  

  data = NULL;
}


MATRIX::MATRIX(int rows, int cols)
{
  m_nRows = rows;
  m_nCols = cols;
  m_nMatrixType = MAT_GENERAL;  

  data = new double[m_nRows*m_nCols];
 
  // Bør nulstille data
  // memset(data, 0, sizeof(data));
  // Sets all variables to 0
  (*this) = 0.0;

//  cout << "End of constructor at address " << data << "\n";
}


MATRIX::~MATRIX()
{
  //Print();

//  cout << "End of destructor at address " << data;

  if (data != NULL)
  {
    delete[] data;
    data = NULL;
  }

//  cout <<  ". New data: " << data << "\n";
}

int MATRIX::Rows(void) const
{
  return m_nRows;
}

int MATRIX::Cols(void) const
{
  return m_nCols;
}


// Matrix scaling
MATRIX MATRIX::operator*(double s)
{
  int i;
  MATRIX res(m_nRows,m_nCols);

  res.m_nMatrixType = m_nMatrixType;

  for (i=0; i<(m_nCols*m_nRows); i++)
  {
    res.data[i] = data[i]*s;
  }

  return res;
}


// Add constant
MATRIX MATRIX::operator+(double s)
{
  int i;
  MATRIX res(m_nRows,m_nCols);

  res.m_nMatrixType = m_nMatrixType;

  for (i=0; i<(m_nCols*m_nRows); i++)
  {
    res.data[i] = data[i]+s;
  }

  return res;
}

// Subtract constant
MATRIX MATRIX::operator-(double s)
{
  int i;
  MATRIX res(m_nRows,m_nCols);

  res.m_nMatrixType = m_nMatrixType;

  for (i=0; i<(m_nCols*m_nRows); i++)
  {
    res.data[i] = data[i]-s;
  }

  return res;
}

/*
MATRIX MATRIX:(double s)*operator
{
  int i;
  MATRIX res(m_nRows,m_nCols);

  res.m_nMatrixType = m_nMatrixType;

  for (i=0; i<(m_nCols*m_nRows); i++)
  {
    res.data[i] = data[i]*s;
  }

  return res;
}
*/


// Assignment operator
// Assigns all elements in a matrix to the value s
MATRIX& MATRIX::operator=(double s)
{
  int i, j;

  for (i = 1; i<=m_nRows; i++)
  {
    for (j = 1; j<=m_nCols; j++)
    {
       (*this)(i,j) = s; 
    }
  }
  return (*this);
}

// Assign a matrix 
MATRIX& MATRIX::operator=(const MATRIX &s)
{
#ifdef CHECK_BOUNDS
  assert(s.m_nCols == m_nCols);
  assert(s.m_nRows == m_nRows);
#endif

  int i, j;

  for (i=1; i<=Rows(); i++) {
    for (j=1; j<=Cols(); j++) {
      (*this)(i,j) = s(i,j);
    }
  }

  return *this;
}

MATRIX::MATRIX(const MATRIX &m)		   // Copy operator
{
  m_nRows = m.Rows();
  m_nCols = m.Cols();
  m_nMatrixType = m.m_nMatrixType;

  int i, j;
  data = new double[m_nCols*m_nRows];

  for (i = 1; i<=m_nRows; i++) {
    for (j = 1; j<=m_nCols; j++) {
       (*this)(i,j) = m(i,j); 
    }
  }
}


// Returns a pointer to a certain element in the matrix
inline double& MATRIX::operator()(int row, int col)
{
#ifdef CHECK_BOUNDS
    assert(row>=1);
    assert(row<=m_nRows);
    assert(col>=1);
    assert(col<=m_nCols);
#endif

  return(data[ m_nCols*(row-1)+col-1]);
}

inline double& MATRIX::operator()(int row, int col) const
{
#ifdef CHECK_BOUNDS
    assert(row>=1);
    assert(row<=m_nRows);
    assert(col>=1);
    assert(col<=m_nCols);
#endif

  return(data[ m_nCols*(row-1)+col-1]);
}


void MATRIX::Print(void) const
{
  for (int i = 1; i<=m_nRows; i++)
  {
    for (int j = 1; j<=m_nCols; j++)
    {
      printf("%-5.3f  ", (*this)(i,j));
    }
    printf("\n");
  }
}

int MATRIX::FPrint(char *fname) const
{

  FILE *F;

  F = fopen(fname, "w");

  if (F != NULL)
  {
    for (int i = 1; i<=m_nRows; i++)
    {
      for (int j = 1; j<=m_nCols; j++)
      {
        fprintf(F, "%-7.4f  ", (*this)(i,j));
      }
      fprintf(F, "\n");
    }
    fclose(F);
    return 0;
  }
  return -1;
}


int MATRIX::FileWriteSymmetric(char *fname)
{
  FILE *F;

  F = fopen(fname, "w");

  if (F != NULL)
  {
    for (int i = 1; i<=m_nRows; i++)
    {
      for (int j = 1; j<=i; j++)
      { 
	if ((*this)(i,j) != 0.0)
	  fprintf(F, "%5d %5d  %9.7f\n", i, j, (*this)(i,j));
      }
    }
    fclose(F);
    return 0;
  }
  return -1;
}

int MATRIX::FileReadSymmetric(char *fname, int nCol, int nRow)
{
  FILE *F;
  int maximum,i, j;
  float res;

  F = fopen(fname, "r");

  if (F != NULL)
  {
    maximum = 0;
    while (fscanf(F, "%5d %5d %9f", &i, &j, &res) != EOF)
    {
      // First figure out the dimensions
      if (i>maximum)
        maximum = i;
    }
    fclose(F);
  }
  else {
    printf("Error: File %s not found\n", fname);
    exit(1);
  }

  F = fopen(fname, "r");

  if (F != NULL && maximum>0)
  {
    //    (*this).Resize(max(nRow, maximum),max(nCol, maximum));
    (*this).Resize(maximum, maximum);

    (*this) = 0.0;

    while (fscanf(F, "%5d %5d  %9f\n", &i, &j, &res) != EOF)
    {
      (*this)(i,j) = res;
      (*this)(j,i) = res;
    }
  }
  fclose(F);

  return 0;
}


MATRIX& MATRIX::Resize(const int rows, const int cols)    // Clears and resizes matrix
{
  delete[] data;

  m_nRows = rows;
  m_nCols = cols;
  m_nMatrixType = MAT_GENERAL;  

  data = new double[m_nRows*m_nCols];
 
  // B†r nulstille data
  // memset(data, 0, sizeof(data));
  // Sets all variables to 0
  (*this) = 0.0;

  return *this;
}


MATRIX& MATRIX::Insert(const int rows, const int cols, const MATRIX& a)
{
  int i, j;

  assert(a.Rows()+rows <= m_nRows);
  assert(a.Cols()+cols <= m_nCols);

  for (i = 0; i< a.Rows(); i++)
  {
    for (j = 0; j< a.Cols(); j++)
    {
      (*this)(rows + i, cols+j) = a(i+1,j+1);
    }
  }
  return *this;
}




MATRIX MATRIX::Subset(const MATRIX& a)
{
  assert(a.Rows() == m_nRows);
  assert(a.Cols() == 1);

  int nDim1, i, j, k, l;

  nDim1 = 0;
  for (i = 1; i<= a.Rows(); i++) {
    if (a(i,1))
      nDim1++;
  }

  MATRIX res(nDim1,m_nCols);

  for (i=1, k=0; i<=a.Rows(); i++) {
    if (!a(i,1)) 
      continue;
    k++;
    for (j=1, l=0; j<=m_nCols; j++) {
      res(k,j) = (*this)(i,j);
    }
  }

  return res;
}





MATRIX MATRIX::Subset(const MATRIX& a, const MATRIX& b)    
{
  assert(a.Rows() == m_nRows);
  assert(b.Rows() == m_nCols);
  assert(a.Cols() == 1);
  assert(b.Cols() == 1);

  int nDim1, nDim2, i, j, k, l;

  nDim1 = 0;
  for (i = 1; i<= a.Rows(); i++) {
    if (a(i,1))
      nDim1++;
  }
  nDim2 = 0;
  for (i = 1; i<= b.Rows(); i++) {
    if (b(i,1))
      nDim2++;
  }

  MATRIX res(nDim1,nDim2);

  for (i=1, k=0; i<=a.Rows(); i++) {
    if (!a(i,1)) 
      continue;
    k++;
    for (j=1, l=0; j<=b.Rows(); j++) {
      if (!b(j,1))
	continue;
      l++;
      res(k,l) = (*this)(i,j);
    }
  }

  return res;
}





double MATRIX::Sum(void)
{
  double res = 0.0;
  long i;

  for (i=0; i<m_nRows*m_nCols; i++)
    res += data[i];

  return(res);

}

double MATRIX::Min(void)
{
  double res = data[0];
  long i;

  for (i=0; i<m_nRows*m_nCols; i++)
    if (data[i]<res)
      res = data[i];

  return(res);

}

double MATRIX::Max(void)
{
  double res = data[0];
  long i;

  for (i=0; i<m_nRows*m_nCols; i++)
    if (data[i]>res)
      res = data[i];

  return(res);

}

double MATRIX::SSD(const double mean)
{
  double res = 0.0;
  long i;

  for (i=0; i<m_nRows*m_nCols; i++)
    res += (data[i]-mean)*(data[i]-mean);

  return(res);

}


// ----------------------------------------------------------------
// Now comes overloading of +, - operators for various combinations
// ----------------------------------------------------------------

// General matrix addition
MATRIX operator+(const MATRIX& m, const MATRIX& n)
{

#ifdef CHECK_BOUNDS
  // Check that matrices conform
  assert(m.Cols() == n.Cols());
  assert(m.Rows() == n.Rows());
#endif


  int i, j;

  MATRIX res(m.Rows(),m.Cols());

  res.m_nMatrixType = m.m_nMatrixType;

  for (i=1; i<=m.Rows(); i++)
  {
    for (j=1; j<=m.Cols(); j++)
      res(i,j) = m(i,j) + n(i,j);
  }

  return res;
}

/*
MATRIX operator+=(MATRIX& m, const MATRIX& n)
{
  int i,j;

  // Check that matrices conform
  assert(m.Cols() == n.Cols());
  assert(m.Rows() == n.Rows());

  for (i=1; i<=m.Rows(); i++)
    for (j=1; j<=m.Cols(); j++)
       m(i,j) += n(i,j);
}

*/

MATRIX MATRIX::operator+=(const MATRIX& m)
{
  int i,j;

  //  (*this).Print();

#ifdef CHECK_BOUNDS
  // Check that matrices conform
  assert(m.Cols() == m_nCols);
  assert(m.Rows() == m_nRows);
#endif

  for (i=1; i<=m_nRows; i++)
    for (j=1; j<=m_nCols; j++)
       (*this)(i,j) += m(i,j); 
  return (*this);
}


MATRIX MATRIX::operator*=(const double a)
{
  int i,j;

  for (i=1; i<=m_nRows; i++)
    for (j=1; j<=m_nCols; j++)
       (*this)(i,j) *= a; 
  return (*this);
}



// General matrix subtraction
MATRIX operator-(const MATRIX& m, const MATRIX& n)
{

  #ifdef CHECK_BOUNDS
  // Check that matrices conform
  assert(m.Cols() == n.Cols());
  assert(m.Rows() == n.Rows());
  #endif

  int i, j;

  MATRIX res(m.Rows(),m.Cols());

  res.m_nMatrixType = m.m_nMatrixType;

  for (i=1; i<=m.Rows(); i++)
  {
    for (j=1; j<=m.Cols(); j++)
      res(i,j) = m(i,j) - n(i,j);
  }

  return res;
}

// General matrix multiplication
MATRIX operator*(const MATRIX& m, const MATRIX& n)
{

  #ifdef CHECK_BOUNDS
  // Check that matrices conform
  assert(m.Cols() == n.Rows());
  #endif

  MATRIX res(m.Rows(),n.Cols());
  double tmp;

  res.m_nMatrixType = m.m_nMatrixType;

  for (int i = 1; i<=m.Rows(); i++)
  {
    for (int j = 1; j<=n.Cols(); j++)
    {
      tmp = 0.0;
      for (int k = 1; k<=m.Cols(); k++)
      {
        tmp += m(i,k)*n(k,j);
      }
      res(i,j) = tmp;
    }
  }
  return res;
}


double Trace(const MATRIX& m)
{
  #ifdef CHECK_BOUNDS
  // Checking square matrix
  assert(m.Cols() == m.Rows());
  #endif

  double res = 0.0;
 
  for (int i=1; i<=m.Cols(); i++)
  {
    res += m(i,i);
  }
  return res;
}

//
//  Returns the LOWER triangular matrix of the cholesky decomposistion
//
MATRIX Cholesky(const MATRIX &m)
{
  #ifdef CHECK_BOUNDS
  // Check that matrix is square
  assert(m.Cols() == m.Rows());
  #endif

  MATRIX Res = m;
  int i, j, k;
  double sum;
  
  for (k = 1; k<=Res.Cols(); k++) {
    // First fixing the diagonal elements
    sum = 0;
    for (j = 1; j<k; j++) {
      sum += Res(k,j)*Res(k,j);
    }
    Res(k,k) = sqrt(Res(k,k)-sum);  

    for (i = k+1; i<=Res.Cols(); i++) {
      sum = 0;
      for (j = 1; j<k; j++) {
	sum += Res(i,j)*Res(k,j);
      }
      Res(i,k) = (Res(i,k)-sum)/(Res(k,k));
    }

  }
  // Sets the upper matrix to 0
  for (k = 1; k<=Res.Cols(); k++) {
    for (i = k+1; i<=Res.Cols(); i++) {
      Res(k,i) = 0;
    }
  }

  return(Res);
}

MATRIX Transpose(const MATRIX &m)
{  
  MATRIX Res(m.Cols(), m.Rows());

  for (int i=1; i<=Res.Rows(); i++)
  {
    for (int j=1; j<=Res.Cols(); j++)
      Res(i,j) = m(j,i);
  }

  return Res;
}

MATRIX SubMatrix(const MATRIX &m, const MATRIX &select, int dimensions)
{
  #ifdef CHECK_BOUNDS
  switch (dimensions) 
  {
  case 1:  assert(m.Rows() == select.Rows()); break;
  case 2:  assert(m.Cols() == select.Rows()); break;
  case 3:  assert(m.Rows() == select.Rows()); assert(m.Cols() == select.Rows()); break;
  } 
  #endif

  int nSize, i, j, row, col;

  MATRIX Res(0,0);

  // Calculate size of submatrix
  nSize = 0;
  for (i=1; i<=select.Rows(); i++)
  {
    if (select(i,1))
      nSize++;    
  }

  switch (dimensions) 
  {
  case 1: Res.Resize(nSize, m.Cols()); 
          row = 0;
          for (i=1; i<=select.Rows(); i++)
          {
            if (select(i,1))
	    {
              row++;
              for (j=1; j<=m.Cols(); j++)
              {
		Res(row,j) = m(i,j);
	      }
	    }
	  }
	  break;
  case 2: Res.Resize(m.Rows(),nSize); 
          col = 0;
          for (i=1; i<=select.Rows(); i++)
          {
            if (select(i,1))
	    {
              col++;
              for (j=1; j<=m.Rows(); j++)
              {
		Res(j,col) = m(i,j);
	      }
	    }
	  }
	  break;
  case 3: Res.Resize(nSize,nSize); 
	  row = 0;
          for (i=1; i<=select.Rows(); i++)
          {
            col = 0;
            if (select(i,1))
	    {
              row++;
              for (j=1; j<=select.Rows(); j++)
              {
		if (select(j,1))
		{
		  col++;
		  Res(row,col) = m(i,j);
		}
	      }
	    }
	  }
	  break;
  default: return Res;
  }

  return Res;
}


MATRIX Sweep(const MATRIX &m, int k)
{
  MATRIX Res = m; // Copies m
  int i, j;

  #ifdef CHECK_BOUNDS
  // Checking square matrix
  assert(m.Cols() == m.Rows());
  #endif

  Res(k,k) = -1/m(k,k);

  for (i = 1; i<=Res.Rows(); i++)
  {
    if (i==k) 
      continue;

    Res(i,k) = -Res(k,k)*m(i,k);
    Res(k,i) = Res(i,k);

  }

  for (i = 1; i<=Res.Rows(); i++)
  {
    if (i==k)
      continue;

    //    Res(i,k) = -Res(k,k)*m(i,k);
    //    Res(k,i) = Res(i,k);

    for (j = i; j<=Res.Cols(); j++)
    {
      if (j==k)
        continue;

      Res(i,j) = m(i,j) - Res(i,k)*m(k,j);
      Res(j,i) = Res(i,j);
    }
  }
  return Res;
}


// NO Pivoting yet
MATRIX GaussJordan(const MATRIX &m)
{
  int i, j, k, n;
  double dummy;
  MATRIX Res = m;

#ifdef CHECK_BOUNDS
  // Checking square matrix
  assert(m.Cols() == m.Rows());
#endif

  n = Res.Cols();
  Res = 0;

  for (i = 1; i <= n ; i++) {
    Res(i,i) = 1;
  }

  for (k = 1; k <= n ; k++) {
    dummy = m(k,k);
    // Scale the k'th row
    for (j = 1; j<=n; j++) {
      Res(k,j) /= dummy;
    }
    for (i = 1; i<= n; i++) {
      if (i == k)
	continue;
      dummy = m(i,k);
      for (j = 1; j<=n; j++) {
	Res(i,j) += -dummy*Res(k,j);
      }
    }
  }
  return(Res);
}


MATRIX Inverse(const MATRIX &m, double *logdet)
{
  // Uses the sweep operator
  // Requires symmetric PD matrix (prevents pivoting)

#ifdef CHECK_BOUNDS
  // Checking square matrix
  assert(m.Cols() == m.Rows());
#endif

  MATRIX Res = m;
  double determinant = 0.0;
  double dummy;

  for (int k = 1; k<=Res.Cols(); k++)
  {
    dummy = Res(k,k);
    determinant += log(fabs(dummy));
    for (int i = 1; i<=Res.Cols(); i++)
    {
      if (i == k)
        continue;
      for (int j = i; j<=Res.Cols(); j++)
      {
        if (j==k)
	  continue;
        Res(i,j) -= Res(i,k)*Res(k,j)/dummy;
        Res(j,i) = Res(i,j);
      }
    }

    for (int i = 1; i<=Res.Cols(); i++)
    {
      if (i == k) 
        continue;
      Res(i,k) /= dummy;
      Res(k,i) = Res(i,k);
    }
    Res(k,k) = -1.0/dummy;
  }

  // If there is a place to hold the log determinant
  if (logdet)
    *logdet = determinant;
  return Res*-1;
}

double LogDeterminant(const MATRIX &m)
{
  m.Print();

  exit(1);
}




// --------------------------------------------------
// Here comes the derived class for vectors
// -------------------------------------------------- 

double ColVector::Sum()
{
  double res = 0;
  for (int i = 1; i<Cols(); i++)
    res += (*this)(i,1);
  return(res);
}

double ColVector::Mean()
{
  double res = 0;
  for (int i = 1; i<Cols(); i++)
    res += (*this)(i,1);
  return(res / Cols());
}

double ColVector::Variance()
{
  double mean = 0, res = 0;
  int i;

  #ifdef CHECK_BOUNDS
  // Checking that at least 1 observation
  assert (Rows()>1);
  #endif

  for (i = 1; i<=Cols(); i++)
    mean += (*this)(i,1);
  mean /= (double)Cols();

  for (i = 1; i<=Cols(); i++)
    res +=  ((*this)(i,1)-mean)*((*this)(i,1)-mean);

  return(res / sqrt(Cols()-1));
}



MATRIX SubMatrix(const MATRIX& b, const MATRIX& a)
{
#ifdef CHECK_BOUNDS
  assert(a.Rows() == b.Rows());
  assert(a.Cols() == 1);
#endif

  int nSize, i, j, ii, jj;
 
  nSize = 0;
  for (i=1; i<=a.Rows(); i++)
  {
    if (a(i,1))
      nSize++;
  }

  MATRIX Res(nSize, nSize);

  ii = 0;
  jj = 0;
  for (i=1; i<=b.Rows(); i++)
  {
    jj=0;
    if (a(i,1))
    {
      ii++;      
      for (j=1; j<=b.Cols(); j++)
      {
        if (a(j,1))
	{
          jj++;	
          Res(ii,jj) = b(i,j);
	}
      }
    }
  }
  return Res;
}


MATRIX Hademard(const MATRIX &m, const MATRIX &n)
{

#ifdef CHECK_BOUNDS
  // Check that matrices conform
  assert(m.Cols() == n.Cols());
  assert(m.Rows() == n.Rows());
#endif

  int i, j;

  MATRIX res(m.Rows(),m.Cols());
  res.m_nMatrixType = m.m_nMatrixType;

  for (i=1; i<=m.Rows(); i++)
  {
    for (j=1; j<=m.Cols(); j++)
      res(i,j) = m(i,j) * n(i,j);
  }

  return res;
}




// XXX laver denne nu. Skal kende størrelsen på familierne, og de skal
// ligge pænt efter hinanden
// Returnerer en liste af matricer
// Bliver ikke pænt slettet
MATRIX* FileReadSymmetricKnownSize(char *fname, int size)
{
  FILE *F;
  int maximum, i, j, nPed, ped;
  float res;
  MATRIX *result = NULL;

  F = fopen(fname, "r");

  if (F != NULL)
  {
    maximum = 0;
    while (fscanf(F, "%5d %5d %9f", &i, &j, &res) != EOF)
    {
      // First figure out the dimensions
      if (i>maximum)
        maximum = i;
    }
    fclose(F);
  }
  else {
    printf("Error: File %s not found\n", fname);
    exit(1);
  }


  F = fopen(fname, "r");

  if (F != NULL && maximum>0)
  {
    //    (*this).Resize(max(nRow, maximum),max(nCol, maximum));

    nPed = maximum/size;

    // Allocates memory
    result = new MATRIX[nPed];

    // Resize all matrices
    for (i = 0; i<nPed; i++) {
      result[i].Resize(size,size);
      result[i] = 0.0;
    }

    // Reads files
    ped = 0;
    while (fscanf(F, "%5d %5d  %9f\n", &i, &j, &res) != EOF)
    {
      if ((i-1) / size != (j-1)/ size) {
	printf("WARNING: Problems in FileReadSymmetricKnownSize: %d %d\n", i, j);
      }

      ped = (i-1) / size;

      result[ped](i - ped*size,j-ped*size) = res;
      result[ped](j - ped*size,i-ped*size) = res;
      
    }
  }
  fclose(F);

  return result;
}




 

//
// This a slow but simple version suitable only for 
// smaller covariance matrices.
//

MATRIX FindGroups(const MATRIX &m)
{
  int n = m.Rows(); 
  int i, j, k, group;
  MATRIX result(n,1), keep(n,1);

  // Initially people are in their own group
  for (i = 1; i<=n ; i++) {
    result(i,1) = i;
    keep(i,1) = 1;
  }

  for (i = 1; i<n ; i++) {
    group = (int) result(i,1);

    for (j = i+1; j<=n ; j++) {      
      if (m(i,j) && (result(j,1) != result(i,1))) {  
	// The pair is in the same group
	// Replace all persons in this group
	for (k = 1; k<=n; k++) {
	  if (result(k,1) == group)
	    result(k,1) = result(j,1);
	}	  
	group = (int) result(j,1);
      }
    }
  }

  // Should normalize the groups
  result = result*-1;
  group = 0;
  for (i = 1; i<=n ; i++) {
    if (result(i,1)<0) {
      group++;
      for (j = i+1; j<=n; j++) {
	if (result(j,1) == result(i,1)) {
	  result(j,1) = group;
	}
      }
      result(i,1) = group;
    }
  }


  return(result);
}


MATRIX CombineMatrices(const MATRIX m[], const int nMat)
{
  int i, currow, curcol, totalrow, totalcol, thiscol, j, k;
  MATRIX res;

  totalrow = 0;
  totalcol = 0;

  for (i = 0; i< nMat; i++) {
    totalrow += m[i].Rows();
    totalcol += m[i].Cols();
  }

  res.Resize(totalrow, totalcol);
  res = 0;

  currow = 0;
  curcol = 0;
  thiscol= 0;
  for (i = 0; i< nMat; i++) {
    for (j = 1; j<=m[i].Rows(); j++) {
      currow++;
      curcol = thiscol;
      for (k = 1; k<=m[i].Cols(); k++) {
	curcol++;
	res(currow, curcol) = m[i](j,k);	
      }      
    }	
    thiscol += m[i].Cols();
  }
  
  return(res);
}



MATRIX AppendMatrices(const MATRIX m[], const int nMat)
{
  int i, currow, curcol, totalrow, totalcol, j, k;
  MATRIX res;

  totalrow = 0;
  totalcol = m[0].Cols();

  for (i = 0; i< nMat; i++) {
    totalrow += m[i].Rows();
    if (m[i].Cols() != totalcol) {
      printf("ERROR: Appending matrices with different number of columns\n");
      exit(1);
    }
  }

  res.Resize(totalrow, totalcol);
  res = 0;

  currow = 0;
  curcol = 0;
  for (i = 0; i< nMat; i++) {
    for (j = 1; j<=m[i].Rows(); j++) {
      currow++;
      curcol = 0;
      for (k = 1; k<=m[i].Cols(); k++) {
	curcol++;
	res(currow, curcol) = m[i](j,k);	
      }      
    }	
  }
  
  return(res);
}
