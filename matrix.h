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


/*
- CE 25/10 1999
  Should probably make some of the functions inline 
  to save time around the code. Especially the matrix operators.
  This could be done by copying the functions from the .cpp file
  and inserting them here.

*/

#ifndef EoFMATRIX
#define EoFMATRIX

#define CHECK_BOUNDS

enum {
  MAT_GENERAL,
  MAT_SYMMETRIC,
  MAXMAT,
};

class MATRIX {

private:
  int m_nRows;                // Rows
  int m_nCols;                // Cols
  //  long m_datalength           // Length of data
  double *data;               // Holds the data

public:
  MATRIX();                   // Constructor
  MATRIX(int rows, int cols);

  virtual ~MATRIX();                  // Destructor

  inline virtual int Rows(void) const; // Number of Rows
  inline virtual int Cols(void) const; // Number of Cols

  MATRIX& operator=(double s);         // Assignment operator
  MATRIX& operator=(const MATRIX &s);  // Assignment operator
  MATRIX(const MATRIX &m);	       // Copy operator

  inline double& operator()(int row, int col);       // Referencing
  inline double& operator()(int row, int col) const; // Referencing

  //  MATRIX& operator()(int row, int col) const; // Create submatrix

  MATRIX& Resize(const int rows, const int cols);    // Clears and resizes matrix
  MATRIX& Insert(const int rows, const int cols, const MATRIX& a);    // Inserts A from row, col

  MATRIX  Subset(const MATRIX& a);                  // Picks subset from indicator vector
  MATRIX  Subset(const MATRIX& a, const MATRIX& b); // Picks subset from indicator vector


  // Simple functions
  MATRIX operator*(double s);        // Scale
  MATRIX operator+(double s);        // Add constant
  MATRIX operator-(double s);        // Subtract constant

  MATRIX operator+=(const MATRIX &m);  // Assignment operator
  MATRIX operator*=(const double a);   // Matrix scaling

  void Print(void) const;           // Prints the matrix
  int FPrint(char *fname) const;    // Prints the matrix to a file

  int FileWriteSymmetric(char *fname);    // Writes the matrix to a file
  int FileReadSymmetric(char *fname, int nCol, int nRow);     // Reads the matrix to a file

  double Sum(void);             // Sum of all elements in the matrix
  double SSD(const double mean);// Sum of Squared Deviations from the mean
  double Min(void);             // Minimum of all elements
  double Max(void);             // Maximum of all elements

public:
  int m_nMatrixType;
  
};




// Column vector class has some additional functions
class ColVector: public MATRIX
{

public:
  double Sum();
  double Mean();
  double Variance();

};



//-----------------------------
// Definition of math functions / relations
//-----------------------------

// NOTE: These should all be inline to prevent the 
// destructor of the result to be called

// Addition
MATRIX operator+(const MATRIX& m, const MATRIX& n);
//MATRIX operator+=(MATRIX& m, const MATRIX& n);

// Subtraction
MATRIX operator-(const MATRIX& m, const MATRIX& n);

// Multiplication
MATRIX operator*(const MATRIX& m, const MATRIX& n);

// Get a submatrix
// select is a n*1 vector, that determines which rows/columns are
// selected
MATRIX SubMatrix(const MATRIX &m, const MATRIX &select, int dimensions);


MATRIX CombineMatrices(const MATRIX m[], const int nMat);
MATRIX AppendMatrices(const MATRIX m[], const int nMat);

// Sweep operator
MATRIX Sweep(const MATRIX &m, int k);

// Inversion
MATRIX Inverse(const MATRIX &m, double *logdet = NULL);
MATRIX GaussJordan(const MATRIX &m);
double LogDeterminant(const MATRIX &m);

// Trace
double Trace(const MATRIX& m);

// Cholesky LU-decomposition
// Should be symmetric and PD
MATRIX Cholesky(const MATRIX &m);


// Transpose
MATRIX Transpose(const MATRIX &m);

MATRIX Hademard(const MATRIX &m, const MATRIX &n);

MATRIX SubMatrix(const MATRIX& b, const MATRIX& a);  // a is a (column-)vector. 
                                // Selects the col/rows of b when a <>0



// Under development

// Returns a vector of correlated groups from a PD matrix
MATRIX FindGroups(const MATRIX &m);

MATRIX *FileReadSymmetricKnownSize(char *fname, int size);     // Reads the matrix to a file

/***************************************************************
 *
 * SYMMATRIX
 *
 *
 *
 ***************************************************************/

// Column vector class has some additional functions
class SYMMATRIX: public MATRIX
{

};



#endif


