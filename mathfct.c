#include <math.h>
#include "pedipet.h"
#include "nrutil.h"

#define SWAP(a,b)      {temp=(a); (a)=(b);(b)=temp;}
#define TINY 1.0e-20;

/* 

  Gauss Jordan elimination 

*/

void gaussj(double **a, int n, double **b, int m)
{
  int *indxc, *indxr, *ipiv;
  int i, icol, irow,j, k, l, ll;
  double big, dum, pivinv, temp;

  icol = 0;
  irow = 0;

  indxc = ivector(1,n);
  indxr = ivector(1,n);
  ipiv  = ivector(1,n);

  for (j=1; j<=n; j++)
    ipiv[j] = 0;
  for (i=1; i<=n; i++)
  {
    big = 0.0;
    for (j=1; j<=n; j++)
    {
      if (ipiv[j] != 1)
      {
        for (k=1; k<=n;k++)
	{
          if (ipiv[k] == 0)
	  {
            if (absolut(a[j][k])>=big)
	    {
              big = absolut(a[j][k]);
              irow = j;
              icol = k;
	    }
	  }
          else if (ipiv[k]>1)
            nrerror("Gauss-Jordan: Singular matrix (1)");
        }
      }
    }
    ++(ipiv[icol]);

    if (irow != icol)
    {
      for (l=1; l<=n; l++) 
        SWAP(a[irow][l],a[icol][l]); 
      for (l=1; l<=m; l++) 
        SWAP(b[irow][l],b[icol][l]); 
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if (a[icol][icol] == 0.0)
      nrerror("Gauss Jordan: Singular matrix (2)");
    pivinv = 1.0/a[icol][icol];
    a[icol][icol]= 1.0;
    for (l=1; l<=n; l++)
      a[icol][l] *= pivinv;
    for (l=1; l<=m; l++)   
      b[icol][l] *= pivinv; 
    for (ll=1; ll<=n; ll++)
    {
      if (ll != icol)
      {
        dum = a[ll][icol];
        a[ll][icol] = 0.0;
        for (l=1; l<=n; l++)    
          a[ll][l] -= a[icol][l]*dum;
        for (l=1; l<=m; l++)    
          b[ll][l] -= b[icol][l]*dum;

      }
    }

  }

  for (l=n; l>=1; l--)
  {
    if (indxr[l] != indxc[l])
      for (k=1; k<=n; k++)
        SWAP(a[k][indxr[l]], a[k][indxc[l]]);
  }
  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);
}

void ludcmp(double **a, int n, int *indx, double *d)
{
        int i,imax,j,k;
        double big,dum,sum,temp;
        double *vv;

	imax = 0;
 
        vv=vector(1,n);
        *d=1.0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if ((temp=fabs(a[i][j])) > big) big=temp;
                if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
                vv[i]=1.0/big;
        }
        for (j=1;j<=n;j++) {
                for (i=1;i<j;i++) {
                        sum=a[i][j];
                        for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                }
                big=0.0;
                for (i=j;i<=n;i++) {
                        sum=a[i][j];
                        for (k=1;k<j;k++)
                                sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                        if ( (dum=vv[i]*fabs(sum)) >= big) {
                                big=dum;
                                imax=i;
                        }
                }
                if (j != imax) {
                        for (k=1;k<=n;k++) {
                                dum=a[imax][k];
                                a[imax][k]=a[j][k];
                                a[j][k]=dum;
                        }
                        *d = -(*d);
                        vv[imax]=vv[j];
                }
                indx[j]=imax;
                if (a[j][j] == 0.0) a[j][j]=TINY;
                if (j != n) {
                        dum=1.0/(a[j][j]);
                        for (i=j+1;i<=n;i++) a[i][j] *= dum;
                }
        }
        free_vector(vv,1,n);
}

double Determinant(double **a, int n)
{
  int j, *indx;
  double d;

  indx = ivector(1,n);

  ludcmp(a, n, indx, &d);
  for (j=1; j<=n; j++)
    d *= a[j][j];

  free_ivector(indx,1,n);
  return d;
}


// Returns log of the determinant, assuming it is positive
double LogDeterminant(double **a, int n)
{
  int j, *indx;
  double d;

  indx = ivector(1,n);

  ludcmp(a, n, indx, &d);
  d = 0.0;
  for (j=1; j<=n; j++)
    d += log(fabs(a[j][j]));

  free_ivector(indx,1,n);
  return d;
}



OLDMATRIX *MtxNew(int n, int m)
{
  OLDMATRIX *newmat;
  
  newmat = cmalloc(sizeof(OLDMATRIX));
  memset (newmat, 0, sizeof(OLDMATRIX));

  newmat->rows = n;
  newmat->cols = m;
  newmat->element = matrix(1, n, 1, m);

  return(newmat);
}

void MtxDel(OLDMATRIX *a)
{
  free_matrix(a->element, 1, a->rows,1,  a->cols);
  free(a);
//   removelist(listofmatrices, a);
}

void MtxFreeList(OLDMATRIX *a)
{
  OLDMATRIX *b, *c;

  for (b = a; b; )
  {
    c = b->next;
    MtxDel(b);
    b = c;
  }
}


// OLDMATRIX multiplication
// Maxtrix mata (dim n1*m1), matb (dim n2*m2)
// Return result in OLDMATRIX res
void MtxMulti(OLDMATRIX *a, OLDMATRIX *b, OLDMATRIX *res)
{
  int i, j, k;
  double tempres, **tmpmatr;

  if (a->cols!=b->rows)
  {
    printf("Matrices does not conform\n");
    return;
  }

  tmpmatr = matrix(1, a->rows, 1, b->cols);

  for (i=1; i<=a->rows; i++)
  {
    for (j=1; j<=b->cols; j++)
    {
      tempres = 0.0;
      for (k=1; k<=a->cols; k++)
      {
        tempres +=a->element[i][k]*b->element[k][j];
      }
      tmpmatr[i][j] = tempres;
    }
  }
  for (i=1; i<=a->rows; i++)
    for (j=1; j<=b->cols; j++)
      res->element[i][j] = tmpmatr[i][j];

  free_matrix(tmpmatr, 1, a->rows, 1, b->cols);
}


// OLDMATRIX dot
// Maxtrix mata (dim n1*m1), matb (dim n2*m2)
// Return result in OLDMATRIX res
double MtxDot(OLDMATRIX *a, OLDMATRIX *b)
{
  int i;
  double tempres;

  if (a->cols!=b->rows)
  {
    printf("Matrices does not conform\n");
    exit(1);
  }

  if (a->rows != 1 || b->cols != 1)
  {
    printf("ERROR: Can only input vectors\n");
    exit(1);
  }

  tempres = 0.0;

  for (i=1; i<=a->cols; i++)
  {
    tempres +=a->element[1][i]*b->element[i][1];
  }
  return tempres;
}


void MtxCopy(OLDMATRIX *a, OLDMATRIX *b)
{
  int i, j;

  for (i=1; i<=a->rows; i++)
    for (j=1; j<=a->cols; j++)
      b->element[i][j] = a->element[i][j];
}

void MtxAdd(OLDMATRIX *a,  OLDMATRIX *b, OLDMATRIX *c)
{
  int i, j;
  double **tmpmatr;

  if (a->cols!=b->cols || a->rows!=b->rows)
  {
    printf("Matrices do not conform\n");
    return;
  }

  tmpmatr = matrix(1, a->rows, 1, b->cols);

  for (i=1; i<=a->rows; i++)
  {
    for (j=1; j<=a->cols; j++)
    {
      tmpmatr[i][j] = a->element[i][j] + b->element[i][j];
    }
  }

  for (i=1; i<=a->rows; i++)
    for (j=1; j<=b->cols; j++)
      c->element[i][j] = tmpmatr[i][j];

  free_matrix(tmpmatr, 1, a->rows, 1, b->cols);

}

void MtxScale(OLDMATRIX *a,  double b)
{
  int i, j;

  for (i=1; i<=a->rows; i++)
  {
    for (j=1; j<=a->cols; j++)
    {
      a->element[i][j] = a->element[i][j]*b;
    }
  }
}

// 
double MtxSum(OLDMATRIX *a)
{
  int i;
  double res;

  res = 0.0;

  if (a->rows != 1 && a->cols !=1)
  {
    printf("Matrix must be a vector to sum\n");
    return 0;
  }

  if (a->rows>1)
  {
    for (i=1; i<=a->rows; i++)
      res += a->element[i][1];
  }
  else if (a->cols>1)
  {
    for (i=1; i<=a->cols; i++)
      res += a->element[1][i];
  }
  return res;
}

int MtxPositive(OLDMATRIX *a)
{
  int i;
  int res;

  res = 0;

  if (a->rows != 1 && a->cols !=1)
  {
    printf("Matrix must be a vector to sum\n");
    return 0;
  }

  if (a->rows>1)
  {
    for (i=1; i<=a->rows; i++)
      if (a->element[i][1])
        res ++;
  }
  else if (a->cols>1)
  {
    for (i=1; i<=a->cols; i++)
      if (a->element[1][i])
        res ++;
  }
  return res;
}




// Keeps the rows/cols in a for which the number is in vect,
// newdim is the number of elements in vect
// c is a n * 1 vector of 0/1 (0 = keep, 1 = remove)
// If input is a column or a row vector the reduce it
OLDMATRIX *MtxReduce(OLDMATRIX *a, OLDMATRIX *c)
{
  OLDMATRIX *b;
  int i, j, x, newdim, k, l;

  

  newdim = c->rows - (int) MtxPositive(c);
  x = 0;
  if (a->rows==1)
    x = 1;
  if (a->cols==1)
    x = 2;
  if (a->cols ==a->rows)
  x = 3;

  if (x==0)
  {
    printf("Matrix must be square or a vector to reduce\n");
    exit(1);
  }

  switch (x)
  {
    case 1:   b = MtxNew(1, newdim); break;
    case 2:   b = MtxNew(newdim, 1); break;
    case 3:   b = MtxNew(newdim, newdim); break;
  default: b = NULL; break;
  }

  k = 1;
  l = 1;
  for (i=1; i<=c->rows; i++)
  {
    l = 1;
    if (!c->element[i][1])
    {
      for (j=1; j<=c->rows; j++)
      {
        if (!c->element[j][1])
        {
          switch (x)
          { 
            case 1: b->element[1][l] = a->element[1][j]; break;
            case 2: b->element[k][1] = a->element[i][1]; break;
            case 3: b->element[k][l] = a->element[i][j]; break;
          }
          l++;
	}
      }
      k++;
    }
  }
  return b;
}



void MtxChangeSign(OLDMATRIX *a)
{
  int i, j;

  for (i=1; i<=a->rows; i++)
    for (j=1; j<=a->cols; j++)
      a->element[i][j] = -a->element[i][j];
}


void MtxInver(OLDMATRIX *a)
{
  double **bsolve;

  if (a->rows != a->cols)
  {
    printf("Matrix is not square. Cannot invert\n"); 
    return;
  }

  bsolve = matrix(1, a->rows, 1, 1); 
  gaussj(a->element, a->rows, bsolve, 1);
  free_matrix(bsolve, 1, a->rows, 1, 1);
}

void MtxTrans(OLDMATRIX *mat, OLDMATRIX *res)
{
  int i, j;

  for (i=1; i<=mat->rows; i++)
  {
    for (j=1; j<=mat->cols; j++)
    {
      res->element[j][i] = mat->element[i][j];
    }
  }
}

void MtxPrint(OLDMATRIX *a)
{
  int i, j;

  for (i=1; i<=a->rows; i++)
  {
    for (j=1; j<=a->cols; j++)
    {
      printf("%-5.2f ", a->element[i][j]);
    }
    printf("\n");
  }
}

void MtxFPrint(char *filename , OLDMATRIX *a)
{
  int i, j;

  cfopen (filename,"w");

  for (i=1; i<=a->rows; i++) {
    for (j=1; j<=a->cols; j++)  {
      fprintf(F, "%-3.2f ", a->element[i][j]);
    }
    fprintf(F, "\n");
  }
  fflush(F);
  fclose(F);
}

void MtxFSymmetricPrint(char *filename , OLDMATRIX *a)
{
  int i, j;

  cfopen (filename,"w");

  for (i=1; i<=a->rows; i++) {
    for (j=1; j<=i; j++) {
      if (a->element[i][j])
	fprintf(F, "%5d %5d  %9.7f\n", i, j, a->element[i][j]);
    }
  }
  fflush(F);
  fclose(F);
}

OLDMATRIX *MtxFRead(char *filename)
{
  int i, j, n;
  OLDMATRIX *res;
  double value;

  printf("start\n");

  cfopen (filename,"r");

  getbuf();

  n = 0;

  if  (strlen(igetstr(buf)))
  {
    n = 1;
    while( strlen(getstr()) )
      n++;
  }

  printf("Size is %d\n", n);

  res = MtxNew(n,n);

  value = atof(igetstr(buf));
  for (i=1; i<=n; i++)
  {
    for (j=1; j<=n; j++)
    {
      //      printf("UGUGU: %d %d\n", i, j);
      res->element[i][j] = value;
      value = atof(getstr());   
    }
    getbuf();
    value = atof(igetstr(buf));
  }
  return(res);

}


double MtxTrace(OLDMATRIX *a)
{
  double res;
  int i;

  if (a->rows != a->cols)
  {
    printf("ERROR: Can only take trace on square matrices\n");
    exit(1);
  }

  res = 0.0;
  for (i=1; i<=a->rows; i++)
    res += a->element[i][i];
  return res;
}

double MtxLogDeterminant(OLDMATRIX *a)
{
  OLDMATRIX *tmpmat;
  int j, *indx;
  double d;

  tmpmat = MtxNew(a->rows, a->cols);
  MtxCopy(a, tmpmat);

  indx = ivector(1,a->rows);

  ludcmp(tmpmat->element, tmpmat->rows, indx, &d);
  d = 0.0;
  for (j=1; j<=a->rows; j++)
    d += log(fabs(tmpmat->element[j][j]));

  free_ivector(indx,1,a->rows);
  MtxDel(tmpmat);
  return d;
}

//
//  Various routines
//

// Calculates the kinship coeff OLDMATRIX
OLDMATRIX *MakeKinshipMatrix(void)
{
  individual *ind, *ind2;
  OLDMATRIX *res;
  int numpers, i, j, done, code;

  numpers = listlen(individuals);
  // Sets tmpstr to 0 for all individual
  // 0 Means not set, 1 means set
  i = 1;
  forind
  {
    strcpy(ind->tmpstr1, "0");
    sprintf(ind->tmpstr2, "%d",i);
    i++;
  }

  // Should figure out the right order. A persons parents should
  // always appear first
  //
  // done is a boolean - are we finished?
  // i is the number to be placed
  // code is a booelan - did we place someone this round?
  done = 0;
  i = 1;
  while (!done)
  {
    // Noone has been placed
    code = 0;
    forind
    {
      // If founder and not previously placed
      if (founder(ind))
      {
        if (atoip(ind->tmpstr1)==0)
        {
          sprintf(ind->tmpstr1, "%d",i);
          i++;
          code = 1;
	}
      }
      else if (atoip(ind->father->tmpstr1) && atoip(ind->mother->tmpstr1) && (atoip(ind->tmpstr1)==0))
      {
        // Non founder with both parents in data already
        sprintf(ind->tmpstr1, "%d",i);
        i++;
        code = 1;
      }
    }
    if (i==numpers+1)
    {
      done = 1;
      code = 1;
    }
    // If noone has been placed this turn
    if (!code)
    {
      printf("Error: Could not create kinship matrix.\nConsistency error in dataset\n");
      return 0;
    }
  }

  // Initialize
  res = MtxNew(numpers, numpers);

  for (i = 1; i<= numpers; i++)
  {
    forind
    {
      if (i==atoi(ind->tmpstr1))
        break;
    }
    if (founder(ind))
      res->element[atoi(ind->tmpstr2)][atoi(ind->tmpstr2)] = 0.5;
    else
    {
      res->element[atoi(ind->tmpstr2)][atoi(ind->tmpstr2)] = 0.5 + 0.5*res->element[atoi(ind->father->tmpstr2)][atoi(ind->mother->tmpstr2)];
      // Should now fix the rest
      for (ind2 = individuals; ind2; ind2 = ind2->next)
      {
        if (atoi(ind2->tmpstr1)<i)
        {
          j = atoi(ind2->tmpstr2);
          res->element[atoi(ind->tmpstr2)][j] = 0.5*res->element[j][atoi(ind->father->tmpstr2)]+0.5*res->element[j][atoi(ind->mother->tmpstr2)];
          res->element[atoi(ind2->tmpstr2)][atoi(ind->tmpstr2)] = res->element[atoi(ind->tmpstr2)][atoi(ind2->tmpstr2)];
        }
      }
    }
  }

  return res;
}

// Only works for non-inbred relatives
// Needs the kinship OLDMATRIX as input
OLDMATRIX *MakeDelta7Matrix(OLDMATRIX *kinship)
{
  individual *ind, *ind2;
  OLDMATRIX *res;
  int i, j, n;

  res = MtxNew(kinship->rows, kinship->cols);
  n = kinship->cols;

  i = 1;
  forind
  {
    sprintf(ind->tmpstr1, "%d", i);
    i++;
  }

  i = 1;
  forind
  {
    res->element[i][i] = 1;
    if (i==n)
      break;

    j = i+1;
    for (ind2 = ind->next; ind2; ind2 = ind2->next)
    {
      if (ind->father && ind->mother && ind2->father && ind2->mother)
      {
        res->element[i][j] = kinship->element[atoi(ind->father->tmpstr1)][atoi(ind2->father->tmpstr1)]*
                             kinship->element[atoi(ind->mother->tmpstr1)][atoi(ind2->mother->tmpstr1)] +
                             kinship->element[atoi(ind->father->tmpstr1)][atoi(ind2->mother->tmpstr1)]*
                             kinship->element[atoi(ind->mother->tmpstr1)][atoi(ind2->father->tmpstr1)];
        res->element[j][i] = res->element[i][j];
      }
      j++;
    }
    i++;
  }

  return res;
}



#undef TINY
