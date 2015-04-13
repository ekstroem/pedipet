#ifndef NRUTIL
#define NRUTIL

void nrerror(char error_text[]);
/* standard error handler */

double *vector(int nl,int nh);
/* allocate a double vector with subscript range v[nl..nh] */

int *ivector(int nl,int nh);
/* allocate an int vector with subscript range v[nl..nh] */

double *dvector(int nl,int nh);
/* allocate a double vector with subscript range v[nl..nh] */

double **matrix(int nrl, int nrh, int ncl, int nch);
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */

double **dmatrix(int nrl,int nrh,int ncl,int nch);

int **imatrix(int nrl, int nrh, int ncl, int nch);
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */

double **submatrix(double **a, int oldrl, int oldrh, int oldcl,
                  int oldch, int newrl, int newcl);
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */

void free_vector(double *v, int nl, int nh);
/* free a double vector allocated with vector() */

void free_ivector(int *v, int nl, int nh);
/* free an int vector allocated with ivector() */

void free_dvector(double *v, int nl, int nh);
/* free a double vector allocated with dvector() */

void free_matrix(double **m,int nrl, int nrh, int ncl, int nch);
/* free a double matrix allocated by matrix() */

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch); 
/* free a double matrix allocated by dmatrix() */

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);
/* free an int matrix allocated by imatrix() */

void free_submatrix(double **b, int nrl, int nrh, int ncl, int nch);
/* free a submatrix allocated by submatrix() */

#endif
