// #include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


#include "nrutil.h"

void nrerror(char error_text[])
/* standard error handler */
{
  //        void _exit();

        fprintf(stderr,"run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        _exit(1);
}

double *vector(int nl,int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double))-nl;
        if (!v) nrerror("allocation failure in vector()");
        return v;
}

int *ivector(int nl,int nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int))-nl;
        if (!v) nrerror("allocation failure in ivector()");
        return v;
}

double *dvector(int nl,int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double))-nl;
        if (!v) nrerror("allocation failure in dvector()");
        return v;
}

double **matrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        int i;
        double **m;

        /* allocate pointers to rows */
        m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*))-nrl;
        if (!m) nrerror("allocation failure 1 in matrix()");

        /* allocate rows and set pointers to them */
        for(i=nrl;i<=nrh;i++) {
                m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double))-ncl;
                if (!m[i]) nrerror("allocation failure 2 in matrix()");
                memset(m[i]+ncl, 0, (nch-ncl+1)*sizeof(double));
        }
        /* return pointer to array of pointers to rows */
        return m;
}

double **dmatrix(int nrl,int nrh,int ncl,int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        int i;
        double **m;

        /* allocate pointers to rows */
        m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*))-nrl;
        if (!m) nrerror("allocation failure 1 in dmatrix()");

        /* allocate rows and set pointers to them */
        for(i=nrl;i<=nrh;i++) {
                m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double))-ncl;
                if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
        }
        /* return pointer to array of pointers to rows */
        return m;
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        int i,**m;

        /* allocate pointers to rows */
        m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*))-nrl;
        if (!m) nrerror("allocation failure 1 in imatrix()");

        /* allocate rows and set pointers to them */
        for(i=nrl;i<=nrh;i++) {
                m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int))-ncl;
                if (!m[i]) nrerror("allocation failure 2 in imatrix()");
        }
        /* return pointer to array of pointers to rows */
        return m;
}

double **submatrix(double **a, int oldrl, int oldrh, int oldcl,
                  int oldch, int newrl, int newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
        int i,j;
        double **m;

        /* allocate array of pointers to rows */
        m=(double **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(double*))-newrl;
        if (!m) nrerror("allocation failure in submatrix()");

        /* set pointers to rows */
        for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i];

        /* return pointer to array of pointers to rows */
        return m;
}

void free_vector(double *v, int nl, int nh)
/* free a double vector allocated with vector() */
{
        free((char*) (v+nl));
}

void free_ivector(int *v, int nl, int nh)
/* free an int vector allocated with ivector() */
{
        free((char*) (v+nl));
}

void free_dvector(double *v, int nl, int nh)
/* free a double vector allocated with dvector() */
{
        free((char*) (v+nl));
}

void free_matrix(double **m,int nrl,int nrh,int ncl,int nch)
/* free a double matrix allocated by matrix() */
{
        int i;

        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
}

void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
/* free a double matrix allocated by dmatrix() */
{
        int i;

        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
}

void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
/* free an int matrix allocated by imatrix() */
{
        int i;

        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
}

void free_submatrix(double **b, int nrl, int nrh, int ncl, int nch)
/* free a submatrix allocated by submatrix() */
{
        free((char*) (b+nrl));
}

