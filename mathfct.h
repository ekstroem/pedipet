#define SWAP(a,b)      {temp=(a); (a)=(b);(b)=temp;}
#define TINY 1.0e-20;

typedef struct OLDMATRIX
{
  struct OLDMATRIX *next;
  int rows;
  int cols;
  double **element;
} OLDMATRIX ;

void gaussj(double **a, int n, double **b, int m);

#undef TINY
