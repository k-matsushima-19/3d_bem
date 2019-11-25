#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include "gauss.h"

typedef unsigned uint;

typedef unsigned longindex;

typedef double real;

typedef real *preal;

typedef const real *pcreal;

typedef double _Complex field;

typedef field *pfield;

typedef const field *pcfield;

typedef struct _amatrix amatrix;

typedef amatrix *pamatrix;

typedef const amatrix *pcamatrix;

struct _amatrix {
  field *a;

  uint ld;

  uint rows;
  uint cols;

  void *owner;
};

typedef struct _tridiag tridiag;

typedef tridiag *ptridiag;

typedef const tridiag *pctridiag;

struct _tridiag {
  preal d;

  preal l;
  
  preal u;

  uint size;

  ptridiag owner;
};


void zsteqr_(const char *compz, const unsigned *n, double *d, double *e,
    double _Complex *z, const unsigned *ldz, double *work, int *info);

void dstev_(const char *jobz, const unsigned *n, double *d, double *e, double *z,
    const unsigned *ldz, double *work, int *info);


static real ABSSQR(field x) {
  real rx = creal(x);
  real ix = cimag(x);

  return rx * rx + ix * ix;
}

static pamatrix init_amatrix(pamatrix a, uint rows, uint cols)
{
  assert(a != NULL);

  a->a = (rows > 0 && cols > 0 ? (field *) malloc(sizeof(field) * rows * cols) : NULL); // (rows > 0 && cols > 0 ? allocmatrix(rows, cols) : NULL);
  a->ld = rows;
  a->rows = rows;
  a->cols = cols;
  a->owner = NULL;

  return a;
}

static pamatrix new_amatrix(uint rows, uint cols)
{
  pamatrix  a;

  a = (pamatrix) malloc(sizeof(amatrix));

  init_amatrix(a, rows, cols);

  return a;
}

static field getentry_amatrix(pcamatrix a, uint row, uint col) {
  longindex lda = a->ld;

  return a->a[row + lda * col];
}

static void uninit_amatrix(pamatrix a)
{
  if (!a->owner)
    free(a->a);
}

static void del_amatrix(pamatrix a)
{
  uninit_amatrix(a);
  free(a);
}

static ptridiag init_tridiag(ptridiag T, uint size)
{
  T->d = (real *) malloc(sizeof(real)*size); // allocreal(size);
  T->u = (size > 1 ? (real *) malloc(sizeof(real)*(size-1)) : NULL); // (size > 1 ? allocreal(size - 1) : NULL);
  T->l = (size > 1 ? (real *) malloc(sizeof(real)*(size-1)) : NULL); // (size > 1 ? allocreal(size - 1) : NULL);
  T->size = size;
  T->owner = NULL;

  return T;
}

static ptridiag new_tridiag(uint size)
{
  ptridiag  T;

  T = (ptridiag) malloc(sizeof(tridiag));

  init_tridiag(T, size);

  return T;
}

static void identity_amatrix(pamatrix a)
{
  uint      rows = a->rows;
  uint      cols = a->cols;
  longindex lda = a->ld;
  uint      i, j;

  for (j = 0; j < cols; j++) {
    for (i = 0; i < rows; i++) {
      a->a[i + j * lda] = 0.0;
    }
  }

  for (i = 0; i < cols && i < rows; i++) {
    a->a[i + i * lda] = 1.0;
  }
}

static uint muleig_tridiag(ptridiag T, pamatrix Q)
{
  preal     work;
  int       info = 0;

  const uint u_one = 1;

  if (T->size > 1) {
    if (Q) {
      work = (real *) malloc(sizeof(real)*(2 * T->size + 2)); // allocreal(2 * T->size + 2);
      zsteqr_("V", &T->size, T->d, T->l, Q->a, &Q->ld, work, &info);
      assert(info >= 0);
      free(work);
    }
    else {
      dstev_("N", &T->size, T->d, T->l, NULL, &u_one, NULL, &info);
      assert(info >= 0);
    }
  }
  return (info != 0);
}

static uint eig_tridiag(ptridiag T, pamatrix Q)
{
  if (Q)
    identity_amatrix(Q);

  return muleig_tridiag(T, Q);
}

static void uninit_tridiag(ptridiag T)
{
  if (T->owner == NULL) {
    free(T->d);
    if (T->size > 1) {
      free(T->l);
      free(T->u);
    }
    else {
      assert(T->l == NULL);
      assert(T->u == NULL);
    }
  }
}

static void del_tridiag(ptridiag T)
{
  uninit_tridiag(T);

  free(T);
}

void generate_gauss(int m, double *x, double *w)
{
  ptridiag  T;
  pamatrix  Q;
  int i;
  uint      info;

  T = new_tridiag(m);
  for (i = 0; i < m; i++)
    T->d[i] = 0.0;
  for (i = 0; i < m - 1; i++)
    T->l[i] = (i + 1.0) / sqrt((2.0 * i + 1.0) * (2.0 * i + 3.0));

  Q = new_amatrix(m, m);

  info = eig_tridiag(T, Q);
  assert(info == 0);

  for (i = 0; i < m; i++) {
    x[i] = T->d[i];
    w[i] = 2.0 * ABSSQR(getentry_amatrix(Q, 0, i));
  }

  del_amatrix(Q);
  del_tridiag(T);
}
