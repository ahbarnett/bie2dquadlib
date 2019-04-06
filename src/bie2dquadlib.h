// minimal definition of interface to bie2dquadlib

#ifndef BIE2DQUADLIB_H
#define BIE2DQUADLIB_H

#include <complex.h>        // C99 type


// Precision-independent real and complex types for interfacing...
#ifdef SINGLE
  typedef float FLT;
  typedef float complex CPX;
  #define EMACH (float)6.0e-8
#else
  typedef double FLT;
  typedef double complex CPX;
  #define EMACH (float)1.1e-16
#endif

typedef struct {  // all of precomputed stuff to allow fast output of A_{I,J}
  int N;          // how many nodes
  FLT *x;         // 2-by-N node locations (Cartesian)
  FLT *nx;        // 2-by-N unit normals at nodes
  FLT *w;         // N vector of "speed weights", ie parameter wgt times speed,
                  // so that line integral uses nodes x and weights w.

  void (*kerfunc)();  // kernel function ptr

  int p;          // number of nodes per panel (same for all)
  int *panind;    // N vector of panel indices associated to the nodes
  int Np;         // number of panels

  // CRS sparse format, panel-wise...
  int *rptr;      // length-Np row ptr array (all rows exist)
  int *cind;      // nonzero column panel indices in each row
  FLT *Aent;      // real kernel p*p block entries
  // how allow CPX Aent vals?
} quadr;


// precomputation routines - return ptr to a quadr struct:
//quadr* quadr_panel_laplace();
//quadr* quadr_panel_helmholtz();
// etc

void getAblk(int *I, int nI, int* J, int nJ, quadr* Q);

void getAblk_proxy(FLT *pxy, int npxy, int* J, int nJ, quadr* Q);

void quadr_free(quadr* Q);


#endif  // BIE2DQUADLIB_H
