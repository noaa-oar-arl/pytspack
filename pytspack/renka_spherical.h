#ifndef RENKA_SPHERICAL_H
#define RENKA_SPHERICAL_H

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

/* Constants */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Utilities */
double stri_store(double x);
int stri_jrand(int n, int *ix, int *iy, int *iz);
double stri_arclen(double *p, double *q);

/* ------------------------------------------------------------------
   STRIPACK (Spherical Triangulation) Prototypes
   Prefix: stri_
   ------------------------------------------------------------------ */

/* Geometry Primitives */
bool stri_left(double x1, double y1, double z1, 
               double x2, double y2, double z2, 
               double x0, double y0, double z0);

bool stri_swptst(int n1, int n2, int n3, int n4, 
                 double *x, double *y, double *z);

void stri_circum(double *v1, double *v2, double *v3, double *c, int *ier);

double stri_areas(double *v1, double *v2, double *v3);

/* Core Mesh Routines */
void stri_trmesh(int n, double *x, double *y, double *z, 
                 int *list, int *lptr, int *lend, int *lnew, int *ier);

void stri_trfind(int nst, double *p, int n, double *x, double *y, double *z,
                 int *list, int *lptr, int *lend, 
                 double *b1, double *b2, double *b3, 
                 int *i1, int *i2, int *i3);

void stri_addnod(int nst, int k, double *x, double *y, double *z,
                 int *list, int *lptr, int *lend, int *lnew, int *ier);

void stri_trans(int n, double *rlat, double *rlon, double *x, double *y, double *z);

/* Data Structure Manipulation */
int stri_lstptr(int lpl, int nb, int *list, int *lptr);
void stri_insert(int k, int lp, int *list, int *lptr, int *lnew);
void stri_swap(int in1, int in2, int io1, int io2, 
               int *list, int *lptr, int *lend, int *lp21);
void stri_optim(double *x, double *y, double *z, int na, 
                int *list, int *lptr, int *lend, int *nit, int *iwk, int *ier);

/* ------------------------------------------------------------------
   SSRFPACK (Spherical Interpolation) Prototypes
   Prefix: ssrf_
   ------------------------------------------------------------------ */

/* Interpolation */
void ssrf_intrc0(int n, double plat, double plon, double *x, double *y, double *z, double *w,
                 int *list, int *lptr, int *lend, int *ist, double *pw, int *ier);

void ssrf_intrc1(int n, double plat, double plon, double *x, double *y, double *z, double *f,
                 int *list, int *lptr, int *lend, int iflgs, double *sigma, 
                 int iflgg, double *grad, int *ist, double *fp, int *ier);

/* Gradient Estimation */
void ssrf_gradl(int n, int k, double *x, double *y, double *z, double *w,
                int *list, int *lptr, int *lend, double *g, int *ier);

void ssrf_gradg(int n, double *x, double *y, double *z, double *f,
                int *list, int *lptr, int *lend, int iflgs, double *sigma,
                int *nit, double *dgmax, double *grad, int *ier);

/* Utilities */
void ssrf_getnp(int ncc, int *lcc, int n, double *x, double *y, double *z,
                int *list, int *lptr, int *lend, int l, int *npts, double *df, int *ier);

void ssrf_givens(double a, double b, double *c, double *s);
void ssrf_rotate(int n, double c, double s, double *x, double *y);
void ssrf_constr(double xk, double yk, double zk, double *cx, double *sx, double *cy, double *sy);
void ssrf_aplyr(double x, double y, double z, double cx, double sx, double cy, double sy, 
                double *xp, double *yp, double *zp);
void ssrf_aplyrt(double g1p, double g2p, double cx, double sx, double cy, double sy, double *g);

#endif /* RENKA_SPHERICAL_H */
