#ifndef RENKA_H
#define RENKA_H

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

/* Constants */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------
   TRIPACK (Planar Triangulation) 
   ------------------------------------------------------------------ */
void tri_trmesh(int n, double *x, double *y, int *list, int *lptr, int *lend, int *lnew, int *ier);
void tri_trfind(int nst, double px, double py, int n, double *x, double *y, 
                int *list, int *lptr, int *lend, int *i1, int *i2, int *i3);
void tri_getnp(int ncc, int *lcc, int n, double *x, double *y, int *list, int *lptr, int *lend, 
               int l, int *npts, double *dist, int *ier);

/* Helper geometric/list functions */
int tri_lstptr(int lpl, int nb, int *list, int *lptr);
void tri_addnod(int k, double xk, double yk, int ist, int ncc, int *lcc, int n, 
                double *x, double *y, int *list, int *lptr, int *lend, int *lnew, int *ier);

/* ------------------------------------------------------------------
   SRFPACK (Planar Interpolation)
   ------------------------------------------------------------------ */
void srf_intrc0(double px, double py, int ncc, int *lcc, int n, double *x, double *y, double *z, 
                int *list, int *lptr, int *lend, int *ist, double *pz, int *ier);
void srf_intrc1(double px, double py, int ncc, int *lcc, int n, double *x, double *y, double *z, 
                int *list, int *lptr, int *lend, int iflgs, double *sigma, double *grad, 
                bool dflag, int *ist, double *pz, double *pzx, double *pzy, int *ier);
void srf_gradl(int k, int ncc, int *lcc, int n, double *x, double *y, double *z, 
               int *list, int *lptr, int *lend, double *dx, double *dy, int *ier);

/* ------------------------------------------------------------------
   STRIPACK (Spherical Triangulation)
   ------------------------------------------------------------------ */
void stri_trmesh(int n, double *x, double *y, double *z, 
                 int *list, int *lptr, int *lend, int *lnew, int *ier);
void stri_trfind(int nst, double *p, int n, double *x, double *y, double *z,
                 int *list, int *lptr, int *lend, 
                 double *b1, double *b2, double *b3, 
                 int *i1, int *i2, int *i3);
void stri_addnod(int nst, int k, double *x, double *y, double *z,
                 int *list, int *lptr, int *lend, int *lnew, int *ier);
void stri_trans(int n, double *rlat, double *rlon, double *x, double *y, double *z);

/* ------------------------------------------------------------------
   SSRFPACK (Spherical Interpolation)
   ------------------------------------------------------------------ */
void ssrf_intrc0(int n, double plat, double plon, double *x, double *y, double *z, double *w,
                 int *list, int *lptr, int *lend, int *ist, double *pw, int *ier);
void ssrf_intrc1(int n, double plat, double plon, double *x, double *y, double *z, double *f,
                 int *list, int *lptr, int *lend, int iflgs, double *sigma, 
                 int iflgg, double *grad, int *ist, double *fp, int *ier);
void ssrf_gradl(int n, int k, double *x, double *y, double *z, double *w,
                int *list, int *lptr, int *lend, double *g, int *ier);

#endif
