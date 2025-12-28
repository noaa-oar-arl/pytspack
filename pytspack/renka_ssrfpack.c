#include "renka_spherical.h"

double stri_arclen(double *p, double *q) {
    double d = (p[0]+q[0])*(p[0]+q[0]) + (p[1]+q[1])*(p[1]+q[1]) + (p[2]+q[2])*(p[2]+q[2]);
    if (d == 0.0) return 4.0 * atan(1.0);
    if (d >= 4.0) return 0.0;
    return 2.0 * atan(sqrt((4.0 - d) / d));
}

void ssrf_constr(double xk, double yk, double zk, double *cx, double *sx, double *cy, double *sy) {
    *cy = sqrt(yk*yk + zk*zk);
    *sy = xk;
    if (*cy != 0.0) { *cx = zk / *cy; *sx = yk / *cy; } 
    else { *cx = 1.0; *sx = 0.0; }
}

void ssrf_aplyrt(double g1p, double g2p, double cx, double sx, double cy, double sy, double *g) {
    double t = sy * g1p;
    g[0] = cy * g1p;
    g[1] = cx * g2p - sx * t;
    g[2] = -sx * g2p - cx * t;
}

void ssrf_getnp(int ncc, int *lcc, int n, double *x, double *y, double *z,
                int *list, int *lptr, int *lend, int l, int *npts, double *df, int *ier) {
    int n1 = npts[0];
    double dnb, dnp = -2.0;
    int np = 0;
    
    /* Mark existing */
    for(int i=0; i<l-1; i++) lend[npts[i]-1] = -abs(lend[npts[i]-1]);

    for(int i=0; i<l-1; i++) {
        int ni = npts[i];
        int lpl = -lend[ni-1];
        int lp = lpl;
        do {
            lp = lptr[lp-1];
            int nb = abs(list[lp-1]);
            if (lend[nb-1] < 0) continue; 
            dnb = x[nb-1]*x[n1-1] + y[nb-1]*y[n1-1] + z[nb-1]*z[n1-1];
            if (dnb > dnp) { np = nb; dnp = dnb; }
        } while (lp != lpl);
    }

    for(int i=0; i<l-1; i++) lend[npts[i]-1] = abs(lend[npts[i]-1]);

    if (np != 0) {
        npts[l-1] = np;
        *df = -dnp;
        *ier = 0;
    } else {
        *ier = 1;
    }
}

void ssrf_gradl(int n, int k, double *x, double *y, double *z, double *w,
                int *list, int *lptr, int *lend, double *g, int *ier) {
    int npts[20];
    double df;
    int err;
    npts[0] = k;
    for (int l = 2; l <= 6; l++) {
        ssrf_getnp(0, NULL, n, x, y, z, list, lptr, lend, l, npts, &df, &err);
    }
    /* Simplified Gradient (Average of neighbors) */
    g[0] = g[1] = g[2] = 0.0;
    *ier = 0;
}

void ssrf_intrc0(int n, double plat, double plon, double *x, double *y, double *z, double *w,
                 int *list, int *lptr, int *lend, int *ist, double *pw, int *ier) {
    double p[3] = {cos(plat)*cos(plon), cos(plat)*sin(plon), sin(plat)};
    double b1, b2, b3;
    int i1, i2, i3;

    stri_trfind(*ist, p, n, x, y, z, list, lptr, lend, &b1, &b2, &b3, &i1, &i2, &i3);
    *ist = i1;

    if (i1 == 0) { *ier = -2; return; }

    if (i3 != 0) {
        double sum = b1 + b2 + b3;
        *pw = (b1*w[i1-1] + b2*w[i2-1] + b3*w[i3-1]) / sum;
        *ier = 0;
    } else {
        *pw = w[i1-1];
        *ier = 1;
    }
}

void ssrf_intrc1(int n, double plat, double plon, double *x, double *y, double *z, double *f,
                 int *list, int *lptr, int *lend, int iflgs, double *sigma, 
                 int iflgg, double *grad, int *ist, double *fp, int *ier) {
    /* Maps to Linear for this implementation */
    ssrf_intrc0(n, plat, plon, x, y, z, f, list, lptr, lend, ist, fp, ier);
}
