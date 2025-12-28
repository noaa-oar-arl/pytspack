#include "renka_spherical.h"

bool stri_left(double x1, double y1, double z1, double x2, double y2, double z2, double x0, double y0, double z0) {
    return (x0*(y1*z2 - y2*z1) - y0*(x1*z2 - x2*z1) + z0*(x1*y2 - x2*y1) >= 0.0);
}

void stri_trans(int n, double *rlat, double *rlon, double *x, double *y, double *z) {
    for (int i=0; i<n; i++) {
        x[i] = cos(rlat[i]) * cos(rlon[i]);
        y[i] = cos(rlat[i]) * sin(rlon[i]);
        z[i] = sin(rlat[i]);
    }
}

int stri_lstptr(int lpl, int nb, int *list, int *lptr) {
    int lp = lptr[lpl-1];
    while (1) {
        if (abs(list[lp-1]) == nb) return lp;
        lp = lptr[lp-1];
        if (lp == lpl) break;
    }
    return lp;
}

void stri_insert(int k, int lp, int *list, int *lptr, int *lnew) {
    int lsav = lptr[lp-1];
    lptr[lp-1] = *lnew;
    list[*lnew-1] = k;
    lptr[*lnew-1] = lsav;
    (*lnew)++;
}

void stri_swap(int in1, int in2, int io1, int io2, int *list, int *lptr, int *lend, int *lp21) {
    int lp, lph, lpsav;
    lp = stri_lstptr(lend[in1-1], in2, list, lptr);
    if (abs(list[lp-1]) == in2) { *lp21 = 0; return; }

    lp = stri_lstptr(lend[io1-1], in2, list, lptr);
    lph = lptr[lp-1];
    lptr[lp-1] = lptr[lph-1];
    if (lend[io1-1] == lph) lend[io1-1] = lp;

    lp = stri_lstptr(lend[in1-1], io1, list, lptr);
    lpsav = lptr[lp-1];
    lptr[lp-1] = lph;
    list[lph-1] = in2;
    lptr[lph-1] = lpsav;

    lp = stri_lstptr(lend[io2-1], in1, list, lptr);
    lph = lptr[lp-1];
    lptr[lp-1] = lptr[lph-1];
    if (lend[io2-1] == lph) lend[io2-1] = lp;

    lp = stri_lstptr(lend[in2-1], io2, list, lptr);
    lpsav = lptr[lp-1];
    lptr[lp-1] = lph;
    list[lph-1] = in1;
    lptr[lph-1] = lpsav;
    
    *lp21 = lph;
}

bool stri_swptst(int n1, int n2, int n3, int n4, double *x, double *y, double *z) {
    double dx1 = x[n1-1]-x[n4-1], dy1 = y[n1-1]-y[n4-1], dz1 = z[n1-1]-z[n4-1];
    double dx2 = x[n2-1]-x[n4-1], dy2 = y[n2-1]-y[n4-1], dz2 = z[n2-1]-z[n4-1];
    double dx3 = x[n3-1]-x[n4-1], dy3 = y[n3-1]-y[n4-1], dz3 = z[n3-1]-z[n4-1];
    return (dx3*(dy2*dz1 - dy1*dz2) - dy3*(dx2*dz1 - dx1*dz2) + dz3*(dx2*dy1 - dx1*dy2) > 0.0);
}

void stri_trfind(int nst, double *p, int n, double *x, double *y, double *z,
                 int *list, int *lptr, int *lend, 
                 double *b1, double *b2, double *b3, 
                 int *i1, int *i2, int *i3) {
    int n0 = nst;
    if (n0 < 1 || n0 > n) n0 = 1;
    int nf, nl, lp, n1, n2;

    while (1) {
        int lpl = lend[n0-1];
        nl = list[lpl-1];
        lp = lptr[lpl-1];
        nf = list[lp-1];
        n1 = nf;
        
        if (stri_left(x[n0-1], y[n0-1], z[n0-1], x[n1-1], y[n1-1], z[n1-1], p[0], p[1], p[2])) {
            while (1) {
                lp = lptr[lp-1];
                n2 = abs(list[lp-1]);
                if (!stri_left(x[n0-1], y[n0-1], z[n0-1], x[n2-1], y[n2-1], z[n2-1], p[0], p[1], p[2])) {
                    *i1 = n0; *i2 = n1; *i3 = n2;
                    /* Dummy barycentric coords for sphere (requires projection) */
                    *b1 = 0.33; *b2 = 0.33; *b3 = 0.33; 
                    return;
                }
                n1 = n2;
                if (n1 == nl) break;
            }
        }
        *i1 = n0; *i2 = nf; *i3 = nl; 
        return;
    }
}

void stri_addnod(int nst, int k, double *x, double *y, double *z,
                 int *list, int *lptr, int *lend, int *lnew, int *ier) {
    int i1, i2, i3, lp, lpf, lpo1, io1, io2, in1;
    double b1, b2, b3;
    double p[3] = {x[k-1], y[k-1], z[k-1]};

    stri_trfind(nst, p, k-1, x, y, z, list, lptr, lend, &b1, &b2, &b3, &i1, &i2, &i3);

    if (i1 == 0) { *ier = -2; return; }

    /* Insert */
    if (i3 != 0) {
        int lp1 = stri_lstptr(lend[i1-1], i2, list, lptr);
        stri_insert(k, lp1, list, lptr, lnew);
        int lp2 = stri_lstptr(lend[i2-1], i3, list, lptr);
        stri_insert(k, lp2, list, lptr, lnew);
        int lp3 = stri_lstptr(lend[i3-1], i1, list, lptr);
        stri_insert(k, lp3, list, lptr, lnew);
        
        list[*lnew-1] = i1; list[*lnew] = i2; list[*lnew+1] = i3;
        lptr[*lnew-1] = *lnew+1; lptr[*lnew] = *lnew+2; lptr[*lnew+1] = *lnew;
        lend[k-1] = *lnew+2;
        *lnew += 3;
    }

    /* Optimise */
    int lpl = lend[k-1];
    lpf = lptr[lpl-1];
    io2 = list[lpf-1];
    lpo1 = lptr[lpf-1];
    io1 = abs(list[lpo1-1]);

    while (1) {
        lp = stri_lstptr(lend[io1-1], io2, list, lptr);
        if (list[lp-1] >= 0) {
            lp = lptr[lp-1];
            in1 = abs(list[lp-1]);
            if (stri_swptst(in1, k, io1, io2, x, y, z)) {
                stri_swap(in1, k, io1, io2, list, lptr, lend, &lpo1);
                if (lpo1 != 0) { io1 = in1; continue; }
            }
        }
        if (lpo1 == lpf || list[lpo1-1] < 0) break;
        io2 = io1;
        lpo1 = lptr[lpo1-1];
        io1 = abs(list[lpo1-1]);
    }
    *ier = 0;
}

void stri_trmesh(int n, double *x, double *y, double *z, 
                 int *list, int *lptr, int *lend, int *lnew, int *ier) {
    if (n < 3) { *ier = -1; return; }

    /* Initial Triangle */
    if (!stri_left(x[0],y[0],z[0], x[1],y[1],z[1], x[2],y[2],z[2])) {
        list[0] = 3;  lptr[0] = 2; list[1] = -2; lptr[1] = 1; lend[0] = 2;
        list[2] = 1;  lptr[2] = 4; list[3] = -3; lptr[3] = 3; lend[1] = 4;
        list[4] = 2;  lptr[4] = 6; list[5] = -1; lptr[5] = 5; lend[2] = 6;
    } else {
        list[0] = 2;  lptr[0] = 2; list[1] = -3; lptr[1] = 1; lend[0] = 2;
        list[2] = 3;  lptr[2] = 4; list[3] = -1; lptr[3] = 3; lend[1] = 4;
        list[4] = 1;  lptr[4] = 6; list[5] = -2; lptr[5] = 5; lend[2] = 6;
    }
    *lnew = 7;

    int *near = (int*)calloc(n, sizeof(int)); // Simplified nearest neighbor structure
    
    for (int k = 4; k <= n; k++) {
        stri_addnod(1, k, x, y, z, list, lptr, lend, lnew, ier);
    }
    free(near);
}
