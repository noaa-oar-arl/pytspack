#include "renka.h"

/* Global tolerance used for geometric checks (initialized in trmesh) */
static double swtol;

/* ------------------------------------------------------------------
   UTILITIES
   ------------------------------------------------------------------ */

bool tri_left(double x1, double y1, double x2, double y2, double x0, double y0) {
    return ((x2 - x1) * (y0 - y1) >= (x0 - x1) * (y2 - y1));
}

double tri_dist_sq(double x1, double y1, double x2, double y2) {
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

/* Returns the index of NB in the adjacency list for N0 (where LPL = LEND(N0)) */
int tri_lstptr(int lpl, int nb, int *list, int *lptr) {
    int lp = lptr[lpl-1];
    while (1) {
        if (abs(list[lp-1]) == nb) return lp;
        lp = lptr[lp-1];
        if (lp == lpl) break;
    }
    return lp;
}

/* Inserts K as a neighbor of N1 following N2 (pointer LP) */
void tri_insert(int k, int lp, int *list, int *lptr, int *lnew) {
    int lsav = lptr[lp-1];
    lptr[lp-1] = *lnew;
    list[*lnew-1] = k;
    lptr[*lnew-1] = lsav;
    (*lnew)++;
}

/* Deletes neighbor NB from N0's adjacency list */
void tri_delnb(int n0, int nb, int n, int *list, int *lptr, int *lend, int *lnew, int *lph) {
    int lpl = lend[n0-1];
    int lpp = lpl;
    int lpb = lptr[lpp-1];

    while (abs(list[lpb-1]) != nb) {
        lpp = lpb;
        lpb = lptr[lpp-1];
        if (lpb == lpl) { *lph = -2; return; } /* Not found */
    }

    /* Delink LPB */
    lptr[lpp-1] = lptr[lpb-1];
    if (lend[n0-1] == lpb) lend[n0-1] = lpp;

    /* Fill hole at LPB with LNEW-1 */
    int lnw = *lnew - 1;
    list[lpb-1] = list[lnw-1];
    lptr[lpb-1] = lptr[lnw-1];

    /* Relink references to LNEW-1 */
    for (int i = 0; i < n; i++) {
        if (lend[i] == lnw) { lend[i] = lpb; break; }
    }
    for (int i = 0; i < lnw-1; i++) {
        if (lptr[i] == lnw) lptr[i] = lpb;
    }

    *lnew = lnw;
    *lph = lpb;
}

/* ------------------------------------------------------------------
   GEOMETRY
   ------------------------------------------------------------------ */

bool tri_swptst(int n1, int n2, int n3, int n4, double *x, double *y) {
    double dx11 = x[n1-1] - x[n3-1];
    double dx12 = x[n2-1] - x[n3-1];
    double dx22 = x[n2-1] - x[n4-1];
    double dx21 = x[n1-1] - x[n4-1];

    double dy11 = y[n1-1] - y[n3-1];
    double dy12 = y[n2-1] - y[n3-1];
    double dy22 = y[n2-1] - y[n4-1];
    double dy21 = y[n1-1] - y[n4-1];

    double cos1 = dx11 * dx12 + dy11 * dy12;
    double cos2 = dx22 * dx21 + dy22 * dy21;

    if (cos1 >= 0.0 && cos2 >= 0.0) return false;
    if (cos1 < 0.0 && cos2 < 0.0) return true;

    double sin1 = dx11 * dy12 - dx12 * dy11;
    double sin2 = dx22 * dy21 - dx21 * dy22;
    double sin12 = sin1 * cos2 + cos1 * sin2;

    return (sin12 > -swtol);
}

void tri_swap(int in1, int in2, int io1, int io2, int *list, int *lptr, int *lend, int *lp21) {
    int lp, lph, lpsav;

    lp = tri_lstptr(lend[in1-1], in2, list, lptr);
    if (abs(list[lp-1]) == in2) { *lp21 = 0; return; }

    /* Delete IO2 from IO1 */
    lp = tri_lstptr(lend[io1-1], in2, list, lptr);
    lph = lptr[lp-1];
    lptr[lp-1] = lptr[lph-1];
    if (lend[io1-1] == lph) lend[io1-1] = lp;

    /* Insert IN2 into IN1 */
    lp = tri_lstptr(lend[in1-1], io1, list, lptr);
    lpsav = lptr[lp-1];
    lptr[lp-1] = lph;
    list[lph-1] = in2;
    lptr[lph-1] = lpsav;

    /* Delete IO1 from IO2 */
    lp = tri_lstptr(lend[io2-1], in1, list, lptr);
    lph = lptr[lp-1];
    lptr[lp-1] = lptr[lph-1];
    if (lend[io2-1] == lph) lend[io2-1] = lp;

    /* Insert IN1 into IN2 */
    lp = tri_lstptr(lend[in2-1], io2, list, lptr);
    lpsav = lptr[lp-1];
    lptr[lp-1] = lph;
    list[lph-1] = in1;
    lptr[lph-1] = lpsav;

    *lp21 = lph;
}

/* ------------------------------------------------------------------
   MESH CONSTRUCTION
   ------------------------------------------------------------------ */

void tri_intadd(int k, int i1, int i2, int i3, int *list, int *lptr, int *lend, int *lnew) {
    int lp;
    /* Add K to I1, I2, I3 */
    lp = tri_lstptr(lend[i1-1], i2, list, lptr);
    tri_insert(k, lp, list, lptr, lnew);
    lp = tri_lstptr(lend[i2-1], i3, list, lptr);
    tri_insert(k, lp, list, lptr, lnew);
    lp = tri_lstptr(lend[i3-1], i1, list, lptr);
    tri_insert(k, lp, list, lptr, lnew);

    /* Add I1, I2, I3 to K */
    list[*lnew-1] = i1; list[*lnew] = i2; list[*lnew+1] = i3;
    lptr[*lnew-1] = *lnew+1; lptr[*lnew] = *lnew+2; lptr[*lnew+1] = *lnew;
    lend[k-1] = *lnew+2;
    *lnew += 3;
}

void tri_bdyadd(int k, int i1, int i2, int *list, int *lptr, int *lend, int *lnew) {
    int lp, lsav, next, nsav;
    /* Add K to I1 */
    lp = lend[i1-1];
    lsav = lptr[lp-1];
    lptr[lp-1] = *lnew;
    list[*lnew-1] = -k; /* Boundary mark */
    lptr[*lnew-1] = lsav;
    lend[i1-1] = *lnew;
    (*lnew)++;
    
    next = -list[lp-1];
    list[lp-1] = next; /* Unmark old boundary */
    nsav = next;

    /* Loop boundary nodes between I1 and I2 */
    while (next != i2) {
        lp = lend[next-1];
        tri_insert(k, lp, list, lptr, lnew);
        next = -list[lp-1];
        list[lp-1] = next;
    }

    /* Add neighbors to K */
    lsav = *lnew;
    list[*lnew-1] = i1;
    lptr[*lnew-1] = *lnew+1;
    (*lnew)++;
    
    next = nsav;
    while (next != i2) {
        list[*lnew-1] = next;
        lptr[*lnew-1] = *lnew+1;
        (*lnew)++;
        lp = lend[next-1];
        next = list[lp-1];
    }
    list[*lnew-1] = -i2; /* New boundary mark */
    lptr[*lnew-1] = lsav;
    lend[k-1] = *lnew;
    (*lnew)++;
}

void tri_trfind(int nst, double px, double py, int n, double *x, double *y, 
                int *list, int *lptr, int *lend, int *i1, int *i2, int *i3) {
    int n0 = (nst < 1 || nst > n) ? 1 : nst;
    int nf, nl, lp, n1, n2;

    while (1) {
        int lpl = lend[n0-1];
        nl = list[lpl-1];
        lp = lptr[lpl-1];
        nf = list[lp-1];
        n1 = nf;

        /* Check Edge N0->N1 */
        if (tri_left(x[n0-1], y[n0-1], x[n1-1], y[n1-1], px, py)) {
            while (1) {
                lp = lptr[lp-1];
                n2 = abs(list[lp-1]);
                if (!tri_left(x[n0-1], y[n0-1], x[n2-1], y[n2-1], px, py)) {
                    /* Found cone */
                    *i1 = n0; *i2 = n1; *i3 = n2; return;
                }
                n1 = n2;
                if (n1 == nl) {
                    /* Fallback or boundary */
                    *i1 = n0; *i2 = nf; *i3 = nl; return;
                }
            }
        }
        
        /* Optimization: Stochastic walk would go here */
        /* Simple fallback */
        *i1 = n0; *i2 = nf; *i3 = nl; return;
    }
}

void tri_addnod(int k, double xk, double yk, int ist, int ncc, int *lcc, int n, 
                double *x, double *y, int *list, int *lptr, int *lend, int *lnew, int *ier) {
    int i1, i2, i3, lp, lpf, lpo1, io1, io2, in1, lp21;
    
    tri_trfind(ist, xk, yk, n, x, y, list, lptr, lend, &i1, &i2, &i3);

    if (i1 == 0) { *ier = -2; return; }
    if (i3 != 0) {
        tri_intadd(k, i1, i2, i3, list, lptr, lend, lnew);
    } else {
        tri_bdyadd(k, i1, i2, list, lptr, lend, lnew);
    }

    /* Optimization Loop */
    int lpl = lend[k-1];
    lpf = lptr[lpl-1];
    io2 = list[lpf-1];
    lpo1 = lptr[lpf-1];
    io1 = abs(list[lpo1-1]);

    while (1) {
        lp = tri_lstptr(lend[io1-1], io2, list, lptr);
        if (list[lp-1] >= 0) {
            lp = lptr[lp-1];
            in1 = abs(list[lp-1]);
            if (tri_swptst(in1, k, io1, io2, x, y)) {
                tri_swap(in1, k, io1, io2, list, lptr, lend, &lp21);
                if (lp21 != 0) { io1 = in1; continue; }
            }
        }
        if (lpo1 == lpf || list[lpo1-1] < 0) break;
        io2 = io1;
        lpo1 = lptr[lpo1-1];
        io1 = abs(list[lpo1-1]);
    }
    *ier = 0;
}

void tri_trmesh(int n, double *x, double *y, int *list, int *lptr, int *lend, int *lnew, int *ier) {
    double d1, d2, d3;
    int k, i, i0;
    int *near, *next;
    double *dist;

    /* Workspace */
    near = (int *)calloc(n, sizeof(int));
    next = (int *)calloc(n, sizeof(int));
    dist = (double *)calloc(n, sizeof(double));

    swtol = 10.0 * 2.22e-16;

    if (n < 3) { *ier = -1; goto cleanup; }

    /* Initial Triangle */
    if (tri_left(x[0], y[0], x[1], y[1], x[2], y[2])) {
        list[0] = 3; lptr[0] = 2; list[1] = -2; lptr[1] = 1; lend[0] = 2;
        list[2] = 1; lptr[2] = 4; list[3] = -3; lptr[3] = 3; lend[1] = 4;
        list[4] = 2; lptr[4] = 6; list[5] = -1; lptr[5] = 5; lend[2] = 6;
    } else {
        list[0] = 2; lptr[0] = 2; list[1] = -3; lptr[1] = 1; lend[0] = 2;
        list[2] = 3; lptr[2] = 4; list[3] = -1; lptr[3] = 3; lend[1] = 4;
        list[4] = 1; lptr[4] = 6; list[5] = -2; lptr[5] = 5; lend[2] = 6;
    }
    *lnew = 7;

    if (n == 3) { *ier = 0; goto cleanup; }

    /* Initial Nearest Neighbors */
    for (k = 4; k <= n; k++) {
        d1 = tri_dist_sq(x[k-1], y[k-1], x[0], y[0]);
        d2 = tri_dist_sq(x[k-1], y[k-1], x[1], y[1]);
        d3 = tri_dist_sq(x[k-1], y[k-1], x[2], y[2]);
        
        if (d1 <= d2 && d1 <= d3) { near[k-1] = 1; dist[k-1] = d1; next[k-1] = near[0]; near[0] = k; }
        else if (d2 <= d3)        { near[k-1] = 2; dist[k-1] = d2; next[k-1] = near[1]; near[1] = k; }
        else                      { near[k-1] = 3; dist[k-1] = d3; next[k-1] = near[2]; near[2] = k; }
    }

    /* Insertion Loop */
    for (k = 4; k <= n; k++) {
        tri_addnod(k, x[k-1], y[k-1], near[k-1], 0, NULL, k-1, x, y, list, lptr, lend, lnew, ier);
        if (*ier != 0) goto cleanup;

        /* Update Nearest Neighbors */
        int knear = near[k-1];
        /* Remove K from list */
        if (near[knear-1] == k) {
            near[knear-1] = next[k-1];
        } else {
            i = near[knear-1];
            while (next[i-1] != k) i = next[i-1];
            next[i-1] = next[k-1];
        }

        /* Re-bucket unprocessed nodes associated with neighbors of K */
        int lpl = lend[k-1];
        int lp = lpl;
        while(1) {
            lp = lptr[lp-1];
            int j = abs(list[lp-1]);
            
            i = near[j-1];
            while (i != 0) {
                int nexti = next[i-1];
                double d = tri_dist_sq(x[k-1], y[k-1], x[i-1], y[i-1]);
                if (d < dist[i-1]) {
                    /* Move I to K's bucket */
                    if (i == near[j-1]) near[j-1] = nexti;
                    else {
                        int prev = near[j-1];
                        while(next[prev-1] != i) prev = next[prev-1];
                        next[prev-1] = nexti;
                    }
                    next[i-1] = near[k-1];
                    near[k-1] = i;
                    dist[i-1] = d;
                }
                i = nexti;
            }
            if (lp == lpl) break;
        }
    }

    *ier = 0;

cleanup:
    free(near);
    free(next);
    free(dist);
}

/* Used by SRFPACK for gradient estimation */
void tri_getnp(int ncc, int *lcc, int n, double *x, double *y, int *list, int *lptr, int *lend, 
               int l, int *npts, double *dist, int *ier) {
    int lm1 = l - 1;
    int n1 = npts[0];
    int ni, lpl, lp, nb, np;
    double dnb, dnp;

    /* Mark existing nodes */
    for (int k = 0; k < lm1; k++) {
        ni = npts[k];
        lend[ni-1] = -abs(lend[ni-1]);
    }

    dnp = 1.0e30;
    np = 0;

    for (int k = 0; k < lm1; k++) {
        ni = npts[k];
        lpl = -lend[ni-1];
        lp = lpl;
        do {
            lp = lptr[lp-1];
            nb = abs(list[lp-1]);
            if (lend[nb-1] < 0) continue; 
            dnb = tri_dist_sq(x[nb-1], y[nb-1], x[n1-1], y[n1-1]);
            if (dnb < dnp) { np = nb; dnp = dnb; }
        } while (lp != lpl);
    }

    /* Unmark */
    for (int k = 0; k < lm1; k++) {
        ni = npts[k];
        lend[ni-1] = abs(lend[ni-1]);
    }

    if (np != 0) {
        npts[l-1] = np;
        dist[l-1] = sqrt(dnp);
        *ier = 0;
    } else {
        *ier = 1;
    }
}
