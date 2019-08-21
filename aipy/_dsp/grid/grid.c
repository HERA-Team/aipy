#include "grid.h"

int grid1D_r(float *buf, long buflen, 
        float *inds, float *data, long datalen, long footprint) {
    long i, j, jmod;
    float find, fdat, fwgt;
    for (i = 0; i < datalen; i++) {
        find = inds[i];
        fdat = data[i];
        for (j = floorf(find-footprint/2); j <= ceilf(find+footprint/2); j++) {
            fwgt = find - j;
            jmod = j % buflen;
            jmod = jmod < 0 ? jmod + buflen : jmod;
            fwgt = 0.79788456080286541 * exp(-2*fwgt*fwgt); // Gaussian, sig=0.5
            buf[jmod] += fwgt * fdat;
        }
    }
    return 0;
}

int grid1D_c(float *buf, long buflen, 
        float *inds, float *data, long datalen, long footprint) {
    long i, j, jmod;
    float find, fdatr, fdati, fwgt;
    for (i = 0; i < datalen; i++) {
        find = inds[i];
        fdatr = data[2*i];
        fdati = data[2*i+1];
        for (j = floorf(find-footprint/2); j <= ceilf(find+footprint/2); j++) {
            fwgt = find - j;
            jmod = j % buflen;
            jmod = jmod < 0 ? jmod + buflen : jmod;
            fwgt = 0.79788456080286541 * exp(-2*fwgt*fwgt); // Gaussian, sig=0.5
            buf[2*jmod]   += fwgt * fdatr;
            buf[2*jmod+1] += fwgt * fdati;
        }
    }
    return 0;
}

int grid2D_c(float *buf, long buflen1, long buflen2,
        float *ind1, float *ind2, float *data, long datalen, long footprint) {
    long i, j1, j2, j1mod, j2mod;
    float find1, find2, fdatr, fdati, fwgt1, fwgt2;
    for (i = 0; i < datalen; i++) {
        find1 = ind1[i];
        find2 = ind2[i];
        fdatr = data[2*i];
        fdati = data[2*i+1];
        for (j1 = floorf(find1-footprint/2); j1 <= ceilf(find1+footprint/2); j1++) {
          j1mod = j1 % buflen1;
          j1mod = j1mod < 0 ? j1mod + buflen1 : j1mod;
          fwgt1 = find1 - j1;
          fwgt1 *= fwgt1;
          for (j2 = floorf(find2-footprint/2); j2 <= ceilf(find2+footprint/2); j2++) {
            j2mod = j2 % buflen2;
            j2mod = j2mod < 0 ? j2mod + buflen2 : j2mod;
            fwgt2 = find2 - j2;
            fwgt2 *= fwgt2;
            fwgt2 = 0.63661977236758149 * exp(-2*(fwgt1+fwgt2)); // 2D Gaussian, sigx,y=0.5
            // XXX should really make sure wgts sum to 1
            buf[2*(j1mod*buflen1+j2mod)]   += fwgt2 * fdatr;
            buf[2*(j1mod*buflen1+j2mod)+1] += fwgt2 * fdati;
          }
        }
    }
    return 0;
}

int degrid2D_c(float *buf, long buflen1, long buflen2,
        float *ind1, float *ind2, float *data, long datalen, long footprint) {
    long i, j1, j2, j1mod, j2mod;
    float find1, find2, fwgt1, fwgt2, tot_wgt;
    for (i = 0; i < datalen; i++) {
        find1 = ind1[i];
        find2 = ind2[i];
        tot_wgt = 0;
        for (j1 = floorf(find1-footprint/2); j1 <= ceilf(find1+footprint/2); j1++) {
          j1mod = j1 % buflen1;
          j1mod = j1mod < 0 ? j1mod + buflen1 : j1mod;
          fwgt1 = find1 - j1;
          fwgt1 *= fwgt1;
          for (j2 = floorf(find2-footprint/2); j2 <= ceilf(find2+footprint/2); j2++) {
            j2mod = j2 % buflen2;
            j2mod = j2mod < 0 ? j2mod + buflen2 : j2mod;
            fwgt2 = find2 - j2;
            fwgt2 *= fwgt2;
            fwgt2 = 0.63661977236758149 * exp(-2*(fwgt1+fwgt2)); // 2D Gaussian, sigx,y=0.5
            tot_wgt += fwgt2;
            data[2*i] += fwgt2 * buf[2*(j1mod*buflen1+j2mod)];
            data[2*i+1] += fwgt2 * buf[2*(j1mod*buflen1+j2mod)+1];
          }
        }
        data[2*i] /= tot_wgt;
        data[2*i+1] /= tot_wgt;
    }
    return 0;
}
