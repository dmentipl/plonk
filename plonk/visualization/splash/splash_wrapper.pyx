# cython: language_level=3
# cython: boundscheck=False

from numpy cimport ndarray
import numpy as np

cdef extern:
    void c_interpolate3d_projection(
        float *x,
        float *y,
        float *z,
        float *hh,
        float *weight,
        float *dat,
        int   *itype,
        int   *npart,
        float *xmin,
        float *ymin,
        float *datsmooth,
        int   *npixx,
        int   *npixy,
        float *pixwidthx,
        float *pixwidthy,
        bint  *normalise,
        float *zobserver,
        float *dscreen,
        bint  *useaccelerate
        )

def interpolate3d_projection(
    float[:] x,
    float[:] y,
    float[:] z,
    float[:] hh,
    float[:] weight,
    float[:] dat,
    int[:]   itype,
    int      npart,
    float    xmin,
    float    ymin,
    int      npixx,
    int      npixy,
    float    pixwidthx,
    float    pixwidthy,
    bint     normalise,
    float    zobserver,
    float    dscreen,
    bint     useaccelerate ):

    cdef float[:, ::1] datsmooth = np.empty((npixx, npixy), dtype=np.single)

    c_interpolate3d_projection(
        &x[0],
        &y[0],
        &z[0],
        &hh[0],
        &weight[0],
        &dat[0],
        &itype[0],
        &npart,
        &xmin,
        &ymin,
        &datsmooth[0,0],
        &npixx,
        &npixy,
        &pixwidthx,
        &pixwidthy,
        &normalise,
        &zobserver,
        &dscreen,
        &useaccelerate )

    return datsmooth
