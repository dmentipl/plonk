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

    void c_interpolate3d_proj_vec(
        float *x,
        float *y,
        float *z,
        float *hh,
        float *weight,
        float *vecx,
        float *vecy,
        int   *itype,
        int   *npart,
        float *xmin,
        float *ymin,
        float *vectsmoothx,
        float *vectsmoothy,
        int   *npixx,
        int   *npixy,
        float *pixwidthx,
        float *pixwidthy,
        bint  *normalise,
        float *zobserver,
        float *dscreen
        )

    void c_interpolate3d_fastxsec(
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
        float *zslice,
        float *datsmooth,
        int   *npixx,
        int   *npixy,
        float *pixwidthx,
        float *pixwidthy,
        bint  *normalise
        )

    void c_interpolate3d_xsec_vec(
        float *x,
        float *y,
        float *z,
        float *hh,
        float *weight,
        float *vecx,
        float *vecy,
        int   *itype,
        int   *npart,
        float *xmin,
        float *ymin,
        float *zslice,
        float *vectsmoothx,
        float *vectsmoothy,
        int   *npixx,
        int   *npixy,
        float *pixwidthx,
        float *pixwidthy,
        bint  *normalise
        )

    void c_interp3d_proj_opacity(
        float *x,
        float *y,
        float *z,
        float *pmass,
        int   *npmass,
        float *hh,
        float *weight,
        float *dat,
        float *zorig,
        int   *itype,
        int   *npart,
        float *xmin,
        float *ymin,
        float *datsmooth,
        float *brightness,
        int   *npixx,
        int   *npixy,
        float *pixwidth,
        float *zobserver,
        float *dscreenfromobserver,
        float *rkappa,
        float *zcut
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

def interpolate3d_proj_vec(
    float[:] x,
    float[:] y,
    float[:] z,
    float[:] hh,
    float[:] weight,
    float[:] vecx,
    float[:] vecy,
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
    float    dscreen ):

    cdef float[:, ::1] vecsmoothx = np.empty((npixx, npixy), dtype=np.single)
    cdef float[:, ::1] vecsmoothy = np.empty((npixx, npixy), dtype=np.single)

    c_interpolate3d_proj_vec(
        &x[0],
        &y[0],
        &z[0],
        &hh[0],
        &weight[0],
        &vecx[0],
        &vecy[0],
        &itype[0],
        &npart,
        &xmin,
        &ymin,
        &vecsmoothx[0,0],
        &vecsmoothy[0,0],
        &npixx,
        &npixy,
        &pixwidthx,
        &pixwidthy,
        &normalise,
        &zobserver,
        &dscreen )

    return vecsmoothx, vecsmoothy

def interpolate3d_fastxsec(
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
    float    zslice,
    int      npixx,
    int      npixy,
    float    pixwidthx,
    float    pixwidthy,
    bint     normalise ):

    cdef float[:, ::1] datsmooth = np.empty((npixx, npixy), dtype=np.single)

    c_interpolate3d_fastxsec(
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
        &zslice,
        &datsmooth[0,0],
        &npixx,
        &npixy,
        &pixwidthx,
        &pixwidthy,
        &normalise )

    return datsmooth

def interpolate3d_xsec_vec(
    float[:] x,
    float[:] y,
    float[:] z,
    float[:] hh,
    float[:] weight,
    float[:] vecx,
    float[:] vecy,
    int[:]   itype,
    int      npart,
    float    xmin,
    float    ymin,
    float    zslice,
    int      npixx,
    int      npixy,
    float    pixwidthx,
    float    pixwidthy,
    bint     normalise ):

    cdef float[:, ::1] vecsmoothx = np.empty((npixx, npixy), dtype=np.single)
    cdef float[:, ::1] vecsmoothy = np.empty((npixx, npixy), dtype=np.single)

    c_interpolate3d_xsec_vec(
        &x[0],
        &y[0],
        &z[0],
        &hh[0],
        &weight[0],
        &vecx[0],
        &vecy[0],
        &itype[0],
        &npart,
        &xmin,
        &ymin,
        &zslice,
        &vecsmoothx[0,0],
        &vecsmoothy[0,0],
        &npixx,
        &npixy,
        &pixwidthx,
        &pixwidthy,
        &normalise )

    return vecsmoothx, vecsmoothy

def interp3d_proj_opacity(
    float[:] x,
    float[:] y,
    float[:] z,
    float[:] pmass,
    int      npmass,
    float[:] hh,
    float[:] weight,
    float[:] dat,
    float[:] zorig,
    int[:]   itype,
    int      npart,
    float    xmin,
    float    ymin,
    int      npixx,
    int      npixy,
    float    pixwidth,
    float    zobserver,
    float    dscreenfromobserver,
    float    rkappa,
    float    zcut ):

    cdef float[:, ::1] datsmooth  = np.empty((npixx, npixy), dtype=np.single)
    cdef float[:, ::1] brightness = np.empty((npixx, npixy), dtype=np.single)

    c_interp3d_proj_opacity(
        &x[0],
        &y[0],
        &z[0],
        &pmass[0],
        &npmass,
        &hh[0],
        &weight[0],
        &dat[0],
        &zorig[0],
        &itype[0],
        &npart,
        &xmin,
        &ymin,
        &datsmooth[0,0],
        &brightness[0,0],
        &npixx,
        &npixy,
        &pixwidth,
        &zobserver,
        &dscreenfromobserver,
        &rkappa,
        &zcut )

    return datsmooth
