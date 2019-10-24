# cython: language_level=3
# cython: boundscheck=False

import numpy as np

cimport libsplash

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
    bint     useaccelerate,
    bint     iverbose
):

    cdef float[:, ::1] datsmooth = np.empty((npixx, npixy), dtype=np.single)

    libsplash.interpolate3d_projection_c(
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
        &useaccelerate,
        &iverbose
    )

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
    float    dscreen,
    bint     iverbose
):

    cdef float[:, ::1] vecsmoothx = np.empty((npixx, npixy), dtype=np.single)
    cdef float[:, ::1] vecsmoothy = np.empty((npixx, npixy), dtype=np.single)

    libsplash.interpolate3d_proj_vec_c(
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
        &dscreen,
        &iverbose
    )

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
    bint     normalise,
    bint     iverbose
):

    cdef float[:, ::1] datsmooth = np.empty((npixx, npixy), dtype=np.single)

    libsplash.interpolate3d_fastxsec_c(
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
        &normalise,
        &iverbose
    )

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
    bint     normalise,
    bint     iverbose
):

    cdef float[:, ::1] vecsmoothx = np.empty((npixx, npixy), dtype=np.single)
    cdef float[:, ::1] vecsmoothy = np.empty((npixx, npixy), dtype=np.single)

    libsplash.interpolate3d_xsec_vec_c(
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
        &normalise,
        &iverbose
    )

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
    float    zcut,
    bint     iverbose
):

    cdef float[:, ::1] datsmooth  = np.empty((npixx, npixy), dtype=np.single)
    cdef float[:, ::1] brightness = np.empty((npixx, npixy), dtype=np.single)

    libsplash.interp3d_proj_opacity_c(
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
        &zcut,
        &iverbose
    )

    return datsmooth

def interpolate3d_proj_geom(
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
    int      igeom,
    int      iplotx,
    int      iploty,
    int      iplotz,
    int[:]   ix,
    float[:] xorigin
):

    cdef float[:, ::1] datsmooth = np.empty((npixx, npixy), dtype=np.single)

    libsplash.interpolate3d_proj_geom_c(
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
        &igeom,
        &iplotx,
        &iploty,
        &iplotz,
        &ix[0],
        &xorigin[0]
    )

    return datsmooth

def interpolate3d_xsec_geom(
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
    bint     normalise,
    int      igeom,
    int      iplotx,
    int      iploty,
    int      iplotz,
    int[:]   ix,
    float[:] xorigin
):

    cdef float[:, ::1] datsmooth = np.empty((npixx, npixy), dtype=np.single)

    libsplash.interpolate3d_xsec_geom_c(
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
        &normalise,
        &igeom,
        &iplotx,
        &iploty,
        &iplotz,
        &ix[0],
        &xorigin[0]
    )

    return datsmooth
