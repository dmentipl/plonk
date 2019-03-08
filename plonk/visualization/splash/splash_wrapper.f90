module splash_wrapper

  use splash,        only: interpolate3d_projection
  use iso_c_binding, only: c_float, c_int, c_bool

  implicit none

contains

subroutine c_interpolate3d_projection(                                        &
    x, y, z, hh, weight, dat, itype, npart, xmin, ymin, datsmooth, npixx,     &
    npixy, pixwidthx, pixwidthy, normalise, zobserver, dscreen, useaccelerate &
    ) bind(c)

  real(c_float),   intent(in)  :: x(npart),      &
                                  y(npart),      &
                                  z(npart),      &
                                  hh(npart),     &
                                  weight(npart), &
                                  dat(npart),    &
                                  xmin,          &
                                  ymin,          &
                                  pixwidthx,     &
                                  pixwidthy,     &
                                  zobserver,     &
                                  dscreen
  integer(c_int),  intent(in)  :: npart,         &
                                  npixx,         &
                                  npixy,         &
                                  itype(npart)
  logical(c_bool), intent(in)  :: normalise,     &
                                  useaccelerate
  real(c_float),   intent(out) :: datsmooth(npixx,npixy)

  logical :: normalise_f, &
             useaccelerate_f

  normalise_f     = normalise
  useaccelerate_f = useaccelerate

  call interpolate3d_projection(x, y, z, hh, weight, dat, itype, npart,     &
     xmin, ymin, datsmooth, npixx, npixy, pixwidthx, pixwidthy, normalise_f, &
     zobserver, dscreen, useaccelerate_f)

end subroutine c_interpolate3d_projection

end module splash_wrapper

