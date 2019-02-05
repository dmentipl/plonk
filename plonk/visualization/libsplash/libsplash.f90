!-----------------------------------------------------------------
!
!  This file is (or was) part of SPLASH, a visualisation tool
!  for Smoothed Particle Hydrodynamics written by Daniel Price:
!
!  http://users.monash.edu.au/~dprice/splash
!
!  SPLASH comes with ABSOLUTELY NO WARRANTY.
!  This is free software; and you are welcome to redistribute
!  it under the terms of the GNU General Public License
!  (see LICENSE file for details) and the provision that
!  this notice remains intact. If you modify this file, please
!  note section 2a) of the GPLv2 states that:
!
!  a) You must cause the modified files to carry prominent notices
!     stating that you changed the files and the date of any change.
!
!  Copyright (C) 2005-2016 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!-----------------------------------------------------------------

module libsplash

 implicit none

 public :: interpolate2D
 public :: interpolate2D_vec
 public :: interpolate2D_xsec
 public :: interpolate_part
 public :: interpolate_part1
 public :: interpolate2D_pixels
 public :: set_interpolation_weights

 private

 private :: w_cubic
 private :: w_quartic
 private :: w_quintic
 private :: w_quartic2h
 private :: w_wendlandc2
 private :: w_wendlandc4
 private :: w_wendlandc6
 private :: pint
 private :: full_2d_mod
 private :: F1_2d
 private :: F2_2d
 private :: F3_2d
 private :: get_sink_type

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_j w_j dat_j W(r-r_j, h_j)
!
!     where _j is the quantity at the neighbouring particle j and
!     W is the smoothing kernel, for which we use the usual cubic spline.
!     For an SPH interpolation the weight for each particle should be
!     the dimensionless quantity
!
!     w_j = m_j / (rho_j * h_j**ndim)
!
!     Other weights can be used (e.g. constants), but in this case the
!     normalisation option should also be set.
!
!     Input: particle coordinates  : x,y    (npart)
!            smoothing lengths     : hh     (npart)
!            interpolation weights : weight (npart)
!            scalar data to smooth : dat    (npart)
!
!            number of pixels in x,y : npixx,npixy
!            pixel width             : pixwidth
!            option to normalise interpolation : normalise (.true. or .false.)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price 2003-2012
!     Exact rendering implemented by Maya Petkova and Daniel Price 2018
!--------------------------------------------------------------------------

interface
subroutine interpolate2D(x,y,hh,weight,dat,itype,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidthx,pixwidthy,&
     normalise,exact,periodicx,periodicy)

  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise,exact,periodicx,periodicy
end subroutine interpolate2D
end interface

!--------------------------------------------------------------------------
!
!     ** this version does vector quantities
!
!     Input: particle coordinates  : x,y   (npart)
!            smoothing lengths     : hh     (npart)
!            interpolation weights : weight (npart)
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!
!     Output: smoothed vector field    : vecsmoothx (npixx,npixy)
!                                      : vecsmoothy (npixx,npixy)
!
!     Daniel Price, University of Exeter, March 2005
!     Exact rendering implemented by Maya Petkova and Daniel Price 2018å
!--------------------------------------------------------------------------

interface
subroutine interpolate2D_vec(x,y,hh,weight,vecx,vecy,itype,npart, &
     xmin,ymin,vecsmoothx,vecsmoothy,npixx,npixy,pixwidthx,pixwidthy,&
     normalise,periodicx,periodicy)

  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,hh,weight,vecx,vecy
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx,vecsmoothy
  logical, intent(in) :: normalise,periodicx,periodicy

end subroutine interpolate2D_vec
end interface

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     this version takes any 1D cross section through a 2D data set
!     the 1D line is specified by two points, (x1,y1) and (x2,y2)
!     (ie. this is for arbitrary oblique cross sections)
!
!     NB: A similar version could be used for 2D oblique cross sections
!         of 3D data. In this case we would need to find the intersection
!         between the smoothing sphere and the cross section plane. However
!         in 3D it is simpler just to rotate the particles first and then take
!         a straight cross section.
!
!     Input: particle coordinates  : x,y   (npart)
!            smoothing lengths     : hh     (npart)
!            interpolation weights : weight (npart)
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx)
!
!     Daniel Price, Institute of Astronomy, Cambridge, Feb 2004
!--------------------------------------------------------------------------

interface
subroutine interpolate2D_xsec(x,y,hh,weight,dat,itype,npart,&
     x1,y1,x2,y2,datsmooth,npixx,normalise)

  integer, intent(in) :: npart,npixx
  real, intent(in), dimension(npart) :: x,y,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: x1,y1,x2,y2
  real, intent(out), dimension(npixx) :: datsmooth
  logical, intent(in) :: normalise
end subroutine interpolate2D_xsec
end interface

!--------------------------------------------------------------------------
!     subroutine to render particles onto a pixel array
!     at the maximum or minimum colour
!
!     Written by Daniel Price 21/7/2008
!--------------------------------------------------------------------------

interface
subroutine interpolate_part(x,y,hh,npart,xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval,brightness)

  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,hh
  real, intent(in) :: xmin,ymin,pixwidth,datval
  real, intent(inout), dimension(npixx,npixy) :: datsmooth
  real, intent(inout), dimension(npixx,npixy), optional :: brightness
end subroutine interpolate_part
end interface

!--------------------------------------------------------------------------
!     subroutine to render a single particle onto a pixel array
!
!     Written by Daniel Price 21/7/2008
!--------------------------------------------------------------------------

interface
subroutine interpolate_part1(xi,yi,hi,xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval,brightness)

  real, intent(in) :: xi,yi,hi,xmin,ymin,pixwidth,datval
  integer, intent(in) :: npixx,npixy
  real, intent(inout), dimension(npixx,npixy) :: datsmooth
  real, intent(inout), dimension(npixx,npixy), optional :: brightness
end subroutine interpolate_part1
end interface

!--------------------------------------------------------------------------
!     subroutine to interpolate from arbitrary data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_j dat_j W(r-r_j, \Delta) / sum_j W(r-r_j, \Delta)
!
!     where _j is the quantity at the neighbouring particle j and
!     W is the smoothing kernel. The interpolation is normalised.
!
!     Input: data points  : x,y    (npart)
!            third scalar to use as weight : dat (npart) (optional)
!
!            number of pixels in x,y : npixx,npixy
!            pixel width             : pixwidth
!            option to normalise interpolation : normalise (.true. or .false.)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price 2017
!--------------------------------------------------------------------------

interface
subroutine interpolate2D_pixels(x,y,itype,npart, &
     xmin,ymin,xmax,ymax,datsmooth,npixx,npixy,&
     normalise,adaptive,dat,datpix2)

  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,xmax,ymax
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise,adaptive
  real, intent(in), dimension(npart), optional :: dat
  real, dimension(npixx,npixy), intent(out), optional :: datpix2
end subroutine interpolate2D_pixels
end interface

!---------------------------------------
!
!  Functional forms of various kernels
!
!--------------------------------------

interface
pure real function w_cubic(q2)
 real, intent(in) :: q2
end function w_cubic
end interface

interface
pure real function w_quartic(q2)
 real, intent(in) :: q2
end function w_quartic
end interface

interface
pure real function w_quintic(q2)
 real, intent(in) :: q2
end function w_quintic
end interface

interface
pure real function w_quartic2h(q2)
 real, intent(in) :: q2
end function w_quartic2h
end interface

interface
pure real function w_wendlandc2(q2)
 real, intent(in) :: q2
end function w_wendlandc2
end interface

interface
pure real function w_wendlandc4(q2)
 real, intent(in) :: q2
end function w_wendlandc4
end interface

interface
pure real function w_wendlandc6(q2)
 real, intent(in) :: q2
end function w_wendlandc6
end interface


!---------------------------------------
!
!  Integrals of the kernel, used
!  when performing exact interpolation
!  See Petkova, Laibe & Bonnell (2018)
!
!---------------------------------------

interface
pure real function pint(r0, d1, d2, hi1)
 real, intent(in) :: r0, d1, d2, hi1
end function pint
end interface

!---------------------------------------
!
! Helper functions for kernel integrals
! See Petkova, Laibe & Bonnell (2018)
!
!---------------------------------------
interface
pure real function full_2d_mod(phi,tphi,q0)
 real, intent(in) :: phi,tphi,q0
end function full_2d_mod
end interface

interface
pure real function F1_2d(phi,tphi,cphi,q0)
 real, intent(in) :: phi,tphi,cphi,q0
end function F1_2d
end interface

interface
pure real function F2_2d(phi, tphi, cphi, q0)
 real, intent(in) :: phi, tphi, cphi, q0
end function F2_2d
end interface

interface
pure real function F3_2d(phi)
 real, intent(in) :: phi
end function F3_2d
end interface

!-------------------------------------------------------------------
! Set interpolation weights for the particles. The weights are
! calculated using:
!
!           w = m/(rho*h**ndim),
!
! where we need to handle a few special scenarios:
!
! 1) Firstly, the weight should be calculated in a consistent
!    set of units. Safest way is to use the data as originally
!    read from the dump file, before any unit scaling was applied.
!
! 2) Particle weights are set to zero for particle types not
!    used in the rendering.
!
! 3) If particle mass not read, it is still possible to perform
!    interpolations, but not using the SPH weights. These
!    interpolations *must* therefore be normalised.
!-------------------------------------------------------------------

interface
subroutine set_interpolation_weights(weighti,dati,iamtypei,usetype, &
           ninterp,npartoftype,masstype,ntypes,ndataplots,irho,ipmass,ih,ndim, &
           iRescale,idensityweighted,inormalise,units,unit_interp,required, &
           rendersinks)

  real, dimension(:), intent(out)              :: weighti
  real, dimension(:,:), intent(in)             :: dati
  integer, dimension(:), intent(in) :: iamtypei
  logical, dimension(:), intent(in)            :: usetype
  logical, dimension(:), intent(in)    :: required
  integer, intent(in)                          :: ih,irho,ipmass,ndim
  integer, intent(in)                          :: ninterp,ntypes,ndataplots
  integer, dimension(:), intent(in)            :: npartoftype
  real,    dimension(:), intent(in)            :: masstype
  logical, intent(in)                          :: iRescale,idensityweighted,rendersinks
  logical, intent(inout)                       :: inormalise
  real, dimension(:), intent(in)       :: units
  real, intent(in)                  :: unit_interp
end subroutine set_interpolation_weights
end interface

!-----------------------------------------------------------------
!
! utility to "guess" which particle type contains sink particles
! from the label
!
!-----------------------------------------------------------------

interface
integer function get_sink_type(ntypes)
 integer, intent(in) :: ntypes
end function get_sink_type
end interface

end module libsplash
