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

module splash

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

 real,    parameter :: pi          = 3.1415926536
 real,    parameter :: radkernel   = 2.
 real,    parameter :: radkernel2  = 4.
 real,    parameter :: cnormk2D    = 0.4547284088
 real,    parameter :: weight_sink = -1.

 integer, parameter :: maxparttypes = 12  ! max # of different particle types

 character(len=20), dimension(maxparttypes) :: labeltype

contains

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

subroutine interpolate2D(x,y,hh,weight,dat,itype,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidthx,pixwidthy,&
     normalise,exact,periodicx,periodicy)

  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise,exact,periodicx,periodicy
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: ipixi,jpixi
  real :: hi,hi1,radkern,q2,wab,const
  real :: term,termnorm,dx,dy,xpix,ypix

! Maya
  real :: pixint,d1,d2,r0

  datsmooth = 0.
  datnorm = 0.
  if (exact) then
     print "(1x,a)",'interpolating from particles to 2D grid (exact)...'
  elseif (normalise) then
     print "(1x,a)",'interpolating from particles to 2D grid (normalised)...'
  else
     print "(1x,a)",'interpolating from particles to 2D grid (non-normalised)...'
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print "(1x,a)",'interpolate2D: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate2D: warning: ignoring some or all particles with h < 0'
  endif
  const = cnormk2D  ! normalisation constant
  !
  !--loop over particles
  !
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts
     !
     !--skip particles with zero weights
     !
     termnorm = const*weight(i)
     if (termnorm.le.0.) cycle over_parts
     !
     !--skip particles with wrong h's
     !
     hi = hh(i)
     if (hi.le.tiny(hi)) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi1 = 1./hi
     radkern = radkernel*hi  ! radius of the smoothing kernel
     term = termnorm*dat(i)
     !
     !--for each particle work out which pixels it contributes to
     !
     ipixmin = int((x(i) - radkern - xmin)/pixwidthx)
     jpixmin = int((y(i) - radkern - ymin)/pixwidthy)
     ipixmax = int((x(i) + radkern - xmin)/pixwidthx) + 1
     jpixmax = int((y(i) + radkern - ymin)/pixwidthy) + 1

     if (.not.periodicx) then
        if (ipixmin.lt.1)     ipixmin = 1      ! make sure they only contribute
        if (ipixmax.gt.npixx) ipixmax = npixx  ! to pixels in the image
     endif
     if (.not.periodicy) then
        if (jpixmin.lt.1)     jpixmin = 1
        if (jpixmax.gt.npixy) jpixmax = npixy
     endif

     if (exact) then
        !
        !--loop over pixels boundaries, adding the contribution from this particle
        !
        !--first pixel row
        !
        if(jpixmax.ge.jpixmin) then
           jpix = jpixmin
           jpixi = jpix
           if (periodicy) then
              if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
              if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1
           endif
           ypix = ymin + (jpix-0.5)*pixwidthy
           dy = ypix - y(i)

           do ipix = ipixmin,ipixmax
              ipixi = ipix
              if (periodicx) then
                 if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
                 if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
              endif
              xpix = xmin + (ipix-0.5)*pixwidthx
              dx = xpix - x(i)

              !--top boundary
              r0 = 0.5*pixwidthy - dy
              d1 = 0.5*pixwidthx + dx
              d2 = 0.5*pixwidthx - dx
              pixint = pint(r0, d1, d2, hi1)

              wab = pixint /pixwidthx/pixwidthy/const*hi**2
              datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
              if (normalise) datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab
           enddo
        end if

        !
        !--first pixel column
        !
        if(ipixmax.ge.ipixmin) then
           ipix = ipixmin
           ipixi = ipix
           if (periodicx) then
              if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
              if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
           endif
           xpix = xmin + (ipix-0.5)*pixwidthx
           dx = xpix - x(i)

           do jpix = jpixmin,jpixmax
              jpixi = jpix
              if (periodicy) then
                 if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
                 if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1
              endif
              ypix = ymin + (jpix-0.5)*pixwidthy
              dy = ypix - y(i)

              !--left boundary
              r0 = 0.5*pixwidthx - dx
              d1 = 0.5*pixwidthy - dy
              d2 = 0.5*pixwidthy + dy
              pixint = pint(r0, d1, d2, hi1)

              wab = pixint /pixwidthx/pixwidthy/const*hi**2
              datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
              if (normalise) datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab
           enddo
        end if

        !
        !--other pixels
        !
        do jpix = jpixmin,jpixmax
           jpixi = jpix
           if (periodicy) then
              if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
              if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1
           endif
           ypix = ymin + (jpix-0.5)*pixwidthy
           dy = ypix - y(i)

           do ipix = ipixmin,ipixmax
              ipixi = ipix
              if (periodicx) then
                 if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
                 if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
              endif
              xpix = xmin + (ipix-0.5)*pixwidthx
              dx = xpix - x(i)
              !
              !--Kernel integral
              !
              !--bottom boundary
              r0 = 0.5*pixwidthy + dy
              d1 = 0.5*pixwidthx - dx
              d2 = 0.5*pixwidthx + dx
              pixint = pint(r0, d1, d2, hi1)

              wab = pixint /pixwidthx/pixwidthy/const*hi**2
              datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
              if (normalise) datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab

              if(jpix < jpixmax) then
                 jpixi = jpix+1
                 if (periodicy) then
                    if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
                    if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1
                 endif

                 datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) - term*wab
                 if (normalise) datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) - termnorm*wab

                 jpixi = jpix
                 if (periodicy) then
                    if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
                    if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1
                 endif
              end if

              !--right boundary
              r0 = 0.5*pixwidthx + dx
              d1 = 0.5*pixwidthy + dy
              d2 = 0.5*pixwidthy - dy
              pixint = pint(r0, d1, d2, hi1)

              wab = pixint /pixwidthx/pixwidthy/const*hi**2
              datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
              if (normalise) datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab

              if(ipix < ipixmax) then
                 ipixi = ipix+1
                 if (periodicx) then
                    if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
                    if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
                 endif

                 datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) - term*wab
                 if (normalise) datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) - termnorm*wab
              end if
           enddo
        enddo
     else
        !
        !--loop over pixels, adding the contribution from this particle
        !
        do jpix = jpixmin,jpixmax
           jpixi = jpix
           if (periodicy) then
              if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
              if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1
           endif
           ypix = ymin + (jpix-0.5)*pixwidthy
           dy = ypix - y(i)
           do ipix = ipixmin,ipixmax
              ipixi = ipix
              if (periodicx) then
                 if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
                 if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
              endif
              xpix = xmin + (ipix-0.5)*pixwidthx
              dx = xpix - x(i)
              q2 = (dx*dx + dy*dy)*hi1*hi1
              !
              !--SPH kernel
              !
              wab = w_cubic(q2)
              !
              !--calculate data value at this pixel using the summation interpolant
              !
              datsmooth(ipixi,jpixi) = datsmooth(ipixi,jpixi) + term*wab
              if (normalise) datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab
           enddo
        enddo
     end if

  enddo over_parts
  if (exact) then
     print*, 'sum of datpix = ', sum(datsmooth)/(npixx*npixy)
     print*, 'max of datpix = ', maxval(datsmooth)
     print*, 'min of datpix = ', minval(datsmooth)
  endif
  !
  !--normalise dat array
  !
  if (normalise) then
     where (datnorm > 0.)
        datsmooth = datsmooth/datnorm
     end where
  endif

  return
end subroutine interpolate2D

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

subroutine interpolate2D_vec(x,y,hh,weight,vecx,vecy,itype,npart, &
     xmin,ymin,vecsmoothx,vecsmoothy,npixx,npixy,pixwidthx,pixwidthy,&
     normalise,periodicx,periodicy)

  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,hh,weight,vecx,vecy
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx,vecsmoothy
  logical, intent(in) :: normalise,periodicx,periodicy
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  integer :: ipixi,jpixi
  real :: hi,hi1,radkern,q2,wab,const
  real :: termnorm,termx,termy,dx,dy,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  datnorm = 0.
  if (normalise) then
     print "(1x,a)",'interpolating vector field from particles to 2D grid (normalised)...'
  else
     print "(1x,a)",'interpolating vector field from particles to 2D grid (non-normalised)...'
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print*,'interpolate2D_vec: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate2D_vec: warning: ignoring some or all particles with h < 0'
  endif
  const = cnormk2D  ! normalisation constant
  !
  !--loop over particles
  !
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts
     !
     !--skip particles with zero weights
     !
     termnorm = const*weight(i)
     if (termnorm.le.0.) cycle over_parts
     !
     !--skip particles with wrong h's
     !
     hi = hh(i)
     if (hi.le.tiny(hi)) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi1 = 1./hi
     radkern = radkernel*hi  ! radius of the smoothing kernel
     termx = termnorm*vecx(i)
     termy = termnorm*vecy(i)
     !
     !--for each particle work out which pixels it contributes to
     !
     ipixmin = int((x(i) - radkern - xmin)/pixwidthx)
     jpixmin = int((y(i) - radkern - ymin)/pixwidthy)
     ipixmax = int((x(i) + radkern - xmin)/pixwidthx) + 1
     jpixmax = int((y(i) + radkern - ymin)/pixwidthy) + 1

     if (.not.periodicx) then
        if (ipixmin.lt.1)     ipixmin = 1
        if (ipixmax.gt.npixx) ipixmax = npixx
     endif
     if (.not.periodicy) then
        if (jpixmin.lt.1)     jpixmin = 1
        if (jpixmax.gt.npixy) jpixmax = npixy
     endif
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        jpixi = jpix
        if (periodicy) then
           if (jpixi.lt.1)     jpixi = mod(jpixi,npixy) + npixy
           if (jpixi.gt.npixy) jpixi = mod(jpixi-1,npixy) + 1
        endif
        ypix = ymin + (jpix-0.5)*pixwidthy
        dy = ypix - y(i)
        do ipix = ipixmin,ipixmax
           ipixi = ipix
           if (periodicx) then
              if (ipixi.lt.1)     ipixi = mod(ipixi,npixx) + npixx
              if (ipixi.gt.npixx) ipixi = mod(ipixi-1,npixx) + 1
           endif
           xpix = xmin + (ipix-0.5)*pixwidthx
           dx = xpix - x(i)
           q2 = (dx*dx + dy*dy)*hi1*hi1
           !
           !--SPH kernel
           !
           wab = w_cubic(q2)
           !
           !--calculate data value at this pixel using the summation interpolant
           !
           vecsmoothx(ipixi,jpixi) = vecsmoothx(ipixi,jpixi) + termx*wab
           vecsmoothy(ipixi,jpixi) = vecsmoothy(ipixi,jpixi) + termy*wab
           if (normalise) datnorm(ipixi,jpixi) = datnorm(ipixi,jpixi) + termnorm*wab

        enddo
     enddo

  enddo over_parts
  !
  !--normalise dat arrays
  !
  if (normalise) then
     where (datnorm > 0.)
        vecsmoothx = vecsmoothx/datnorm
        vecsmoothy = vecsmoothy/datnorm
     end where
  endif

  return

end subroutine interpolate2D_vec

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

subroutine interpolate2D_xsec(x,y,hh,weight,dat,itype,npart,&
     x1,y1,x2,y2,datsmooth,npixx,normalise)

  integer, intent(in) :: npart,npixx
  real, intent(in), dimension(npart) :: x,y,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: x1,y1,x2,y2
  real, intent(out), dimension(npixx) :: datsmooth
  logical, intent(in) :: normalise
  real, dimension(npixx) :: datnorm

  integer :: i,ipix,ipixmin,ipixmax
  real :: hi,hi1,radkern,q2,wab,const
  real :: term,termnorm,dx,dy,xpix,ypix,pixwidth,xpixwidth,xlength
  real :: gradient,yintercept,aa,bb,cc,determinant,det
  real :: xstart,xend,ystart,yend,rstart,rend
  real :: tol
  logical :: xsame, ysame, debug

  debug = .false.
  !
  !--check for errors in input
  !
  tol = 1.e-3
  ysame = (abs(y2 - y1).lt.tol)
  xsame = (abs(x2 - x1).lt.tol)
  if (xsame.and.ysame) then
     print*,'error: interpolate2D_xsec: zero length cross section'
     return
  endif
  if (npixx.eq.0) then
     print*,'error: interpolate2D_xsec: npix = 0 '
     return
  endif
  print*,'oblique 1D cross section through 2D data: npix =',npixx
  !
  !--work out the equation of the line y = mx + c from the two points input
  !
  gradient = 0.
  if (.not.xsame) gradient = (y2-y1)/(x2-x1)
  yintercept = y2 - gradient*x2
  print*,'cross section line: y = ',gradient,'x + ',yintercept
  !
  !--work out length of line and divide into pixels
  !
  xlength = sqrt((x2-x1)**2 + (y2-y1)**2)
  pixwidth = xlength/real(npixx)
  xpixwidth = (x2 - x1)/real(npixx)
  if (debug) then
     print*,'length of line = ',xlength
     print*,'pixel width = ',pixwidth, ' in x direction = ',xpixwidth
  endif
  !
  !--now interpolate to the line of pixels
  !
  datsmooth = 0.
  datnorm = 0.
  const = cnormk2D   ! normalisation constant
  !
  !--loop over particles
  !
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts
     !
     !--skip particles with zero weights
     !
     termnorm = const*weight(i)
     if (termnorm.le.0.) cycle over_parts
     !
     !--skip particles with wrong h's
     !
     hi = hh(i)
     if (hi.le.tiny(hi)) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi1 = 1./hi
     radkern = radkernel*hi    ! radius of the smoothing kernel
     term = termnorm*dat(i)
     !
     !--for each particle work out which pixels it contributes to
     !  to do this we need to work out the two points at which the line
     !  intersects the particles smoothing circle
     !  given by the equation (x-xi)^2 + (y-yi)^2 = (2h)^2.
     !  The x co-ordinates of these points are the solutions to a
     !  quadratic with co-efficients:

     aa = 1. + gradient**2
     bb = 2.*gradient*(yintercept - y(i)) - 2.*x(i)
     cc = x(i)**2 + y(i)**2 - 2.*yintercept*y(i) + yintercept**2 &
          - radkern**2
     !
     !--work out whether there are any real solutions and find them
     !
     determinant = bb**2 - 4.*aa*cc
     if (determinant < 0) then
        !!print*,' particle ',i,': does not contribute ',x(i),y(i)
     else
        det = sqrt(determinant)
        xstart = (-bb - det)/(2.*aa)
        xend =  (-bb + det)/(2.*aa)
        if (xstart.lt.x1) xstart = x1
        if (xstart.gt.x2) xstart = x2
        if (xend.gt.x2) xend = x2
        if (xend.lt.x1) xend = x1
        ystart = gradient*xstart + yintercept
        yend = gradient*xend + yintercept
        !
        !--work out position in terms of distance (no. of pixels) along the line
        !
        rstart = sqrt((xstart-x1)**2 + (ystart-y1)**2)
        rend = sqrt((xend-x1)**2 + (yend-y1)**2)

        ipixmin = int(rstart/pixwidth)
        ipixmax = int(rend/pixwidth) + 1

        if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
        if (ipixmax.lt.1) ipixmax = 1
        if (ipixmax.gt.npixx) ipixmax = npixx
        if (ipixmin.gt.npixx) ipixmax = npixx
               !
        !--loop over pixels, adding the contribution from this particle
        !
        !if (debug) print*,' particle ',i,': ',ipixmin,ipixmax,xstart,x(i),xend
        do ipix = ipixmin,ipixmax

           xpix = x1 + (ipix-0.5)*xpixwidth
           ypix = gradient*xpix + yintercept
           dy = ypix - y(i)
           dx = xpix - x(i)
           q2 = (dx*dx + dy*dy)*hi1*hi1
           !
           !--SPH kernel
           !
           wab = w_cubic(q2)
           !
           !--calculate data value at this pixel using the summation interpolant
           !
           datsmooth(ipix) = datsmooth(ipix) + term*wab
           if (normalise) datnorm(ipix) = datnorm(ipix) + termnorm*wab

        enddo

     endif

  enddo over_parts
  !
  !--normalise dat array
  !
  if (normalise) then
     where (datnorm > 0.)
        datsmooth = datsmooth/datnorm
     end where
  endif

  return

end subroutine interpolate2D_xsec

!--------------------------------------------------------------------------
!     subroutine to render particles onto a pixel array
!     at the maximum or minimum colour
!
!     Written by Daniel Price 21/7/2008
!--------------------------------------------------------------------------

subroutine interpolate_part(x,y,hh,npart,xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval,brightness)

  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,hh
  real, intent(in) :: xmin,ymin,pixwidth,datval
  real, intent(inout), dimension(npixx,npixy) :: datsmooth
  real, intent(inout), dimension(npixx,npixy), optional :: brightness
  integer :: i

  if (pixwidth.le.0.) then
     print "(1x,a)",'interpolate_part: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate_part: warning: ignoring some or all particles with h < 0'
  endif
  !
  !--loop over particles
  !
  if (present(brightness)) then
     do i=1,npart
        call interpolate_part1(x(i),y(i),hh(i),xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval)
     enddo
  else
     do i=1,npart
        call interpolate_part1(x(i),y(i),hh(i),xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval,brightness)
     enddo
  endif
  return

end subroutine interpolate_part

!--------------------------------------------------------------------------
!     subroutine to render a single particle onto a pixel array
!
!     Written by Daniel Price 21/7/2008
!--------------------------------------------------------------------------

subroutine interpolate_part1(xi,yi,hi,xmin,ymin,datsmooth,npixx,npixy,pixwidth,datval,brightness)

  real, intent(in) :: xi,yi,hi,xmin,ymin,pixwidth,datval
  integer, intent(in) :: npixx,npixy
  real, intent(inout), dimension(npixx,npixy) :: datsmooth
  real, intent(inout), dimension(npixx,npixy), optional :: brightness
  integer :: ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: radkern,radkern2,rab2
  real :: dx,dy2,xpix,ypix

  !
  !--skip particles with wrong h's
  !
  if (hi.le.tiny(hi)) return
  !
  !--set kernel related quantities
  !
  radkern = max(hi,2.*pixwidth)
  radkern2 = radkern*radkern  ! radius of the smoothing kernel
  !
  !--for each particle work out which pixels it contributes to
  !
  ipixmin = int((xi - radkern - xmin)/pixwidth)
  jpixmin = int((yi - radkern - ymin)/pixwidth)
  ipixmax = int((xi + radkern - xmin)/pixwidth) + 1
  jpixmax = int((yi + radkern - ymin)/pixwidth) + 1

  if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
  if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
  if (ipixmax.gt.npixx) ipixmax = npixx
  if (jpixmax.gt.npixy) jpixmax = npixy
  !
  !--loop over pixels, adding the contribution from this particle
  !
  do jpix = jpixmin,jpixmax
     ypix = ymin + (jpix-0.5)*pixwidth
     dy2 = (ypix - yi)**2
     do ipix = ipixmin,ipixmax
        xpix = xmin + (ipix-0.5)*pixwidth
        dx = xpix - xi
        rab2 = dx**2 + dy2
        !
        !--set data value at this pixel to maximum
        !
        if (rab2.lt.radkern2) then
           datsmooth(ipix,jpix) = datval
           if (present(brightness)) then
              brightness(ipix,jpix) = 1.0
           endif
        endif
     enddo
  enddo

  return
end subroutine interpolate_part1

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

  real, dimension(npixx,npixy) :: datnorm,datold
  real, dimension(npixx) :: dx2i

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax,its,itsmax
  real :: hi,hi1,radkernx,radkerny,q2,wab,const
  real :: term,termnorm,dy,xpix,ypix,ddx,ddy
  real :: xi,yi,pixwidthx,pixwidthy,dy2

  if (adaptive) then
     print "(1x,a)",'interpolating from particles to 2D pixels (adaptive)...'
  else
     print "(1x,a)",'interpolating from particles to 2D pixels...'
  endif

  if (adaptive) then
     itsmax = 3
  else
     itsmax = 1
  endif

  iterations: do its=1,itsmax
  datsmooth = 0.
  datnorm = 0.

  const = cnormk2D  ! normalisation constant
  !
  !--loop over particles
  !
  ddx = npixx/(xmax - xmin)
  ddy = npixy/(ymax - ymin)

  pixwidthx = 1. !/npixx
  pixwidthy = 1. !/npixy
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print "(1x,a)",'interpolate2D: error: pixel width <= 0'
     return
  endif

  !$omp parallel do default(none) &
  !$omp shared(npart,itype,x,y,xmin,ymin,ddx,ddy,its) &
  !$omp shared(datold,datsmooth,datnorm,npixx,npixy,const,radkernel,radkernel2,dat) &
  !$omp shared(pixwidthx,pixwidthy,normalise) &
  !$omp private(i,xi,yi,ipix,jpix,hi,hi1) &
  !$omp private(radkernx,radkerny,ipixmin,ipixmax,jpixmin,jpixmax) &
  !$omp private(dx2i,xpix,ypix,dy,dy2,q2,wab,term,termnorm)
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts

     !
     !--scale particle positions into viewport coordinates
     !
     xi = (x(i) - xmin)*ddx
     yi = (y(i) - ymin)*ddy
     hi = 1.0*pixwidthx  ! in units of pixel spacing

     ipix = int(xi)
     jpix = int(yi)
     if (its > 1 .and. ipix.ge.1 .and. ipix.le.npixx.and. jpix.ge.1 .and. jpix.le.npixy) then
        hi = min(1.5/sqrt(datold(ipix,jpix)),20.) !) !,1.)
     endif
     hi1 = 1./hi
     termnorm = const*hi1*hi1
     !
     !--set kernel related quantities
     !
     radkernx = radkernel*hi  ! radius of the smoothing kernel
     radkerny = radkernel*hi  ! radius of the smoothing kernel
     if (present(dat)) then
        term = termnorm*dat(i)
     else
        term = termnorm
     endif
     !
     !--for each particle work out which pixels it contributes to
     !
     ipixmin = int((xi - radkernx))
     jpixmin = int((yi - radkerny))
     ipixmax = int((xi + radkernx)) + 1
     jpixmax = int((yi + radkerny)) + 1

     if (ipixmin.lt.1)     ipixmin = 1      ! make sure they only contribute
     if (ipixmax.gt.npixx) ipixmax = npixx  ! to pixels in the image
     if (jpixmin.lt.1)     jpixmin = 1
     if (jpixmax.gt.npixy) jpixmax = npixy
     !
     !--precalculate an array of dx2 for this particle (optimisation)
     !
     do ipix=ipixmin,ipixmax
        xpix = (ipix-0.5)*pixwidthx
        dx2i(ipix) = ((xpix - xi)**2)*hi1*hi1
     enddo
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        ypix = (jpix-0.5)*pixwidthy
        dy = ypix - yi
        dy2 = dy*dy*hi1*hi1

        do ipix = ipixmin,ipixmax
           q2 = dx2i(ipix) + dy2
           !
           !--SPH kernel
           !
           if (q2 < radkernel2) then
              wab = w_cubic(q2)
              !
              !--calculate data value at this pixel using the summation interpolant
              !
              !$omp atomic
              datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab
              if (normalise) then
                 !$omp atomic
                 datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab
              endif
           endif

        enddo
     enddo

  enddo over_parts
  !$omp end parallel do

  if (present(dat)) then
     datold = datnorm
  else
     datold = datsmooth
  endif
  !
  !--normalise dat array
  !
  if (normalise) then
     where (datnorm > 0.)
        datsmooth = datsmooth/datnorm
     end where
  endif

  enddo iterations

  if (present(datpix2)) datpix2 = datnorm

  return

end subroutine interpolate2D_pixels

!---------------------------------------
!
!  Functional forms of various kernels
!
!--------------------------------------

pure real function w_cubic(q2)

 real, intent(in) :: q2
 real :: q

 if (q2.lt.1.0) then
    q = sqrt(q2)
    w_cubic = 1.-1.5*q2 + 0.75*q2*q
 elseif (q2.lt.4.0) then
    q = sqrt(q2)
    w_cubic = 0.25*(2.-q)**3
 else
    w_cubic = 0.
 endif

end function w_cubic

pure real function w_quartic(q2)

 real, intent(in) :: q2
 real :: q

 q = sqrt(q2)
 if (q.lt.0.5) then
    w_quartic = (2.5-q)**4 - 5.*(1.5-q)**4 + 10.*(0.5-q)**4
 elseif (q.lt.1.5) then
    w_quartic = (2.5-q)**4 - 5.*(1.5-q)**4
 elseif (q.lt.2.5) then
    w_quartic = (2.5-q)**4
 else
    w_quartic = 0.
 endif

end function w_quartic

pure real function w_quintic(q2)

 real, intent(in) :: q2
 real :: q,q4

 if (q2.lt.1.0) then
    q = sqrt(q2)
    q4 = q2*q2
    w_quintic = 66.-60.*q2 + 30.*q4 - 10.*q4*q
 elseif ((q2.ge.1.0).and.(q2.lt.4.0)) then
    q = sqrt(q2)
    w_quintic = (3.-q)**5 - 6.*(2.-q)**5
 elseif ((q2.ge.4.0).and.(q2.lt.9.0)) then
    q = sqrt(q2)
    w_quintic = (3.-q)**5
 else
    w_quintic = 0.0
 endif

end function w_quintic

pure real function w_quartic2h(q2)

 real, intent(in) :: q2
 real :: q

 q = sqrt(q2)
 if (q.lt.0.4) then
    w_quartic2h = (2.-q)**4 - 5.*(1.2-q)**4 + 10.*(0.4-q)**4
 elseif (q.lt.1.2) then
    w_quartic2h = (2.-q)**4 - 5.*(1.2-q)**4
 elseif (q.lt.2.) then
    w_quartic2h = (2.-q)**4
 else
    w_quartic2h = 0.
 endif

end function w_quartic2h

pure real function w_wendlandc2(q2)

 real, intent(in) :: q2
 real :: q

 if (q2.lt.4.) then
    q = sqrt(q2)
    w_wendlandc2 = (1. - 0.5*q)**4*(2.*q + 1.)
 else
    w_wendlandc2 = 0.
 endif

end function w_wendlandc2

pure real function w_wendlandc4(q2)

 real, intent(in) :: q2
 real :: q

 if (q2.lt.4.) then
    q = sqrt(q2)
    w_wendlandc4 = (1. - 0.5*q)**6*(35./12.*q2 + 3.*q + 1.)
 else
    w_wendlandc4 = 0.
 endif

end function w_wendlandc4

pure real function w_wendlandc6(q2)

 real, intent(in) :: q2
 real :: q

 if (q2.lt.4.) then
    q = sqrt(q2)
    w_wendlandc6 = (1. - 0.5*q)**8*(4.*q2*q + 25./4.*q2 + 4.*q + 1.)
 else
    w_wendlandc6 = 0.
 endif

end function w_wendlandc6

!---------------------------------------
!
!  Integrals of the kernel, used
!  when performing exact interpolation
!  See Petkova, Laibe & Bonnell (2018)
!
!---------------------------------------
pure real function pint(r0, d1, d2, hi1)
 real, intent(in) :: r0, d1, d2, hi1
 real :: ar0,q0,tphi1,tphi2,phi1,phi2

 if (abs(r0) < tiny(0.)) then
    pint = 0.
    return
 elseif (r0 > 0.) then
    pint = 1.
    ar0 = r0
 else
    pint = -1.
    ar0 = -r0
 endif
 q0 = ar0*hi1
 tphi1 = abs(d1)/ar0
 tphi2 = abs(d2)/ar0
 phi1 = atan(tphi1)
 phi2 = atan(tphi2)

 if (d1*d2 >= 0.) then
    pint = pint*(full_2d_mod(phi1,tphi1,q0) + full_2d_mod(phi2,tphi2,q0))
 else
    if (abs(d1) < abs(d2)) then
       pint = pint*(full_2d_mod(phi2,tphi2,q0) - full_2d_mod(phi1,tphi1,q0))
    else
       pint = pint*(full_2d_mod(phi1,tphi1,q0) - full_2d_mod(phi2,tphi2,q0))
    endif
 endif

end function pint

!---------------------------------------
!
! Helper functions for kernel integrals
! See Petkova, Laibe & Bonnell (2018)
!
!---------------------------------------
pure real function full_2d_mod(phi,tphi,q0)
 real, intent(in) :: phi,tphi,q0
 real :: q, phi1, phi2, tphi1, tphi2, cphi, cphi1, cphi2

 if (q0 <= 1.0) then
    cphi = cos(phi)
    q = q0/cphi

    if (q <= 1.0) then
       full_2d_mod = F1_2d(phi,tphi,cphi,q0)
    elseif (q <= 2.0) then
       cphi1 = q0
       phi1 = acos(q0)
       tphi1 = tan(phi1)
       full_2d_mod = F2_2d(phi,tphi,cphi,q0) - F2_2d(phi1,tphi1,cphi1,q0) &
                   + F1_2d(phi1,tphi1,cphi1,q0)
    else
       cphi1 = q0
       phi1 = acos(q0)
       cphi2 = 0.5*q0
       phi2 = acos(0.5*q0)
       tphi1 = tan(phi1)
       tphi2 = tan(phi2)
       full_2d_mod = F3_2d(phi) - F3_2d(phi2) + F2_2d(phi2,tphi2,cphi2,q0) &
                                - F2_2d(phi1,tphi1,cphi1,q0) &
                                + F1_2d(phi1,tphi1,cphi1,q0)
    endif
 elseif (q0 <= 2.0) then
    cphi = cos(phi)
    q = q0/cphi

    if (q <= 2.0) then
       full_2d_mod = F2_2d(phi,tphi,cphi,q0)
    else
       cphi2 = 0.5*q0
       phi2 = acos(0.5*q0)
       tphi2 = tan(phi2)
       full_2d_mod = F3_2d(phi) - F3_2d(phi2) + F2_2d(phi2,tphi2,cphi2,q0)
    endif
 else ! q0 > 2
    full_2d_mod = F3_2d(phi)
 endif

end function full_2d_mod

pure real function F1_2d(phi,tphi,cphi,q0)
 real, intent(in) :: phi,tphi,cphi,q0
 real :: I2, I4, I5, logs, cphi2, q02, q03

 !tphi = tan(phi)
 !cphi = cos(phi)
 cphi2 = cphi*cphi

 q02 = q0*q0
 q03 = q02*q0

 logs = log(tan(phi/2.+pi/4.))

 I2 = tphi
 I4 = 1./3. * tphi*(2. + 1./cphi2)

 I5 = 1./16. * (0.5*(11.*sin(phi) + 3.*sin(3.*phi))/cphi2/cphi2 + 6.*logs)

 F1_2d =  5./7.*q02/pi  * (I2 - 3./4.*q02 *I4 + 0.3*q03 *I5)

end function F1_2d

pure real function F2_2d(phi, tphi, cphi, q0)
 real, intent(in) :: phi, tphi, cphi, q0
 real :: I0, I2, I3, I4, I5, logs, cphi2, q02, q03

! tphi = tan(phi)
! cphi = cos(phi)
 cphi2 = cphi*cphi

 q02 = q0*q0
 q03 = q02*q0

 logs = log(tan(phi/2.+pi/4.))

 I0 = phi
 I2 = tphi
 I4 = 1./3. * tphi*(2. + 1./(cphi2))

 I3 = 1./2.  * (tphi/cphi + logs)
 I5 = 1./16. * (0.5*(11.*sin(phi) + 3.*sin(3.*phi))/cphi2/cphi2 + 6.*logs)

 F2_2d =  5./7.*q02/pi  * (2.*I2 - 2.*q0 *I3 + 3./4.*q02 *I4 - 1./10.*q03 *I5 &
                           - 1./10./q02 *I0)

end function F2_2d

pure real function F3_2d(phi)
 real, intent(in) :: phi
 real :: I0

 I0 = phi
 F3_2d = 0.5/pi  *I0

end function F3_2d

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
  integer         :: i2,i1,itype,ipart,isinktype
  real :: dunitspmass,dunitsrho,dunitsh

  !
  !-- unit_interp is a multiplication factor that can
  !   be used to scale the "weight" in case the
  !   units read from the read_data routine
  !   are inconsistent (this is the case for the SEREN read)
  !
  dunitspmass = 1.d0
  dunitsrho   = 1.d0
  dunitsh     = 1.d0
  if (iRescale) then
     if (ipmass.gt.0) dunitspmass = 1.d0/units(ipmass)
     if (ih.gt.0)     dunitsh     = 1.d0/units(ih)
     if (irho.gt.0)   dunitsrho   = 1.d0/units(irho)
  endif
  dunitspmass = dunitspmass * unit_interp

  isinktype = get_sink_type(ntypes)

  if (ipmass.gt.0 .and. ipmass.le.ndataplots .and. &
      irho.gt.0 .and. irho.le.ndataplots .and. &
      ih .gt. 0 .and. ih.le.ndataplots .and. &
      required(ipmass) .and. required(irho) .and. required(ih)) then

     if (size(iamtypei) > 1) then
        !
        !--particles with mixed types
        !
        !$omp parallel do default(none) &
        !$omp shared(ninterp,iamtypei,weighti,dati,rendersinks,isinktype) &
        !$omp shared(usetype,idensityweighted,dunitsrho)    &
        !$omp shared(ipmass,ih,irho,dunitspmass,dunitsh,ndim) &
        !$omp private(ipart,itype)
        do ipart=1,ninterp
           itype = iamtypei(ipart)
           if (.not.usetype(itype)) then
              if (rendersinks .and. itype.eq.isinktype) then
                 weighti(ipart) = weight_sink
              else
                 weighti(ipart) = 0.
              endif
           elseif (idensityweighted) then
              if (dati(ipart,ih) > tiny(dati)) then
                 weighti(ipart) = (dati(ipart,ipmass)*dunitspmass)/ &
                                 ((dati(ipart,ih)*dunitsh)**ndim)
              else
                 weighti(ipart) = 0.
              endif
           else
              if (dati(ipart,irho) > tiny(dati) .and. dati(ipart,ih) > tiny(dati)) then
                 weighti(ipart) = (dati(ipart,ipmass)*dunitspmass)/ &
                                 ((dati(ipart,irho)*dunitsrho)*(dati(ipart,ih)*dunitsh)**ndim)
              else
                 weighti(ipart) = 0.
              endif
           endif
        enddo
        !$omp end parallel do
     else
        !
        !--particles ordered by type
        !
        i2 = 0
        over_types: do itype=1,ntypes
           i1 = i2 + 1
           i2 = i2 + npartoftype(itype)
           i2 = min(i2,ninterp)
           if (i1 > i2) cycle over_types
           !--set weights to zero for particle types not used in the rendering
           if (.not.usetype(itype)) then
              if (rendersinks .and. itype.eq.isinktype) then
                 weighti(i1:i2) = weight_sink
              else
                 weighti(i1:i2) = 0.
              endif
           elseif (idensityweighted) then
           !--for density weighted interpolation use m/h**ndim
              where(dati(i1:i2,ih) > tiny(dati))
                 weighti(i1:i2) = (dati(i1:i2,ipmass)*dunitspmass)/ &
                                  ((dati(i1:i2,ih)*dunitsh)**ndim)
              elsewhere
                 weighti(i1:i2) = 0.
              endwhere
           else
           !--usual interpolation use m/(rho h**ndim)
              where(dati(i1:i2,irho) > tiny(dati) .and. dati(i1:i2,ih) > tiny(dati))
                 weighti(i1:i2) = (dati(i1:i2,ipmass)*dunitspmass)/ &
                                  ((dati(i1:i2,irho)*dunitsrho)*(dati(i1:i2,ih)*dunitsh)**ndim)
              elsewhere
                 weighti(i1:i2) = 0.
              endwhere
           endif
        enddo over_types
     endif

     if (idensityweighted) then
        print "(a)",' USING DENSITY WEIGHTED INTERPOLATION '
        inormalise = .true.
     endif

  elseif (any(masstype(1:ntypes).gt.0.) .and. &
          irho.gt.0 .and. irho.le.ndataplots .and. &
          ih .gt. 0 .and. ih.le.ndataplots .and. &
          required(irho) .and. required(ih)) then

     if (size(iamtypei) > 1) then
        !
        !--particles with mixed types
        !
        !$omp parallel do default(none) &
        !$omp shared(ninterp,iamtypei,weighti,dati,rendersinks,isinktype) &
        !$omp shared(usetype,idensityweighted,dunitsrho,masstype) &
        !$omp shared(ih,irho,dunitspmass,dunitsh,ndim)    &
        !$omp private(ipart,itype)
        do ipart=1,ninterp
           itype = iamtypei(ipart)
           if (.not.usetype(itype)) then
              if (rendersinks .and. itype.eq.isinktype) then
                 weighti(ipart) = weight_sink
              else
                 weighti(ipart) = 0.
              endif
           elseif (idensityweighted) then
              if (dati(ipart,ih) > tiny(dati)) then
                 weighti(ipart) = (masstype(itype)*dunitspmass)/ &
                                 ((dati(ipart,ih)*dunitsh)**ndim)
              else
                 weighti(ipart) = 0.
              endif
           else
              if (dati(ipart,irho) > tiny(dati) .and. dati(ipart,ih) > tiny(dati)) then
                 weighti(ipart) = (masstype(itype)*dunitspmass)/ &
                                 ((dati(ipart,irho)*dunitsrho)*(dati(ipart,ih)*dunitsh)**ndim)
              else
                 weighti(ipart) = 0.
              endif
           endif
        enddo
        !$omp end parallel do
     else
        !
        !--particles ordered by type
        !
        i2 = 0
        over_types2: do itype=1,ntypes
           i1 = i2 + 1
           i2 = i2 + npartoftype(itype)
           i2 = min(i2,ninterp)
           if (i1 > i2) cycle over_types2
           !--set weights to zero for particle types not used in the rendering
           if (.not.usetype(itype)) then
              if (rendersinks .and. itype.eq.isinktype) then
                 weighti(i1:i2) = weight_sink
              else
                 weighti(i1:i2) = 0.
              endif
           else
              where(dati(i1:i2,irho) > tiny(dati) .and. dati(i1:i2,ih) > tiny(dati))
                 weighti(i1:i2) = masstype(itype)/ &
                                ((dati(i1:i2,irho)*dunitsrho)*(dati(i1:i2,ih)*dunitsh)**ndim)
              elsewhere
                 weighti(i1:i2) = 0.
              endwhere
           endif
        enddo over_types2
     endif

     if (idensityweighted) then
        print "(a)",' USING DENSITY WEIGHTED INTERPOLATION '
        inormalise = .true.
     endif
  else
     if (required(ih) .and. required(irho) .and. ih.gt.0 .and. irho.gt.0) then
        print "(a)",' WARNING: particle mass not set: using normalised interpolations'
     endif
     weighti(1:ninterp) = 1.0

  !--if particle mass has not been set, then must use normalised interpolations
     inormalise = .true.

     if (size(iamtypei) > 1) then
        !
        !--particles with mixed types
        !
        !$omp parallel do default(none) &
        !$omp shared(ninterp,iamtypei,weighti,usetype) &
        !$omp private(ipart,itype)
        do ipart=1,ninterp
           itype = iamtypei(ipart)
           if (.not.usetype(itype)) weighti(ipart) = 0.
        enddo
        !$omp end parallel do
     endif
  endif

end subroutine set_interpolation_weights

!-----------------------------------------------------------------
!
! utility to "guess" which particle type contains sink particles
! from the label
!
!-----------------------------------------------------------------

integer function get_sink_type(ntypes)

 integer, intent(in) :: ntypes
 integer :: i

 get_sink_type = 0
 do i=1,ntypes
    if (get_sink_type.eq.0 .and. index(labeltype(i),'sink').ne.0) get_sink_type = i
 enddo

end function get_sink_type

end module splash
