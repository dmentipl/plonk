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
!  Copyright (C) 2005-2017 Daniel Price. All rights reserved.
!  Contact: daniel.price@monash.edu
!
!  CHANGES made by Daniel Mentiplay, starting on 25-01-2019.
!
!-----------------------------------------------------------------

module splash
 implicit none

 public :: interpolate3D_projection
 public :: interpolate3D_proj_vec
 public :: interpolate3D_fastxsec
 public :: interpolate3D_xsec_vec
 public :: interp3D_proj_opacity

 private

 integer, parameter :: maxcoltable = 1000
 real :: coltable(maxcoltable)
 real :: dq2table = 4./maxcoltable
 real :: ddq2table = maxcoltable/4.

 real :: radkernel   = 2.
 real :: radkernel2  = 4.
 real :: cnormk3D    = 0.3183098862
 real :: weight_sink = -1

contains

!-------------------------------------------------------------
! tabulates the integral through the cubic spline kernel
! tabulated in (r/h)**2 so that sqrt is not necessary
!-------------------------------------------------------------
subroutine setup_integratedkernel

 integer :: i,j
 real :: rxy2,deltaz,dz,z,q2,wkern,coldens
 integer, parameter :: npts = 100

 dq2table = radkernel2/maxcoltable
 ddq2table = 1./dq2table

 do i=0,maxcoltable-1
!
!--tabulate for (cylindrical) r**2 between 0 and radkernel**2
!
    rxy2 = i*dq2table
!
!--integrate z between 0 and sqrt(radkernel^2 - rxy^2)
!
    deltaz = sqrt(radkernel2 - rxy2)
    dz = deltaz/real(npts-1)
    coldens = 0.
    if (deltaz.ne.deltaz) print "(a)",'WARNING: NaN in kernel table setup'
    do j=1,npts
       z = (j-1)*dz
       q2 = rxy2 + z*z
#ifdef KERNEL
#if KERNEL==1
       wkern = w_cubic(q2)
#elif KERNEL==2
       wkern = w_quartic(q2)
#elif KERNEL==3
       wkern = w_quintic(q2)
#elif KERNEL==4
       wkern = w_quartic2h(q2)
#elif KERNEL==5
       wkern = w_wendlandc2(q2)
#elif KERNEL==6
       wkern = w_wendlandc4(q2)
#elif KERNEL==7
       wkern = w_wendlandc6(q2)
#else
       wkern = w_cubic(q2)
#endif
#else
       wkern = w_cubic(q2)
#endif
       if (j.eq.1 .or. j.eq.npts) then
          coldens = coldens + 0.5*wkern*dz ! trapezoidal rule
       else
          coldens = coldens + wkern*dz
       endif
    enddo
    coltable(i)=2.0*coldens*cnormk3D
 end do
 coltable(maxcoltable) = 0.

 return
end subroutine setup_integratedkernel

!-------------------------------------------------------------
! This function interpolates from the table of integrated kernel values
! to give w(q)
!-------------------------------------------------------------
real function wfromtable(q2)

 real, intent(in) :: q2
 real :: dxx,dwdx
 integer :: index, index1
 !
 !--find nearest index in table
 !
 index = max(int(q2*ddq2table),0) ! the max prevents seg faults on NaNs for q2
 !index = min(index,maxcoltable) ! should be unnecessary if q2 < radkernel checked
 index1 = min(index + 1,maxcoltable)
 !
 !--find increment along from this index
 !
 dxx = q2 - index*dq2table
 !
 !--find gradient
 !
 dwdx = (coltable(index1) - coltable(index))*ddq2table
 !
 !--compute value of integrated kernel
 !
 wfromtable = coltable(index) + dwdx*dxx

end function wfromtable

!--------------------------------------------------------------------------
!     subroutine to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b weight_b dat_b W(r-r_b, h_b)
!
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     ** In this version 3D data is interpolated to a 2D grid by use of an
!     ** integrated form of the kernel (that is W_ab in this case is
!     ** the integral through the 3D kernel to give a 2D kernel)
!     ** This results in a column density map of the interpolated quantity
!     ** From a similar routine by Matthew Bate.
!
!     The (dimensionless) weight for each particle should be
!
!     weight = pmass/(rho*h^3)
!
!     the interface is written in this form to avoid floating exceptions
!     on physically scaled data.
!
!     Input: particle coordinates  : x,y,z (npart) - note that z is only required for perspective
!            smoothing lengths     : hh    (npart)
!            weight for each particle : weight (npart)
!            scalar data to smooth : dat   (npart)
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Written by Daniel Price September 2003
!     3D perspective added Nov 2005
!--------------------------------------------------------------------------
subroutine interpolate3D_projection(x,y,z,hh,weight,dat,itype,npart, &
     xmin,ymin,datsmooth,npixx,npixy,pixwidthx,pixwidthy,normalise,zobserver,dscreen, &
     useaccelerate)

  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy,zobserver,dscreen
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise
  real, dimension(npixx,npixy) :: datnorm
  logical, intent(in) :: useaccelerate
  real :: row(npixx)

  integer :: ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax,npixpartx,npixparty
  integer :: iprintinterval, iprintnext,ipixi,jpixi,jpixcopy
  integer :: nsubgrid,nfull,nok,ncpus
#ifdef _OPENMP
  integer :: omp_get_num_threads,i
#else
  integer(kind=selected_int_kind(10)) :: iprogress,i  ! up to 10 digits
#endif
  real :: hi,hi1,hi21,radkern,wab,q2,xi,yi,xminpix,yminpix
  real :: term,termnorm,dy,dy2,ypix,zfrac,hsmooth,horigi
  real :: xpixmin,xpixmax,xmax,ypixmin,ypixmax,ymax
  real :: hmin,fac,hminall !,dhmin3
  real, dimension(npixx) :: xpix,dx2i
  logical :: iprintprogress,use3Dperspective,accelerate
  character(len=32) :: string

  datsmooth = 0.
  term = 0.
  string = 'projecting'
  if (normalise) then
     string = trim(string)//' (normalised)'
     datnorm = 0.
  elseif (useaccelerate) then
     string = trim(string)//' (fast)'
  endif

  ncpus = 0
  !$omp parallel
  !$omp master
  !$ ncpus = omp_get_num_threads()
  !$omp end master
  !$omp end parallel

  if (ncpus > 0) then
     write (*,"(1x,a,': ',i4,' x ',i4,' on ',i3,' cpus')") trim(string),npixx,npixy,ncpus
  else
     write (*,"(1x,a,': ',i4,' x ',i4)") trim(string),npixx,npixy
  endif

  if (pixwidthx.le.0. .or. pixwidthy.le.0) then
     print "(1x,a)",'ERROR: pixel width <= 0'
     return
  endif
  !nout = count(hh(1:npart).le.0.)
  !if (nout.gt.0) then
  !   print*,'interpolate3D_projection: warning: ignoring ',nout,' particles with h <= 0'
  !endif
  !
  !--check column density table has actually been setup
  !
  if (abs(coltable(1)).le.1.e-5) then
     call setup_integratedkernel
  endif
  !
  !--print a progress report if it is going to take a long time
  !  (a "long time" is, however, somewhat system dependent)
  !
  iprintprogress = (npart .ge. 100000) .or. (npixx*npixy .gt.100000)
  !
  !--loop over particles
  !
  iprintinterval = 25
  if (npart.ge.1e6) iprintinterval = 10
  iprintnext = iprintinterval
  use3Dperspective = abs(dscreen).gt.tiny(dscreen)

  xminpix = xmin - 0.5*pixwidthx
  yminpix = ymin - 0.5*pixwidthy
  xmax = xmin + npixx*pixwidthx
  ymax = ymin + npixy*pixwidthy
!
!--use a minimum smoothing length on the grid to make
!  sure that particles contribute to at least one pixel
!
  hmin = 0.5*max(pixwidthx,pixwidthy)
  !dhmin3 = 1./(hmin*hmin*hmin)
!
!--store x value for each pixel (for optimisation)
!
  do ipix=1,npixx
     xpix(ipix) = xminpix + ipix*pixwidthx
  enddo
  nsubgrid = 0
  nok = 0
  hminall = huge(hminall)

!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,dat,itype,npart) &
!$omp shared(xmin,ymin,xmax,ymax,xminpix,yminpix,xpix,pixwidthx,pixwidthy) &
!$omp shared(npixx,npixy,dscreen,zobserver,use3dperspective,useaccelerate) &
!$omp shared(normalise,radkernel,radkernel2,datsmooth,datnorm) &
!$omp firstprivate(hmin) & !,dhmin3) &
!$omp private(hi,zfrac,xi,yi,radkern,xpixmin,xpixmax,ypixmin,ypixmax) &
!$omp private(hsmooth,horigi,hi1,hi21,term,termnorm,npixpartx,npixparty,jpixi,ipixi) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax,accelerate) &
!$omp private(dx2i,row,q2,ypix,dy,dy2,wab) &
!$omp private(i,ipix,jpix,jpixcopy,fac) &
!$omp reduction(+:nsubgrid,nok) &
!$omp reduction(min:hminall)
!$omp do schedule (guided, 2)
  over_particles: do i=1,npart
     !
     !--report on progress
     !
#ifndef _OPENMP
     if (iprintprogress) then
        iprogress = 100*i/npart
        if (iprogress.ge.iprintnext) then
           write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,i
           iprintnext = iprintnext + iprintinterval
        endif
     endif
#endif
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_particles
     !
     !--set h related quantities
     !
     hi = hh(i)
     horigi = hi
     if (hi.le.0.) then
        cycle over_particles
     elseif (use3Dperspective) then
        if (z(i).gt.zobserver) cycle over_particles
        zfrac = abs(dscreen/(z(i)-zobserver))
        hi = hi*zfrac
     endif

     radkern = radkernel*hi ! radius of the smoothing kernel

     !--cycle as soon as we know the particle does not contribute
     xi = x(i)
     xpixmin = xi - radkern
     if (xpixmin.gt.xmax) cycle over_particles
     xpixmax = xi + radkern
     if (xpixmax.lt.xmin) cycle over_particles

     yi = y(i)
     ypixmin = yi - radkern
     if (ypixmin.gt.ymax) cycle over_particles
     ypixmax = yi + radkern
     if (ypixmax.lt.ymin) cycle over_particles

     !--take resolution length as max of h and 1/2 pixel width
     if (hi.lt.hmin) then
        hminall = min(hi,hminall)
        nsubgrid = nsubgrid + 1
        hsmooth = hmin
        fac = 1. !(horigi*horigi*horigi)*dhmin3      ! factor by which to adjust the weight
     else
        fac = 1.
        hsmooth = hi
        nok = nok + 1
     endif
     radkern = radkernel*hsmooth
     !
     !--set kernel related quantities
     !
     hi1 = 1./hsmooth
     hi21 = hi1*hi1
     termnorm = weight(i)*fac*horigi
     term = termnorm*dat(i) ! h gives the z length scale (NB: no perspective)
     !
     !--for each particle work out which pixels it contributes to
     !
     npixpartx = int(radkern/pixwidthx) + 1
     npixparty = int(radkern/pixwidthy) + 1
     jpixi = int((yi-ymin)/pixwidthy) + 1
     ipixi = int((xi-xmin)/pixwidthx) + 1
     ipixmin = ipixi - npixpartx
     ipixmax = ipixi + npixpartx
     jpixmin = jpixi - npixparty
     jpixmax = jpixi + npixparty

!     ipixmin = int((xi - radkern - xmin)/pixwidth)
!     jpixmin = int((yi - radkern - ymin)/pixwidth)
!     ipixmax = ipixmin + npixpart !!int((xi + radkern - xmin)/pixwidth) + 1
!     jpixmax = jpixmin + npixpart !!int((yi + radkern - ymin)/pixwidth) + 1
     !
     !--loop over pixels, adding the contribution from this particle
     !  copy by quarters if all pixels within domain
     !
     accelerate = useaccelerate .and. (.not.normalise) &
                 .and. npixpartx.gt.5 .and. npixparty.gt.5 &
                 .and. ipixmin.ge.1 .and. ipixmax.le.npixx &
                 .and. jpixmin.ge.1 .and. jpixmax.le.npixy

     if (accelerate) then
        !--adjust xi, yi to centre of pixel
        xi = xminpix + ipixi*pixwidthx
        yi = yminpix + jpixi*pixwidthy
        !
        !--precalculate an array of dx2 for this particle (optimisation)
        !
        do ipix=ipixmin,ipixmax
           dx2i(ipix) = ((xpix(ipix) - xi)**2)*hi21
        enddo
        do jpix = jpixi,jpixmax
           ypix = yminpix + jpix*pixwidthy
           dy = ypix - yi
           dy2 = dy*dy*hi21
           do ipix = ipixi,ipixmax
              q2 = dx2i(ipix) + dy2
              !
              !--SPH kernel - integral through cubic spline
              !  interpolate from a pre-calculated table
              !
              if (q2.lt.radkernel2) then
                 wab = wfromtable(q2)
                 !
                 !--calculate data value at this pixel using the summation interpolant
                 !
                 !$omp atomic
                 datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab
                 row(ipix) = term*wab
              else
                 row(ipix) = 0.
              endif
           enddo
           !--NB: the following actions can and should be vectorized (but I don't know how...)
           !--copy top right -> top left
           do ipix=ipixmin,ipixi-1
              !$omp atomic
              datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + row(ipixmax-(ipix-ipixmin))
           enddo
           if (jpix.ne.jpixi) then
              jpixcopy = jpixi - (jpix-jpixi)
              !--copy top right -> bottom left
              do ipix=ipixmin,ipixi-1
                 !$omp atomic
                 datsmooth(ipix,jpixcopy) = datsmooth(ipix,jpixcopy) + row(ipixmax-(ipix-ipixmin))
              enddo
              !--copy top right -> bottom right
              do ipix=ipixi,ipixmax
                 !$omp atomic
                 datsmooth(ipix,jpixcopy) = datsmooth(ipix,jpixcopy) + row(ipix)
              enddo
           endif
        enddo

     else
        ipixmin = int((xi - radkern - xmin)/pixwidthx)
        ipixmax = int((xi + radkern - xmin)/pixwidthx)
        jpixmin = int((yi - radkern - ymin)/pixwidthy)
        jpixmax = int((yi + radkern - ymin)/pixwidthy)

        if (ipixmin.lt.1) ipixmin = 1  ! make sure they only contribute
        if (jpixmin.lt.1) jpixmin = 1  ! to pixels in the image
        if (ipixmax.gt.npixx) ipixmax = npixx ! (note that this optimises
        if (jpixmax.gt.npixy) jpixmax = npixy !  much better than using min/max)
        !
        !--precalculate an array of dx2 for this particle (optimisation)
        !
        do ipix=ipixmin,ipixmax
           dx2i(ipix) = ((xpix(ipix) - xi)**2)*hi21
        enddo

        do jpix = jpixmin,jpixmax
           ypix = yminpix + jpix*pixwidthy
           dy = ypix - yi
           dy2 = dy*dy*hi21
           do ipix = ipixmin,ipixmax
              !xpix = xminpix + ipix*pixwidthx
              !dx = xpix - xi
              !rab2 = (xminpix + ipix*pixwidthx - xi)**2 + dy2
              q2 = dx2i(ipix) + dy2 ! dx2 pre-calculated; dy2 pre-multiplied by hi21
              !
              !--SPH kernel - integral through cubic spline
              !  interpolate from a pre-calculated table
              !
              if (q2.lt.radkernel2) then
                 wab = wfromtable(q2)
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
     endif

  enddo over_particles
!$omp end do
!$omp end parallel
!
!--normalise dat array
!
  if (normalise) then
     !--normalise everywhere (required if not using SPH weighting)
     where (datnorm > tiny(datnorm))
        datsmooth = datsmooth/datnorm
     end where
  endif
!
!--warn about subgrid interpolation
!
  if (nsubgrid.gt.1) then
     nfull = int((xmax-xmin)/(hminall)) + 1
     if (nsubgrid.gt.0.1*nok) &
     print "(a,i9,a,/,a,i6,a)",' Warning: pixel size > 2h for ',nsubgrid,' particles', &
                               '          need',nfull,' pixels for full resolution'
  endif

  return

end subroutine interpolate3D_projection

!--------------------------------------------------------------------------
!
!     Same as previous but for a vector quantity
!
!     Input: particle coordinates  : x,y   (npart)
!            smoothing lengths     : hh    (npart)
!            weight for each particle : weight (npart)
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!
!     Output: smoothed vector field   : vecsmoothx (npixx,npixy)
!                                     : vecsmoothy (npixx,npixy)
!
!     Daniel Price 23/12/04
!--------------------------------------------------------------------------
subroutine interpolate3D_proj_vec(x,y,z,hh,weight,vecx,vecy,itype,npart,&
     xmin,ymin,vecsmoothx,vecsmoothy,npixx,npixy,pixwidthx,pixwidthy,normalise,zobserver,dscreen)

  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,vecx,vecy
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy,zobserver,dscreen
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx, vecsmoothy
  logical, intent(in) :: normalise
  real, dimension(:,:), allocatable :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax,ierr
  real :: hi,hi1,hi21,radkern,q2,wab,rab2,const,zfrac,hsmooth
  real :: termx,termy,termnorm,dx,dy,dy2,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  termx = 0.
  termy = 0.
  if (normalise) then
     print "(1x,a)",'projecting vector (normalised) from particles to pixels...'
     allocate(datnorm(npixx,npixy),stat=ierr)
     if (ierr /= 0) then
        print "(a)",'interpolate3D_proj_vec: error allocating memory'
        return
     endif
     datnorm = 0.
  else
     print "(1x,a)",'projecting vector from particles to pixels...'
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print "(a)",'interpolate3D_proj_vec: error: pixel width <= 0'
     return
  endif
  !
  !--loop over particles
  !
!$omp parallel default(none) &
!$omp shared(hh,z,x,y,weight,vecx,vecy,itype,vecsmoothx,vecsmoothy,npart) &
!$omp shared(xmin,ymin,pixwidthx,pixwidthy,zobserver,dscreen,datnorm) &
!$omp shared(npixx,npixy,normalise,radkernel,radkernel2) &
!$omp private(hi,radkern,const,zfrac,ypix,xpix) &
!$omp private(hsmooth,hi1,hi21,termx,termy,termnorm) &
!$omp private(ipixmin,ipixmax,jpixmin,jpixmax) &
!$omp private(dy,dy2,dx,rab2,q2,wab) &
!$omp private(i,ipix,jpix)

!$omp do schedule(guided, 2)
  over_particles: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_particles
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     const = weight(i)*hi ! h gives the z length scale (NB: no perspective)
     if (hi.le.0.) then
        cycle over_particles
     elseif (abs(dscreen).gt.tiny(dscreen)) then
        if (z(i).gt.zobserver) cycle over_particles
        zfrac = abs(dscreen/(z(i)-zobserver))
        hi = hi*zfrac
     endif

     !--take resolution length as max of h and 1/2 pixel width
     hsmooth = max(hi,0.5*min(pixwidthx,pixwidthy))

     radkern = radkernel*hsmooth    ! radius of the smoothing kernel
     hi1 = 1./hsmooth
     hi21 = hi1*hi1

     termx = const*vecx(i)
     termy = const*vecy(i)
     termnorm = const
     !
     !--for each particle work out which pixels it contributes to
     !
     ipixmin = int((x(i) - radkern - xmin)/pixwidthx)
     jpixmin = int((y(i) - radkern - ymin)/pixwidthy)
     ipixmax = int((x(i) + radkern - xmin)/pixwidthx) + 1
     jpixmax = int((y(i) + radkern - ymin)/pixwidthy) + 1

     ! PRINT*,'particle ',i,' x, y, z = ',x(i),y(i),z(i),dat(i),rho(i),hi
     ! PRINT*,'pixels = ',ipixmin,ipixmax,jpixmin,jpixmax

     if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
     if (ipixmax.gt.npixx) ipixmax = npixx
     if (jpixmax.gt.npixy) jpixmax = npixy
     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        ypix = ymin + (jpix-0.5)*pixwidthy
        dy = ypix - y(i)
        dy2 = dy*dy
        do ipix = ipixmin,ipixmax
           xpix = xmin + (ipix-0.5)*pixwidthx
           dx = xpix - x(i)
           rab2 = dx**2 + dy2
           q2 = rab2*hi21
           !
           !--SPH kernel - integral through cubic spline
           !  interpolate from a pre-calculated table
           !
           if (q2.lt.radkernel2) then
              wab = wfromtable(q2)
              !
              !--calculate data value at this pixel using the summation interpolant
              !
              vecsmoothx(ipix,jpix) = vecsmoothx(ipix,jpix) + termx*wab
              vecsmoothy(ipix,jpix) = vecsmoothy(ipix,jpix) + termy*wab
              if (normalise) datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab
           endif

        enddo
     enddo

  enddo over_particles
!$omp end do
!$omp end parallel

  if (normalise .and. allocated(datnorm)) then
     !--normalise everywhere
     where (datnorm > tiny(datnorm))
        vecsmoothx = vecsmoothx/datnorm
        vecsmoothy = vecsmoothy/datnorm
     end where
  endif
  if (allocated(datnorm)) deallocate(datnorm)

  return

end subroutine interpolate3D_proj_vec

!--------------------------------------------------------------------------
!
!     ** In this version 3D data is interpolated to a single 2D cross section
!     ** This is much faster than interpolating to a 3D grid
!     ** and is efficient if only one or two cross sections are needed.
!
!     ** Note that the cross section is always taken in the z co-ordinate
!     ** so should submit the appropriate arrays as x, y and z.
!
!     Input: particle coordinates  : x,y,z (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            scalar data to smooth : dat   (npart)
!            cross section location: zslice
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!
!     Daniel Price, Institute of Astronomy, Cambridge, 23/9/03
!--------------------------------------------------------------------------
subroutine interpolate3D_fastxsec(x,y,z,hh,weight,dat,itype,npart,&
     xmin,ymin,zslice,datsmooth,npixx,npixy,pixwidthx,pixwidthy,normalise)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy,zslice
  real, intent(out), dimension(npixx,npixy) :: datsmooth
  logical, intent(in) :: normalise
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,radkern,q2,wab,const,xi,yi,hi21
  real :: termnorm,term,dy,dy2,dz,dz2,ypix,rescalefac
  real, dimension(npixx) :: dx2i

  datsmooth = 0.
  datnorm = 0.
  if (normalise) then
     print*,'taking fast cross section (normalised)...',zslice
  else
     print*,'taking fast cross section (non-normalised)...',zslice
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print*,'interpolate3D_xsec: error: pixel width <= 0'
     return
  elseif (npart.le.0) then
     print*,'interpolate3D_xsec: error: npart = 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D_xsec: WARNING: ignoring some or all particles with h < 0'
  endif
  const = cnormk3D
  !
  !--renormalise dat array by first element to speed things up
  !
  if (dat(1).gt.tiny(dat)) then
     rescalefac = dat(1)
  else
     rescalefac = 1.0
  endif
  !
  !--loop over particles
  !
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) cycle over_parts
     hi1 = 1./hi
     hi21 = hi1*hi1
     radkern = radkernel*hi    ! radius of the smoothing kernel
     !
     !--for each particle, work out distance from the cross section slice.
     !
     dz = zslice - z(i)
     dz2 = dz**2*hi21
     !
     !--if this is < 2h then add the particle's contribution to the pixels
     !  otherwise skip all this and start on the next particle
     !
     if (dz2 .lt. radkernel2) then

        xi = x(i)
        yi = y(i)
        termnorm = const*weight(i)
        term = termnorm*dat(i)/rescalefac
        !
        !--for each particle work out which pixels it contributes to
        !
        ipixmin = int((xi - radkern - xmin)/pixwidthx)
        jpixmin = int((yi - radkern - ymin)/pixwidthy)
        ipixmax = int((xi + radkern - xmin)/pixwidthx) + 1
        jpixmax = int((yi + radkern - ymin)/pixwidthy) + 1

        if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
        if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
        if (ipixmax.gt.npixx) ipixmax = npixx
        if (jpixmax.gt.npixy) jpixmax = npixy
        !
        !--precalculate an array of dx2 for this particle (optimisation)
        !
        do ipix=ipixmin,ipixmax
           dx2i(ipix) = ((xmin + (ipix-0.5)*pixwidthx - xi)**2)*hi21 + dz2
        enddo
        !
        !--loop over pixels, adding the contribution from this particle
        !
        do jpix = jpixmin,jpixmax
           ypix = ymin + (jpix-0.5)*pixwidthy
           dy = ypix - yi
           dy2 = dy*dy*hi21
           do ipix = ipixmin,ipixmax
              q2 = dx2i(ipix) + dy2
              !
              !--SPH kernel - standard cubic spline
              !
              if (q2.lt.radkernel2) then
                 wab = w_cubic(q2)
                 !
                 !--calculate data value at this pixel using the summation interpolant
                 !
                 datsmooth(ipix,jpix) = datsmooth(ipix,jpix) + term*wab
                 if (normalise) datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab

              endif

           enddo
        enddo

     endif                  ! if particle within 2h of slice
  enddo over_parts                    ! over particles
  !
  !--normalise dat array
  !
  if (normalise) then
     !--normalise everywhere (required if not using SPH weighting)
     where (datnorm > tiny(datnorm))
        datsmooth = datsmooth/datnorm
     end where
  endif
  datsmooth = datsmooth*rescalefac

  return

end subroutine interpolate3D_fastxsec

!--------------------------------------------------------------------------
!     program to interpolate from particle data to even grid of pixels
!
!     The data is smoothed using the SPH summation interpolant,
!     that is, we compute the smoothed array according to
!
!     datsmooth(pixel) = sum_b m_b dat_b/rho_b W(r-r_b, h_b)
!
!     where _b is the quantity at the neighbouring particle b and
!     W is the smoothing kernel, for which we use the usual cubic spline
!
!     ** In this version 3D data is interpolated to a single 2D cross section
!     ** This is much faster than interpolating to a 3D grid
!     ** and is efficient if only one or two cross sections are needed.
!
!     ** Note that the cross section is always taken in the z co-ordinate
!     ** so should submit the appropriate arrays as x, y and z.
!
!     Input: particle coordinates  : x,y,z (npart)
!            particle masses       : pmass (npart)
!            density on particles  : rho   (npart) - must be computed separately
!            smoothing lengths     : hh    (npart) - could be computed from density
!            vector data to smooth : vecx  (npart)
!                                    vecy  (npart)
!            cross section location: zslice
!
!     Output: smoothed vector field   : vecsmoothx (npixx,npixy)
!                                     : vecsmoothy (npixx,npixy)
!
!     Daniel Price, Institute of Astronomy, Cambridge, 23/9/03
!--------------------------------------------------------------------------
subroutine interpolate3D_xsec_vec(x,y,z,hh,weight,vecx,vecy,itype,npart,&
     xmin,ymin,zslice,vecsmoothx,vecsmoothy,npixx,npixy,pixwidthx,pixwidthy,normalise)

  implicit none
  integer, intent(in) :: npart,npixx,npixy
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,vecx,vecy
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidthx,pixwidthy,zslice
  real, intent(out), dimension(npixx,npixy) :: vecsmoothx, vecsmoothy
  logical, intent(in) :: normalise
  real, dimension(npixx,npixy) :: datnorm

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax
  real :: hi,hi1,radkern,q2,wab,const
  real :: termx,termy,termnorm,dx,dy,dz,dz2,xpix,ypix

  vecsmoothx = 0.
  vecsmoothy = 0.
  datnorm = 0.
  if (normalise) then
     print*,'taking fast cross section (normalised)...',zslice
  else
     print*,'taking fast cross section (non-normalised)...',zslice
  endif
  if (pixwidthx.le.0. .or. pixwidthy.le.0.) then
     print*,'interpolate3D_xsec_vec: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D_xsec_vec: WARNING: ignoring some or all particles with h < 0'
  endif
  const = cnormk3D ! normalisation constant (3D)
  !
  !--loop over particles
  !
  over_parts: do i=1,npart
     !
     !--skip particles with itype < 0
     !
     if (itype(i).lt.0) cycle over_parts
     !
     !--set kernel related quantities
     !
     hi = hh(i)
     if (hi.le.0.) cycle over_parts
     hi1 = 1./hi
     radkern = radkernel*hi    ! radius of the smoothing kernel
     !
     !--for each particle, work out distance from the cross section slice.
     !
     dz = zslice - z(i)
     dz2 = dz**2
     !
     !--if this is < 2h then add the particle's contribution to the pixels
     !  otherwise skip all this and start on the next particle
     !
     if (abs(dz) .lt. radkern) then
        termnorm = const*weight(i)
        termx = termnorm*vecx(i)
        termy = termnorm*vecy(i)
        !
        !--for each particle work out which pixels it contributes to
        !
        ipixmin = int((x(i) - radkern - xmin)/pixwidthx)
        jpixmin = int((y(i) - radkern - ymin)/pixwidthy)
        ipixmax = int((x(i) + radkern - xmin)/pixwidthx) + 1
        jpixmax = int((y(i) + radkern - ymin)/pixwidthy) + 1

        if (ipixmin.lt.1) ipixmin = 1 ! make sure they only contribute
        if (jpixmin.lt.1) jpixmin = 1 ! to pixels in the image
        if (ipixmax.gt.npixx) ipixmax = npixx
        if (jpixmax.gt.npixy) jpixmax = npixy
        !
        !--loop over pixels, adding the contribution from this particle
        !
        do jpix = jpixmin,jpixmax
           ypix = ymin + (jpix-0.5)*pixwidthy
           dy = ypix - y(i)
           do ipix = ipixmin,ipixmax
              xpix = xmin + (ipix-0.5)*pixwidthx
              dx = xpix - x(i)
              q2 = (dx*dx + dy*dy + dz2)*hi1*hi1
              !
              !--SPH kernel - standard cubic spline
              !
              if (q2.lt.radkernel2) then
                 wab = w_cubic(q2)
                 !
                 !--calculate data value at this pixel using the summation interpolant
                 !
                 vecsmoothx(ipix,jpix) = vecsmoothx(ipix,jpix) + termx*wab
                 vecsmoothy(ipix,jpix) = vecsmoothy(ipix,jpix) + termy*wab
                 if (normalise) datnorm(ipix,jpix) = datnorm(ipix,jpix) + termnorm*wab
              endif

           enddo
        enddo

     endif                  ! if particle within 2h of slice
  enddo over_parts                    ! over particles
  !
  !--normalise dat array(s)
  !
  if (normalise) then
     where (datnorm > tiny(datnorm))
        vecsmoothx = vecsmoothx/datnorm
        vecsmoothy = vecsmoothy/datnorm
     end where
  endif

  return

end subroutine interpolate3D_xsec_vec

!--------------------------------------------------------------------------
! $Id: interpolate3D_opacity.f90,v 1.16 2007/11/20 17:05:35 dprice Exp $
!
!     subroutine to do a ray trace through the particle data
!
!     we use the radiation transport equation along a ray, that is
!     the change in intensity from one side of a particle to the other is
!     given by:
!
!     I_nu = I_nu(0) exp(-tau_i) + S_nu (1 - exp(-tau_i))
!
!     where tau_i is the integrated optical depth through the particle,
!     and S_nu is the colour calculated from a colour table for the rendered data.
!     We calculate an intensity in red, green and blue for colour plots.
!
!     tau_i = kappa \int rho dz
!
!     this is calculated using the SPH kernel for rho, so for each pixel
!     the optical depth is incremented as the sum
!
!     tau_i = kappa \sum_j m_j \int W dz
!
!     where \int W dz is the SPH kernel integrated along one spatial dimension.
!     This is interpolated from a pre-calculated table (see module projections3D for this).
!
!     kappa is the monochromatic mass extinction coefficient
!     (particle cross section per unit mass) and is a constant for all particles
!     which must be given as input (although see below for calculations of a
!     meaningful values for kappa in terms of "surface depth in units of smoothing lengths")
!
!     Input: particle coordinates  : x,y,z (npart) - note that z is only required for perspective
!            particle masses       : pmass (npmass)
!            smoothing lengths     : hh    (npart)
!            weight                : m/(h^3 rho) (not used, but skips particles with w <= 0)
!            scalar data to smooth : dat   (npart)
!
!     Particle masses can be sent in as either a single scalar (npmass = 1)
!      or as an array of length npart (npmass=npart)
!
!     Settings: zobs, dz1 : settings for 3D projection
!               rkappa    : particle cross section per unit mass
!
!     Output: smoothed data            : datsmooth (npixx,npixy)
!             brightness array         : brightness (npixx,npixy)
!
!--------------------------------------------------------------------------
subroutine interp3D_proj_opacity(x,y,z,pmass,npmass,hh,weight,dat,zorig,itype,npart, &
     xmin,ymin,datsmooth,brightness,npixx,npixy,pixwidth,zobserver,dscreenfromobserver, &
     rkappa,zcut)

  implicit none
  real, parameter :: pi=3.1415926536
  integer, intent(in) :: npart,npixx,npixy,npmass
  real, intent(in), dimension(npart) :: x,y,z,hh,weight,dat,zorig
  real, intent(in), dimension(npmass) :: pmass
  integer, intent(in), dimension(npart) :: itype
  real, intent(in) :: xmin,ymin,pixwidth,zobserver,dscreenfromobserver, &
                      zcut,rkappa
  real, dimension(npixx,npixy), intent(out) :: datsmooth, brightness

  integer :: i,ipix,jpix,ipixmin,ipixmax,jpixmin,jpixmax,nused,nsink
  integer :: iprintinterval,iprintnext
  integer, dimension(npart) :: iorder
  integer(kind=selected_int_kind(12)) :: ipart
  real :: hi,hi1,hi21,radkern,q2,wab,pmassav
  real :: term,dy,dy2,ypix,zfrac,hav,zcutoff
  real :: fopacity,tau,rkappatemp,termi,xi,yi
  logical :: iprintprogress,adjustzperspective,rendersink
  real, dimension(npixx) :: xpix,dx2i
  real :: xminpix,yminpix
!#ifdef _OPENMP
!  integer :: OMP_GET_NUM_THREADS
!#else
  integer(kind=selected_int_kind(12)) :: iprogress
!#endif

  datsmooth = 0.
  term = 0.
  brightness = 0.
  print "(1x,a)",'ray tracing from particles to pixels...'
  if (pixwidth.le.0.) then
     print "(a)",'interpolate3D_opacity: error: pixel width <= 0'
     return
  endif
  if (any(hh(1:npart).le.tiny(hh))) then
     print*,'interpolate3D_opacity: warning: ignoring some or all particles with h < 0'
  endif
  !--check that npmass is sensible
  if (npmass.lt.1 .or. npmass.gt.npart) then
     print*,'interpolate3D_opacity: ERROR in input number of particle masses '
     return
  endif
  !--these values for npmass are not sensible but the routine will still work
  if (npmass.ne.1 .and. npmass.ne.npart) then
     print*,'WARNING: interpolate3D_opacity: number of particle masses input =',npmass
  endif

  if (abs(dscreenfromobserver).gt.tiny(dscreenfromobserver)) then
     adjustzperspective = .true.
     zcutoff = zobserver
  else
     adjustzperspective = .false.
     zcutoff = huge(zobserver)
  endif

!
!--kappa is the opacity in units of length^2/mass
!  sent as an input parameter as it should be kept constant throughout the simulation
!
!  However we compute a reasonable estimate below based on the current plot so that
!  we can give the "actual" optical depth for the current frame in terms of number of
!  smoothing lengths. This is purely for diagnostic purposes only.
!
!--calculate average h
  hav = sum(hh(1:npart))/real(npart)
!--average particle mass
  pmassav = sum(pmass(1:npmass))/real(npmass)
  rkappatemp = pi*hav*hav/(pmassav*coltable(0))
  print*,'average h = ',hav,' average mass = ',pmassav
  print "(1x,a,f6.2,a)",'typical surface optical depth is ~',rkappatemp/rkappa,' smoothing lengths'
  !
  !--print a progress report if it is going to take a long time
  !  (a "long time" is, however, somewhat system dependent)
  !
  iprintprogress = (npart .ge. 100000) .or. (npixx*npixy .gt.100000)
  !
  !--loop over particles
  !
  iprintinterval = 25
  if (npart.ge.1e6) iprintinterval = 10
  iprintnext = iprintinterval

!
!--first sort the particles in z so that we do the opacity in the correct order
!
  call indexx(npart,z,iorder)
!
!--store x value for each pixel (for optimisation)
!
  xminpix = xmin - 0.5*pixwidth
  yminpix = ymin - 0.5*pixwidth
  do ipix=1,npixx
     xpix(ipix) = xminpix + ipix*pixwidth
  enddo

  nused = 0
  nsink = 0

!!$OMP PARALLEL default(none) &
!!$OMP SHARED(hh,z,x,y,zorig,pmass,dat,itype,datsmooth,npmass,npart) &
!!$OMP SHARED(xmin,ymin,xminpix,yminpix,xpix,pixwidth) &
!!$OMP SHARED(npixx,npixy,dscreenfromobserver,zobserver,adjustzperspective) &
!!$OMP SHARED(zcut,zcutoff,iorder,rkappa,brightness) &
!!$OMP PRIVATE(hi,zfrac,xi,yi,radkern) &
!!$OMP PRIVATE(hi1,hi21,term,termi) &
!!$OMP PRIVATE(ipixmin,ipixmax,jpixmin,jpixmax) &
!!$OMP PRIVATE(dx2i,q2,ypix,dy,dy2,wab) &
!!$OMP PRIVATE(ipart,i,ipix,jpix,tau,fopacity) &
!!$OMP REDUCTION(+:nused)
!!$OMP MASTER
!#ifdef _OPENMP
!  print "(1x,a,i3,a)",'Using ',OMP_GET_NUM_THREADS(),' cpus'
!#endif
!!$OMP END MASTER

!!$OMP DO ORDERED SCHEDULE(dynamic)
  over_particles: do ipart=1,npart
     !
     !--report on progress
     !
!#ifndef _OPENMP
     if (iprintprogress) then
        iprogress = 100*(ipart/npart)
        if (iprogress.ge.iprintnext) then
           write(*,"('(',i3,'% -',i12,' particles done)')") iprogress,ipart
           iprintnext = iprintnext + iprintinterval
        endif
     endif
!#endif
     !
     !--render in order from back to front
     !
     i = iorder(ipart)
     !
     !--skip particles with itype < 0
     !
     if (itype(i) < 0) cycle over_particles

     !
     !--skip particles with weight < 0
     !  but not if weight == weight_sink (=-1)
     !
     rendersink = .false.
     if (abs(weight(i) - weight_sink) < tiny(0.)) then
        rendersink = .true.
     elseif (weight(i) <= 0.) then
        cycle over_particles
     endif

     !
     !--allow slicing [take only particles with z(unrotated) < zcut]
     !
     particle_within_zcut: if (zorig(i).lt.zcut .and. z(i).lt.zcutoff) then

     !  count particles within slice
     nused = nused + 1
     !
     !--adjust h according to 3D perspective
     !  need to be careful -- the kernel quantities
     !  change with z (e.g. radkern, r^2/h^2)
     !  but *not* the 1/h^2 in tau (because the change in 1/h^2 in tau
     !  would be cancelled by the corresponding change to h^2 in kappa)
     !
     hi = hh(i)
     if (hi.le.0.) then
        cycle over_particles
     elseif (adjustzperspective) then
        zfrac = abs(dscreenfromobserver/(z(i)-zobserver))
        hi = hi*zfrac
     endif

     !--these are the quantities used in the kernel r^2/h^2
     radkern = 2.*hi
     hi1 = 1./hi
     hi21 = hi1*hi1
     !--this is the term which multiplies tau
     if (npmass.eq.npart) then
        term = pmass(i)/(hh(i)*hh(i))
     else
        term = pmass(1)/(hh(i)*hh(i))
     endif
     !
     !--determine colour contribution of current point
     !  (work out position in colour table)
     !
!     dati = dat(i)
     xi = x(i)
     yi = y(i)
     termi = dat(i)
     !
     !--sink particles can have weight set to -1
     !  indicating that we should include them in the rendering
     !
     if (rendersink) then
        termi = pmass(i)/(4./3.*pi*hh(i)**3)  ! define "density" of a sink
        nsink = nsink + 1
     endif
     !
     !--for each particle work out which pixels it contributes to
     !
     ipixmin = int((xi - radkern - xmin)/pixwidth)
     jpixmin = int((yi - radkern - ymin)/pixwidth)
     ipixmax = int((xi + radkern - xmin)/pixwidth) + 1
     jpixmax = int((yi + radkern - ymin)/pixwidth) + 1

     if (ipixmin.lt.1) ipixmin = 1  ! make sure they only contribute
     if (jpixmin.lt.1) jpixmin = 1  ! to pixels in the image
     if (ipixmax.gt.npixx) ipixmax = npixx ! (note that this optimises
     if (jpixmax.gt.npixy) jpixmax = npixy !  much better than using min/max)
     !
     !--precalculate an array of dx2 for this particle (optimisation)
     !
     do ipix=ipixmin,ipixmax
        dx2i(ipix) = ((xpix(ipix) - xi)**2)*hi21
     enddo

     !
     !--loop over pixels, adding the contribution from this particle
     !
     do jpix = jpixmin,jpixmax
        ypix = yminpix + jpix*pixwidth
        dy = ypix - yi
        dy2 = dy*dy*hi21
        do ipix = ipixmin,ipixmax
           q2 = dx2i(ipix) + dy2
           !
           !--SPH kernel - integral through cubic spline
           !  interpolate from a pre-calculated table
           !
           if (q2.lt.radkernel2) then
              wab = wfromtable(q2)
              !
              !--get incremental tau for this pixel from the integrated SPH kernel
              !
              tau = rkappa*wab*term
              fopacity = 1. - exp(-tau)
              !
              !--render, obscuring previously drawn pixels by relevant amount
              !  also calculate total brightness (`transparency') of each pixel
              !
              datsmooth(ipix,jpix) = (1.-fopacity)*datsmooth(ipix,jpix) + fopacity*termi
              brightness(ipix,jpix) = brightness(ipix,jpix) + fopacity
           endif

        enddo
     enddo

     endif particle_within_zcut

  enddo over_particles
!!$OMP END DO
!!$OMP END PARALLEL

  if (nsink > 99) then
     print*,'rendered ',nsink,' sink particles'
  elseif (nsink > 0) then
     print "(1x,a,i2,a)",'rendered ',nsink,' sink particles'
  endif
  if (zcut.lt.huge(zcut)) print*,'slice contains ',nused,' of ',npart,' particles'

  return

end subroutine interp3D_proj_opacity

!---------------------------------------
!
!  Functional forms of various kernels
!
!--------------------------------------
pure real function w_cubic(q2)
 implicit none
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
 implicit none
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
 implicit none
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
 implicit none
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
 implicit none
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
 implicit none
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
 implicit none
 real, intent(in) :: q2
 real :: q

 if (q2.lt.4.) then
    q = sqrt(q2)
    w_wendlandc6 = (1. - 0.5*q)**8*(4.*q2*q + 25./4.*q2 + 4.*q + 1.)
 else
    w_wendlandc6 = 0.
 endif

end function w_wendlandc6

subroutine indexx(n, arr, indx)
!************************************************************
!                                                           *
!  This is INDEXX using the quicksort algorithm.            *
!                                                           *
!************************************************************
 implicit none
 integer, parameter :: m=7, nstack=500
 integer, intent(in) :: n
 real, dimension(n), intent(in) :: arr
 integer, dimension(n), intent(out) :: indx

 integer :: i,j,k,l,ir,jstack,indxt,itemp
 integer, dimension(nstack) :: istack
 real :: a

 do j = 1, n
    indx(j) = j
 enddo
 jstack = 0
 l = 1
 ir = n

1 if (ir - l.lt.m) then
     do j = l + 1, ir
        indxt = indx(j)
        a = arr(indxt)
        do i = j - 1, 1, -1
           if (arr(indx(i)).le.a) goto 2
           indx(i + 1) = indx(i)
        end do
        i = 0
2       indx(i + 1) = indxt
     end do
     if (jstack.eq.0) return
     ir = istack(jstack)
     l = istack(jstack - 1)
     jstack = jstack - 2
  else
     k = (l + ir)/2
     itemp = indx(k)
     indx(k) = indx(l + 1)
     indx(l + 1) = itemp
     if (arr(indx(l + 1)).gt.arr(indx(ir))) then
        itemp = indx(l + 1)
        indx(l + 1) = indx(ir)
        indx(ir) = itemp
     endif
     if (arr(indx(l)).gt.arr(indx(ir))) then
        itemp = indx(l)
        indx(l) = indx(ir)
        indx(ir) = itemp
     endif
     if (arr(indx(l + 1)).gt.arr(indx(l))) then
        itemp = indx(l + 1)
        indx(l + 1) = indx(l)
        indx(l) = itemp
     endif
     i = l + 1
     j = ir
     indxt = indx(l)
     a = arr(indxt)
3    continue
     i = i + 1
     if (arr(indx(i)).lt.a) goto 3
4    continue
     j = j - 1
     if (arr(indx(j)).gt.a) goto 4
     if (j.lt.i) goto 5
     itemp = indx(i)
     indx(i) = indx(j)
     indx(j) = itemp
     goto 3

5    indx(l) = indx(j)
     indx(j) = indxt
     jstack = jstack + 2
     if (jstack.gt.nstack) then
        print*,'fatal error!!! stacksize exceeded in sort'
        print*,'(need to set parameter nstack higher in subroutine indexx '
        print*,' this is in the file sort.f90)'
        stop
     endif
     if (ir - i + 1.ge.j - l) then
        istack(jstack) = ir
        istack(jstack - 1) = i
        ir = j - 1
     else
        istack(jstack) = j - 1
        istack(jstack - 1) = l
        l = i
     endif
  endif

goto 1
end subroutine indexx
end module splash
