!-------------------------------------------------------------------------------

!                      Code_Saturne version 2.1.4
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine uscfth &
!================

 ( nvar   , nscal  ,                                              &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   sorti1 , sorti2 , gamagr , xmasmr )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Define thermodynamic laws (especially for the compressible flow scheme).

!    This user subroutine is mandatory for the compressible flow scheme.


! Introduction
! ============

! This user subroutine allows to define all physical properties and
! variables, through the implementation of thermodynamic laws.


! Avalable thermodynamic laws
! ===========================

!  1. Perfect gas (the molar mass 'xmasml' must be provided)
!  2. Perfect gas with non constant Gamma (example to be adapted)
!  3. Van Der Waals (not yet implemented)


! Implemented calculations
! ========================

! This user subroutine implements the computation of several quantities.
! Each calculation has to be explicitly implemented in the appropriate
! section below (already done for perfect gas).


! Selection of the quantity to return
! ===================================

! When calling the user subroutine, the integer 'iccfth' specifies which
! calculation has to be performed (and which quantity has to be returned).
! The values for 'iccfth' for each case are provided below.
! For some configurations, two systems of references are used for 'iccfth'
! (this is useful to make tests easier to implement in the calling
! subroutines): both systems are explained hereafter for information.

! First system:

!   the variables are referred to using an index i:
!     Variable  P  rho  T   e   h   s  'internal energy - CvT'
!        Index  1   2   3   4   5   6              7

!   iccfth is as follows, depending on which quantity needs to be computed:
!     - compute all variables at cell centers from variable i
!                                              and variable j (i<j):
!               => iccfth = 10*i+j
!     - compute all variables at boundary faces from variable i
!                                                and variable j (i<j):
!               => iccfth = 10*i+j+900

! Second system:

!   the variables are referred to using a different index i:
!     Variable  P  rho  T  e  s
!        Index  2   3   5  7  13

!   iccfth is as follows, depending on which quantity needs to be computed:
!     - compute all variables at cell centers from variable i
!                                              and variable j (i<j):
!               => iccfth = i*j*10000
!     - compute all variables at boundary faces from variable i
!                                                and variable j (i<j):
!               => iccfth = i*j*10000+900

! Other quantities:

!   the variables are referred to using the index of the first system.
!   iccfth is defined as follows:
!     - compute variable i at cell centers (for s and 'internal energy-CvT')
!               => iccfth = i
!                                   \partial(variable i)|
!     - compute partial derivative  --------------------|
!                                   \partial(variable j)|variable k
!               => iccfth = 100*i+10*j+k
!     - compute boundary conditions, resp. symmetry, wall, inlet, outlet:
!               => iccfth = 91, 92, 93, 94


! Values of iccfth
! ================

! To summarize, the values for iccfth are as follows:

!   Values at the cell centers:

!   -> set calculation options (cst/variable cp)   : iccfth = -1
!   -> set default initialization                  : iccfth =  0
!   -> calculate gamma                             : iccfth =  1
!   -> verification of the density                 : iccfth = -2
!   -> verification of the energy                  : iccfth = -4
!   -> calculation of temperature and energy
!                     from pressure and density    : iccfth =  12 or  60000
!   -> calculation of density and energy
!                     from pressure and temperature: iccfth =  13 or 100000
!   -> calculation of density and temperature
!                     from pressure and energy     : iccfth =  14 or 140000
!   -> calculation of pressure and energy
!                     from density and temperature : iccfth =  23 or 150000
!   -> calculation of pressure and temperature
!                     from density and energy      : iccfth =  24 or 210000
!
!                      2    dP |
!   -> calculation of c  = ----|                   : iccfth = 126
!                          drho|s
!
!                            dP|
!   -> calculation of beta = --|                   : iccfth = 162
!                            ds|rho
!
!                          de|
!   -> calculation of Cv = --|                     : iccfth = 432
!                          dT|rho
!
!   -> calculation of entropie                     : iccfth =   6
!
!
!   Values at the boundary faces
!
!   -> calculation of the boundary conditions:
!     - symmetry                                   : iccfth =  90
!     - wall                                       : iccfth =  91
!     - inlet                                      : iccfth =  92
!     - outlet                                     : iccfth =  93
!     - different outlet,not implemented yet       : iccfth =  94
!
!   -> calculation of the variables at the faces for boundary conditions:
!     - temperature and energy
!         from pressure and density                : iccfth = 912 ou  60900
!     - density and energy
!         from pressure and temperature            : iccfth = 913 ou 100900
!     - density and temperature
!         from pressure and energy                 : iccfth = 914 ou 140900
!     - pressure and energy
!         from density and temperature             : iccfth = 923 ou 150900
!     - pressure and temperature
!         from density and energy                  : iccfth = 924 ou 210900


!   Values at the cell centers and at the boundary faces

!   -> calculation of 'internal energy - Cv.T'     : iccfth =   7

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! sorti1,2(*)      ! ra ! --> ! output variable (unused if iccfth.lt.0)        !
! gamagr(*)        ! ra ! --> ! equivalent "gamma" constant of the gas         !
!                  !    !     !   (unused if iccfth.lt.0)                      !
!                  !    !     !   (first value only used for perfect gas)      !
! xmasmr(*)        ! ra ! --> ! molar mass of the components of the gas        !
!                  !    !     !   (unused if iccfth.lt.0)                      !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use parall
use pointe
use entsor
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          iccfth   , imodif

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*),propfa(nfac,*),propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision sorti1(*), sorti2(*), gamagr(*), xmasmr(*)

! Local variables

integer          ifac0
integer          ierr
integer          iel    , ifac   , ivar
integer          irh    , itk    , ien
integer          iclp   , iclr   , iclt   , icle
integer          iclu   , iclv   , iclw
integer          iutile
double precision gamagp , xmasml , enint
double precision xmach  , xmachi , xmache , dxmach

integer          npmax
parameter (npmax = 1000)
double precision cstgr(npmax)

!===============================================================================


!===============================================================================
! 0. Initialization.
!    No user input required.
!===============================================================================

! Error indicator (stop if non zero)
ierr   = 0

! Rank of the variables in their associated arrays
if(iccfth.ge.0.or.iccfth.le.-2) then
  irh = isca(irho  )
  itk = isca(itempk)
  ien = isca(ienerg)
  iclp = iclrtp(ipr,icoef)
  iclr = iclrtp(irh,icoef)
  iclt = iclrtp(itk,icoef)
  icle = iclrtp(ien,icoef)
  iclu = iclrtp(iu,icoef)
  iclv = iclrtp(iv,icoef)
  iclw = iclrtp(iw,icoef)
endif

! For calculation of values at the cell centers,
!   ifac0 > indicates that the array rtp must be modified
! For calculation of values at the cell faces,
!   ifac0 is the number of the current face
ifac0 = imodif

!===============================================================================
! 1. Thermodynamic law choice (for each phase)
!    User input required.
!===============================================================================

! --> ieos = 1: Perfect gas with constant Gamma
! --> ieos = 2: Perfect gas with variable Gamma
! --> ieos = 3: Van Der Waals

ieos = 1


! Warning: once the thermodynamic law has been chosen,
! =======  the remainder of the user subroutine must be modified


!===============================================================================
! 2. Perfect gas
!===============================================================================

if(ieos.eq.1) then

!===============================================================================
! 2.1. Parameters to be completed by the user
!===============================================================================

! For each phase

! --- Molar mass of the gas (kg/mol)

  if(iccfth.ge.0) then
    xmasml = 28.8d-3
  endif

!===============================================================================
! 2.2. Default laws
!      No user input required.
!===============================================================================

! --- Calculation of the constant gamagp

  if(iccfth.gt.0) then

    ! Gamagp is supposed to be superior or equal to 1.
    ! It is computed at each call, even if this may seem costly,
    !   to be coherent with the "constant gamma" case for which this
    !   constant is not saved. A ''save'' instruction and a test would
    !   be sufficient to avoid computing gamagp at each call if necessary.

    gamagp = 1.d0 + rr/(xmasml*cp0-rr)

    if(gamagp.lt.1.d0) then
      write(nfecra,1010) gamagp
      call csexit (1)
    endif

    ! Gamma is returned if required

    if(iccfth.eq.1) then
      gamagr(1) = gamagp
    endif

  endif


! --- Calculation options: constant Cp and Cv (perfect gas)

  if(iccfth.eq.-1) then

    ! The value for the isobaric specific heat Cp0 must be provided in
    !   the user subroutine ''usini1''. The value for the isochoric
    !   specific heat Cv0 is calculated in a subsequent section (from Cp0)

    icp = 0
    icv = 0


! --- Default initializations (before uscfxi)

!     T0 is positive (this assumption has been checked in
!       the user programme 'verini')

  elseif(iccfth.eq.0) then

    cv0 = cp0 - rr/xmasml

    if ( isuite .eq. 0 ) then
      do iel = 1, ncel
        rtp(iel,irh) = p0*xmasml/(rr*t0)
        rtp(iel,ien) = cv0*t0
      enddo
    endif


! --- Verification of the density

  elseif(iccfth.eq.-2) then

    ! If the density is lower or equal to zero: clipping, write and stop.
    !   Indeed, if this is the case, the thermodynamic computations will
    !   most probably fail.
    ! This call is done at the end of the density calculation (after
    !   a classical clipping and before parallel communications).

    ierr = 0
    do iel = 1, ncel
      if(rtp(iel,irh).le.0.d0) then
        rtp(iel,irh) = epzero
        ierr = ierr + 1
      endif
    enddo
    if(irangp.ge.0) then
      call parcpt (ierr)
    endif
    if(ierr.gt.0) then
      ntmabs = ntcabs
      write(nfecra,8000)ierr, epzero
    endif


! --- Verification of the energy

  elseif(iccfth.eq.-4) then

    ! If the total energy <= zero: clipping, write and stop
    !   Indeed, if this is the case, the thermodynamic computations will
    !   most probably fail.

    ierr = 0
    do iel = 1, ncel
      enint = rtp(iel,ien)                                     &
               - 0.5d0*( rtp(iel,iu)**2                        &
                       + rtp(iel,iv)**2                        &
                       + rtp(iel,iw)**2 )
      if(enint.le.0.d0) then
        rtp(iel,ien) = epzero                                  &
               + 0.5d0*( rtp(iel,iu)**2                        &
                       + rtp(iel,iv)**2                        &
                       + rtp(iel,iw)**2 )
        ierr = ierr + 1
      endif
    enddo
    if(irangp.ge.0) then
      call parcpt (ierr)
    endif
    if(ierr.gt.0) then
      ntmabs = ntcabs
      write(nfecra,8100)ierr, epzero
    endif


! --- Calculation of temperature and energy from pressure and density

  elseif(iccfth.eq.12.or.iccfth.eq.60000) then

    ! Verification of the values of the density
    ierr = 0
    do iel = 1, ncel
      if(rtp(iel,irh).le.0.d0) then
        write(nfecra,3010)rtp(iel,irh),iel
      endif
    enddo
    ! Stop if a negative value is detected (since the density has been
    ! provided by the user, one potential cause is a wrong user
    ! initialization)
    if(ierr.eq.1) then
      call csexit (1)
    endif

    do iel = 1, ncel
      ! Temperature
      sorti1(iel) = xmasml*rtp(iel,ipr)/(rr*rtp(iel,irh))
      ! Total energy
      sorti2(iel) = cv0*sorti1(iel)                        &
           + 0.5d0*( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 )
    enddo

    ! Transfer to the array rtp
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,itk) = sorti1(iel)
        rtp(iel,ien) = sorti2(iel)
      enddo
    endif


! --- Calculation of density and energy from pressure and temperature:

  elseif(iccfth.eq.13.or.iccfth.eq.100000) then

    ! Verification of the values of the temperature
    ierr = 0
    do iel = 1, ncel
      if(rtp(iel,itk).le.0.d0) then
        write(nfecra,2010)rtp(iel,itk),iel
      endif
    enddo
    ! Stop if a negative value is detected (since the temperature has been
    ! provided by the user, one potential cause is a wrong user
    ! initialization: a value not provided in Kelvin for example)
    if(ierr.eq.1) then
      call csexit (1)
    endif

    do iel = 1, ncel
      ! Density
      sorti1(iel) = xmasml*rtp(iel,ipr)/(rr*rtp(iel,itk))
      ! Total energy
      sorti2(iel) = cv0*rtp(iel,itk)                    &
           + 0.5d0*( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 )
    enddo

    ! Transfer to the array rtp
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,irh) = sorti1(iel)
        rtp(iel,ien) = sorti2(iel)
      enddo
    endif


! --- Calculation of density and temperature from pressure and energy

  elseif(iccfth.eq.14.or.iccfth.eq.140000) then

    do iel = 1, ncel
      ! Internal energy (to avoid the need to divide by the temperature
      ! to compute density)
      enint = rtp(iel,ien)                                     &
               - 0.5d0*( rtp(iel,iu)**2                        &
                       + rtp(iel,iv)**2                        &
                       + rtp(iel,iw)**2 )
      ! Density
      sorti1(iel) = rtp(iel,ipr) / ( (gamagp-1.d0) * enint )
      ! Temperature
      sorti2(iel) = xmasml * (gamagp-1.d0) * enint / rr
    enddo

    ! Transfer to the array rtp
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,irh) = sorti1(iel)
        rtp(iel,itk) = sorti2(iel)
      enddo
    endif


! --- Calculation of pressure and energy from density and temperature

  elseif(iccfth.eq.23.or.iccfth.eq.150000) then

    do iel = 1, ncel
      ! Pressure
      sorti1(iel) = rtp(iel,irh)*rtp(iel,itk)*rr/xmasml
      ! Total energy
      sorti2(iel) = cv0*rtp(iel,itk)                    &
           + 0.5d0*( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 )
    enddo

    ! Transfer to the array rtp
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,ipr) = sorti1(iel)
        rtp(iel,ien) = sorti2(iel)
      enddo
    endif


! --- Calculation of pressure and temperature from density and energy

  elseif(iccfth.eq.24.or.iccfth.eq.210000) then

    do iel = 1, ncel
      ! Internal energy (to avoid the need to divide by the temperature
      ! to compute density)
      enint = rtp(iel,ien)                                     &
               - 0.5d0*( rtp(iel,iu)**2                        &
                       + rtp(iel,iv)**2                        &
                       + rtp(iel,iw)**2 )
      ! Pressure
      sorti1(iel) = (gamagp-1.d0) * rtp(iel,irh) * enint
      ! Temperature
      sorti2(iel) = xmasml * (gamagp-1.d0) * enint / rr
    enddo

    ! Transfer to the array rtp
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,ipr) = sorti1(iel)
        rtp(iel,itk) = sorti2(iel)
      enddo
    endif


!                     2                            2         P
! --- Calculation of c from pressure and density: c = gamma*---
!                                                           rho

  elseif(iccfth.eq.126) then

    ! Verification of the values of the density
    !   This test can be discarded to reduce the CPU time (if
    !     density is <= 0, the calculation will simply fail)
    !   It is discarded here with iutile = 0
    iutile = 0
    if(iutile.eq.1) then
      ierr = 0
      do iel = 1, ncel
        if(rtp(iel,irh).le.0.d0) then
          write(nfecra,4010)rtp(iel,irh),iel
        endif
      enddo
      if(ierr.eq.1) then
        call csexit (1)
      endif
    endif

    do iel = 1, ncel
      sorti1(iel) = gamagp * rtp(iel,ipr) / rtp(iel,irh)
    enddo


!                                                              gamma
! --- Calculation of beta from pressure and density: beta = rho

  elseif(iccfth.eq.162) then

    ! Verification of the values of the density
    !   This test can be discarded to reduce the CPU time (if
    !     density is <= 0, the calculation will simply fail)
    !   It is discarded here with iutile = 0
    iutile = 0
    if(iutile.eq.1) then
      ierr = 0
      do iel = 1, ncel
        if(rtp(iel,irh).lt.0.d0) then
          write(nfecra,4020)rtp(iel,irh),iel
        endif
      enddo
      if(ierr.eq.1) then
        call csexit (1)
      endif
    endif

    do iel = 1, ncel
      sorti1(iel) = rtp(iel,irh)**gamagp
    enddo


! --- Calculation of the isochoric specific heat

    ! It is a constant: nothing to do


!                                                                  P
! --- Calculation of the entropy from pressure and density: s = --------
!                                                                  gamma
!                                                               rho

  elseif(iccfth.eq.6) then

    ! Verification of the values of the density
    !   This test can be discarded to reduce the CPU time (if
    !     density is <= 0, the calculation will simply fail)
    ierr = 0
    do iel = 1, ncel
      if(rtp(iel,irh).le.0.d0) then
        write(nfecra,4030)rtp(iel,irh),iel
      endif
    enddo
    if(ierr.eq.1) then
      call csexit (1)
    endif

    do iel = 1, ncel
      sorti1(iel) = rtp(iel,ipr) / (rtp(iel,irh)**gamagp)
    enddo


! --- Calculation of 'internal energy - Cv.T'

  elseif(iccfth.eq.7) then

    ! It is zero for a perfect gas

    !   At the cell centers
    do iel = 1, ncel
      sorti1(iel) = 0.d0
    enddo

    !   On the boundary faces
    do ifac = 1, nfabor
      sorti2(ifac) = 0.d0
    enddo


! --- Calculation of the boundary conditions on the face ifac = ifac0

!  -- Wall

  elseif(iccfth.eq.91) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Calculation of the Mach number at the boundary face, using the
    !   cell center velocity projected on the vector normal to the boundary
    xmach =                                                       &
         ( rtp(iel,iu)*surfbo(1,ifac)                          &
         + rtp(iel,iv)*surfbo(2,ifac)                          &
         + rtp(iel,iw)*surfbo(3,ifac) ) / surfbn(ifac)         &
         / sqrt( gamagp * rtp(iel,ipr) / rtp(iel,irh) )

    ! Pressure

    !   A Neumann boundary condition is used. This does not allow to use
    !     the Rusanov scheme, but some stabilization effect is expected.
    !     A test based on the value of coefb at the previous time step
    !     is implemented to avoid oscillating between a rarefaction
    !     situation and a shock configuration from one time step to the
    !     next.

    !   Rarefaction
    if(xmach.lt.0.d0.and.coefb(ifac,iclp).le.1.d0) then

      if(xmach.gt.2.d0/(1.d0-gamagp)) then
        coefb(ifac,iclp) = (1.d0 + (gamagp-1.d0)/2.d0 * xmach)    &
             ** (2.d0*gamagp/(gamagp-1.d0))
      else
        ! In case the rarefaction is too strong, a zero Dirichlet value
        !   is used for pressure (the value of coefb is used here as an
        !   indicator and will be modified later in cfxtcl)
        coefb(ifac,iclp) = rinfin
      endif

      !  Shock
    elseif(xmach.gt.0.d0.and.coefb(ifac,iclp).ge.1.d0) then

      coefb(ifac,iclp) = 1.d0 + gamagp*xmach                      &
            *( (gamagp+1.d0)/4.d0*xmach                           &
                + sqrt(1.d0 + (gamagp+1.d0)**2/16.d0*xmach**2) )

      !  Oscillation between rarefaction and shock or zero Mach number
    else
      coefb(ifac,iclp) = 1.d0
    endif


!  -- Symmetry

  elseif(iccfth.eq.90) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! A zero flux condition (homogeneous Neumann condition) is
    !   prescribed by default.
    ! No user input required


!  -- Subsonic inlet with prescribed density and velocity

    ! The subsonic nature of the inlet is postulated.

    ! Further testing may be required here. Contrary to the initial
    !   development, an explicit Dirichlet condition is prescribed for
    !   pressure instead of a Neumann condition (however, the same
    !   physical value for pressure is used).
    ! The advantage of this approach is to allow the use of the Rusanov
    !   scheme to stabilize the user defined inlet conditions.
    ! Moreover, with this approach, coefb does not have to be filled in
    !   here (it is not a major point, since coefb has to be filled in
    !   for the wall boundary condition anyway)
    ! Shall an oscillatory behavior (in time) be observed, it might be
    !   worth trying to add a test to avoid switching between
    !   rarefaction and shock from one time step to the other (just as
    !   for the wall boundary condition).
    ! The relevance of this approach remains to be demonstrated.

  elseif(iccfth.eq.92) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Calculation of the Mach number at the boundary face, using the
    !   cell center velocity projected on the vector normal to the boundary
    xmachi =                                                      &
         ( rtp(iel,iu)*surfbo(1,ifac)                          &
         + rtp(iel,iv)*surfbo(2,ifac)                          &
         + rtp(iel,iw)*surfbo(3,ifac) ) / surfbn(ifac)         &
         / sqrt( gamagp * rtp(iel,ipr) / rtp(iel,irh) )
    xmache =                                                      &
         ( coefa(ifac,iclu)*surfbo(1,ifac)                        &
         + coefa(ifac,iclv)*surfbo(2,ifac)                        &
         + coefa(ifac,iclw)*surfbo(3,ifac) ) /surfbn(ifac)        &
         / sqrt( gamagp * rtp(iel,ipr) / rtp(iel,irh) )
    dxmach = xmachi - xmache

    ! Pressure: rarefaction wave (Rusanov)
    if(dxmach.le.0.d0) then

      if(dxmach.gt.2.d0/(1.d0-gamagp)) then
        coefa(ifac,iclp) = rtp(iel,ipr)*                       &
             ( (1.d0 + (gamagp-1.d0)*0.50d0*dxmach)               &
               ** (2.d0*gamagp/(gamagp-1.d0))    )
      elseif(dxmach.le.2.d0/(1.d0-gamagp) ) then
        coefa(ifac,iclp) = 0.d0
      endif

      ! Pressure: shock (Rusanov)
    else
      coefa(ifac,iclp) = rtp(iel,ipr)*                         &
           (  1.d0 + gamagp*dxmach                                &
           *( (gamagp+1.d0)*0.25d0*dxmach                         &
           + sqrt(1.d0 + (gamagp+1.d0)**2/16.d0*dxmach**2) )  )
    endif

    ! This choice overrides the previous Rusanov choice
    coefa(ifac,iclp) = rtp(iel,ipr)

    ! Total energy
    coefa(ifac,icle) =                                            &
         coefa(ifac,iclp)/((gamagp-1.d0)*coefa(ifac,iclr))        &
         + 0.5d0*(coefa(ifac,iclu)**2                             &
                + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)


!  -- Subsonic inlet with prescribed mass and enthalpy flow rates
    ! The quantities prescribed are rho*u and rho*u*h

    ! The subsonic nature of the inlet is postulated.

    ! This section remains to be implemented: stop for the moment

    ! One may proceed as follows:
    !   Pressure computed with a Newton method
    !   Velocity and density computed from pressure
    !   Total energy computed from enthalpy
    !   (written on paper, to be implemented: contact the user support)

  elseif(iccfth.eq.94) then

    ifac = ifac0
    iel  = ifabor(ifac)

    write(nfecra,7000)

    call csexit (1)
    !==========


!  -- Subsonic outlet

    ! The subsonic nature of the inlet is postulated.

  elseif(iccfth.eq.93) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Rarefaction case
    if(coefa(ifac,iclp).le.rtp(iel,ipr)) then

      ! Density
      coefa(ifac,iclr) = rtp(iel,irh)                          &
           * (coefa(ifac,iclp)/rtp(iel,ipr))**(1.d0/gamagp)

      ! Velocity
      coefa(ifac,iclu) = rtp(iel,iu)                           &
           + 2.d0/(gamagp-1.d0)                                   &
           * sqrt(gamagp*rtp(iel,ipr)/rtp(iel,irh))         &
           * (1.d0-(coefa(ifac,iclp)/rtp(iel,ipr)              &
                        )**((gamagp-1.d0)/(2.d0*gamagp)))         &
           * surfbo(1,ifac)/surfbn(ifac)

      coefa(ifac,iclv) = rtp(iel,iv)                           &
           + 2.d0/(gamagp-1.d0)                                   &
           * sqrt( gamagp*rtp(iel,ipr)/rtp(iel,irh))        &
           * (1.d0-(coefa(ifac,iclp)/rtp(iel,ipr)              &
                        )**((gamagp-1.d0)/(2.d0*gamagp)))         &
           * surfbo(2,ifac)/surfbn(ifac)

      coefa(ifac,iclw) = rtp(iel,iw)                           &
           + 2.d0/(gamagp-1.d0)                                   &
           * sqrt( gamagp*rtp(iel,ipr)/rtp(iel,irh))        &
           * (1.d0-(coefa(ifac,iclp)/rtp(iel,ipr)              &
                        )**((gamagp-1.d0)/(2.d0/gamagp)))         &
           * surfbo(3,ifac)/surfbn(ifac)

      ! Total energy
      coefa(ifac,icle) =                                          &
           coefa(ifac,iclp)/((gamagp-1.d0)*coefa(ifac,iclr))      &
           + 0.5d0*(coefa(ifac,iclu)**2                           &
                  + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)

    ! Shock
    else

      ! Density
      coefa(ifac,iclr) = rtp(iel,irh)                          &
           * ( (gamagp+1.d0)*coefa(ifac,iclp)                     &
             + (gamagp-1.d0)*rtp(iel,ipr) )                    &
           / ( (gamagp-1.d0)*coefa(ifac,iclp)                     &
             + (gamagp+1.d0)*rtp(iel,ipr) )

      ! Velocity
      coefa(ifac,iclu) = rtp(iel,iu)                           &
           - (coefa(ifac,iclp)-rtp(iel,ipr))                   &
           * sqrt(2.d0/                                           &
                  (rtp(iel,irh)                                &
                   *((gamagp+1.d0)*coefa(ifac,iclp)               &
                    +(gamagp-1.d0)*rtp(iel,ipr) )))            &
           * surfbo(1,ifac)/surfbn(ifac)

      coefa(ifac,iclv) = rtp(iel,iv)                           &
           - (coefa(ifac,iclp)-rtp(iel,ipr))                   &
           * sqrt(2.d0/                                           &
                  (rtp(iel,irh)                                &
                   *((gamagp+1.d0)*coefa(ifac,iclp)               &
                    +(gamagp-1.d0)*rtp(iel,ipr) )))            &
           * surfbo(2,ifac)/surfbn(ifac)

      coefa(ifac,iclw) = rtp(iel,iw)                           &
           - (coefa(ifac,iclp)-rtp(iel,ipr))                   &
           * sqrt(2.d0/                                           &
                  (rtp(iel,irh)                                &
                   *((gamagp+1.d0)*coefa(ifac,iclp)               &
                    +(gamagp-1.d0)*rtp(iel,ipr) )))            &
           * surfbo(3,ifac)/surfbn(ifac)

      ! Total energy
      coefa(ifac,icle) =                                          &
           coefa(ifac,iclp)/((gamagp-1.d0)*coefa(ifac,iclr))      &
           + 0.5d0*(coefa(ifac,iclu)**2                           &
                  + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)

    endif


! --- Calculation of temperature and energy from pressure and density

    ! It is postulated that the pressure and density values are
    !   strictly positive

  elseif(iccfth.eq.912.or.iccfth.eq.60900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Temperature
    coefa(ifac,iclt) =                                            &
         xmasml*coefa(ifac,iclp)/(rr*coefa(ifac,iclr))

    ! Energie totale
    coefa(ifac,icle) =                                            &
         cv0*coefa(ifac,iclt)                              &
         + 0.5d0*( coefa(ifac,iclu)**2                            &
                 + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2 )


! --- Calculation of density and energy from pressure and temperature

  elseif(iccfth.eq.913.or.iccfth.eq.100900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Density
    coefa(ifac,iclr) =                                            &
         xmasml*coefa(ifac,iclp)/(rr*coefa(ifac,iclt))

    ! Total energy
    coefa(ifac,icle) =                                            &
         cv0*coefa(ifac,iclt)                              &
         + 0.5d0*( coefa(ifac,iclu)**2                            &
                 + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2 )


! --- Calculation of density and temperature from pressure and total energy

  elseif(iccfth.eq.914.or.iccfth.eq.140900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Density
    coefa(ifac,iclr) = coefa(ifac,iclp)/( (gamagp-1.d0)*          &
         (coefa(ifac,icle)                                        &
         - 0.5d0*( coefa(ifac,iclu)**2                            &
                 + coefa(ifac,iclv)**2                            &
                 + coefa(ifac,iclw)**2 ) ) )

    ! Temperature
    coefa(ifac,iclt)=                                             &
         xmasml*coefa(ifac,iclp)/(rr*coefa(ifac,iclr))


! --- Calculation of pressure and energy from density and temperature

  elseif(iccfth.eq.923.or.iccfth.eq.150900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Pressure
    coefa(ifac,iclp) = coefa(ifac,iclr)*rr/xmasml                 &
                                       *coefa(ifac,iclt)

    ! Total energy
    coefa(ifac,icle) = cv0 * coefa(ifac,iclt)              &
         + 0.5d0*( coefa(ifac,iclu)**2                            &
                 + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2 )


! --- Calculation of pressure and temperature from density and energy

  elseif(iccfth.eq.924.or.iccfth.eq.210900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Pressure
    coefa(ifac,iclp) = (gamagp-1.d0)*coefa(ifac,iclr)             &
          *( coefa(ifac,icle)                                     &
            - 0.5d0*( coefa(ifac,iclu)**2                         &
                    + coefa(ifac,iclv)**2                         &
                    + coefa(ifac,iclw)**2 ) )


    ! Temperature
    coefa(ifac,iclt)=                                             &
         xmasml*coefa(ifac,iclp)/(rr*coefa(ifac,iclr))


! --- End of the treatment of the perfect gas
  endif


!===============================================================================
! 3. Perfect gas with variable gamma
!===============================================================================

! This section requires further checking and testing

elseif(ieos.eq.2) then

!===============================================================================

!===============================================================================
! 3.1. Parameters to be completed by the user
!===============================================================================



! --- Examples (to be copied and adapted in section ''3.1. Parameters ...''

!-------------------------------------------------------------------------------
! This test allows the user to ensure that the version of this subroutine
!   used is that from his case definition, and not that from the library.

  if(0.eq.1) then

! --- Ex. 1: Perfect gas containing 3 components
!     Molar mass, gamma

    ! Molar mass of the components (kg/mol)
    cstgr(1)  = 18.d-3
    cstgr(2)  = 32.d-3
    cstgr(3)  = 28.d-3

    if(iccfth.gt.0) then

      ! Calculation of the molar mass of the mixture at cell centers
      do iel = 1, ncel
        xmasmr(iel) = 1.d0 / ( rtp(iel,isca(1))/cstgr(1)          &
                             + rtp(iel,isca(2))/cstgr(2)          &
                             + rtp(iel,isca(3))/cstgr(3) )
      enddo

      ! Calculation of the equivalent gamma of the mixture at cell centers
      do iel = 1, ncel
        gamagr(iel) = propce(iel,ipproc(icp))              &
           / ( propce(iel,ipproc(icp)) - rr/xmasmr(iel) )
      enddo

    endif

  endif

!-------------------------------------------------------------------------------

! End of the examples


! Verification of the values of gamagr: gamagr >= 1., otherwise stop

  ierr = 0

  do iel = 1, ncel
    if(iccfth.gt.0 .and. gamagr(iel).lt.1.d0) then
      ierr = 1
      write(nfecra,1020) iel, gamagr(iel)
    endif
  enddo

  if(ierr.eq.1) then
    call csexit (1)
  endif


! --- Calculation options: variable Cp and Cv
!     (isobaric and isochoric specific heat)

  if(iccfth.eq.-1) then

    icp = 1
    cp0 = epzero
    icv = 1
    cv0 = epzero


! Default initializations

  elseif(iccfth.eq.0) then

    do iel = 1, ncel
      propce(iel,ipproc(icp)) = cp0
      propce(iel,ipproc(icv)) =                            &
           cp0 - rr/xmasmr(iel)
      rtp(iel,irh) = p0*xmasmr(iel)/rr/t0
      rtp(iel,ien) = propce(iel,ipproc(icv))*t0
    enddo


! --- Calculation of temperature and energy from pressure and density

  elseif(iccfth.eq.12) then

    do iel = 1, ncel

      ! Temperature
      sorti1(iel) =                                               &
           xmasmr(iel)/rr*rtp(iel,ipr)/rtp(iel,irh)

      ! Total energy
      sorti2(iel) = propce(iel,ipproc(icv))*sorti1(iel)    &
    + 0.5d0*( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 )

    enddo

    ! Transfer to the array rtp
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,itk) = sorti1(iel)
        rtp(iel,ien) = sorti2(iel)
      enddo
    endif


! --- Calculation of density and energy from pressure and temperature:

  elseif(iccfth.eq.13) then

    do iel = 1, ncel

      ! Density
      sorti1(iel) =                                               &
           xmasmr(iel)/rr*rtp(iel,ipr)/rtp(iel,itk)

      ! Total energy
      sorti2(iel) =                                               &
           propce(iel,ipproc(icv))*rtp(iel,itk)         &
    + 0.5d0*( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 )

    enddo

    ! Transfer to the array rtp
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,irh) = sorti1(iel)
        rtp(iel,ien) = sorti2(iel)
      enddo
    endif


! --- Calculation of density and temperature from pressure and energy

  elseif(iccfth.eq.14) then

    do iel = 1, ncel

      ! Density
      sorti1(iel) =                                               &
           rtp(iel,ipr)/(gamagr(iel)-1.d0)/( rtp(iel,ien)   &
  - 0.5d0*( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 ) )

      ! Temperature
      sorti2(iel) = xmasmr(iel)/rr*rtp(iel,ipr)/sorti1(iel)

    enddo

    ! Transfer to the array rtp
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,irh) = sorti1(iel)
        rtp(iel,itk) = sorti2(iel)
      enddo
    endif


! --- Calculation of pressure and energy from density and temperature

  elseif(iccfth.eq.23) then

    do iel = 1, ncel

      ! Pressure
      sorti1(iel) =                                               &
           rtp(iel,irh)*rr/xmasmr(iel)*rtp(iel,itk)

      ! Total energy
      sorti2(iel) =                                               &
           propce(iel,ipproc(icv))*rtp(iel,itk)         &
    + 0.5d0*( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 )

    enddo

    ! Transfer to the array rtp
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,ipr) = sorti1(iel)
        rtp(iel,ien) = sorti2(iel)
      enddo
    endif


! --- Calculation of pressure and temperature from density and energy

  elseif(iccfth.eq.24) then

    do iel = 1, ncel

      ! Pressure
      sorti1(iel) =                                               &
           (gamagr(iel)-1.d0)*rtp(iel,irh)*( rtp(iel,ien)   &
  - 0.5d0*( rtp(iel,iu)**2 + rtp(iel,iv)**2 + rtp(iel,iw)**2 ) )

      ! Temperature
      sorti2(iel) = xmasmr(iel)/rr*sorti1(iel)/rtp(iel,irh)

    enddo

    ! Transfer to the array rtp
    if(imodif.gt.0) then
      do iel = 1, ncel
        rtp(iel,ipr) = sorti1(iel)
        rtp(iel,itk) = sorti2(iel)
      enddo
    endif

!                     2                            2         P
! --- Calculation of c from pressure and density: c = gamma*---
!                                                           rho

  elseif(iccfth.eq.126) then

    do iel = 1, ncel

      ! Verification of the positivity of the pressure
      if(rtp(iel,ipr).lt.0.d0) then
        write(nfecra,1110) iel , rtp(iel,ipr)
        ierr = 1

      ! Verification of the positivity of the density
      elseif(rtp(iel,irh).le.0.d0) then
        write(nfecra,1120) iel , rtp(iel,irh)
        ierr = 1

      else

        ! Computation
        sorti1(iel) =                                             &
             gamagr(iel) * rtp(iel,ipr) / rtp(iel,irh)

      endif

    enddo

    ! Stop if error detected
    if(ierr.eq.1) call csexit (1)


!                                                              gamma
! --- Calculation of beta from pressure and density: beta = rho

  elseif(iccfth.eq.162) then

    do iel = 1, ncel

      ! Verification of the positivity of the density
      if(rtp(iel,irh).lt.0.d0) then
        write(nfecra,1220) iel , rtp(iel,irh)
        ierr = 1

      else

        ! Computation
        sorti1(iel) = rtp(iel,irh)**gamagr(iel)

      endif

    enddo

    ! Stop if error detected
    if(ierr.eq.1) call csexit (1)


! --- Calculation of the isochoric specific heat: Cv = Cp - R/M

  elseif(iccfth.eq.432) then

    do iel = 1, ncel

      sorti1(iel) = propce(iel,ipproc(icp))-rr/xmasmr(iel)

    enddo

    ! Stop if error detected (kept by consistance with other sections)
    if(ierr.eq.1) call csexit (1)

!                                                                  P
! --- Calculation of the entropy from pressure and density: s = --------
!                                                                  gamma
!                                                               rho

  elseif(iccfth.eq.6) then

    do iel = 1, ncel

      ! Verification of the positivity of the pressure
      if(rtp(iel,ipr).lt.0.d0) then
        write(nfecra,1310) iel , rtp(iel,ipr)
        ierr = 1

      ! Verification of the positivity of the density
      elseif(rtp(iel,irh).le.0.d0) then
        write(nfecra,1320) iel , rtp(iel,irh)
        ierr = 1

      else

        ! Computation
        sorti1(iel) =                                             &
             rtp(iel,ipr) / (rtp(iel,irh)**gamagr(iel))

      endif

    enddo

    ! Stop if error detected
    if(ierr.eq.1) call csexit (1)


! --- Calculation of 'internal energy - Cv.T'

  elseif(iccfth.eq.7) then

    ! It is zero for a perfect gas

    !   At the cell centers
    do iel = 1, ncel
      sorti1(iel) = 0.d0
    enddo

    !   On the boundary faces
    do ifac = 1, nfabor
      sorti2(ifac) = 0.d0
    enddo

    ! Stop if error detected (kept by consistance with other sections)
    if(ierr.eq.1) call csexit (1)


! --- Calculation of the boundary conditions on the face ifac = ifac0

!  -- Wall/symmetry

  elseif(iccfth.eq.91) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Calculation of the Mach number at the boundary face, using the
    !   cell center velocity projected on the vector normal to the boundary
    xmach = ( rtp(iel,iu)*surfbo(1,ifac)                       &
           + rtp(iel,iv)*surfbo(2,ifac)                        &
           + rtp(iel,iw)*surfbo(3,ifac) ) / surfbn(ifac)       &
         / sqrt( gamagr(iel)*rtp(iel,ipr)/rtp(iel,irh) )

    coefa(ifac,iclp) = 0.d0

    ! Pression and entropy: rarefaction

    if(xmach.le.0.d0 .and. xmach.gt.2.d0/(1.d0-gamagr(iel))) then
      coefb(ifac,iclp) = (1.d0 + (gamagr(iel)-1.d0)/2.d0 * xmach) &
           ** (2.d0*gamagr(iel)/(gamagr(iel)-1.d0))
      coefb(ifac,iclt) = 1.d0

    elseif(xmach.le.2.d0/(1.d0-gamagr(iel)) ) then
      coefb(ifac,iclp) = 0.d0
      coefb(ifac,iclt) = 1.d0

      ! Pressure and entropy: shock

    else
      coefb(ifac,iclp) = 1.d0 + gamagr(iel)*xmach                 &
            *( (gamagr(iel)+1.d0)/4.d0*xmach                      &
           + sqrt(1.d0 + (gamagr(iel)+1.d0)**2/16.d0*xmach**2) )
      coefb(ifac,iclt) = coefb(ifac,iclp)/(1.d0-coefb(ifac,iclp)) &
          / rtp(iel,ipr) * ( rtp(iel,irh)                         &
              * (rtp(iel,iu)**2+rtp(iel,iv)**2+rtp(iel,iw)**2)    &
              + rtp(iel,ipr) *(1.d0-coefb(ifac,iclp)) )
    endif

    ! Total energy: 'internal energy - Cv T'

    coefa(ifac,icle) = 0.d0

    ! Stop if error detected
    if(ierr.eq.1) call csexit (1)


!  -- Inlet

  elseif(iccfth.eq.92) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Calculation of the Mach number at the boundary face, using the
    !   cell center velocity projected on the vector normal to the boundary
    xmachi = ( rtp(iel,iu)*surfbo(1,ifac)                      &
         + rtp(iel,iv)*surfbo(2,ifac)                          &
         + rtp(iel,iw)*surfbo(3,ifac) )/surfbn(ifac)           &
         / sqrt(gamagr(iel)*rtp(iel,ipr)/rtp(iel,irh))
    xmache = ( coefa(ifac,iclu)*surfbo(1,ifac)                    &
         + coefa(ifac,iclv)*surfbo(2,ifac)                        &
         + coefa(ifac,iclw)*surfbo(3,ifac) )/surfbn(ifac)         &
         / sqrt(gamagr(iel)*rtp(iel,ipr)/rtp(iel,irh))
    dxmach = xmachi - xmache

    ! Pressure: rarefaction wave
    if(dxmach.le.0.d0) then

      if(dxmach.gt.2.d0/(1.d0-gamagr(iel))) then
        coefa(ifac,iclp) = rtp(iel,ipr)*                       &
             ( (1.d0 + (gamagr(iel)-1.d0)*0.50d0*dxmach)          &
               ** (2.d0*gamagr(iel)/(gamagr(iel)-1.d0))  )
      elseif(dxmach.le.2.d0/(1.d0-gamagr(iel)) ) then
        coefa(ifac,iclp) = 0.d0
      endif

    ! Pressure: shock
    else
      coefa(ifac,iclp) = rtp(iel,ipr)*                         &
           (  1.d0 + gamagr(iel)*dxmach                           &
           *( (gamagr(iel)+1.d0)*0.25d0*dxmach                    &
           + sqrt(1.d0 + (gamagr(iel)+1.d0)**2/16.d0              &
                                           *dxmach**2) )  )
    endif

    ! This choice overrides the previous Rusanov choice
    coefa(ifac,iclp) = rtp(iel,ipr)

    ! Total energy
    coefa(ifac,icle) =                                            &
         coefa(ifac,iclp)/((gamagr(iel)-1.d0)*coefa(ifac,iclr))   &
         + 0.5d0*(coefa(ifac,iclu)**2                             &
                + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)

!  -- Outlet

  elseif(iccfth.eq.93) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Calculation of the Mach number at the boundary face, using the
    !   cell center velocity projected on the vector normal to the boundary
    xmach = ( rtp(iel,iu)*surfbo(1,ifac)                       &
           + rtp(iel,iv)*surfbo(2,ifac)                        &
           + rtp(iel,iw)*surfbo(3,ifac) ) / surfbn(ifac)       &
         / sqrt(gamagr(iel)*rtp(iel,ipr)/rtp(iel,irh))

    ! Supersonic outlet: Dirichlet for all variables
    if(xmach.ge.1.d0) then
      do ivar = 1, nvar
        coefa(ifac,iclrtp(ivar,icoef)) = rtp(iel,ivar)
      enddo

      ! Entropy
      coefa(ifac,iclt) =                                          &
           rtp(iel,ipr)/rtp(iel,irh)**gamagr(iel)

    ! Subsonic outlet
    elseif(xmach.lt.1.d0 .and. xmach.ge.0.d0) then

      ! Rarefaction:
      if(coefa(ifac,iclp).le.rtp(iel,ipr)) then

        ! Density
        coefa(ifac,iclr) = rtp(iel,irh)                        &
             * (coefa(ifac,iclp)/rtp(iel,ipr))                 &
                **(1.d0/gamagr(iel))

        ! Velocity
        coefa(ifac,iclu) = rtp(iel,iu)                         &
             + 2.d0/(gamagr(iel)-1.d0)                            &
 * sqrt( gamagr(iel) * rtp(iel,ipr) / rtp(iel,irh) )        &
             * ( 1.d0                                             &
 - (coefa(ifac,iclp)/rtp(iel,ipr))                             &
               **((gamagr(iel)-1.d0)/2.d0/gamagr(iel)) )          &
 * surfbo(1,ifac) / surfbn(ifac)

        coefa(ifac,iclv) = rtp(iel,iv)                         &
             + 2.d0/(gamagr(iel)-1.d0)                            &
 * sqrt( gamagr(iel) * rtp(iel,ipr) / rtp(iel,irh) )        &
             * ( 1.d0                                             &
 - (coefa(ifac,iclp)/rtp(iel,ipr))                             &
               **((gamagr(iel)-1.d0)/2.d0/gamagr(iel)) )          &
 * surfbo(2,ifac) / surfbn(ifac)

        coefa(ifac,iclw) = rtp(iel,iw)                         &
             + 2.d0/(gamagr(iel)-1.d0)                            &
 * sqrt( gamagr(iel) * rtp(iel,ipr) / rtp(iel,irh) )        &
             * ( 1.d0                                             &
 - (coefa(ifac,iclp)/rtp(iel,ipr))                             &
               **((gamagr(iel)-1.d0)/2.d0/gamagr(iel)) )          &
 * surfbo(3,ifac) / surfbn(ifac)

        ! Total energy
        coefa(ifac,icle) = coefa(ifac,iclp)                       &
 /( (gamagr(iel)-1.d0)*coefa(ifac,iclr) )                         &
              + 0.5d0*(coefa(ifac,iclu)**2                        &
                     + coefa(ifac,iclv)**2                        &
                     + coefa(ifac,iclw)**2)

        ! Entropy
        coefa(ifac,iclt) = coefa(ifac,iclp)                       &
                             /coefa(ifac,iclr)**gamagr(iel)

      ! Shock:
      else

        ! Density
        coefa(ifac,iclr) = rtp(iel,irh)                        &
 * ( (gamagr(iel)+1.d0)*coefa(ifac,iclp)                          &
   + (gamagr(iel)-1.d0)*rtp(iel,ipr) )                         &
 / ( (gamagr(iel)-1.d0)*coefa(ifac,iclp)                          &
   + (gamagr(iel)+1.d0)*rtp(iel,ipr) )

        ! Velocity
        coefa(ifac,iclu) = rtp(iel,iu)                         &
 - (coefa(ifac,iclp)-rtp(iel,ipr))*sqrt(2.d0/rtp(iel,irh)   &
 / ( (gamagr(iel)+1.d0)*coefa(ifac,iclp)                          &
   + (gamagr(iel)-1.d0)*rtp(iel,ipr) ))                        &
 * surfbo(1,ifac) / surfbn(ifac)

        coefa(ifac,iclv) = rtp(iel,iv)                         &
 - (coefa(ifac,iclp)-rtp(iel,ipr))*sqrt(2.d0/rtp(iel,irh)   &
 / ( (gamagr(iel)+1.d0)*coefa(ifac,iclp)                          &
   + (gamagr(iel)-1.d0)*rtp(iel,ipr) ))                        &
 * surfbo(2,ifac) / surfbn(ifac)

        coefa(ifac,iclw) = rtp(iel,iw)                         &
 - (coefa(ifac,iclp)-rtp(iel,ipr))*sqrt(2.d0/rtp(iel,irh)   &
 / ( (gamagr(iel)+1.d0)*coefa(ifac,iclp)                          &
   + (gamagr(iel)-1.d0)*rtp(iel,ipr) ))                        &
 * surfbo(3,ifac) / surfbn(ifac)

        ! Total energy
        coefa(ifac,icle) = coefa(ifac,iclp)                       &
 /( (gamagr(iel)-1.d0)*coefa(ifac,iclr) )                         &
     + 0.5d0*(coefa(ifac,iclu)**2                                 &
            + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)

        ! Entropy
        coefa(ifac,iclt) = coefa(ifac,iclp)                       &
                             /coefa(ifac,iclr)**gamagr(iel)

      endif

    else
      write(nfecra,*) 'iccfth = ',iccfth,'  Mach = ',xmach
      ierr = 1
    endif

    if(ierr.eq.1) call csexit (1)


! --- Calculation of temperature and energy from pressure and density

  elseif(iccfth.eq.912.or.iccfth.eq.60900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Temperature
    coefa(ifac,iclt) = xmasmr(iel)/rr*coefa(ifac,iclp)            &
                                        /coefa(ifac,iclr)

    ! Total energy
    coefa(ifac,icle) = propce(iel,ipproc(icv))             &
               * coefa(ifac,iclt) + 0.5d0*( coefa(ifac,iclu)**2   &
                    + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)


! --- Calculation of density and energy from pressure and temperature

  elseif(iccfth.eq.913.or.iccfth.eq.100900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Density
    coefa(ifac,iclr) = xmasmr(iel)/rr*coefa(ifac,iclp)            &
                                       /coefa(ifac,iclt)

    ! Total energy
    coefa(ifac,icle) = propce(iel,ipproc(icv))             &
               * coefa(ifac,iclt) + 0.5d0*( coefa(ifac,iclu)**2   &
                    + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)


! --- Calculation of density and temperature from pressure and total energy

  elseif(iccfth.eq.914.or.iccfth.eq.140900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Density
    coefa(ifac,iclr) = coefa(ifac,iclp)/(gamagr(iel)-1.d0)        &
           / (coefa(ifac,icle) - 0.5d0*( coefa(ifac,iclu)**2      &
                 + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2 ))

    ! Temperature
    coefa(ifac,iclt)= xmasmr(iel)/rr*coefa(ifac,iclp)             &
                                       /coefa(ifac,iclr)


! --- Calculation of pressure and energy from density and temperature

  elseif(iccfth.eq.923.or.iccfth.eq.150900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Pressure
    coefa(ifac,iclp) = coefa(ifac,iclr)*rr/xmasmr(iel)            &
                                       *coefa(ifac,iclt)

    ! Total energy
    coefa(ifac,icle) = propce(iel,ipproc(icv))             &
               * coefa(ifac,iclt) + 0.5d0*( coefa(ifac,iclu)**2   &
                    + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2)


! --- Calculation of pressure and temperature from density and energy

  elseif(iccfth.eq.924.or.iccfth.eq.210900) then

    ifac = ifac0
    iel  = ifabor(ifac)

    ! Pressure
    coefa(ifac,iclp) = (gamagr(iel)-1.d0)*coefa(ifac,iclr)        &
          *( coefa(ifac,icle) - 0.5d0*( coefa(ifac,iclu)**2       &
                + coefa(ifac,iclv)**2 + coefa(ifac,iclw)**2 ) )


    ! Temperature
    coefa(ifac,iclt)= xmasmr(iel)/rr*coefa(ifac,iclp)             &
                                       /coefa(ifac,iclr)


! --- End of perfect gas with variable gamma
  endif

! --- End of test on the thermodynamic laws
endif


!--------
! Formats
!--------

 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     Gamma = ',e12.4   ,/,                                      &
'@     Gamma must be a real number greater or equal to 1.',/,     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1020 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     In cell ',i10   ,', Gamma = ',e12.4   ,/,                  &
'@     Gamma must be a real number greater or equal to 1.',/,     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 2010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The computation of density failed.',/,                     &
'@',/,                                                            &
'@     Temperature = ',e12.4   ,' in cell ',i10  ,/,              &
'@     Temperature must be strictly positive.',/,                 &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 3010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The computation of temperature failed.',/,                 &
'@',/,                                                            &
'@     Density = ',e12.4   ,' in cell ',i10  ,/,                  &
'@     Density must be strictly positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 4010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the squared speed of sound failed.',/,  &
'@',/,                                                            &
'@     Density = ',e12.4   ,' in cell ',i10  ,/,                  &
'@     Density must be strictly positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 4020 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the variable beta failed.',/,           &
'@',/,                                                            &
'@     Density = ',e12.4   ,' in cell ',i10  ,/,                  &
'@     Density must be strictly positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 4030 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the entropy failed.',/,                 &
'@',/,                                                            &
'@     Density = ',e12.4   ,' in cell ',i10  ,/,                  &
'@     Density must be strictly positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 7000 format (                                                    &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     The boundary condition of the type ''prescribed mass',/,   &
'@     and enthalpy flow rates '' is not available in the ',/,    &
'@     current release.',/,                                       &
'@',/,                                                            &
'@     Modify the user subroutine ''uscfth''.',/,                 &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 8000 format (                                                    &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     Negative values of the density were encountered ',/,       &
'@     in ',i10   ,' cells.',/,                                   &
'@     The density was clipped at ',e12.4  ,/                     &
'@     The run was stopped.',/,                                   &
'@',/,                                                            &
'@     If it is desired to continue the run in spite of this ',/, &
'@     behavior, it is possible to force a standard clipping ',/, &
'@     by setting a minimum value for the density variable in',/, &
'@     the GUI or in the user subroutine ''usini1'' (set the ',/, &
'@     scamin value associated to the variable ',/,               &
'@     isca(irho).',/,                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 8100 format (                                                    &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with constant gamma.',/,                 &
'@',/,                                                            &
'@     Negative values of the internal energy were encountered',/,&
'@     in ',i10   ,' cells.',/,                                   &
'@     The internal energy  was clipped at ',e12.4  ,/            &
'@     The run was stopped.',/,                                   &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)


! The following formats may be discarded if or when the
! gamma variable option will have been fixed


 1110 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with variable gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the squared speed of sound failed.',/,  &
'@',/,                                                            &
'@     In cell ',i10   ,' Pressure = ',e12.4   ,/,                &
'@     Pressure must be positive.',/,                             &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1120 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with variable gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the squared speed of sound failed.',/,  &
'@',/,                                                            &
'@     In cell ',i10   ,' Density = ',e12.4   ,/,                 &
'@     Density must be strictly positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1220 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with variable gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the variable beta failed.',/,           &
'@',/,                                                            &
'@     In cell ',i10   ,' Density = ',e12.4   ,/,                 &
'@     Density must be strictly positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1310 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with variable gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the entropy failed.',/,                 &
'@',/,                                                            &
'@     In cell ',i10   ,' Pressure = ',e12.4   ,/,                &
'@     Pressure must be positive.',/,                             &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1320 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in thermodynamics computations',/,         &
'@    =======',/,                                                 &
'@     Error encountered in the user subroutine ''uscfth'', ',/,  &
'@       for perfect gas with variable gamma.',/,                 &
'@',/,                                                            &
'@     The computation of the entropy failed.',/,                 &
'@',/,                                                            &
'@     In cell ',i10   ,' Density = ',e12.4   ,/,                 &
'@     Density must be striclty positive.',/,                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)


!----
! End
!----

return
end subroutine
