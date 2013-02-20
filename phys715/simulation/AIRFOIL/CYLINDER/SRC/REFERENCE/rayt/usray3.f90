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

subroutine usray3 &
!================

 ( nvar   , nscal  , iappel ,                                     &
   itypfb ,                                                       &
   izfrdp ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   ck     )

!===============================================================================
! Purpose:
! --------

! Absorption coefficient for radiative module
! ----------------------

! It is necessary to define the value of the fluid's absorption
! coefficient Ck.

! For a transparent medium, the coefficient should be set to 0.d0.

! In the case of the P-1 model, we check that the optical length is at
! least of the order of 1.

! Caution:
! ========
!   For specific physics (Combustion, coal, ...),

!   it is Forbidden to define the absorption coefficient here.
!         =========

!   See subroutine ppcabs.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! izfrdp(nfabor    ! ia ! <-- ! zone number for boundary faces                 !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! ck(ncelet)       ! ra ! --> ! medium's absorption coefficient                !
!                  !    !     ! (zero if transparent)                          !
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
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use cpincl
use ppincl
use radiat
use ihmpre
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal  , iappel

integer          itypfb(nfabor)
integer          izfrdp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)


double precision ck(ncelet)

! Local variables

integer          iel, ifac, iok
double precision vv, sf, xlc, xkmin, pp

!===============================================================================



!===============================================================================
! 0 - Memory management
!===============================================================================


! Stop flag (to determine if faces were forgotten)
iok = 0

!===============================================================================
! Absorption coefficient of the medium (m-1)

! In the case of specific physics (gas/coal/fuel combustion, elec)

! Ck must not be defined here
!============
! (it is determined automatically, possibly from the parametric file)

! In other cases, Ck must be defined (it is zero by default)
!===============================================================================

if (ippmod(iphpar).le.1) then

  do iel = 1, ncel
    ck(iel) = 0.d0
  enddo

  !--> P1 model: standard control of absorption coefficient values.
  !              this coefficient must ensure an optical length
  !              at least of the order of 1.

  if (iirayo.eq.2) then
    sf = 0.d0
    vv = 0.d0

    ! Compute characteristic length of calculation domain

    do ifac = 1,nfabor
      sf = sf + sqrt(surfbo(1,ifac)**2 + surfbo(2,ifac)**2 + surfbo(3,ifac)**2)
    enddo
    if (irangp.ge.0) then
      call parsom(sf)
      !==========
    endif

    do iel = 1,ncel
      vv = vv + volume(iel)
    enddo
    if (irangp.ge.0) then
      call parsom(vv)
      !==========
    endif

    xlc = 3.6d0 * vv / sf

    ! Clipping for variable CK

    xkmin = 1.d0 / xlc

    iok = 0

    do iel = 1,ncel
      if (ck(iel).lt.xkmin) then
        iok = iok + 1
      endif
    enddo

    ! Stop at the end of the time step if the optical thickness is too big
    ! (istpp1 = 1 allows stopping cleanly at the end of the current time step).
    pp = xnp1mx/100.0d0
    if (dble(iok).gt.pp*dble(ncel)) then
      write(nfecra,3000) xkmin, dble(iok)/dble(ncel)*100.d0, xnp1mx
      istpp1 = 1
      ! call csexit(1)
      !==========
    endif
  endif

endif

! -------
! Formats
! -------

 3000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:    P1 radiation approximation (usray3)         ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    The optical length of the semi-transparent medium       ',/,&
'@      must be at least of the order of one to be in the     ',/,&
'@      domain of validity of the P-1 approximation.          ',/,&
'@    This does not seem to be the case here.                 ',/,&
'@                                                            ',/,&
'@    The minimum absorption coefficient to ensure this       ',/,&
'@      optical length is XKmin = ', e10.4                     ,/,&
'@    This value is not reached for ', e10.4,'%               ',/,&
'@      of the meshe''s cells.                                ',/,&
'@    The percentage of mesh cells for which we allow this    ',/,&
'@      condition not to be rspected is set by default in     ',/,&
'@      usini1 to xnp1mx = ', e10.4,'%                        ',/,&
'@                                                            ',/,&
'@    The calculation is interrupted.                         ',/,&
'@                                                            ',/,&
'@    Check the values of the absorption coefficient Ck       ',/,&
'@      in ppcabs, usray3 or the thermochemistry data file.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine
