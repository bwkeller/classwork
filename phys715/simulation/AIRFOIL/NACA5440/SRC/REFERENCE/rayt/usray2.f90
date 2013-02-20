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

subroutine usray2 &
!================

 ( nvar   , nscal  ,                                              &
   itypfb ,                                                       &
   icodcl , izfrdp , isothp ,                                     &
   tmin   , tmax   , tx     ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  ,                                              &
   tparop , qincid , hfcnvp , flcnvp ,                            &
   xlamp  , epap   , epsp   , textp  , tintp  )

!===============================================================================
! Purpose:
! --------

! User subroutine for input of radiative transfer parameters: boundary conditions

! Sketch of thermal flux in boundary wall


!  Outside      Inside wall       Fluid domain
!            |               |
!            |               |
!            |               |
!            | EPSP Qincid<--|<----- Qincid             [Incident flux]
!            |               |
!            | [absorption]  |
!            |               |
!            |               |-----> (1-EPSP) Qincid    [reflexion]
!            |               |
!            |               |
!            |               |                       4
!            |               |-----> EPSP sigma Tparop  [emission]
!            |               |
!            |               |
!            |               |
!            | [conduction]  |     .....o Tfluide       [convection]
!            |               |    .
!            |     XLAMP     |   .    Hfluide
!            |               |  .
!            |            ...o..
!            |          ..   | Tparop (Inside wall temperature)
!            |         .     |
!            |       ..      |
!            |      .        |
!            |    ..         |
!            |   .           |
!            o...            |
!       Textp|               |
!(Outside    |               |
! wall       |<------------->|
! temprature)|     EPAP      |
!            |               |


!  The radiative boundary condition ios based on the calculation of
!  a new wall temperature.
!  This temparature is computed with a thermal flux balance:

!  Q           =  Q           + (Qrayt           - Qrayt        )
!   conduction     convection         absorption        emission

!  therfore:

!   XLAMP
!   -----(Tparop-Textp) =
!   EPAP
!                                                             4
!         Hfluide (Tfluide-Tparop) + EPSP (QINCID - SIGMA Tparop )

!  Note: in Code_Saturne the flux is positive when it is oriented
!  from inside to outside.
!
!
!      CORPS                     Emissivity
!      ------------------------------------
!      polished steel               0,06
!      oxidized steel               0,80
!      steel rough                  0,94
!      polished aluminium           0,04
!      oxidized aluminium (inside)  0,09
!      oxidized aluminium (wet air) 0,90
!      brick                        0,93
!      concrete                     0,93
!      paper                        0,8 to 0,9
!      water                        0,96


! Boundary faces identification
! =============================

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.

! Note: these usefull constants are defined
!       TKELVI = 273.16D0
!       SIG    = 5.6703D-8

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypfb           ! ia ! <-- ! boundary face types                            !
! icodcl           ! ia ! <-- ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1  -> Dirichlet                              !
!                  !    !     ! = 2  -> flux density                           !
!                  !    !     ! = 4  -> sliding wall and u.n=0 (velocity)      !
!                  !    !     ! = 5  -> friction and u.n=0 (velocity)          !
!                  !    !     ! = 6  -> roughness and u.n=0 (velocity)         !
!                  !    !     ! = 9  -> free inlet/outlet (velocity)           !
!                  !    !     !         inflowing possibly blocked             !
! izfrdp(nfabor)   ! ia ! --> ! boundary faces -> zone number                  !
! isothp(nfabor)   ! ia ! --> ! boundary face type for radative transfer       !
!                  !    !     ! = itpimp -> Gray wall with fixed inside temp   !
!                  !    !     ! = ipgrno -> Gray wall with fixed outside temp  !
!                  !    !     ! = iprefl -> Reflecting wall with fixed         !
!                  !    !     !                                  outside temp  !
!                  !    !     ! = ifgrno -> Gray wall with fixed               !
!                  !    !     !                               conduction flux  !
!                  !    !     ! = ifrefl -> Reflecting wall with fixed         !
!                  !    !     !                               conduction flux  !
! tmin             ! r  !     ! min value of the wall temperature              !
! tmax             ! r  !     ! max value of the wall temperature              !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! rcodcl           ! ra ! --> ! boundary condition values                      !
!                  !    !     ! rcodcl(3) = flux density value                 !
!                  !    !     !  (negative for gain) in w/m2                   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! tparop(nfabor)   ! ra ! <-- ! inside current wall temperature (K)            !
! qincid(nfabor)   ! ra ! <-- ! radiative incident flux  (W/m2)                !
! hfcnvp(nfabor)   ! ra ! <-- ! convective exchange coefficient (W/m2/K)       !
! flcnvp(nfabor)   ! ra ! <-- ! convective flux (W/m2)                         !
! xlamp(nfabor)    ! ra ! --> ! conductivity (W/m/K)                           !
! epap(nfabor)     ! ra ! --> ! thickness (m)                                  !
! epsp(nfabor)     ! ra ! --> ! emissivity (>0)                                !
! textp(nfabor)    ! ra ! --> ! outside temperature (K)                        !
! tintp(nfabor)    ! ra ! --> ! initial inside temperature (K)                 !
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
use radiat
use ihmpre
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          itypfb(nfabor)

integer          icodcl(nfabor,nvar)
integer          izfrdp(nfabor), isothp(nfabor)

double precision tmin , tmax , tx

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision coefa(nfabor,*), coefb(nfabor,*)


double precision tparop(nfabor), qincid(nfabor)
double precision hfcnvp(nfabor),flcnvp(nfabor)
double precision xlamp(nfabor), epap(nfabor)
double precision epsp(nfabor)
double precision textp(nfabor), tintp(nfabor)

! Local variables

integer          ifac , ivar, iok
integer          ilelt, nlelt

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

!===============================================================================
! 0. Initialization
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))


!===============================================================================
! 1. IVAR: number of the thermal variable
!===============================================================================

ivar = isca(iscalt)

!===============================================================================
!  2. Min and max values for the wall temperatures (clipping otherwise)
!   TMIN and TMAX are given in Kelvin.
!===============================================================================

tmin = 0.d0
tmax = grand + tkelvi

!===============================================================================
! 3. Assign boundary conditions to boundary wall
!===============================================================================

!     ZONES DEFINITION
!     ================

!     We define zones of wall boundary, and we assign a type.
!       This allows to apply the boundary conditions and realize
!       balance sheets by treating them separately for each zone.

!     For each boundary face ifac (not just the faces of wall)
!       the user defines his own choice by a number of zone
!       IZFRDP(ifac) from color of the boundary face
!         or more generally, their properties (color, groups ...),
!         or boundary conditions specified in usclim,
!         or even of their coordinates.
!     Warning: it is essential that ALL boundary faces
!       have been assigned to a zone.
!     The number of zones (the value of IZFRDP(ifac)) is
!       arbitrarily chosen by the user, but must be a
!       positive integer and less than or equal to NBZRDM
!       (value set in parameter radiat.h).



!     WALL CARACTERISTICS
!     ===================

!      WARNING: the unity of the temperature is the Kelvin
!      -------

!      Mandatory data:
!      ---------------
!      isothp(ifac) boundary face type
!                  = itpimp -> Gray wall with fixed inside temperature
!                  = ipgrno -> Gray wall with fixed outside temperature
!                  = iprefl -> Reflecting wall with fixed outside temperature
!                  = ifgrno -> Gray wall with fixed conduction flux
!                  = ifrefl -> Reflecting wall with fixed conduction flux

!      tintp(ifac) inside wall temperature (Kelvin)
!                  initialize tparop at the first time step.
!                  If isothp = itpimp, the value of tparop is fixed to tintp
!                  In the other case, tintp is only for initialization.


!      Other data (depend of the isothp):
!      ----------------------------------

!      rcodcl = conduction flux
!      epsp   = emissivity
!      xlamp  = conductivity (W/m/K)
!      epap   = thickness (m)
!      textp  = outside temperature (K)


!     EXAMPLE
!     =======

!        Wall boundary faces (IPAROI and IPARUG), are devided into 5 zones
!          located with IFRFAC(IFAC) in the range of number from 51 to 55.
!          For each location a different radiative boundary condition is applied.
!        For all other boundary that are not wall (i.e. inlet, oulet, symetry)
!          the user can define arbritay new zone using the array IFRFAC(IFAC),
!          for wich a value can be arbitrarily choosen between 1 and NBZRDM.
!
!     Warning: it is forbidden to modify tparop and qincid in this subroutine
!     ========

!    Indicator for forgotten faces.
iok = 0

!   -------------------------------------------------------------------
!-->  Example 1:
!       For wall boundary faces, selection criteria: color 1
!       Gray or black wall with profil of fixed inside temperature
!       ------------------------------------

CALL GETFBR('1',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if ( itypfb(ifac).eq.iparoi ) then

!      zone number
    izfrdp(ifac) = 51

!      Type of condition: gray or black wall with fixed inside temperature
    isothp(ifac) = itpimp

!      Emissivity
    epsp  (ifac) = 0.1d0

!      Profil of fixed inside temperature
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo

!   -------------------------------------------------------------------
!-->  Example 2 :
!       For wall boundary faces, selection criteria: color 2
!       Gray or black wall with fixed outside temperature TEXTP
!       ------------------------------------

CALL GETFBR('2',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if ( itypfb(ifac).eq.iparug ) then

!      zone number
    izfrdp(ifac) = 52

!      Type of condition: gray or black wall with fixed outside temperature TEXTP
    isothp(ifac) = ipgrno

!        Emissivity
    epsp  (ifac) = 0.9d0
!        Conductivity (W/m/K)
    xlamp (ifac) = 3.0d0
!        Thickness    (m)
    epap  (ifac) = 0.1d0
!        Fixed outside temperature: 473.16 K
    textp (ifac) = 200.d0 + tkelvi
!        Initial inside temperature: 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo

!   -------------------------------------------------------------------
!-->  Exemple 3 :
!       For wall boundary faces, selection criteria: color 3
!       Reflecting wall (EPSP = 0) with fixed outside temperature TEXTP
!       ------------------------------------

CALL GETFBR('3',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if ( itypfb(ifac).eq.iparoi ) then

!      zone number
    izfrdp(ifac) = 53

!      Type of condition: reflecting wall with fixed outside temperature TEXTP
    isothp(ifac) = iprefl

!        Conductivity (W/m/K)
    xlamp (ifac) = 3.0d0
!        Thickness    (m)
    epap  (ifac) = 0.1d0
!        Fixed outside temperature: 473.16 K
    textp (ifac) = 200.d0 + tkelvi
!        Initial inside temperature: 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo

!   -------------------------------------------------------------------
!-->  Example 4 :
!      For wall boundary faces which have the color 4:
!           gray or black wall and fixed conduction flux through the wall

!        XLAMP
!        -----(Tparop-Textp) = fixed conduction flux     (W/m2)
!        EPAP
!                         = RODCL(IFAC,IVAR,3)

!       If the conduction flux is zero then the wall is adiabatic.
!       The array RCODCL(IFAC,IVAR,3) has the value of the flux.
!       Flux density (< 0 if gain for the fluid)
!         For temperatures T,    in Watt/m2:
!            RCODCL(IFAC,IVAR,3) = CP*(VISCLS+VISCT/SIGMAS) * GRAD T
!         For enthalpies H,      in Watt/m2:
!            RCODCL(IFAC,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD H
!       ------------------------------------

CALL GETFBR('4',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if ( itypfb(ifac).eq.iparoi ) then

!      zone number
    izfrdp(ifac) = 54

!      Type of condition: gray or black wall with fixed conduction flux through the wall
    isothp(ifac) = ifgrno

!      Emissivity
    epsp  (ifac) = 0.9d0
!      Conduction flux (W/m2)
    rcodcl(ifac,ivar,3) = 0.d0
!      Initial inside temperature: 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo

!   -------------------------------------------------------------------
!-->  Example 5 :
!      Pour les faces PAROI de couleur 5 :
!           Paroi reflechissante (EPSP = 0) et
!           Flux de conduction impose dans la paroi
!      For wall boundary faces which have the color 5:
!           reflecting wall and fixed conduction flux through the wall

!      Equivalent a imposer une condition de flux au fluide

!        XLAMP
!        -----(Tparop-Textp) = fixed conduction flux and EPSP = 0
!        EPAP
!                         = RODCL(IFAC,IVAR,3)

!       If the conduction flux is zero then the wall is adiabatic.
!       Flux density (< 0 if gain for the fluid)
!         For temperatures T,    in Watt/m2:
!            RCODCL(IFAC,IVAR,3) = CP*(VISCLS+VISCT/SIGMAS) * GRAD T
!         For enthalpies H,      in Watt/m2:
!            RCODCL(IFAC,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD H
!       ------------------------------------

CALL GETFBR('5',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if ( itypfb(ifac).eq.iparoi ) then

!      zone number
    izfrdp(ifac) = 55

!      Type of condition: reflecting wall with fixed conduction flux through the wall
    isothp(ifac) = ifrefl

!      Conduction flux (W/m2)
    rcodcl(ifac,ivar,3) = 0.d0
!      Initial inside temperature: 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo

!                           WARNING

!   -------------------------------------------------------------------
!-->   For all boundary faces that are not wall it is MANDATORY to
!      impose a number of zone in the array izfrdp.
!      For each zone, informations will be displayed in the listing.
!       ------------------------------------

do ifac = 1, nfabor

  if     ( itypfb(ifac).eq.isolib                  ) then
    izfrdp(ifac) = 61
  elseif ( itypfb(ifac).eq.ientre.and.                      &
           cdgfbo(2,ifac)    .gt.0.d0                    ) then
    izfrdp(ifac) = 62
  elseif ( itypfb(ifac).eq.ientre.and.                      &
           cdgfbo(2,ifac)    .le.0.d0                    ) then
    izfrdp(ifac) = 63
  elseif ( itypfb(ifac).eq.isymet                  ) then
    izfrdp(ifac) = 64



!   -------------------------------------------------------------------
!-->  Example 7 :
!      Verification that all boundary faces have been treated.
!       ------------------------------------

  elseif ( itypfb(ifac).eq.iparoi .or.                      &
           itypfb(ifac).eq.iparug     ) then
    if (izfrdp(ifac) .eq. -1) then
      write(nfecra,1000)ifac
      iok = iok + 1
    endif
  endif




!     End of the loop on the boundary faces
!     -------------------------------------

enddo

! Stop if there are forgotten faces
if(iok.ne.0) then
  call csexit (1)
  !==========
endif

! -------
! FORMAT
! -------


 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in definition of boundary conditions',/,   &
'@    =======',/,                                                 &
'@   Radiative data are missing for face: ',I10,/,                &
'@',/,                                                            &
'@     The user subroutine ''usray2'' must be completed.',/, &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

! ---
! END
! ---

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
