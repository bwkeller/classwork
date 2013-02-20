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

subroutine user_fuel_bcond &
!================

 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl )

!===============================================================================
! PURPOSE  :
! --------

!    USER SUBROUTINE for extended physic
!                Combsution of heavy Fuel oil
!    Allocation of boundary conditions (ICODCL, RCODCL)
!    for variables unknowns during USCLIM


! Introduction
! ============

! Here we define boundary conditions on a per-face basis.

! Boundary faces may be identified using the 'getfbr' subroutine.

!  getfbr(string, nelts, eltlst) :
!  - string is a user-supplied character string containing
!    selection criteria;
!  - nelts is set by the subroutine. It is an integer value
!    corresponding to the number of boundary faces verifying the
!    selection criteria;
!  - lstelt is set by the subroutine. It is an integer array of
!    size nelts containing the list of boundary faces verifying
!    the selection criteria.

!  string may contain:
!  - references to colors (ex.: 1, 8, 26, ...
!  - references to groups (ex.: inlet, group1, ...)
!  - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!  These criteria may be combined using logical operators
!  ('and', 'or') and parentheses.
!  Example: '1 and (group2 or group3) and y < 1' will select boundary
!  faces of color 1, belonging to groups 'group2' or 'group3' and
!  with face center coordinate y less than 1.



! Boundary condition types
! ========================

! Boundary conditions may be assigned in two ways.


!    For "standard" boundary conditions:
!    -----------------------------------

!     (inlet, free outlet, wall, symmetry), we define a code
!     in the 'itypfb' array (of dimensions number of boundary faces,
!     number of phases). This code will then be used by a non-user
!     subroutine to assign the following conditions (scalars in
!     particular will receive the conditions of the phase to which
!     they are assigned). Thus:

!     Code      |  Boundary type
!     --------------------------
!      ientre   |   Inlet
!      isolib   |   Free outlet
!      isymet   |   Symmetry
!      iparoi   |   Wall (smooth)
!      iparug   |   Rough wall

!     Integers ientre, isolib, isymet, iparoi, iparug
!     are defined elsewhere (param.h). Their value is greater than
!     or equal to 1 and less than or equal to ntypmx
!     (value fixed in paramx.h)


!     In addition, some values must be defined:


!     - Inlet (more precisely, inlet/outlet with prescribed flow, as
!              the flow may be prescribed as an outflow):

!       -> Dirichlet conditions on variables
!         other than pressure are mandatory if the flow is incoming,
!         optional if the flow is outgoing (the code assigns 0 flux
!         if no Dirichlet is specified); thus,
!         at face 'ifac', for the variable 'ivar': rcodcl(ifac, ivar, 1)


!     - Smooth wall: (= impermeable solid, with smooth friction)

!       -> Velocity value for sliding wall if applicable
!         at face ifac, rcodcl(ifac, iu, 1)
!                       rcodcl(ifac, iv, 1)
!                       rcodcl(ifac, iw, 1)
!       -> Specific code and prescribed temperature value
!         at wall, if applicable:
!         at face ifac, icodcl(ifac, ivar)    = 5
!                       rcodcl(ifac, ivar, 1) = prescribed temperature
!       -> Specific code and prescribed flux value
!         at wall, if applicable:
!         at face ifac, icodcl(ifac, ivar)    = 3
!                       rcodcl(ifac, ivar, 3) = prescribed flux
!                                        =
!        Note that the default condition for scalars
!         (other than k and epsilon) is homogeneous Neumann.


!     - Rough wall: (= impermeable solid, with rough friction)

!       -> Velocity value for sliding wall if applicable
!         at face ifac, rcodcl(ifac, iu, 1)
!                       rcodcl(ifac, iv, 1)
!                       rcodcl(ifac, iw, 1)
!       -> Value of the dynamic roughness height to specify in
!                       rcodcl(ifac, iu, 3) (value for iv et iw not used)
!       -> Specific code and prescribed temperature value
!         at rough wall, if applicable:
!         at face ifac, icodcl(ifac, ivar)    = 6
!                       rcodcl(ifac, ivar, 1) = prescribed temperature
!                       rcodcl(ifac, ivar, 3) = dynamic roughness height
!       -> Specific code and prescribed flux value
!         at rough wall, if applicable:
!         at face ifac, icodcl(ifac, ivar)    = 3
!                       rcodcl(ifac, ivar, 3) = prescribed flux
!                                        =
!        Note that the default condition for scalars
!         (other than k and epsilon) is homogeneous Neumann.

!     - Symmetry (= impermeable frictionless wall):

!       -> Nothing to specify


!     - Free outlet (more precisely free inlet/outlet with prescribed pressure)

!       -> Nothing to prescribe for pressure and velocity
!          For scalars and turbulent values, a Dirichlet value may optionally
!            be specified. The behavior is as follows:
!              * pressure is always handled as a Dirichlet condition
!              * if the mass flow is inflowing:
!                  we retain the velocity at infinity
!                  Dirichlet condition for scalars and turbulent values
!                    (or zero flux if the user has not specified a
!                    Dirichlet value)
!                if the mass flow is outflowing:
!                  we prescribe zero flux on the velocity, the scalars,
!                  and turbulent values

!       Note that the pressure will be reset to P0
!           on the first free outlet face found


!    For "non-standard" conditions:
!    ------------------------------

!     Other than (inlet, free outlet, wall, symmetry), we define
!      - on one hand, for each face:
!        -> an admissible 'itypfb' value
!           (i.e. greater than or equal to 1 and less than or equal to
!            ntypmx; see its value in paramx.h).
!           The values predefined in paramx.h:
!           'ientre', 'isolib', 'isymet', 'iparoi', 'iparug' are in
!           this range, and it is preferable not to assign one of these
!           integers to 'itypfb' randomly or in an inconsiderate manner.
!           To avoid this, we may use 'iindef' if we wish to avoid
!           checking values in paramx.h. 'iindef' is an admissible
!           value to which no predefined boundary condition is attached.
!           Note that the 'itypfb' array is reinitialized at each time
!           step to the non-admissible value of 0. If we forget to
!           modify 'typfb' for a given face, the code will stop.

!      - and on the other hand, for each face and each variable:
!        -> a code             icodcl(ifac, ivar)
!        -> three real values  rcodcl(ifac, ivar, 1)
!                              rcodcl(ifac, ivar, 2)
!                              rcodcl(ifac, ivar, 3)
!     The value of 'icodcl' is taken from the following:
!       1: Dirichlet      (usable for any variable)
!       3: Neumann        (usable for any variable)
!       4: Symmetry       (usable only for the velocity and
!                          components of the Rij tensor)
!       5: Smooth wall    (usable for any variable except for pressure)
!       6: Rough wall     (usable for any variable except for pressure)
!       9: Free outlet    (usable only for velocity)
!     The values of the 3 'rcodcl' components are
!      rcodcl(ifac, ivar, 1):
!         Dirichlet for the variable          if icodcl(ifac, ivar) =  1
!         wall value (sliding velocity, temp) if icodcl(ifac, ivar) =  5
!         The dimension of rcodcl(ifac, ivar, 1) is that of the
!           resolved variable: ex U (velocity in m/s),
!                                 T (temperature in degrees)
!                                 H (enthalpy in J/kg)
!                                 F (passive scalar in -)
!      rcodcl(ifac, ivar, 2):
!         "exterior" exchange coefficient (between the prescribed value
!                          and the value at the domain boundary)
!                          rinfin = infinite by default
!         For velocities U,                in kg/(m2 s):
!           rcodcl(ifac, ivar, 2) =          (viscl+visct) / d
!         For the pressure P,              in  s/m:
!           rcodcl(ifac, ivar, 2) =                     dt / d
!         For temperatures T,              in Watt/(m2 degres):
!           rcodcl(ifac, ivar, 2) = Cp*(viscls+visct/sigmas) / d
!         For enthalpies H,                in kg /(m2 s):
!           rcodcl(ifac, ivar, 2) =    (viscls+visct/sigmas) / d
!         For other scalars F              in:
!           rcodcl(ifac, ivar, 2) =    (viscls+visct/sigmas) / d
!              (d has the dimension of a distance in m)
!
!      rcodcl(ifac, ivar, 3) if icodcl(ifac, ivar) <> 6:
!        Flux density (< 0 if gain, n outwards-facing normal)
!                         if icodcl(ifac, ivar)= 3
!         For velocities U,                in kg/(m s2) = J:
!           rcodcl(ifac, ivar, 3) =         -(viscl+visct) * (grad U).n
!         For pressure P,                  en kg/(m2 s):
!           rcodcl(ifac, ivar, 3) =                    -dt * (grad P).n
!         For temperatures T,              in Watt/m2:
!           rcodcl(ifac, ivar, 3) = -Cp*(viscls+visct/sigmas) * (grad T).n
!         For enthalpies H,                in Watt/m2:
!           rcodcl(ifac, ivar, 3) = -(viscls+visct/sigmas) * (grad H).n
!         For other scalars F in :
!           rcodcl(ifac, ivar, 3) = -(viscls+visct/sigmas) * (grad F).n

!      rcodcl(ifac, ivar, 3) if icodcl(ifac, ivar) = 6:
!        Roughness for the rough wall law
!         For velocities U, dynamic roughness
!           rcodcl(ifac, ivar, 3) = rugd
!         For other scalars, thermal roughness
!           rcodcl(ifac, ivar, 3) = rugt


!      Note that if the user assigns a value to itypfb equal to
!       ientre, isolib, isymet, iparoi, or iparug
!       and does not modify icodcl (zero value by default),
!       itypfb will define the boundary condition type.

!      To the contrary, if the user prescribes
!        icodcl(ifac, ivar) (nonzero),
!        the values assigned to rcodcl will be used for the considered
!        face and variable (if rcodcl values are not set, the default
!        values will be used for the face and variable, so:
!                                 rcodcl(ifac, ivar, 1) = 0.d0
!                                 rcodcl(ifac, ivar, 2) = rinfin
!                                 rcodcl(ifac, ivar, 3) = 0.d0)
!        Especially, we may have for example:
!        -> set itypfb(ifac) = iparoi
!        which prescribes default wall conditions for all variables at
!        face ifac,
!        -> and define IN ADDITION for variable ivar on this face
!        specific conditions by specifying
!        icodcl(ifac, ivar) and the 3 rcodcl values.


!      The user may also assign to itypfb a value not equal to
!       ientre, isolib, isymet, iparoi, iparug, iindef
!       but greater than or equal to 1 and less than or equal to
!       ntypmx (see values in param.h) to distinguish
!       groups or colors in other subroutines which are specific
!       to the case and in which itypfb is accessible.
!       In this case though it will be necessary to
!       prescribe boundary conditions by assigning values to
!       icodcl and to the 3 rcodcl fields (as the value of itypfb
!       will not be predefined in the code).


! Consistency rules
! =================

!       A few consistency rules between 'icodcl' codes for
!         variables with non-standard boundary conditions:

!           Codes for velocity components must be identical
!           Codes for Rij components must be identical
!           If code (velocity or Rij) = 4
!             we must have code (velocity and Rij) = 4
!           If code (velocity or turbulence) = 5
!             we must have code (velocity and turbulence) = 5
!           If code (velocity or turbulence) = 6
!             we must have code (velocity and turbulence) = 6
!           If scalar code (except pressure or fluctuations) = 5
!             we must have velocity code = 5
!           If scalar code (except pressure or fluctuations) = 6
!             we must have velocity code = 6


! Remarks
! =======

!       Caution: to prescribe a flux (nonzero) to Rij,
!                the viscosity to take into account is viscl
!                even if visct exists (visct=rho cmu k2/epsilon)

!       We have the ordering array for boundary faces from the
!           previous time step (except for the fist time step,
!           where 'itrifb' has not been set yet).
!       The array of boundary face types 'itypfb' has been
!           reset before entering the subroutine.


!       Note how to access some variables:

! Cell values
!               Let         iel = ifabor(ifac)

! * Density                         cell iel:
!                  propce(iel, ipproc(irom))
! * Dynamic molecular viscosity     cell iel:
!                  propce(iel, ipproc(iviscl))
! * Turbulent viscosity   dynamique cell iel:
!                  propce(iel, ipproc(ivisct))
! * Specific heat                   cell iel:
!                  propce(iel, ipproc(icp))
! * Diffusivity: lambda          scalaire iscal, cell iel:
!                  propce(iel, ipproc(ivisls(iscal)))

! Boundary face values

! * Density                        boundary face ifac :
!                  propfb(ifac, ipprob(irom))
! * Mass flow relative to variable ivar, boundary face ifac:
!      (i.e. the mass flow used for convecting ivar)
!                  propfb(ifac, pprob(ifluma(ivar )))
! * For other values                  at boundary face ifac:
!      take as an approximation the value in the adjacent cell iel
!      i.e. as above with iel = ifabor(ifac).

!-------------------------------------------------------------------------------
!   All cells can be identified by using the subroutine 'getcel'.
!    Syntax of getcel:
!     getcel(string, nelts, eltlst) :
!     - string is a user-supplied character string containing
!       selection criteria;
!     - nelts is set by the subroutine. It is an integer value
!       corresponding to the number of boundary faces verifying the
!       selection criteria;
!     - lstelt is set by the subroutine. It is an integer array of
!       size nelts containing the list of boundary faces verifying
!       the selection criteria.

!       string may contain:
!       - references to colors (ex.: 1, 8, 26, ...
!       - references to groups (ex.: inlet, group1, ...)
!       - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!
!       These criteria may be combined using logical operators
!       ('and', 'or') and parentheses.
!       Example: '1 and (group2 or group3) and y < 1' will select boundary
!       faces of color 1, belonging to groups 'group2' or 'group3' and
!       with face center coordinate y less than 1.
!
!   All boundary faces may be identified using the 'getfbr' subroutine.
!    Syntax of getfbr:
!     getfbr(string, nelts, eltlst) :
!     - string is a user-supplied character string containing
!       selection criteria;
!     - nelts is set by the subroutine. It is an integer value
!       corresponding to the number of boundary faces verifying the
!       selection criteria;
!     - lstelt is set by the subroutine. It is an integer array of
!       size nelts containing the list of boundary faces verifying
!       the selection criteria.
!
!     string may contain:
!     - references to colors (ex.: 1, 8, 26, ...
!     - references to groups (ex.: inlet, group1, ...)
!     - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!
!     These criteria may be combined using logical operators
!     ('and', 'or') and parentheses.
!
!   All internam faces may be identified using the 'getfac' subroutine.
!    Syntax of getfac:
!     getfac(string, nelts, eltlst) :
!     - string is a user-supplied character string containing
!       selection criteria;
!     - nelts is set by the subroutine. It is an integer value
!       corresponding to the number of boundary faces verifying the
!       selection criteria;
!     - lstelt is set by the subroutine. It is an integer array of
!       size nelts containing the list of boundary faces verifying
!       the selection criteria.
!
!     string may contain:
!     - references to colors (ex.: 1, 8, 26, ...
!     - references to groups (ex.: inlet, group1, ...)
!     - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!
!     These criteria may be combined using logical operators
!     ('and', 'or') and parentheses.
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! icodcl           ! ia ! --> ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1  -> Dirichlet                              !
!                  !    !     ! = 2  -> flux density                           !
!                  !    !     ! = 4  -> sliding wall and u.n=0 (velocity)      !
!                  !    !     ! = 5  -> friction and u.n=0 (velocity)          !
!                  !    !     ! = 6  -> roughness and u.n=0 (velocity)         !
!                  !    !     ! = 9  -> free inlet/outlet (velocity)           !
!                  !    !     !         inflowing possibly blocked             !
! itrifb(nfabor)   ! ia ! <-- ! indirection for boundary faces ordering)       !
! itypfb(nfabor)   ! ia ! --> ! boundary face types                            !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! rcodcl           ! ra ! --> ! boundary condition values                      !
!                  !    !     ! rcodcl(1) = Dirichlet value                    !
!                  !    !     ! rcodcl(2) = exterior exchange coefficient      !
!                  !    !     !  (infinite if no exchange)                     !
!                  !    !     ! rcodcl(3) = flux density value                 !
!                  !    !     !  (negative for gain) in w/m2 or                !
!                  !    !     !  roughness height (m) if icodcl=6              !
!                  !    !     ! for velocities           ( vistl+visct)*gradu  !
!                  !    !     ! for pressure                         dt*gradp  !
!                  !    !     ! for scalars    cp*(viscls+visct/sigmas)*gradt  !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use ihmpre
use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use ppcpfu
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)

! LOCAL VARIABLES


integer          ifac, ii
integer          izone
integer          iclafu
integer          ilelt, nlelt

double precision uref2, d2s3
double precision xkent, xeent

integer, allocatable, dimension(:) :: lstelt

!===============================================================================



!===============================================================================
! 1.  INITIALISATIONS

!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))

d2s3 = 2.d0/3.d0

!===============================================================================
! 2.  ALLOCATION OF BOUNDARY CONDITIONS VECTOR
!       LOOP ON BOUNDARY FACE
!         FAMILY AND PROPERTIES ARE DETERMINED
!         BOUNDARY CONDITION ALLOCATED
!
!     BOUNDARY CONDITIONS ON BONDARY FACE HAVE TO BE ALLOCATED HERE
!
!     USER'S WORK TO DO
!
!===============================================================================

!  Each kind of condition for extended physic is allocated with a number
!  IZONE ( 0<IZONE<= NOZPPM ; NOZPPM allocated in ppppar.h)

! ---- The 12 color is a pure air inlet

CALL GETFBR('12',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   kind of boundary condition for standard variables
  itypfb(ifac) = ientre

!   Zone number allocation
  izone = 1

!   Zone number storage
  izfppp(ifac) = izone

! ------ This inlet have a fixed mass flux

  ientat(izone) = 1
  iqimp(izone)  = 1
  inmoxy(izone) = 1
!     - Air mass flow rate kg/s
  qimpat(izone) = 1.46d-03
!     - Air inlet temperature K
  timpat(izone) = 400.d0 + tkelvi
!     - Fuel mass flow rate kg/s
  qimpfl(izone) = 0.d0

! ----- The 12 color is a fixed mass flow rate inlet
!        user gives only the speed vector direction
!        (spedd vector norm is irrelevent)
!
  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 5.d0

!   Boundary conditions of turbulence
  icalke(izone) = 1
!
!    - If ICALKE = 0 the boundary conditions of turbulence at
!      the inlet are calculated as follows:

  if(icalke(izone).eq.0) then

    uref2 = rcodcl(ifac,iu,1)**2                           &
           +rcodcl(ifac,iv,1)**2                           &
           +rcodcl(ifac,iw,1)**2
    uref2 = max(uref2,1.d-12)
    xkent  = epzero
    xeent  = epzero

    call keenin                                                   &
    !==========
      ( uref2, xintur(izone), dh(izone), cmu, xkappa,             &
        xkent, xeent )

    if    (itytur.eq.2) then

      rcodcl(ifac,ik,1)  = xkent
      rcodcl(ifac,iep,1) = xeent

    elseif(itytur.eq.3) then

      rcodcl(ifac,ir11,1) = d2s3*xkent
      rcodcl(ifac,ir22,1) = d2s3*xkent
      rcodcl(ifac,ir33,1) = d2s3*xkent
      rcodcl(ifac,ir12,1) = 0.d0
      rcodcl(ifac,ir13,1) = 0.d0
      rcodcl(ifac,ir23,1) = 0.d0
      rcodcl(ifac,iep,1)  = xeent

    elseif (iturb.eq.50) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iep,1)  = xeent
      rcodcl(ifac,iphi,1) = d2s3
      rcodcl(ifac,ifb,1)  = 0.d0

    elseif (iturb.eq.60) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iomg,1) = xeent/cmu/xkent

    endif

  endif
!
!    - If ICALKE = 1 the boundary conditions of turbulence at
!      the inlet refer to both, a hydraulic diameter and a
!      reference velocity given in usini1.f90.
!
  dh(izone)     = 0.032d0
!
!    - If ICALKE = 2 the boundary conditions of turbulence at
!      the inlet refer to a turbulence intensity.
!
  xintur(izone) = 0.d0

! ------ Treatment of user's scalars

  if ( (nscal-nscapp).gt.0 ) then
    do ii = 1, (nscal-nscapp)
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo

! ---- Inlet of both primary Air and Fuel

CALL GETFBR('11',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   kind of boundary condition for standard variables
  itypfb(ifac) = ientre

!   Zone number allocation
  izone = 2

!   Zone number storage
  izfppp(ifac) = izone

! ------ This inlet have a fixed mass flux

  ientfl(izone) = 1
  iqimp(izone)  = 1
  inmoxy(izone) = 1
!     - Air mass flow rate in kg/s
  qimpat(izone) = 1.46d-03
!     - Air Temperature at inlet in K
  timpat(izone) = 800.d0  + tkelvi

!     - Fuel mass flow rate in kg/s
  qimpfl(izone) = 1.46d-04/360.d0

!     - PERCENTAGE mass fraction of each granulometric class
!       ICLAFU (1 < ICLAFU < NCLAFU )
  iclafu = 1
  distfu(izone,iclafu) = 100.d0

!     - Fuel Temperature at inlet in K
  timpfl(izone) = 100.d0  + tkelvi

! ----- The 11 color is a fixed mass flow rate inlet
!        user gives only the speed vector direction
!        (spedd vector norm is irrelevent)
!

  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 5.d0

!   Boundary conditions of turbulence
  icalke(izone) = 1
!
!    - If ICALKE = 0 the boundary conditions of turbulence at
!      the inlet are calculated as follows:

  if(icalke(izone).eq.0) then

    uref2 = rcodcl(ifac,iu,1)**2                           &
           +rcodcl(ifac,iv,1)**2                           &
           +rcodcl(ifac,iw,1)**2
    uref2 = max(uref2,1.d-12)
    xkent  = epzero
    xeent  = epzero

    call keenin                                                   &
    !==========
      ( uref2, xintur(izone), dh(izone), cmu, xkappa,             &
        xkent, xeent )

    if    (itytur.eq.2) then

      rcodcl(ifac,ik,1)  = xkent
      rcodcl(ifac,iep,1) = xeent

    elseif(itytur.eq.3) then

      rcodcl(ifac,ir11,1) = d2s3*xkent
      rcodcl(ifac,ir22,1) = d2s3*xkent
      rcodcl(ifac,ir33,1) = d2s3*xkent
      rcodcl(ifac,ir12,1) = 0.d0
      rcodcl(ifac,ir13,1) = 0.d0
      rcodcl(ifac,ir23,1) = 0.d0
      rcodcl(ifac,iep,1)  = xeent

    elseif (iturb.eq.50) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iep,1)  = xeent
      rcodcl(ifac,iphi,1) = d2s3
      rcodcl(ifac,ifb,1)  = 0.d0

    elseif (iturb.eq.60) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iomg,1) = xeent/cmu/xkent

    endif

  endif
!
!    - If ICALKE = 1 the boundary conditions of turbulence at
!      the inlet refer to both, a hydraulic diameter and a
!      reference velocity given in usini1.f90.
!
  dh(izone)     = 0.032d0
!
!    - If ICALKE = 2 the boundary conditions of turbulence at
!      the inlet refer to a turbulence intensity.
!
  xintur(izone) = 0.d0
enddo

! --- Color 15 is a wall

CALL GETFBR('15',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          WALL  : nul mass flow rate (nul pressure flux)
!                  rubbing for speed (and turbulence)
!                  nul scalar fluxes

!   kind of boundary condition for standard variables
  itypfb(ifac)   = iparoi


!   Zone number allocation
  izone = 3

!   Zone number storage
  izfppp(ifac) = izone

enddo

! --- Color 19 is an outlet

CALL GETFBR('19',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          OUTLET : nul fluxes for speed and scalar
!                   pressure fixed

!   kind of boundary condition for standard variables
    itypfb(ifac)   = isolib

!   Zone number allocation
    izone = 4

!   Zone number storage
    izfppp(ifac) = izone

  enddo

! --- 14 and 4 are symetry

CALL GETFBR('14 or 4',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SYMETRY

!   kind of boundary condition for standard variables
  itypfb(ifac)   = isymet

!   Zone number allocation
  izone = 5

!   Zone number storage
  izfppp(ifac) = izone

enddo

! Deallocate the temporary array
deallocate(lstelt)

!----
! END
!----

return
end subroutine
