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

subroutine uslag2 &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa  , ifrlag ,                            &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   )

!===============================================================================
! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module:
!   -------------------------------------

!    User subroutine (Mandatory intervention)

!    User subroutine for the boundary conditions associated to the particles
!    (inlet and treatment of the other boundaries)


! Boundary faces identification
! =============================

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! i  ! <-- ! maximum number of particles allowed            !
! nvp              ! i  ! <-- ! number of particle variables                   !
! nvp1             ! i  ! <-- ! nvp minus position, fluid and part. velocities !
! nvep             ! i  ! <-- ! number of particle properties (integer)        !
! nivep            ! i  ! <-- ! number of particle properties (integer)        !
! ntersl           ! i  ! <-- ! number of source terms of return coupling      !
! nvlsta           ! i  ! <-- ! nb of Lagrangian statistical variables         !
! nvisbr           ! i  ! <-- ! number of boundary statistics                  !
! itrifb(nfabor)   ! ia ! <-- ! indirection for the sorting of the             !
! itypfb(nfabor)   ! ia ! <-- ! type of the boundary faces                     !
! ifrlag(nfabor    ! ia ! --> ! type of the Lagrangian boundary faces          !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! transported variables at the previous          !
! (ncelet,*)       !    !     ! time step                                      !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! coefa, coefb     ! ra ! <-- ! boundary conditions at the boundary faces      !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
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
use cstnum
use cstphy
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use cpincl
use ihmpre
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itypfb(nfabor) , itrifb(nfabor)
integer          itepa(nbpmax,nivep) , ifrlag(nfabor)

double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)

! Local variables

integer          ifac , izone, nbclas, iclas
integer          icha
integer          ilelt, nlelt

double precision pis6 , mp0 , temp

integer, allocatable, dimension(:) :: lstelt

!===============================================================================



!===============================================================================
! 1.  Memory management
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))


!===============================================================================
! 2. Initialization
!===============================================================================

pis6 = pi / 6.d0

!===============================================================================
! 3. Construction of the boundary zones
!===============================================================================


!     Definition of the boundary zones
!     --------------------------------

!     For the Lagrangian module, the user defines nfrlag boundary zones
!     from the color of the boundary faces, or more generally from their properties
!     (colors, groups..) or from the boundary conditions prescribed in usclim, or
!     even from their coordinates. To do that, we fill the ifrlag(nfabor) array
!     which gives for every boundary face the number of the zone to which it belongs ifrlag(ifac)
!
!     Be careful, all the boundary faces must have been affected.
!
!     The number of the zones (thus the values of ifrlag(ifac)) is arbitrarily
!     chosen by the user, but must be a positive integer and inferior or equal to
!     nflagm (parameter prescribed in lagpar.h).
!
!     Afterwards, we assign to every zone a type named itylag that will be used
!     to impose global boundary conditions.

izone = -1

! ---> First zone numbered izone=1 ( = color 10)
CALL GETFBR('10',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone        = 1
  ifrlag(ifac) = izone

enddo

! ---> Second zone numbered izone=2 ( = part of color 4)
CALL GETFBR('4 and Y < 1.0',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone        = 2
  ifrlag(ifac) = izone

enddo

! ---> Third zone numbered izone=3 ( = inlet phase 1)
do ifac = 1, nfabor
  if(itypfb(ifac).eq.ientre) then
    izone        = 4
    ifrlag(ifac) = izone
  endif
enddo

! ---> Nth zone numbered izone=5 (= color 3)
CALL GETFBR('3',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone        = 5
  ifrlag(ifac) = izone

enddo


!===============================================================================
! 4. Injection per particle class into the calculation domain
!===============================================================================

!   TO PROVIDE INFORMATION ABOUT THE PARTICLE CLASSES,
!   WE FOLLOW A TWO-STEP PROCEDURE:

!   1) FIRST, THE NUMBER OF PARTICLE CLASSES IS PRESCRIBED
!      FOR EACH BOUNDARY ZONE: IUSNCL (by default, this parameter is equal to zero)

!   2) AFTERWARDS, FOR EACH ZONE AND FOR EACH CLASS, WE PRESCRIBE
!      THE PHYSICAL PROPERTIES OF THE PARTICLES
!




! --> Number of particle classes entering the domain
!   We assign here the number of classes for each zone previously identified.
!
!   This number is zero by default.
!   The maximal number of classes is nclagm (defined in lagpar.h)

! ---> First zone numbered izone = 1: 1 class injected
izone     = 1
nbclas    = 1
iusncl(izone) = nbclas

! ---> Second zone numbered izone = 2: 0 class injected
izone     = 2
nbclas    = 0
iusncl(izone) = nbclas

! ---> Third zone numbered izone = 4 : 0 class injected
izone     = 4
nbclas    = 0
iusncl(izone) = nbclas

! ---> Zone numbered izone = 5 : 0 class injected
izone     = 5
nbclas    = 0
iusncl(izone) = nbclas


! --> For every class associated with a zone,
!     we give the followong information.


!     iusncl number of classes per zone
!     iusclb boundary conditions for the particles
!     = ientrl -> zone of particle inlet
!     = isortl -> particle outlet
!     = irebol -> rebound of the particles
!     = idepo1 -> definitive deposition
!     = idepo2 -> definitive deposition, but the particle remains in memory
!                 (useful only if iensi2 = 1)
!     = idepo3 -> deposition and resuspension possible
!                 following the conditions of the flow
!     = idepfa -> deposition of the particle with attachment force,
!                 the velocity is conserved, and resuspension is possible
!                 (Possible if ladlvo = 1 )
!     = iencrl -> fouling (coal only iphyla = 2)
!     = jbord1 -> user-defined particle/boundary interaction (cf. uslabo)
!     = jbord2 -> user-defined particle/boundary interaction (cf. uslabo)
!     = jbord3 -> user-defined particle/boundary interaction (cf. uslabo)
!     = jbord4 -> user-defined particle/boundary interaction (cf. uslabo)
!     = jbord5 -> user-defined particle/boundary interaction (cf. uslabo)



!     Array iuslag :
!     ================
!        ijnbp : number of particles per class and per zone
!        ijfre : injection frequency. If ijfre = 0, then the injection
!                occurs only at the first absolute iteration.
!        iclst : number of the group to which the particle belongs
!                (only if one wishes to calculate statistics per group)
!        ijuvw : type of condition on the velocity
!                  = -1 imposed flow velocity
!                  =  0 imposed velocity along the normal direction of the
!                      boundary face, with norm equal to RUSLAG(ICLAS,IZONE,IUNO)
!                  =  1 imposed velocity: we prescribe   RUSLAG(ICLAS,IZONE,IUPT)
!                                                        RUSLAG(ICLAS,IZONE,IVPT)
!                                                        RUSLAG(ICLAS,IZONE,IWPT)
!                  =  2 user-defined profile
!        ijprtp : type of temperature condition
!                  =  1 imposed temperature: we prescribe RUSLAG(ICLAS,IZONE,ITPT)
!                  =  2 user-defined profile
!        ijprdp : type of diameter condition
!                  =  1 imposed diameter: we prescribe  RUSLAG(ICLAS,IZONE,IDPT)
!                                                       RUSLAG(ICLAS,IZONE,IVDPT)
!                  =  2 user-defined profile
!        inuchl : number of the coal of the particle (only if iphyla = 2)

!     Array ruslag :
!     ===============
!        iuno  : Norm of the velocity (m/s)
!        iupt  : Velocity along the X axis, for each class and for each zone (m/s)
!        ivpt  : Velocity along the Y axis, for each class and for each zone (m/s)
!        iwpt  : Velocity along the Z axis, for each class and for each zone (m/s)
!        idebt : Mass flow rate (kg/s)
!        ipoit : Statistical weight (number of samples) associated
!                to the particle (automatically computed to respect a mass
!                flow rate if it is defined)

!        Physical characteristics:
!          idpt   : diameter (m)
!          ivdpt  : standard deviation of the diameter (m)
!          itpt   : temperature in Celsius degress (no enthalpy)
!          icpt   : specific heat (J/kg/K)
!          iepsi  : emissivity (if =0 then no radiative effect is taken into account)
!          iropt  : density (kg/m3)

!         If coal (iphyla=2)
!            ihpt  : temperature in Celsius degress (no enthalpy)
!            imcht : mass of reactive coal (kg)
!            imckt : masse of coke (kg)


! ---> EXAMPLE : First zone, numbered IZONE = 1 (NBCLAS classes)
!        IUSCLB : adherence of the particle to a boundary face
!        IJNBP  : 10 particles for each class,
!        IJFRE  : injection every other time step
!        IJUVW, IUPT, IVPT, IWPT : imposed velocity on 1.1D0, 0.0D0, 0.0D0
!        ICPT   : cp equal to 10000
!        ITPT   : temperature equal to 25 Celsius degress
!        IDPT   : diameter equal to 50.E-6 m
!        IEPSI  : emissivity equal to 0.7
!        IVDPT  : constant diameter ==> standard deviation null
!        IROPT  : density
!        IPOIT  : statistical weight (number of physical particles
!                 represented by one statistical particle)
!        IDEBT  : mass flow rate


izone     = 1
nbclas    = iusncl(izone)
iusclb (izone)         =  ientrl
do iclas  = 1, nbclas

  iuslag (iclas,izone,ijnbp) = 10
  iuslag (iclas,izone,ijfre) = 2

  if (nbclst.gt.0) then
    iuslag(iclas,izone,iclst) = 1
  endif

  iuslag (iclas,izone,ijuvw) = -1
  ruslag (iclas,izone,iupt)  = 1.1d0
  ruslag (iclas,izone,ivpt)  = 0.0d0
  ruslag (iclas,izone,iwpt)  = 0.0d0
  iuslag (iclas,izone,ijprpd)= 1
  ruslag (iclas,izone,ipoit) = 1.d0
  ruslag (iclas,izone,idebt) = 0.d0

!    if the physics is " simple"

  if ( iphyla.eq.0 .or. iphyla.eq.1 ) then

!        Mean value and standard deviation of the diameter

    iuslag (iclas,izone,ijprdp)= 1
    ruslag (iclas,izone,idpt)  = 50.d-6
    ruslag (iclas,izone,ivdpt) = 0.d0

!        Density

    ruslag(iclas,izone,iropt) = 2500.d0

    if ( iphyla.eq.1 ) then

!        Temperature and Cp

      if ( itpvar.eq.1 ) then
        iuslag (iclas,izone,ijprtp) = 1
        ruslag(iclas,izone,itpt)    = 20.d0

        ruslag(iclas,izone,icpt)    = 1400.d0
        ruslag(iclas,izone,iepsi)   = 0.7d0
      endif

    endif

!    Coal

  else if ( iphyla.eq.2 ) then

!    CAUTION :   1) To transport and burn coal particles with the Lagrangian
!                   module, a specific physics for the dispersed phase must
!                   be activated for the carrier phase.
!
!                2) The physical properties of the coal particles are known
!                   from the thermo-chemical file: dp_FCP
!
!                3) For the current phase ICLAS, and for the current boundary zone
!                   NB, we assign to the coal particles the properties of the coal ICHA
!                   of the icha class taken from the file dp_FCP.
!
!                4) icha : number of the coal between 1 and ncharb defined by the user
!                   in the file dp_FCP.
!


    icha = ichcor(iclas)
    temp = 800.d0

!        Number of the coal

    iuslag(iclas,izone,inuchl) = icha

!        Temperature and Cp

    ruslag(iclas,izone,ihpt) = temp
    ruslag(iclas,izone,icpt) = cp2ch(icha)

!        Mean value and standard deviation of the diameter

    ruslag (iclas,izone,idpt)  = diam20(iclas)
    ruslag (iclas,izone,ivdpt) = 0.d0

!        Density

    ruslag(iclas,izone,iropt) =  rho0ch(icha)

!        Mass of reactive coal and
!        mass of coke (null if the coal has never burnt)

    mp0 = pis6 * ( ruslag(iclas,izone,idpt)**3 )                  &
               * ruslag(iclas,izone,iropt)
    ruslag(iclas,izone,imcht) = mp0 * (1.d0-xashch(icha))
    ruslag(iclas,izone,imckt) = 0.d0

  endif

enddo

! ---> Second zone, numbered izone = 2
!        IUSCLB : rebound of the particle

izone     = 2
iusclb (izone)         =  irebol


! same procedure for the other zones...

!===============================================================================

!--------
! Formats
!--------

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
