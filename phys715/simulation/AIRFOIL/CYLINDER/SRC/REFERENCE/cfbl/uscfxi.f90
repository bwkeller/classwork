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

subroutine uscfxi &
!================

 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Initialize the unknown variables for the compressible flow scheme.


! Description
! ===========

! This subroutine is similar to the user subroutine 'usiniv', but
! is dedicated to the compressible flow scheme.
! It is called at the beginning of the computation (only if it is
! not a restart), just before the time marching loop starts.
! It allows to initialize all the unknown variables.

! The standard initialization has been reproduced here as an example.

! More examples can be found in 'usiniv'.


! Physical properties
! ===================

! The physical properties (viscosity, specific heat, thermal
! conductivity, Schmidt number) that are stored in the arrays propce,
! propfa and propfb must not be modified here: if it is necessary to
! do so, it must be done in the dedicated user programme 'uscfpv'.


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp(ncelet,*)    ! ra ! <-> ! calculated variables at cell centers           !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use mesh
use ppincl

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

! Local variables

integer          iel

integer          iccfth, iscal, imodif, iutile

integer, allocatable, dimension(:) :: lstelt

double precision, allocatable, dimension(:) :: w1, w2, w3, w4

!===============================================================================


!===============================================================================
! 1. Control print
!===============================================================================

write(nfecra,9001)

!===============================================================================
! 2.  Initialization of local variables
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet), w4(ncelet))

imodif = 1

!===============================================================================
! 3. Unknown variable initialization
!      for initial calculations (not in case of restart)
!===============================================================================

if ( isuite.eq.0 ) then

! --- Velocity components

  do iel = 1, ncel
    rtp(iel,iu) = 0.d0
    rtp(iel,iv) = 0.d0
    rtp(iel,iw) = 0.d0
  enddo


! --- User defined scalars

  ! If there are user defined scalars
  if(nscaus.gt.0) then
    ! For each scalar
    do iscal = 1, nscaus
      ! Initialize each cell value
      do iel = 1, ncel
        rtp(iel,isca(iscal)) = 0.d0
      enddo
    enddo
  endif


! --- Pressure, Density, Temperature, Total Energy

  ! Only 2 out of these 4 variables are independent: one may choose to
  ! initialize any pair of variables picked out of these 4, except
  ! (Temperature-Energy). The remaining 2 variables will be deduced
  ! automatically.


  ! Initialize 2 and only 2 variables

  !   To do so, set iutile=1 for each of the 2 selected variables
  !             and iutile=0 for each of the 2 others

  !   In the example provided below, Pressure and Temperature are
  !   initialized.


  ! iccfth indicates which variables have been set:
  !   it is completed automatically for each variable and
!     it must not be modified.
  iccfth = 10000


  ! 1. Pressure (Pa)
  iutile = 1
  if(iutile.eq.1) then
    iccfth = iccfth*2
    do iel = 1, ncel
      rtp(iel,ipr         ) = p0
    enddo
  endif

  ! 2. Density (kg/m3)
  iutile = 0
  if(iutile.eq.1) then
    iccfth = iccfth*3
    do iel = 1, ncel
      rtp(iel,isca(irho  )) = ro0
    enddo
  endif

  ! 3. Temperature (K -- Warning: Kelvin)
  iutile = 1
  if(iutile.eq.1) then
    iccfth = iccfth*5
    do iel = 1, ncel
      rtp(iel,isca(itempk)) = t0
    enddo
  endif

  ! 4. Total Energy (J/kg)
  iutile = 0
  if(iutile.eq.1) then
    iccfth = iccfth*7
    do iel = 1, ncel
      rtp(iel,isca(ienerg)) = cv0*t0
    enddo
  endif


  ! ** The following subroutine returns automatically the values for the
  ! two remaining variables that need to be computed, using the
  ! indicator iccfth.

  call uscfth                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     )


endif

! Free memory
deallocate(w1, w2, w3, w4)

!----
! Formats
!----

 9001 format(                                                     &
/,                                                                &
'  uscfxi: User defined initialization of the variables.',/,      &
/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
