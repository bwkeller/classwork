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

subroutine usvima &
!================

 ( nvar   , nscal  ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   viscmx , viscmy , viscmz )

!===============================================================================
! Purpose:
! -------

!    User subroutine dedicated the use of ALE (Arbitrary Lagrangian Eulerian Method) :
!                 fills mesh viscosity arrays.


!    This subroutine is called at the beginning of each time step.

!    Here one can modify mesh viscosity value to prevent cells and nodes
!    from huge displacements in awkward areas, such as boundary layer for example.

!    IF variable IORTVM = 0, mesh viscosity modeling is isotropic therefore VISCMX
!    array only needs to be filled.
!    IF variable IORTVM = 1, mesh viscosity modeling is orthotropic therefore
!    all arrays VISCMX, VISCMY and VISCMZ need to be filled.

!    Note that VISCMX, VISCMY and VISCMZ arrays are initialized at the first time step to the value of 1.

! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
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
! viscmx(ncelet)    ! ra ! <-- ! mesh viscosity in X direction                 !
! viscmy(ncelet)    ! ra ! <-- ! mesh viscosity in Y direction                 !
! viscmz(ncelet)    ! ra ! <-- ! mesh viscosity in Z direction                 !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use albase
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision viscmx(ncelet), viscmy(ncelet), viscmz(ncelet)

! Local variables

integer          iel
double precision rad, xr2, xcen, ycen, zcen

!===============================================================================



!===============================================================================
!  1. Example :
!       One gives a huge mesh viscosity value to the cells included in a circle
!       with a radius 'rad' and a center '(xcen, ycen, zcen)'

!     In general it appears quite much easier to fill mesh viscosity arrays at
!     the beginning of the calculations basing on the initial geometry.

if (ntcabs.eq.0) then
  rad = (1.d-3)**2
  xcen  = 1.d0
  ycen  = 0.d0
  zcen  = 0.d0

  do iel = 1, ncel
    xr2 = (xyzcen(1,iel)-xcen)**2 + (xyzcen(2,iel)-ycen)**2       &
         + (xyzcen(3,iel)-zcen)**2
    if (xr2.lt.rad) viscmx(iel) = 1.d10
  enddo

! 2. In case of orthotropic mesh viscosity modeling one can choose
!    to submit nodes to a lower stress in Z direction
  if (iortvm.eq.1) then
    do iel = 1, ncel
      viscmy(iel) = viscmx(iel)
      viscmz(iel) = 1.d0
    enddo
  endif

endif

!----
! Formats
!----

!----
! End
!----

return
end subroutine
