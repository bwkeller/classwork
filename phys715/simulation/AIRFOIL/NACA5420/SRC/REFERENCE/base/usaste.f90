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

subroutine usaste &
!================

 ( idfstr )

!===============================================================================
! Purpose:
! -------

!    User subroutine dedicated the Fluid-Structure external coupling
!    with Code_Aster :
!          Here one defines the boundary faces coupled
!          with Code_Aster and the fluid forces components
!          which are given to structural calculations

!    Boundary faces identification
!    =============================

!    Boundary faces may be identified using the 'getfbr' subroutine.
!    The syntax of this subroutine is described in the 'usclim' subroutine,
!    but a more thorough description can be found in the user guide.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idfstr(nfabor)   ! ia ! <-- ! boundary faces -> structure definition         !
! ia(*)            ! ia ! --- ! main integer work array                        !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use optcal
use entsor
use albase
use parall
use period
use alaste
use mesh

!===============================================================================

implicit none

! Arguments

integer          nbstru

integer          idfstr(nfabor)

! Local variables

integer          ifac
integer          ilelt, nlelt

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

!===============================================================================
! 1.  INITIALIZATION
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))


!===============================================================================
! 2.  Definition of external structures
!===============================================================================

!    Here one fills array IDFSTR(NFABOR)
!    For each boundary face IFAC, IDFSTR(IFAC) is the number of the structure
!    the face belongs to (if IDFSTR(IFAC) = 0, the face IFAC doesn't
!    belong to any structure.)
!    When using external coupling with Code_Aster, structure number necessarily
!    needs to be negative (as shown in following examples).

!    The number of "external" structures with Code Aster is automatically
!    defined with the maximum absolute value of IDFSTR table, meaning that
!    external structure numbers must be defined sequentially with negative values
!    beginning with integer value '-1'.


!    In following example, boundary faces with color 2 and which abscissa X < 2.0
!    belong to external structure '-1'.
!    Boundary faces with color 2 and which abscissa X > 2.0 belong to external
!    structure '-2'. The total number of external structures coupled with Code_Aster
!    equals 2.

!    Boundary faces identification
!    =============================

!    Boundary faces may be identified using the 'getfbr' subroutine.
!    The syntax of this subroutine is described in the 'usclim' subroutine,
!    but a more thorough description can be found in the user guide.

!================================================================================

CALL GETFBR('2 and X < 2.0',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = -1

enddo


CALL GETFBR('2 and X > 2.0',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = -2

enddo

! --- The movement of external structure called '-1' is blocked in Z direction.

asddlf(3,1) = 0

! --- The movement of external structure called '-2' is blocked in Z direction.

asddlf(3,2) = 0

!----
! Formats
!----

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
