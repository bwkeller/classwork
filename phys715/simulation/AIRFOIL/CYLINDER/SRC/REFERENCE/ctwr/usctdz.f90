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

subroutine usctdz
!================

!===============================================================================
! FONCTION :
! --------

! ROUTINE UTILISATEUR POUR LA DEFINITION DES ZONES D'ECHANGE
! D'UN AEROREFRIGERANT

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use entsor
use cstphy
use parall
use period
use ppppar
use ppthch
use ppincl
use ctincl
use mesh

!===============================================================================

implicit none

! Local variables

integer          imzech,ntypze,idimze,neleze

double precision teaueze,qeaueze,deltat
double precision xap,xnp,surface,dgout


!===============================================================================




! La modelisation choisie doit etre coherente avec IPPMOD

! IMZECH = 0 - pas de modele
!          1 - modele de Merkel
!          2 - modele de Poppe

imzech  = ippmod(iaeros)

! Definition de la zone d'echange

idimze  = 2
ntypze  = 2
neleze  = 20
deltat  = 10.0d0
teaueze = 36d0
qeaueze = 33737.d0
xap     = 0.2d0
xnp     = 0.5d0
surface = 5183.6d0
dgout   = 0.005d0

call defct &
!=========
 ( idimze, '2 or 3', imzech, ntypze, neleze,          &
   deltat, teaueze, qeaueze, xap, xnp, surface, dgout )


!===============================================================================
!     FORMATS
!===============================================================================

return
end subroutine

