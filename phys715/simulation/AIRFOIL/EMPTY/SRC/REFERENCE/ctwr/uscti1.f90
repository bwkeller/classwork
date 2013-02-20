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

subroutine uscti1
!================


!===============================================================================
!  FONCTION  :
!  ---------


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
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use ctincl

!===============================================================================

implicit none

double precision cpa,cpe,cpv,hv0,rhoe,visc,conduc

!===============================================================================

!===============================================================================


!===============================================================================
! 1.  PARAMETRES POUR L'ECART DE TEMPERATURE IMPOSE
!===============================================================================

!     ACTIVATION
iaeeri = 0

!     ECART DE REFRIGERATION A IMPOSER
vaeeri = 13.d0

!     FREQUENCE DE MODIFICATION DE LA TEMPERATURE
iaeerp = 5

!     PAS DE TEMPERATURE POUR LE CALCUL DE LA PENTE DE ECARTREF(TEAU)
paseri = 0.015d0

!     MAXIMUM DE LA TEMPERATURE D'EAU CHAUDE MOYENNE PONDEREE
aetemx = 80.d0

!     MINIMUM DE LA TEMPERATURE D'EAU REFROIDIE MOYENNE PONDEREE
aetemn = 10.d0

!     NOMBRE DE ZONES D'ECHANGES AYANT UNE FRONTIERE ENTREE EAU

nbzsup = 2

!     LISTE DES NBZSUP ZONES D'ECHANGES EN BORD DE L'ENTREE EAU

lizsup(1) = 1
lizsup(2) = 2

!     NOMBRE DE ZONES D'ECHANGES AYANT UNE FRONTIERE SORTIE EAU
nbzinf = 2

!     LISTE DES NBZINF ZONES D'ECHANGES EN BORD DE LA SORTIE EAU

lizinf(1) = 1
lizinf(2) = 2

!     INSTANT ACTIVATION ECART IMPOSE

inbaei = 1000.D0

!===============================================================================
! 2.  POST-PROCESSING DES ZONES D'ECHANGES
!===============================================================================

ichrze = 1

!===============================================================================
! 3.  CALCUL SUITE AEROREFRIGERANT
!===============================================================================

isuict = isuite

!===============================================================================
! 4.  PROPRIETES DE L'AIR
!===============================================================================

! Il est deconseille de modifier ici ces proprietes

cpa    = 1006.0d0
cpv    = 1831.0d0
cpe    = 4179.0d0
hv0    = 2501600.0d0
rhoe   = 997.85615d0
visc   = 1.765d-5
conduc = 0.02493d0

call ctprof &
!==========
( cpa, cpv, cpe, hv0, rhoe, visc, conduc, gx, gy, gz )

!----
! FIN
!----

return
end subroutine
