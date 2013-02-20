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

subroutine usctcl &
!================

 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl )

!===============================================================================
! FONCTION :
! --------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE
!            MODELISATION AEROREFRIGERANTS
!    REMPLISSAGE DU TABLEAU DE CONDITIONS AUX LIMITES
!    (ICODCL,RCODCL) POUR LES VARIABLES INCONNUES
!    PENDANT DE USCLIM.F



!    CE SOUS PROGRAMME UTILISATEUR EST OBLIGATOIRE
!    =============================================


! Introduction
! ============

! Here one defines boundary conditions on a per-face basis.

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


! Boundary condition types
! ========================

! Boundary conditions setup for standard variables (pressure, velocity,
! turbulence, scalars) is described precisely in the 'usclim' subroutine.

! Detailed explanation will be found in the theory guide.


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
! itrifb           ! ia ! <-- ! indirection for boundary faces ordering        !
! itypfb           ! ia ! --> ! boundary face types                            !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! rcodcl           ! ra ! --> ! boundary condition values                      !
!  (nfabor,nvar,3) !    !     ! rcodcl(1) = Dirichlet value                    !
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
use ppincl
use ctincl
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

! Local variables

integer          ifac, iel, ii, ivar
integer          ilelt, nlelt, izone

integer, allocatable, dimension(:) :: lstelt

!===============================================================================


!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))


!===============================================================================
! 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR LES FACES DE BORD
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LA CONDITION LIMITE

!          IMPOSER ICI LES CONDITIONS LIMITES SUR LES FACES DE BORD

!===============================================================================

!   Pour chaque type de condition relative aux physiques particulieres
!       on affecte un numero de zone de maniere a pouvoir donner les
!       conditions aux limites par zone physique et non par face de maillage
!       Un numero de zone est un entier arbitraire strictement positif
!         et inferieur ou egal a NOZPPM (dont la valeur est fixee en
!         parametre dans ppppar.h)


! --- On impose en couleur 2 une sortie libre pour toutes les phases
CALL GETFBR('2',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone = 1
  izfppp(ifac) = izone

!          SORTIE : FLUX NUL VITESSE ET TEMPERATURE, PRESSION IMPOSEE
!            Noter que la pression sera recalee a P0
!                sur la premiere face de sortie libre (ISOLIB)


  itypfb(ifac) = isolib

!       On impose un profil de pression sur toutes les faces
  icodcl(ifac,ipr) = 1
  rcodcl(ifac,ipr,1) =                                     &
    ro0*(  gx*(cdgfbo(1,ifac)-xyzp0(1))              &
                + gy*(cdgfbo(2,ifac)-xyzp0(2))              &
                + gz*(cdgfbo(3,ifac)-xyzp0(3)))

enddo

! --- On impose en couleur 4 une paroi pour toutes les phases
CALL GETFBR('4',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone = 2
  izfppp(ifac) = izone

!          PAROI : DEBIT NUL (FLUX NUL POUR LA PRESSION)
!                  FROTTEMENT POUR LES VITESSES (+GRANDEURS TURB)
!                  FLUX NUL SUR LES SCALAIRES

  itypfb(ifac)   = iparoi

enddo

! --- On impose une symetrie en couleur 5 pour toutes les phases
CALL GETFBR('5',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone = 3
  izfppp(ifac) = izone

!          SYMETRIES

  itypfb(ifac) = isymet

enddo

!----
! FORMATS
!----

!----
! FIN
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
