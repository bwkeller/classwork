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

subroutine uskpdc &
!================

 ( nvar   , nscal  ,                                              &
   ncepdp , iappel ,                                              &
   icepdc , izcpdc ,                                              &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc )

!===============================================================================
! FONCTION :
! ----------

!                    PERTES DE CHARGE (PDC)

! IAPPEL = 1 :
!             CALCUL DU NOMBRE DE CELLULES OU L'ON IMPOSE UNE PDC
! IAPPEL = 2 :
!             REPERAGE DES CELLULES OU L'ON IMPOSE UNE PDC
! IAPPEL = 3 :
!             CALCUL DES VALEURS DES COEFS DE PDC


! CKUPDC EST LE COEFF DE PDC CALCULE.

!  IL INTERVIENT DANS LA QDM COMME SUIT :
!    RHO DU/DT = - GRAD P + TSPDC        (+ AUTRES TERMES)
!                      AVEC TSPDC = - RHO CKUPDC U ( en kg/(m2 s))


!  POUR UNE PDC REPARTIE,

!    SOIT KSIL = DHL/(0.5 RHO U**2) DONNE DANS LA LITTERATURE
!    (DHL EST LA PERTE DE CHARGE PAR UNITE DE LONGUEUR)

!    LE TERME SOURCE TSPDC VAUT DHL = - KSIL *(0.5 RHO U**2)

!    ON A CKUPDC = 0.5 KSIL ABS(U)


!  POUR UNE PDC SINGULIERE,

!    SOIT KSIS = DHS/(0.5 RHO U**2) DONNE DANS LA LITTERATURE
!    (DHS EST LA PERTE DE CHARGE SINGULIERE)

!    LE TERME SOURCE TSPDC VAUT DHS/L = - KSIS/L *(0.5 RHO U**2)

!    ON A CKUPDC = 0.5 KSIS/L ABS(U)

!    OU L DESIGNE LA LONGUEUR SUR LAQUELLE
!           ON A CHOISI DE REPRESENTER LA ZONE DE PDC SINGULIERE


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! iappel           ! e  ! <-- ! indique les donnes a renvoyer                  !
! icepdc(ncepdp    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! izcpdc(ncelet)   ! ia ! <-- ! cells zone for head loss definition            !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
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
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp
integer          iappel

integer          icepdc(ncepdp)
integer          izcpdc(ncel)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6)

! Local variables

integer          iel, ielpdc, ikpdc
integer          ilelt, nlelt
integer          izone
integer          iutile

double precision alpha, cosalp, sinalp, vit, ck1, ck2

integer, allocatable, dimension(:) :: lstelt

!===============================================================================



!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))


if(iappel.eq.1.or.iappel.eq.2) then

!===============================================================================

! 1. POUR CHAQUE PHASE : UN OU DEUX APPELS

!      PREMIER APPEL :

!        IAPPEL = 1 : NCEPDP : CALCUL DU NOMBRE DE CELLULES
!                                AVEC PERTES DE CHARGE


!      DEUXIEME APPEL (POUR LES PHASES AVEC NCEPDP > 0) :

!        IAPPEL = 2 : ICEPDC : REPERAGE DU NUMERO DES CELLULES
!                                AVEC PERTES DE CHARGE

! REMARQUES :

!        Ne pas utiliser CKUPDC dans cette section
!          (il est rempli au troisieme appel, IAPPEL = 3)

!        Ne pas utiliser ICEPDC dans cette section
!           au premier appel (IAPPEL = 1)

!        On passe ici a chaque pas de temps
!           (ATTENTION au cout calcul de vos developpements)

!===============================================================================


!  1.1 A completer par l'utilisateur : selection des cellules
!  -----------------------------------------------------------

! --- Exemple 1 : Aucune pdc (defaut)

  ielpdc = 0


! --- Exemple 2 : Pdc definies par coordonnees pour la phase 1
!                 Pas de pertes de charge pour la phase 2
!                 Le traitement etant different pour les phases
!                   un test est necessaire.

!       Ce test permet de desactiver l'exemple
  if(1.eq.0) then

    izone = 0
    ielpdc = 0

    call getcel('X <= 6.0 and X >= 4.0 and Y >= 2.0 and'//      &
         'Y <= 8.0',nlelt,lstelt)

    izone = izone + 1

    do ilelt = 1, nlelt
      iel = lstelt(ilelt)
      izcpdc(iel) = izone
      ielpdc = ielpdc + 1
      if (iappel.eq.2) icepdc(ielpdc) = iel
    enddo

  endif


!  1.2 Sous section generique a ne pas modifier
!  ---------------------------------------------

! --- Pour IAPPEL = 1,
!      Renseigner NCEPDP, nombre de cellules avec pdc
!      Le bloc ci dessous est valable pourles 2 exemples ci dessus

  if (iappel.eq.1) then
    ncepdp = ielpdc
  endif

!-------------------------------------------------------------------------------

elseif(iappel.eq.3) then

!===============================================================================

! 2. POUR CHAQUE PHASE AVEC NCEPDP > 0 , TROISIEME APPEL

!      TROISIEME APPEL (POUR LES PHASES AVEC NCEPDP > 0) :

!       IAPPEL = 3 : CKUPDC : CALCUL DES COEFFICIENTS DE PERTE DE CHARGE
!                             DANS LE REPERE DE CALCUL
!                             STOCKES DANS L'ORDRE
!                             K11, K22, K33, K12, K13, K23


!    REMARQUE :

!        Veillez a ce que les coefs diagonaux soient positifs.

!        Vous risquez un PLANTAGE si ce n'est pas le cas.

!        AUCUN controle ulterieur ne sera effectue.

!      ===========================================================


!  2.1 A completer par l'utilisateur : valeur des coefs
!  -----------------------------------------------------

! --- Attention
!   Il est important que les CKUPDC soient completes (par des valeurs
!     nulles eventuellement) dans la mesure ou ils seront utilises pour
!     calculer un terme source dans les cellules identifiees precedemment.

!   On les initialise tous par des valeurs nulles.
!       Et on demande a l'utilisateur de conserver cette initialisation.
!                                        =========

  do ikpdc = 1, 6
    do ielpdc = 1, ncepdp
      ckupdc(ielpdc,ikpdc) = 0.d0
    enddo
  enddo

  ! --- Tenseur diagonal
  !   Exemple de pertes de charges dans la direction x

  iutile = 0
  if (iutile.eq.0) return

  do ielpdc = 1, ncepdp
    iel=icepdc(ielpdc)
    vit = sqrt(  rtpa(iel,iu)**2                         &
         + rtpa(iel,iv)**2                         &
         + rtpa(iel,iw)**2)
    ckupdc(ielpdc,1) = 10.d0*vit
    ckupdc(ielpdc,2) =  0.d0*vit
    ckupdc(ielpdc,3) =  0.d0*vit
  enddo

  ! --- Tenseur 3x3
  !   Exemple de pertes de charges a ALPHA = 45 degres x,y
  !      la direction x resiste par ck1 et y par ck2
  !      ck2 nul represente des ailettes comme ceci :  ///////
  !      dans le repere de calcul X Y

  !                 Y|    /y
  !                  |  /
  !                  |/
  !                  \--------------- X
  !                   \ / ALPHA
  !                    \
  !                     \ x

  iutile = 0
  if (iutile.eq.0) return

  alpha  = pi/4.d0
  cosalp = cos(alpha)
  sinalp = sin(alpha)
  ck1 = 10.d0
  ck2 =  0.d0

  do ielpdc = 1, ncepdp
    iel=icepdc(ielpdc)
    vit = sqrt(  rtpa(iel,iu)**2                         &
         + rtpa(iel,iv)**2                         &
         + rtpa(iel,iw)**2)
    ckupdc(ielpdc,1) = (cosalp**2*ck1 + sinalp**2*ck2)*vit
    ckupdc(ielpdc,2) = (sinalp**2*ck1 + cosalp**2*ck2)*vit
    ckupdc(ielpdc,3) =  0.d0
    ckupdc(ielpdc,4) = cosalp*sinalp*(-ck1+ck2)*vit
    ckupdc(ielpdc,5) =  0.d0
    ckupdc(ielpdc,6) =  0.d0
  enddo

!-------------------------------------------------------------------------------

endif

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
