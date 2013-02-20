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

subroutine usativ &
!================

 ( nvar   , nscal  ,                                              &
   nbmetd , nbmett , nbmetm ,                                     &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  , &
   tmprom , ztprom , zdprom , xmet   , ymet   , pmer   ,          &
   ttprom , qvprom , uprom  , vprom  , ekprom , epprom ,          &
   rprom  , tpprom , phprom )

!===============================================================================
! FONCTION :
! --------

! ROUTINE UTILISATEUR : INITIALISATION DES VARIABLES DE CALCUL
!     POUR LA PHYSIQUE PARTICULIERE: VERSION ATMOSPHERIQUE
!     PENDANT DE USINIV

! Cette routine est appelee en debut de calcul (suite ou non)
!     avant le debut de la boucle en temps

! Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
!     les variables de calcul,
!     les valeurs du pas de temps


! On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
!     ou relues d'un fichier suite
! On ne dispose des variables VISCLS, CP (quand elles sont
!     definies) que si elles ont pu etre relues dans un fichier
!     suite de calcul

! Les proprietes physiaues sont accessibles dans le tableau
!     PROPCE (prop au centre), PROPFA (aux faces internes),
!     PROPFB (prop aux faces de bord)
!     Ainsi,
!      PROPCE(IEL,IPPROC(IROM  )) designe ROM   (IEL)
!      PROPCE(IEL,IPPROC(IVISCL)) designe VISCL (IEL)
!      PROPCE(IEL,IPPROC(ICP   )) designe CP    (IEL)
!      PROPCE(IEL,IPPROC(IVISLS(ISCAL))) designe VISLS (IEL ,ISCAL)

!      PROPFA(IFAC,IPPROF(IFLUMA(IVAR ))) designe FLUMAS(IFAC,IVAR)

!      PROPFB(IFAC,IPPROB(IROM  )) designe ROMB  (IFAC)
!      PROPFB(IFAC,IPPROB(IFLUMA(IVAR ))) designe FLUMAB(IFAC,IVAR)





! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME USPHYV
!     ET PAS ICI


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa coefb      ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
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
use entsor
use parall
use period
use mesh

!===============================================================================

implicit none

integer          nvar   , nscal
integer          nbmetd , nbmett , nbmetm

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision tmprom(nbmetm)
double precision ztprom(nbmett) , zdprom(nbmetd)
double precision xmet(nbmetm)   , ymet(nbmetm)  , pmer(nbmetm)
double precision ttprom(nbmett,nbmetm) , qvprom(nbmett,nbmetm)
double precision uprom(nbmetd,nbmetm)  , vprom(nbmetd,nbmetm)
double precision ekprom(nbmetd,nbmetm) , epprom(nbmetd,nbmetm)
double precision rprom(nbmett,nbmetm)  , tpprom(nbmett,nbmetm)
double precision phprom(nbmett,nbmetm)

! Local variables

integer          iel, iutile
double precision d2s3
double precision zent,xuent,xvent,xkent,xeent,tpent

integer, allocatable, dimension(:) :: lstelt

!===============================================================================



!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))

d2s3 = 2.d0/3.d0

!===============================================================================
! 2. INITIALISATION DES INCONNUES :
!      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
!===============================================================================

! --- INCONNUES
!     Exemple : Initialisation a partir des profils meteo


if (isuite.eq.0) then

! Initialisation with the input meteo profile

  do iel = 1, ncel

    zent=xyzcen(3,iel)

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, uprom , zent  , ttcabs, xuent )

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, vprom , zent  , ttcabs, xvent )

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, ekprom, zent  , ttcabs, xkent )

    call intprf                                                   &
    !==========
   (nbmetd, nbmetm,                                               &
    zdprom, tmprom, epprom, zent  , ttcabs, xeent )

    rtp(iel,iu)=xuent
    rtp(iel,iv)=xvent
    rtp(iel,iw)=0.d0

!     ITYTUR est un indicateur qui vaut ITURB/10
    if    (itytur.eq.2) then

      rtp(iel,ik)  = xkent
      rtp(iel,iep) = xeent

    elseif(itytur.eq.3) then

      rtp(iel,ir11) = d2s3*xkent
      rtp(iel,ir22) = d2s3*xkent
      rtp(iel,ir33) = d2s3*xkent
      rtp(iel,ir12) = 0.d0
      rtp(iel,ir13) = 0.d0
      rtp(iel,ir23) = 0.d0
      rtp(iel,iep)  = xeent

    elseif(iturb.eq.50) then

      rtp(iel,ik)   = xkent
      rtp(iel,iep)  = xeent
      rtp(iel,iphi) = d2s3
      rtp(iel,ifb)  = 0.d0

    elseif(iturb.eq.60) then

      rtp(iel,ik)   = xkent
      rtp(iel,iomg) = xeent/cmu/xkent

    elseif(iturb.eq.70) then

      rtp(iel,inusa) = cmu*xkent**2/xeent

    endif

    if (iscalt.ge.0) then
! On suppose que le scalaire est la temperature potentielle :
      call intprf                                                 &
      !==========
   (nbmett, nbmetm,                                               &
    ztprom, tmprom, tpprom, zent  , ttcabs, tpent )

      rtp(iel,isca(iscalt)) = tpent

    endif
  enddo

endif

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
