! Adapted from init_vloc

!----------------------------------------------------------------------
SUBROUTINE init_vloc()
!----------------------------------------------------------------------
  !! This routine computes the fourier coefficient of the local
  !! potential vloc(ig,it) for each type of atom.
  !
  USE atom,           ONLY : msh, rgrid
  USE m_gth,          ONLY : vloc_gth
  USE kinds,          ONLY : DP
  USE uspp_param,     ONLY : upf
  USE ions_base,      ONLY : ntyp => nsp
  USE cell_base,      ONLY : omega, tpiba2
  USE vlocal,         ONLY : vloc
  USE gvect,          ONLY : ngl, gl
  !
  IMPLICIT NONE
  !
  INTEGER :: nt ! counter on atomic types
  
  !
  ! Cases of GTH pspot and Coulomb potential are removed
  !

  !
  vloc(:,:) = 0._DP
  !
  DO nt = 1, ntyp
    !
    ! normal case
    !
    CALL my_vloc_of_g( rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r, &
                    upf(nt)%vloc(:), upf(nt)%zp, tpiba2, ngl, gl, omega, &
                    vloc(:,nt) )
    write(*,*) 'init_vloc: nt = ', nt, ' sum vloc(:,nt) in Ha = ', sum(vloc(:,nt))*0.5d0
  ENDDO

END SUBROUTINE init_vloc

