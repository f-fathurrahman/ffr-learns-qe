! ffr: 23 Juni 2016
!

!------------------------------------------------------------------------------
SUBROUTINE test_gvectors()
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE gvecw, ONLY : gcutw
  USE wvfct, ONLY : npwx, npw, g2kin
  USE gvect, ONLY : ngm, g, gg, gl, mill
  USE klist, ONLY : init_igk
  USE klist, ONLY : ngk, nks, xk, igk_k
  IMPLICIT NONE
  ! Local
  INTEGER :: ik, ig, jg

  ALLOCATE( ngk(nks) )
  CALL n_plane_waves( gcutw, nks, xk, g, ngm, npwx, ngk )

  ALLOCATE( g2kin(npwx) )

  WRITE(*,*) 'size(mill) = ', size(mill)
  CALL init_igk ( npwx, ngm, g, gcutw )
  DO ik = 1, nks
    IF(ionode) WRITE(stdout,*) 'ik, npw, ngk(ik) = ', ik, npw, ngk(ik)
    DO ig=1,npw
      DO jg=ig,npw
        WRITE(*,'(1x,2I5,3F13.5)') ig, jg, g(1:3, igk_k(ig,ik) ) - g(1:3, igk_k(jg,ik) )
      ENDDO
      !WRITE(*,'(1x,4I5,3F13.5)') ig, mill(1:3,igk(ig)), g(1:3,igk_k(ig,ik))
    ENDDO
  ENDDO

  WRITE(*,*) 'size(gg) = ', size(gg)
  WRITE(*,*) 'size(gl) = ', size(gl)

END SUBROUTINE


!=-----------------------------------------------------------------------------
PROGRAM main
!=-----------------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE mp_global, ONLY : mp_startup, mp_global_end
  IMPLICIT NONE

  CALL mp_startup()

  CALL setup_structure_H()

  CALL setup_fft( 30_DP, 120_DP )

  CALL setup_symmetry()

  CALL setup_kpoints_000()

  CALL setup_gvect( .FALSE. )

  !CALL test_gvectors()  ! not yet working
  !CALL t_import_gvect()
  !CALL t_import_gvecs()
  !CALL t_import_gvecw()

  CALL gvect_info() 

  CALL mp_global_end()

END PROGRAM
