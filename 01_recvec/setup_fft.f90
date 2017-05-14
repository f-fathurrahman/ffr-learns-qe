!------------------------------------------------------------------------------
SUBROUTINE setup_fft( ecutwfc_, ecutrho_ )
!------------------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE io_global, ONLY: stdout, ionode
  USE cell_base, ONLY: tpiba2, at, bg
  USE gvecw, ONLY: ecutwfc, gcutw
  USE gvect, ONLY: ecutrho, gcutm
  USE gvecs, ONLY: ecuts, gcutms, dual, doublegrid
  USE fft_base, ONLY: dfftp, dffts
  USE fft_types, ONLY: fft_type_allocate, fft_type_init
  USE mp_bands, ONLY: intra_bgrp_comm
  USE paw_variables, ONLY: okpaw
  USE uspp, ONLY: okvan
  IMPLICIT NONE
  REAL(DP) :: ecutwfc_, ecutrho_

  ! FIXME
  okvan = .FALSE.
  okpaw = .FALSE.

  ecutwfc = ecutwfc_
  ecutrho = ecutrho_
  dual    = ecutrho/ecutwfc

  doublegrid = ( dual > 4.D0 )
  IF ( doublegrid .AND. (.NOT.okvan .AND. .not.okpaw) ) &
     CALL infomsg ( 'setup', 'no reason to have ecutrho>4*ecutwfc' )

  gcutm = dual * ecutwfc / tpiba2
  gcutw = ecutwfc / tpiba2
  !
  IF ( doublegrid ) THEN
    gcutms = 4.D0 * ecutwfc / tpiba2
  ELSE
    gcutms = gcutm
  END IF

  IF( ionode ) THEN
    WRITE(stdout,'(/,1x,A,F18.5)') 'dual = ', dual
    WRITE(stdout,'(1x,A,2F18.5)') 'gcutm, gcutms = ', gcutm, gcutms
  ENDIF

  ! ... calculate dimensions of the FFT grid
  !
  ! ... if the smooth and dense grid must coincide, ensure that they do
  ! ... also if dense grid is set from input and smooth grid is not
  !
  IF ( ( dfftp%nr1 /= 0 .AND. dfftp%nr2 /= 0 .AND. dfftp%nr3 /= 0 ) .AND. &
       ( dffts%nr1 == 0 .AND. dffts%nr2 == 0 .AND. dffts%nr3 == 0 ) .AND. &
       .NOT. doublegrid ) THEN
     dffts%nr1 = dfftp%nr1
     dffts%nr2 = dfftp%nr2
     dffts%nr3 = dfftp%nr3
  END IF
  CALL fft_type_allocate( dfftp, at, bg, gcutm, intra_bgrp_comm )
  CALL fft_type_allocate( dffts, at, bg, gcutms, intra_bgrp_comm )

  IF( ionode ) THEN
    WRITE(stdout,'(/,1x,A,3I8)') 'Dense : nr1,nr2,nr3', dfftp%nr1, dfftp%nr2, dfftp%nr3
    WRITE(stdout,'(1x,A,3I8)') 'Smooth: nr1,nr2,nr3', dffts%nr1, dffts%nr2, dffts%nr3
  ENDIF
END SUBROUTINE
