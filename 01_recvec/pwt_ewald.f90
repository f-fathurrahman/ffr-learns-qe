SUBROUTINE ewald_debug()
  USE kinds, ONLY : DP
  USE cell_base, ONLY : alat, at, bg, omega
  USE ions_base, ONLY : zv, nat, nsp, ityp, tau
  USE gvect, ONLY : ngm, gstart, g, gg, gcutm
  USE cell_base, ONLY : tpiba2
  USE gvecs, ONLY : dual
  USE gvecw, ONLY : ecutwfc
  USE vlocal, ONLY : strf
  USE control_flags, ONLY : gamma_only
  IMPLICIT NONE
  REAL(DP), EXTERNAL :: ewald
  REAL(DP) :: Ene_Ewald

  WRITE(*,*)
  WRITE(*,*) 'sum(strf)  = ', sum(strf)
  WRITE(*,*) 'gstart     = ', gstart
  WRITE(*,*) 'gamma_only = ', gamma_only
  WRITE(*,*) 'dual       = ', dual
  WRITE(*,*) 'tpiba2     = ', tpiba2
  WRITE(*,*) 'ecutwf     = ', ecutwfc
  WRITE(*,*) 'gcutm      = ', gcutm
  WRITE(*,*) 'omega      = ', omega
  WRITE(*,*) 'alat       = ', alat
  WRITE(*,*) 'shape(ityp) = ', shape(ityp)
  WRITE(*,*) 'shape(zv)   = ', shape(zv)
  WRITE(*,*) 'zv = ', zv

  Ene_Ewald = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )

  WRITE(*,'(/,1x,A,F18.10)') 'Ene_Ewald (in Ha) = ', Ene_Ewald*0.5_DP

END SUBROUTINE


PROGRAM pwt_ewald
  USE mp_global, ONLY: mp_startup
  USE environment, ONLY: environment_start
  USE read_input, ONLY: read_input_file
  USE check_stop, ONLY: check_stop_init
  USE parameters, ONLY: ntypx, npk, lmaxx
  ! for v-5.1 and later
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  ! Local
  INTEGER :: exit_status

  CALL mp_startup( )

  CALL environment_start('PWSCF')

  CALL read_input_file('PW', input_file_)

  CALL iosys()

  ! If we don't include the following two subroutine calls, there will
  ! be error before the program stops regarding find_unit.
  CALL setup()
  CALL init_run()

  CALL ewald_debug()

  CALL stop_run( exit_status )
  CALL do_stop( exit_status )

END PROGRAM
