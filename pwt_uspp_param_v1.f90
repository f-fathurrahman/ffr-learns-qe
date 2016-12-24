! efefer, October 1st, 2014

PROGRAM t_uspp_param
  USE io_global, ONLY: stdout
  USE mp_global, ONLY: mp_startup
  USE environment, ONLY: environment_start
  USE read_input, ONLY: read_input_file
  USE check_stop, ONLY: check_stop_init
  USE parameters, ONLY: ntypx, npk, lmaxx
  USE uspp_param, ONLY: nvb, newpseudo, upf
  ! for v-5.1 and later
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  ! Local
  INTEGER :: ips
  INTEGER :: exit_status

  CALL mp_startup( )

  CALL environment_start('PWSCF')

  CALL read_input_file('PW', input_file_)

  CALL iosys()
  ! If we don't include the following two subroutine calls, there will
  ! be error before the program stops regarding find_unit.
  CALL setup()
  CALL init_run()

  WRITE(stdout,*) 'nvb = ', nvb
  WRITE(stdout,*) 'newpseudo = ', newpseudo

  DO ips=1,size(upf)
    WRITE(stdout,*)
    WRITE(stdout,'(1x,A,I2)') 'UPF information for species #', ips
    WRITE(stdout,*) '*******************************'
    WRITE(stdout,*)
    WRITE(stdout,'(1x,A,F4.1)') 'Valence = ', upf(ips)%zp
    WRITE(stdout,*) 'lmax = ', upf(ips)%lmax
    WRITE(stdout,*) 'lloc = ', upf(ips)%lloc
    WRITE(stdout,*) 'no. of projectors = ', upf(ips)%nbeta
    WRITE(stdout,*) 'shape(upf%vnl) = ', shape( upf(ips)%vnl )
    WRITE(stdout,*) 'size(vloc) = ', size( upf(ips)%vloc )
    WRITE(stdout,*) 'shape(dion) = ', shape( upf(ips)%dion )
    WRITE(stdout,*) 'upf%dion = ', upf(ips)%dion
    WRITE(stdout,*) 'shape(pswfc) = ', shape( upf(ips)%pswfc )
    WRITE(stdout,*) 'shape(aewfc) = ', shape( upf(ips)%aewfc )
    WRITE(stdout,*) 'shape(beta) = ', shape( upf(ips)%beta )
  ENDDO

  CALL stop_run( exit_status )
  CALL do_stop( exit_status )

END PROGRAM

