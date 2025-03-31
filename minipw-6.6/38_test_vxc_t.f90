program test_vxc_t
  USE mp_global,         ONLY : mp_startup, mp_global_end
  USE environment,       ONLY : environment_start
  USE command_line_options, ONLY: input_file_
  !
  IMPLICIT NONE
  CHARACTER(LEN=9) :: code = 'LD1'
  integer  :: lsd ! 1 in the LSDA case, 0 otherwise
  real(8) :: rho(2), rhoc ! the system density
  real(8) :: exc(1), vxc(2)

  CALL mp_startup()
  CALL environment_start(code)
  CALL ld1_readin(input_file_)
  CALL ld1_setup()

  lsd = 0
  rho(1) = 2.d0
  rho(2) = 0.d0  
  rhoc = 0.d0

  call vxc_t(lsd, rho, rhoc, exc, vxc)

  write(*,*) 'exc (in Ha) = ', exc*0.5d0
  write(*,*) 'vxc (in Ha) = ', vxc*0.5d0

end program

