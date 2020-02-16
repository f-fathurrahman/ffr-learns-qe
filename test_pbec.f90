PROGRAM test_pbec
  IMPLICIT NONE 
  REAL(8) :: rho, grho
  REAL(8) :: sc, v1c, v2c
  INTEGER :: iflag

  rho = 1.2d0
  grho = 0.1d0
  iflag = 1
  CALL pbec( rho, grho, iflag, sc, v1c, v2c )
  WRITE(*,'(1x,A,F18.10)') 'QE:    ', sc

  CALL pbec_libxc( rho, grho, iflag, sc, v1c, v2c )
  WRITE(*,'(1x,A,F18.10)') 'Libxc: ', sc

END PROGRAM 


SUBROUTINE pbec_libxc(rho, grho, iflag, sc, v1c, v2c)
  !
  ! PBE correlation (without LDA part)
  ! iflag=1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
  ! iflag=2: J.P.Perdew et al., PRL 100, 136406 (2008).
  ! iflag=3: L. Chiodo et al, PRL 108, 126402 (2012)  (PBEQ2D)
  !
  use xc_f90_types_m
  use xc_f90_lib_m
  IMPLICIT NONE 
  INTEGER, PARAMETER :: DP=8
  INTEGER, INTENT(IN) :: iflag
  REAL(DP), INTENT(IN) :: rho, grho
  REAL(DP), INTENT(OUT):: sc, v1c, v2c
  ! local variables
  integer :: func_id = -1 ! not set
  integer :: N = 1
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info
  real(dp) :: exc, ec_lda = 0.0d0 , vc_lda = 0.0d0

  IF( iflag == 1)  func_id = 130
  IF( iflag == 2)  func_id = 133
  IF( iflag == 3)  CALL errore('pbec', 'case not implemented with libxc', iflag)

  CALL xc_f90_func_init( xc_func, xc_info, func_id, XC_UNPOLARIZED )
  call xc_f90_gga_exc_vxc( xc_func, N, rho, grho, exc, v1c, v2c )
  call xc_f90_func_end( xc_func )

  ! remove PW correlation for compatibility with QE  
  call xc_f90_func_init( xc_func, xc_info, 12, XC_UNPOLARIZED )    
  call xc_f90_lda_exc_vxc( xc_func, N, rho, ec_lda, vc_lda )
  call xc_f90_func_end( xc_func ) 
  
  exc = exc - ec_lda
  v1c = v1c - vc_lda
  
  sc = exc * rho  ! e_x = rho * \epsilon_x
  v2c = v2c*2.0_dp

END SUBROUTINE 



