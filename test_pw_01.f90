PROGRAM my_pwscf
  
  USE mp_global, ONLY : mp_startup
  USE mp_world, ONLY : world_comm
  USE mp_pools, ONLY : intra_pool_comm
  USE mp_diag, ONLY : mp_start_diag
  USE mp_bands, ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE command_line_options, ONLY : ndiag_, input_file_
  USE environment, ONLY : environment_start
  USE read_input, ONLY : read_input_file

  IMPLICIT NONE 
  
  INTEGER :: exit_status

  CALL mp_startup(start_images=.true.)

  CALL mp_start_diag( ndiag_, world_comm, intra_bgrp_comm, &
                      do_distr_diag_inside_bgrp_=.true. )

  CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, inter_bgrp_comm )

  CALL environment_start( 'PWSCF' )

  CALL read_input_file( 'PW', input_file_ )

  CALL my_run_pwscf( exit_status )

  !CALL unset_mpi_comm_4_solvers()
  !CALL stop_run( exit_status )
  !CALL do_stop( exit_status )

  WRITE(*,*) 'Program ended normally'

END PROGRAM 


SUBROUTINE my_run_pwscf( exit_status )

  IMPLICIT NONE 
  INTEGER :: exit_status

  exit_status = 0

  CALL iosys()
  
  CALL setup()

  CALL init_run()

  CALL test_uspp_param()

  CALL test_vlocal()

  CALL test_v_hartree()

END SUBROUTINE 


SUBROUTINE test_uspp_param()
  USE uspp_param, ONLY : upf
  
  IMPLICIT NONE 

  WRITE(*,*) 'size upf = ', size(upf)
END SUBROUTINE 

SUBROUTINE test_vlocal()

  USE vlocal, ONLY : vloc, strf
  USE gvect, ONLY : ngl, ngm
  USE ions_base, ONLY : ntyp => nsp
  USE fft_base,  ONLY : dfftp
  USE gvect, ONLY : ngm
  USE gvect, ONLY : igtongl, gcutm
  USE scf, ONLY : v_of_0
  USE fft_interfaces, ONLY : invfft
  USE cell_base, ONLY : tpiba2

  IMPLICIT NONE 

  REAL(8), PARAMETER :: PI=4.d0*atan(1.d0)
  INTEGER :: nt, ng
  
  COMPLEX(8), ALLOCATABLE :: aux(:)
  REAL(8), ALLOCATABLE :: vloc_r(:)

  ALLOCATE( aux(dfftp%nnr) )
  ALLOCATE( vloc_r(dfftp%nnr) )

  aux(:) = (0.d0,0.d0)
  vloc_r(:) = 0.d0


  WRITE(*,*) 'dfftp%nnr = ', dfftp%nnr

  WRITE(*,*) 'ngl = ', ngl  ! no. of g-shells
  WRITE(*,*) 'ngm = ', ngm
  WRITE(*,*) 'size vloc = ', size(vloc)

  DO nt = 1, ntyp
      DO ng = 1, ngm
          aux(dfftp%nl(ng)) = aux(dfftp%nl(ng)) + vloc(igtongl(ng),nt)*strf(ng, nt)
      END DO
  END DO

  WRITE(*,*) 'v_of_0 = ', v_of_0

  CALL invfft('Rho', aux, dfftp)
  vloc_r = dble(aux)

  WRITE(*,*) 'sum vloc_r in Ha', sum(vloc_r)/2.d0
  WRITE(*,*) 'vloc_r', vloc_r(1:5)

  WRITE(*,*) 'gcutm = ', gcutm*tpiba2
  WRITE(*,*) 'PI    = ', PI
  WRITE(*,*) 4.d0*15.d0/2.d0/PI**2

END SUBROUTINE 


SUBROUTINE test_v_hartree()

  USE gvect, ONLY : ngm
  USE fft_base, ONLY : dfftp
  USE lsda_mod,  ONLY : nspin
  IMPLICIT NONE
  
  COMPLEX(8), ALLOCATABLE :: rhog(:)
  REAL(8), ALLOCATABLE :: rhor(:)
  REAL(8), ALLOCATABLE :: v(:,:)
  REAL(8) :: ehart, charge

  ALLOCATE( rhor(dfftp%nnr) )
  ALLOCATE( rhog(ngm) )
  ALLOCATE( v(dfftp%nnr,nspin) )

  rhog(:) = (1.d0,0.d0)
  rhog(4) = (2.d0,0.d0)
  rhog(3) = (4.d0,0.d0)
  v(:,:) = 0.d0

  call v_h( rhog, ehart, charge, v )
  WRITE(*,*)
  WRITE(*,*) 'From reciprocal space rhoe'
  WRITE(*,*) 'ehart  = ', ehart
  WRITE(*,*) 'charge = ', charge

  call my_v_h( rhog, ehart, charge, v )
  WRITE(*,*)
  WRITE(*,*) 'From reciprocal space rhoe my version'
  WRITE(*,*) 'ehart  = ', ehart
  WRITE(*,*) 'charge = ', charge

  
  rhor(:) = 1.d0
  rhor(4) = 2.d0
  rhor(3) = 4.d0
  call v_h_of_rho_r( rhor, ehart, charge, v(:,1) )
  WRITE(*,*)
  WRITE(*,*) 'From real space rhoe'
  WRITE(*,*) 'ehart  = ', ehart
  WRITE(*,*) 'charge = ', charge

  DEALLOCATE( rhor )
  DEALLOCATE( rhog )
  DEALLOCATE( v )
END SUBROUTINE


SUBROUTINE my_v_h( rhog, ehart, charge, v )
  ! ... Hartree potential VH(r) from n(G)
  USE constants, ONLY : fpi, e2
  USE kinds,     ONLY : DP
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces,ONLY : invfft
  USE gvect,     ONLY : ngm, gg, gstart
  USE lsda_mod,  ONLY : nspin
  USE cell_base, ONLY : omega, tpiba2

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN)  :: rhog(ngm)
  REAL(DP),  INTENT(INOUT) :: v(dfftp%nnr,nspin)
  REAL(DP),    INTENT(OUT) :: ehart, charge
  !
  REAL(DP)              :: fac
  REAL(DP), ALLOCATABLE :: aux1(:,:)
  REAL(DP) :: rgtot_re, rgtot_im
  INTEGER :: is, ig
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  !
  ALLOCATE( aux(dfftp%nnr), aux1(2,ngm) )


  WRITE(*,*) '----------------'
  WRITE(*,*) 'Starting my_v_h:'
  WRITE(*,*) '----------------'

  charge = 0.D0
  IF ( gstart == 2 ) THEN
     charge = omega*REAL( rhog(1) )
     WRITE(*,*) 'From my_v_h: charge = ', charge
  END IF
  
  ehart     = 0.D0
  aux1(:,:) = 0.D0

  WRITE(*,*) 'gstart = ', gstart
  DO ig = gstart, ngm
    fac = 1.D0 / gg(ig) 
    rgtot_re = REAL(  rhog(ig) )
    rgtot_im = AIMAG( rhog(ig) )
    ehart = ehart + ( rgtot_re**2 + rgtot_im**2 ) * fac
    aux1(1,ig) = rgtot_re * fac
    aux1(2,ig) = rgtot_im * fac
  ENDDO

  fac = e2 * fpi / tpiba2
  ehart = ehart * fac
  aux1 = aux1 * fac
  ehart = ehart * 0.5D0 * omega
  aux(:) = 0.D0
  aux(dfftp%nl(1:ngm)) = CMPLX ( aux1(1,1:ngm), aux1(2,1:ngm), KIND=8 )

  ! transform hartree potential to real space
  CALL invfft ('Rho', aux, dfftp)
  
  ! add hartree potential to the xc potential
  DO is = 1, nspin
    v(:,is) = v(:,is) + DBLE (aux(:))
  END DO
  !
  DEALLOCATE( aux, aux1 )
  RETURN
END SUBROUTINE my_v_h


