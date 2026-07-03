!------------------------------------------------------------------------
SUBROUTINE my_exx_fft_create()
!------------------------------------------------------------------------
  !! Initialise the custom grid that allows to put the wavefunction
  !! onto the new (smaller) grid for \rho=\psi_{k+q}\psi^*_k and vice versa.  
  !! Set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  use kinds, only: dp
  !
  USE gvecw,          ONLY : ecutwfc
  USE gvect,          ONLY : ecutrho, ngm, g, gg, gstart, mill
  USE cell_base,      ONLY : at, bg, tpiba2
  USE recvec_subs,    ONLY : ggen, ggens
  USE fft_base,       ONLY : smap
  USE fft_types,      ONLY : fft_type_init
  USE symm_base,      ONLY : fft_fact
  USE mp_exx,         ONLY : nproc_egrp, negrp, intra_egrp_comm
  USE mp_bands,       ONLY : nproc_bgrp, intra_bgrp_comm, nyfft
  !
  USE klist,          ONLY : nks, xk
  USE mp_pools,       ONLY : inter_pool_comm
  USE mp,             ONLY : mp_max, mp_sum
  !
  USE control_flags,  ONLY : tqr
  USE realus,         ONLY : qpointlist, tabxx, tabp
  USE command_line_options, ONLY : nmany_
  !
  ! from exx module
  USE control_flags,        ONLY : gamma_only, tqr
  USE fft_types,            ONLY : fft_type_descriptor
  USE stick_base,           ONLY : sticks_map, sticks_map_deallocate
  !
  use exx, only: exx_fft_initialized, ecutfock, dfftt, gstart_t, ngmt_g, npwt, gt, ggt

  !
  IMPLICIT NONE
  !
  ! local variables
  !
  INTEGER :: ik, ngmt
  INTEGER, ALLOCATABLE :: ig_l2gt(:), millt(:,:)
  INTEGER, EXTERNAL :: n_plane_waves
  REAL(DP) :: gkcut, gcutmt
  LOGICAL :: lpara

  write(*,*)
  write(*,*) '<div> ENTER my_exx_fft_create'
  write(*,*)

  IF( exx_fft_initialized ) then
    write(*,*) 'No need to do anythin, exx_fft is already initialized'
    write(*,*) 'Early return in my_exx_fft_create'
    write(*,*)
    write(*,*) '</div> EXIT my_exx_fft_create'
    write(*,*)
    RETURN
  ENDIF
  !
  ! Initialise the custom grid that allows us to put the wavefunction
  ! onto the new (smaller) grid for \rho=\psi_{k+q}\psi^*_k and vice versa
  !
  ! gkcut is such that all |k+G|^2 < gkcut (in units of (2pi/a)^2)
  ! Note that with k-points, gkcut > ecutwfc/(2pi/a)^2
  ! gcutmt is such that |q+G|^2 < gcutmt
  !
  IF( gamma_only ) THEN
    gkcut = ecutwfc / tpiba2
    gcutmt = ecutfock / tpiba2
  ELSE
    gkcut = 0.0_DP
    DO ik = 1, nks
      gkcut = MAX(gkcut, SQRT(SUM(xk(:,ik)**2)))
    ENDDO
    CALL mp_max( gkcut, inter_pool_comm )
    ! Alternatively, variable "qnorm" earlier computed in "exx_grid_init"
    ! could be used as follows:
    ! gkcut = ( SQRT(ecutwfc/tpiba2) + qnorm )**2
    gkcut = ( SQRT(ecutwfc/tpiba2) + gkcut )**2
    ! 
    ! the following instruction may be needed if ecutfock \simeq ecutwfc
    ! and guarantees that all k+G are included
    !
    gcutmt = MAX(ecutfock/tpiba2, gkcut)
  ENDIF
  write(*,*) 'ecutwfc = ', ecutwfc
  write(*,*) 'ecutfock = ', ecutfock
  write(*,*) 'gcutmt = ', gcutmt
  !
  ! set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  ! ffr: negrp must be 1, band parallelization is removed
  !
  ! no band parallelization: exx grid is a subgrid of general grid
  !
  lpara = ( nproc_bgrp > 1 )
  CALL fft_type_init( dfftt, smap, "rho", gamma_only, lpara,         &
                    intra_bgrp_comm, at, bg, gcutmt, gcutmt/gkcut, &
                    fft_fact=fft_fact, nyfft=nyfft, nmany=nmany_ )
  CALL ggens( dfftt, gamma_only, at, g, gg, mill, gcutmt, ngmt, gt, ggt )
  gstart_t = gstart
  npwt = n_plane_waves(ecutwfc/tpiba2, nks, xk, gt, ngmt)
  ngmt_g = ngmt
  CALL mp_sum( ngmt_g, intra_bgrp_comm )
  !
  ! define clock labels (this enables the corresponding fft too)
  dfftt%rho_clock_label = 'fftc'
  dfftt%wave_clock_label = 'fftcw' 
  !
  WRITE(*, '(/5x,"EXX grid: ",i8," G-vectors", 5x,       &
    &   "FFT dimensions: (",i4,",",i4,",",i4,")")') ngmt_g, &
    &   dfftt%nr1, dfftt%nr2, dfftt%nr3
  !
  exx_fft_initialized = .TRUE.
  !
  IF (tqr) THEN
    IF (ecutfock == ecutrho) THEN
      WRITE(*, '(5x,"Real-space augmentation: EXX grid -> DENSE grid")' )
      tabxx => tabp
    ELSE
      WRITE(*, '(5x,"Real-space augmentation: initializing EXX grid")' )
      CALL qpointlist( dfftt, tabxx )
    ENDIF
  ENDIF

  write(*,*)
  write(*,*) '</div> EXIT my_exx_fft_create'
  write(*,*)
  
  RETURN
  !
END SUBROUTINE
