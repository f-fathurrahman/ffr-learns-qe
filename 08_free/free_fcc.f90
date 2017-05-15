! ffr: 7 Dec 2015
!
! An example of how to work with gvect: calculation of empty-lattice
! band structure of FCC
!


!------------------------------------------------------------------------------
SUBROUTINE setup_structure()
!------------------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE io_global, ONLY: stdout, ionode
  USE cell_base, ONLY: cell_base_init
  USE cell_base, ONLY: omega, at, bg, alat
  USE ions_base, ONLY : if_pos, ityp, tau, extfor, fixatom, nsp, &
                        atm, nat, amass, tau_format
  IMPLICIT NONE
  !
  INTEGER :: ibrav_
  REAL(8) :: celldm_(6)
  REAL(8) :: a_, b_, c_, cosab_, cosac_, cosbc_
  LOGICAL :: trd_ht
  REAL(8) :: rd_ht(3,3)
  CHARACTER(len=10) :: cell_units_
  !
  INTEGER :: ii, ia
  INTEGER :: iuxsf

  ! Parameters required to initialize unit cell
  ibrav_ = 1
  celldm_(:) = 0.0_DP  ! in bohr
  a_ = 5_DP   ! in angstrom
  b_ = 0.0_DP
  c_ = 0.0_DP
  cosab_ = 0.0_DP
  cosac_ = 0.0_DP
  cosbc_ = 0.0_DP
  trd_ht = .FALSE.
  rd_ht(:,:) = 0.0_DP  ! unit cell
  cell_units_ = 'bohr'

  ! The variables defined in module cell_base are initialized here
  CALL cell_base_init( ibrav_, celldm_, a_, b_, c_, cosab_, cosac_, &
          & cosbc_, trd_ht, rd_ht, cell_units_ )
  
  ! Atomic positions
  !
  nsp = 1
  nat = 1
  ! We need to allocate several allocatable arrays that depend on nsp and nat
  ! Note that, we don't need to allocate `atm` and amass` because
  ! they are static arrays with size `ntypmax`.
  ALLOCATE( ityp(nat) )
  ALLOCATE( tau(3,nat) )
  ALLOCATE( if_pos(3,nat) )
  ALLOCATE( extfor(3,nat) )

  ityp(1) = 1
  !
  atm(1) = 'Si'
  !
  tau_format = 'crystal'
  !
  tau(:,1) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
  !
  if_pos(:,:) = 1
  extfor(:,:) = 0.0_DP
  !
  CALL convert_tau(tau_format, nat, tau)
  !
  IF(ionode) THEN
    WRITE(stdout, '(/,5x,"site n.     atom                  positions (alat units)")')
    WRITE(stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
               (ia, atm(ityp(ia)), ia, (tau(ii,ia), ii=1,3), ia=1,nat)
  ENDIF

  iuxsf = 111
  IF(ionode) THEN
    OPEN(unit=iuxsf, file='STRUCT.xsf', action='write')
    CALL xsf_struct(alat, at, nat, tau, atm, ityp, iuxsf)
    CLOSE(iuxsf)
  ENDIF

END SUBROUTINE



!------------------------------------------------------------------------------
SUBROUTINE setup_fft()
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE cell_base, ONLY : tpiba2, at, bg
  USE gvecw, ONLY : ecutwfc, gcutw
  USE gvect, ONLY : ecutrho, gcutm
  USE gvecs, ONLY : ecuts, gcutms, dual
  USE fft_base, ONLY : dfftp, dffts
  USE fft_types, ONLY : fft_type_allocate
  USE mp_bands, ONLY : intra_bgrp_comm
  IMPLICIT NONE
  
  ecutwfc = 10_DP
  ecutrho = 40_DP
  dual    = ecutrho/ecutwfc
  gcutms  = 4_DP*ecutwfc/tpiba2
  gcutm   = dual*ecutwfc/tpiba2
  gcutw   = ecutwfc / tpiba2

  IF( ionode ) THEN
    WRITE(stdout,'(/,1x,A,F18.5)') 'dual = ', dual
    WRITE(stdout,'(1x,A,2F18.5)') 'gcutm, gcutms = ', gcutm, gcutms
  ENDIF

  CALL fft_type_allocate( dfftp, at, bg, gcutm, intra_bgrp_comm )
  CALL fft_type_allocate( dffts, at, bg, gcutms, intra_bgrp_comm )

  IF( ionode ) THEN
    WRITE(stdout,'(/,1x,A,3I8)') 'Dense : nr1,nr2,nr3', dfftp%nr1, dfftp%nr2, dfftp%nr3
    WRITE(stdout,'(1x,A,3I8)') 'Smooth: nr1,nr2,nr3', dffts%nr1, dffts%nr2, dffts%nr3
  ENDIF
END SUBROUTINE



!------------------------------------------------------------------------------
SUBROUTINE setup_symmetry()
!------------------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE io_global, ONLY : stdout, ionode
  USE symm_base, ONLY : set_sym_bl, find_sym
  USE ions_base, ONLY : nat, tau, ityp
  USE fft_base, ONLY : dfftp
  USE symm_base, ONLY : nrot, nsym
  IMPLICIT NONE
  ! Local
  LOGICAL :: magnetic_sym, monopole
  REAL(DP) :: m_loc(3,nat)

  ! We ignore the magnetization
  magnetic_sym = .FALSE.
  !
  monopole = .FALSE.
  m_loc(:,:)   = 0_DP

  CALL set_sym_bl()
  CALL find_sym( nat, tau, ityp, magnetic_sym, m_loc, monopole )

  IF(ionode) THEN
    WRITE(stdout,*) 'nat, nsym, nrot = ', nat, nsym, nrot
  ENDIF

END SUBROUTINE



! An example of initializing k-points for band-structure calculation
!------------------------------------------------------------------------------
SUBROUTINE setup_kpoints()
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE parameters, ONLY : npk
  USE cell_base, ONLY : at, bg
  USE klist, ONLY : nks, wk, xk, nkstot
  IMPLICIT NONE
  !
  INTEGER :: ik, ii

  ! k-points for band structure calculation
  ! TODO: read from file
  nks = 60
  nkstot = nks
  wk(1:nks) = 1.0_DP  ! it is not really used in band structure calculation
  xk(:,1:nks) = reshape( (/ &
    0.00000000, 0.00000000, 0.50000000, &
    0.00000000, 0.00000000, 0.45454545, &
    0.00000000, 0.00000000, 0.40909091, &
    0.00000000, 0.00000000, 0.36363636, &
    0.00000000, 0.00000000, 0.31818182, &
    0.00000000, 0.00000000, 0.27272727, &
    0.00000000, 0.00000000, 0.22727273, &
    0.00000000, 0.00000000, 0.18181818, &
    0.00000000, 0.00000000, 0.13636364, &
    0.00000000, 0.00000000, 0.09090909, &
    0.00000000, 0.00000000, 0.04545455, &
    0.00000000, 0.00000000, 0.00000000, &
    0.02500000, 0.02500000, 0.02500000, &
    0.05000000, 0.05000000, 0.05000000, &
    0.07500000, 0.07500000, 0.07500000, &
    0.10000000, 0.10000000, 0.10000000, &
    0.12500000, 0.12500000, 0.12500000, &
    0.15000000, 0.15000000, 0.15000000, &
    0.17500000, 0.17500000, 0.17500000, &
    0.20000000, 0.20000000, 0.20000000, &
    0.22500000, 0.22500000, 0.22500000, &
    0.25000000, 0.25000000, 0.25000000, &
    0.27500000, 0.27500000, 0.27500000, &
    0.30000000, 0.30000000, 0.30000000, &
    0.32500000, 0.32500000, 0.32500000, &
    0.35000000, 0.35000000, 0.35000000, &
    0.37500000, 0.37500000, 0.37500000, &
    0.40000000, 0.40000000, 0.40000000, &
    0.42500000, 0.42500000, 0.42500000, &
    0.45000000, 0.45000000, 0.45000000, &
    0.47500000, 0.47500000, 0.47500000, &
    0.50000000, 0.50000000, 0.50000000, &
    0.45833333, 0.50000000, 0.50000000, &
    0.41666667, 0.50000000, 0.50000000, &
    0.37500000, 0.50000000, 0.50000000, &
    0.33333333, 0.50000000, 0.50000000, &
    0.29166667, 0.50000000, 0.50000000, &
    0.25000000, 0.50000000, 0.50000000, &
    0.20833333, 0.50000000, 0.50000000, &
    0.16666667, 0.50000000, 0.50000000, &
    0.12500000, 0.50000000, 0.50000000, &
    0.08333333, 0.50000000, 0.50000000, &
    0.04166667, 0.50000000, 0.50000000, &
    0.00000000, 0.50000000, 0.50000000, &
    0.00000000, 0.46875000, 0.46875000, &
    0.00000000, 0.43750000, 0.43750000, &
    0.00000000, 0.40625000, 0.40625000, &
    0.00000000, 0.37500000, 0.37500000, &
    0.00000000, 0.34375000, 0.34375000, &
    0.00000000, 0.31250000, 0.31250000, &
    0.00000000, 0.28125000, 0.28125000, &
    0.00000000, 0.25000000, 0.25000000, &
    0.00000000, 0.21875000, 0.21875000, &
    0.00000000, 0.18750000, 0.18750000, &
    0.00000000, 0.15625000, 0.15625000, &
    0.00000000, 0.12500000, 0.12500000, &
    0.00000000, 0.09375000, 0.09375000, &
    0.00000000, 0.06250000, 0.06250000, &
    0.00000000, 0.03125000, 0.03125000, &
    0.00000000, 0.00000000, 0.00000000 /), (/ 3,nks /) )

  CALL cryst_to_cart( nkstot, xk, bg, 1 )
  IF(ionode) THEN
    WRITE(stdout,*) 'npk, nks = ', npk, nks
    WRITE(stdout, '(23x,"cart. coord. in units 2pi/alat")')
    DO ik=1,nkstot
      WRITE(stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') ik, &
        ( xk(ii,ik), ii=1,3), wk(ik)
    ENDDO
  ENDIF

END SUBROUTINE


!------------------------------------------------------------------------------
SUBROUTINE setup_gvect()
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE recvec_subs, ONLY : ggen
  USE cell_base, ONLY : at, bg
  IMPLICIT NONE
  !
  LOGICAL :: gamma_only
  
  gamma_only = .FALSE.
  CALL data_structure( gamma_only )
  CALL ggen( gamma_only, at, bg )
END SUBROUTINE


!------------------------------------------------------------------------------
SUBROUTINE test_hamiltonian1()
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE cell_base, ONLY : tpiba2
  USE gvecw, ONLY : ecutwfc, gcutw
  USE wvfct, ONLY : npwx, npw, current_k, g2kin
  USE gvect, ONLY : ngm, g
  USE klist, ONLY : ngk, nks, xk, igk_k, init_igk
  IMPLICIT NONE
  ! Local
  INTEGER :: ik, ig, ib
  ! Full Hamiltonian and eigenvectors (column-wise)
  COMPLEX(8), ALLOCATABLE :: H_k(:,:), psi_k(:,:)
  !
  INTEGER :: nbnd_plot = 5
  REAL(8), ALLOCATABLE :: e_k(:), e_bnd(:,:), distk(:)
  !
  INTEGER :: n_plane_waves

  ALLOCATE( ngk(nks) )
  npwx = n_plane_waves( gcutw, nks, xk, g, ngm )
  WRITE(*,*) 'npwx = ', npwx

  ALLOCATE( igk_k(npwx,nks) )
  ALLOCATE( g2kin(npwx) )

  ALLOCATE( H_k(npwx,npwx) )
  ALLOCATE( psi_k(npwx,npwx) )
  ALLOCATE( e_k(npwx) )
  ALLOCATE( e_bnd(nks,nbnd_plot) )
  ALLOCATE( distk(nks) )

  CALL init_igk( npwx, ngm, g, gcutw )

  DO ik=1,nks
    ! generate Gk and sort them
!    WRITE(*,*) 'npw, ngm = ', npw, ngm
!    CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk_k(:,ik), g2kin )
!    CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk_k(:,ik), g2kin )
    npw = ngk(ik)
    IF(ionode) WRITE(stdout,*) 'ik, npw, ngk(ik) = ', ik, npw, ngk(ik)
    ! Build the Hamiltonian
    H_k(:,:) = cmplx(0_DP, 0_DP)
    DO ig=1,npw
      H_k(ig,ig) = ( ( xk(1,ik) + g(1,igk_k(ig,ik)) )**2 + &
                     ( xk(2,ik) + g(2,igk_k(ig,ik)) )**2 + &
                     ( xk(3,ik) + g(3,igk_k(ig,ik)) )**2 ) * tpiba2
    ENDDO
    WRITE(*,*) 'sum(H_k) = ', sum(H_k)
!    STOP 
    CALL cdiagh( npw, H_k, npwx, e_k, psi_k )
    ! For plotting purpose
    e_bnd(ik,:) = e_k(1:nbnd_plot)
  ENDDO

!  STOP 

  distk(1) = 0.0_DP
  DO ik=2,nks
    !distk(ik) = distk(ik-1) + sqrt( sum( (xk(:,ik) - xk(:,ik-1) )**2 ) )
    distk(ik) = distk(ik-1) + sqrt( (xk(1,ik)-xk(1,ik-1))**2 + &
      (xk(2,ik)-xk(2,ik-1))**2 + (xk(3,ik)-xk(3,ik-1))**2 )
  ENDDO

  IF( ionode ) THEN
    DO ib=1,nbnd_plot
      DO ik=1,nks
        WRITE(111,'(1x,2F18.10)') distk(ik), e_bnd(ik,ib)
      ENDDO
      WRITE(111,*)
    ENDDO
  ENDIF

  DEALLOCATE( H_k, psi_k, e_k, e_bnd, distk )
END SUBROUTINE


!=-----------------------------------------------------------------------------
PROGRAM t_gvect
!=-----------------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE io_global, ONLY: stdout
  USE mp_global, ONLY : mp_startup, mp_global_end
  IMPLICIT NONE

  CALL mp_startup()

  CALL setup_structure()

  CALL setup_fft()

  CALL setup_symmetry()

  CALL setup_kpoints()

  CALL setup_gvect()

  CALL test_hamiltonian1()

  CALL mp_global_end()

END PROGRAM




