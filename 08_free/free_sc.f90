! ffr: 8 Dec 2015
!
! An example of how to work with gvect: calculation of empty-lattice
! Gamma-only


!------------------------------------------------------------------------------
SUBROUTINE setup_structure()
!------------------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE constants, ONLY : PI, BOHR_RADIUS_ANGS
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
  a_ = PI*BOHR_RADIUS_ANGS   ! in angstrom
  b_ = 0.0_DP ! we only need to specify a_ for ibrav_ = 1
  c_ = 0.0_DP
  cosab_ = 0.0_DP
  cosac_ = 0.0_DP
  cosbc_ = 0.0_DP
  trd_ht = .FALSE.
  rd_ht(:,:) = 0.0_DP  ! unit cell
  cell_units_ = 'bohr'  ! only take effect for trd_ht = .true.

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
  atm(1) = 'H'
  !
  tau_format = 'crystal'
  !
  tau(:,1) = (/ 0.5_DP, 0.5_DP, 0.5_DP /)
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
  USE kinds, ONLY: DP
  USE io_global, ONLY: stdout, ionode
  USE cell_base, ONLY: tpiba2, at, bg
  USE gvecw, ONLY: ecutwfc, gcutw
  USE gvect, ONLY: ecutrho, gcutm
  USE gvecs, ONLY: ecuts, gcutms, dual
  USE fft_base, ONLY: dfftp, dffts
  USE grid_subroutines, ONLY: realspace_grid_init
  IMPLICIT NONE
  
  ecutwfc = 40_DP
  ecutrho = 4_DP*ecutwfc
  dual    = ecutrho/ecutwfc
  gcutms  = 4_DP*ecutwfc/tpiba2
  gcutm   = dual*ecutwfc/tpiba2
  gcutw   = ecutwfc / tpiba2

  IF( ionode ) THEN
    WRITE(stdout,'(/,1x,A,F18.5)') 'dual = ', dual
    WRITE(stdout,'(1x,A,3F18.5)') 'gcutm, gcutms, gcutw = ', gcutm, gcutms, gcutw
  ENDIF

  CALL realspace_grid_init( dfftp, at, bg, gcutm )
  CALL realspace_grid_init( dffts, at, bg, gcutms )

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
  LOGICAL :: magnetic_sym
  REAL(DP) :: m_loc(3,nat)

  ! We ignore the magnetization
  magnetic_sym = .FALSE.
  m_loc(:,:)   = 0_DP

  CALL set_sym_bl()
  CALL find_sym( nat, tau, ityp, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
    magnetic_sym, m_loc )

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
  nks = 1
  nkstot = nks
  wk(1:nks) = 1.0_DP  ! it is not really used in band structure calculation
  xk(:,1:nks) = reshape( (/ 0.00000000, 0.00000000, 0.00000000 /), (/ 3,nks /) )

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
  
  gamma_only = .false.
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
  USE wvfct, ONLY : npwx, igk, npw, current_k, g2kin
  USE gvect, ONLY : ngm, g
  USE klist, ONLY : ngk, nks, xk
  IMPLICIT NONE
  ! Local
  INTEGER :: ik, ig, ib
  ! Full Hamiltonian and eigenvectors (column-wise)
  COMPLEX(8), ALLOCATABLE :: H_k(:,:), psi_k(:,:)
  !
  INTEGER :: nbnd_plot = 10
  REAL(8), ALLOCATABLE :: e_k(:), e_bnd(:,:), distk(:)

  ALLOCATE( ngk(nks) )
  CALL n_plane_waves( gcutw, nks, xk, g, ngm, npwx, ngk )

  ALLOCATE( igk(npwx) )
  ALLOCATE( g2kin(npwx) )

  ALLOCATE( H_k(npwx,npwx) )
  ALLOCATE( psi_k(npwx,npwx) )
  ALLOCATE( e_k(npwx) )
  ALLOCATE( e_bnd(nks,nbnd_plot) )
  ALLOCATE( distk(nks) )

  DO ik=1,nks
    ! generate Gk and sort them
    CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
    IF(ionode) WRITE(stdout,*) 'ik, npw, ngk(ik) = ', ik, npw, ngk(ik)
    ! Build the Hamiltonian
    H_k(:,:) = cmplx(0_DP, 0_DP)
    DO ig=1,npw
      H_k(ig,ig) = ( ( xk(1,ik) + g(1,igk(ig)) )**2 + &
                     ( xk(2,ik) + g(2,igk(ig)) )**2 + &
                     ( xk(3,ik) + g(3,igk(ig)) )**2 ) * tpiba2
    ENDDO
    CALL cdiagh( npw, H_k, npwx, e_k, psi_k )
    ! For plotting purpose
    e_bnd(ik,:) = e_k(1:nbnd_plot)
  ENDDO

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




