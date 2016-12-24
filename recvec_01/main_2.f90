! ffr: 23 Juni 2016
!
! An example of how to work with gvect: calculation of empty-lattice


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
  !
  ! Atomic positions
  nsp = 1
  nat = 1
  !
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
  USE gvecs, ONLY: ecuts, gcutms, dual, doublegrid
  USE fft_base, ONLY: dfftp, dffts
  USE fft_types, ONLY: fft_type_allocate, fft_type_init
  USE mp_bands, ONLY: intra_bgrp_comm
  USE paw_variables, ONLY: okpaw
  USE uspp, ONLY: okvan
  IMPLICIT NONE

  ! FIXME
  okvan = .FALSE.
  okpaw = .FALSE.

  ! FIXME: This should become arguments of this subroutine
  ecutwfc = 30_DP
  ecutrho = 120_DP

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
  CALL gshells( .false. )
END SUBROUTINE


!------------------------------------------------------------------------------
SUBROUTINE test_gvectors()
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE cell_base, ONLY : tpiba2
  USE gvecw, ONLY : ecutwfc, gcutw
  USE wvfct, ONLY : npwx, npw, current_k, g2kin
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
  USE io_global, ONLY: stdout
  USE mp_global, ONLY : mp_startup, mp_global_end
  IMPLICIT NONE

  CALL mp_startup()

  CALL setup_structure()

  CALL setup_fft()

  CALL setup_symmetry()

  CALL setup_kpoints()

  CALL setup_gvect()

  !CALL test_gvectors()
  CALL t_import_gvect()
  CALL t_import_gvecs()
  CALL t_import_gvecw()

  CALL mp_global_end()

END PROGRAM
