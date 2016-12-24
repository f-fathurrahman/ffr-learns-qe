! ffr: 5 Dec 2015
!
! An example of how to initialize G-vectors


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

  !
  ! ... Unit cell
  !
  ibrav_ = 1
  celldm_(:) = 0.0_DP
  a_ = 5.4307_DP
  b_ = 0.0_DP
  c_ = 0.0_DP
  cosab_ = 0.0_DP
  cosac_ = 0.0_DP
  cosbc_ = 0.0_DP
  trd_ht = .FALSE.
  rd_ht(:,:) = 0.0_DP
  cell_units_ = 'bohr'

  ! The variables defined in module cell_base are initialized here
  CALL cell_base_init( ibrav_, celldm_, a_, b_, c_, cosab_, cosac_, &
          & cosbc_, trd_ht, rd_ht, cell_units_ )

  !
  ! ... Atomic positions
  !
  nsp = 1
  nat = 8
  ! We need to allocate several allocatable arrays that depend on nsp and nat
  ! Note that, we don't need to allocate `atm` and amass` because
  ! they are static arrays with size `ntypmax`.
  ALLOCATE( ityp(nat) )
  ALLOCATE( tau(3,nat) )
  ALLOCATE( if_pos(3,nat) )
  ALLOCATE( extfor(3,nat) )

  ityp(:) = 1
  !
  atm(1) = 'Si'
  !
  tau_format = 'crystal'
  !
  tau = reshape( (/ &
         0.00_DP, 0.00_DP, 0.00_DP, &
         0.00_DP, 0.50_DP, 0.50_DP, &
         0.25_DP, 0.25_DP, 0.25_DP, &
         0.50_DP, 0.00_DP, 0.50_DP, &
         0.25_DP, 0.75_DP, 0.75_DP, &
         0.75_DP, 0.75_DP, 0.25_DP, &
         0.75_DP, 0.25_DP, 0.75_DP, &
         0.50_DP, 0.50_DP, 0.00_DP /), (/ 3, nat /) )
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



!------------------------------------------------------------------------------
SUBROUTINE setup_kpoints()
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE parameters, ONLY : npk
  USE cell_base, ONLY : at, bg
  USE klist, ONLY : nks, wk, xk, nkstot
  USE symm_base, ONLY : nrot, s, t_rev
  ! This is actually can be bypassed. However, we include it and set the
  ! relevant variables in anticipation that they will be accessed by
  ! other subroutines in PWSCF
  USE start_k, ONLY : nk1,nk2,nk3, k1,k2,k3
  IMPLICIT NONE
  !
  LOGICAL :: time_reversal, skip_equivalence
  INTEGER :: ik, ii

  time_reversal    = .TRUE.
  skip_equivalence = .FALSE. ! set to FALSE to reduce number of k-points

  ! k1,k2,k3 are the offset, the actual grid sampling are given by nk1,nk2,nk3
  k1 = 0; nk1 = 3
  k2 = 0; nk2 = 3
  k3 = 0; nk3 = 3

  CALL kpoint_grid( 1, time_reversal, skip_equivalence, s, t_rev, bg, &
    npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)

  nkstot = nks

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

  CALL mp_global_end()

END PROGRAM
