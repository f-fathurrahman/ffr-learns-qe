!------------------------------------------------------------------------------
SUBROUTINE setup_structure_H()
!------------------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE constants, ONLY : PI, BOHR_RADIUS_ANGS
  USE io_global, ONLY : stdout, ionode
  USE cell_base, ONLY : cell_base_init
  USE cell_base, ONLY : at, alat
  USE ions_base, ONLY : if_pos, ityp, tau, extfor, nsp, &
                        atm, nat, tau_format
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
  a_ = 2_DP*PI*BOHR_RADIUS_ANGS   ! in angstrom
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
