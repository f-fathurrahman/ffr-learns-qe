!------------------------------------------------------------------------------
SUBROUTINE setup_structure_Si8()
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
