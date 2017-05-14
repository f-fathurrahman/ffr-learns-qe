! ffr: 4 Dec 2015
!
! An example of how to initialize a crystalline or molecular structure using
! several routines available in Quantum Espresso.


PROGRAM t_struct
  USE kinds, ONLY: DP
  USE constants, ONLY : pi
  USE io_global, ONLY: stdout
  USE cell_base, ONLY: cell_base_init
  USE cell_base, ONLY: omega, at, bg, alat
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
  REAL(8) :: mag(3,3)
  !
  INTEGER :: ii, ia

  ibrav_ = 0
  celldm_(:) = 0.0_DP
  a_ = 0.0_DP
  b_ = 0.0_DP
  c_ = 0.0_DP
  cosab_ = 0.0_DP
  cosac_ = 0.0_DP
  cosbc_ = 0.0_DP
  trd_ht = .TRUE.
  rd_ht(1,:) = (/ 10.0_DP, 0.0_DP, 0.0_DP /)
  rd_ht(2,:) = (/ 8.0_DP, 5.0_DP, 0.0_DP /)
  rd_ht(3,:) = (/ 0.0_DP, 0.0_DP, 5.0_DP /)
  cell_units_ = 'bohr'

  ! The variables defined in module cell_base are initialized here
  CALL cell_base_init( ibrav_, celldm_, a_, b_, c_, cosab_, cosac_, &
          & cosbc_, trd_ht, rd_ht, cell_units_ )

  WRITE(stdout,*) 'alat = ', alat
  WRITE(stdout,*) 'Unit cell volume:', omega
  !
  WRITE(stdout,*) 'at = '
  DO ii=1,3
    WRITE(stdout,'(1x,3F18.10)') at(ii,:)
  ENDDO
  !
  WRITE(stdout,*) 'bg = '
  DO ii=1,3
    WRITE(stdout,'(1x,3F18.10)') bg(ii,:)
  ENDDO

  mag = matmul( at, transpose(bg) )
  WRITE(stdout,*) 'at*bg'' = '
  DO ii=1,3
    WRITE(stdout,'(1x,3F18.10)') mag(ii,:)
  ENDDO

  nsp = 2
  nat = 3
  ! We need to allocate several allocatable arrays that depend on nsp and nat
  ! Note that, we don't need to allocate `atm` and amass` because
  ! they are static arrays with size `ntypmax`.
  ALLOCATE( ityp(nat) )
  ALLOCATE( tau(3,nat) )
  ALLOCATE( if_pos(3,nat) )
  ALLOCATE( extfor(3,nat) )

  ityp(:) = (/ 1, 2, 2 /)
  !
  atm(1) = 'C'
  atm(2) = 'O'
  !
  tau_format = 'angstrom'
  !
  tau(:,1) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
  tau(:,2) = (/ 1.1_DP, 0.0_DP, 0.0_DP /)
  tau(:,3) = (/ -1.1_DP, 0.0_DP, 0.0_DP /)
  !
  if_pos(:,:) = 1
  extfor(:,:) = 0.0_DP
  !
  CALL convert_tau(tau_format, nat, tau)
  !
  !WRITE(stdout,*)
  WRITE(stdout, '(/,5x,"site n.     atom                  positions (alat units)")')
  WRITE(stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
             (ia, atm(ityp(ia)), ia, (tau(ii,ia), ii=1,3), ia=1,nat)

  WRITE(stdout, '(/,5x,"site n.     atom                  positions (bohr)")')
  WRITE(stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
             (ia, atm(ityp(ia)), ia, (tau(ii,ia)*alat, ii=1,3), ia=1,nat)

  CALL xsf_struct(alat, at, nat, tau, atm, ityp, 111)

END PROGRAM



