! ffr: 8 Dec 2015
!
! An example of how to initialize and work with pseudopotential data
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


! using readpp to read pseudopotential
SUBROUTINE setup_pseudo()
  USE io_global, ONLY : ionode, stdout
  USE io_files, ONLY : pseudo_dir, pseudo_dir_cur, psfile
  USE read_pseudo_mod, ONLY : readpp
  USE uspp_param, ONLY : upf
  IMPLICIT NONE
  CHARACTER(len=256) :: input_dft, dft_name
  LOGICAL :: printout

  input_dft = 'none'
  printout = .TRUE.
  pseudo_dir = '/home/efefer/pseudo_old/'
  !psfile(1)  = 'Si.pbe-rrkj.UPF'
  psfile(1)  = 'Si.pbe-n-van.UPF'
!  IF(ionode) THEN
!    WRITE(stdout,*) trim(pseudo_dir)
!    write(pseudo_dir_cur, len(psfile)
!  ENDIF
  CALL readpp( input_dft, printout )
  dft_name = upf(1)%dft
  IF( ionode ) THEN
    WRITE(stdout,*) 'dft_name = ', trim(dft_name)
    WRITE(stdout,*) 'zp = ', upf(1)%zp
    WRITE(stdout,*) 'typ, mesh = ', upf(1)%typ, upf(1)%mesh
  ENDIF


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
  CALL setup_pseudo()

!  CALL setup_fft()
!  CALL setup_symmetry()
!  CALL setup_kpoints()
!  CALL setup_gvect()
!  CALL test_hamiltonian1()

  CALL mp_global_end()

END PROGRAM




