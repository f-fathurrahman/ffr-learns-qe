INCLUDE 'prepare_all.f90'


PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  CALL test_PAW_rad2lm()
END PROGRAM

!---------------------------
subroutine test_PAW_rad2lm()
!---------------------------
  use paw_onecenter, only: PAW_rad2lm
  USE atom, ONLY : g => rgrid
  USE ions_base, ONLY : ityp
  USE uspp_param, ONLY : upf
  USE paw_variables, ONLY : paw_info, rad
  !
  IMPLICIT NONE
  !
  TYPE(paw_info) :: i
  INTEGER :: ia
  INTEGER :: l2
  REAL(8), ALLOCATABLE :: gc_rad(:,:,:)
  REAL(8), ALLOCATABLE :: gc_lm(:,:,:)

  ! Choose atom index and which partial waves to be used (AE or PS)
  ia = 1

  i%a = ia   ! atom's index
  i%t = ityp(ia) ! type of atom ia
  i%m = g(i%t)%mesh ! radial mesh size for atom i%t
  i%b = upf(i%t)%nbeta  ! number of beta functions for i%t
  i%l = upf(i%t)%lmax_rho + 1  ! max ang.mom. in augmentation for ia
  l2  = i%l**2

  write(*,*) 'Nrmesh = ', i%m

  ALLOCATE( gc_rad(i%m, rad(i%t)%nx, 1) )
  ALLOCATE( gc_lm(i%m, i%l**2, 1) )


  gc_rad(:,:,:) = 1.1d0

  CALL PAW_rad2lm(i, gc_rad, gc_lm, i%l, 1)

  write(*,*) 'lmax_loc = ', i%l
  write(*,*) 'sum wwylm = ', sum(rad(i%t)%wwylm)
  write(*,*) 'sum gc_rad = ', sum(gc_rad)
  write(*,*) 'sum gc_lm = ', sum(gc_lm)

  deallocate(gc_rad, gc_lm)

end subroutine


!--------------------------------------------------------------------------------
SUBROUTINE my_PAW_rad2lm( i, F_rad, F_lm, lmax_loc, nspin )
!------------------------------------------------------------------------------
  use kinds, only : dp
  use paw_variables, only : paw_info, rad
  implicit none
  !! Computes:
  !! \[ F_{lm}(r) = \int d \Omega\ F(r,\text{th},\text{ph})\ Y_{lm}(\text{th},
  !! \text{ph}) \]
  !
  TYPE(paw_info), INTENT(IN) :: i
  !! atom's minimal info
  INTEGER, INTENT(IN) :: nspin
  !! spin configuration label
  INTEGER,  INTENT(IN) :: lmax_loc
  !! In some cases I have to keep higher angular components
  !! than the default ones (=lmaxq =the ones present in rho)
  REAL(DP), INTENT(OUT):: F_lm(i%m, lmax_loc**2, nspin)
  !! lm component of F up to lmax_loc
  REAL(DP), INTENT(IN) :: F_rad(i%m, rad(i%t)%nx, nspin)
  !! radial samples of F
  !
  ! ... local variables
  !
  INTEGER :: ix    ! counter for integration
  INTEGER :: lm    ! counter for angmom
  INTEGER :: ispin ! counter for spin
  INTEGER :: j

  !write(*,*) 'shape F_lm = ', shape(F_lm)
  !write(*,*) 'shape F_rad = ', shape(F_rad)
  !write(*,*) 'sum F_lm = ', sum(F_lm)

  write(*,*)
  write(*,*) '*** Using my_PAW_rad2lm'
  write(*,*)

  DO ispin = 1, nspin
    DO lm = 1, lmax_loc**2
      F_lm(:,lm,ispin) = 0.d0
      DO ix = 1, rad(i%t)%nx
        DO j = 1, i%m
          F_lm(j, lm, ispin) = F_lm(j, lm, ispin) + F_rad(j,ix,ispin) * rad(i%t)%wwylm(ix,lm)
          !write(111,'(1x,3I5,2F18.10)') lm, ix, j, rad(i%t)%wwylm(ix,lm), F_lm(j,lm,ispin)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  !do ispin = 1, nspin
  !  do lm = 1, lmax_loc**2
  !    do j = 1, i%m
  !      write(222,'(1x,2I5,F18.10)') j, lm, F_lm(j,lm,ispin)
  !    enddo
  !  enddo
  !enddo

  return
END SUBROUTINE