INCLUDE 'prepare_all.f90'

PROGRAM main
  !
  USE kinds, ONLY: DP
  USE gvect, ONLY: ngm
  !
  IMPLICIT NONE 

  COMPLEX(DP), ALLOCATABLE :: rhocg(:,:)
  INTEGER :: nspina

  CALL prepare_all()

  nspina = 1
  ALLOCATE( rhocg(ngm,nspina) )

  CALL my_atomic_rho_g(rhocg, nspina)

  DEALLOCATE( rhocg )

END PROGRAM 



!-----------------------------------------------------------------------
SUBROUTINE my_atomic_rho_g( rhocg, nspina )
!-----------------------------------------------------------------------
  !! Compute superposition of atomic charges in reciprocal space.
  !
  !! Three cases:
  !
  !! * if \(\text{nspina}=1\) the total atomic charge density is calculated;
  !! * if \(\text{nspina}=2\) collinear case. The total density is calculated
  !!               in the first component and the magnetization in 
  !!               the second;
  !! * if \(\text{nspina}=4\) noncollinear case. Total density in the first
  !!               component and magnetization vector in the
  !!               other three.
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : eps8
  USE atom,                 ONLY : rgrid, msh
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : tpiba, omega
  USE gvect,                ONLY : ngm, ngl, gstart, gl, igtongl
  USE lsda_mod,             ONLY : starting_magnetization
  USE vlocal,               ONLY : starting_charge, strf
  USE noncollin_module,     ONLY : angle1, angle2
  USE uspp_param,           ONLY : upf
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nspina
  !! number of spin components to be calculated. It may differ from
  !! nspin because in some cases the total charge only is needed, 
  !! even in a LSDA calculation.
  COMPLEX(DP), INTENT(OUT) :: rhocg(ngm,nspina)
  !! contains G-space components of the superposition of atomic charges
  !! contained in the array upf%rho_at (read from pseudopotential files).
  !
  ! ... local variables
  !
  REAL(DP) :: rhoscale, gx
  REAL(DP), ALLOCATABLE :: rhocgnt(:), aux(:)
  REAL(DP) :: angular(nspina)
  INTEGER :: ir, is, ig, igl, nt, ndm
  !
  ! allocate work space 
  !
  ndm = MAXVAL ( msh(1:ntyp) )
  ALLOCATE (rhocgnt( ngl))
  ALLOCATE (aux(ndm))
  !
  DO nt = 1, ntyp
    !
    ! Here we compute the G=0 term
    !
    IF( gstart == 2 ) THEN 
      DO ir = 1, msh (nt)
        aux (ir) = upf(nt)%rho_at (ir)
      ENDDO
      CALL simpson( msh(nt), aux, rgrid(nt)%rab, rhocgnt(1) )
    ENDIF
    !
    ! Here we compute the G<>0 term
    !
    DO igl = gstart, ngl
      gx = sqrt (gl (igl) ) * tpiba
      DO ir = 1, msh (nt)
        IF( rgrid(nt)%r(ir) < eps8 ) THEN 
          aux(ir) = upf(nt)%rho_at(ir)
        ELSE 
          aux(ir) = upf(nt)%rho_at(ir) * sin(gx*rgrid(nt)%r(ir)) / (rgrid(nt)%r(ir)*gx)
        ENDIF
      ENDDO
      CALL simpson( msh(nt), aux, rgrid(nt)%rab, rhocgnt(igl) )
    ENDDO
    !
    ! we compute the 3D atomic charge in reciprocal space
    !
    IF( upf(nt)%zp > eps8 ) THEN
      rhoscale = MAX(0.0_dp, upf(nt)%zp - starting_charge(nt)) / upf(nt)%zp
      WRITE(*,*) 'Pass here 104 in my_atomic_rho'
    ELSE
      rhoscale = 1.0_dp
    ENDIF
    WRITE(*,*) 'rhoscale = ', rhoscale
    WRITE(*,*) 'starting_charge(nt) = ', starting_charge(nt)
    !
    !
    DO ig = 1, ngm
      rhocg(ig,1) = rhocg(ig,1) + strf(ig,nt) * rhoscale * rhocgnt(igtongl(ig)) / omega
    ENDDO
    WRITE(*,*) 'sum rhocg = ', sum(rhocg(:,1))
    !
    IF( nspina >= 2 ) THEN
      !
      angular(1) = 1._dp
      IF( nspina == 4 ) THEN
        angular(1) = sin(angle1(nt))*cos(angle2(nt))
        angular(2) = sin(angle1(nt))*sin(angle2(nt))
        angular(3) = cos(angle1(nt))
      ENDIF
      !
      DO is = 2, nspina
        DO ig = 1, ngm
          rhocg(ig,is) = rhocg(ig,is) + &
                         starting_magnetization(nt) * angular(is-1) * &
                         strf(ig,nt) * rhoscale * rhocgnt(igtongl(ig)) / omega
        ENDDO
      ENDDO
      !
    ENDIF 
    ! must complete the computation of rhocg before updating rhocgnt
    ! for the next type
  ENDDO

  DEALLOCATE(aux)
  DEALLOCATE(rhocgnt)

END SUBROUTINE 

