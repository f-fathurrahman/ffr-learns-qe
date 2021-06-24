MODULE my_symme
  
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : s, sname, ft, nrot, nsym, t_rev, time_reversal, &
                         irt, invs, invsym

  !
  LOGICAL :: no_rho_sym=.false.
  INTEGER :: ngs              ! number of symmetry-related G-vector shells
  
  TYPE shell_type
     INTEGER, POINTER :: vect(:)
  END TYPE shell_type
  
  ! shell contains a list of symmetry-related G-vectors for each shell
  TYPE(shell_type), ALLOCATABLE :: shell(:)
  
  ! Arrays used for parallel symmetrization
  INTEGER, ALLOCATABLE :: sendcnt(:), recvcnt(:), sdispls(:), rdispls(:)

CONTAINS 

  SUBROUTINE my_sym_rho_init_shells ( ngm_, g_ )
    !-----------------------------------------------------------------------
    !
    !  Initialize G-vector shells needed for symmetrization
    ! 
    USE constants, ONLY : eps8
    USE mp_bands,  ONLY : nproc_bgrp
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ngm_
    REAL(DP), INTENT(IN) :: g_(3,ngm_)
    !
    LOGICAL, ALLOCATABLE :: done(:)
    INTEGER, ALLOCATABLE :: n(:,:), igsort(:)
    REAL(DP), ALLOCATABLE :: g2sort_g(:)
    INTEGER :: i,j,is,ig, iig, jg, ng, sn(3), gshell(3,48)
    LOGICAL :: found
    !
    ngs = 0
    ! shell should be allocated to the number of symmetry shells
    ! since this is unknown, we use the number of all G-vectors
    ALLOCATE ( shell(ngm_) )
    ALLOCATE ( done(ngm_), n(3,ngm_) )
    ALLOCATE ( igsort (ngm_))
    DO ig=1,ngm_
       !
       done(ig) = .false.
       ! G-vectors are stored as integer indices in crystallographic axis:
       !    G = n(1)*at(1) + n(2)*at(2) + n(3)*at(3)
       n(:,ig) = nint ( at(1,:)*g_(1,ig) + at(2,:)*g_(2,ig) + at(3,:)*g_(3,ig) )
       !
       NULLIFY(shell(ig)%vect)
       !
    END DO
!
!   The following algorithm can become very slow if ngm_ is large and
!   g vectors are not ordered in increasing order. This happens 
!   in the parallel case.
!
    IF (nproc_bgrp > 1 .AND. ngm_ > 20000) THEN
       ALLOCATE ( g2sort_g(ngm_))
       g2sort_g(:)=g_(1,:)*g_(1,:)+g_(2,:)*g_(2,:)+g_(3,:)*g_(3,:)
       igsort(1) = 0
       CALL hpsort_eps( ngm_, g2sort_g, igsort, eps8 )
       DEALLOCATE( g2sort_g)
    ELSE
       DO ig=1,ngm_
          igsort(ig)=ig
       ENDDO
    ENDIF
    !
    DO iig=1,ngm_
       !
       ig=igsort(iig)
       IF ( done(ig) ) CYCLE
       !
       ! we start a new shell of symmetry-equivalent G-vectors
       ngs = ngs+1
       ! ng: counter on G-vectors in this shell
       ng  = 0
       DO is=1,nsym
          ! integer indices for rotated G-vector
          sn(:)=s(:,1,is)*n(1,ig)+s(:,2,is)*n(2,ig)+s(:,3,is)*n(3,ig)
          found = .false.
          ! check if this rotated G-vector is equivalent to any other
          ! vector already present in this shell
shelloop: DO i=1,ng
             found = ( sn(1)==gshell(1,i) .and. &
                       sn(2)==gshell(2,i) .and. &
                       sn(3)==gshell(3,i) )
             if (found) exit shelloop
          END DO shelloop
          IF ( .not. found ) THEN
             ! add rotated G-vector to this shell
             ng = ng + 1
             IF (ng > 48) CALL errore('sym_rho_init_shell','internal error',48)
             gshell(:,ng) = sn(:)
          END IF
       END DO
       ! there are ng vectors gshell in shell ngs
       ! now we have to locate them in the list of G-vectors
       ALLOCATE ( shell(ngs)%vect(ng))
       DO i=1,ng
gloop:    DO jg=iig,ngm_
             j=igsort(jg)
             IF (done(j)) CYCLE gloop
                found = ( gshell(1,i)==n(1,j) .and. &
                          gshell(2,i)==n(2,j) .and. &
                          gshell(3,i)==n(3,j) )
             IF ( found ) THEN
                done(j)=.true.
                shell(ngs)%vect(i) = j
                EXIT gloop
             END IF
          END DO gloop
          IF (.not. found) CALL errore('sym_rho_init_shell','lone vector',i)
       END DO
       !
    END DO
    DEALLOCATE ( n, done ) 
    DEALLOCATE( igsort)

  END SUBROUTINE


END MODULE