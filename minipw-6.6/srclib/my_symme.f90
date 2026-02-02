
!--------------------------------------------------------------------------
MODULE my_symme
!------------------------------------------------------------------------
  !! This module contains routines used for symmetrization.
  !
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : s, sname, ft, nrot, nsym, t_rev, time_reversal, &
                         irt, invs, invsym
  !
  SAVE
  !PRIVATE
  !
  ! General-purpose symmetrizaton routines
  !
  !PUBLIC ::  symscalar, symvector, symtensor, symmatrix, symv, &
  !           symtensor3, symmatrix3, crys_to_cart, cart_to_crys
  ! For symmetrization in reciprocal space (all variables are private)
  !
  !PUBLIC :: sym_rho_init, sym_rho, sym_rho_deallocate
  !
  LOGICAL :: no_rho_sym=.true.      ! do not perform symetrization of charge density
  INTEGER :: ngs              ! number of symmetry-related G-vector shells

  TYPE shell_type
     INTEGER, POINTER :: vect(:)
  END TYPE shell_type

  ! shell contains a list of symmetry-related G-vectors for each shell
  TYPE(shell_type), ALLOCATABLE :: shell(:)

  ! Arrays used for parallel symmetrization
  INTEGER, ALLOCATABLE :: sendcnt(:), recvcnt(:), sdispls(:), rdispls(:)
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  LOGICAL FUNCTION rho_sym_needed( )
  !-----------------------------------------------------------------------
    !! TRUE if rho symmetrization is needed.
    !
    rho_sym_needed = .NOT. no_rho_sym
    !
  END FUNCTION rho_sym_needed

  !-------------------------------------------------------------------------
  SUBROUTINE symscalar( nat, scalar )
  !-----------------------------------------------------------------------
     !! Symmetrize a scalar function \(f(na)\), where na is the atom index.
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     !! number of atoms
     REAL(DP), INTENT(INOUT) :: scalar(nat)
     !! function to symmetrize
     !
     INTEGER :: isym
     REAL(DP), ALLOCATABLE :: work (:)
     !
     IF (nsym == 1) RETURN
     !
     ALLOCATE (work(nat))
     work(:) = 0.0_dp
     DO isym = 1, nsym
        work (:) = work (:) +  scalar(irt(isym,:))
     END DO
     scalar(:) = work(:) / DBLE(nsym)
     DEALLOCATE (work)
     !
   END SUBROUTINE symscalar
   !
   !--------------------------------------------------------------------------
   SUBROUTINE symvector( nat, vect )
     !-----------------------------------------------------------------------
     !! Symmetrize a function \(f(i,na)\) (e.g. the forces in cartesian axis),
     !! where \(i\) is the cartesian component, \(na\) the atom index.
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     !! number of atoms
     REAL(DP), INTENT(INOUT) :: vect(3,nat)
     !! vector function to symmetrize
     !
     ! ... local variables
     !
     INTEGER :: na, isym, nar
     REAL(DP), ALLOCATABLE :: work(:,:)
     !
     IF (nsym == 1) RETURN
     !
     ALLOCATE (work(3,nat))
     !
     ! bring vector to crystal axis
     !
     DO na = 1, nat
        work(:,na) = vect(1,na)*at(1,:) + &
                     vect(2,na)*at(2,:) + &
                     vect(3,na)*at(3,:)
     END DO
     !
     ! symmetrize in crystal axis
     !
     vect (:,:) = 0.0_dp
     DO na = 1, nat
        DO isym = 1, nsym
           nar = irt (isym, na)
           vect (:, na) = vect (:, na) + &
                          s (:, 1, isym) * work (1, nar) + &
                          s (:, 2, isym) * work (2, nar) + &
                          s (:, 3, isym) * work (3, nar)
        END DO
     END DO
     work (:,:) = vect (:,:) / DBLE(nsym)
     !
     ! bring vector back to cartesian axis
     !
     DO na = 1, nat
        vect(:,na) = work(1,na)*bg(:,1) + &
                     work(2,na)*bg(:,2) + &
                     work(3,na)*bg(:,3)
     END DO
     !
     DEALLOCATE (work)
     !
   END SUBROUTINE symvector
   !
   !--------------------------------------------------------------------------
   SUBROUTINE symtensor( nat, tens )
     !-----------------------------------------------------------------------
     !! Symmetrize a function \(f(i,j,na)\) (e.g. the effective charges in
     !! cartesian axis), where \(i,j\) are the cartesian components and \(na\)
     !! is the atom index.
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     !! number of atoms
     REAL(DP), INTENT(INOUT) :: tens(3,3,nat)
     !! tensor function to symmetrize
     !
     ! ... local variables
     !
     INTEGER :: na, isym, nar, i,j,k,l
     REAL(DP), ALLOCATABLE :: work (:,:,:)
     !
     IF (nsym == 1) RETURN
     !
     ! bring tensor to crystal axis
     !
     DO na=1,nat
        CALL cart_to_crys ( tens (:,:,na) )
     END DO
     !
     ! symmetrize in crystal axis
     !
     ALLOCATE (work(3,3,nat))
     work (:,:,:) = 0.0_dp
     DO na = 1, nat
        DO isym = 1, nsym
           nar = irt (isym, na)
           DO i = 1, 3
              DO j = 1, 3
                 DO k = 1, 3
                    DO l = 1, 3
                       work (i,j,na) = work (i,j,na) + &
                          s (i,k,isym) * s (j,l,isym) * tens (k,l,nar)
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO
     tens (:,:,:) = work (:,:,:) / DBLE(nsym)
     DEALLOCATE (work)
     !
     ! bring tensor back to cartesian axis
     !
     DO na=1,nat
        CALL crys_to_cart ( tens (:,:,na) )
     END DO
     !
     !
   END SUBROUTINE symtensor
   !
   !-----------------------------------------------------------------------
   SUBROUTINE symv( vect )
     !--------------------------------------------------------------------
     !! Symmetrize a vector \(f(i)\), i=cartesian components
     !! The vector is supposed to be axial: inversion does not change it.
     !! Time reversal changes its sign. Note that only groups compatible with
     !! a finite magnetization give a nonzero output vector.
     !
     IMPLICIT NONE
     !
     REAL (DP), INTENT(inout) :: vect(3)
     !! the vector to rotate
     !
     ! ... local variables
     !
     INTEGER :: isym
     REAL(DP) :: work(3), segno
     !
     IF (nsym == 1) RETURN
     !
     ! bring vector to crystal axis
     !
     work(:) = vect(1)*at(1,:) + vect(2)*at(2,:) + vect(3)*at(3,:)
     vect = work
     work=0.0_DP
     do isym = 1, nsym
        segno=1.0_DP
        IF (sname(isym)(1:3)=='inv') segno=-1.0_DP
        IF (t_rev(isym)==1) segno=-1.0_DP*segno
        work (:) = work (:) + segno * ( &
                       s (:, 1, isym) * vect (1) + &
                       s (:, 2, isym) * vect (2) + &
                       s (:, 3, isym) * vect (3) )
     enddo
     work=work/nsym
   !
   !  And back in cartesian coordinates.
   !
   vect(:) = work(1) * bg(:,1) + work(2) * bg(:,2) + work(3) * bg(:,3)
   !
   end subroutine symv
   !
   !-------------------------------------------------------------------------
   SUBROUTINE symmatrix( matr )
     !-----------------------------------------------------------------------
     !! Symmetrize a function \(f(i,j)\) (e.g. stress, dielectric tensor in
     !! cartesian axis), where \(i,j\) are the cartesian components.
     !
     IMPLICIT NONE
     !
     REAL(DP), INTENT(INOUT) :: matr(3,3)
     !! the function \(f(i,j)\) to symmetrize
     !
     ! ... local variables
     !
     INTEGER :: isym, i,j,k,l
     REAL(DP) :: work (3,3)
     !
     IF (nsym == 1) RETURN
     !
     ! bring matrix to crystal axis
     !
     CALL cart_to_crys ( matr )
     !
     ! symmetrize in crystal axis
     !
     work (:,:) = 0.0_dp
     DO isym = 1, nsym
        DO i = 1, 3
           DO j = 1, 3
              DO k = 1, 3
                 DO l = 1, 3
                    work (i,j) = work (i,j) + &
                       s (i,k,isym) * s (j,l,isym) * matr (k,l)
                 END DO
              END DO
           END DO
        END DO
     END DO
     matr (:,:) = work (:,:) / DBLE(nsym)
     !
     ! bring matrix back to cartesian axis
     !
     CALL crys_to_cart ( matr )
     !
   END SUBROUTINE symmatrix
   !
   !------------------------------------------------------------------------
   SUBROUTINE symmatrix3( mat3 )
     !-----------------------------------------------------------------------
     !! Symmetrize a function \(f(i,j,k)\) (e.g. nonlinear susceptibility),
     !! where \(i,j,k\) are the cartesian components.
     !! BEWARE: input in crystal axis, output in cartesian axis.
     !
     IMPLICIT NONE
     !
     REAL(DP), INTENT(INOUT) :: mat3(3,3,3)
     !! function f(i,j,k) to symmetrize
     !
     ! ... local variables
     !
     INTEGER :: isym, i,j,k,l,m,n
     REAL(DP) :: work(3,3,3)
     !
     IF (nsym > 1) THEN
        !
        work (:,:,:) = 0.0_dp
        DO isym = 1, nsym
           DO i = 1, 3
              DO j = 1, 3
                 DO k = 1, 3
                    DO l = 1, 3
                       DO m = 1, 3
                          DO n = 1, 3
                             work (i, j, k) = work (i, j, k) + &
                                s (i, l, isym) * s (j, m, isym) * &
                                s (k, n, isym) * mat3 (l, m, n)
                          END DO
                       END DO
                    END DO
                 END DO
              END DO
           END DO
        END DO
        mat3 = work/ DBLE(nsym)
        !
     END IF
     !
     ! Bring to cartesian axis
     !
     CALL crys_to_cart_mat3 ( mat3 )
     !
   END SUBROUTINE symmatrix3
   !
   !-------------------------------------------------------------------------
   SUBROUTINE symtensor3( nat, tens3 )
     !-----------------------------------------------------------------------
     !! Symmetrize a function \(f(i,j,k, na)\) (e.g. the Raman tensor), where
     !! \(i,j,k\) are the cartesian axes, \(na\) is the atom index.
     !! BEWARE: input in crystal axis, output in cartesian axis
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     !! number of atoms
     REAL(DP), INTENT(INOUT) :: tens3(3,3,3,nat)
     !! the function f(i,j,k, na) to symmetrize
     !
     ! ... local variables
     !
     INTEGER :: na, isym, nar, i,j,k,l,n,m
     REAL(DP), ALLOCATABLE :: work (:,:,:,:)
     !
     IF (nsym > 1) THEN
        !
        ! symmetrize in crystal axis
        !
        ALLOCATE (work(3,3,3,nat))
        work (:,:,:,:) = 0.0_dp
        DO na = 1, nat
           DO isym = 1, nsym
              nar = irt (isym, na)
              DO i = 1, 3
                 DO j = 1, 3
                    DO k = 1, 3
                       DO l = 1, 3
                          DO m =1, 3
                             DO n =1, 3
                                work (i, j, k, na) = work (i, j, k, na) + &
                                  s (i, l, isym) * s (j, m, isym) *    &
                                  s (k, n, isym) * tens3 (l, m, n, nar)
                             END DO
                          END DO
                       END DO
                    END DO
                 END DO
              END DO
           END DO
        END DO
        tens3 (:,:,:,:) =   work(:,:,:,:) / DBLE (nsym)
        DEALLOCATE (work)
        !
     END IF
     !
     ! Bring to cartesian axis
     !
     DO na = 1, nat
        CALL crys_to_cart_mat3 ( tens3(:,:,:,na) )
     END DO
     !
   END SUBROUTINE symtensor3
   !
   ! Routines for crystal to cartesian axis conversion
   !
   !INTERFACE cart_to_crys
   !  MODULE PROCEDURE cart_to_crys_mat, cart_to_crys_mat3
   !END INTERFACE
   !INTERFACE crys_to_cart
   !  MODULE PROCEDURE crys_to_cart
   !END INTERFACE
   !
   !-------------------------------------------------------------------------
   SUBROUTINE cart_to_crys( matr )
     !-----------------------------------------------------------------------
     !! Cartesian to crystal axis conversion.
     !
     IMPLICIT NONE
     !
     REAL(DP), INTENT(INOUT) :: matr(3,3)
     !! Axis conversion matrix
     !
     ! ... local variables
     !
     REAL(DP) :: work(3,3)
     INTEGER :: i,j,k,l
     !
     work(:,:) = 0.0_dp
     DO i = 1, 3
        DO j = 1, 3
           DO k = 1, 3
              DO l = 1, 3
                 work(i,j) = work(i,j) + matr(k,l) * at(k,i) * at(l,j)
              END DO
           END DO
        END DO
     END DO
     !
     matr(:,:) = work(:,:)
     !
   END SUBROUTINE cart_to_crys
   !
   !-------------------------------------------------------------------------
   SUBROUTINE crys_to_cart( matr )
     !-----------------------------------------------------------------------
     !! Crystal to cartesian axis conversion.
     !
     IMPLICIT NONE
     !
     REAL(DP), INTENT(INOUT) :: matr(3,3)
     !! Axis conversion matrix
     !
     ! ... local variables
     !
     REAL(DP) :: work(3,3)
     INTEGER :: i,j,k,l
     !
     work(:,:) = 0.0_dp
     DO i = 1, 3
        DO j = 1, 3
           DO k = 1, 3
              DO l = 1, 3
                 work(i,j) = work(i,j) + &
                             matr(k,l) * bg(i,k) * bg(j,l)
              END DO
           END DO
        END DO
     END DO
     matr(:,:) = work(:,:)
     !
   END SUBROUTINE crys_to_cart
   !
   !------------------------------------------------------------------------
   SUBROUTINE crys_to_cart_mat3( mat3 )
     !-----------------------------------------------------------------------
     !! Crystal to cartesian axis conversion for \(f(i,j,k)\) objects.
     !
     IMPLICIT NONE
     !
     REAL(DP), INTENT(INOUT) :: mat3(3,3,3)
     !! Axis conversion tensor
     !
     REAL(DP) :: work(3,3,3)
     INTEGER :: i,j,k,l,m,n
     !
     ! ... local variables
     !
     work(:,:,:) = 0.0_dp
     DO i = 1, 3
        DO j = 1, 3
           DO k = 1, 3
              DO l = 1, 3
                 DO m = 1, 3
                    DO n = 1, 3
                       work (i, j, k) = work (i, j, k) +  &
                          mat3 (l, m, n) * bg (i, l) * bg (j, m) * bg (k, n)
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO
     mat3(:,:,:) = work (:,:,:)
     !
   END SUBROUTINE crys_to_cart_mat3
   !
   ! G-space symmetrization
   !
   !------------------------------------------------------------------------
   SUBROUTINE sym_rho_init( gamma_only )
    !-----------------------------------------------------------------------
    !! Initialize arrays needed for symmetrization in reciprocal space.
    !
    USE gvect, ONLY : ngm, g
    !
    LOGICAL, INTENT(IN) :: gamma_only
    !
    no_rho_sym = gamma_only .OR. (nsym==1)
    IF (no_rho_sym) RETURN
    CALL sym_rho_init_shells( ngm, g )
    !
  END SUBROUTINE sym_rho_init

  !-----------------------------------------------------------------------
  SUBROUTINE sym_rho_init_shells( ngm_, g_ )
  !-----------------------------------------------------------------------
    !! Initialize G-vector shells needed for symmetrization.
    !
    USE constants, ONLY : eps8
    USE mp_bands,  ONLY : nproc_bgrp
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ngm_
    !! number of g-points
    REAL(DP), INTENT(IN) :: g_(3,ngm_)
    !! G-vectors
    !
    ! ... local variables
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
    ALLOCATE( shell(ngm_) )
    ALLOCATE( done(ngm_), n(3,ngm_) )
    ALLOCATE( igsort (ngm_) )
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

  END SUBROUTINE sym_rho_init_shells
  !
  !-----------------------------------------------------------------------
  SUBROUTINE sym_rho (nspin, rhog)
    !-----------------------------------------------------------------------
    !! Symmetrize the charge density rho in reciprocal space.
    !
    !! Distributed parallel algorithm: collects entire shells of G-vectors
    !! and corresponding rho(G), calls sym_rho_serial to perform the
    !! symmetrization, re-distributed rho(G) into original ordering.
    !
    USE constants,            ONLY : eps8, eps6
    USE gvect,                ONLY : ngm, g
    USE parallel_include
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nspin
    !! nspin=1,2,4 \(\rightarrow\) unpolarized, LSDA,
    !! non-colinear magnetism
    COMPLEX(DP), INTENT(INOUT) :: rhog(ngm,nspin)
    !! components of rho: rhog(ig) = rho(G(:,ig)).
    !! Unsymmetrized on input, symmetrized on output
    !
    !
    IF ( no_rho_sym) RETURN
    CALL sym_rho_serial ( ngm, g, nspin, rhog )
    RETURN
  END SUBROUTINE sym_rho

  !---------------------------------------------------------------------
  SUBROUTINE sym_rho_deallocate( )
    !-------------------------------------------------------------------
    !! Deallocates symmetrization objects.
    !
    IMPLICIT NONE
    !
    INTEGER :: i
    !
    IF ( ALLOCATED (rdispls) ) DEALLOCATE (rdispls)
    IF ( ALLOCATED (recvcnt) ) DEALLOCATE (recvcnt)
    IF ( ALLOCATED (sdispls) ) DEALLOCATE (sdispls)
    IF ( ALLOCATED (sendcnt) ) DEALLOCATE (sendcnt)
    IF ( ALLOCATED (shell) ) THEN
       DO i=1,SIZE(shell)
          IF ( ASSOCIATED(shell(i)%vect) ) DEALLOCATE (shell(i)%vect)
       END DO
       DEALLOCATE (shell)
    END IF
    !
  END SUBROUTINE sym_rho_deallocate
  !
END MODULE my_symme
