! efefer, 19 Dec 2015
!
! Experimental codes for implementation of CheFSI method for
! solving Kohn-Sham equations
!
! Refs:
! (1) Y. Zhou, et. al. J. Comput. Phys 219 (2006) pp. 172-184
! (2) Y. Zhou, et. al. Phys. Rev. E 74 (2006) pp. 066704
! (3) Y. Zhou, et. al. J. Comput. Phys 274 (2014) pp. 770-782
!


!------------------------------------------------------------------------------
SUBROUTINE cheFSI_scaled( ik, degree, lb, ub, scal )
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode, stdout
  USE wvfct, ONLY : npwx, npw, nbnd, et
  USE wavefunctions_module, ONLY : evc
  IMPLICIT NONE
  ! Arguments
  INTEGER, INTENT(in) :: ik
  INTEGER, INTENT(in) :: degree
  REAL(DP), INTENT(in) :: lb, ub, scal
  !
  COMPLEX(DP), ALLOCATABLE :: Y(:,:), Ynew(:,:)
  COMPLEX(DP), ALLOCATABLE :: T(:,:)
  REAL(DP), ALLOCATABLE :: eval_T(:) ! TODO: use et directly? FIXME: How about pools?
  !
  REAL(DP) :: e, c, sigma, sigmanew, tau
  !
  INTEGER :: i

  ALLOCATE( Y(npwx,nbnd) )
  ALLOCATE( Ynew(npwx,nbnd) )
  ALLOCATE( T(nbnd,nbnd) )
  ALLOCATE( eval_T(nbnd) )

  e = (ub - lb)/2.0_DP
  c = (ub + lb)/2.0_DP
  sigma = e/(c - scal)
  tau = 2.0_DP/sigma

  ! Y <-- H*evc
  CALL h_psi( npwx, npw, nbnd, evc, Y )
  ! Y <-- ( Y - c*evc)*sigma/e
  Y(1:npw,1:nbnd) = ( Y(1:npw,1:nbnd) - c*evc(1:npw,1:nbnd) ) * sigma/e

  DO i=2,degree
    sigmanew = 1.0_DP/( tau - sigma )

    ! Ynew <-- H*Y
    CALL h_psi( npwx, npw, nbnd, Y, Ynew )
    ! Ynew <-- ( Ynew - c*Y )*2*sigmanew/e - sigma*sigmanew*evc
    Ynew(1:npw,1:nbnd) = ( Ynew(1:npw,1:nbnd) - c*Y(1:npw,1:nbnd) )*2.0_DP*sigmanew/e - &
      sigma*sigmanew*evc(1:npw,1:nbnd)
    !
    evc(1:npw,1:nbnd) = Y(1:npw,1:nbnd)
    !
    Y(1:npw,1:nbnd) = Ynew(1:npw,1:nbnd)
    !
    sigma = sigmanew

  ENDDO

  !IF(ionode) WRITE(stdout,*) 'Before ortho:'
  !CALL test_orthonorm( Y, npwx, npw, nbnd )
  ! The filtered vector are in Y
  ! Orthonormalize Y
  CALL ortho_gram_schmidt( Y, npwx, npw, nbnd )
  !IF(ionode) WRITE(stdout,*) 'After ortho:'
  !CALL test_orthonorm( Y, npwx, npw, nbnd )

  !
  ! Compute Rayleigh quotient to obtain the eigenvectors
  !
  ! FIXME Use rotate_wfc instead?
  ! Ynew <-- H*Y
  CALL h_psi( npwx, npw, nbnd, Y, Ynew)
  ! T <-- <Y|H*Y>
  CALL zgemm( 'Conjg', 'No', nbnd, nbnd, npw, (1.d0,0.d0), Y,npwx, Ynew,npwx, &
    (0.d0,0.d0), T,nbnd )
  ! FIXME Symmetrize T?
  T = (T + conjg(transpose(T)))*0.5_DP  ! symmetrize
  CALL eig_zheevd( T,nbnd, eval_T,nbnd )
  ! evc <-- Y*T
  CALL zgemm('N','N', npw,nbnd,nbnd, (1.0_DP,0.0_DP), Y,npwx, T,nbnd, (0.d0,0.d0),evc,npwx)
  ! Copy eigenvalues
  DO i=1,nbnd
    et(i,ik) = eval_T(i)
  ENDDO
  !IF(ionode) WRITE(stdout,*) 'Orthonormalizing evc'
  !CALL ortho_gram_schmidt(evc, npwx, npw, nbnd)

  !
  DEALLOCATE( Y, Ynew, T, eval_T )

END SUBROUTINE





!------------------------------------------------------------------------------
SUBROUTINE cheFSI( ik, degree, lb, ub )
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode, stdout
  USE wvfct, ONLY : npwx, npw, nbnd, et
  USE wavefunctions_module, ONLY : evc
  IMPLICIT NONE
  ! Arguments
  INTEGER, INTENT(in) :: ik
  INTEGER, INTENT(in) :: degree
  REAL(DP), INTENT(in) :: lb, ub
  !
  COMPLEX(DP), ALLOCATABLE :: Y(:,:), Ynew(:,:)
  COMPLEX(DP), ALLOCATABLE :: T(:,:)
  REAL(DP), ALLOCATABLE :: eval_T(:) ! TODO: use et directly? FIXME: How about pools?
  !
  REAL(DP) :: e, c, sigma, sigma1, sigma2
  !
  INTEGER :: i

  ALLOCATE( Y(npwx,nbnd) )
  ALLOCATE( Ynew(npwx,nbnd) )
  ALLOCATE( T(nbnd,nbnd) )
  ALLOCATE( eval_T(nbnd) )

  Y(:,:)    = (0.0_DP, 0.0_DP)  ! TODO: need this?
  Ynew(:,:) = (0.0_DP, 0.0_DP)

  e = (ub - lb)/2.0_DP
  c = (ub + lb)/2.0_DP
  sigma = e/(lb - c)
  sigma1 = sigma

  ! Y <-- H*evc
  CALL h_psi( npwx, npw, nbnd, evc, Y )
  ! Y <-- ( Y - c*evc)*sigma1/e
  Y(1:npw,1:nbnd) = ( Y(1:npw,1:nbnd) - c*evc(1:npw,1:nbnd) ) * sigma1/e

  DO i=2,degree
    sigma2 = 1.0_DP/( 2.0_DP/sigma1 - sigma )

    ! Ynew <-- H*Y
    CALL h_psi( npwx, npw, nbnd, Y, Ynew )
    ! Ynew <-- 2*( Ynew - c*Y )*sigma2/e - sigma*sigma2*evc
    Ynew(1:npw,1:nbnd) = 2.0_DP* ( Ynew(1:npw,1:nbnd) - c*Y(1:npw,1:nbnd) )*sigma2/e - &
      sigma*sigma2*evc(1:npw,1:nbnd)
    !
    evc(1:npw,1:nbnd) = Y(1:npw,1:nbnd)
    !
    Y(1:npw,1:nbnd) = Ynew(1:npw,1:nbnd)
    !
    sigma = sigma2

  ENDDO

  !IF(ionode) WRITE(stdout,*) 'Before ortho:'
  !CALL test_orthonorm( Y, npwx, npw, nbnd )
  ! The filtered vector are in Y
  ! Orthonormalize Y
  CALL ortho_gram_schmidt( Y, npwx, npw, nbnd )
  !CALL ortho_qr( Y, npwx, npw, nbnd )
  !
  !IF(ionode) WRITE(stdout,*) 'After ortho:'
  !CALL test_orthonorm( Y, npwx, npw, nbnd )


  !
  ! Compute Rayleigh quotient to obtain the eigenvectors
  !
  ! Ynew <-- H*Y
  CALL h_psi( npwx, npw, nbnd, Y, Ynew)
  ! T <-- <Y|H*Y>
  CALL zgemm( 'Conjg', 'No', nbnd, nbnd, npw, (1.d0,0.d0), Y,npwx, Ynew,npwx, &
    (0.d0,0.d0), T,nbnd )
  ! FIXME Symmetrize T?
  CALL eig_zheevd( T,nbnd, eval_T,nbnd )
  ! evc <-- Y*T
  CALL zgemm('N','N', npw,nbnd,nbnd, (1.0_DP,0.0_DP), Y,npwx, T,nbnd, (0.d0,0.d0),evc,npwx)
  ! Copy eigenvalues
  DO i=1,nbnd
    et(i,ik) = eval_T(i)
  ENDDO
  ! Is this required?
  !IF(ionode) WRITE(stdout,*) 'Orthonormalizing evc'
  !CALL ortho_gram_schmidt( evc, npwx, npw, nbnd )
  !CALL test_orthonorm( evc, npwx, npw, nbnd )

  !
  DEALLOCATE( Y, Ynew, T, eval_T )

END SUBROUTINE



!----------------------------------------------
SUBROUTINE test_orthonorm(vec, lda, nrow, ncol)
!----------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode, stdout
  IMPLICIT NONE
  INTEGER, INTENT(in) :: lda, nrow, ncol
  COMPLEX(8), INTENT(in) :: vec(lda,ncol)
  !
  REAL(DP) :: norm1
  INTEGER :: ic
  !
  COMPLEX(DP) :: zdotc
  
  IF(ionode) WRITE(stdout,*) '>>>>> test_orthonorm <<<<<'

  IF(ionode) WRITE(stdout,*) 'Norms2:'
  DO ic=1,ncol
    norm1 = real( zdotc( nrow, vec(1,ic),1, vec(1,ic),1 ), kind=DP )
    IF(ionode) WRITE(stdout,*) ic, norm1
  ENDDO
  !
  IF(ionode) WRITE(stdout,*) 'Norms2 wrt to state1'
  DO ic=2,ncol
    norm1 = real( zdotc( nrow, vec(1,ic),1, vec(1,1),1 ), kind=DP )
    IF(ionode) WRITE(stdout,*) ic, norm1
  ENDDO


END SUBROUTINE


!------------------------------
SUBROUTINE test_orthonorm_evc()
!------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode, stdout
  USE wavefunctions_module, ONLY : evc
  USE wvfct, ONLY : npwx, npw, current_k, nbnd
  !
  REAL(DP) :: norm1
  INTEGER :: ib
  !
  COMPLEX(DP) :: zdotc

  IF(ionode) THEN
    WRITE(stdout,*) '>>>> test_orthonorm_evc <<<<'
    WRITE(stdout,*) 'current_k = ', current_k
    WRITE(stdout,*) 'npwx, npw = ', npwx, npw
  ENDIF
  IF(ionode) WRITE(stdout,*) 'evc-Norms:'
  DO ib=1,nbnd
    norm1 = sqrt( real( zdotc( npw, evc(1:npw,ib),1, evc(1:npw,ib),1 ), kind=DP ) )
    IF(ionode) WRITE(stdout,*) ib, norm1
  ENDDO
!  IF(ionode) WRITE(stdout,*) 'evc-Norms wrt to state1'
!  DO ib=2,nbnd
!    norm1 = sqrt( real( zdotc( npw, evc(1:npw,ib),1, evc(1:npw,1),1 ), kind=DP ) )
!    IF(ionode) WRITE(stdout,*) ib, norm1
!  ENDDO

END SUBROUTINE



!------------------------------------------------
SUBROUTINE ortho_gram_schmidt(v, ldv, nrow, ncol)
!------------------------------------------------
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ldv, nrow, ncol
  COMPLEX(DP), INTENT(inout) :: v(ldv,ncol)
  !
  INTEGER :: ii, jj
  COMPLEX(DP) :: zz, puv
  !
  COMPLEX(DP) :: zdotc

  DO ii = 1, ncol
    zz = zdotc( nrow, v(1:nrow,ii),1, v(1:nrow,ii),1 )
    v(1:nrow,ii) = v(1:nrow,ii)/sqrt( REAL( zz,kind=DP ) )
    DO jj = ii+1, ncol
      puv = prj( nrow, v(1:nrow,ii), v(1:nrow,jj) )
      v(1:nrow,jj) = v(1:nrow,jj) - puv*v(1:nrow,ii)
    ENDDO
  ENDDO

  CONTAINS

    ! compute prj = <v|u>/<u|u>
    FUNCTION prj(N,u,v)
      IMPLICIT NONE
      !
      COMPLEX(DP) :: prj
      INTEGER :: N
      COMPLEX(DP) :: u(N), v(N)
      !
      COMPLEX(DP) :: vu, uu
      COMPLEX(DP) :: zdotc
      !
      ! FIXME: I got the vectors to be orthogonal when I reverse the arguments
      ! for zdotc
      vu = zdotc( N, u,1, v,1 )
      uu = zdotc( N, u,1, u,1 )
      prj = vu/uu
    END FUNCTION

END SUBROUTINE


!------------------------------------------------
SUBROUTINE ortho_gram_schmidt_old(v, ldv, nrow, ncol)
!------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ldv, nrow, ncol
  COMPLEX(DP), INTENT(inout) :: v(ldv,ncol)
  !
  COMPLEX(DP), ALLOCATABLE :: w(:,:), ww(:), wvn(:)
  INTEGER :: ic, icc
  ! BLAS
  COMPLEX(DP) :: zdotc

  ALLOCATE( w(nrow,ncol) )
  ALLOCATE( ww(ncol) )
  ALLOCATE( wvn(ncol-1) )

  w(:,:) = (0.0_DP, 0.0_DP)
  ww(:) = (0.0_DP, 0.0_DP)
  
  w(1:nrow,1) = v(1:nrow,1)

  DO ic=2,ncol
    ww(ic-1) = zdotc( nrow, w(1:nrow,ic-1),1, w(1:nrow,ic-1),1 )
    IF(ionode) WRITE(stdout,*) 'Gram-Schmidt: ic-1, ww(ic-1)', ic-1, real(ww(ic-1),kind=DP)
    DO icc=1,ic-1
      wvn(icc) = zdotc( nrow, w(1:nrow,icc),1, v(1:nrow,ic),1 )
      w(1:nrow,ic) = w(1:nrow,ic) - wvn(icc)/real(ww(ic-1),kind=DP)*w(:,icc)
    ENDDO
    w(1:nrow,ic) = v(1:nrow,ic) + w(1:nrow,ic)   ! FIXME beware of this
  ENDDO

  ww(ncol) = zdotc( nrow, w(1:nrow,ncol),1, w(1:nrow,ncol),1 )
  DO ic=1,ncol
    v(1:nrow,ic) = w(1:nrow,ic)/sqrt( real( ww(ic), kind=8 ) )
  ENDDO

  DEALLOCATE( w, ww, wvn )
END SUBROUTINE


!
! FIXME: Test the performance of this subroutine
!
SUBROUTINE ortho_qr(X, npwx, npw, nbnd)
  USE kinds, ONLY: DP
  USE io_global, ONLY: stdout
  IMPLICIT NONE
  ! arguments
  INTEGER :: npwx, npw, nbnd
  COMPLEX(DP) :: X(npwx,nbnd)
  !
  ! Local arrays
  !
  COMPLEX(DP), ALLOCATABLE :: tau(:), work(:)
  !
  ! Local scalars
  !
  INTEGER :: lwork
  INTEGER :: info

  lwork = nbnd

  ALLOCATE( tau(nbnd) )  ! double the space for tau
  ALLOCATE( work(lwork) )

  !
  ! Compute QR factorization
  !
  CALL zgeqrf(npw, nbnd, X,npwx, tau, work,lwork, info)
  IF( info /= 0 ) THEN
    WRITE(stdout,*)
    WRITE(stdout,*) 'Error calling zgeqrf'
    RETURN
  ENDIF
  !
  ! Construct Q explicitly
  !
  CALL zungqr( npw,nbnd,nbnd, X,npwx, tau, work,lwork, info)
  IF( info /= 0 ) THEN
    WRITE(stdout,*)
    WRITE(stdout,*) 'Error calling zungqr'
    RETURN
  ENDIF

  DEALLOCATE( tau )
  DEALLOCATE( work )
END SUBROUTINE



!----------------------------------------
subroutine eig_zheevd(A,lda, lambda,N)
!----------------------------------------
  IMPLICIT NONE
  ! Arguments
  INTEGER :: ldA,N
  REAL(8) :: lambda(N)
  COMPLEX(8) :: A(ldA,N)
  ! Local
  INTEGER :: lwork, lrwork, liwork
  COMPLEX(8), ALLOCATABLE :: work(:)
  REAL(8), ALLOCATABLE :: rwork(:)
  INTEGER, ALLOCATABLE :: iwork(:)
  INTEGER :: info

  lwork = N*N + 2*N
  lrwork = 2*N*N + 5*N + 1
  liwork = 5*N + 3
  
  ALLOCATE(work(lwork))
  ALLOCATE(rwork(lrwork))
  ALLOCATE(iwork(liwork))

  !WRITE(*,*) 'sum(A) = ', sum(A)
  call zheevd('V','U', N, A,ldA, lambda, work,lwork, rwork,lrwork, iwork,liwork, info)
  IF(info /= 0) THEN
    WRITE(*,'(1x,A,I4)') 'ERROR calling ZHEEVD : info = ', info
    RETURN
  ENDIF

  DEALLOCATE(work)
  DEALLOCATE(rwork)
  DEALLOCATE(iwork)
END SUBROUTINE



!
! Set lower bound and upper bound for CheFSI procedure
!
! Algorithm 4.4 of Ref. (1)
!
!------------------------------------------------------------------------------
SUBROUTINE kstep_lanczos(ik, npw, nlancz, lb, ub)
!------------------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode, stdout
  USE wvfct, ONLY : npwx, nbnd, et
  IMPLICIT NONE
  INTEGER, INTENT(in) :: ik  ! current kpoint
  INTEGER, INTENT(in) :: npw  ! no. of PW of this kpoint
  INTEGER, INTENT(in) :: nlancz ! no. of Lanczos iteration
  REAL(DP), INTENT(out) :: lb, ub  ! lower and upper bound
  ! Local
  COMPLEX(DP), ALLOCATABLE :: v(:), f(:), v0(:)
  REAL(DP), ALLOCATABLE :: T(:,:), eval_T(:)
  !
  REAL(DP) :: rnd1, rnd2, alpha, beta
  COMPLEX(DP) :: tz
  !
  INTEGER :: ipw, j
  ! For dsyev
  INTEGER :: lwork, info
  REAL(DP), ALLOCATABLE :: work(:)
  !
  COMPLEX(DP) :: zdotc
  REAL(DP) :: induced_norm

  ALLOCATE( v(npwx), v0(npwx) )
  ALLOCATE( f(npwx) )
  ALLOCATE( T(nlancz,nlancz) )
  ALLOCATE( eval_T(nlancz) )

  lwork = 3*nlancz - 1
  ALLOCATE( work(lwork) )

  v(:)   = (0.0_DP, 0.0_DP)  ! FIXME: Need this?
  v0(:)  = (0.0_DP, 0.0_DP)  ! FIXME: Need this?
  f(:)   = (0.0_DP, 0.0_DP)  ! FIXME: Need this?
  T(:,:) = 0.0_DP

  !
  ! Random initial
  !
  DO ipw=1,npw
    CALL random_number( rnd1 )
    CALL random_number( rnd2 )
    v(ipw) = cmplx(rnd1, rnd2, kind=DP)
  ENDDO
  tz = zdotc( npw, v(:),1, v(:),1 )
  v(:) = v(:)/sqrt( real(tz,kind=DP) )  ! FIXME use BLAS 

  ! 
  ! Step no. 2
  !
  CALL h_psi( npwx, npw,1, v(:), f(:) ) ! f <-- H v
  alpha = real( zdotc(npw, f(:),1, v(:),1 ), kind=DP )    ! alpha <-- <f|v>
  f(1:npw) = f(1:npw) - alpha*v(1:npw)
  T(1,1) = alpha

  ! Step. 3
  DO j=2,nlancz
    ! Step 4
    beta = sqrt( real( zdotc( npw, f(:),1, f(:),1 ), kind=DP ) )
    IF( beta < 1.d-06 ) THEN
      IF(ionode) WRITE(stdout,*) 'WARNING: small beta in kstep_lanczos'
    ENDIF
    ! Step 5
    v0(1:npw) = v(1:npw)
    v(1:npw) = f(1:npw)/beta
    ! Step 6
    CALL h_psi( npwx, npw,1, v(:), f(:) )
    f(1:npw) = f(1:npw) - beta*v0(1:npw)
    ! Step 7
    alpha = real( zdotc(npw, f(:),1, v(:),1 ), kind=DP )    ! alpha <-- <f|v>
    f(1:npw) = f(1:npw) - alpha*v(1:npw)
    ! Step 8
    T(j,j-1) = beta
    T(j-1,j) = beta
    T(j,j)   = alpha
  ENDDO

  lb = et(nbnd,ik) + 1.0_DP
  alpha = sqrt( real( zdotc(npw, f(:),1, f(:),1 ), kind=DP ) ) ! now alpha has different meaning
  !ub = induced_norm(T, nlancz, nlancz) + alpha
  CALL dsyev('N', 'U', nlancz, T,nlancz, eval_T, work, lwork, info)
  ub = eval_T(nlancz) + alpha

  !IF(ionode) THEN
  !  WRITE(stdout,*) '>>>>> kstep_lanczos <<<<<<'
  !  WRITE(stdout,*) 'lb = ', lb
  !  WRITE(stdout,*) 'up = ', ub
  !ENDIF

  DEALLOCATE( v, f, v0 )
  DEALLOCATE( T, work, eval_T )

END SUBROUTINE


!------------------------------------------------------------------------------
FUNCTION induced_norm( A, ldA, N) RESULT(norm)
!------------------------------------------------------------------------------
  USE kinds, ONLY: DP
  IMPLICIT NONE
  ! Arguments
  INTEGER :: ldA
  INTEGER :: N
  REAL(DP) :: A(ldA,N)
  ! Result
  REAL(DP) :: norm
  ! Local array
  REAL(DP), ALLOCATABLE :: v1(:), v(:)
  ! Local scalar
  INTEGER :: i
  REAL(DP) :: d
  ! BLAS
  REAL(DP) :: ddot

  ALLOCATE( v1(N), v(N) )

  ! Create a unit norm vector
  d = 1.d0/sqrt( dble(N) )
  DO i=1,N
    v1(i) = d
  ENDDO
  ! Multiply
  CALL dgemv( 'No', N,N, 1.d0,A,ldA, v1,1, 0.d0,v,1 )
  ! Calculate norm
  norm = sqrt( ddot( N, v,1, v,1 ) )

  DEALLOCATE( v1, v )

END FUNCTION


