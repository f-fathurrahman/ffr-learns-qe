INCLUDE 'prepare_all.f90'

PROGRAM main
  USE uspp_param, ONLY: lmaxkb
  IMPLICIT NONE 
  CALL prepare_all()
  CALL my_aainit(lmaxkb + 1)
END PROGRAM

!-----------------------------------------------------------------------
SUBROUTINE my_aainit(lli)
!-----------------------------------------------------------------------
!
  ! this routine computes the coefficients of the expansion of the product
  ! of two real spherical harmonics into real spherical harmonics.
  !
  !     Y_limi(r) * Y_ljmj(r) = \sum_LM  ap(LM,limi,ljmj)  Y_LM(r)
  !
  ! On output:
  ! ap     the expansion coefficients
  ! lpx    for each input limi,ljmj is the number of LM in the sum
  ! lpl    for each input limi,ljmj points to the allowed LM
  !
  ! The indices limi,ljmj and LM assume the order for real spherical
  ! harmonics given in routine ylmr2
  !
  USE upf_kinds, ONLY : DP
  USE upf_invmat
  USE upf_params, ONLY: lmaxx, lqmax
  USE uspp_param, ONLY: lmaxkb
  IMPLICIT NONE 
  !
  ! input: the maximum li considered
  !  
  INTEGER :: lli
  !
  ! local variables
  !
  INTEGER :: llx, l, li, lj
  REAL(DP) , ALLOCATABLE :: r(:,:), rr(:), ylm(:,:), mly(:,:)
  ! an array of random vectors: r(3,llx)
  ! the norm of r: rr(llx)
  ! the real spherical harmonics for array r: ylm(llx,llx)
  ! the inverse of ylm considered as a matrix: mly(llx,llx)
  !
  ! Some variables in uspp module
  INTEGER :: nlx, mx
  INTEGER, ALLOCATABLE :: lpx(:,:), lpl(:,:,:)
  REAL(DP), ALLOCATABLE :: ap(:,:,:)
  !
  ! Function
  REAL(DP) :: my_compute_ap

  ! maximum number of combined angular momentum
  nlx = (lmaxx+1)**2
  ! maximum magnetic angular momentum of Q
  mx = 2*lqmax-1
  ! for each pair of combined momenta lm(1),lm(2): 
  ALLOCATE( lpx(nlx,nlx) ) ! maximum combined angular momentum LM
  ALLOCATE( lpl(nlx,nlx,mx) )  ! list of combined angular momenta  LM
  
  WRITE(*,*) 'lmaxkb = ', lmaxkb
  WRITE(*,*) 'lli (maximum li considered) = ', lli

  WRITE(*,*)
  WRITE(*,*) 'Some hardcoded parameters'
  WRITE(*,*) '-------------------------'
  WRITE(*,*) 'lmaxx (max nonlocal proj angmom) = ', lmaxx
  WRITE(*,*) 'nlx (max no. of combined angmom) = ', nlx
  WRITE(*,*) 'lqmax (max no. of angmom of Q)   = ', lqmax
  WRITE(*,*) 'mx (max magnetic angmom of Q)    = ', mx

  ! Clebsch-Gordan coefficients for spherical harmonics
  ALLOCATE( ap(lqmax*lqmax,nlx,nlx) )

  IF( lli < 0 ) CALL upf_error('aainit', 'lli not allowed',lli)

  IF( lli*lli > nlx ) CALL upf_error('aainit','nlx is too small ',lli*lli)

  llx = (2*lli-1)**2
  IF( 2*lli-1 > lqmax ) CALL upf_error('aainit', 'ap leading dimension is too small', llx)

  WRITE(*,*) 'llx = ', llx

  ALLOCATE( r(3, llx) )
  ALLOCATE( rr(llx) )
  ALLOCATE( ylm(llx,llx) )
  ALLOCATE( mly(llx,llx) )

  r(:,:)   = 0.0_DP
  ylm(:,:) = 0.0_DP
  mly(:,:) = 0.0_DP
  ap(:,:,:)= 0.0_DP

  ! generate an array of random vectors (uniform deviate on unitary sphere)
  CALL my_gen_rndm_r(llx, r, rr)
  WRITE(*,*)
  WRITE(*,*) 'Randomly generated unit vectors'
  DO li = 1, llx
    WRITE(*,'(1x,I4,4F18.10)') li, r(1,li), r(2,li), r(3,li), rr(li)
  ENDDO 

  ! generate the real spherical harmonics for the array: ylm(ir,lm)
  CALL ylmr2(llx, llx, r, rr, ylm)

  ! store the inverse of ylm(ir,lm) in mly(lm,ir)
  CALL invmat(llx, ylm, mly)

  ! for each li,lj compute ap(l,li,lj) and the indices, lpx and lpl
  WRITE(*,*)
  DO li = 1, lli*lli
    WRITE(*,*)
    DO lj = 1, lli*lli
      lpx(li,lj) = 0
      DO l = 1, llx
        ap(l,li,lj) = my_compute_ap(l,li,lj,llx,ylm,mly)
        IF( abs(ap(l,li,lj) ) > 1.d-3 ) THEN 
          lpx(li,lj) = lpx(li,lj) + 1 ! increment
          IF(lpx(li,lj) > mx) CALL upf_error('aainit', 'mx dimension too small', lpx(li,lj))
          lpl(li,lj,lpx(li,lj)) = l
          WRITE(*,'(1x,3I4,F18.10)') l, li, lj, ap(l, li, lj)          
        ENDIF 
      ENDDO 
    ENDDO 
  ENDDO 

  DEALLOCATE(mly)
  DEALLOCATE(ylm)
  DEALLOCATE(rr)
  DEALLOCATE(r)
  DEALLOCATE(lpl)
  DEALLOCATE(lpx)

  RETURN 
END SUBROUTINE 


!-----------------------------------------------------------------------
FUNCTION my_compute_ap(l,li,lj,llx,ylm,mly)
!-----------------------------------------------------------------------
  !-  given an l and a li,lj pair compute ap(l,li,lj)
  USE upf_kinds, ONLY: DP
  implicit none
  !
  ! first the I/O variables
  !
  INTEGER :: llx ! the dimension of ylm and mly
  INTEGER :: l,li,lj       ! the arguments of the array ap
  
  REAL(DP) :: my_compute_ap ! this function
  REAL(DP) :: ylm(llx,llx)  ! the real spherical harmonics for array r
  REAL(DP) :: mly(llx,llx)  ! the inverse of ylm considered as a matrix
  !
  ! here the local variables
  !
  INTEGER :: ir
    
  my_compute_ap = 0.0_DP
  DO ir = 1,llx
    my_compute_ap = my_compute_ap + mly(l,ir)*ylm(ir,li)*ylm(ir,lj)
  ENDDO
    
  RETURN 
END FUNCTION 



!-----------------------------------------------------------------------
SUBROUTINE my_gen_rndm_r(llx, r, rr)
!-----------------------------------------------------------------------
  ! - generate an array of random vectors (uniform deviate on unitary sphere)
  !
  USE upf_kinds, ONLY : DP
  USE upf_const,  ONLY: tpi
  IMPLICIT NONE 
  
  !
  ! first the I/O variables
  !
  INTEGER :: llx   ! input: the dimension of r and rr
  REAL(DP) :: r(3,llx) ! output: an array of random vectors
  REAL(DP) :: rr(llx)     ! output: the norm of r
  !
  ! here the local variables
  !
  INTEGER :: ir
  REAL(DP) :: costheta, sintheta, phi
  REAL(DP) :: my_randy ! FUNCTION

  DO ir = 1, llx
    costheta = 2.0_DP * my_randy() - 1.0_DP
    sintheta = SQRT ( 1.0_DP - costheta*costheta)
    phi = tpi * my_randy()
    r(1,ir) = sintheta * cos(phi)
    r(2,ir) = sintheta * sin(phi)
    r(3,ir) = costheta
    rr(ir)   = 1.0_DP
  ENDDO 
  RETURN 
END SUBROUTINE 

!------------------------------------------------------------------------
FUNCTION my_randy()
!------------------------------------------------------------------------
  !   
  ! x=randy(n): reseed with initial seed idum=n ( 0 <= n <= ic, see below)
  !             if randy is not explicitly initialized, it will be
  !             initialized with seed idum=0 the first time it is called
  ! x=randy() : generate uniform real(DP) numbers x in [0,1]
  !   
  USE upf_kinds, only : DP
  IMPLICIT NONE 
  !
  REAL(DP) :: my_randy
  !INTEGER, optional    :: irand
  !   
  INTEGER , PARAMETER  :: m    = 714025, &
                          ia   = 1366, &
                          ic   = 150889, &
                          ntab = 97
  REAL(DP), PARAMETER  :: rm = 1.0_DP / m 
  INTEGER              :: j
  INTEGER, SAVE        :: ir(ntab), iy, idum=0
  LOGICAL, SAVE        :: first=.true.
  !   
  
  !IF( present(irand) ) THEN
  !  idum = MIN( ABS(irand), ic) 
  !  first=.true.
  !ENDIF

  IF ( first ) THEN
    first = .false.
    idum = MOD( ic - idum, m ) 
    DO j=1,ntab
      idum = mod(ia*idum+ic,m)
      ir(j) = idum
    ENDDO
    idum = mod(ia*idum+ic,m)
    iy = idum
  ENDIF

  j = 1 + (ntab*iy)/m
  IF( j > ntab .OR. j <  1 ) CALL upf_error('my_randy', 'j out of range', ABS(j)+1)
  iy = ir(j)
  my_randy = iy*rm
  idum = mod(ia*idum+ic,m)
  ir(j) = idum
  RETURN
END FUNCTION 

