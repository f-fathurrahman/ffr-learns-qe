!------------------------------------
SUBROUTINE my_pseudo_q(qfunc, qfuncl)
!------------------------------------       
  USE kinds, ONLY : DP
  USE ld1_parameters, ONLY : nwfsx
  USE ld1inc, ONLY : rcut, lls,  grid, ndmx, lmx2, nbeta, ikk, ecutrho, &
              rmatch_augfun, rmatch_augfun_nc
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: qfunc(ndmx,nwfsx,nwfsx)
  REAL(DP), INTENT(OUT) :: qfuncl(ndmx,nwfsx,nwfsx,0:lmx2) 
  REAL(DP),  EXTERNAL :: int_0_inf_dr
  !
  ! variables for aug. functions generation
  ! 
  INTEGER  :: irc, ns, ns1, l1, l2, l3, lll, mesh, n, ik
  INTEGER  :: l1_e, l2_e
  REAL(DP) :: aux(ndmx)
  REAL(DP) :: augmom, ecutrhoq, rmatch

  write(*,*)
  write(*,*) '<div> ENTER my_pseudo_q'
  write(*,*)
  write(*,*) 'rmatch_augfunc_nc = ', rmatch_augfun_nc

  ecutrho = 0.0_DP
  mesh = grid%mesh
  qfuncl = 0.0_DP
  do ns = 1,nbeta
    l1 = lls(ns)
    do ns1 = ns,nbeta
      l2 = lls(ns1)
      !
      ! Find the matching point
      !
      ik = 0
      IF (rmatch_augfun_nc) THEN
        rmatch = min(rcut(ns),rcut(ns1))
      ELSE
        rmatch = rmatch_augfun
      ENDIF
      !
      do n = 1,mesh
        if( grid%r(n) > rmatch ) then
          ik = n
          exit
        endif
      enddo
      IF( ik==0 .or. ik > mesh-20) THEN
        call errore('pseudo_q', 'wrong rmatch_augfun', 1)
      ENDIF
      !
      ! Do the pseudization
      !
      do l3 = abs(l1-l2), l1+l2, 2
        CALL my_compute_q_3bess(l3, l1+l2, ik, qfunc(1,ns,ns1), qfuncl(1,ns,ns1,l3), ecutrhoq)
        IF( ecutrhoq > ecutrho) then
          ecutrho = ecutrhoq
          l1_e = l1
          l2_e = l2
        ENDIF
        qfuncl(1:mesh,ns1,ns,l3) = qfuncl(1:mesh,ns,ns1,l3)
      enddo
    enddo
  enddo
  !
  !  Check that multipoles have not changed
  !
  irc = maxval(ikk(1:nbeta)) + 8
  augmom = 0.0_DP
  DO ns = 1,nbeta
    l1 = lls(ns)
    DO ns1 = ns,nbeta
      l2 = lls(ns1)
      DO l3 = abs(l1-l2), l1+l2, 2
        aux(1:irc) = (qfuncl(1:irc,ns,ns1,l3) - qfunc(1:irc,ns,ns1)) * grid%r(1:irc)**l3
        lll = l1 + l2 + 2 + l3
        augmom = int_0_inf_dr(aux(1:irc),grid,irc,lll)
        IF( abs(augmom) > 1.d-5) THEN
          WRITE (*,'(5x,a,2i3,a,2i3,a,i3,f15.7)') " Problem with multipole",ns,l1,":",ns1,l2, " l3=",l3, augmom
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  WRITE(*,'(/,5x, "Q pseudized with Bessel functions")')
  WRITE(*,'(5x,"Expected ecutrho= ",f12.4," due to l1=",i3,"   l2=",i3)') ecutrho, l1_e, l2_e
  
  write(*,*)
  write(*,*) '</div> EXIT my_pseudo_q'
  write(*,*)
  
  RETURN
END SUBROUTINE


!--------------------------------------------------------------------------
subroutine my_compute_q_3bess(ldip, lam, ik, chir, phi_out, ecutrho)
!--------------------------------------------------------------------------
  !
  ! This routine computes the phi_out function by pseudizing the
  ! chir function with a linear combination of three Bessel functions
  ! multiplied by r**2. In input it receives the point
  ! ik where the cut is done, the angular momentum lam of the 
  ! bessel functions and the function chir.
  ! Phi_out has the same ldip dipole moment of chir.
  !
  use kinds, only : DP
  use radial_grids, only: ndmx
  use ld1inc, only: grid
  implicit none

  integer ::    &
       ldip,    & ! input: the order of the dipole
       lam,     & ! input: the angular momentum
       ik         ! input: the point corresponding to rc

  real(DP) :: &
       xc(8)      ! output: the coefficients of the Bessel functions

  real(DP) ::         &
       chir(ndmx),    &   ! input: the all-electron function
       phi_out(ndmx)      ! output: the phi function
  !
  real(DP) ::  &
       ecutrho,& ! the expected cut-off on the charge density for this q
       fae,    & ! the value of the all-electron function
       f1ae,   & ! its first derivative
       f2ae,   & ! the second derivative
       dip       ! the norm of the function

  integer ::    &
       n, nst, nc

  real(DP) :: &
       gi(ndmx), j1(ndmx,3), jnor, cm(3), bm(3), delta, gam

  real(DP), external :: deriv_7pts, deriv2_7pts, int_0_inf_dr

  integer ::  &
       iok,   &  ! flag
       nbes      ! number of Bessel functions to be used

  nbes = 3
  !
  nst = lam + 2 + ldip
  !
  ! compute the first and second derivative of input function at r(ik)
  !
  fae = chir(ik)
  f1ae = deriv_7pts(chir, ik, grid%r(ik), grid%dx)
  f2ae = deriv2_7pts(chir, ik, grid%r(ik), grid%dx)
  !
  ! compute the ldip dipole moment of the input function
  !
  do n=1,ik+1
    gi(n) = chir(n) * grid%r(n)**ldip  
  enddo
  dip = int_0_inf_dr(gi, grid, ik, nst)
  !
  ! RRKJ: the pseudo-wavefunction is written as an expansion into 3  
  !       spherical Bessel functions for r < r(ik)
  ! find q_i with the correct log derivatives
  !   
  call find_qi(f1ae/fae, xc(nbes+1), ik, ldip, nbes, 2, iok)
  if( iok .ne. 0 ) then
    call errore('compute_q_3bess', 'problem with the q_i coefficients', 1)
  endif
  !
  !   compute the Bessel functions and multiply by r**2
  !
  do nc = 1,nbes
    call sph_bes(ik + 5, grid%r, xc(nbes+nc), ldip, j1(1,nc))
    jnor = j1(ik,nc)*grid%r2(ik)
    do n = 1,ik+5
      j1(n,nc) = j1(n,nc)*grid%r2(n)*chir(ik)/jnor
    enddo
  enddo
  !
  ! compute the bm functions (second derivative of the j1)
  ! and the ldip dipole moment of the Bessel function (cm)
  do nc = 1, nbes
    bm(nc) = deriv2_7pts(j1(1,nc), ik, grid%r(ik), grid%dx)
    do n = 1,ik
      gi(n) = j1(n,nc)*grid%r(n)**ldip
    enddo
    cm(nc) = int_0_inf_dr(gi,grid,ik,nst)
  enddo
  !
  ! solve the linear system to find the coefficients
  gam = (bm(3) - bm(1))/(bm(2) - bm(1))
  delta = (f2ae - bm(1))/(bm(2) - bm(1))
  !   
  xc(3) = (dip - cm(1) + delta*(cm(1) - cm(2)))/(gam*(cm(1)-cm(2))+cm(3)-cm(1))
  xc(2) = -xc(3)*gam + delta
  xc(1) = 1.0_dp - xc(2) - xc(3)
  !
  ! Set the function for r<=r(ik)
  do n = 1,ik
     phi_out(n) = xc(1)*j1(n,1) + xc(2)*j1(n,2) + xc(3)*j1(n,3)
  enddo
  !
  ! for r > r(ik) the function does not change
  !
  do n = ik+1,grid%mesh
    phi_out(n) = chir(n)
  enddo
  ecutrho = 2.0_dp*xc(6)**2

  return
end subroutine
