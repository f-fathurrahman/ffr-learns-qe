!--------------------------------------------------------------------------
subroutine my_descreening()
!--------------------------------------------------------------------------
  !
  ! This routine descreens the local potential and the ddd
  ! coefficients (the latter only in the US case)
  ! The charge density is computed with the test configuration,
  ! not the one used to generate the pseudopotential
  !      
  use kinds, only: dp
  use mp,        only : mp_bcast
  use radial_grids, only: ndmx
  use ld1_parameters, only: nwfsx
  use ld1inc, only: grid, nlcc, vxt, lsd, vpstot, vpsloc, file_screen, &
                    vh, enne, rhoc, latt, rhos, enl, &
                    nbeta, bmat, qvan, qvanl, jjs, lls, ikk, pseudotype, &
                    nwfts, enlts, octs, llts, jjts, phits, nstoaets, &
                    which_augfun
  implicit none

  integer ::  &
       ns,    &  ! counter on pseudo functions
       ns1,   &  ! counter on pseudo functions
       ib,jb, &  ! counter on beta functions
       lam       ! the angular momentum

  real(DP) :: &
       vaux(ndmx,2)     ! work space

  real(DP), external :: int_0_inf_dr ! the integral function

  real(DP), parameter :: &
       thresh= 1.e-12_dp          ! threshold for selfconsistency

  integer  :: &
       n, nst, iwork(nwfsx), ios, nerr
  
  write(*,*)
  write(*,*) '<div> ENTER my_descreening'
  write(*,*)
  
  ! descreening the local potential: NB: this descreening is done with
  ! the occupation of the test configuration. This is required
  ! for pseudopotentials with semicore states. In the other cases
  ! a test configuration equal to the one used for pseudopotential
  ! generation is strongly suggested
  !
  do n = 1,nwfts
    enlts(n) = enl(nstoaets(n))
  enddo
  !
  ! compute the pseudowavefunctions in the test configuration
  !
  call my_ascheqps_drv(vpsloc, 1, thresh, .false., nerr)
  !
  ! descreening the D coefficients
  !
  if (pseudotype == 3) then
    do ib = 1,nbeta
      do jb = 1,ib
        if( lls(ib) == lls(jb) .and. abs(jjs(ib)-jjs(jb)) < 1.e-7_dp ) then
          lam = lls(ib)
          nst = (lam+1)*2
          IF( which_augfun == 'PSQ' ) then
            do n = 1,ikk(ib)
              vaux(n,1) = qvanl(n,ib,jb,0)*vpsloc(n)
            enddo
          ELSE
            do n = 1,ikk(ib)
              vaux(n,1) = qvan(n,ib,jb)*vpsloc(n)
            enddo
          ENDIF
          bmat(ib,jb) = bmat(ib,jb) - int_0_inf_dr(vaux(1,1),grid,ikk(ib),nst)
        endif
        bmat(jb,ib) = bmat(ib,jb)
      enddo ! do jb
    enddo ! do ib 
    write(*,'(/5x,'' The ddd matrix'')')
    do ns1 = 1,nbeta
      write(*,'(6f12.5)') (bmat(ns1,ns),ns=1,nbeta)
    enddo
  endif
  !
  ! descreening the local pseudopotential
  iwork = 1
  call my_chargeps(rhos, phits, nwfts, llts, jjts, octs, iwork)
  !
  call new_potential(ndmx, grid%mesh, grid, 0.0_dp, vxt, lsd, nlcc, latt, enne,&
       rhoc, rhos, vh, vaux, 1)

  do n = 1,grid%mesh
    vpstot(n,1) = vpsloc(n)
    vpsloc(n) = vpsloc(n) - vaux(n,1)
  enddo

  if (file_screen .ne.' ') then
    open(unit=20, file=file_screen, status='unknown', iostat=ios)
    do n = 1,grid%mesh
      write(20,'(i5,7e12.4)') n, grid%r(n), vpsloc(n)+vaux(n,1), vpsloc(n), vaux(n,1), rhos(n,1)
    enddo
    close(20)
  endif

  write(*,*)
  write(*,*) '</div> EXIT my_descreening'
  write(*,*)

  return
end subroutine


!---------------------------------------------------------------
subroutine my_chargeps(rho_i, phi_i, nwf_i, ll_i, jj_i, oc_i, iswf_i)
!---------------------------------------------------------------
  !
  ! calculate the (spherical) pseudo charge density 
  !
  use kinds, only: dp
  use ld1_parameters, only: nwfsx
  use radial_grids, only: ndmx
  use ld1inc, only: grid, pseudotype, qvan, nbeta, betas, lls, jjs, ikk,  &
                    which_augfun, qvanl, nspin
  implicit none

  integer :: &
       nwf_i,        & ! input: the number of wavefunctions
       iswf_i(nwfsx),& ! input: their spin
       ll_i(nwfsx)     ! input: their angular momentum

  real(DP) ::  &
       jj_i(nwfsx), & ! input: their total angular momentum
       oc_i(nwfsx), & ! input: the occupation
       phi_i(ndmx,nwfsx), & ! input: the functions to add
       rho_i(ndmx,2)   ! output: the (nspin) components of the charge

  integer ::     &
       is,     &   ! counter on spin
       n,n1,n2,&   ! counters on beta and mesh function
       ns,nst,ikl  ! counter on wavefunctions

  real(DP) ::    &
       work(nwfsx), & ! auxiliary variable for becp
       int_0_inf_dr,& ! integration function
       gi(ndmx)        ! used to compute the integrals


  rho_i = 0.0_dp
  !
  ! compute the square modulus of the eigenfunctions
  !
  do ns = 1,nwf_i
    if(oc_i(ns) > 0.0_dp) then
      is = iswf_i(ns)
      do n = 1,grid%mesh
        rho_i(n,is) = rho_i(n,is) + oc_i(ns)*phi_i(n,ns)**2
      end do
     endif
  enddo
  !
  ! if US pseudopotential compute the augmentation part
  !
  if( pseudotype == 3 ) then
    do ns = 1,nwf_i
      if( oc_i(ns) > 0.0_dp) then
        is = iswf_i(ns)
        do n1 = 1,nbeta
          if( ll_i(ns) == lls(n1).and. abs(jj_i(ns)-jjs(n1)) < 1.e-7_dp) then
            nst = (ll_i(ns) + 1)*2
            ikl = ikk(n1)
            do n = 1,ikl
              gi(n) = betas(n,n1)*phi_i(n,ns)
            enddo
            work(n1) = int_0_inf_dr(gi, grid, ikl, nst)
          else
            work(n1) = 0.0_dp
          endif
        enddo
        !
        ! and adding to the charge density
        !
        do n1 = 1,nbeta
          do n2 = 1,nbeta
            IF( which_augfun == 'PSQ' ) then
              do n = 1,grid%mesh
                rho_i(n,is) = rho_i(n,is) + qvanl(n,n1,n2,0)*oc_i(ns)*work(n1)*work(n2)
              enddo
            ELSE
              do n = 1,grid%mesh
                rho_i(n,is) = rho_i(n,is) + qvan(n,n1,n2)*oc_i(ns)*work(n1)*work(n2)
              enddo
            ENDIF
          enddo ! do n2
        enddo ! do n1
      endif
    enddo ! do ns
  endif
  !
  ! Check for negative charge
  !
  do is = 1,nspin
    do n = 2,grid%mesh !ffr: why start from 2 ?
      if( rho_i(n,is) < -1.d-12) then 
        call errore('chargeps','negative rho',1)
      endif
    enddo
  enddo

  return
end subroutine
