!--------------------------------------------------------------------------
subroutine my_ascheqps_drv(veff, ncom, thresh, flag_all, nerr)
!--------------------------------------------------------------------------
  ! This routine is a driver that calculates for the test
  ! configuration the solutions of the Kohn and Sham equation
  ! with a fixed pseudo-potential. The potentials are assumed
  ! to be screened. The effective potential veff is given in input.
  ! The output wavefunctions are written in phits and are normalized.
  ! If flag is .true. compute all wavefunctions, otherwise only
  ! the wavefunctions with positive occupation.
  !      
  use kinds, only: dp
  use ld1_parameters, only: nwfsx
  use radial_grids, only: ndmx
  use ld1inc, only: grid, pseudotype, rel, &
                    lls, jjs, qq, ikk, ddd, betas, nbeta, vnl, &
                    nwfts, iswts, octs, llts, jjts, nnts, enlts, phits, nbeta
  implicit none

  integer ::    &
          nerr, &     ! control the errors of the routine ascheqps
          ncom        ! number of components of the pseudopotential

  real(DP) :: &
       veff(ndmx,ncom)    ! work space for writing the potential 

  logical :: flag_all    ! if true calculates all the wavefunctions

  integer ::  &
       ns,    &  ! counter on pseudo functions
       is,    &  ! counter on spin
       nbf,   &  ! auxiliary nbeta
       n,     &  ! index on r point
       nstop, &  ! errors in each wavefunction
       ind

  real(DP) :: &
       vaux(ndmx,2)     ! work space for writing the potential 

  real(DP) :: thresh         ! threshold for selfconsistency
  integer :: i, j
  
  write(*,*)
  write(*,*) '<div> ENTER my_ascheqps_drv'
  write(*,*)
  write(*,*) 'thresh = ', thresh

  !
  ! compute the pseudowavefunctions in the test configuration
  !
  if (pseudotype == 1) then
    nbf = 0
  else
    nbf = nbeta
  endif

  nerr = 0
  ! Loop over all states for test
  do ns = 1,nwfts
    if( octs(ns) > 0.0_dp .or. ( octs(ns) > -1.0_dp .and. flag_all ) ) then
      write(*,*)
      write(*,*) 'Calling my_ascheqps for input configuration'
      write(*,*) 'ns, nnts, llts, jjts'
      write(*,'(1x,3I3,F5.1)') ns, nnts(ns), llts(ns), jjts(ns)
      write(*,*)
      write(*,'(1x,A,F18.10)') 'At input: energy (in Ha) = ', enlts(ns)*0.5d0
      !
      is = iswts(ns)
      if( ncom==1 .and. is==2) then
        call errore('ascheqps_drv','incompatible spin',1)
      endif
      !
      if( pseudotype == 1 ) then
        !
        if( rel < 2 .or. llts(ns) == 0 .or. &
          & abs(jjts(ns)-llts(ns)+0.5_dp) < 0.001_dp) then
          ind = 1
        !
        elseif( rel == 2 .and. llts(ns) > 0 .and. &
              & abs(jjts(ns)-llts(ns)-0.5_dp) < 0.001_dp) then
          ind = 2
        else
          call errore('my_ascheqps_drv', 'unexpected case', 1)
        endif
        !
        do n = 1,grid%mesh
          vaux(n,is) = veff(n,is) + vnl(n,llts(ns),ind)
        enddo
      else
        !ffr: for other pseudotypes vaux is veff (V_Ps_loc)
        do n = 1,grid%mesh
          vaux(n,is) = veff(n,is)
        enddo
      endif
      !
      ! ddd is input here, should not be modified
      !write(*,*) 'Before my_ascheqps: ddd matrix = ' ! should be spin-dependent 
      !do i = 1,nbeta
      !  write(*,'(6f12.5)') (ddd(i,j,is), j=1, nbeta)
      !enddo
      !
      call my_ascheqps( nnts(ns), llts(ns), jjts(ns), enlts(ns), grid%mesh, ndmx, &
                    &   grid, vaux(1,is), thresh, phits(1,ns), betas, ddd(1,1,is), qq, nbf, &
                    &   nwfsx, lls, jjs, ikk, nstop)
      write(*,*)
      write(*,'(1x,A15,F18.10)') 'After my_ascheqps: energy (in Ha) = ', enlts(ns)*0.5d0
      !write(*,*)
      !write(*,*) 'After my_ascheqps: ddd matrix = ' ! should be spin-dependent 
      !do i = 1,nbeta
      !  write(*,'(6f12.5)') (ddd(i,j,is), j=1, nbeta)
      !enddo
      !
      ! normalize the wavefunctions 
      call normalize(phits(1,ns), llts(ns), jjts(ns), ns)
      !
      ! not sure whether the "best" error code should be like this:
      ! IF ( octs(ns) > 0.0_dp ) nerr = nerr + nstop
      !   i.e. only for occupied states, or like this:
      nerr = nerr + nstop
    endif ! if octs is larger than zero
  enddo

  write(*,*)
  write(*,*) '</div> EXIT my_ascheqps_drv'
  write(*,*)

  return
end subroutine

