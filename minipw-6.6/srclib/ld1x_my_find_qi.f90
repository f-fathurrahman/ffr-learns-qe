!--------------------------------------------------------------------------
subroutine my_find_qi(logderae, xc, ik, lam, ncn, flag, iok)
!--------------------------------------------------------------------------
  !
  ! This routine finds three values of q such that the
  ! functions f_l have a logarithmic derivative equal to
  ! logderae at the point ik
  !
  !  f_l = j_l(r) * r**flag
  !
  use kinds, only: dp
  use ld1inc, only: grid
  implicit none
  integer, parameter :: ncmax=10   ! maximum allowed nc
  integer ::      &
       ik,    & ! input: the point corresponding to rcut
       ncn,   & ! input: the number of qi to compute
       flag,  & ! input: the type of function
       iok,   & ! output: if 0 the calculation in this routine is ok
       lam      ! input: the angular momentum

  real(DP) :: &
       xc(ncn), & ! output: the values of qi
       logderae  ! input: the logarithmic derivative

  real(DP) ::   &
       j1(ncmax), & ! the bessel function in three points
       qmax,qmin, & ! the limits of the q search
       logdermax, logdermin, & ! the maximum and minimum logder
       logder, & ! the actual logder
       jlmin, jlmax, & ! the value of jl in qmin and qmax
       my_compute_log, &! function for log derivative
       dq, dq_0 ! the step to braket the q

  integer ::    &
       nc,  &    ! counter on the q found
       icount, &  ! too many iterations
       icount1, & ! too many iterations
       imax,&   ! maximum number of iteration to braket
       iq      ! counter on iteration

  iok = 0
  if(ncn > ncmax) then
    call errore('find_qi','ncn is too large',1)
  endif

  if( flag == 0 .and. lam /= 0) then
    call errore('find_qi','lam too large for this iflag',1)
  endif

  if( lam > 6) then
    call errore('find_qi','l not programmed',1)
  endif
  !
  ! fix deltaq and the maximum step number
  dq_0 = 0.05d0
  imax = 600
  !
  ! prepare for the first iteration. For too small q the function could
  ! have noise.
  !
  qmax = 0.5_dp
  call sph_bes(7, grid%r(ik-3), qmax, lam, j1)
  j1(1:7) = j1(1:7)*grid%r(ik-3:ik+3)**flag
  logdermax = my_compute_log(j1, grid%r(ik), grid%dx) - logderae
  jlmax = j1(4)

  icount = 0
  do nc = 1,ncn
    !
    ! bracket the zero
    !
    icount1 = 0
    200 continue
    dq = dq_0
    qmin = qmax
    logdermin = logdermax
    jlmin = jlmax
    do iq = 1,imax
      !
      qmax = qmin + dq
      !
      call sph_bes(7, grid%r(ik-3), qmax, lam, j1)
      !
      j1(1:7) = j1(1:7)*grid%r(ik-3:ik+3)**flag
      !
      logdermax = my_compute_log(j1, grid%r(ik), grid%dx) - logderae
      !
      jlmax = j1(4)
      !
      ! the zero has been bracketed?
      !
      if( jlmin * jlmax > 0.d0 ) then ! far from an asintote
        !
        if( logdermax*logdermin < 0.d0 ) goto 100
        !
        qmin = qmax
        logdermin = logdermax
        jlmin = jlmax
      else
        if( logdermax*logdermin < 0.d0) then
          qmin = qmax
          logdermin=logdermax
          jlmin = jlmax
        else
          dq = 0.5d0 * dq
        endif
      endif
    enddo ! iq

    call infomsg ('find_qi','qmax not found ')
    iok = 1
    return

100 continue
    !
    !      start bisection loop
    !
    xc(nc) = qmin + (qmax-qmin)/2.d0
    !
    call sph_bes(7, grid%r(ik-3), xc(nc), lam, j1)
    j1(1:7) = j1(1:7)*grid%r(ik-3:ik+3)**flag
    !
    logder = my_compute_log(j1, grid%r(ik), grid%dx) - logderae
    !
    if( logder*logdermin < 0.d0) then
      qmax = xc(nc)
      logdermax = logder
    else
      qmin = xc(nc)
      logdermin = logder
    endif
    !
    ! avoid the asintotes
    !
    if( abs(logdermin-logdermax) > 1.e3_dp ) then 
      qmax = xc(nc)  
      logdermax = logder
      icount1 = icount1 + 1
      if( icount1 < 20) then
        goto 200
      else
        call errore('find_q','problem finding q',1)
      endif
    endif
    !
    ! check for convergence
    !
    icount = icount + 1
    if( icount > 1000 ) call errore('find_q','too many iterations',1)
    if( abs(logdermax-logdermin) > 1.e-8_dp) goto 100 ! bisection
  enddo

  return
end subroutine

!--------------------------------
function my_compute_log(j1,rj,dx)
!--------------------------------
  use kinds, only : DP
  implicit none
  real(DP) ::   &
       my_compute_log, &
       deriv_7pts,  &
       dx,          &
       j1(7),       &
       rj
  my_compute_log = deriv_7pts(j1, 4, rj, dx)/j1(4)

  return
end function
