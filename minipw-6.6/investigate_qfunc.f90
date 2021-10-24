subroutine investigate_qfunc()
  use ions_base, only: nat, ntyp => nsp, ityp, tau
  use uspp_param, only: upf, lmaxq, nh, nhm, lmaxkb
  use uspp, only: indv_ijkb0
  implicit none
  integer :: nt

  do nt = 1, ntyp
    write(*,*)
    write(*,*) 'nt = ', nt
    if(.not. upf(nt)%q_with_l) then
      write(*,*) 'qfunc = ', upf(nt)%qfunc(5,1)
    endif
    write(*,*) 'q_with_l = ', upf(nt)%q_with_l
    write(*,*) 'qfuncl = ', upf(nt)%qfuncl(10,1,:)
    write(*,*) 'shape qfuncl = ', shape(upf(nt)%qfuncl)
    write(*,*) upf(nt)%qfuncl(1:5,1,0)
    !write(*,*) upf(nt)%qfuncl(1,1,1)
  enddo

  write(*,*) 'nhm = ', nhm
  write(*,*) 'nh = ', nh
  write(*,*) 'indv_ijkb0 = ', indv_ijkb0

end