subroutine info_upf()

  USE ions_base,  ONLY : ntyp => nsp
  USE uspp_param, ONLY: upf

  implicit none
  integer :: nt

  write(*,*)
  write(*,*) '      info_upf:'
  write(*,*) '      ---------'

  do nt = 1,ntyp
    write(*,*) 'upf ', nt, ' multiproj = ', upf(nt)%is_multiproj
  enddo

  write(*,*)

  return

end subroutine