INCLUDE 'prepare_all.f90'


PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  CALL do_write_arrays()
END PROGRAM


subroutine do_write_arrays()
  use paw_variables, only : rad
  USE ions_base, ONLY : nat, nsp
  implicit none
  integer :: isp
  character(30) :: filename

  if(nsp > 9) then
    stop 'nsp is harcoded to be not larger than 9'
  endif

  do isp = 1,nsp
    ! wwylm
    write(filename,"(A5,I1,A4)") 'wwylm', isp, '.dat'
    call write_array2_r8(filename, size(rad(1)%wwylm, 1), size(rad(1)%wwylm, 2), rad(1)%wwylm)
    !
    write(filename,"(A2,I1,A4)") 'ww', isp, '.dat'
    call write_array1_r8(filename, size(rad(1)%ww, 1), rad(1)%ww)
  enddo

  return

end subroutine