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

  if(nsp /= 2) then
    stop 'nsp is hardcoded to 2'
  endif

  call write_array2_r8('wwylm1.dat', size(rad(1)%wwylm, 1), size(rad(1)%wwylm, 2), rad(1)%wwylm)
  call write_array1_r8('ww1.dat', size(rad(1)%ww, 1), rad(1)%ww)

  call write_array2_r8('wwylm2.dat', size(rad(2)%wwylm, 1), size(rad(2)%wwylm, 2), rad(2)%wwylm)
  call write_array1_r8('ww2.dat', size(rad(2)%ww, 1), rad(2)%ww)

  return

end subroutine