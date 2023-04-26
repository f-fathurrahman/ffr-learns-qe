subroutine write_array1_r8(filename, n1, arr)
  implicit none
  character(*) :: filename
  integer :: n1, i
  real(8) :: arr(n1)
  integer :: iu

  write(*,*) 'write_array1_r8 to filename = ', trim(filename)

  iu = 1234
  open(unit=iu, file=trim(filename))
  do i = 1,n1
    write(iu,*) arr(i)
  enddo
  close(iu)

  return
end subroutine


subroutine write_array2_r8(filename, n1, n2, arr)
  implicit none
  character(*) :: filename
  integer :: n1, n2, i, j
  real(8) :: arr(n1,n2)
  integer :: iu

  write(*,*) 'write_array2_r8 to filename = ', trim(filename)

  iu = 1234
  open(unit=iu, file=trim(filename))
  do j = 1,n2
    do i = 1,n1
      write(iu,*) arr(i,j)
    enddo
  enddo
  close(iu)

  return
end subroutine