subroutine ld1x_print_variables()
  use ld1inc
  implicit none
  integer :: n

  write(*,*)

  write(*,*) 'nwfs = ', nwfs ! no of pseudowavefunc
  write(*,*) 'nwf = ', nwf

  write(*,'(1x,A)', advance='no') 'el (wavefunction labels) = '
  do n = 1,nwf
    write(*,'(1x,2A,A)', advance='no') el(n), ' '
  enddo
  write(*,*)

  write(*,'(1x,A)', advance='no') 'nn = '
  do n = 1,nwf
    write(*,'(1x,I2,A)',advance='no') nn(n), ' '
  enddo
  write(*,*)

  write(*,'(1x,A)', advance='no') 'll = '
  do n = 1,nwf
    write(*,'(1x,I2,A)',advance='no') ll(n), ' '
  enddo
  write(*,*)


  write(*,'(1x,A)', advance='no') 'isw = '
  do n = 1,nwfs
    write(*,'(1x,I2,A)',advance='no') isw(n), ' '
  enddo
  write(*,*)


  write(*,'(1x,A)',advance='no') 'core_state = '
  do n = 1,nwf
    write(*,'(1x,1L,A)',advance='no') core_state(n), ' '
  enddo
  write(*,*)

  write(*,'(1x,A)', advance='no') 'nstoae = '
  do n = 1,nwfs
    write(*,'(1x,I2,A)',advance='no') nstoae(n), ' '
  enddo
  write(*,*)

  write(*,'(1x,A)', advance='no') 'els (pseudo wavefunction labels) = '
  do n = 1,nwfs
    write(*,'(1x,2A,A)', advance='no') els(n), ' '
  enddo
  write(*,*)

  write(*,'(1x,A)', advance='no') 'nns = '
  do n = 1,nwfs
    write(*,'(1x,I2,A)',advance='no') nns(n), ' '
  enddo
  write(*,*)

  write(*,'(1x,A)', advance='no') 'lls = '
  do n = 1,nwfs
    write(*,'(1x,I2,A)',advance='no') lls(n), ' '
  enddo
  write(*,*)


  write(*,'(1x,A)', advance='no') 'isws = '
  do n = 1,nwfs
    write(*,'(1x,I2,A)',advance='no') isws(n), ' '
  enddo
  write(*,*)

  write(*,*) 'lloc = ', lloc
  write(*,*) 'nbeta = ', nbeta
  write(*,*) 'nsloc = ', nsloc



end subroutine