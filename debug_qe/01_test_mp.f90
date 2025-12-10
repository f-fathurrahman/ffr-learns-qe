program test_mp
  use io_global, only: ionode
  use mp_global, only: mp_startup, mp_global_end
  !USE environment, only : environment_start, environment_end
  implicit none

  !call mp_startup(start_images=.true., images_only=.true.)
  call mp_startup()
  !call environment_start('my_env')

  if(ionode) then
    write(*,*) 'Pass here'
  endif

  !CALL environment_end('my_env')
  call mp_global_end()
  stop

end program
