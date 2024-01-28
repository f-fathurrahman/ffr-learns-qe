!=========================================================================
!
! Routines:
!
! (1) write_program_header Originally by DAS Last Modified 10/15/2010 (DAS)
!
! Write a header for the beginning of each program.
!
!=========================================================================

module write_program_header_m
  use global_m
  use version_m
  implicit none
  private
  public :: &
    write_program_header
contains
subroutine write_program_header(aname, write_to_7)
  character(len=*), intent(in) :: aname
  logical, intent(in) :: write_to_7
  character :: ainfo*256, aversion*32, adate*11, atime*14
  character(len=256) :: flags_str
 
  aversion = ' code, ' + trim("Real") + ' version, run on '
  call get_version_info(ainfo)
  call date_time(adate,atime)
  if (peinf%inode==0) then
    write(6,10) trim(ainfo)
    if (write_to_7) write(7,10) trim(ainfo)
 10 format(/' BerkeleyGW ',a)
    write(6,20) trim(aname), trim(aversion), trim(adate), ' at ', trim(atime)
    if (write_to_7) write(7,20) trim(aname), trim(aversion), trim(adate), ' at ', trim(atime)
 20 format(1x,2a,1x,3a)
    write(6,'()')
    call write_logo()
    write(6,'(1x,a)') '--------------------------------------------------------------------------------'
    write(6,'(1x,a)') 'Please cite the main BerkeleyGW paper when using results from this calculation:'
    write(6,'()')
    write(6,'(1x,a)') 'Jack Deslippe, Georgy Samsonidze, David A. Strubbe, Manish Jain, Marvin L.'
    write(6,'(1x,a)') 'Cohen, and Steven G. Louie, "BerkeleyGW: A Massively Parallel Computer Package'
    write(6,'(1x,a)') 'for the Calculation of the Quasiparticle and Optical Properties of Materials'
    write(6,'(1x,a)') 'and Nanostructures," Comput. Phys. Commun. 183, 1269 (2012)'
    write(6,'()')
    write(6,'(1X,a)') 'Aditionally, PLEASE CITE ALL PAPERS REPORTED AT THE END OF THIS CALCULATION.'
    write(6,'(1x,a)') '--------------------------------------------------------------------------------'
    write(6,30)
    if (write_to_7) write(7,30)
 30 format(1x,'Running serial version (no MPI)')
    write(6,'()')
    write(6,'(1x,a)') 'Compilation flags:'
    write(6,'(1x,2a)') '- Compiler: ', "GNU"
    flags_str = ''
    call write_flag("MPI", "MPI")
    call write_flag("OMP", "OMP")
    write(6,'(1x,2a)') '- Para. flags: ', trim(flags_str)
    flags_str = ''
    call write_flag("USESCALAPACK", "USESCALAPACK")
    call write_flag("UNPACKED", "UNPACKED")
    call write_flag("USEFFTW3", "1")
    call write_flag("HDF5", "HDF5")
    write(6,'(1x,2a)') '- Math  flags: ', trim(flags_str)
    flags_str = ''
    call write_flag("DEBUG", "DEBUG")
    call write_flag("VERBOSE", "VERBOSE")
    write(6,'(1x,2a)') '- Debug flags: ', trim(flags_str)
    !call write_flag("X", "X")
    write(6,'()')
  endif
 
  return
contains
  subroutine write_flag(flag_name, flag_tostring)
    character(len=*), intent(in) :: flag_name, flag_tostring
   
    if (flag_name/=flag_tostring) then
      if (len(trim(flags_str))>0) then
        flags_str = trim(flags_str) + ', ' + flag_name
      else
        flags_str = flag_name
      endif
    endif
   
  end subroutine write_flag
  subroutine write_berkeleygw()
   
   
  end subroutine write_berkeleygw
  subroutine write_logo()
   
    ! FHJ: Logo generated with jp2a utility. See LOGO/donkey2ascii.sh.
    ! FHJ: Banner was manually generated. It is based on the output of Figlet with font "Varsity"
    write(6,'(a)') '                                                                 ..o.          '
    write(6,'(a)') '                                                                .oxxo.         '
    write(6,'(a)') '                                                               .oxxxxo...      '
    write(6,'(a)') '                                                               oxxxxxxo.       '
    write(6,'(a)') '                                                              .oxxxxxxx.       '
    write(6,'(a)') '                                                              .ooooooxxo..     '
    write(6,'(a)') '                                                              .oooooooxo..     '
    write(6,'(a)') '                                                              .oooooxxo...     '
    write(6,'(a)') '                                                       .........oxooo......    '
    write(6,'(a)') '                                                 ............................  '
    write(6,'(a)') '                                             ................................. '
    write(6,'(a)') '                                          .................................... '
    write(6,'(a)') '           .          ..oo. ....  .................................oooxxxxxxxo.'
    write(6,'(a)') '     .............oxxxx@ox@@@x@x.....................o...........ooooooooooxx. '
    write(6,'(a)') '    .o.........oox@x.oo........xxx@@............ooxxxxo..........ooooxxxxxoxo  '
    write(6,'(a)') '    .x........x@xxo...............o@xxo........oxxx@@@xoooooooooooooooxxxo...  '
    write(6,'(a)') '    .o......ox@@o..................oox@o.....ooxxx@xoooxxxxxxxoooooooooooo.... '
    write(6,'(a)') '    o..ooooo@@xoooo....ooo...........x@o.....ooxxxxo   .oxxxxxxxxxxooooooo.... '
    write(6,'(a)') '    . .oooo@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....ooooox.     .oxx@@@@@xxoo........ '
    write(6,'(a)') '      .ooooxxxxxxxooooooxooooooooooooooxo...oooooxx.        ..ox@xxxoo.........'
    write(6,'(a)') '      .ooooooxxxx@xoooooooooooooxxoooooooooooxxxxx.            .oxxooooooooxxo.'
    write(6,'(a)') '     .oooooxxxxx@@@xxxxxxxxxxxxxxxxxxxxxxxxoxxxxo.                .oxxxxxxxxo. '
    write(6,'(a)') '   ....oooxxxxx@@@xo..oxxx@@@@@@xxxxoxxoooooooxx.                  .oxxxoo..   '
    write(6,'(a)') '  .....ooxxxx@@xo.       ........    .ooooooooxxo                              '
    write(6,'(a)') '  ..oooxxxx@@@o                       .oooooooxoo.                             '
    write(6,'(a)') '  ....oxooxxxxx.                       .ooo..oooo.                             '
    write(6,'(a)') '  .....o.ooxxxxxo.                     .oooooooxo.                             '
    write(6,'(a)') ' ......ooooxxxxxxo.                     .ooooooxoo..                           '
    write(6,'(a)') '........ooxxxxxxxo..                     .o....oxoo...                         '
    write(6,'(a)') '.......ooooxxxxxxxo.                     ........oooo.                         '
    write(6,'(a)') '.ooooooo..ooxxxxoooo.                    .........ooo...                       '
    write(6,'(a)') '..oxo...ooooxxxoooo..                    .ooo......oooo...                     '
    write(6,'(a)') '  .ooooo....o.                            .oxxxoo....ooo....                   '
    write(6,'(a)') '    .oooooo...                              ...ooooo...ooo..                   '
    write(6,'(a)') '       ...                                     .oo.......                      '
    write(6,'(a)') '                                               ....ooo...                      '
    write(6,'(a)') "                      __              __                                       "
    write(6,'(a)') " ______              [  |            [  |                 ._____  _          _ "
    write(6,'(a)') "|_   _ \              | |  _          | |                /  ___ \| |        | |"
    write(6,'(a)') "  | |_) | .---.  _. _.| | / |  .---.  | | .---.  _    _ / /    \_|\ \  /\  / / "
    write(6,'(a)') "  |  __'./ /__\\[ /`\_| '' <  / /__\\ | |/ /__\\| \  | | |   _____ \ \/  \/ /  "
    write(6,'(a)') " _| |__| | \__. | |   | |`\ \ | \___. | || \___. \ \/ / \ \.___| |  \  /\  /   "
    write(6,'(a)') "|_______/ \.__./[_]  [__|  \_] \.__./[___]\.__./  \  /   \.____./    \/  \/    "
    write(6,'(a)') "                                                  / /                          "
    write(6,'(a)') "                                                 /_/                           "
   
  end subroutine write_logo
end subroutine write_program_header
end module write_program_header_m
