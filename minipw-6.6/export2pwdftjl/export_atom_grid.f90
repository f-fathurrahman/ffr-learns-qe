!----------------------------
SUBROUTINE export_atom_grid()
!----------------------------
  USE ions_base, ONLY: nsp
  USE atom
  !
  USE json_module
  !
  IMPLICIT NONE
  !
  INTEGER :: isp
  CHARACTER(20) :: filename
  ! JSON stuffs
  TYPE(json_core) :: json
  TYPE(json_value), POINTER :: p, inp

  CALL json%initialize()

  DO isp = 1,nsp

    CALL json%create_object(p, '')

    CALL json%create_object(inp, 'atom_grid')
    CALL json%add(p, inp)

    CALL json%add(inp, 'msh', msh(isp))

    CALL json%add(inp, 'r', rgrid(isp)%r)
    CALL json%add(inp, 'rab', rgrid(isp)%rab)

    ! Set the filename
    WRITE(filename,"(A10,I1,A5)") 'atom_grid_', isp, '.json'

    ! write to file
    CALL json%print(p, trim(filename))
  
    ! Free the object
    CALL json%destroy(p)

    IF( json%failed() ) STOP 'Failed destroying JSON object'

  ENDDO

  RETURN

END SUBROUTINE