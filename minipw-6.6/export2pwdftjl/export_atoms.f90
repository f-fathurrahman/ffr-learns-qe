!------------------------
SUBROUTINE export_atoms()
!------------------------
  USE cell_base, ONLY: tpiba2, at, bg, omega
  USE ions_base, ONLY: nsp, tau, ityp, atm
  USE json_module
  !
  IMPLICIT NONE
  !
  ! JSON stuffs
  TYPE(json_core) :: json
  TYPE(json_value), pointer :: p, inp
  !
  CHARACTER(20) :: filename

  filename = 'atoms.json'


  CALL json%initialize()
  CALL json%create_object(p, '')

  CALL json%create_object(inp, 'atoms')
  CALL json%add(p, inp)

  CALL json%add(inp, 'tpiba2', tpiba2)
  CALL json%add(inp, 'omega', omega)

  CALL json%add(inp, 'shape_at', shape(at))
  CALL json%add(inp, 'at', reshape(at, [size(at)]))

  CALL json%add(inp, 'shape_bg', shape(bg))
  CALL json%add(inp, 'bg', reshape(bg, [size(bg)]))

  CALL json%add(inp, 'shape_tau', shape(tau))
  CALL json%add(inp, 'tau', reshape(tau, [size(tau)]))

  CALL json%add(inp, 'ityp', ityp)
  CALL json%add(inp, 'atm', atm(1:nsp))

  nullify(inp)

  ! write to file
  CALL json%print(p, trim(filename))

  CALL json%destroy(p)

  IF( json%failed() ) STOP 'Failed destroying JSON object'

  RETURN

END SUBROUTINE