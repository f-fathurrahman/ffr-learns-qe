!-------------------------------
subroutine export_vars_scf_mod()
!-------------------------------
  use scf
  !
  USE json_module
  !
  IMPLICIT NONE
  !
  CHARACTER(20) :: filename
  ! JSON stuffs
  TYPE(json_core) :: json
  TYPE(json_value), POINTER :: p, inp

  CALL json%initialize()

  CALL json%create_object(p, '')

  CALL json%create_object(inp, 'scf_mod')
  CALL json%add(p, inp)

  CALL json%add(inp, 'v_of_0', v_of_0)

  if( allocated(rho_core) ) then
    CALL json%add(inp, 'rho_core', rho_core)
  endif

  if( allocated(vltot) ) then
    CALL json%add(inp, 'vltot', vltot)
  endif

  if( allocated(vrs) ) then
    CALL json%add(inp, 'shape_vrs', shape(vrs))
    CALL json%add(inp, 'vrs', reshape(vrs, [size(vrs)]))
  endif

  if( allocated(v%of_r) ) then
    CALL json%add(inp, 'shape_v_of_r', shape(v%of_r))
    CALL json%add(inp, 'v_of_r', reshape(v%of_r, [size(v%of_r)]))
  endif

  if( allocated(rho%of_r) ) then
    CALL json%add(inp, 'shape_rho_of_r', shape(rho%of_r))
    CALL json%add(inp, 'rho_of_r', reshape(rho%of_r, [size(rho%of_r)]))
  endif


  ! Set the filename
  filename = 'scf_mod.json'

  ! write to file
  CALL json%print(p, trim(filename))
  
  ! Free the object
  CALL json%destroy(p)

  IF( json%failed() ) STOP 'Failed destroying JSON object'

  RETURN


end subroutine