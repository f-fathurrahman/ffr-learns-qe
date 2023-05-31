!-------------------------
SUBROUTINE export_vlocal()
!-------------------------
  USE vlocal
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

  CALL json%create_object(inp, 'vlocal')
  CALL json%add(p, inp)

  CALL json%add(inp, 'shape_vloc', shape(vloc))
  CALL json%add(inp, 'vloc', reshape(vloc, [size(vloc)]))

  CALL json%add(inp, 'shape_strf', shape(strf))
  CALL json%add(inp, 'strf_real', reshape(real(strf), [size(strf)]))
  CALL json%add(inp, 'strf_imag', reshape(imag(strf), [size(strf)]))

  ! Set the filename
  filename = 'vlocal_mod.json'

  ! write to file
  CALL json%print(p, trim(filename))
  
  ! Free the object
  CALL json%destroy(p)

  IF( json%failed() ) STOP 'Failed destroying JSON object'

  RETURN

END SUBROUTINE
