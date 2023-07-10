!--------------------------------
SUBROUTINE export_paw_variables()
!--------------------------------
  USE ions_base, ONLY: nsp
  USE paw_variables
  !
  USE json_module
  
  IMPLICIT NONE 
  INTEGER :: isp
  CHARACTER(30) :: filename
  ! JSON stuffs
  TYPE(json_core) :: json
  TYPE(json_value), POINTER :: p, inp

  IF( .not. okpaw ) THEN 
    WRITE(*,*) 'paw_variables: early return'
    RETURN 
  ENDIF 

  CALL json%initialize()

  CALL json%create_object(p, '')

  CALL json%create_object(inp, 'paw_variables')
  CALL json%add(p, inp)


  CALL json%add(inp, 'shape_ddd_paw', shape(ddd_paw))
  CALL json%add(inp, 'ddd_paw', reshape(ddd_paw, [size(ddd_paw)]))
 
  ! Set the filename
  WRITE(filename,"(A18,I1,A5)") 'paw_variables.json'

  ! write to file
  CALL json%print(p, trim(filename))
  
  ! Free the object
  CALL json%destroy(p)

  IF( json%failed() ) STOP 'Failed destroying JSON object'

  RETURN 

END SUBROUTINE