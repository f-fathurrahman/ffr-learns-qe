!----------------------------
SUBROUTINE export_pspot_upf()
!----------------------------
  USE ions_base, ONLY: nsp
  USE uspp_param, ONLY: upf
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

    CALL json%create_object(inp, 'pseudo_upf')
    CALL json%add(p, inp)

    CALL json%add(inp, 'dft', trim(upf(isp)%dft))

    CALL json%add(inp, 'lmax', upf(isp)%lmax)
    CALL json%add(inp, 'lmax_rho', upf(isp)%lmax_rho)

    CALL json%add(inp, 'r', upf(isp)%r)
    CALL json%add(inp, 'rab', upf(isp)%rab)
    CALL json%add(inp, 'vloc', upf(isp)%vloc)

    CALL json%add(inp, 'shape_dion', shape(upf(isp)%dion))
    CALL json%add(inp, 'dion', reshape(upf(isp)%dion, [size(upf(isp)%dion)]))


    IF( allocated(upf(isp)%qfuncl) ) THEN
      CALL json%add(inp, 'shape_qfuncl', shape(upf(isp)%qfuncl))
      CALL json%add(inp, 'qfuncl', reshape(upf(isp)%qfuncl, [size(upf(isp)%qfuncl)]))
    ENDIF

    IF( allocated(upf(isp)%qfunc) ) THEN
      CALL json%add(inp, 'shape_qfunc', shape(upf(isp)%qfunc))
      CALL json%add(inp, 'qfunc', reshape(upf(isp)%qfunc, [size(upf(isp)%qfunc)]))
    ENDIF

    ! Set the filename
    WRITE(filename,"(A11,I1,A5)") 'pseudo_upf_', isp, '.json'

    ! write to file
    CALL json%print(p, trim(filename))
  
    ! Free the object
    CALL json%destroy(p)

    IF( json%failed() ) STOP 'Failed destroying JSON object'

  ENDDO

  RETURN

END SUBROUTINE