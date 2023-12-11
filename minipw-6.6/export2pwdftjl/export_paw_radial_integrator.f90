!----------------------------------------
SUBROUTINE export_paw_radial_integrator()
!----------------------------------------
  USE ions_base, ONLY: nsp
  USE paw_variables, ONLY : okpaw, rad
  USE uspp_param, ONLY: upf
  !
  USE json_module
  
  IMPLICIT NONE 
  INTEGER :: isp
  CHARACTER(30) :: filename
  ! JSON stuffs
  TYPE(json_core) :: json
  TYPE(json_value), POINTER :: p, inp

  IF( .not. okpaw ) THEN 
    WRITE(*,*) 'paw_radial_integrator: early return'
    RETURN 
  ENDIF 

  CALL json%initialize()

  DO isp = 1,nsp

    ! Skip if no PAW data for this species
    IF( .not. upf(isp)%tpawp ) CYCLE

    CALL json%create_object(p, '')

    CALL json%create_object(inp, 'paw_radial_integrator')
    CALL json%add(p, inp)

    CALL json%add(inp, 'lmax', rad(isp)%lmax)
    CALL json%add(inp, 'ladd', rad(isp)%ladd)
    CALL json%add(inp, 'lm_max', rad(isp)%lm_max)
 
    CALL json%add(inp, 'nx', rad(isp)%nx)

    CALL json%add(inp, 'ww', rad(isp)%ww)

    CALL json%add(inp, 'shape_ylm', shape(rad(isp)%ylm))
    CALL json%add(inp, 'ylm', reshape(rad(isp)%ylm, [size(rad(isp)%ylm)]))
 
    CALL json%add(inp, 'shape_wwylm', shape(rad(isp)%wwylm))
    CALL json%add(inp, 'wwylm', reshape(rad(isp)%wwylm, [size(rad(isp)%wwylm)]))
 
    if( associated(rad(isp)%dylmt) ) then
      CALL json%add(inp, 'shape_dylmt', shape(rad(isp)%dylmt))
      CALL json%add(inp, 'dylmt', reshape(rad(isp)%dylmt, [size(rad(isp)%dylmt)]))
    endif

    if( associated(rad(isp)%dylmp) ) then
      CALL json%add(inp, 'shape_dylmp', shape(rad(isp)%dylmp))
      CALL json%add(inp, 'dylmp', reshape(rad(isp)%dylmp, [size(rad(isp)%dylmp)]))
    endif

    if( associated(rad(isp)%cotg_th) ) then
      CALL json%add(inp, 'shape_cotg_th', shape(rad(isp)%cotg_th))
      CALL json%add(inp, 'cotg_th', reshape(rad(isp)%cotg_th, [size(rad(isp)%cotg_th)]))
    endif

    ! Set the filename
    WRITE(filename,"(A22,I1,A5)") 'paw_radial_integrator_', isp, '.json'

    ! write to file
    CALL json%print(p, trim(filename))
  
    ! Free the object
    CALL json%destroy(p)

    IF( json%failed() ) STOP 'Failed destroying JSON object'

  ENDDO 

  RETURN 

END SUBROUTINE