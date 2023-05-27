! Using json-fortran module to export various data from PWSCF
! 
! Intended to be used for data that cannot be accessed easily from
! Julia, especially, the ones that defined using derived data types.
! Other varibles such as scalars and arrays that are exposed directly
! in the .so file should be accessed directly using cglobal and
! other C-interface functions of Julia.
!
! Current limitations of json-fortran:
!
! - Multidimensional arrays are not supported.
!   Work-around: using reshape to one-dimensional array and
!   also its shape
!
! - Complex scalars and arrays are not supported
!   Work-around: write real and imag parts separately


INCLUDE 'prepare_all.f90'

PROGRAM main

  CALL prepare_all()
  CALL export_atoms()
  CALL export_pspot_upf()
  CALL export_paw_radial_integrator()

END PROGRAM main


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
    WRITE(*,*) 'No PAW Data: early return'
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


