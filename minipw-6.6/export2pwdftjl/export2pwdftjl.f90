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
  call export_atoms()
  CALL export_pspot_upf()

END PROGRAM main


subroutine export_atoms()
  use cell_base
  use ions_base
  use json_module
  !
  implicit none
  !
  ! JSON stuffs
  type(json_core) :: json
  type(json_value), pointer :: p, inp
  !
  character(20) :: filename

  filename = 'atoms.json'


  call json%initialize()
  call json%create_object(p, '')

  call json%create_object(inp, 'atoms')
  call json%add(p, inp)

  call json%add(inp, 'tpiba2', tpiba2)
  call json%add(inp, 'omega', omega)

  call json%add(inp, 'shape_at', shape(at))
  call json%add(inp, 'at', reshape(at, [size(at)]))

  call json%add(inp, 'shape_bg', shape(bg))
  call json%add(inp, 'bg', reshape(bg, [size(bg)]))

  call json%add(inp, 'shape_tau', shape(tau))
  call json%add(inp, 'tau', reshape(tau, [size(tau)]))

  call json%add(inp, 'ityp', ityp)
  call json%add(inp, 'atm', atm)

  nullify(inp)

  ! write to file
  call json%print(p, trim(filename))

  call json%destroy(p)

  if(json%failed()) stop 'Failed destroying JSON object' 

end subroutine




!----------------------------
subroutine export_pspot_upf()
!----------------------------
  use ions_base, only: nsp
  use uspp_param, only: upf
  !
  use json_module
  !
  implicit none
  !
  integer :: isp
  character(20) :: filename
  ! JSON stuffs
  type(json_core) :: json
  type(json_value), pointer :: p, inp


  isp = 1


  call json%initialize()
  call json%create_object(p, '')

  call json%create_object(inp, 'pseudo_upf')
  call json%add(p, inp)

  call json%add(inp, 'dft', trim(upf(isp)%dft))

  call json%add(inp, 'lmax', upf(isp)%lmax)
  call json%add(inp, 'lmax_rho', upf(isp)%lmax_rho)


  call json%add(inp, 'shape_dion', shape(upf(isp)%dion))
  call json%add(inp, 'dion', reshape(upf(isp)%dion, [size(upf(isp)%dion)]))

  if(allocated(upf(isp)%qfuncl)) then
    call json%add(inp, 'shape_qfuncl', shape(upf(isp)%qfuncl))
    call json%add(inp, 'qfuncl', reshape(upf(isp)%qfuncl, [size(upf(isp)%qfuncl)]))
  endif

  if(allocated(upf(isp)%qfunc)) then
    call json%add(inp, 'shape_qfunc', shape(upf(isp)%qfunc))
    call json%add(inp, 'qfunc', reshape(upf(isp)%qfunc, [size(upf(isp)%qfunc)]))
  endif


  write(filename,"(A11,I1,A5)") 'pseudo_upf_', isp, '.json'

  ! write to file
  call json%print(p, trim(filename))

  call json%destroy(p)

  if(json%failed()) stop 'Failed destroying JSON object' 

end subroutine