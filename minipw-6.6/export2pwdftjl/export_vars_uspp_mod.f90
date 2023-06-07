!--------------------------------
SUBROUTINE export_vars_uspp_mod()
!--------------------------------
  USE uspp
  USE us  ! in file pwcom.f90
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

  CALL json%create_object(inp, 'uspp')
  CALL json%add(p, inp)

  if( okvan ) then
    CALL json%add(inp, 'shape_lpx', shape(lpx))
    CALL json%add(inp, 'lpx', reshape(lpx, [size(lpx)]))

    CALL json%add(inp, 'shape_lpl', shape(lpl))
    CALL json%add(inp, 'lpl', reshape(lpl, [size(lpl)]))

    CALL json%add(inp, 'shape_ap', shape(ap))
    CALL json%add(inp, 'ap', reshape(ap, [size(ap)]))
  endif

  CALL json%add(inp, 'nkb', nkb)
  CALL json%add(inp, 'nkbus', nkbus)

  if( allocated(indv) ) then
    CALL json%add(inp, 'shape_indv', shape(indv))
    CALL json%add(inp, 'indv', reshape(indv, [size(indv)]))
  endif

  if( allocated(nhtol) ) then
    CALL json%add(inp, 'shape_nhtol', shape(nhtol))
    CALL json%add(inp, 'nhtol', reshape(nhtol, [size(nhtol)]))
  endif

  if( allocated(nhtolm) ) then
    CALL json%add(inp, 'shape_nhtolm', shape(nhtolm))
    CALL json%add(inp, 'nhtolm', reshape(nhtolm, [size(nhtolm)]))
  endif

  if( allocated(ijtoh) ) then
    CALL json%add(inp, 'shape_ijtoh', shape(ijtoh))
    CALL json%add(inp, 'ijtoh', reshape(ijtoh, [size(ijtoh)]))
  endif

  if( allocated(indv_ijkb0) ) then
    CALL json%add(inp, 'shape_indv_ijkb0', shape(indv_ijkb0))
    CALL json%add(inp, 'indv_ijkb0', reshape(indv_ijkb0, [size(indv_ijkb0)]))
  endif

  if( allocated(becsum) ) then
    CALL json%add(inp, 'shape_becsum', shape(becsum))
    CALL json%add(inp, 'becsum', reshape(becsum, [size(becsum)]))
  endif


  if( allocated(vkb) ) then
    CALL json%add(inp, 'shape_vkb', shape(vkb))
    CALL json%add(inp, 'vkb_real', reshape(real(vkb), [size(vkb)]))
    CALL json%add(inp, 'vkb_imag', reshape(imag(vkb), [size(vkb)]))
  endif


  if( allocated(qq_nt) ) then
    CALL json%add(inp, 'shape_qq_nt', shape(qq_nt))
    CALL json%add(inp, 'qq_nt', reshape(qq_nt, [size(qq_nt)]))
  endif

  if( allocated(qq_at) ) then
    CALL json%add(inp, 'shape_qq_at', shape(qq_at))
    CALL json%add(inp, 'qq_at', reshape(qq_at, [size(qq_at)]))
  endif


  if( allocated(dvan) ) then
    CALL json%add(inp, 'shape_dvan', shape(dvan))
    CALL json%add(inp, 'dvan', reshape(dvan, [size(dvan)]))
  endif


  if( allocated(deeq) ) then
    CALL json%add(inp, 'shape_deeq', shape(deeq))
    CALL json%add(inp, 'deeq', reshape(deeq, [size(deeq)]))
  endif



  ! From file pwcom.f90
  if( allocated(qrad) ) then
    CALL json%add(inp, 'shape_qrad', shape(qrad))
    CALL json%add(inp, 'qrad', reshape(qrad, [size(qrad)]))
  endif

  ! Set the filename
  filename = 'uspp_mod.json'

  ! write to file
  CALL json%print(p, trim(filename))
  
  ! Free the object
  CALL json%destroy(p)

  IF( json%failed() ) STOP 'Failed destroying JSON object'

  RETURN

end subroutine