subroutine my_vpack( ndim, ndimx, nspin, vin, vout, iflag )
  use kinds, ONLY: DP
  !
  implicit none
  !
  integer :: ndim, ndimx, nspin, iflag, n
  real(DP) :: vin(ndimx * nspin), vout(ndimx * nspin)
  
  ! only useful for nspin > 1
  if( nspin == 1 .or. ndim == ndimx ) then
    return
  endif

  if( iflag == 1 ) then
    do n = 1, ndim
      vin(n + ndim) = vin(n + ndimx)
      vout(n + ndim) = vout(n + ndimx)
    enddo
  elseif( iflag ==  - 1 ) then
    do n = ndim,1,-1
      vin(n + ndimx) = vin(n + ndim)
      vout(n + ndimx) = vout(n + ndim)
    enddo
    do n = ndim + 1, ndimx
      vin(n) = 0.d0
      vout(n) = 0.d0
    enddo
  else
     call errore('vpack', ' wrong flag ', 1)
  endif
  return
end subroutine
