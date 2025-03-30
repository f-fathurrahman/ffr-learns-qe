subroutine my_vpack (ndim, ndimx, nspin, vin, vout, iflag)
  use kinds, ONLY: DP
  implicit none
  integer :: ndim, ndimx, nspin, iflag, n

  real(DP) :: vin (ndimx * nspin), vout (ndimx * nspin)
  if (nspin.eq.1.or.ndim.eq.ndimx) return
  if (iflag.eq.1) then
     do n = 1, ndim
        vin (n + ndim) = vin (n + ndimx)
        vout (n + ndim) = vout (n + ndimx)
     enddo
  elseif (iflag.eq. - 1) then
     do n = ndim, 1, - 1
        vin (n + ndimx) = vin (n + ndim)
        vout (n + ndimx) = vout (n + ndim)
     enddo
     do n = ndim + 1, ndimx
        vin (n) = 0.0_DP
        vout (n) = 0.0_DP
     enddo
  else
     call errore ('vpack', ' wrong flag ', 1)

  endif
  return
end subroutine
