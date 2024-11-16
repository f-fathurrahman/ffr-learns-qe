
! ffr: use global variables directly

subroutine ld1x_driver_starting_potential()
  USE radial_grids, ONLY: ndmx
  USE ld1inc, ONLY: grid, enne, vpot, vxt, enl, &
                    nspin, nn, ll, oc, nwf, &
                     zed, zval, v0

  CALL starting_potential( ndmx, grid%mesh, zval, zed, nwf, oc, nn, ll,&
                           grid%r, enl, v0, vxt, vpot, enne, nspin )

  return
end subroutine
