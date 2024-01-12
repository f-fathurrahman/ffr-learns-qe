!----------------------------------------------------------------------
subroutine my_dvloc_of_g(mesh, msh, rab, r, vloc_at, zp, tpiba2, ngl, gl, &
     omega, dvloc)
  !----------------------------------------------------------------------
  !
  ! dvloc = D Vloc (g^2) / D g^2 = (1/2g) * D Vloc(g) / D g
  !
  USE kinds
  USE constants , ONLY : pi, fpi, e2, eps8
  USE Coul_cut_2D, ONLY: do_cutoff_2D
  USE esm, ONLY : do_comp_esm, esm_bc
  implicit none
  !
  !    first the dummy variables
  !
  integer, intent(in) :: ngl, mesh, msh
  ! the number of shell of G vectors
  ! max number of mesh points
  ! number of mesh points for radial integration
  real(DP), intent(in) :: zp, rab(mesh), r(mesh), vloc_at(mesh), &
                          tpiba2, omega, gl(ngl)
  ! valence pseudocharge
  ! the derivative of the radial grid
  ! the radial grid
  ! the pseudo on the radial grid
  ! 2 pi / alat
  ! the volume of the unit cell
  ! the moduli of g vectors for each s
  !
  real(DP), intent(out) ::  dvloc(ngl)
  ! the fourier transform dVloc/dG
  !
  real(DP) :: vlcp, g2a, gx
  real(DP), allocatable :: aux(:), aux1(:)
  real(DP), external ::  qe_erf

  integer :: i, igl, igl0
  ! counter on erf functions or gaussians
  ! counter on g shells vectors
  ! first shell with g != 0

  ! the  G=0 component is not computed
  if( gl(1) < eps8 ) then
    dvloc(1) = 0.0d0
    igl0 = 2
  else
    igl0 = 1
  endif

  ! Pseudopotentials in numerical form (Vloc contains the local part)
  ! In order to perform the Fourier transform, a term erf(r)/r is
  ! subtracted in real space and added again in G space
  allocate(aux1(mesh))
  !
  ! This is the part of the integrand function
  ! indipendent of |G| in real space
  !
  do i = 1, msh
    aux1(i) = r(i)*vloc_at(i) + zp * e2 * qe_erf(r(i))
  enddo
  !
  allocate(aux(mesh))
  !
  do igl = igl0, ngl
    gx = sqrt(gl(igl)*tpiba2)
    !
    ! and here we perform the integral, after multiplying for the |G|
    ! dependent  part
    !
    ! DV(g)/Dg = Integral of r (Dj_0(gr)/Dg) V(r) dr
    do i = 1, msh
      aux(i) = aux1(i) * ( r(i) * cos(gx*r(i) ) / gx - sin(gx*r(i)) / gx**2 )
    enddo
    call simpson( msh, aux, rab, vlcp )
    ! DV(g^2)/Dg^2 = (DV(g)/Dg)/2g
    vlcp = fpi / omega / 2.0d0 / gx * vlcp
    ! for ESM stress
    ! In ESM, vloc and dvloc have only short term.
    IF( ( .not. do_comp_esm ) .or. ( esm_bc .eq. 'pbc' ) ) THEN
      ! subtract the long-range term
      IF( .not. do_cutoff_2D ) then ! 2D cutoff: do not re-add LR part here (re-added later in stres_loc)
       g2a = gl(igl) * tpiba2 / 4.d0
       vlcp = vlcp + fpi/omega * zp * e2 * exp(-g2a) * (g2a + 1.d0) / ( gl(igl) * tpiba2 )**2
      ENDIF
    ENDIF
    dvloc(igl) = vlcp
  enddo

  deallocate(aux)
  deallocate(aux1)

  return
end subroutine my_dvloc_of_g