!----------------------------------------------------------------------
SUBROUTINE my_vloc_of_g( mesh, msh, rab, r, vloc_at, zp, tpiba2, ngl, &
                      gl, omega, vloc )
  !----------------------------------------------------------------------
  !! This routine computes the Fourier transform of the local
  !! part of an atomic pseudopotential, given in numerical form.
  !! A term erf(r)/r is subtracted in real space (thus making the
  !! function short-ramged) and added again in G space (for G<>0)
  !! The G=0 term contains \int (V_loc(r)+ Ze^2/r) 4pi r^2 dr.
  !! This is the "alpha" in the so-called "alpha Z" term of the energy.
  !! Atomic Ry units everywhere.
  !
  USE kinds
  USE constants, ONLY : pi, fpi, e2, eps8
  USE Coul_cut_2D, ONLY: do_cutoff_2D
  !
  IMPLICIT NONE
  !
  integer, intent(in) :: ngl
  !! the number of shells of G vectors
  INTEGER, INTENT(IN) :: mesh
  !! number of grid points in the radial grid
  INTEGER, INTENT(IN) :: msh
  !! as above, used for radial integration
  REAL(DP), INTENT(IN) :: zp
  !! valence pseudocharge
  REAL(DP), INTENT(IN) :: rab(mesh)
  !! the derivative of mesh points
  REAL(DP), INTENT(IN) :: r(mesh)
  !! the mesh points
  REAL(DP), INTENT(IN) :: vloc_at(mesh)
  !! local part of the atomic pseudopotential on the radial mesh
  REAL(DP), INTENT(IN) :: tpiba2
  !! 2 pi / alat
  REAL(DP), INTENT(IN) :: omega
  !! the volume of the unit cell
  REAL(DP), INTENT(IN) :: gl(ngl)
  !! the moduli of g vectors for each shell 
  REAL(DP), INTENT(OUT) :: vloc(ngl)
  !! the fourier transform of the potential
  !
  ! ... local variables
  !
  REAL(DP) :: vlcp, fac, gx
  REAL(DP), ALLOCATABLE :: aux(:), aux1(:)
  INTEGER :: igl, igl0, ir
  ! igl :counter on g shells vectors
  ! igl0:first shell with g != 0
  ! ir  :counter on mesh points
  !
  REAL(DP), EXTERNAL :: qe_erf
  !
  allocate ( aux(msh), aux1(msh) )
  if( gl(1) < eps8 ) then
    !
    ! first the G=0 term
    !
    ! do_comp_esp and do_cutoff_2D are removed
    !
    do ir = 1, msh
      aux(ir) = r(ir) * ( r(ir) * vloc_at(ir) + zp * e2 )
    enddo
    call simpson( msh, aux, rab, vlcp )
    vloc(1) = vlcp        
    igl0 = 2
  else
    igl0 = 1
  endif
  !
  ! here the G<>0 terms, we first compute the part of the integrand 
  ! function independent of |G| in real space
  !
  do ir = 1, msh
    aux1(ir) = r(ir) * vloc_at(ir) + zp * e2 * qe_erf(r(ir))
  enddo
  fac = zp * e2 / tpiba2
  !
  ! and here we perform the integral, after multiplying for the |G|
  ! dependent part
  !
  do igl = igl0, ngl
    gx = sqrt(gl(igl) * tpiba2)
    do ir = 1, msh
      aux(ir) = aux1(ir) * sin(gx*r(ir))/gx
    enddo
    call simpson( msh, aux, rab, vlcp )
    !
    ! here we re-add the analytic fourier transform of the erf function
    !
    IF( .not. do_cutoff_2D ) THEN
      vlcp = vlcp - fac*exp(-gl(igl)*tpiba2*0.25d0) / gl (igl)
    ENDIF
    vloc(igl) = vlcp
  enddo
  vloc(:) = vloc(:)*fpi/omega
  deallocate (aux, aux1)

  return
END SUBROUTINE
