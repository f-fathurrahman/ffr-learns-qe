!
! Copyright (C) 2004-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
function int_0_inf_dr(f,grid,mesh,nst)
  !---------------------------------------------------------------
  !
  ! integral of f from 0 to infinity
  ! f is given on a logarithmic mesh. 
  ! f(r) is assumed to be proportional to r**nst for small r
  !
  USE kinds, only : DP
  USE radial_grids, only: radial_grid_type, series
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER, INTENT(in) :: mesh, nst
  REAL(DP), INTENT(in):: f(mesh)
  TYPE(radial_grid_type), INTENT(in) :: grid
  REAL(DP) :: int_0_inf_dr
  !
  ! local variables
  !
  REAL(DP):: fs(4), b(4), sum1
  INTEGER :: i
  !
  ! series development: contribution for small r
  !
  IF( mesh > grid%mesh ) &
     CALL errore('int_0_inf_dr','value of mesh is larger than expected',mesh)

  DO i = 1,4
    fs(i) = f(i)/grid%r(i)**nst
  ENDDO 
  CALL series( fs, grid%r, grid%r2, b )
  int_0_inf_dr = ( b(1)/(nst+1) + grid%r(1)* &
              &  ( b(2)/(nst+2) + grid%r(1)*b(3)/(nst+3)) ) * grid%r(1)**(nst+1)
  !
  ! simpson integration (logarithmic mesh: dr ==> r dx)
  !
  sum1 = 0.0_DP
  DO i = 1,mesh-2,2
     sum1 = sum1 + f(i)*grid%r(i) + 4.0_dp*f(i+1)*grid%r(i+1) + f(i+2)*grid%r(i+2)
  ENDDO 
  int_0_inf_dr = int_0_inf_dr + sum1*grid%dx/3.0_DP
  RETURN
END FUNCTION int_0_inf_dr

