!this subroutines calls solvers for eigen-value problem of the excitions
SUBROUTINE simple_eigen(sinp)
  USE input_simple_exc
  USE simple_objects
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE derived_objects

  IMPLICIT NONE 

  TYPE(input_options) :: sinp
  TYPE(exc), POINTER :: a(:) !excitons to be found
  TYPE(bands) :: bd
  TYPE(product) :: pd
  TYPE(potential) :: pt
  TYPE(prod_proj) :: pp
  TYPE(prod_mix) :: pm

  COMPLEX(kind=DP), ALLOCATABLE :: ene(:)!their energies
  INTEGER :: i

  ! read in excitonic Hamiltonian stuff   
  CALL read_bands(sinp,bd)
   
  !read in product stuff
  CALL initialize_product(pd)
  CALL read_product(sinp,pd)

  !read in potential stuff
  CALL initialize_potential(pt)
  CALL read_potential(sinp,pt)
   
  !set up product contractions
  CALL initialize_prod_proj(pp)
  CALL build_prod_proj(bd,pd,pp)
   
  !set up mixed contractions
  CALL initialize_prod_mix(pm)
  CALL build_prod_mix(sinp,bd,pd,pm,pt)

  ALLOCATE(a(sinp%nvec))
  DO i=1,sinp%nvec
    CALL setup_exc(bd,a(i))
  ENDDO 

  ALLOCATE(ene(sinp%nvec))
  !CALL solver
  SELECT CASE(sinp%diago)
  CASE(0) !steepest descent
    CALL diago_exc_sd(sinp,bd,pp,pt,pm,a)
  CASE(1) !CG
    CALL diago_exc_cg(sinp,bd,pp,pt,pm,a)
  END SELECT

  DO i=1,sinp%nvec
     CALL nice_write_exc(bd,sinp,a(i),i)
  ENDDO 

  DO i=1,sinp%nvec
     CALL deallocate_exc(a(i))
  ENDDO 
  DEALLOCATE(a)
  DEALLOCATE(ene)
  
  CALL deallocate_bands(bd)
  CALL deallocate_product(pd)
  CALL deallocate_potential(pt)
  CALL deallocate_prod_proj(pp)
  CALL deallocate_prod_mix(pm)

  RETURN 
END SUBROUTINE simple_eigen
