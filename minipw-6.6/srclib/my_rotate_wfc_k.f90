!----------------------------------------------------------------------------
SUBROUTINE my_rotate_wfc_k( h_psi, s_psi, overlap, &
                         npwx, npw, nstart, nbnd, npol, psi, evc, e )
!----------------------------------------------------------------------------
  !
  ! Serial version of rotate_wfc for colinear, k-point calculations
  !
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id, &
                            me_bgrp, root_bgrp
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INCLUDE 'laxlib.fh'
  !
  ! I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! number of spin polarizations
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP) :: psi(npwx*npol,nstart), evc(npwx*npol,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e(nbnd)
    ! eigenvalues
  !
  ! local variables
  !
  INTEGER :: kdim, kdmx
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
  REAL(DP),    ALLOCATABLE :: en(:)
  INTEGER :: n_start, n_end, my_n
  !
  EXTERNAL  h_psi,    s_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)

  IF( npol == 1 ) THEN
    !
    kdim = npw
    kdmx = npwx
    !
  ELSE
    !
    kdim = npwx*npol
    kdmx = npwx*npol
    !
  ENDIF
  !
  ALLOCATE( aux(kdmx, nstart ) )    
  ALLOCATE( hc( nstart, nstart) )    
  ALLOCATE( sc( nstart, nstart) )    
  ALLOCATE( vc( nstart, nstart) )    
  ALLOCATE( en( nstart ) )
  !
  ! Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  !      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  CALL h_psi( npwx, npw, nstart, psi, aux )
  !
  hc = (0.D0,0.D0)
  CALL divide(inter_bgrp_comm, nstart, n_start, n_end)
  my_n = n_end - n_start + 1;
  write(*,*) 'nstart, n_start, n_end =', nstart, n_start, n_end
  !
  if( n_start <= n_end ) then
    call ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi, kdmx, &
              & aux(1,n_start), kdmx, (0.D0,0.D0), hc(1,n_start), nstart )
  endif
  CALL mp_sum( hc, inter_bgrp_comm )
  CALL mp_sum( hc, intra_bgrp_comm )
  !
  sc = (0.D0,0.D0)
  IF( overlap ) THEN
    ! overlap is needed
    !
    CALL s_psi( npwx, npw, nstart, psi, aux )
    if(n_start <= n_end) then
      CALL ZGEMM( 'C','N', nstart, my_n, kdim, (1.D0,0.D0), psi, kdmx, &
                &  aux(1,n_start), kdmx, (0.D0,0.D0), sc(1,n_start), nstart )
    endif
    !
  ELSE
    ! overlap is not needed (it is identity?)
    if( n_start <= n_end ) then
      CALL ZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0,0.D0), psi, kdmx, &
                & psi(1,n_start), kdmx, (0.D0,0.D0), sc(1,n_start), nstart )
    endif
  ENDIF
  !
  CALL mp_sum( sc, inter_bgrp_comm )
  CALL mp_sum( sc, intra_bgrp_comm )
  !
  ! Diagonalize
  !
  CALL diaghg( nstart, nbnd, hc, sc, nstart, en, vc, me_bgrp, root_bgrp, intra_bgrp_comm )
  !
  e(:) = en(1:nbnd)
  !
  !  update the basis set
  !  
  aux = (0.D0,0.D0)
  if( n_start <= n_end ) then
    CALL ZGEMM( 'N','N', kdim, nbnd, my_n, (1.D0,0.D0), psi(1,n_start), &
              &  kdmx, vc(n_start,1), nstart, (0.D0,0.D0), aux, kdmx )
  endif
  CALL mp_sum( aux, inter_bgrp_comm )
  !     
  evc(:,:) = aux(:,1:nbnd)
  !
  DEALLOCATE( en )
  DEALLOCATE( vc )
  DEALLOCATE( sc )
  DEALLOCATE( hc )
  DEALLOCATE( aux )
  RETURN
  !
END SUBROUTINE my_rotate_wfc_k
