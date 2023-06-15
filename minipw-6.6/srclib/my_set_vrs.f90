!------------------------------------------------------------------------------
SUBROUTINE my_set_vrs( vrs, vltot, vr, kedtau, kedtaur, nrxx, nspin, doublegrid )
!-----------------------------------------------------------------------------
  !! Set the total local potential vrs on the smooth mesh to be used in 
  !! \(\texttt{h_psi}\), adding the (spin dependent) scf (H+xc) part and
  !! the sum of all the local pseudopotential contributions.
  !
  USE kinds
  USE funct,    ONLY : dft_is_meta
  USE fft_base, ONLY : dffts 
  !
  IMPLICIT NONE

  !! input: number of spin components: 1 if lda, 2 if lsd, 4 if noncolinear
  INTEGER :: nspin

  !! input: the fft grid dimension
  INTEGER :: nrxx

  !! output: total local potential on the smooth grid vrs=vltot+vr
  REAL(DP) :: vrs(nrxx,nspin)
  
  !! input: the total local pseudopotential
  REAL(DP) :: vltot(nrxx)
  
  !! input: the scf(H+xc) part of the local potential
  REAL(DP) :: vr(nrxx,nspin)

  !! position dependent kinetic energy enhancement factor
  REAL(DP) :: kedtau(dffts%nnr,nspin)

  !! the kinetic energy density in R-space
  REAL(DP) :: kedtaur(nrxx,nspin)
  
  ! input: true if a doublegrid is used
  LOGICAL :: doublegrid

  INTEGER :: i

  write(*,*)
  write(*,*) '------------------------------------------------------------'
  write(*,*) 'ENTER my_set_vrs'
  write(*,*) '------------------------------------------------------------'

  write(*,*)
  write(*,*) 'my_set_vrs: nrxx = ', nrxx
  write(*,*) 'my_set_vrs: shape vltot = ', shape(vltot)
  write(*,*) 'my_set_vrs: shape vrs = ', shape(vrs)
  write(*,*)
  write(*,*) 'my_set_vrs: values given in Ha'
  write(*,*)
  write(*,*) 'my_set_vrs: sum vltot before = ', sum(vltot)*0.5d0
  write(*,*) 'my_set_vrs: sum vr (dense) before = ', sum(vr)*0.5d0
  write(*,*) 'my_set_vrs: sum vrs (smooth) before = ', sum(vrs)*0.5d0

  CALL my_sum_vrs( nrxx, nspin, vltot, vr, vrs )
  CALL my_interpolate_vrs( nrxx, nspin, doublegrid, kedtau, kedtaur, vrs )

  write(*,*)
  write(*,*) 'my_set_vrs: sum vltot after = ', sum(vltot)*0.5d0
  write(*,*) 'my_set_vrs: sum vr (dense) after = ', sum(vr)*0.5d0
  write(*,*) 'my_set_vrs: sum vrs (smooth) after = ', sum(vrs)*0.5d0

  write(*,*)
  write(*,*) '------------------------------------------------------------'
  write(*,*) 'EXIT my_set_vrs'
  write(*,*) '------------------------------------------------------------'

  RETURN
  !
END SUBROUTINE my_set_vrs
!
!
!--------------------------------------------------------------------
SUBROUTINE my_sum_vrs( nrxx, nspin, vltot, vr, vrs )
!--------------------------------------------------------------------
  !! Accumulates local potential contributions into vrs (the total local potential).
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER :: nspin
  !! input: number of spin components: 1 if lda, 2 if lsd, 4 if noncolinear
  
  !! input: the fft grid dimension
  INTEGER :: nrxx

  !! output: total local potential on the smooth grid:  
  !! \(\text{vrs\}=\text{vltot}+\text{vr}\)
  REAL(DP) :: vrs(nrxx,nspin)

  !! input: the total local pseudopotential  
  REAL(DP) :: vltot(nrxx)

  !! input: the scf(H+xc) part of the local potential  
  REAL(DP) :: vr(nrxx,nspin)

  ! local variable
  INTEGER :: is
  !
  !
  DO is = 1, nspin
    !
    ! define the total local potential (external + scf) for each spin ...
    !
    IF (is > 1 .AND. nspin == 4) THEN
      !
      ! noncolinear case: only the first component contains vltot
      ! ffr: why?
      !
      vrs(:,is) = vr(:,is)
    ELSE
      vrs(:,is) = vltot(:) + vr(:,is)
    ENDIF
    !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE my_sum_vrs


!--------------------------------------------------------------------------
SUBROUTINE my_interpolate_vrs( nrxx, nspin, doublegrid, kedtau, kedtaur, vrs )
!--------------------------------------------------------------------------
  !! Interpolates local potential on the smooth mesh if necessary.
  !
  USE kinds
  USE funct,           ONLY : dft_is_meta
  USE fft_base,        ONLY : dffts, dfftp
  USE fft_interfaces,  ONLY : fft_interpolate
  !
  IMPLICIT NONE

  !! input: number of spin components: 1 if lda, 2 if lsd, 4 if noncolinear
  INTEGER :: nspin
  
  !! input: the fft grid dimension
  INTEGER :: nrxx
  
  !! output: total local potential interpolated on the smooth grid
  REAL(DP) :: vrs(nrxx,nspin)
  
  !! position dependent kinetic energy enhancement factor
  REAL(DP) :: kedtau(dffts%nnr,nspin)

  !! the kinetic energy density in R-space
  REAL(DP) :: kedtaur(nrxx,nspin)

  !! input: true if a doublegrid is used
  LOGICAL :: doublegrid

  ! local variable
  INTEGER :: is

  ! interpolate it on the smooth mesh if necessary
  write(*,*) 'my_interpolate_vrs: sum vrs before interpolate: ', sum(vrs)*0.5d0

  if( dft_is_meta() ) then
    write(*,*)
    write(*,*) 'Before fft_interpolate'
    write(*,*) 'sum(kedtau) = ', sum(kedtau)
    write(*,*) 'sum(kedtaur) = ', sum(kedtaur)
  endif

  DO is = 1, nspin
    IF( doublegrid ) CALL fft_interpolate( dfftp, vrs(:, is), dffts, vrs(:,is) )
    IF( dft_is_meta() ) CALL fft_interpolate( dfftp, kedtaur(:,is), dffts, kedtau(:,is) )
  ENDDO

  if( dft_is_meta() ) then
    write(*,*)
    write(*,*) 'After fft_interpolate'
    write(*,*) 'sum(kedtau) = ', sum(kedtau)
    write(*,*) 'sum(kedtaur) = ', sum(kedtaur)
  endif

  write(*,*) 'my_interpolate_vrs: sum vrs after interpolate: ', sum(vrs)*0.5d0

  RETURN

END SUBROUTINE my_interpolate_vrs
