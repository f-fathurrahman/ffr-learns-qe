! This will be called in setup.f90

!------------------------------------------------------------------------
SUBROUTINE my_exx_base_exx_grid_init( reinit )
!------------------------------------------------------------------------
  use kinds, only : dp
  !
  USE symm_base,         ONLY : nsym, s
  USE cell_base,         ONLY : bg, at, tpiba
  USE spin_orb,          ONLY : domag
  USE noncollin_module,  ONLY : nspin_lsda, noncolin
  USE klist,             ONLY : xk, nkstot, nks, qnorm
  USE start_k,           ONLY : nk1,nk2,nk3
  USE control_flags,     ONLY : iverbosity
  !
  use exx_base, only : index_sym, index_xk, index_xkq, xkq_collect, &
                     & x_gamma_extrapolation, nqs, nq1, nq2, nq3, nkqs, &
                     & grid_factor, eps, exx_grid_initialized
  !
  IMPLICIT NONE
  !
  LOGICAL, OPTIONAL :: reinit
  !! reinitialize exx if .TRUE.
  !
  ! local variables
  !
  INTEGER :: isym, ik, ikq, iq, max_nk, temp_nkqs, idx, sign_
  INTEGER :: nqx(3)
  INTEGER, ALLOCATABLE :: temp_index_xk(:), temp_index_sym(:)
  INTEGER, ALLOCATABLE :: temp_index_ikq(:)
  REAL(DP), ALLOCATABLE :: temp_xkq(:,:), xk_collect(:,:)
  LOGICAL :: xk_not_found
  REAL(DP) :: sxk(3), dxk(3), xk_cryst(3)
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  CALL start_clock ( 'exx_grid' )
  !
  IF ( PRESENT (reinit) ) THEN
    IF ( reinit ) THEN
      IF (ALLOCATED(xkq_collect))  DEALLOCATE( xkq_collect )
      IF (ALLOCATED(index_xk)   )  DEALLOCATE( index_xk    )
      IF (ALLOCATED(index_sym)  )  DEALLOCATE( index_sym   )
      exx_grid_initialized = .FALSE.
      nkqs = 0
    ENDIF
  ENDIF
  !
  !ffr: this is the default value
  IF (nq1 <= 0) nq1 = nk1
  IF (nq2 <= 0) nq2 = nk2
  IF (nq3 <= 0) nq3 = nk3
  !
  !ffr: this is for molecule
  IF (nkstot == nspin_lsda) THEN
    nq1=1; nq2=1; nq3=1
  ENDIF
  !
  IF (ANY( (/nq1,nq2,nq3/) <=0 )) CALL errore( 'exx_grid_init', "wrong EXX q grid", 1 )
  !
  IF (exx_grid_initialized) CALL errore( 'exx_grid_init', "grid already initialized", 1 )
  exx_grid_initialized = .TRUE.
  !
  ! definitions and checks
  !
  grid_factor = 1._dp
  IF (x_gamma_extrapolation) grid_factor = 8.d0/7.d0
  !
  nqs = nq1 * nq2 * nq3
  !
  ! all processors on all pools need to have access to all k+q points
  !
  ALLOCATE( xk_collect(3,nkstot) )
  !
  CALL poolcollect( 3, nks, xk, nkstot, xk_collect )
  !
  ! set a safe limit as the maximum number of auxiliary points we may need
  ! and allocate auxiliary arrays
  max_nk = nkstot * MIN(48, 2 * nsym)
  ALLOCATE( temp_index_xk(max_nk), temp_index_sym(max_nk) )
  ALLOCATE( temp_index_ikq(max_nk) )
  ALLOCATE( temp_xkq(3,max_nk) )
  !
  ! find all k-points equivalent by symmetry to the points in the k-list
  !
  temp_nkqs = 0
  !
  !ffr: loop over symmetry elements
  DO isym = 1, nsym
    DO ik = 1, nkstot
      !
      ! go to crystalline coordinates
      xk_cryst(:) = xk_collect(:,ik)
      CALL cryst_to_cart( 1, xk_cryst, at, -1 )
      ! rotate with this sym.op.
      sxk(:) = s(:,1,isym)*xk_cryst(1) + s(:,2,isym)*xk_cryst(2) + s(:,3,isym)*xk_cryst(3)
      ! add sxk to the auxiliary list IF it is not already present
      xk_not_found = .TRUE.
      ! *** do-loop skipped the first time because temp_nkqs == 0
      DO ikq = 1, temp_nkqs
        IF (xk_not_found ) THEN
          dxk(:) = sxk(:)-temp_xkq(:,ikq) - NINT(sxk(:)-temp_xkq(:,ikq))
          IF ( ABS(dxk(1))<=eps .AND. ABS(dxk(2))<=eps .AND. ABS(dxk(3))<=eps ) THEN    
            xk_not_found = .FALSE.
          ENDIF
        ENDIF
      ENDDO
      !
      IF (xk_not_found) THEN
        temp_nkqs = temp_nkqs + 1
        temp_xkq(:,temp_nkqs) = sxk(:)
        temp_index_xk(temp_nkqs) = ik
        temp_index_sym(temp_nkqs) = isym
      ENDIF
      !
      sxk(:) = -sxk(:)
      xk_not_found = .TRUE.
      DO ikq = 1, temp_nkqs
        IF (xk_not_found ) THEN
          dxk(:) = sxk(:) - temp_xkq(:,ikq) - NINT(sxk(:) - temp_xkq(:,ikq))
          IF ( ABS(dxk(1))<=eps .AND. ABS(dxk(2))<=eps .AND. ABS(dxk(3))<=eps ) THEN
            xk_not_found = .FALSE.
          ENDIF
        ENDIF
      ENDDO
      !
      IF (xk_not_found .AND. .NOT. (noncolin.AND.domag) ) THEN
        temp_nkqs = temp_nkqs + 1
        temp_xkq(:,temp_nkqs) = sxk(:)
        temp_index_xk(temp_nkqs)  = ik
        temp_index_sym(temp_nkqs) =-isym
      ENDIF
      !
    ENDDO
  ENDDO
  !
  ! Find good q-point grid. Decrease the nqX until a good grid is found or
  ! until it is 1 x 1 x 1 (always good)
  idx = 1
  sign_ = -1
  nqx = (/nq1, nq2, nq3/)
  DO WHILE(.TRUE.)
    CALL my_exx_qgrid_init(temp_nkqs, xk_collect, temp_xkq, &
                        nkqs, temp_index_ikq, dxk)

    ! Good q-point mesh
    IF (ALL(ABS(dxk) < eps ) ) THEN
      !
      IF (idx > 1) &
        WRITE(*, '(5x,a)') "EXX: WARNING: q-point mesh has been updated!"
      !
      WRITE(*, '(5x,a,3i5)') "EXX: q-point mesh: ", nq1, nq2, nq3
      EXIT ! DO WHILE
    ENDIF
    !
    ! Try q-points around the input mesh, prioritizing smaller mesh
    !
    nq1 = nqx(1) + idx * sign_
    nq2 = nqx(2) + idx * sign_
    nq3 = nqx(3) + idx * sign_
    !
    ! Ensure no values smaller than 1
    IF (nq1 < 1) nq1 = 1
    IF (nq2 < 1) nq2 = 1
    IF (nq3 < 1) nq3 = 1
    !
    ! Enforce nqX <= nkX. This is important for surfaces to keep the
    ! Z q-point 1.
    !
    IF (nq1 > nk1) nq1 = nk1
    IF (nq2 > nk2) nq2 = nk2
    IF (nq3 > nk3) nq3 = nk3
    !
    nqs = nq1 * nq2 * nq3
    !
    sign_ = -1 * sign_
    !
    ! Increase idx every other time sign is changed
    IF (sign_ < 0) idx = idx + 1
    !
  ENDDO
  !
  ! allocate and fill the arrays xkq(3,nkqs), index_xk(nkqs) and index_sym(nkqs)
  ! NOTE: nkqs will be redefined as nspin_lsda*nkqs later 
  !
  ALLOCATE( xkq_collect(3,nspin_lsda*nkqs), index_xk(nspin_lsda*nkqs),  &
            index_sym(nspin_lsda*nkqs) )
  !
  DO ik = 1, nkqs
    ikq               = temp_index_ikq(ik)
    xkq_collect(:,ik) = temp_xkq(:,ikq)
    index_xk(ik)      = temp_index_xk(ikq)
    index_sym(ik)     = temp_index_sym(ikq)
  ENDDO
  !
  CALL cryst_to_cart( nkqs, xkq_collect, bg, +1 )
  !
  IF (nkqs > 1) THEN
    WRITE(*, '(5x,3a)') "EXX: setup a grid of "//TRIM(int_to_char(nkqs))&
                          //" q-points centered on each k-point"
    IF ( nkqs < 100 .OR. iverbosity > 0 ) THEN
        WRITE( *, '(5x,a)' ) '(k+q)-points:'
        DO ik = 1, nkqs
          WRITE( *, '(3f12.7,5x,2i5)') (xkq_collect(ikq,ik), ikq=1,3), &
                index_xk(ik), index_sym(ik)
        ENDDO
    ELSE
        WRITE( *, '(5x,a)' ) "(set verbosity='high' to see the list)"
    END IF
  ELSE
    WRITE(*, '(5X,"EXX: grid of k+q points same as grid of k-points")')
  ENDIF
  !
  ! if nspin == 2, the kpoints are repeated in couples (spin up, spin down)
  IF (nspin_lsda == 2) THEN
    DO ik = 1, nkstot/2
        DO iq = 1, nqs
          index_xkq(nkstot/2+ik,iq) = index_xkq(ik,iq) + nkqs
        ENDDO
    ENDDO
    !
    DO ikq = 1, nkqs
      xkq_collect(:,ikq+nkqs) = xkq_collect(:,ikq)
      index_xk(ikq + nkqs)  = index_xk(ikq) + nkstot/2
      index_sym(ikq + nkqs) = index_sym(ikq)
    ENDDO
    nkqs = 2 * nkqs
  ENDIF
  !
  ! clean up
  DEALLOCATE( temp_index_xk, temp_index_sym, temp_index_ikq, temp_xkq )
  !
  ! check that everything is what it should be
  CALL exx_grid_check( xk_collect(:,:) )
  DEALLOCATE( xk_collect )
  !
  ! qnorm = max |q|, used in allocate_nlpot to compute the maximum size
  !         of some arrays (e.g. qrad) - beware: needed for US/PAW+EXX
  !
  qnorm = 0.0_dp
  DO iq = 1, nkqs
    DO ik = 1, nks
      qnorm = MAX(qnorm, SQRT( SUM((xk(:,ik)-xkq_collect(:,iq))**2) ))
    ENDDO
  ENDDO
  qnorm = qnorm * tpiba
  !
  CALL stop_clock( 'exx_grid' )
  !
  RETURN
  !
END SUBROUTINE my_exx_base_exx_grid_init