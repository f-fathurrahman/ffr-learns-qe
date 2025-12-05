!----------------------------------------------------------------------------
SUBROUTINE my_wfcinit()
!----------------------------------------------------------------------------
  ! This routine computes an estimate of the starting wavefunctions
  ! from superposition of atomic wavefunctions and/or random wavefunctions.
  ! It also open needed files or memory buffers
  USE io_global,            ONLY : stdout
  USE basis,                ONLY : natomwfc, starting_wfc
  USE bp,                   ONLY : lelfield
  USE klist,                ONLY : xk, nks, ngk, igk_k
  USE control_flags,        ONLY : io_level, lscf
  USE fixed_occ,            ONLY : one_atom_occupations
  USE ldaU,                 ONLY : lda_plus_u, U_projection, wfcU, lda_plus_u_kind
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE io_files,             ONLY : nwordwfc, nwordwfcU, iunhub, iunwfc,&
                                   diropn, xmlfile, restart_dir
  USE buffers,              ONLY : open_buffer, close_buffer, get_buffer, save_buffer
  USE uspp,                 ONLY : nkb, vkb
  USE wavefunctions,        ONLY : evc
  USE wvfct,                ONLY : nbnd, current_k
  USE mp,                   ONLY : mp_bcast, mp_sum
  IMPLICIT NONE
  !
  INTEGER :: ik
  LOGICAL :: exst_mem, exst_file

  write(*,*)
  write(*,*) '<div> ENTER my_wfcinit'
  write(*,*)

  !
  ! Orthogonalized atomic functions needed for DFT+U and other cases
  !
  !IF ( use_wannier .OR. one_atom_occupations ) CALL orthoatwfc ( use_wannier )
  IF ( one_atom_occupations ) CALL orthoatwfc ( .false. ) ! ffr
  IF ( lda_plus_u ) CALL orthoUwfc()
  !
  ! open files/buffer for wavefunctions (nwordwfc set in openfil)
  ! io_level > 1 : open file, otherwise: open buffer
  !

  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst_file )
  !
  IF( TRIM(starting_wfc) == 'file') THEN
    stop 'disabled in my_wfcinit'
  ENDIF

  !
  ! state what will happen
  !
  IF ( TRIM(starting_wfc) == 'file' ) THEN
     !
     WRITE( stdout, '(5X,"Starting wfcs from file")' )
     !
  ELSE IF ( starting_wfc == 'atomic' ) THEN
     !
     IF ( natomwfc >= nbnd ) THEN
        WRITE( stdout, '(5X,"Starting wfcs are ",I4," atomic wfcs")' ) natomwfc
     ELSE
        WRITE( stdout, '(5X,"Starting wfcs are ",I4," atomic + ", &
             &           I4," random wfcs")' ) natomwfc, nbnd-natomwfc
     END IF
     !
  ELSE IF ( TRIM(starting_wfc) == 'atomic+random' .AND. natomwfc > 0) THEN
     !
     IF ( natomwfc >= nbnd ) THEN
        WRITE( stdout, '(5X,"Starting wfcs are ",I4," randomized atomic wfcs")')&
             natomwfc
     ELSE
        WRITE( stdout, '(5X,"Starting wfcs are ",I4," randomized atomic wfcs + "&
             &          ,I4," random wfcs")' ) natomwfc, nbnd-natomwfc
     END IF
     !
  ELSE
     !
     WRITE( stdout, '(5X,"Starting wfcs are random")' )
     !
  END IF
  !
  ! exit here if starting from file or for non-scf calculations.
  ! In the latter case the starting wavefunctions are not 
  ! calculated here but just before diagonalization (to reduce I/O)
  !
  IF (  ( .NOT. lscf .AND. .NOT. lelfield ) .OR. TRIM(starting_wfc) == 'file' ) THEN
    RETURN
  ENDIF

  ! calculate and write all starting wavefunctions to buffer
  DO ik = 1, nks
    !
    ! Hpsi initialization: k-point index, spin, kinetic energy
    !
    current_k = ik
    IF( lsda ) then
      current_spin = isk(ik)
      write(*,*) 'current_spin index = ', current_spin
    endif
    call g2_kin(ik)  ! calculate (G+k)^2
    !
    ! More Hpsi initialization: nonlocal pseudopotential projectors |beta>
    !
    IF ( nkb > 0 ) then
      ! XXX why need this?
      ! XXX This will compute projectors
      CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
    endif
    !
    ! Needed for DFT+U
    !
    IF( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) then
      CALL get_buffer( wfcU, nwordwfcU, iunhub, ik )
    endif
    !
    ! DFT+U+V: calculate the phase factor at a given k point
    !
    IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) then
      CALL phase_factor(ik)
    endif
    !
    ! calculate starting wavefunctions (calls Hpsi)
    !
    CALL my_init_wfc( ik )
    !
    ! write  starting wavefunctions to file
    !
    IF( nks > 1 .OR. (io_level > 1) .OR. lelfield ) then
      CALL save_buffer( evc, nwordwfc, iunwfc, ik )
    endif
    !
  ENDDO

  write(*,*)
  write(*,*) '</div> EXIT my_wfcinit'
  write(*,*)

  RETURN
  !
END SUBROUTINE my_wfcinit


