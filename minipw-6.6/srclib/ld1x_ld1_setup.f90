!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
SUBROUTINE ld1_setup
!---------------------------------------------------------------
  !
  !     this routine computes the general variables needed in
  !     the atomic calculation
  !
  !
  USE kinds, ONLY : dp
  USE ld1inc, ONLY : iswitch, lsd, rel,sl3, lmx, nconf, nn, ll, &
                nbeta, lloc, nsloc, etot, etot0, etots0, vxt, tau, vtau, &
                file_wavefunctions, file_logder,  &
                file_wavefunctionsps, file_logderps, file_pawexp, &
                prefix, new, core_state, enls, enl,  &
                nwf,                  jj,    el,    isw,    oc, nstoae, &
                nwfs,          lls,   jjs,   els,   isws,   ocs, &
                nwfts,  nnts,  llts,  jjts,  elts,  iswts,  octs, nstoaets, &
                nwftsc, nntsc, lltsc, jjtsc, eltsc, iswtsc, octsc, nstoaec, lpaw
  USE funct, ONLY : get_iexch, dft_is_meta, start_exx, dft_is_nonlocc 
  IMPLICIT NONE

  INTEGER :: n, n1, nc
  LOGICAL :: hf, oep, meta, non_locc
  real(DP) :: ocs_loc
  
  write(*,*)
  write(*,*) '<div> ENTER ld1_setup'
  write(*,*)
  
  write(*,*) 'nwfs = ', nwfs
  write(*,*) 'nwf = ', nwf

  write(*,*) 'Label for all states'
  do n = 1,nwf
    write(*,'(1x,A2,A)',advance="no") el(n), ' '
  enddo
  !
  write(*,*)
  write(*,*) 'Label for valence states'
  do n = 1,nwfs
    write(*,'(1x,A2,A)',advance="no") els(n), ' '
  enddo
  write(*,*)

  !
  !
  meta = dft_is_meta()
  IF ( meta .and. rel > 1 ) THEN
    CALL errore('ld1_setup','meta-GGA not implemented for fully-relativistic case', 1)
  ENDIF
  !
  IF ( meta .and. lsd == 1 ) THEN
    CALL errore('ld1_setup','meta-GGA not implemented for LSDA', 2)
  ENDIF
  !
  IF( meta .and. iswitch > 1 ) THEN
    CALL errore('ld1_setup','meta-GGA implemented only for all-electron case', 3)
  ENDIF
  !
  hf = ( get_iexch() == 5 )
  IF( hf ) THEN 
    CALL errore('ld1_setup','HF not implemented yet',1)
  ENDIF
  !
  oep = ( get_iexch() == 4 )
  IF( oep .and. iswitch > 1 ) THEN
    CALL errore('ld1_setup', 'OEP is implemented only for all-electron calc.', 1)
  ENDIF
  !
  IF( oep .and. rel > 0) THEN
    CALL errore('ld1_setup', 'OEP is implemented only for non-relativistic calc.', 1)
  ENDIF
  IF( oep ) THEN 
    CALL start_exx()
  ENDIF
  ! 
  non_locc = dft_is_nonlocc()
  IF ( non_locc ) THEN
    CALL errore('ld1_setup','non-local functionals not implemented yet', 1)
  ENDIF
  !
  CALL set_sl3(sl3, lmx)
  !
  ! make the correspondence all-electron pseudopotential
  !
  WRITE(*,*) 'lsd = ', lsd
  IF (iswitch >= 3) THEN
    DO nc = 1, nconf
      DO n = 1, nwftsc(nc)
        nstoaec(n,nc) = 0
        DO n1 = 1,nwf
          IF( lsd == 1 ) THEN
            IF( eltsc(n,nc)==el(n1) .and. iswtsc(n,nc)==isw(n1) ) THEN
              nstoaec(n,nc) = n1
            ENDIF
          ELSE
            ! This is the usual case
            IF( eltsc(n,nc)==el(n1) .and. abs(jjtsc(n,nc)-jj(n1)) < 1.d-7 ) THEN
              nstoaec(n,nc) = n1
            ENDIF
          ENDIF
        ENDDO
        !
        IF( nstoaec(n,nc)==0 ) THEN
          CALL errore('ld1_setup', 'all electron wfc corresponding to pseudo-state '//eltsc(n,nc)//' not found', nc)
        ENDIF
      ENDDO
    ENDDO
    !
    ! set the test configuration for descreening
    !
    IF( iswitch == 3 ) THEN
      nwfts = nwftsc(1)
      DO n = 1,nwfts
        nnts(n) = nntsc(n,1)
        llts(n) = lltsc(n,1)
        elts(n) = eltsc(n,1)
        jjts(n) = jjtsc(n,1)
        iswts(n) = iswtsc(n,1)
        octs(n) = octsc(n,1)
        nstoaets(n) = nstoaec(n,1)
      ENDDO
    ENDIF
  ENDIF
  !
  new(:) = .false.
  !
  ! divide the core and valence states
  !
  gen_and_test:&
  IF( iswitch == 3 ) THEN
    isws = 1
    !
    DO n = 1,nwf
      core_state(n) = .true.
    ENDDO
    !
    DO n = 1,nwfs
      nstoae(n) = 0
      DO n1 = 1,nwf
        IF( rel==2 ) THEN
          IF( els(n)==el(n1) .and. jjs(n)==jj(n1) ) THEN
            nstoae(n) = n1
            core_state(n1) = .false.
          ENDIF
        ELSE
          IF(els(n)==el(n1) ) THEN
            nstoae(n) = n1
            core_state(n1) = .false.
          ENDIF
        ENDIF
      ENDDO
      !
      IF( nstoae(n)==0 ) THEN
        CALL errore('ld1_setup', 'no all electron for this ps',n)
      ENDIF
      !
      IF( enls(n) /= 0.0_dp) THEN
        new(n) = .true.
      ENDIF
    ENDDO
    !
    IF( lloc > -1 ) THEN
      nsloc = nwfs
      nbeta = nwfs - 1
      IF( rel==2 .and. lloc/=0 ) THEN
        nsloc = nwfs - 1
        nbeta = nwfs - 2
        IF( lls(nsloc+1) /= lloc ) THEN
          CALL errore('ld1_setup', 'mismatch between lloc and l of spin-orbit split wfc chosen for local potential', nsloc)
        ENDIF
      ENDIF
      !
      IF( lls(nsloc) /= lloc) THEN
        IF (rel==2) THEN
          CALL errore('ld1_setup', 'mismatch between lloc and l of spin-orbit split wfc chosen for local potential', nsloc)
        ELSE
          CALL errore('ld1_setup', 'mismatch between lloc and l of the wavefunction chosen for local potential', nsloc)
        ENDIF
      ENDIF
      !
      ocs_loc = ocs(nsloc)
      IF( rel==2 .and. lloc > 0 ) THEN 
        ocs_loc = ocs_loc + ocs(nsloc+1)
      ENDIF
      IF( lpaw .and. ocs_loc > 0.0_DP) THEN 
        CALL errore('ld1_setup','Paw generation with electrons in the local channel is not available',1)
      ENDIF
    ELSE
      nsloc = -1
      nbeta = nwfs
    ENDIF
    !
    ! test the occupations: for pseudopotential generation
    ! all-electron and pseudopotential occupations must match
    DO n = 1,nwfs
      IF( .not. new(n) ) THEN
        IF( abs(oc(nstoae(n)) - ocs(n)) > 1.0d-8 ) THEN 
          CALL errore('ld1_setup','mismatched all-electron/pseudo occupations',n)
        ENDIF
      ENDIF
    ENDDO
  ENDIF gen_and_test
  !
  ! zero the external potential, metaGGA potential, total energies
  !
  vxt = 0.0_dp
  vtau = 0.0_dp
  tau = 0.0_dp
  etot0 = 0.0_dp
  etots0 = 0.0_dp
  enl = 0.0_dp
  !
  ! initialize file names (used in all_electron)
  !
  file_wavefunctions = trim(prefix)//'.wfc'
  file_wavefunctionsps = trim(prefix)//'ps.wfc'
  file_logder = trim(prefix)//'.dlog'
  file_logderps = trim(prefix)//'ps.dlog'
  file_pawexp = trim(prefix)//'.pwe'

  write(*,*) 'nwf = ', nwf
  write(*,*) 'nn = ', nn(1:nwf)
  write(*,*) 'll = ', ll(1:nwf)
  write(*,*) 'oc = ', oc(1:nwf)
  write(*,*) 'core_state = ', core_state(1:nwf)

  write(*,*)
  write(*,*) '</div> EXIT ld1_setup'
  write(*,*)

  RETURN
END SUBROUTINE ld1_setup

