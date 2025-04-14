subroutine ld1x_debug_setup_gen_and_test()
  USE kinds, ONLY : dp
  USE ld1inc
  !
  implicit none
  !
  INTEGER :: n, n1
  real(DP) :: ocs_loc

  if( iswitch /= 3 ) then
    stop 'Need iswitch==3'
  endif

  isws = 1
  !
  DO n = 1,nwf
    core_state(n)=.true.
  ENDDO
  !
  DO n = 1,nwfs
    nstoae(n) = 0
    DO n1 = 1,nwf
      IF( rel==2 ) THEN
        IF( els(n) == el(n1) .and. jjs(n)==jj(n1) ) THEN
          nstoae(n)=n1
          core_state(n1)=.false.
        ENDIF
      ELSE
        IF( els(n) == el(n1) ) THEN
          nstoae(n) = n1
          core_state(n1) = .false.
        ENDIF
      ENDIF
    ENDDO
    IF( nstoae(n) == 0 ) then 
      CALL errore('ld1_setup', 'no all electron for this ps',n)
    endif
    IF( enls(n) /= 0.0_dp ) then
      new(n) = .true.
    endif
  ENDDO
  IF( lloc > -1 ) THEN
    nsloc = nwfs
    nbeta = nwfs-1
    IF( rel == 2 .and. lloc /=0 ) THEN
      nsloc = nwfs - 1
      nbeta = nwfs - 2
      IF( lls(nsloc+1) /= lloc ) &
         CALL errore('ld1_setup','mismatch between lloc and l of ' // &
       &           'spin-orbit split wfc chosen for local potential',nsloc)
    ENDIF
    IF( lls(nsloc) /= lloc) THEN
       IF(rel==2) THEN
          CALL errore('ld1_setup','mismatch between lloc and l of ' // &
       &            'spin-orbit split wfc chosen for local potential',nsloc)
       ELSE
          CALL errore('ld1_setup','mismatch between lloc and l of ' // &
       &            'the wavefunction chosen for local potential',nsloc)
       ENDIF
    ENDIF
    !
    ocs_loc = ocs(nsloc)
    IF(rel==2.and.lloc>0) ocs_loc=ocs_loc+ocs(nsloc+1)
    IF (lpaw .and. ocs_loc>0.0_DP) &
        CALL errore('ld1_setup','Paw generation with electrons' // &
        &    'in the local channel is not available',1)
  ELSE
    nsloc=-1
    nbeta=nwfs
  ENDIF
  !
  ! test the occupations: for pseudopotential generation
  ! all-electron and pseudopotential occupations must match
  !
  DO n=1,nwfs
    IF(.not. new(n)) THEN
      IF(abs(oc(nstoae(n))-ocs(n)) > 1.0d-8 ) CALL errore &
         ('ld1_setup','mismatched all-electron/pseudo occupations',n)
    ENDIF
  ENDDO

  call ld1x_print_variables()

end subroutine