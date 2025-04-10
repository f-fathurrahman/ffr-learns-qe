!--------------------------------------------------------------------------
SUBROUTINE my_pseudovloc()
!--------------------------------------------------------------------------
  !
  ! This routine generate a local pseudopotential 
  ! The output of the routine are:
  ! vpsloc: the local pseudopotential
  !      
  USE kinds, ONLY : DP
  USE radial_grids, ONLY : ndmx
  USE io_global, ONLY : stdout
  USE ld1inc, ONLY : lloc, rcloc, grid, vpot, vpsloc, rel, nsloc, &
                     phis, els, chis, psipsus, &
                     jjs, nstoae, enls, new, psi, enl, rcut, psipaw, &
                     psipaw_rel
  IMPLICIT NONE

  INTEGER :: &
       nwf0, &  ! used to specify the all electron function
       nst,  &  ! auxiliary
       ik       ! the point corresponding to rc

  REAL(DP) ::             &
       xc(8),              &  ! the coefficients of the fit
       vaux(ndmx,2),        &  ! keeps the potential
       psi_in(ndmx)            ! auxiliary

  INTEGER ::         &
       n,        &  ! counter on mesh points
       ns,       &  ! auxiliary
       indi,rep, &  ! auxiliary
       indns(0:1)    ! auxiliary

  write(*,*)
  write(*,*) '**** ENTER my_pseudovloc'
  write(*,*)


  IF( lloc < 0 ) THEN
    !
    ! Compute the potential by smoothing the AE potential
    !
    ! Compute the ik which correspond to this cutoff radius
    !
    WRITE(stdout, &
         "(/,5x,' Generating local potential from pseudized AE potential:',&
           &  /,5x,' Matching radius rcloc = ',f8.4)") rcloc
    ik = 0
    ! find largest index n for which r(n) < rcloc
    DO n = 1,grid%mesh
      IF( grid%r(n) < rcloc ) ik = n
    ENDDO
    WRITE(*,*) 'ik = ', ik
    !
    IF( mod(ik,2) == 0 ) THEN
      WRITE(*,*) 'Found ik is even number, making it odd'
      ik = ik + 1
    ENDIF
    !
    IF( ik <= 1 .or. ik > grid%mesh ) THEN
      CALL errore('pseudovloc', 'wrong matching point', 1)
    ENDIF
    WRITE(*,*) 'grid%mesh = ', grid%mesh
    WRITE(*,*) 'my_pseudovloc: Found ik = ', ik
    !
    ! smooth the potential before ik.
    !
    ! with the original recipe
    IF( lloc==-1 ) THEN
      ! used for example Si_paw
      WRITE(*,*) 'Calling my_compute_potps'
      CALL my_compute_potps(ik, vpot, vpsloc, xc)
    ENDIF
    !
    ! or with a modified recipe that enforce V''(0)=0 as suggested by TM
    IF( lloc==-2 ) THEN
      WRITE(stdout,"(5x,' Enforcing V''''(0)=0 (lloc=-2)')")
    ENDIF
    ! XXX: merge this if statement?
    IF( lloc==-2 ) THEN
      CALL compute_potps_new(ik, vpot, vpsloc, xc)
    ENDIF
    
    WRITE(stdout, 110) grid%r(ik), xc(5)**2 
110  FORMAT(/5x, ' Local pseudo, rcloc=',f6.3,' Estimated cut-off energy= ', f8.2,' Ry')
  
    !
    !
  ELSE
    !
    ! if a given angular momentum gives the local component this is done here
    !
    WRITE(*,*) 'lloc is given: lloc = ', lloc
    nst = (lloc + 1)*2
    WRITE(*,*) 'nst = ', nst
    !
    IF( rel==2 .and. lloc > 0 ) THEN
      rep = 1
      indns(0) = nsloc
      indns(1) = nsloc + 1
      IF( jjs(nsloc) > jjs(nsloc+1) ) then
        indns(0) = nsloc + 1
        indns(1) = nsloc
      ENDIF
    ELSE
      rep = 0
      indns(0) = nsloc
    ENDIF
    vpsloc = 0.0_dp
    vaux = 0.0_dp
    !
    WRITE(*,*) 'my_pseudovloc: rep = ', rep
    !
    DO indi = 0,rep
      nwf0 = nstoae(nsloc+indi)
      IF( enls(nsloc+indi) == 0.0_dp ) then
        enls(nsloc+indi) = enl(nwf0)
      ENDIF
      !
      ! compute the ik closer to r_cut
      !
      ik = 0
      DO n = 1,grid%mesh
        IF( grid%r(n) < rcut(nsloc+indi) ) THEN
          ik = n
        ENDIF
      ENDDO
      IF( mod(ik,2) == 0) THEN
        ik = ik + 1
      ENDIF
      IF( ik <= 1 .or. ik > grid%mesh ) THEN
        CALL errore('pseudovloc','wrong matching point',1)
      ENDIF
      !
      rcloc = rcut(nsloc+indi)
      !
      IF( rep == 0 ) THEN
        WRITE(stdout,"(/,5x,' Generating local pot.: lloc=',i1, &
                & ', matching radius rcloc = ',f8.4)") lloc, rcloc
      ELSE
        IF( rel==2 ) THEN
          WRITE(stdout,"(/,5x,' Generating local pot.: lloc=',i1, &
               &', j=',f5.2,', matching radius rcloc = ',f8.4)") &
               lloc, lloc-0.5d0+indi, rcloc
        ELSE
          WRITE(stdout,"(/,5x,' Generating local pot.: lloc=',i1, &
               &', spin=',i1,', matching radius rcloc = ',f8.4)") &
               lloc, indi+1, rcloc
        ENDIF
      ENDIF
      !
      ! compute the phi functions
      !
      ns = indns(indi)
      WRITE(*,*) 'ns = ', ns
      IF( new(ns) ) then
        WRITE(*,*) 'Calling set_psi_in'
        CALL my_set_psi_in(ik,lloc,jjs(ns),enls(ns),psi_in,psipaw_rel)
      ELSE
        psi_in(:) = psi(:,1,nwf0)
      ENDIF
      psipaw(:,ns) = psi_in(:)
      !
      ! compute the phi and chi functions
      !
      WRITE(*,*) 'Calling compute_phi_tm and compute_chi_tm'
      CALL compute_phi_tm( lloc, ik, psi_in, phis(1,ns), 0, xc,enls(ns), els(ns) )
      CALL compute_chi_tm( lloc, ik, ik+10, phis(1,ns), chis(1,ns), xc, enls(ns) )
      !
      ! set the local potential equal to the all-electron one at large r
      !
      DO n=1,grid%mesh
        IF (grid%r(n) > rcloc) then
          vaux(n,indi+1)=vpot(n,1)
        ELSE
          vaux(n,indi+1)=chis(n,ns)/phis(n,ns)
        ENDIF
      ENDDO
      !
      psipsus(:,ns)=phis(:,ns)
      !
    ENDDO  ! indi=0,rep
    !
    IF( rep==0 ) THEN
      DO n=1,grid%mesh
        vpsloc(n) = vaux(n,1)
      ENDDO
    ELSE
      DO n = 1,grid%mesh
        vpsloc(n) = (lloc*vaux(n,1)+(lloc+1.0_dp)*vaux(n,2)) / (2.0_dp*lloc + 1.0_dp)
      ENDDO
    ENDIF
  
  ENDIF ! lloc

  WRITE(*,*)
  WRITE(*,*) '**** EXIT my_pseudovloc'
  WRITE(*,*)

  RETURN
END SUBROUTINE


