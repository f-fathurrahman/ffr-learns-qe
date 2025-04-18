
!--------------------------------------------------------------------------
subroutine my_gener_pseudo()
!--------------------------------------------------------------------------
  !
  ! This routine generate a pseudopotential in separable form
  ! It can be of NC type or of US type
  ! Multiple projections are allowed.
  ! Spin-orbit split pseudopotentials are also available.
  ! NB: bmat indices are as in the Vanderbilt paper PRB (1990)
  !
  ! The output of the routine are:
  !
  ! phis: the pseudo wavefunctions
  ! betas: the nonlocal projectors
  ! bmat:  the pseudopotential coefficients
  ! qq:    the integrals of the q functions
  ! qvan:  the augmentation functions
  ! vpsloc: the local pseudopotential
  ! chis:   auxiliary functions
  !
  ! The construction of a PAW dataset can also be done (experimental)
  !      
  use kinds, only: dp
  use radial_grids, only: ndmx
  use ld1_parameters, only: nwfsx
  use io_global, only : stdout
  use mp,        only : mp_bcast
  use ld1inc, only: grid, lls, jjs, new, psipaw, psi, ecutwfc, ecutrho, &
                    psipsus, tm, ocs, phis, els, nwfs, nspin, rel, nlcc, &
                    file_chi, file_beta, file_wfcncgen, file_qvan, &
                    file_wfcusgen, file_wfcaegen, psccharge, aeccharge, &
                    qvan, qvanl, qq, bmat, ddd, betas, nbeta, ikk, pseudotype, &
                    pawsetup, zval, vpsloc, vpot, vnl, lpaw, rcloc, rcutus, &
                    enl, enls, rcut, chis, nstoae, rmatch_augfun,&
                    lnc2paw, rhos, which_augfun, psipaw_rel, &
                    cau_fact, ik, ikus
  use atomic_paw, only : us2paw, paw2us
  use matrix_inversion
  implicit none

  integer ::   &
       ikloc, &  ! the point corresponding to rc local
       ns,    &  ! counter on pseudo functions
       ns1,   &  ! counter on pseudo functions
       nnode, &  ! the number of nodes of phi
       lam,   &  ! the angular momentum
       irc       ! the maximum radius of the sphere

  real(DP) ::    &
       xc(8),         &  ! parameters of bessel functions
       psi_in(ndmx),  &  ! the all_electron wavefunction
       gi(ndmx),      &  ! auxiliary to compute the integrals
       occ, norm1, norm2, &  !
       speed        ! an estimate of electron speed at rcutus

  real(DP), allocatable :: &
       b(:,:), binv(:,:) ! the B matrix and its inverse

  real(DP) ::    &
       aekin(nwfsx,nwfsx),  & ! AE kinetic energies
       pskin(nwfsx,nwfsx),  & ! PS kinetic energies
       kindiff(nwfsx,nwfsx)   ! AE-PS k.e.

  real(DP), external :: int_0_inf_dr    ! the function calculating the integral 

  character(len=6), external :: int_to_char

  integer :: n, nwf0, nst, ikl, is, ind

  logical :: lbes4     ! use 4 Bessel functions expansion

  ! additional vars for paw
  real(DP) :: vpotpaw (ndmx) ! total potential to be used for PAW 
                             ! generation (normally the AE potential)


  WRITE(*,*)
  WRITE(*,*) '***** ENTER my_gener_pseudo *****'
  WRITE(*,*)


  if (lpaw) then
    write(stdout, '(/,5x,21(''-''),'' Generating PAW atomic setup '',20(''-''),/)')
  elseif (pseudotype == 1.or.pseudotype == 2) then
    write(stdout, '(/,5x,21(''-''),'' Generating NC pseudopotential '',21(''-''),/)')
  elseif (pseudotype == 3) then
    write(stdout, '(/,5x,21(''-''),'' Generating US pseudopotential '',21(''-''),/)')
  else
    call errore('gener_pseudo', 'pseudotype not programmed',1)
  endif

  if( pseudotype == 1.and.rel == 2 ) then
    call errore('gener_pseudo', 'not programmed', 2)
  endif

  if( pseudotype /= 3.and. lpaw ) then
    call errore('gener_pseudo', 'please start from a US for generating a PAW dataset' ,pseudotype)
  endif

  psipaw = 0.0_dp
  phis = 0.0_dp
  !
  ! compute the local potential from the all-electron potential
  !
  call my_pseudovloc()
  !
  ! initialize total potential for PAW generation
  if( lpaw ) then
    if( .not. lnc2paw ) then
      write(*,*) 'Not using lnc2paw, setting vpotpaw to vpot'
      vpotpaw(1:grid%mesh) = vpot(1:grid%mesh,1)
    else
      vpotpaw(1:grid%mesh) = vpsloc(1:grid%mesh)
    endif
  endif
  !
  ! if nlcc is true compute here the core charge
  ! the core charge is needed also for the PAW dataset
  !
  if( nlcc .or. lpaw ) then
    call my_set_rho_core()
  endif
  !
  ! set the appropriate energies and the correspondence all-electron
  ! pseudo
  !
  do n = 1,nwfs
    if( enls(n) == 0.0_dp ) then
      enls(n) = enl(nstoae(n))
    endif
  enddo
  !
  ik = 0
  ikus = 0
  ikloc = 0
  ! find ik, ikus, and ikloc
  do ns = 1,nbeta
    !
    do n=1,grid%mesh
      !
      if( grid%r(n) < rcut(ns) ) then
        ik(ns) = n
      endif
      !
      if( grid%r(n) < rcutus(ns) ) then
        ikus(ns) = n
      endif
      !
      if( grid%r(n) < rcloc) then
        ikloc = n
      endif
    enddo
    ! make ik odd
    if( mod(ik(ns),2) == 0) then
      ik(ns) = ik(ns) + 1
    endif
    ! make ikus odd
    if( mod(ikus(ns),2) == 0) then
      ikus(ns) = ikus(ns) + 1
    endif
    ! make ikloc odd
    if( mod(ikloc,2) == 0) then
      ikloc = ikloc + 1
    endif
    ! check validity of ik, ikus, and ikloc
    if( ik(ns) > grid%mesh) then
      call errore('gener_pseudo','ik is wrong ',ns)
    endif
    ! why?
    if( ikus(ns) + 10 > grid%mesh) then
      call errore('gener_pseudo','ikus is wrong ',ns)
    endif
    ! why?
    if( ikloc + 5 > grid%mesh) then
      call errore('gener_pseudo','ikloc is wrong ',ns)
    endif
    !
    if( pseudotype == 3 ) then
      ikk(ns) = max(ikus(ns)+10, ikloc+5)
    else
      ikk(ns) = max(ik(ns)+10, ikloc+5)
    endif
    !
    if( ikk(ns) > grid%mesh) then
      call errore('gener_pseudo', 'ikk is wrong ',ns)
    endif
  enddo ! ns
  write(*,*) 'nbeta = ', nbeta
  write(*,*) 'ik = ', ik(1:nbeta)
  write(*,*) 'ikus = ', ikus(1:nbeta)
  write(*,*) 'ikloc = ', ikloc
  !
  do ns = 1,nbeta
    do ns1 = 1,nbeta
      if( lls(ns) == lls(ns1) .and. ikk(ns1) > ikk(ns)) then
        ikk(ns) = ikk(ns1)
      endif
    enddo
  enddo
  !
  irc = maxval(ikk(1:nbeta)) + 8
  IF(mod(irc,2) == 0) then
    irc = irc + 1
  endif
  !
  IF(irc > grid%mesh) then
    CALL errore('gener_pseudo','irc is too large',1)
  endif
  !
  ! Set the all-electron wavefunctions, calculating those at user supplied
  ! energies. The wavefunctions are written on file at this point, so
  ! the user can check them also when the pseudopotential generation is  
  ! unsuccessful
  ! 
  IF (rel==2) WRITE(stdout,'(/,5x,"Fully relativistic calculation: |psi|^2")')
  !
  do ns = 1,nbeta
    nwf0 = nstoae(ns)
    if( new(ns) ) then
      write(*,*) 'new(ns) is true'
      call my_set_psi_in( ikus(ns), lls(ns), jjs(ns), enls(ns), psipaw(1,ns), psipaw_rel(1,ns) )
    else
      lam = lls(ns)
      nst = (lam+1)*2
      psipaw(:,ns) = psi(:,1,nwf0)
      !
      do n = 1,grid%mesh
        gi(n) = psipaw(n,ns)*psipaw(n,ns)
      enddo
      !
      norm1 = sqrt(int_0_inf_dr(gi,grid,grid%mesh,nst))
      !
      IF( rel==2 ) THEN
        psipaw_rel(:,ns)=psi(:,2,nwf0)
        DO n = 1,irc
          gi(n) = psipaw_rel(n,ns)*psipaw_rel(n,ns)
        ENDDO
        norm2 = sqrt(int_0_inf_dr(gi,grid,irc,nst))
        WRITE(stdout,'(/,5x,"Wfc ",a2," LC norm =",f12.8," SC norm ="&
                      &,f12.8," missing",f12.8)') els(ns), norm1**2, &
                      norm2**2, 1.0_DP-(norm1**2+norm2**2)
        norm1 = sqrt(norm1**2+norm2**2)
        psipaw_rel(:,ns)=psipaw_rel(:,ns)/norm1
        IF( lpaw ) THEN
          speed = psipaw(ikus(ns),ns)*(enls(ns)-vpotpaw(ikus(ns)))*psipaw(ikus(ns),ns)
          WRITE(stdout,'(5x,"Speed at rcutus ",f10.3, " a.u.,  &
               &v/c =",f10.6,",  (v/c)^2 =",f12.8 )') sqrt(abs(speed)), &
                            sqrt(abs(speed)/cau_fact**2),speed/cau_fact**2
        ENDIF
      ENDIF
      ! normalize
      psipaw(:,ns) = psipaw(:,ns)/norm1
     ENDIF
  ENDDO
  !
  call write_wfcfile(file_wfcaegen,psipaw,els,nwfs)
  !
  ! compute the pseudowavefunctions by expansion in spherical
  ! bessel function before r_c
  !
  ecutrho = 0.0_dp
  ecutwfc = 0.0_dp
  do ns = 1,nbeta
    lam = lls(ns)
    nst = (lam + 1)*2
    nwf0 = nstoae(ns)
    !
    if( new(ns) ) then
      occ = 1.0_DP
    else
      occ = ocs(ns)
    endif 
    !
    ! save the all-electron function for the PAW setup
    !
    psi_in(1:grid%mesh) = psipaw(1:grid%mesh,ns) 
    !
    ! compute the phi functions
    !
    write(*,*) 'lnc2paw = ', lnc2paw
    !
    if( lpaw .and. lnc2paw ) then
      ! first compute possibly harder NC pseudowfcs to be
      ! used as AE reference for PAW generation
      nnode = 0
      if( tm ) then 
        write(*,*) 'Calling compute_phi_tm'
        call compute_phi_tm(lam, ik(ns), psi_in, phis(1,ns), 1, xc, enls(ns), els(ns) )
      else
        write(*,*) 'Calling compute_phi'
        call compute_phi(lam, ik(ns), psi_in, phis(1,ns), xc, 1, occ, enls(ns), els(ns) )
      endif
      psipaw(1:grid%mesh,ns) = phis(1:grid%mesh,ns)
    endif
    !
    IF( (which_augfun=='PSQ' .AND. ik(ns) /= ikus(ns)) .OR. &
         (lpaw .AND. .NOT. lnc2paw) ) THEN
      psipsus(:,ns)=psi_in(:)
    ELSE
      if( tm ) then
        call compute_phi_tm(lam,ik(ns), psi_in,phis(1,ns), 1, xc, enls(ns), els(ns))
      else
        call compute_phi(lam,ik(ns), psi_in, phis(1,ns),xc,1,occ,enls(ns),els(ns))
        ecutrho = max(ecutrho, 8.0_dp*xc(6)**2 )
      endif
      !
      ! US only on the components where ikus <> ik
      ! 
      psipsus(:,ns) = phis(:,ns) 
    ENDIF
    !
    if( ikus(ns) /= ik(ns) ) then
      write(*,*) 'Calling compute_phius'
      call compute_phius(lam,ikus(ns),psipsus(1,ns),phis(1,ns),xc,1,els(ns))
      ecutwfc = max(ecutwfc,2.0_dp*xc(5)**2)
      lbes4 = .true.
    else
      lbes4 = .false.
      if(.not. tm) then
        ecutwfc = max(ecutwfc, 2.0_dp*xc(6)**2)
      endif
    endif
    !
    if(tm .and. ik(ns)==ikus(ns) ) then
      write(*,*) 'Calling compute_chi_tm'
      call compute_chi_tm(lam,ik(ns),ikk(ns),phis(1,ns),chis(1,ns),xc,enls(ns))
    else
      write(*,*) 'Calling compute_chi'
      call compute_chi(lam, ikk(ns), phis(1,ns), chis(1,ns), xc, enls(ns), lbes4)
    endif
  enddo
  !
  ! for each angular momentum take the same integration point
  !
  !
  ! construct B_{ij}
  !
  bmat = 0.0_dp
  do ns = 1,nbeta
    do ns1 = 1,nbeta
      if( lls(ns) == lls(ns1) .and. abs(jjs(ns)-jjs(ns1)) < 1.e-7_dp ) then
        nst = (lls(ns)+1)*2
        ikl = ikk(ns1)
        do n=1,grid%mesh
          gi(n) = phis(n,ns)*chis(n,ns1)
        enddo
        bmat(ns,ns1) = int_0_inf_dr(gi,grid,ikl,nst)
      endif
    enddo
  enddo

  allocate( b(nbeta, nbeta), binv(nbeta, nbeta) )

  if( pseudotype == 1 ) then
    !
    ! NC single-projector PP: construct the semilocal potential 
    !
    vnl = 0.0_dp
    do ns = 1,nbeta
      !
      lam = lls(ns)
      !
      if( rel < 2 .or. lls(ns) == 0 .or. abs(jjs(ns)-lls(ns)+0.5_dp) < 0.001_dp ) then
        ind = 1
      elseif( rel == 2 .and. lls(ns) > 0 .and. abs(jjs(ns)-lls(ns)-0.5_dp) < 0.001_dp ) then
        ind = 2
      endif
      !XXX ind may be uninitialized here
      !
      do n = 1,ikk(ns)
        vnl(n,lam,ind) = chis(n,ns)/phis(n,ns)
      enddo
    enddo
    !
  elseif( pseudotype == 2 ) then
    !
    ! symmetrize the B matrix
    !
    do ns = 1,nbeta
      do ns1 = 1,ns-1
        bmat(ns,ns1) = 0.5_dp*(bmat(ns,ns1) + bmat(ns1,ns))
        bmat(ns1,ns) = bmat(ns,ns1)
      enddo
    enddo
  endif
  !
  do ns = 1,nbeta
    do ns1 = 1,nbeta
      b(ns,ns1) = bmat(ns,ns1)
    enddo
  enddo
  !
  ! compute the inverse of the matrix B_{ij}:  B_{ij}^-1
  !
  write(stdout,'(/5x,'' The bmat matrix'')')
  do ns1=1,nbeta
    write(stdout,'(6f12.5)') (bmat(ns1,ns),ns=1,nbeta)
  enddo
  if( nbeta > 0 ) then
    call invmat(nbeta, b, binv)
  endif
  !
  ! compute the beta functions
  !
  betas = 0.0_dp
  do ns = 1,nbeta
    do ns1 = 1,nbeta
      do n = 1,grid%mesh
        betas(n,ns) = betas(n,ns) + binv(ns1,ns)*chis(n,ns1)
      enddo
    enddo
  enddo
  deallocate(b, binv)
  !
  qq = 0.0_dp
  write(*,*) 'pseudotype = ', pseudotype
  if( pseudotype == 3 ) then
    !
    ! compute the Q functions
    !
    do ns = 1,nbeta
      do ns1 = 1,ns
        ikl = max(ikk(ns),ikk(ns1))
        if( which_augfun=='PSQ' .and. rel==2) then
          do n=1, ikl
            qvan(n,ns,ns1) = psipsus(n,ns) * psipsus(n,ns1) &
                         & + psipaw_rel(n,ns) * psipaw_rel(n,ns1) &
                         & - phis(n,ns) * phis(n,ns1)
            gi(n) = qvan(n,ns,ns1)
          enddo
        else
          do n=1, ikl
            qvan(n,ns,ns1) = psipsus(n,ns) * psipsus(n,ns1) - phis(n,ns) * phis(n,ns1)
            gi(n) = qvan(n,ns,ns1)
          enddo
        endif
        !
        do n = ikl+1,grid%mesh
          qvan(n,ns,ns1) = 0.0_dp
        enddo
        !
        ! and puts its integral in qq
        !
        if( lls(ns) == lls(ns1) .and. abs(jjs(ns)-jjs(ns1)) < 1.e-8_dp) then
          nst = (lls(ns)+1)*2
          qq(ns,ns1) = int_0_inf_dr(gi,grid,ikl,nst)
        endif
        !
        ! set the bmat with the eigenvalue part
        !
        bmat(ns,ns1) = bmat(ns,ns1)+enls(ns1)*qq(ns,ns1)
        !
        ! Use symmetry of the n,ns1 indices to set qvan and qq and bmat
        !
        if(ns /= ns1) then
          do n = 1,grid%mesh
            qvan(n,ns1,ns) = qvan(n,ns,ns1)
          enddo
          qq(ns1,ns) = qq(ns,ns1)
          bmat(ns1,ns) = bmat(ns1,ns) + enls(ns)*qq(ns1,ns)
        endif
        !
      enddo ! ns1: 1..ns
      !
    enddo ! ns: 1..beta
    !
    write(stdout,'(/5x,'' The bmat + epsilon qq matrix'')')
    do ns1=1,nbeta
      write(stdout,'(6f12.5)') (bmat(ns1,ns),ns=1,nbeta)
    enddo
    write(stdout,'(/5x,'' The qq matrix'')')
    do ns1 = 1,nbeta
      write(stdout,'(6f12.5)') (qq(ns1,ns),ns=1,nbeta)
    enddo
  endif
  !
  do is=1,nspin
    ddd(:,:,is) = bmat(:,:)
  enddo
  !
  ! generate a PAW dataset if required
  !
  if( lpaw ) then
    if (lnc2paw) then
      write (stdout,'(/5x,''WARNING: __PAW_FROM_NC__'')')
    endif
    !
    !symbol=atom_name(nint(zed))
    !
    ! compute kinetic energy differences, using:
    ! AE:   T |psi> = (e - Vae) |psi>
    ! PS:   T |phi> = (e - Vps) |phi> - |chi>
    do ns=1,nbeta
      do ns1=1,ns
        if( lls(ns)==lls(ns1) .and. jjs(ns)==jjs(ns1) ) then
          !
          ikl = max(ikk(ns),ikk(ns1))
          nst = 2*(lls(ns)+1)
          !
          do n=1,ikl
            gi(n) = psipaw(n,ns)*(enls(ns1) - vpotpaw(n))*psipaw(n,ns1)
          enddo
          !
          aekin(ns,ns1) = int_0_inf_dr(gi(1:grid%mesh),grid,ikl,nst)
          !
          IF (rel==2) THEN
            do n=1,irc
              gi(n) = psipaw_rel(n,ns)*(enls(ns1)-vpotpaw(n))*psipaw_rel(n,ns1)
            enddo
            aekin(ns,ns1) = aekin(ns,ns1) + int_0_inf_dr(gi(1:grid%mesh),grid,irc,nst)
          ENDIF
          !
          do n=1,ikl
            gi(n)=phis(n,ns)*( (enls(ns1)-vpsloc(n))*phis(n,ns1) - chis(n,ns1) )
          end do
          pskin(ns,ns1)=int_0_inf_dr(gi(1:grid%mesh),grid,ikl,nst)
        !
        else
          ! zeros
          aekin(ns,ns1)=0._dp
          pskin(ns,ns1)=0._dp
        endif
        !
        kindiff(ns,ns1) = aekin(ns,ns1) - pskin(ns,ns1)
        kindiff(ns1,ns) = aekin(ns,ns1) - pskin(ns,ns1)
      enddo
    end do
    !
    ! create the 'pawsetup' object containing the atomic setup for PAW
    write(*,*) 'Calling us2paw'
    call us2paw( pawsetup, &
         & zval, grid, rmatch_augfun, ikk,  &
         & nbeta, lls, jjs, ocs, enls, els, rcutus, psipaw, psipaw_rel, &
         & phis, betas, qvan, kindiff, nlcc, aeccharge, psccharge, vpotpaw, &
         & vpsloc, which_augfun, rel )
    !
    ! reread augmentation functions and descreened potentials from PAW
    write(*,*) 'Calling paw2us'
    call paw2us( pawsetup, zval, grid, nbeta, lls, jjs, ikk, betas, &
               & qq, qvan, vpsloc, bmat, rhos, els, rcutus, pseudotype, &
               & psipaw_rel )
    !
  else
    write(*,*) 'Calling pseudo_q and descreening'
    !
    !  Pseudize the Q functions if required. This might be needed for
    !  pseudo-potentials with semicore states. In this case the cut-off radius
    !  for the norm conserving wavefunctions is quite small and without
    !  the Q pseudization the augmentation charges are very hard making the
    !  ASR in phonon calculation very difficult to converge.
    ! 
    IF (which_augfun=='PSQ') CALL pseudo_q(qvan,qvanl)
    !
    ! unscreen the local potential and the D coefficients
    !
    call descreening()
  endif
  !
  ! write the main functions on files
  !
  ! The beta functions  
  !
  call write_wfcfile(file_beta,betas,els,nbeta)
  call write_wfcfile_ft(file_beta,betas,nbeta)
  !
  ! The chi functions
  !
  call write_wfcfile(file_chi,chis,els,nbeta)
  call write_wfcfile_ft(file_chi,chis,nbeta)
  !
  ! The augmentations functions
  !
  if (file_qvan /= ' ') then
    do ns1=1,nbeta
      call write_wfcfile(TRIM(file_qvan)//TRIM( int_to_char(ns1)), qvan(1,1,ns1), els, ns1)
    enddo
  endif
  !
  ! The norm conserving wavefunctions
  !
  call write_wfcfile(file_wfcncgen,psipsus,els,nwfs)
  call write_wfcfile_ft(file_wfcncgen,psipsus,nwfs)
  !
  ! The us wavefunctions
  !
  call write_wfcfile(file_wfcusgen,phis,els,nwfs)
  call write_wfcfile_ft(file_wfcusgen,phis,nwfs)

  write(stdout,"(/,5x,19('-'),' End of pseudopotential generation ',19('-'),/)")

  WRITE(*,*)
  WRITE(*,*) '***** EXIT my_gener_pseudo *****'
  WRITE(*,*)

  return
end subroutine
