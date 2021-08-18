!this SUBROUTINE calculates 
!the Coulombian operator on the optimal product basis
SUBROUTINE v_product()
  USE kinds, ONLY : DP
  USE wvfct, ONLY : npw,npwx,et
  USE mp, ONLY : mp_sum,mp_barrier
  USE klist, ONLY : nks,ngk,xk
  USE noncollin_module, ONLY: npol, noncolin
  USE mp_world, ONLY : world_comm,mpime
  USE spin_orb, ONLY: lspinorb
  USE io_global, ONLY : stdout, ionode
  USE input_simple
  USE gvect, ONLY : ngm, gstart,gg, g
  USE constants, ONLY : e2, fpi
  USE cell_base, ONLY: tpiba,tpiba2,omega,bg,at
  USE io_files,  ONLY : prefix, tmp_dir
  USE klist, ONLY : nks,xk
  USE fft_base, ONLY : dfftp, dffts
  USE fft_interfaces, ONLY : fwfft, invfft
  USE io_files, ONLY : prefix, tmp_dir, diropn
  USE wavefunctions, ONLY : psic
  USE polarization
  USE cell_base, ONLY : alat

  IMPLICIT NONE 

  
  INTEGER, EXTERNAL :: find_free_unit
  INTEGER :: ig,ii,iun,jj,kk
  REAL(kind=DP), ALLOCATABLE :: fac(:)
  REAL(kind=DP) :: qq
  COMPLEX(kind=DP), ALLOCATABLE :: prod_g(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: vmat(:,:)

  INTEGER :: ik,jk,ijk(3),ll
  REAL(kind=DP) :: qk(3),gk(3),sca

  INTEGER :: iunp
  LOGICAL :: exst
  REAL(kind=DP), ALLOCATABLE :: p_basis_r(:)
  COMPLEX(kind=DP), ALLOCATABLE :: p_basis(:,:)
  TYPE(polaw) :: pw
  INTEGER :: ix,iy,iz,ipol,iw
  INTEGER, PARAMETER :: n_int=20
  REAL(kind=DP) :: qx(3),qy(3),qz(3),qt(3),qq_fact
  COMPLEX(kind=DP), ALLOCATABLE :: amat(:,:),tmp_mat(:,:),p_mat(:,:)
  REAL(kind=DP), ALLOCATABLE :: facg(:)
  INTEGER, parameter :: n_int_loc = 20*50
  REAL(kind=DP) :: model
  INTEGER :: n_trovato

  WRITE(*,*)
  WRITE(stdout,*) 'Routine v_product'

  ALLOCATE(fac(npw_max*npol))
   
  IF(l_truncated_coulomb) THEN 
    DO ig=1,npw_max
      qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
      IF(qq > 1.d-8) THEN 
        fac(ig) = (e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
      ELSE
        fac(ig) = e2*fpi*(truncation_radius**2.d0/2.d0)
      ENDIF 
    ENDDO 
  ELSE
    DO ig=1,npw_max
      qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
      IF( qq > 1.d-8 ) THEN 
        fac(ig) = e2*fpi/(tpiba2*qq)
      ELSE
        fac(ig) = 0.d0
      ENDIF 
    ENDDO 
  ENDIF 

  fac=fac/omega/nks
  IF(npol>1) fac(npw_max+1:npw_max*npol)=fac(1:npw_max)

  ALLOCATE(prod_g(npw_max*npol,nprod_e))
  prod_g(1:npw_max*npol,1:nprod_e)=prod_e(1:npw_max*npol,1:nprod_e)
  DO ii=1,nprod_e
    prod_g(1:npw_max*npol,ii)=fac(1:npw_max*npol)*prod_g(1:npw_max*npol,ii)
  ENDDO 
  ALLOCATE(vmat(nprod_e,nprod_e))
  CALL ZGEMM('C','N',nprod_e,nprod_e,npw_max*npol,(1.d0,0.d0),prod_g, &
      npw_max*npol,prod_e,npw_max*npol,(0.d0,0.d0),vmat,nprod_e)
  CALL mp_sum(vmat,world_comm)
  
  !write v_mat on disk
  WRITE(*,*) 'Writing vmat to file'
  IF(ionode) THEN 
    iun=find_free_unit()
    open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.v_mat0', &
          status='unknown',form='unformatted')
    !product number  
    WRITE(iun) nprod_e
    DO ii=1,nprod_e
      WRITE(iun) vmat(1:nprod_e,ii)
    ENDDO 
  ENDIF 
  WRITE(*,*) 'End of writing vmat to file'

  IF(nks /= nkpoints(1)*nkpoints(2)*nkpoints(3))  THEN 
    WRITE(stdout,*) 'K points ill-defined'
    stop
  ENDIF 
   
  IF(ionode) WRITE(iun) nkpoints(1:3)
  DO ik=1,nks !k'
    DO jk=1,nks !k
      qk(1:3)=xk(1:3,ik)-xk(1:3,jk)
      DO ii=1,3
        sca = qk(1)*at(1,ii) + qk(2)*at(2,ii) + qk(3)*at(3,ii)
        sca = sca*nkpoints(ii)
        ijk(ii) = nint(sca)
        IF(ijk(ii)<0) ijk(ii) = ijk(ii) + nkpoints(ii)
        IF(ijk(ii)>= nkpoints(ii)) ijk(ii) = ijk(ii) - nkpoints(ii)
      ENDDO 
      IF(ionode) WRITE(iun) ijk(1:3)
    ENDDO 
  ENDDO 

  !obtain k'-k --> ijk table
  !write k'-k --> ijk table on disk
  !loop on i,j,k 
  n_trovato=0 ! n_found
  DO ii=0,nkpoints(1)-1
    DO jj=0,nkpoints(2)-1
      DO kk=0,nkpoints(3)-1
        
        WRITE(stdout,*) 'v_product: ', ii, jj, kk

        qk(1:3) = bg(1:3,1)*real(ii)/real(nkpoints(1)) + &
                  bg(1:3,2)*real(jj)/real(nkpoints(2)) + &
                  bg(1:3,3)*real(kk)/real(nkpoints(3))
        
        qk = -qk !(-) dovrebbe esser giusto da analisi G,k
        
        write(*,*) 'qk = ', qk
        
        DO ig=1,npw_max
          gk(1:3) = qk(1:3) + g(1:3,ig)
          qq = gk(1)**2.d0 + gk(2)**2.d0 + gk(3)**2.d0
          IF(qq > 1.d-8) THEN 
            !WRITE(*,*) 'qq not small'
            !fac(ig)=e2*fpi/(tpiba2*qq)
            fac(ig)=0.d0
            DO ix=-n_int+1,n_int
            DO iy=-n_int+1,n_int
            DO iz=-n_int+1,n_int
              qx(:) = 0.5d0*(1.d0/dble(n_int*nkpoints(1))*(dble(ix-1)) + &
                      0.5d0/dble(n_int*nkpoints(1)))*bg(:,1)
              qy(:) = 0.5d0*(1.d0/dble(n_int*nkpoints(2))*(dble(iy-1)) + &
                      0.5d0/dble(n_int*nkpoints(2)))*bg(:,2)
              qz(:) = 0.5d0*(1.d0/dble(n_int*nkpoints(3))*(dble(iz-1)) + &
                      0.5d0/dble(n_int*nkpoints(3)))*bg(:,3)
              qt(1:3) = qx(1:3) + qy(1:3) + qz(1:3) + gk(1:3)
              qq_fact = qt(1)**2 + qt(2)**2 + qt(3)**2
              fac(ig) = fac(ig) + 1.d0/qq_fact
            ENDDO 
            ENDDO 
            ENDDO 
            !
            fac(ig)=fac(ig)*e2*fpi/(8.d0*(dble(n_int))**3.d0)/tpiba2 
            !
          ELSE
            fac(ig)=0.d0
            !WRITE(*,'(1x,A,6F10.5)') 'TROVATO', qk(1:3),g(1:3,ig)
            WRITE(*,*) 'Handling small qq'
            n_trovato=mpime+1
            DO ix=-n_int_loc+1,n_int_loc
            DO iy=-n_int_loc+1,n_int_loc
            DO iz=-n_int_loc+1,n_int_loc
              qx(:) = 0.5d0*(1.d0/dble(n_int_loc*nkpoints(1))*(dble(ix-1)) + &
                      0.5d0/dble(n_int_loc*nkpoints(1)))*bg(:,1)
              qy(:) = 0.5d0*(1.d0/dble(n_int_loc*nkpoints(2))*(dble(iy-1)) + &
                      0.5d0/dble(n_int_loc*nkpoints(2)))*bg(:,2)
              qz(:) = 0.5d0*(1.d0/dble(n_int_loc*nkpoints(3))*(dble(iz-1)) + &
                      0.5d0/dble(n_int_loc*nkpoints(3)))*bg(:,3)
              qt(1:3) = qx(1:3) + qy(1:3) + qz(1:3) + qk(1:3) + g(1:3,ig)
              qq_fact = qt(1)**2 + qt(2)**2 + qt(3)**2
              fac(ig) = fac(ig) + 1.d0/qq_fact 
            ENDDO 
            ENDDO 
            ENDDO 
            WRITE(*,*) 'End of handling small qq'
            fac(ig)=fac(ig)*e2*fpi/(8.d0*(dble(n_int_loc))**3.d0)/tpiba2
          ENDIF 
        ENDDO 
        
        fac = fac/omega/nks
          
        IF( npol > 1) fac(npw_max+1:npw_max*npol)=fac(1:npw_max)

        prod_g(1:npw_max*npol,1:nprod_e)=prod_e(1:npw_max*npol,1:nprod_e)
        DO ll=1,nprod_e
          prod_g(1:npw_max*npol,ll)=fac(1:npw_max*npol)*prod_g(1:npw_max*npol,ll)
        ENDDO 
            
        CALL ZGEMM('C','N',nprod_e,nprod_e,npw_max*npol,(1.d0,0.d0),prod_g,npw_max*npol,&
                 &prod_e,npw_max*npol,(0.d0,0.d0),vmat,nprod_e)
        CALL mp_sum(vmat,world_comm)
        IF(ionode) THEN 
          DO ll=1,nprod_e
            WRITE(iun) vmat(1:nprod_e,ll)
          ENDDO 
        ENDIF 

      ENDDO 
    ENDDO 
  ENDDO 

  CALL mp_sum(n_trovato, world_comm)
  !WRITE(stdout,*) 'MPI TROVATO', n_trovato

  IF(w_type==1) THEN 

    WRITE(*,*) 'wtype = 1'
    
    DO ii=0,nkpoints(1)-1
    DO jj=0,nkpoints(2)-1
    DO kk=0,nkpoints(3)-1
      qk(1:3) = bg(1:3,1)*real(ii)/real(nkpoints(1)) + &
                bg(1:3,2)*real(jj)/real(nkpoints(2)) + &
                bg(1:3,3)*real(kk)/real(nkpoints(3))
      qk = -qk !era -
      DO ig=1,npw_max
        gk(1:3) = qk(1:3) + g(1:3,ig)
        qq = gk(1)**2.d0 + gk(2)**2.d0 + gk(3)**2.d0
        !HERE
        qq = qq*tpiba2/0.529**2.d0
        model = (1.d0/epsm-1.d0)*exp(-3.1415926d0*qq/2.d0/lambdam**2.d0)
        !
        IF( qq > 1.d-8 ) THEN 
          fac(ig) = 0.d0
          DO ix = -n_int+1,n_int
          DO iy = -n_int+1,n_int
          DO iz = -n_int+1,n_int
            qx(:) = 0.5d0*(1.d0/dble(n_int*nkpoints(1))*(dble(ix-1)) + &
                    0.5d0/dble(n_int*nkpoints(1)))*bg(:,1)
            qy(:) = 0.5d0*(1.d0/dble(n_int*nkpoints(2))*(dble(iy-1)) + &
                    0.5d0/dble(n_int*nkpoints(2)))*bg(:,2)
            qz(:) = 0.5d0*(1.d0/dble(n_int*nkpoints(3))*(dble(iz-1)) + &
                    0.5d0/dble(n_int*nkpoints(3)))*bg(:,3)
            qt(1:3) = qx(1:3) + qy(1:3) + qz(1:3) + gk(1:3)
            qq_fact = qt(1)**2 + qt(2)**2 + qt(3)**2
            fac(ig)=fac(ig)+1.d0/qq_fact
          ENDDO 
          ENDDO 
          ENDDO 
          fac(ig)=fac(ig)*model
        
        ELSE
      
          fac(ig)=0.d0
          DO ix = -n_int_loc+1,n_int_loc
          DO iy = -n_int_loc+1,n_int_loc
          DO iz = -n_int_loc+1,n_int_loc
            qx(:) = 0.5d0*(1.d0/dble(n_int_loc*nkpoints(1))*(dble(ix-1)) + &
                    0.5d0/dble(n_int_loc*nkpoints(1)))*bg(:,1)
            qy(:) = 0.5d0*(1.d0/dble(n_int_loc*nkpoints(2))*(dble(iy-1)) + &
                    0.5d0/dble(n_int_loc*nkpoints(2)))*bg(:,2)
            qz(:) = 0.5d0*(1.d0/dble(n_int_loc*nkpoints(3))*(dble(iz-1)) + &
                    0.5d0/dble(n_int_loc*nkpoints(3)))*bg(:,3)
            qt(1:3)=qx(1:3)+qy(1:3)+qz(1:3)+qk(1:3)+g(1:3,ig)
            qq_fact=qt(1)**2+qt(2)**2+qt(3)**2
            fac(ig)=fac(ig)+1.d0/qq_fact
          ENDDO 
          ENDDO 
          ENDDO 
          fac(ig)=fac(ig)*e2*fpi/(8.d0*(dble(n_int_loc))**3.d0)/tpiba2
          fac(ig)=fac(ig)*model
        
        ENDIF 
      
      ENDDO 

      fac=fac/omega/nks
      IF(npol>1) fac(npw_max+1:npw_max*npol)=fac(1:npw_max)
      prod_g(1:npw_max*npol,1:nprod_e)=prod_e(1:npw_max*npol,1:nprod_e)
      
      DO ll=1,nprod_e
        prod_g(1:npw_max*npol,ll)=fac(1:npw_max*npol)*prod_g(1:npw_max*npol,ll)
      ENDDO 
      
      CALL ZGEMM('C','N',nprod_e,nprod_e,npw_max*npol,(1.d0,0.d0),prod_g,npw_max*npol,&
           prod_e,npw_max*npol,(0.d0,0.d0),vmat,nprod_e)
      CALL mp_sum(vmat,world_comm)
      
      IF(ionode) THEN 
        DO ll=1,nprod_e
          WRITE(iun) vmat(1:nprod_e,ll)
        ENDDO 
      ENDIF 
    
    ENDDO ! kk
    ENDDO ! jj
    ENDDO ! ii

  ELSE 

    !
    ! This part need GWW
    !

    !calculate v(i,j,k) matrix
    !write it on disk

    ! part relative to W_c
    ! read in polarisability basis and put it on G space
    iunp=find_free_unit()
    ALLOCATE(p_basis_r(dfftp%nnr))
    ALLOCATE(p_basis(npw_max,numpw))
    CALL diropn( iunp, 'basis2simple',dfftp%nnr, exst )
    WRITE(stdout,*) 'exst = ', exst
    DO ii=1,numpw
      CALL davcio(p_basis_r,dfftp%nnr,iunp,ii,-1)
      !DO fft
      psic(1:dfftp%nnr)=p_basis_r(1:dfftp%nnr)
      CALL fwfft('Rho', psic, dfftp)
      DO ig=1,npw_max
        p_basis(ig,ii)=psic(dfftp%nl(ig))
      ENDDO 
    ENDDO 
    close(iunp)
    DEALLOCATE(p_basis_r)

    !read in irreducible polarisability operator
    !read P
    CALL initialize_polaw(pw)
    CALL read_polaw_global(0, pw)
    ALLOCATE(amat(nprod_e,numpw),p_mat(numpw,numpw))
    ALLOCATE(tmp_mat(nprod_e,numpw))
      
    p_mat(1:numpw,1:numpw)=pw%pw(1:numpw,1:numpw)
    !calculate |G|
    ALLOCATE(facg(npw_max*npol))
   
    DO ig=1,npw_max
      gk(1:3) = g(1:3,ig)
      qq = gk(1)**2.d0 + gk(2)**2.d0 + gk(3)**2.d0
      IF(qq > 1.d-8) THEN 
        !facg(ig)=e2*fpi/tpiba2*qq
        facg(ig)=0
        DO ix=-n_int+1,n_int
        DO iy=-n_int+1,n_int
        DO iz=-n_int+1,n_int           
          qx(:)=0.5d0*(1.d0/dble(n_int)*(dble(ix-1))+0.5d0/dble(n_int))*bg(:,1)
          qy(:)=0.5d0*(1.d0/dble(n_int)*(dble(iy-1))+0.5d0/dble(n_int))*bg(:,2)
          qz(:)=0.5d0*(1.d0/dble(n_int)*(dble(iz-1))+0.5d0/dble(n_int))*bg(:,3)
          qt(:)=qx(1:3)+qy(1:3)+qz(1:3)+g(1:3,ig)
          qq_fact=qt(1)**2+qt(2)**2+qt(3)**2
          facg(ig)=facg(ig)+1.d0/qq_fact
        ENDDO 
        ENDDO 
        ENDDO
        ! 
        facg(ig)=facg(ig)*e2*fpi/(8.d0*(dble(n_int))**3.d0)/tpiba
        !
      ELSE
        !DO integrate                                   
        facg(ig)=0.d0 
        DO ix=-n_int_loc+1,n_int_loc
        DO iy=-n_int_loc+1,n_int_loc
        DO iz=-n_int_loc+1,n_int_loc             
          qx(:) = 0.5d0*(1.d0/dble(n_int_loc)*(dble(ix-1))+0.5d0/dble(n_int_loc))*bg(:,1)
          qy(:) = 0.5d0*(1.d0/dble(n_int_loc)*(dble(iy-1))+0.5d0/dble(n_int_loc))*bg(:,2)
          qz(:) = 0.5d0*(1.d0/dble(n_int_loc)*(dble(iz-1))+0.5d0/dble(n_int_loc))*bg(:,3)
          qt(:)=qx(1:3)+qy(1:3)+qz(1:3)+g(1:3,ig)
          qq_fact=qt(1)**2+qt(2)**2+qt(3)**2
          facg(ig)=facg(ig)+1.d0/qq_fact    
        ENDDO 
        ENDDO 
        ENDDO 
        !
        facg(ig)=facg(ig)*e2*fpi/(8.d0*(dble(n_int_loc))**3.d0)/tpiba2
      ENDIF 
    
    ENDDO 
    
    IF(npol>1) facg(1+(npol-1)*npw_max:npol*npw_max)=facg(1:npw_max)

    facg=facg/omega
    
    !loop on q
    !calculate fac terms and in case DO integrate
    DO ii=0,nkpoints(1)-1
    DO jj=0,nkpoints(2)-1
    DO kk=0,nkpoints(3)-1
      qk(1:3) = bg(1:3,1)*real(ii)/real(nkpoints(1)) + &
                bg(1:3,2)*real(jj)/real(nkpoints(2)) + &
                bg(1:3,3)*real(kk)/real(nkpoints(3))            
      qk = -qk !!- dovrebbe esser giusto da analisi G,k       
      DO ig=1,npw_max
        gk(1:3) = qk(1:3) + g(1:3,ig)
        qq = gk(1)**2.d0 + gk(2)**2.d0 + gk(3)**2.d0
        IF(qq > 1.d-8) THEN 
          !fac(ig)=e2*fpi/(tpiba2*qq)
          fac(ig)=0.d0
          DO ix = -n_int+1,n_int
          DO iy = -n_int+1,n_int
          DO iz = -n_int+1,n_int
            qx(:) = 0.5d0*(1.d0/dble(n_int*nkpoints(1))*(dble(ix-1)) + &
                    0.5d0/dble(n_int*nkpoints(1)))*bg(:,1)
            qy(:) = 0.5d0*(1.d0/dble(n_int*nkpoints(2))*(dble(iy-1)) + &
                    0.5d0/dble(n_int*nkpoints(2)))*bg(:,2)
            qz(:) = 0.5d0*(1.d0/dble(n_int*nkpoints(3))*(dble(iz-1)) + &
                    0.5d0/dble(n_int*nkpoints(3)))*bg(:,3)
            qt(:) = qx(1:3) + qy(1:3) + qz(1:3) + gk(1:3)
            qq_fact = qt(1)**2 + qt(2)**2 + qt(3)**2
            fac(ig)=fac(ig)+1.d0/qq_fact
          ENDDO 
          ENDDO 
          ENDDO            
          !
          fac(ig)=fac(ig)*e2*fpi/(8.d0*(dble(n_int))**3.d0)/tpiba2          
          !
        ELSE 
          !DO integrate
          fac(ig) = 0.d0
          !WRITE(stdout,*) ' TROVATO',qk(1:3),g(1:3,ig)
          DO ix=-n_int_loc+1,n_int_loc
          DO iy=-n_int_loc+1,n_int_loc
          DO iz=-n_int_loc+1,n_int_loc                    
            qx(:) = 0.5d0*(1.d0/dble(n_int_loc*nkpoints(1))*(dble(ix-1)) + &
                    0.5d0/dble(n_int_loc*nkpoints(1)))*bg(:,1)
            qy(:) = 0.5d0*(1.d0/dble(n_int_loc*nkpoints(2))*(dble(iy-1)) + &
                    0.5d0/dble(n_int_loc*nkpoints(2)))*bg(:,2)
            qz(:) = 0.5d0*(1.d0/dble(n_int_loc*nkpoints(3))*(dble(iz-1)) + &
                    0.5d0/dble(n_int_loc*nkpoints(3)))*bg(:,3)
            qt(:) = qx(1:3) + qy(1:3) + qz(1:3) + qk(1:3) + g(1:3,ig)
            qq_fact=qt(1)**2+qt(2)**2+qt(3)**2
            fac(ig)=fac(ig)+1.d0/qq_fact
          ENDDO 
          ENDDO 
          ENDDO 
          !
          fac(ig)=fac(ig)*e2*fpi/(8.d0*(dble(n_int_loc))**3.d0)/tpiba2 
        ENDIF 
      ENDDO 
      fac=fac/omega/nks
      fac=fac*sqrt(facg)/sqrt(fac)
            
      !calculate <\mathcal{E_\alpha}|(V(q)\Phi_\mu>
      !calcuate W_c(q)_\alpha\beta
      vmat=0.d0
      DO ipol=0,npol-1            
        DO iw=1,nprod_e
          prod_g(1:npw_max,iw)=fac(1:npw_max)*prod_e(1+ipol*npw_max:npw_max+npw_max*ipol,iw)
        ENDDO 
        CALL ZGEMM('C','N',nprod_e,numpw,npw_max,(1.d0,0.d0),prod_g,npw_max*npol,p_basis,npw_max,(0.d0,0.d0),amat,nprod_e)
        CALL mp_sum(amat, world_comm)
        CALL ZGEMM('N','N',nprod_e,numpw,numpw,(1.d0,0.d0),amat,nprod_e,p_mat,numpw,(0.d0,0.d0),tmp_mat,nprod_e)
        CALL ZGEMM('N','C',nprod_e,nprod_e,numpw,(1.d0,0.d0),tmp_mat,nprod_e,amat,nprod_e,(1.d0,0.d0),vmat,nprod_e)
      ENDDO 
      
      !writes on disk
      IF(ionode) THEN 
        DO ll=1,nprod_e
           WRITE(iun) vmat(1:nprod_e,ll)
        ENDDO 
      ENDIF 
    ENDDO 
    ENDDO 
    ENDDO 

    DEALLOCATE(p_basis)
    DEALLOCATE(amat)
    DEALLOCATE(tmp_mat)
    DEALLOCATE(facg)
    CALL free_memory_polaw(pw)
  ENDIF 

  IF(ionode) close(iun)
  DEALLOCATE(vmat)
  DEALLOCATE(prod_g)
  DEALLOCATE(fac)
  WRITE(stdout,*) 'Out of v_product'

END SUBROUTINE  v_product
