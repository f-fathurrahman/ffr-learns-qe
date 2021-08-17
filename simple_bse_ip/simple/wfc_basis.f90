!this SUBROUTINE creates a optimal basis for the periodic functions |u_{nk}>

SUBROUTINE wfc_basis()

  USE kinds, ONLY : DP
  USE io_files, ONLY : prefix, iunwfc, nwordwfc,tmp_dir
  USE wavefunctions, ONLY : evc,psic
  USE wvfct, ONLY : nbnd,npw,npwx,et 
  USE mp, ONLY : mp_sum,mp_barrier
  USE klist, ONLY : nks,ngk,xk, igk_k
  USE becmod, ONLY : bec_type, becp, calbec,allocate_bec_type, deallocate_bec_type
  USE uspp, ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param, ONLY : upf, nh
  USE noncollin_module, ONLY: npol, noncolin
  USE mp_world, ONLY : world_comm
  USE spin_orb, ONLY: lspinorb
  USE ions_base,  ONLY : nat, nsp, ityp
  USE io_global, ONLY : stdout, ionode
  USE input_simple
  USE wannier_gw, ONLY : num_nbndv
  USE fft_base, ONLY : dffts,dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  USE gvect, ONLY : ngm, gstart,gg, g

  IMPLICIT NONE 

  COMPLEX(kind=DP), ALLOCATABLE :: omat(:,:),omatj(:),omatij(:,:)
  INTEGER ::  ipol
  TYPE(BEC_TYPE), ALLOCATABLE :: bec0(:),bec1(:),bec2(:)
  INTEGER :: ijkb0,nt,na,jh,ih,i,j,ikb,jkb
  INTEGER :: ik,ig,nfound
  INTEGER :: npw1,ntot,isca
  COMPLEX(kind=DP), ALLOCATABLE :: wfc0(:,:),wfc1(:,:),wfc3(:,:)
  COMPLEX(kind=DP) :: csca
  COMPLEX(kind=DP), EXTERNAL :: ZDOTC

  COMPLEX(kind=DP), ALLOCATABLE :: wfc_t(:,:)
  INTEGER, EXTERNAL :: find_free_unit
  INTEGER :: iun
  COMPLEX(kind=DP), ALLOCATABLE :: rwfc(:,:), cwfc(:),rprod(:),prod_g(:,:)
  real(kind=DP), ALLOCATABLE :: grids(:),gridd(:)
  INTEGER :: ii,jj
  COMPLEX(kind=DP), ALLOCATABLE :: valbec(:),zmat(:,:),valbecj(:,:,:)
  LOGICAL :: debug=.false.
  INTEGER :: ikk, jkk, kk

  COMPLEX(kind=DP), ALLOCATABLE :: wfc_e_tmp1(:),wfc_e_tmp2(:)
  TYPE(BEC_TYPE), ALLOCATABLE :: bec_e_tmp1(:),bec_e_tmp2(:)

  INTEGER :: iq(3)
  !
  CALL start_clock('optimal_basis')
  
  ! determine npw_max
  IF( nks > 1 ) THEN 
   !rewind (unit = iunigk)
   !READ( iunigk ) igk
   npw = ngk (1)
   npw_max=maxval(igk_k(1:npw,1))
   DO ik=2,nks
     !READ(iunigk) igk
     npw = ngk(ik)
     isca=maxval(igk_k(1:npw,ik))
     IF( isca > npw_max ) npw_max=isca
   enddo
  ELSE
    npw = ngk (1)
    npw_max=maxval(igk_k(1:npw,1))
  endif
  
  WRITE(stdout,*) 'NPWX NPW_MAX', npwx, npw_max

  !IF(nks>1) rewind (unit = iunigk)
  npw = ngk(1)
  
  !IF ( nks > 1 ) READ( iunigk ) igk
   
  CALL davcio(evc, 2*nwordwfc, iunwfc, 1, -1)
      

  ! allocate basis 
  ntot = num_val + num_cond
  ntot_e = num_val + num_cond
  ALLOCATE( wfc_t(npol*npwx,ntot_e) )

  !put k=1 wfcs on basis
  wfc_t=(0.d0,0.d0)
  DO ipol=0,npol-1
     wfc_t(1+ipol*npwx:npw+ipol*npwx,1:ntot_e)=evc(1+ipol*npwx:npw+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
  enddo
   
!check orthonormality
!   ALLOCATE(omat(ntot_e,ntot_e))
!   CALL ZGEMM('C','N',ntot_e,ntot_e,npol*npwx,(1.d0,0.d0),wfc_t,npol*npwx,wfc_t,npol*npwx,(0.d0,0.d0),omat,ntot_e)
!   CALL mp_sum(omat,world_comm)


!   if(ionode) THEN 
!      DO i=1,ntot_e
!         !DO j=1,ntot_e
!            WRITE(stdout,*) 'K POINT 1, ORTHONORMALITY CHECK: ', i, i, omat(i,i)
         !enddo
!      ENDDO 
!   endif
!   DEALLOCATE(omat)

  ALLOCATE(wfc0(npw_max*npol,ntot_e),wfc1(npw_max*npol,ntot_e))
  wfc0 = (0.d0, 0.d0)
  DO ipol=0,npol-1
    DO ig=1,npw
      wfc0(igk_k(ig,1)+ipol*npw_max,1:ntot_e) = &
        evc(ig+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
    ENDDO 
  ENDDO 


  !loop on k
  DO ik = 2,nks !DEBUG should start from 2
    
    !WRITE(stdout,*) 'K POINT: ', ik
    
    CALL start_clock('wfc_loop')
      
    npw1 = ngk(ik)
    npw = ngk(ik)
      
    CALL davcio(evc, 2*nwordwfc, iunwfc, ik, -1)
      

    wfc_t = (0.d0,0.d0)
    DO ipol = 0,npol-1
      wfc_t(1+ipol*npwx:npw+ipol*npwx,1:ntot) = &
        evc(1+ipol*npwx:npw+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
    ENDDO

    wfc1 = (0.d0,0.d0)
    DO ipol = 0,npol-1
      DO ig=1,npw1
        wfc1(igk_k(ig,ik)+ipol*npw_max,1:ntot) = &
          evc(ig+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
      ENDDO 
    ENDDO 


    ! we need the bec factors only at k-point ik
 
       
    ! allocate  overlap matrix
    ALLOCATE(omat(ntot_e,ntot))

    ! at this point wfc_e are expressed by k-states up to ik-1

    ! calculate overlap
    CALL start_clock('zgemm')
    CALL ZGEMM('C','N',ntot_e,ntot,npol*npw_max, (1.d0,0.d0), wfc0, &
      npol*npw_max, wfc1, npol*npw_max, (0.d0,0.d0), omat, ntot_e)
    CALL stop_clock('zgemm')
    CALL mp_sum(omat,world_comm)
        
    CALL mp_barrier(world_comm)
    IF( ionode .and. debug ) THEN 
      DO i = 1,ntot_e
        DO j = 1,ntot
          WRITE(stdout,*) 'K POINT I J, ORTHONORMALITY CHECK: ', i,j, omat(i,j)
        ENDDO 
      ENDDO 
    ENDIF 
       
    ! project out
    CALL start_clock('zgemm')
    CALL ZGEMM('N','N',npw_max*npol,ntot,ntot_e,(-1.d0,0.d0),wfc0, &
      npw_max*npol,omat,ntot_e,(1.d0,0.d0),wfc1,npw_max*npol )
    CALL stop_clock('zgemm')

    ! updates bec1
    ! loop on k-points up to ik

    !DEBUG part check if they are really orthogonal
    
    !!caluclate overlap
    IF( debug ) THEN 
      CALL start_clock('zgemm')
      CALL ZGEMM('C','N',ntot_e,ntot,npol*npw_max,(1.d0,0.d0),wfc0, &
        npol*npw_max,wfc1,npol*npw_max,(0.d0,0.d0),omat,ntot_e)
      CALL stop_clock('zgemm')
      CALL mp_sum(omat,world_comm)
          
      if(ionode.and.debug) THEN 
        DO i=1,ntot_e
          DO j=1,ntot
            WRITE(stdout,*) 'K POINT I J, ORTHONORMALITY CHECK2: ', i,j, omat(i,j)
          ENDDO 
        ENDDO 
      ENDIF 
    ENDIF  ! debug
       
    !!orhtonormalize them
    !!add to basis and updates arrays and counters
    DEALLOCATE(omat)

    CALL start_clock('wfc_optimal')   
    CALL optimal_gram_schmidt_z(ntot,wfc1,1,s_bands,nfound)
    CALL stop_clock('wfc_optimal')

    ! update wfc
    !
    ALLOCATE(wfc3(npw_max*npol,ntot_e+nfound))
    wfc3 = (0.d0,0.d0)
    wfc3(1:npw_max*npol,1:ntot_e) = wfc0(1:npw_max*npol,1:ntot_e)
    wfc3(1:npw_max*npol,ntot_e+1:ntot_e+nfound) = wfc1(1:npw_max*npol,1:nfound)
    DEALLOCATE(wfc0)
    ntot_e = ntot_e + nfound
    ALLOCATE(wfc0(npw_max*npol,ntot_e))
    wfc0 = (0.d0,0.0)
    wfc0(1:npw_max*npol,1:ntot_e) = wfc3(1:npw_max*npol,1:ntot_e)
    DEALLOCATE(wfc3)
    
    !DO the same copying for the bec factors
     
    !check orthonormality of basis
    IF( debug ) THEN 
      ALLOCATE(omat(ntot_e,ntot_e))
      CALL ZGEMM('C','N',ntot_e,ntot_e,npol*npw_max,(1.d0,0.d0), &
        wfc0,npol*npw_max,wfc0,npol*npw_max,(0.d0,0.d0),omat,ntot_e )
      CALL mp_sum(omat,world_comm)
      DO ii=1,ntot_e
        DO jj=1,ntot_e
          WRITE(stdout,*) 'ORTHONORMALITY :', ik, omat(ii,jj)
        ENDDO 
      ENDDO
      DEALLOCATE(omat)
    ENDIF 

    !recalculate bec's
    WRITE(stdout,*) 'DIMENSION OF BASIS', ntot_e
    CALL stop_clock('wfc_loop')
    !CALL print_clock('wfc_loop')
    !CALL print_clock('wfc_optimal')
    !CALL print_clock('zgemm')
  ENDDO 

  ! copy results to common variable (sigh..)
  ALLOCATE(wfc_e(npw_max*npol,ntot_e))
  wfc_e(1:npw_max*npol,1:ntot_e)=wfc0(1:npw_max*npol,1:ntot_e)
 
  !write data on file
  IF( ionode ) THEN 
    iun = find_free_unit()
    OPEN( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.wfc_basis', &
          status='unknown',form='unformatted')
    ! number of k points
    WRITE(iun) nks
    ! number of valence states to be used
    WRITE(iun) num_val
    ! number of conduction states to be used
    WRITE(iun) num_cond
    ! ntot_e
    WRITE(iun) ntot_e
  ENDIF 
  
  !loop on k points
  !if(nks>1)  rewind (unit = iunigk)
  ALLOCATE(omat(ntot_e,ntot))
  DO ik=1,nks
    
    WRITE(stdout,*) ' IK', ik
    
    !!k-point xyz
    IF( ionode ) WRITE(iun) xk(1:3,ik) !!weight
    npw = ngk(ik)

    CALL davcio(evc, 2*nwordwfc, iunwfc, ik, -1)
            
    !!energies
    IF( ionode ) WRITE(iun) et(num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond,ik) 
    
    !!matrix <wfc_e|wfc_ik>
    wfc1 = (0.d0,0.d0)
    DO ipol = 0,npol-1
      DO ig = 1,npw
        wfc1(igk_k(ig,ik)+ipol*npw_max,1:ntot) = &
            evc(ig+ipol*npwx,num_nbndv(1)-num_val+1:num_nbndv(1)+num_cond)
      ENDDO
    ENDDO
    
    CALL start_clock('zgemm')
    CALL ZGEMM('C','N',ntot_e,ntot,npol*npw_max,(1.d0,0.d0),wfc_e, &
       npol*npw_max,wfc1,npol*npw_max,(0.d0,0.d0),omat,ntot_e)
    CALL stop_clock('zgemm')
    CALL mp_sum(omat,world_comm)
     
    !check for completness
    IF( debug ) THEN 
      ALLOCATE(zmat(ntot,ntot))
      CALL start_clock('zgemm')
      CALL ZGEMM('C','N',ntot,ntot,ntot_e,(1.d0,0.d0),omat, &
        ntot_e,omat,ntot_e,(0.d0,0.d0),zmat,ntot)
      CALL stop_clock('zgemm')
      DO ii=1,ntot
        DO jj=1,ntot
          WRITE(stdout,*) 'BASIS CHECK ik ii jj', ik,ii,zmat(ii,jj)
        ENDDO 
      ENDDO
      DEALLOCATE(zmat)
    ENDIF 
    
    !write omat on disk
    IF( ionode ) THEN 
      DO i=1,ntot
        WRITE(iun) omat(1:ntot_e,i)
      ENDDO
    ENDIF
    
    !if required  calculates product terms for being checked
    IF( debug ) THEN 
      ALLOCATE(rwfc(dfftp%nnr*npol,ntot))
      ALLOCATE(cwfc(1:dffts%nnr*npol))
      ALLOCATE(grids(dffts%nnr),gridd(dfftp%nnr))
      ALLOCATE(rprod(dfftp%nnr*npol))
      ALLOCATE(prod_g(npw_max*npol,ntot))
      ALLOCATE(zmat(ntot,ntot))
      DO ii=1,ntot
        psic(1:dffts%nnr) = 0.d0
        psic(dffts%nl(1:npw_max)) = wfc1(1:npw_max,ii)
        CALL invfft('Wave', psic, dffts)
        cwfc(1:dffts%nnr) = psic(1:dffts%nnr)

        IF( npol > 1 ) THEN 
          psic(1:dffts%nnr)=0.d0
          psic(dffts%nl(1:npw_max))=wfc1(npw_max+1:npw_max+npw_max,ii)
          CALL invfft ('Wave', psic, dffts)
          cwfc(dffts%nnr+1:2*dffts%nnr)=psic(1:dffts%nnr)
        ENDIF


        DO ipol=0,npol-1
          grids(1:dffts%nnr)=dble(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
          CALL interpolate(gridd, grids, 1)
          rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)=gridd(1:dfftp%nnr)
          grids(1:dffts%nnr)=dimag(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
          CALL interpolate(gridd, grids, 1)
          rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii) = &
            rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii) + &
            (0.d0,1.d0)*gridd(1:dfftp%nnr)
        ENDDO

      ENDDO

      DEALLOCATE(cwfc,grids,gridd)
      
      DO ii=1,ntot
        DO jj=1,ntot
          !form product in real space
          rprod(1:dfftp%nnr*npol) = conjg(rwfc(1:dfftp%nnr*npol,ii))*&
                    rwfc(1:dfftp%nnr*npol,jj)
          DO ipol=0,npol-1
            psic(1:dfftp%nnr)=rprod( 1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)
            CALL fwfft ('Rho', psic, dfftp)
            DO ig=1,npw_max
              prod_g(ig+npw_max*ipol,jj)=psic(dfftp%nl(ig))
            ENDDO
          ENDDO
        ENDDO
        !calculate squared modulus of product
        CALL start_clock('zgemm')
        CALL ZGEMM('C','N', ntot, ntot, npw_max*npol, (1.d0,0.d0), prod_g, &
          npw_max*npol, prod_g, npw_max*npol, (0.d0,0.d0), zmat, ntot)
        CALL stop_clock('zgemm')
        CALL mp_sum(zmat, world_comm)
        IF( debug ) THEN 
          DO jj=1,ntot
            WRITE(stdout,*) 'DEBUG PRODUCTS ik, i ,j', ik,ii,jj,zmat(jj,jj)
          ENDDO
        ENDIF
      ENDDO

      !DO the same starting from e vectors
      ALLOCATE(wfc_e_tmp1(npw_max*npol))
      ALLOCATE(wfc_e_tmp2(npw_max*npol))
      ALLOCATE(bec_e_tmp1(nks),bec_e_tmp2(nks))
      DO ikk=1,nks
        CALL allocate_bec_type (nkb,1,bec_e_tmp1(ikk))
        CALL allocate_bec_type (nkb,1,bec_e_tmp2(ikk))
      ENDDO 

      DO ii=1,ntot
        wfc_e_tmp1=0.d0
        DO kk=1,ntot_e
          wfc_e_tmp1(1:npw_max*npol)=wfc_e_tmp1(1:npw_max*npol) + &
             wfc_e(1:npw_max*npol,kk)*omat(kk,ii)
        ENDDO
             
        DO jj=1,ntot
          wfc_e_tmp2=0.d0
          DO kk=1,ntot_e
            wfc_e_tmp2(1:npw_max*npol) = wfc_e_tmp2(1:npw_max*npol) + &
              wfc_e(1:npw_max*npol,kk)*omat(kk,jj)
          ENDDO
          !
          ALLOCATE(cwfc(1:dffts%nnr*npol))
          ALLOCATE(grids(dffts%nnr),gridd(dfftp%nnr))
          psic(1:dffts%nnr)=0.d0
          psic(dffts%nl(1:npw_max))=wfc_e_tmp1(1:npw_max)
          CALL invfft ('Wave', psic, dffts)
          cwfc(1:dffts%nnr)=psic(1:dffts%nnr)

          IF( npol > 1 ) THEN 
            psic(1:dffts%nnr) = 0.d0
            psic(dffts%nl(1:npw_max)) = wfc_e_tmp1(npw_max+1:npw_max+npw_max)
            CALL invfft('Wave', psic, dffts)
            cwfc(dffts%nnr+1:2*dffts%nnr) = psic(1:dffts%nnr)
          ENDIF
      
          DO ipol=0,npol-1
            grids(1:dffts%nnr)=dble(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
            CALL interpolate(gridd, grids, 1)
            rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii)=gridd(1:dfftp%nnr)
            grids(1:dffts%nnr)=dimag(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
            CALL interpolate(gridd, grids, 1)
            rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii) = &
              rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii) + &
              (0.d0,1.d0)*gridd(1:dfftp%nnr)
          ENDDO
          psic(1:dffts%nnr)=0.d0
          psic(dffts%nl(1:npw_max))=wfc_e_tmp2(1:npw_max)
          CALL invfft('Wave', psic, dffts)
          cwfc(1:dffts%nnr)=psic(1:dffts%nnr)

          IF( npol > 1 ) THEN 
            psic(1:dffts%nnr)=0.d0
            psic(dffts%nl(1:npw_max)) = wfc_e_tmp2(npw_max+1:npw_max+npw_max)
            CALL invfft ('Wave', psic, dffts)
            cwfc(dffts%nnr+1:2*dffts%nnr) = psic(1:dffts%nnr)
          ENDIF

          DO ipol=0,npol-1
            grids(1:dffts%nnr) = dble(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
            CALL interpolate(gridd, grids, 1)
            rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,jj)=gridd(1:dfftp%nnr)
            grids(1:dffts%nnr) = dimag(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
            CALL interpolate(gridd, grids, 1)
            rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,jj) = &
              rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,jj) + &
              (0.d0,1.d0)*gridd(1:dfftp%nnr)
          ENDDO
                
          rprod(1:dfftp%nnr*npol)=conjg(rwfc(1:dfftp%nnr*npol,ii)) * &
            rwfc(1:dfftp%nnr*npol,jj)

          DO ipol=0,npol-1
            psic(1:dfftp%nnr)=rprod( 1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)
            CALL fwfft('Rho', psic, dfftp)
            DO ig=1,npw_max
              prod_g(ig+npw_max*ipol,1)=psic(dfftp%nl(ig))
            ENDDO
          ENDDO

          csca=ZDOTC(npw_max*npol,prod_g,1,prod_g,1)
          CALL mp_sum(csca,world_comm)

          if(debug) WRITE(stdout,*) 'DEBUG PRODUCTS2 ik, i ,j', ik,ii,jj,csca
          DEALLOCATE(cwfc,grids,gridd)
        ENDDO
      ENDDO

      DEALLOCATE(wfc_e_tmp1,wfc_e_tmp2)
      DEALLOCATE(bec_e_tmp1,bec_e_tmp2)
      DEALLOCATE(rwfc,rprod,prod_g,zmat)
    
    ENDIF !debug
  
  ENDDO !ik
    
  WRITE(stdout,*) ' '
  WRITE(stdout,*) 'TOTAL NUMBER OF OPTIMAL BASIS VECTORS :', ntot_e
  WRITE(stdout,*) ' '
    
  ! check for completness of basis
  IF( debug ) THEN 
    ALLOCATE(zmat(ntot_e,ntot_e))
    CALL start_clock('zgemm')
    CALL ZGEMM('C','N',ntot_e,ntot_e,npol*npw_max, &
      (1.d0,0.d0), wfc_e, npol*npw_max, &
      wfc_e, npol*npw_max, (0.d0,0.d0), zmat, ntot_e)
    CALL stop_clock('zgemm')
    CALL mp_sum(zmat,world_comm)
   
    DO ii=1,ntot_e
      DO jj=1,ntot_e
        WRITE(stdout,*) 'CHECK OPTIMAL BASIS:', ii,jj, zmat(ii,jj)
      ENDDO 
    ENDDO 
    DEALLOCATE(zmat)
  ENDIF !debug
   
  DEALLOCATE(omat)
  if(ionode) close(iun)

  DEALLOCATE(wfc0,wfc1)

  CALL stop_clock('optimal_basis')

  WRITE(*,*)
  WRITE(*,*) 'End of wfc_basis'
  WRITE(*,*)

  RETURN 
END SUBROUTINE  wfc_basis




SUBROUTINE optimal_gram_schmidt_z(num_in,wfcs,ithres,thres,num_out)
  !this SUBROUTINE performs a gram_schmidt orthonormalization and retains
  !vectors which are above the give threshold

  USE kinds,                ONLY : DP
  USE mp_world, ONLY : world_comm, mpime, nproc
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout, ionode,ionode_id
  USE noncollin_module, ONLY: npol, noncolin
  USE input_simple, ONLY : npw_max,vkb_max
  USE becmod,        ONLY : bec_type,calbec,allocate_bec_type, deallocate_bec_type
  USE uspp,     ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param, ONLY : upf, nh
  USE noncollin_module, ONLY: npol, noncolin
  USE spin_orb, ONLY: lspinorb
  USE ions_base,  ONLY : nat, nsp, ityp

  IMPLICIT NONE 

  INTEGER, INTENT(in) :: num_in !number of initial vectors
  !
  !in input non-orthonormal in output optimal basis
  COMPLEX(kind=DP), INTENT(inout) :: wfcs(npw_max*npol,num_in)
  !
  INTEGER, INTENT(in) :: ithres ! kind of threshold
  REAL(kind=DP), INTENT(in) :: thres ! thrshold for the optimal basis
  INTEGER, INTENT(out) :: num_out ! final number of orthonormal basis functions

  INTEGER :: i,j
  COMPLEX(kind=DP), ALLOCATABLE :: prod(:)
  COMPLEX(kind=DP) :: csca
  COMPLEX(kind=DP), EXTERNAL :: zdotc
  REAL(kind=DP) :: sca

  TYPE(BEC_TYPE) :: bec0
  TYPE(BEC_TYPE), ALLOCATABLE :: bec1(:)
  INTEGER :: ijkb0,nt,na,jh,ih,ikb,jkb,ipol

  ALLOCATE(prod(num_in))
  num_out=0

  DO i=1,num_in
    IF( num_out > 0 ) THEN 
      CALL zgemv('C',npw_max*npol,num_out,(1.d0,0.d0), wfcs,npw_max*npol,wfcs(1,i),1,(0.d0,0.d0),prod,1)
      CALL mp_sum(prod(1:num_out),world_comm)
      CALL start_clock('zgemm')
      CALL zgemm('N','N',npw_max*npol,1,num_out,(-1.d0,0.d0),wfcs,npw_max*npol,prod,num_in,(1.d0,0.d0),wfcs(1,i),npw_max*npol)
      CALL stop_clock('zgemm')
    ENDIF 
    csca = zdotc(npw_max*npol,wfcs(1,i),1, wfcs(1,i),1)
    CALL mp_sum(csca,world_comm)

    IF( dble(csca) >= thres ) THEN 
      num_out = num_out + 1
      sca = dsqrt(dble(csca))
      wfcs(1:npw_max*npol,num_out) = wfcs(1:npw_max*npol,i)/sca
    ENDIF
  ENDDO 

  DEALLOCATE(prod)
  RETURN 
END SUBROUTINE optimal_gram_schmidt_z

SUBROUTINE debug_wfc(npw)
  INTEGER :: npw
  CALL errore('debug_wfc','not implemented',abs(npw))
END SUBROUTINE debug_wfc
