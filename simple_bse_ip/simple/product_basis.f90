!this routine forms a common basis for products of wave functions
!it write it to disk
!it forms product with basis for wave-functions and write it to disk

SUBROUTINE product_basis
  USE kinds, ONLY : DP
  USE io_files,             ONLY : prefix, iunwfc, nwordwfc, tmp_dir
  USE wavefunctions, ONLY : psic
  USE wvfct, ONLY : npw,npwx,et
  USE mp, ONLY : mp_sum
  USE klist, ONLY : nks,ngk,xk
  USE becmod,        ONLY : bec_type, becp, calbec,allocate_bec_type, deallocate_bec_type
  USE uspp,     ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param, ONLY : upf, nh
  USE noncollin_module, ONLY: npol, noncolin
  USE mp_world, ONLY : world_comm
  USE spin_orb, ONLY: lspinorb
  USE ions_base,  ONLY : nat, nsp, ityp
  USE io_global, ONLY : stdout, ionode
  USE input_simple
  USE gvect, ONLY : ngm, gstart,gg, g
  USE fft_base,         ONLY : dffts,dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  



  IMPLICIT NONE 

  INTEGER :: ii,jj,ipol,ig,ir,num,num_out,iun
  COMPLEX(kind=DP), ALLOCATABLE :: rwfc(:,:), cwfc(:),rprod(:),prod_g(:,:)
  real(kind=DP), ALLOCATABLE :: grids(:),gridd(:)
  COMPLEX(kind=DP), ALLOCATABLE :: valbec(:),gprod(:),sprod(:),valbecm(:,:,:,:,:,:)
  LOGICAL, PARAMETER :: ldebug=.false.
  COMPLEX(kind=DP):: csca
  COMPLEX(kind=DP), ALLOCATABLE :: omat(:,:),prod_tmp(:,:),zmat(:,:),omatij(:,:),omatm(:,:)
  INTEGER, EXTERNAL :: find_free_unit
  TYPE(BEC_TYPE) :: bec1
  INTEGER :: ikk,jkk,jjj
  INTEGER :: iq(3)
  INTEGER, PARAMETER :: block_states=1
  INTEGER :: nblock,nks_prod
  WRITE(stdout,*) 'Routine product_basis'

  !put wfc in real space without augmentation part

  ALLOCATE(rwfc(dfftp%nnr*npol,ntot_e))
  ALLOCATE(cwfc(1:dffts%nnr*npol))
  ALLOCATE(grids(dffts%nnr),gridd(dfftp%nnr))
  ALLOCATE(rprod(dfftp%nnr*npol))

  nks_prod = (2*nkpoints(1)+1)*(2*nkpoints(2)+1)*(2*nkpoints(3)+1)
  rwfc = (0.d0,0.d0)
  DO ii = 1,ntot_e
    psic(1:dffts%nnr) = 0.d0
    psic(dffts%nl(1:npw_max)) = wfc_e(1:npw_max,ii)
    CALL invfft('Wave', psic, dffts)
    cwfc(1:dffts%nnr) = psic(1:dffts%nnr)
    IF( npol > 1 ) THEN 
      psic(1:dffts%nnr) = 0.d0
      psic(dffts%nl(1:npw_max)) = wfc_e(npw_max+1:npw_max+npw_max,ii)
      CALL invfft('Wave', psic, dffts)
      cwfc(dffts%nnr+1:2*dffts%nnr) = psic(1:dffts%nnr)
    ENDIF 
    DO ipol=0,npol-1
      rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr,ii) = &
          cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr)
    ENDDO 
  ENDDO

  DEALLOCATE(cwfc,grids,gridd)

  !calculate bec's for US case
    
  !IMPORTANT BEC_E MUST BE ALREADY ALLOCATED AND PROPERLY INITIALISED

  ! loop on wfc states
  nprod_e = 0
  ALLOCATE(prod_g(npw_max*npol,ntot_e))
    
  DO ii=1,ntot_e
    CALL start_clock('Product_ciclo')
    WRITE(stdout,*) 'Forming products with ', ii
    DO jjj=1,ntot_e,block_states
      !!form products with all the remaining i>=j
      DO jj=jjj,min(jjj+block_states-1,ntot_e)
        rprod(1:dfftp%nnr*npol) = conjg(rwfc(1:dfftp%nnr*npol,ii)) * &
            rwfc(1:dfftp%nnr*npol,jj)
        DO ipol=0,npol-1
          psic(1:dfftp%nnr)=rprod( 1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)
          CALL fwfft ('Rho', psic, dfftp)
          DO ig=1,npw_max
            prod_g(ig+npw_max*ipol,jj)=psic(dfftp%nl(ig))
          ENDDO 
        ENDDO 
      ENDDO  
    ENDDO 
  
    !!project out basis already built
    !!gram schmidt with threshold

    !check for consistency
    IF(ldebug) THEN 
      ALLOCATE(cwfc(1:dffts%nnr*npol))
      ALLOCATE(sprod(dffts%nnr*npol))
      DO jj=1,ntot_e
        psic(1:dffts%nnr)=0.d0
        psic(dffts%nl(1:npw_max))=prod_g(1:npw_max,jj)
        CALL invfft ('Wave', psic, dffts)
        cwfc(1:dffts%nnr)=psic(1:dffts%nnr)
        IF(npol>1) THEN 
           psic(1:dffts%nnr)=0.d0
           psic(dffts%nl(1:npw_max))=prod_g(npw_max+1:npw_max+npw_max,jj)
           CALL invfft ('Wave', psic, dffts)
           cwfc(dffts%nnr+1:2*dffts%nnr)=psic(1:dffts%nnr)
        ENDIF
        sprod(1:dffts%nnr*npol)=cwfc(1:dffts%nnr*npol)
        csca=(0.d0,0.d0)
        DO ir=1,dffts%nnr*npol
          csca=csca+sprod(ir)
        ENDDO 
        CALL mp_sum(csca, world_comm)
        csca=csca/(dffts%nr1*dffts%nr2*dffts%nr3)
        WRITE(stdout,*) 'VERIFICA :' , ii, jj, csca
      ENDDO 
      DEALLOCATE(sprod,cwfc)
    ENDIF 
       
    !WRITE(stdout,*) 'Projecting out'
    !!project out basis already built
    num=ntot_e
    IF(nprod_e>0) THEN 
      ALLOCATE(omat(nprod_e,num))
      CALL start_clock('Product_zgemm')
      CALL ZGEMM('C','N',nprod_e,num,npw_max*npol,(1.0d0,0.d0),prod_e, &
        npw_max*npol,prod_g,npw_max*npol,(0.d0,0.d0),omat,nprod_e)
      CALL stop_clock('Product_zgemm')
      CALL mp_sum(omat,world_comm)
      CALL start_clock('Product_zgemm')
      CALL ZGEMM('N','N',npw_max*npol,num,nprod_e,(-1.d0,0.d0),prod_e, &
        npw_max*npol,omat,nprod_e,(1.d0,0.d0),prod_g,npw_max*npol)
      CALL stop_clock('Product_zgemm')
      DEALLOCATE(omat)
    ENDIF 
       
    !WRITE(stdout,*) 'Calling Gram Schmidt'

    !gram schmidt with threshold
    CALL start_clock('Product_gram')
    CALL optimal_gram_schmidt_nc(num,prod_g,s_product,num_out)
    CALL stop_clock('Product_gram')
       
    IF(num_out>0 ) THEN 
      IF(nprod_e>0) THEN  
        ALLOCATE(prod_tmp(npw_max*npol,nprod_e))
        prod_tmp(1:npw_max*npol,1:nprod_e)=prod_e(1:npw_max*npol,1:nprod_e)
        DEALLOCATE(prod_e)
        ALLOCATE(prod_e(npw_max*npol,nprod_e+num_out))
        prod_e(1:npw_max*npol,1:nprod_e)=prod_tmp(1:npw_max*npol,1:nprod_e)
        prod_e(1:npw_max*npol,nprod_e+1:nprod_e+num_out)=prod_g(1:npw_max*npol,1:num_out)
        nprod_e=nprod_e+num_out
        DEALLOCATE(prod_tmp)
      ELSE
        ALLOCATE(prod_e(npw_max*npol,nprod_e+num_out))
        prod_e(1:npw_max*npol,1:num_out)=prod_g(1:npw_max*npol,1:num_out)
        nprod_e=num_out
      ENDIF 
    ENDIF 
      
    !WRITE(stdout, *) 'PRODUCT BASIS: TOT NEW ',nprod_e,num_out
    CALL stop_clock('Product_ciclo')
  ENDDO ! on ii

  !check for orthonormality
  IF( ldebug ) THEN 
    ALLOCATE(zmat(nprod_e,nprod_e))
    CALL ZGEMM('C','N',nprod_e,nprod_e,npol*npw_max,(1.d0,0.d0),prod_e, &
      npol*npw_max,prod_e,npol*npw_max, (0.d0,0.d0),zmat,nprod_e)
    CALL mp_sum(zmat,world_comm)
    DO ii=1,nprod_e
      DO jj=1,nprod_e
        WRITE(stdout,*) 'PRODUCT BASIS ORHTONORMALITY:', ii,jj,zmat(ii,jj)
      ENDDO 
    ENDDO 
    DEALLOCATE(zmat)
  ENDIF 

  !form products of product basis with couple of wfc basis
  !consider cubic matrices for simplicity
  !and write them on disk
  IF(ionode) THEN 
    iun=find_free_unit()
    open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.product_basis', &
        status='unknown',form='unformatted')
    WRITE(iun) nprod_e !product number
    WRITE(iun) ntot_e !wfc number
  ENDIF 

  ALLOCATE(omat(nprod_e,ntot_e))

  DO ii=1,ntot_e
    DO jjj=1,ntot_e,block_states
      !!form products with all the remaining i>=j 
      DO jj=jjj,min(jjj+block_states-1,ntot_e)  
        rprod(1:dfftp%nnr*npol)=conjg(rwfc(1:dfftp%nnr*npol,ii))*rwfc(1:dfftp%nnr*npol,jj)
        DO ipol=0,npol-1
          psic(1:dfftp%nnr)=rprod( 1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)
          CALL fwfft ('Rho', psic, dfftp)
          DO ig=1,npw_max
            prod_g(ig+npw_max*ipol,jj)=psic(dfftp%nl(ig))
          ENDDO 
        ENDDO 
      ENDDO   ! jj
    ENDDO ! jjj

    IF(ldebug) THEN 
      DO jj=1,ntot_e
        csca=0.d0
        DO ig=1,npw_max*npol
           csca=csca+conjg(prod_g(ig,jj))*prod_g(ig,jj)
        ENDDO 
        CALL mp_sum(csca, world_comm)
        WRITE(stdout,*) 'DEBUG E* E PRODUCT', ii, jj, csca
      ENDDO 
    ENDIF 

    CALL ZGEMM('C','N',nprod_e,ntot_e,npw_max*npol,(1.d0,0.d0),prod_e, &
        npw_max*npol,prod_g,npw_max*npol,(0.d0,0.d0),omat,nprod_e)
    CALL mp_sum(omat,world_comm)
    
    IF(ionode) THEN 
      DO jj=1,ntot_e
        WRITE(iun) omat(1:nprod_e,jj)
      ENDDO 
    ENDIF 
  
  ENDDO 

  IF(ionode) close(iun)

  DEALLOCATE(omat)
  DEALLOCATE(rprod,prod_g)
  DEALLOCATE(rwfc)

 END SUBROUTINE  product_basis



!!To be updated in a norm-conserving way

SUBROUTINE optimal_gram_schmidt_nc(num_in,wfcs,thres,num_out)
!this SUBROUTINE performs a gram_schmidt orthonormalization and retains
!vectors which are above the give threshold

  USE kinds,                ONLY : DP
  USE mp_world, ONLY : world_comm, mpime, nproc
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout, ionode,ionode_id
  USE noncollin_module, ONLY: npol, noncolin
  USE input_simple, ONLY : npw_max,vkb_max
  USE noncollin_module, ONLY: npol, noncolin
  USE spin_orb, ONLY: lspinorb
  USE ions_base,  ONLY : nat, nsp, ityp



 IMPLICIT NONE 

  INTEGER, INTENT(in) :: num_in!number of initial vectors
  COMPLEX(kind=DP), INTENT(inout) :: wfcs(npw_max*npol,num_in)!in input non-orthonormal in output optimal basis
  REAL(kind=DP), INTENT(in) :: thres!threshold for the optimal basis          
  INTEGER, INTENT(out) :: num_out!final number of orthonormal basis functions


  INTEGER :: i,j
  COMPLEX(kind=DP), ALLOCATABLE :: prod(:)
  COMPLEX(kind=DP) :: csca
  COMPLEX(kind=DP), EXTERNAL :: zdotc
  REAL(kind=DP) :: sca


  num_out=0


  ALLOCATE(prod(num_in))
  DO i=1,num_in
     IF(num_out >0) THEN 

        CALL zgemv('C',npw_max*npol,num_out,(1.d0,0.d0), wfcs,npw_max*npol,wfcs(1,i),1,(0.d0,0.d0),prod,1)
        CALL mp_sum(prod(1:num_out),world_comm)

        CALL zgemm('N','N',npw_max*npol,1,num_out,(-1.d0,0.d0),wfcs,npw_max*npol,prod,num_in,(1.d0,0.d0),wfcs(1,i),npw_max*npol)
     ENDIF 
     csca = zdotc(npw_max*npol,wfcs(1,i),1,wfcs(1,i),1)
     CALL mp_sum(csca,world_comm)


     IF(dble(csca) >= thres) THEN 
        num_out=num_out+1
        sca=dsqrt(dble(csca))
        wfcs(1:npw_max*npol,num_out)=wfcs(1:npw_max*npol,i)/sca

     ENDIF 
  ENDDO 


  DEALLOCATE(prod)
  RETURN 
END SUBROUTINE optimal_gram_schmidt_nc



