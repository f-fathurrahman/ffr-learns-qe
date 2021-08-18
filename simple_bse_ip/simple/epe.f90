!SUBROUTINE for the matrix elements <e_i|p|e_j>
SUBROUTINE epe()

  USE kinds, only : DP
  USE input_simple
  USE mp, only : mp_sum,mp_barrier
  USE mp_world, only : world_comm
  USE io_global, only : ionode,stdout
  USE io_files, only : prefix, tmp_dir
  USE klist, ONLY : nks,ngk,xk
  USE becmod,        ONLY : bec_type, becp, calbec,allocate_bec_type, deallocate_bec_type
  USE uspp,     ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param, ONLY : upf, nh
  USE noncollin_module, ONLY: npol, noncolin
  USE wavefunctions, ONLY : psic
  USE fft_base,         ONLY : dffts,dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  USE gvect, ONLY : ngm, gstart,gg, g
  USE cell_base, ONLY : tpiba

  IMPLICIT NONE 
  !wfc_e(npw_max,ntot_e)
  INTEGER, EXTERNAL :: find_free_unit
  INTEGER :: i,j,a,iun
  complex(kind=DP), ALLOCATABLE :: omat(:,:,:),ge(:,:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: rwfc(:), cwfc(:),rprod(:),prod_g(:)
  real(kind=DP), ALLOCATABLE :: grids(:),gridd(:)
  INTEGER :: ipol, ii, ig, ikk, jkk, ir
  COMPLEX(kind=DP), ALLOCATABLE ::  phase(:),wfc_us_e(:,:)
  REAL(kind=DP) :: sca
  INTEGER :: iq(3),nks_prod
 
  WRITE(stdout,*) 'Enter Subroutine epe'

  ALLOCATE(omat(ntot_e,ntot_e,3))
  ALLOCATE(ge(npw_max*npol,ntot_e,3))

  !Neglect non-local terms
  !applico G ad e
  IF( .not. okvan ) THEN 
    WRITE(*,*) 'PASS HERE: .NOT. OKVAN'
    ge=0.d0
    DO a=1,3 !3 direzioni di G
       DO ipol=1,npol
          DO ig=1,npw_max
            DO i=1,ntot_e
              ge(ig+(ipol-1)*npw_max,i,a) = g(a,ig)*wfc_e(ig+(ipol-1)*npw_max,i)*2.d0*tpiba
            ENDDO  
          ENDDO  
       ENDDO 
    ENDDO  
    DO a=1,3
       CALL zgemm('C','N',ntot_e,ntot_e,npw_max*npol,(1.0d0,0.0d0),wfc_e,npw_max*npol,ge(:,:,a),&
                     &  npw_max*npol,(0.0d0,0.0d0),omat(:,:,a),ntot_e)
    ENDDO  
    CALL mp_sum(omat,world_comm)
     
  ELSE 

!we use the following approximation also for SOC
!first we calculate |\psi_i(r)|^2=\rho_i(r)
!THEN  \psi_i'(r)=sqrt(rho(r))*eps(i\phi(r)) with the same phase of the smooth part of \psi_i
!COULD BE USED EVEN FOR NC 

    ALLOCATE(rwfc(dfftp%nnr*npol))
    ALLOCATE(cwfc(1:dffts%nnr*npol))
    ALLOCATE(grids(dffts%nnr),gridd(dfftp%nnr))
    ALLOCATE(rprod(dfftp%nnr*npol))
    ALLOCATE(prod_g(npw_max*npol))
    ALLOCATE(phase(dfftp%nnr*npol))
    ALLOCATE(wfc_us_e(npw_max*npol,ntot_e))
    
    nks_prod=(2*nkpoints(1)+1)*(2*nkpoints(2)+1)*(2*nkpoints(3)+1)

    rwfc=(0.d0,0.d0)
    DO ii=1,ntot_e
      !
      psic(1:dffts%nnr)=0.d0
      psic(dffts%nl(1:npw_max))=wfc_e(1:npw_max,ii)
      CALL invfft ('Wave', psic, dffts)
      cwfc(1:dffts%nnr)=psic(1:dffts%nnr)
      !
      IF(npol > 1) THEN 
         psic(1:dffts%nnr)=0.d0
         psic(dffts%nl(1:npw_max))=wfc_e(npw_max+1:npw_max+npw_max,ii)
         CALL invfft ('Wave', psic, dffts)
         cwfc(dffts%nnr+1:2*dffts%nnr)=psic(1:dffts%nnr)
      ENDIF 
      !  
      DO ipol=0,npol-1
        grids(1:dffts%nnr)=dble(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
        !CALL interpolate (gridd, grids, 1) DEBUG ATTENZIONE DA COMPLETARE
        rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)=gridd(1:dfftp%nnr)
        grids(1:dffts%nnr)=dimag(cwfc(1+ipol*dffts%nnr:dffts%nnr+ipol*dffts%nnr))
        !CALL interpolate (gridd, grids, 1) DEBUG ATTENZIONE DA COMPLETARE     
        rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr) = &
           rwfc(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr) + &
           (0.d0,1.d0)*gridd(1:dfftp%nnr)
      ENDDO 
      !
      rprod(1:dfftp%nnr*npol) = conjg(rwfc(1:dfftp%nnr*npol)) * &
           rwfc(1:dfftp%nnr*npol)
      !
      DO ipol=0,npol-1
        psic(1:dfftp%nnr)=rprod( 1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)
        CALL fwfft ('Rho', psic, dfftp)
        DO ig=1,npw_max
          prod_g(ig+npw_max*ipol)=psic(dfftp%nl(ig))
        ENDDO 
      ENDDO 
      !
      !determine phase
      !
      DO ipol=0,npol-1
        DO ir=1,dfftp%nnr
          sca=sqrt(dble(conjg(rwfc(ir+ipol*dfftp%nnr))*rwfc(ir+ipol*dfftp%nnr)))
          IF( sca >= 1d-10 ) THEN 
            phase(ir+ipol*dfftp%nnr) = rwfc(ir+ipol*dfftp%nnr)/sca
          ELSE 
            phase(ir+ipol*dfftp%nnr) = (1.d0,0.d0)
          ENDIF 
        ENDDO 
      ENDDO 
      !
      !fft back 
      !
      DO ipol=0,npol-1
        psic(1:dfftp%nnr)=0.d0
        psic(dfftp%nl(1:npw_max))=prod_g(npw_max+1:npw_max+npw_max)
        CALL invfft ('Wave', psic, dfftp)
        rprod(dfftp%nnr*ipol+1:dfftp%nnr+dfftp%nnr*ipol)=psic(1:dfftp%nnr)  
      ENDDO 
      !
      !add phase
      DO ipol=0,npol-1
        DO ir=1,dfftp%nnr
          sca = sqrt(abs(dble(rprod(ir+ipol*dfftp%nnr))))
          rprod(ir+ipol*dfftp%nnr) = phase(ir+ipol*dfftp%nnr)*sca
        ENDDO 
      ENDDO 
      !fft forth
      DO ipol=0,npol-1
        psic(1:dfftp%nnr) = rprod(1+ipol*dfftp%nnr:dfftp%nnr+ipol*dfftp%nnr)
        CALL fwfft('Rho', psic, dfftp)
        DO ig=1,npw_max
          wfc_us_e(ig+npw_max*ipol,ii) = psic(dfftp%nl(ig))
        ENDDO 
      ENDDO 
      !
    ENDDO 

    DO a=1,3!3 direzioni di G
      DO i=1,ntot_e
        DO ipol=0,npol-1
          ge(1+ipol*npw_max:npw_max+ipol*npw_max,i,a) = &
            g(a,1:npw_max)*wfc_us_e(1+ipol*npw_max:npw_max+ipol*npw_max,i)*2.d0*tpiba
        ENDDO 
      ENDDO 
    ENDDO 
       
    DO a=1,3
      CALL zgemm('C','N',ntot_e,ntot_e,npw_max*npol,(1.0d0,0.0d0),wfc_us_e, &
          npw_max*npol,ge(1,1,a),npw_max*npol,(0.0d0,0.0d0), omat(1,1,a),ntot_e)
    ENDDO
    CALL mp_sum(omat,world_comm)

    DEALLOCATE(cwfc,grids,gridd)
    DEALLOCATE(rwfc,prod_g,rprod)
    DEALLOCATE(phase,wfc_us_e)

  ENDIF ! okvan

  IF(ionode) THEN 
    iun=find_free_unit()
    open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.epe', &
            &status='unknown',form='unformatted')
    WRITE(iun) ntot_e
    DO a=1,3
      DO i=1,ntot_e
        WRITE(iun) omat(1:ntot_e,i,a)
      ENDDO  
    ENDDO   
    close(iun)
  ENDIF  

  DEALLOCATE(ge)
  DEALLOCATE(omat)
  RETURN 
END SUBROUTINE  epe
