! eFeFeR, March 2014

PROGRAM get_evc
  USE kinds, ONLY : DP
  !
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp_global, ONLY : mp_startup, mp_global_end
  USE environment, ONLY : environment_start
  !
  USE io_files,  ONLY : tmp_dir, prefix, iunwfc, nwordwfc
  USE mp, ONLY : mp_bcast
  USE mp_world, ONLY : world_comm
  !
  USE cell_base, ONLY : omega, tpiba2, at, alat
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau, atm
  USE lsda_mod, ONLY : lsda, isk, nspin, current_spin
  USE klist, ONLY : nks, xk,  nkstot
  USE gvect, ONLY : nl, ngm, g
  USE gvecs, ONLY : nls
  USE gvecw, ONLY : ecutwfc
  USE wvfct, ONLY: igk, g2kin, npw, npwx, nbnd
  USE wavefunctions_module, ONLY : evc, psic
  USE uspp, ONLY : vkb
  !
  USE fft_base, ONLY : dfftp, dffts
  USE fft_interfaces, ONLY : fwfft, invfft

  IMPLICIT NONE
  
  ! Local
  INTEGER :: ik, ig
  INTEGER :: ibnd
  INTEGER :: ounit, ios
  !
  INTEGER :: ndata, idata
  INTEGER, ALLOCATABLE :: k_evc(:)
  INTEGER, ALLOCATABLE :: band_evc(:)
  CHARACTER(256) :: xsfname, str1, str2

  !
  ! Initialize environment
  !
  CALL mp_startup()  ! I will assume that MPI is used as default
  !
  CALL environment_start('get_evc')
  !
  IF(ionode) CALL input_from_file()


  ! Open
  OPEN(UNIT=52,FILE='EVC_LIST',ACTION='READ',STATUS='OLD',FORM='FORMATTED',IOSTAT=IOS)
  IF(IOS /= 0) THEN
    WRITE(stdout,*) 'Error opening file EVC_LIST'
    STOP
  ENDIF
  READ(52,*) ndata
  ALLOCATE( k_evc(ndata) )
  ALLOCATE( band_evc(ndata) )
  DO idata=1,ndata
    READ(52,*) k_evc(idata), band_evc(idata)
  ENDDO
  CLOSE(52)

  WRITE(stdout,*) 'Read ndata = ', ndata
  WRITE(stdout,*) 'k_evc = ', k_evc
  WRITE(stdout,*) 'band_evc = ', band_evc

  ! Hard-coded FIXME
  prefix = 'pwscf'
  tmp_dir = '../tmp/'  ! the last '/' is important
  !
  ! Probably not needed because we only use one processor
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )

  !
  ! Here we go ....
  !
  CALL read_file()

  CALL openfil_pp()   ! some subroutines need this

  !
  ! Main loop
  !
  DO idata=1,ndata
    ik = k_evc(idata)
    ibnd = band_evc(idata)

    WRITE(stdout,*) 'ik = ', ik
    WRITE(stdout,*) 'ibnd = ', ibnd

    IF(lsda) current_spin = isk(ik)
   
    WRITE(stdout,*) 'Calling gk_sort'
    CALL gk_sort (xk(1, ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    ! read eigenvectors
    WRITE(stdout,*) 'Calling davcio'
    CALL davcio(evc, nwordwfc, iunwfc, ik, -1)

    WRITE(stdout,*) 'Calling init_us_2'
    CALL init_us_2 (npw, igk, xk(1,ik), vkb)
    !CALL callbec(npw, vkb, evc, becp)  ! need this?

    psic(1:dffts%nnr) = (0.d0,0.d0)
    DO ig=1,npw
      psic(nls(igk(ig))) = evc(ig,ibnd)
    ENDDO
    CALL invfft('Wave', psic, dffts)  ! use capital 'W' in Wave

    ounit = 33   ! somewhat arbitarily chosen
    IF(ik<10) THEN
      WRITE(str1, '(I1)') ik
    ELSE
      WRITE(str1, '(I2)') ik
    ENDIF
    IF(ibnd<10) THEN
      WRITE(str2, '(I1)') ibnd
    ELSEIF(ibnd<100) THEN
      WRITE(str2, '(I2)') ibnd
    ELSEIF(ibnd<1000) THEN
      WRITE(str2, '(I3)') ibnd
    ENDIF
    xsfname = 'PSIC_'//trim(str1)//'_'//trim(str2)//'.xsf'
    !WRITE(stdout,*) trim(xsfname)
    OPEN(unit=ounit, file=xsfname, action='WRITE')
    CALL xsf_struct (alat, at, nat, tau, atm, ityp, ounit)
    CALL xsf_fast_datagrid_3d( real(psic), dffts%nr1, dffts%nr2, dffts%nr3, &
       dffts%nr1x, dffts%nr2x, dffts%nr3x, at, alat, ounit)
    CLOSE(unit=ounit)

  ENDDO  ! end loop over ndata

  ! Free memory
  DEALLOCATE( k_evc )
  DEALLOCATE( band_evc )

  !
  ! Stop the program
  !
  !CALL stop_pp()
  CALL mp_global_end()


END PROGRAM


