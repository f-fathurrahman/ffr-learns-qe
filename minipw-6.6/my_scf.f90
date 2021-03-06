PROGRAM main

  IMPLICIT NONE 

  CALL prepare_all()

  CALL my_scf()

END PROGRAM 


SUBROUTINE my_scf()

  USE kinds, ONLY : DP

  USE cell_base, ONLY : alat, at, bg, omega

  USE ions_base, ONLY : zv, nat, nsp, ityp, tau

  USE fft_base, ONLY : dfftp

  USE gvect, ONLY : ngm, g, gg, gstart, gcutm
  USE gvecw, ONLY : ecutwfc

  USE vlocal, ONLY : strf

  USE klist, ONLY : nelec

  USE ener, ONLY : ewld, deband, hwf_energy, eband, ehart, demet, &
                   etxc, etxcc, vtxc, etot

  USE scf, ONLY : scf_type, create_scf_type, destroy_scf_type, &
                  scf_type_COPY, &
                  rho, v, vnew, &
                  vltot, vrs, &
                  rho_core, rhog_core, &
                  open_mix_file, close_mix_file

  USE lsda_mod, ONLY : nspin

  USE control_flags, ONLY : gamma_only, niter, mixing_beta, ethr, tr2, &
                            conv_elec, nmix, scf_must_converge
  
  USE io_files, ONLY : iunmix

  USE ldaU, ONLY : eth

  USE extfield, ONLY : etotefield

  USE io_global, ONLY : stdout

  IMPLICIT NONE 

  ! External function
  REAL(DP) :: ewald

  TYPE(scf_type) :: rhoin

  INTEGER :: idum, iter, i
  REAL(DP) :: dr2
  LOGICAL :: first, exst
  INTEGER :: Npoints
  REAL(DP) :: deband_hwf
  REAL(DP) :: tr2_min, tr2_final
  REAL(DP) :: charge
  REAL(DP) :: descf


  WRITE(*,*)
  WRITE(*,*) 'my_scf is starting'
  WRITE(*,*)

  Npoints = dfftp%nr1 * dfftp%nr2 * dfftp%nr3

  ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, gamma_only, strf )

  WRITE(*,*) 'Ewald energy (in Ha) = ', 0.5d0*ewld


  WRITE(*,*)
  WRITE(*,'(1x,A,F18.10)') 'integ rhoe = ', sum(rho%of_r(:,:))*omega/Npoints
  WRITE(*,*)

  CALL create_scf_type( rhoin )


  CALL open_mix_file( iunmix, 'mix', exst )

  WRITE(*,*) 'exst = ', exst


  dr2 = 0.d0
  iter = 0

  tr2_final = tr2
  WRITE(*,*) 'tr2 = ', tr2

  DO idum = 1, niter
    
    iter = iter + 1

    WRITE(*,'(1x,A,I4,A,F9.2,A,F5.2)') &
         'iter = ', iter, ' ecut = ', ecutwfc, ' betamix = ', mixing_beta


    IF( iter > 1 ) THEN 
      IF( iter == 2 ) ethr = 1.0d-2
      ethr = min( ethr, 0.1d0*dr2/max(1.d0, nelec) )
      ethr = max( ethr, 1.0d-13 )
    ENDIF 

    first = ( iter == 1 )

    deband_hwf = my_delta_e()

    WRITE(*,*) 'deband_hwf = ', deband_hwf

    CALL scf_type_COPY( rho, rhoin )

    scf_step: DO 
      
      tr2_min = 0.d0

      IF( first ) tr2_min = ethr*max(1.d0, nelec)

      CALL c_bands( iter )

      !CALL poolrecover( et, nbnd, nkstot, nks )

      CALL sum_band()

      hwf_energy = eband + deband_hwf + (etxc - etxcc) + ewld + ehart + demet

      deband = my_delta_e()

      CALL mix_rho( rho, rhoin, mixing_beta, dr2, tr2_min, iter, nmix, iunmix, conv_elec )
      !CALL bcast_scf_type( rhoin, root_pool, inter_pool_comm )
      !CALL mp_bcast( dr2, root_pool, inter_pool_comm )
      !CALL mp_bcast( conv_elec, root_pool, inter_pool_comm )

      IF( .NOT. scf_must_converge .AND. idum == niter ) conv_elec = .true.


      ! rediagonalization is skipped


    
      IF( .NOT. conv_elec ) THEN 
        
        CALL v_of_rho( rhoin, rho_core, rhog_core, &
                       ehart, etxc, vtxc, eth, etotefield, charge, v )
        WRITE(*,*) 'ehart = ', ehart/2

        !
        descf = my_delta_escf()

        ! rhoin <- rho
        CALL scf_type_COPY( rhoin, rho )

    
      ELSE ! convergence is reached

        vnew%of_r(:,:) = v%of_r(:,:)
        CALL v_of_rho( rho, rho_core, rhog_core, &
                       ehart, etxc, vtxc, eth, etotefield, charge, v )

        vnew%of_r(:,:) = v%of_r(:,:) - vnew%of_r(:,:)

        descf = 0.d0

      ENDIF 

      EXIT scf_step
    ENDDO scf_step


    CALL sum_vrs( dfftp%nnr, nspin, vltot, v%of_r, vrs )

    etot = eband + (etxc - etxcc) + ewld + ehart + deband + demet + descf

    WRITE(*,*)
    WRITE(*,*) 'Total energy components in Ha'
    WRITE(*,*)
    WRITE(*,'(1x,A,F18.10)') 'eband        = ', eband/2
    WRITE(*,'(1x,A,F18.10)') 'etxc         = ', etxc/2
    WRITE(*,'(1x,A,F18.10)') 'etxcc        = ', etxcc/2
    WRITE(*,'(1x,A,F18.10)') 'ehart        = ', ehart/2
    WRITE(*,'(1x,A,F18.10)') 'deband       = ', deband/2
    WRITE(*,'(1x,A,F18.10)') 'demet        = ', demet/2
    WRITE(*,'(1x,A,F18.10)') 'descf        = ', descf/2
    WRITE(*,*) '---------------------------------'
    WRITE(*,'(1x,A,F18.10)') 'Etot in (Ha) = ', etot/2
    WRITE(*,*) '---------------------------------'



    IF( conv_elec ) THEN 
      WRITE(*,*)
      WRITE(*,*) 'Convergence is achieved'
      GOTO 10
    ENDIF 


  ENDDO 


10 flush(stdout)

  WRITE(*,*)
  WRITE(*,*) 'my_scf is finished'
  WRITE(*,*)



CONTAINS 

!====================
FUNCTION my_delta_e()
!====================
  IMPLICIT NONE 
  REAL(DP) :: my_delta_e

  my_delta_e = 0.d0

  ! nspin == 1
  my_delta_e = -sum( rho%of_r(:,:) * v%of_r(:,:) )
  my_delta_e = omega*my_delta_e / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)

END FUNCTION 


FUNCTION my_delta_escf()
  IMPLICIT NONE 
  REAL(DP) :: my_delta_escf
  ! local
  !REAL(8) :: rho_dif(2)

  my_delta_escf = 0.d0

  ! must have nspin==1

  my_delta_escf = -sum( (rhoin%of_r(:,:) - rho%of_r(:,:))*v%of_r(:,:) )

  my_delta_escf = my_delta_escf*omega/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )

  !call mp_sum( delta_escf, intra_bgrp_comm )


END FUNCTION 




END SUBROUTINE 

