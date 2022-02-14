INCLUDE 'prepare_all.f90'

PROGRAM main
  IMPLICIT NONE 
  CALL prepare_all()
  call my_compute_qrad()
END PROGRAM


!
! Compute interpolation table qrad(i,nm,l+1,nt) = Q^{(L)}_{nm,nt}(q_i)
! of angular momentum L, for atom of type nt, on grid q_i, where
! nm = combined index for n,m=1,nh(nt)
!
!----------------------------------------------------------------------
SUBROUTINE my_compute_qrad()
!----------------------------------------------------------------------
  USE kinds,        ONLY : dp
  USE constants,    ONLY : fpi
  USE ions_base,    ONLY : ntyp => nsp
  USE cell_base,    ONLY : omega
  USE atom,         ONLY : rgrid
  USE uspp_param,   ONLY : upf, lmaxq, nbetam
  USE us,           ONLY : nqxq, dq
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,           ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER :: ndm, startq, lastq, nt, l, nb, mb, ijv, iq, ir
  ! various indices
  REAL(DP) :: prefr
  ! the prefactor of the Q functions
  REAL(DP) :: q
  REAL(DP), ALLOCATABLE :: qrad(:,:,:,:) ! this local to this function
  !
  ! various work space
  REAL(DP), ALLOCATABLE :: aux(:), besr(:)
  !
  prefr = fpi / omega
  ndm = MAXVAL( upf(:)%kkbeta )

  write(*,*)
  write(*,*) 'nqxq (size of interpolation table) = ', nqxq
  write(*,*) 'nbetam = ', nbetam
  write(*,*) 'ndm = ', ndm
  write(*,*) 'lmaxq = ', lmaxq

  IF( lmaxq > 0 ) THEN
    ALLOCATE( qrad(nqxq, nbetam*(nbetam+1)/2, lmaxq, ntyp) )
  ELSE
    STOP 'ERROR: lmaxq <= 0'
  ENDIF
  
  ALLOCATE( aux(ndm) )
  ALLOCATE( besr(ndm) )

  CALL divide(intra_bgrp_comm, nqxq, startq, lastq)
  !
  qrad(:,:,:,:) = 0.d0
  !
  DO nt = 1, ntyp

    IF( upf(nt)%tvanp ) then

      DO l = 0, upf(nt)%nqlc -1
        ! note that l is the true (combined) angular momentum
        ! and that the arrays have dimensions 0..l (no more 1..l+1)
        DO iq = startq, lastq
          !
          q = (iq - 1) * dq
          ! here we compute the spherical bessel function for each q_i
          CALL sph_bes( upf(nt)%kkbeta, rgrid(nt)%r, q, l, besr)
          !
          DO nb = 1, upf(nt)%nbeta
            !  the Q are symmetric with respect to indices
            DO mb = nb, upf(nt)%nbeta
              ijv = mb * (mb - 1) / 2 + nb
              IF ( ( l >= abs(upf(nt)%lll(nb) - upf(nt)%lll(mb)) ) .AND. &
                   ( l <=     upf(nt)%lll(nb) + upf(nt)%lll(mb)  ) .AND. &
                   ( mod(l+upf(nt)%lll(nb) + upf(nt)%lll(mb),2)==0 ) ) THEN
                DO ir = 1, upf(nt)%kkbeta
                  aux(ir) = besr(ir) * upf(nt)%qfuncl(ir,ijv,l)
                ENDDO
                !
                ! and then we integrate with all the Q functions
                CALL simpson( upf(nt)%kkbeta, aux, rgrid(nt)%rab, qrad(iq,ijv,l+1,nt) )
              ENDIF
            ENDDO
          ENDDO
        ! igl
        ENDDO
        ! l
      ENDDO
      qrad(:,:,:,nt) = qrad(:,:,:,nt)*prefr
      CALL mp_sum( qrad(:,:,:,nt), intra_bgrp_comm )
    ENDIF 
  ENDDO ! ntyp

  ! Test loop
  DO nt = 1, ntyp
    IF( upf(nt)%tvanp ) then
      write(*,*)
      write(*,'(1x,3(A,I4))') 'Species: ', nt, ' nqlc = ', upf(nt)%nqlc, ' nbeta = ', upf(nt)%nbeta
      write(*,*) 'kkbeta  = ', upf(nt)%kkbeta
      write(*,*) 'shape qfuncl = ', shape(upf(nt)%qfuncl)
      write(*,*) 'lmax = ', upf(nt)%lmax
      ! 3rd dim of qfuncl is 0:nql-1 or 0:2*upf%lmax
      DO l = 0, upf(nt)%nqlc -1
        write(*,*) 'l = ', l
        DO nb = 1, upf(nt)%nbeta
          !  the Q are symmetric with respect to indices
          DO mb = nb, upf(nt)%nbeta
            ijv = mb * (mb - 1) / 2 + nb
            !
            write(*,'(1x,A,I4)', advance='no') 'ijv = ', ijv
            write(*,'(1x,A,2I4,A)', advance='no') '(', nb, upf(nt)%lll(nb), ') '
            write(*,'(1x,A,2I4,A)', advance='no') '(', mb, upf(nt)%lll(mb), ') '
            write(*,'(A,I4)', advance='no') 'is even number = ', l + upf(nt)%lll(nb) + upf(nt)%lll(mb)
            !
            IF ( ( l >= abs(upf(nt)%lll(nb) - upf(nt)%lll(mb)) ) .AND. &
                 ( l <=     upf(nt)%lll(nb) + upf(nt)%lll(mb)  ) .AND. &
                 (mod(l+upf(nt)%lll(nb) + upf(nt)%lll(mb),2)==0) ) THEN
              
              write(*,'(1x,A,I4,A)') '**** ijv = ', ijv, ' is accessed'
            else
              write(*,*)  ! new line
            endif
          enddo
        enddo
      enddo
    endif
  enddo



  write(*,*)
  write(*,*) 'nbetam = ', nbetam


  write(*,*) 'sum qrad = ', sum(qrad)





  !
  DEALLOCATE(besr)
  DEALLOCATE(aux)
  !
END SUBROUTINE 
