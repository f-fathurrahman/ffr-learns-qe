!----------------------------------------------------------------------------
SUBROUTINE my_average_pp( ntyp ) 
!----------------------------------------------------------------------------
  !! Spin-orbit pseudopotentials transformed into standard pseudopotentials.
  !
  USE kinds,            ONLY : DP
  USE atom,             ONLY : rgrid
  USE uspp_param,       ONLY : upf
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ntyp
  !! number of species
  !
  ! ... local variables
  !
  INTEGER :: nt, nb, nbe, ind, ind1, l
  REAL(DP) :: vionl
  !
  !  
  DO nt = 1, ntyp
    !
    IF( upf(nt)%has_so ) THEN
      !
      write(*,*) 'Average PP will be executed'   
      !
      IF ( upf(nt)%tvanp ) THEN
        CALL errore('my_average_pp', 'FR-PP please use lspinorb=.true.', 1 )
      ENDIF
      nbe = 0
      DO nb = 1, upf(nt)%nbeta
        nbe = nbe + 1
        IF ( upf(nt)%lll(nb) /= 0 .AND. &
             ABS( upf(nt)%jjj(nb) - upf(nt)%lll(nb) - 0.5D0 ) < 1.D-7 ) THEN
          nbe = nbe - 1
        ENDIF
      ENDDO
      !
      upf(nt)%nbeta = nbe
      nbe = 0
      DO nb = 1, upf(nt)%nbeta
        nbe = nbe + 1
        l = upf(nt)%lll(nbe)
        IF ( l /= 0 ) THEN
          IF (ABS(upf(nt)%jjj(nbe)-upf(nt)%lll(nbe)+0.5d0) < 1.d-7) THEN
            IF ( ABS( upf(nt)%jjj(nbe+1)-upf(nt)%lll(nbe+1)-0.5d0 ) > 1.d-7 ) THEN
              CALL errore( 'average_pp', 'wrong beta functions', 1 )
            ENDIF
            ind = nbe + 1
            ind1 = nbe
          ELSE
            IF (ABS(upf(nt)%jjj(nbe+1)-upf(nt)%lll(nbe+1)+0.5d0) > 1.d-7) THEN
              CALL errore( 'average_pp', 'wrong beta functions', 2 )
            ENDIF
            ind = nbe
            ind1 = nbe + 1
          ENDIF
          !
          vionl = ( ( l + 1.D0 ) * upf(nt)%dion(ind,ind) + &
                  l * upf(nt)%dion(ind1,ind1) ) / ( 2.D0 * l + 1.D0 )
          !
          upf(nt)%beta(1:rgrid(nt)%mesh,nb) = 1.D0 / ( 2.D0 * l + 1.D0 ) * &
               ( ( l + 1.D0 ) * SQRT( upf(nt)%dion(ind,ind) / vionl ) *    &
               upf(nt)%beta(1:rgrid(nt)%mesh,ind) +          &
               l * SQRT( upf(nt)%dion(ind1,ind1) / vionl ) * &
               upf(nt)%beta(1:rgrid(nt)%mesh,ind1) )
          !
          upf(nt)%dion(nb,nb) = vionl
          !
          nbe = nbe + 1
          !
        ELSE
           !
           upf(nt)%beta(1:rgrid(nt)%mesh,nb) = upf(nt)%beta(1:rgrid(nt)%mesh,nbe)
           !
           upf(nt)%dion(nb,nb) = upf(nt)%dion(nbe,nbe)
           !
        ENDIF
        !
        upf(nt)%lll(nb) = upf(nt)%lll(nbe)
        !
      ENDDO ! nb
      !
      nbe = 0
      !
      DO nb = 1, upf(nt)%nwfc
         !
         nbe = nbe + 1
         !
         IF ( upf(nt)%lchi(nb) /= 0 .AND. &
              ABS(upf(nt)%jchi(nb)-upf(nt)%lchi(nb)-0.5D0 ) < 1.D-7 ) &
            nbe = nbe - 1
         !
      ENDDO
      !
      upf(nt)%nwfc = nbe
      nbe = 0
      DO nb = 1, upf(nt)%nwfc
        nbe = nbe + 1
        l = upf(nt)%lchi(nbe)
        IF ( l /= 0 ) THEN
          IF (ABS(upf(nt)%jchi(nbe)-upf(nt)%lchi(nbe)+0.5d0) < 1.d-7) THEN
            IF ( ABS(upf(nt)%jchi(nbe+1)-upf(nt)%lchi(nbe+1)-0.5d0) > 1.d-7) then
              CALL errore( 'average_pp', 'wrong chi functions', 3 )
            ENDIF
            ind = nbe + 1
            ind1 = nbe
          ELSE
            IF( ABS(upf(nt)%jchi(nbe+1)-upf(nt)%lchi(nbe+1)+0.5d0) > 1.d-7) THEN
              CALL errore( 'average_pp', 'wrong chi functions', 4 )
            ENDIF
            ind = nbe
            ind1 = nbe+1
          ENDIF
          !
          upf(nt)%chi(1:rgrid(nt)%mesh,nb) = &
             ((l+1.D0) * upf(nt)%chi(1:rgrid(nt)%mesh,ind)+ &
               l * upf(nt)%chi(1:rgrid(nt)%mesh,ind1)) / ( 2.D0 * l + 1.D0 )
          !   
          nbe = nbe + 1
          !
        ELSE
          !
          upf(nt)%chi(1:rgrid(nt)%mesh,nb) = upf(nt)%chi(1:rgrid(nt)%mesh,nbe)
          !
        ENDIF
        !
        upf(nt)%lchi(nb) = upf(nt)%lchi(nbe)
        !
      ENDDO
      !
    ENDIF
    !
    WRITE(*,*) 'before: upf(nt)%has_so = ', upf(nt)%has_so
    upf(nt)%has_so = .FALSE.
    !
  ENDDO
  !
END SUBROUTINE my_average_pp
