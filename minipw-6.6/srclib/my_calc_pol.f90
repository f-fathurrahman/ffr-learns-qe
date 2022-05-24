! originally inner function in electrons_scf

!-----------------------------------------------------------------------
FUNCTION my_calc_pol( ) RESULT ( en_el )
!-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  use io_global, only : stdout
  use cell_base, only : at
  USE constants, ONLY : pi
  USE bp,        ONLY : ion_pol, el_pol, fc_pol, l_el_pol_old, &
                         el_pol_acc, el_pol_old, efield, l3dstring, gdir, &
                         transform_el, efield_cart
  !
  IMPLICIT NONE
  REAL (DP) :: en_el
  !
  INTEGER :: i, j 
  REAL(DP):: sca, el_pol_cart(3),  el_pol_acc_cart(3)


  IF (.NOT.l3dstring) THEN
    !
    CALL c_phase_field( el_pol(gdir), ion_pol(gdir), fc_pol(gdir), gdir )
    !
    IF (.NOT.l_el_pol_old) THEN
       l_el_pol_old = .TRUE.
       el_pol_old(gdir) = el_pol(gdir)
       en_el = -efield*(el_pol(gdir)+ion_pol(gdir))
       el_pol_acc(gdir) = 0.d0
    ELSE
       sca = (el_pol(gdir)-el_pol_old(gdir))/fc_pol(gdir)
       IF (sca < - pi) THEN
          el_pol_acc(gdir) = el_pol_acc(gdir)+2.d0*pi*fc_pol(gdir)
       ELSEIF (sca > pi) THEN
          el_pol_acc(gdir) = el_pol_acc(gdir)-2.d0*pi*fc_pol(gdir)
       ENDIF
       en_el = -efield*(el_pol(gdir)+ion_pol(gdir)+el_pol_acc(gdir))
       el_pol_old = el_pol
    ENDIF
    !
  ELSE
    !
    DO i = 1, 3
      CALL c_phase_field( el_pol(i), ion_pol(i), fc_pol(i), i )
    ENDDO
    el_pol_cart(:) = 0.d0
    DO i = 1, 3
       DO j = 1, 3
          !el_pol_cart(i)=el_pol_cart(i)+transform_el(j,i)*el_pol(j)
          el_pol_cart(i) = el_pol_cart(i)+at(i,j)*el_pol(j) / &
                           (SQRT(at(1,j)**2.d0+at(2,j)**2.d0+at(3,j)**2.d0))
       ENDDO
    ENDDO
    !
    WRITE( stdout,'( "Electronic Dipole on Cartesian axes" )' )
    DO i = 1, 3
       WRITE(stdout,*) i, el_pol_cart(i)
    ENDDO
    !
    WRITE( stdout,'( "Ionic Dipole on Cartesian axes" )' )
    DO i = 1, 3
       WRITE(stdout,*) i, ion_pol(i)
    ENDDO
    !
    IF (.NOT.l_el_pol_old) THEN
       l_el_pol_old = .TRUE.
       el_pol_old(:) = el_pol(:)
       en_el = 0.d0
       DO i = 1, 3
          en_el = en_el-efield_cart(i)*(el_pol_cart(i)+ion_pol(i))
       ENDDO
       el_pol_acc(:) = 0.d0
    ELSE
       DO i = 1, 3
          sca = (el_pol(i)-el_pol_old(i))/fc_pol(i)
          IF (sca < - pi) THEN
             el_pol_acc(i) = el_pol_acc(i)+2.d0*pi*fc_pol(i)
          ELSEIF (sca > pi) THEN
             el_pol_acc(i) = el_pol_acc(i)-2.d0*pi*fc_pol(i)
          ENDIF
       ENDDO
       el_pol_acc_cart(:) = 0.d0
       DO i = 1, 3
          DO j = 1, 3
             el_pol_acc_cart(i) = el_pol_acc_cart(i)+transform_el(j,i)*el_pol_acc(j)
          ENDDO
       ENDDO
       en_el = 0.d0
       DO i = 1, 3
          en_el = en_el-efield_cart(i)*(el_pol_cart(i)+ion_pol(i)+el_pol_acc_cart(i))
       ENDDO
       el_pol_old(:) = el_pol(:)
    ENDIF
    !
  ENDIF
  !
END FUNCTION