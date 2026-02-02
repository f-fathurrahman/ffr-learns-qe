!-----------------------------------------------------------------------
SUBROUTINE my_sym_rho_serial( ngm_, g_, nspin_, rhog_ )
!-----------------------------------------------------------------------
  !! Symmetrize the charge density rho in reciprocal space.    
  !
  USE kinds, only: DP
  USE constants, ONLY : tpi
  USE cell_base,  ONLY : at, bg
  USE symm_base,  ONLY : s, sname, ft, nsym, t_rev, invs
  use my_symme, only: ngs, shell
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngm_
  !! number of g points
  INTEGER, INTENT(IN) :: nspin_
  !! number of spin components to be symmetrized
  REAL(DP), INTENT(IN) :: g_( 3, ngm_ )
  !! list of G-vectors
  COMPLEX(DP), INTENT(INOUT) :: rhog_( ngm_, nspin_ )
  !! rho in reciprocal space: rhog_(ig) = rho(G(:,ig)). 
  !! Unsymmetrized on input, symmetrized on output
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: g0(:,:)
  REAL(DP) :: sg(3), ft_(3,48), arg
  COMPLEX(DP) :: fact, rhosum(2), mag(3), magrot(3), magsum(3)
  INTEGER :: irot(48), ig, isg, igl, ng, ns, nspin_lsda, is
  LOGICAL, ALLOCATABLE :: done(:)
  LOGICAL :: non_symmorphic(48)
  !
  ! convert fractional translations to cartesian, in a0 units
  !
  DO ns=1,nsym
    non_symmorphic(ns) = ( ft(1,ns) /= 0.0_dp .OR. &
                          ft(2,ns) /= 0.0_dp .OR. &
                          ft(3,ns) /= 0.0_dp )
    IF ( non_symmorphic(ns) ) then 
      ft_(:,ns) = at(:,1)*ft(1,ns) + at(:,2)*ft(2,ns) + at(:,3)*ft(3,ns)
    ENDIF
  ENDDO
  !
  IF( nspin_ == 4 ) THEN
    nspin_lsda = 1
  ELSEIF ( nspin_ == 1 .OR. nspin_ == 2 ) THEN
    nspin_lsda = nspin_
  ELSE
    CALL errore('sym_rho_serial','incorrect value of nspin', nspin_)
  ENDIF
  !
  ! scan shells of G-vectors
  !
  DO igl=1,ngs
    !
    ! symmetrize: \rho_sym(G) = \sum_S rho(SG) for all G-vectors in the star
    !
    ng = SIZE ( shell(igl)%vect )
    allocate ( g0(3,ng), done(ng) )
    IF ( ng < 1 ) CALL errore('sym_rho_serial','internal error',1)
    !
    !  bring G-vectors to crystal axis
    !
    DO ig=1,ng
      g0(:,ig) = g_(:,shell(igl)%vect(ig) )
    ENDDO
    CALL cryst_to_cart (ng, g0, at,-1)
    !
    !  rotate G-vectors
    !
    done(1:ng) = .false.
    DO ig=1,ng
      IF ( .NOT. done(ig)) THEN
        rhosum(:) = (0.0_dp, 0.0_dp)
        magsum(:) = (0.0_dp, 0.0_dp)
        ! S^{-1} are needed here
        DO ns=1,nsym
          sg(:) = s(:,1,invs(ns)) * g0(1,ig) + &
                  s(:,2,invs(ns)) * g0(2,ig) + &
                  s(:,3,invs(ns)) * g0(3,ig)
          irot(ns) = 0
          DO isg=1,ng
            IF ( ABS ( sg(1)-g0(1,isg) ) < 1.0D-5 .AND. &
                ABS ( sg(2)-g0(2,isg) ) < 1.0D-5 .AND. &
                ABS ( sg(3)-g0(3,isg) ) < 1.0D-5 ) THEN
              irot(ns) = isg
              EXIT
            ENDIF
          ENDDO
          IF ( irot(ns) < 1 .OR. irot(ns) > ng ) then
            CALL errore('sym_rho_serial','internal error',2)
          ENDIF
              ! isg is the index of rotated G-vector
              isg = shell(igl)%vect(irot(ns))
              !
              ! non-spin-polarized case: component 1 is the charge
              ! LSDA case: components 1,2 are spin-up and spin-down charge
              ! non colinear case: component  1 is the charge density,
              !                    components 2,3,4 are the magnetization
              ! non colinear case: components 2,3,4 are the magnetization
              !
              IF ( nspin_ == 4 ) THEN
                  ! bring magnetization to crystal axis
                  mag(:) = rhog_(isg, 2) * bg(1,:) + &
                          rhog_(isg, 3) * bg(2,:) + &
                          rhog_(isg, 4) * bg(3,:)
                  ! rotate and add magnetization
                  magrot(:) = s(1,:,invs(ns)) * mag(1) + &
                              s(2,:,invs(ns)) * mag(2) + &
                              s(3,:,invs(ns)) * mag(3)
                  IF (sname(invs(ns))(1:3)=='inv') magrot(:)=-magrot(:)
                  IF (t_rev(invs(ns)) == 1)        magrot(:)=-magrot(:)
              ENDIF
              IF ( non_symmorphic (ns) )  THEN
                  arg = tpi * ( g_(1,isg) * ft_(1,ns) + &
                                g_(2,isg) * ft_(2,ns) + &
                                g_(3,isg) * ft_(3,ns) )
                  fact = CMPLX ( COS(arg), -SIN(arg), KIND=dp )
                  DO is=1,nspin_lsda
                    rhosum(is) = rhosum(is) + rhog_(isg, is) * fact
                  ENDDO
                  IF ( nspin_ == 4 ) &
                      magsum(:) = magsum(:) + magrot(:) * fact
              ELSE
                  DO is=1,nspin_lsda
                    rhosum(is) = rhosum(is) + rhog_(isg, is)
                  ENDDO
                  IF ( nspin_ == 4 ) &
                      magsum(:) = magsum(:) + magrot(:)
              ENDIF
            ENDDO
            !
            DO is=1,nspin_lsda
              rhosum(is) = rhosum(is) / nsym
            ENDDO
            IF ( nspin_ == 4 ) magsum(:) = magsum(:) / nsym
            !
            !  now fill the shell of G-vectors with the symmetrized value
            !
            DO ns=1,nsym
              isg = shell(igl)%vect(irot(ns))
              IF ( nspin_ == 4 ) THEN
                  ! rotate magnetization
                  magrot(:) = s(1,:,ns) * magsum(1) + &
                              s(2,:,ns) * magsum(2) + &
                              s(3,:,ns) * magsum(3)
                  IF (sname(ns)(1:3)=='inv') magrot(:)=-magrot(:)
                  IF (t_rev(ns) == 1)        magrot(:)=-magrot(:)
                  ! back to cartesian coordinates
                  mag(:) = magrot(1) * at(:,1) + &
                          magrot(2) * at(:,2) + &
                          magrot(3) * at(:,3)
              ENDIF
              IF ( non_symmorphic (ns) )  THEN
                  arg = tpi * ( g_(1,isg) * ft_(1,ns) + &
                                g_(2,isg) * ft_(2,ns) + &
                                g_(3,isg) * ft_(3,ns) )
                  fact = CMPLX ( COS(arg), SIN(arg), KIND=dp )
                  DO is=1,nspin_lsda
                    rhog_(isg,is) = rhosum(is) * fact
                  ENDDO
                  IF ( nspin_ == 4 ) THEN
                    DO is=2,nspin_
                        rhog_(isg, is) = mag(is-1)*fact
                    ENDDO
                  ENDIF
              ELSE
                  DO is=1,nspin_lsda
                    rhog_(isg,is) = rhosum(is)
                  ENDDO
                  IF ( nspin_ == 4 ) THEN
                    DO is=2,nspin_
                        rhog_(isg, is) = mag(is-1)
                    ENDDO
                  ENDIF
              ENDIF
              done(irot(ns)) =.true.
            ENDDO
        ENDIF
      ENDDO
      DEALLOCATE ( done, g0 )
  ENDDO
  !
  RETURN
END SUBROUTINE my_sym_rho_serial