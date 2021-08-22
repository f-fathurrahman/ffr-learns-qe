!this module describes the most important objects

!==============================================================================
MODULE simple_objects
!==============================================================================

USE kinds, ONLY : DP

  !data regarding bands and wavefunctions
  !----------------------------
  TYPE bands
  !----------------------------
    INTEGER :: numv !number of valence states (those considered for excitons only)
    INTEGER :: numc !number of conduction states
    INTEGER :: num!numv+numc
    INTEGER :: ntot_e!dimension of global to all k, basis for KS states
    INTEGER :: nk!total number of k, points
    INTEGER :: nk_loc!local number of k points
    INTEGER :: ik_first!first local k point
    INTEGER :: ik_last!last local k point
    REAL(kind=DP) :: scissor!scissor in Ry a.u.
    REAL(kind=DP) :: ene(4)!tot,diagonal,exchange,direct
    REAL(kind=DP), DIMENSION(:,:), POINTER :: k!coordinates of local k points in atomic units 
    COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: omat !overlap matrix (ntot_e,num,nk_loc)
    !KS or GW  energies of valence states (numv,nk_loc) 
    REAL(kind=DP), DIMENSION(:,:), POINTER :: en_v 
    !KS or GW  energies of valence states (numc,nk_loc)
    REAL(kind=DP), DIMENSION(:,:), POINTER :: en_c
  END TYPE bands


  ! this describes an electron-hole excitons
  !--------------------------
  TYPE exc
  !--------------------------
    INTEGER :: numv !number of valence states (those considered for excitons only)
    INTEGER :: numc !number of conduction states
    INTEGER :: num!numv+numc           
    INTEGER :: nk!total number of k, points            
    INTEGER :: nk_loc!local number of k points 
    INTEGER :: ik_first!first local k point 
    INTEGER :: ik_last!last local k point
    REAL(kind=DP) :: ene(4)!tot,diagonal,exchange,direct 
    COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER  :: avc!A_vck terms (numv,numc,nk_loc)
  END TYPE exc


  !this type describes the basis for the products of global wave-function basis vectors
  !this is not distributed with k-points parallelization
  !----------------------
  TYPE product
  !----------------------
    INTEGER :: nprod_e !number of product terms
    INTEGER :: ntot_e !dimension of global to all k, basis for KS states
    !<\mathcal{E_\alpha}|(e_{i}^*e_{j}> terms 
    COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: fij
  END TYPE product


  !this object described a potential matrix in the space of the common basis for products
  !-----------------------------
  TYPE potential
  !-----------------------------
    INTEGER :: nprod_e !number of product terms 
    COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: vpot
    INTEGER :: nkp(3) !k-points grid equally spaced
    INTEGER :: nk !total number of k, points
    INTEGER, DIMENSION(:,:,:), POINTER :: ijk ! correspondence k'-k-->i,j,k (3,k,k')
    ! V_{k'-k} (nprod_e,nprod_e,i,j,k), ijk from 1 to nkp() IMPORTANT
    COMPLEX(kind=DP), DIMENSION(:,:,:,:,:), POINTER :: vpotq 
    ! Wc_{k'-k} (nprod_e,nprod_e,i,j,k), ijk from 1 to nkp() IMPORTANT
    COMPLEX(kind=DP), DIMENSION(:,:,:,:,:), POINTER :: wpotq 
  END TYPE potential

  
!-------------------
TYPE epe
!-------------------
  INTEGER :: ntot_e
  COMPLEX(KIND=DP), DIMENSION(:,:,:), POINTER :: mat
END TYPE epe
  

INTERFACE OPERATOR(+)
   MODULE PROCEDURE sum_exc
END INTERFACE


INTERFACE OPERATOR(*)
  MODULE PROCEDURE prod_exc
  MODULE PROCEDURE prod_c_exc
END INTERFACE


INTERFACE ASSIGNMENT(=)
   MODULE PROCEDURE assign_exc 
END INTERFACE


INTERFACE OPERATOR(.hd.)
   MODULE PROCEDURE h_diag
END INTERFACE


CONTAINS


SUBROUTINE nice_write_exc(bd,simple_in,a,label)
  !write in characters an exciton, with label 
  USE input_simple_exc, ONLY : input_options
  USE mp, ONLY : mp_sum
  USE mp_world, ONLY : world_comm
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : ionode_id, ionode, stdout
  USE mp_world,  ONLY : mpime, nproc

  IMPLICIT NONE 
  TYPE(bands) :: bd
  TYPE(input_options) :: simple_in
  TYPE(exc) :: a
  INTEGER :: label

  INTEGER, EXTERNAL :: find_free_unit
  INTEGER :: iun
  CHARACTER(4) :: nfile
  INTEGER :: ik, ikl,iv,ic
  REAL(kind=DP) :: kappa(3)
  COMPLEX(kind=DP), ALLOCATABLE :: avc(:,:)

  ALLOCATE(avc(a%numv,a%numc))

  WRITE(nfile,'(4i1)') label/1000,mod(label,1000)/100,mod(label,100)/10,mod(label,10)

  IF(ionode) THEN 
    iun = find_free_unit()
    open( unit=iun, file=trim(tmp_dir)//trim(simple_in%prefix)//'.exciton.'//nfile, status='unknown')
    WRITE(iun,*) a%numv,a%numc,a%nk
    WRITE(iun,*) a%ene(1:4)
  ENDIF  
  
  DO ik= 1,a%nk
    kappa=0.d0
    DO ikl=a%ik_first,a%ik_last
      IF(ikl==ik) THEN 
        kappa(1:3)=bd%k(1:3,ik-a%ik_first+1)
      ENDIF 
    ENDDO 
    CALL mp_sum(kappa,world_comm)
    IF(ionode) WRITE(iun,*) kappa(1:3)
  ENDDO 
    
  DO ik= 1,a%nk
    avc=(0.d0,0.d0)
    DO ikl=a%ik_first,a%ik_last
      IF(ikl==ik) THEN 
        avc(1:a%numv,1:a%numc)=a%avc(1:a%numv,1:a%numc,ik-a%ik_first+1)
      ENDIF 
    ENDDO 
    CALL mp_sum(avc,world_comm)
    IF(ionode) THEN 
      DO ic=1,a%numc
        DO iv=1,a%numv
          WRITE(iun,*) avc(iv,ic)
        ENDDO 
      ENDDO 
    ENDIF 
  ENDDO 
    
  IF(ionode) close(iun)
  DEALLOCATE(avc)
  RETURN 
END SUBROUTINE nice_write_exc

SUBROUTINE initialize_epe(element)
  IMPLICIT NONE 
  TYPE(epe) :: element
  nullify(element%mat)
  RETURN 
END SUBROUTINE  initialize_epe
    
SUBROUTINE deallocate_epe(element)
  IMPLICIT NONE 
  TYPE(epe) :: element
  IF(associated(element%mat)) DEALLOCATE(element%mat)
  nullify(element%mat)
  RETURN 
END SUBROUTINE  deallocate_epe

SUBROUTINE read_epe(simple_in,element)
  USE input_simple_exc, ONLY : input_options
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : ionode_id, ionode, stdout
  USE mp_world,  ONLY : mpime, nproc
      
  IMPLICIT NONE 

  TYPE(input_options) :: simple_in
  TYPE(epe) :: element
  INTEGER, EXTERNAL :: find_free_unit
  INTEGER :: iun,i,a
  WRITE(*,*)'epe:opening file'
  IF(ionode) THEN 
    iun = find_free_unit()
    open( unit=iun, file=trim(tmp_dir)//trim(simple_in%prefix)//'.epe', status='old',form='unformatted')
    WRITE(*,*)'file opened'
    read(iun) element%ntot_e
  ENDIF  
  CALL mp_bcast(element%ntot_e,ionode_id,world_comm)
  WRITE(*,*)element%ntot_e
  ALLOCATE(element%mat(element%ntot_e,element%ntot_e,3))
  IF(ionode) THEN 
     DO a=1,3
        DO i=1,element%ntot_e
           read(iun) element%mat(1:element%ntot_e,i,a)
        ENDDO  
     ENDDO  
     close(iun)
  ENDIF  
  CALL mp_bcast(element%mat,ionode_id,world_comm)
  WRITE(*,*)'epe:read all'
  RETURN 
END SUBROUTINE  read_epe


SUBROUTINE deallocate_bands(bd)
  IMPLICIT NONE 
  TYPE(bands) :: bd
  IF(associated(bd%k)) DEALLOCATE(bd%k)
  IF(associated(bd%omat)) DEALLOCATE(bd%omat)
  IF(associated(bd%en_v)) DEALLOCATE(bd%en_v)
  IF(associated(bd%en_c)) DEALLOCATE(bd%en_c)
  nullify(bd%k)
  nullify(bd%omat)
  nullify(bd%en_v)
  nullify(bd%en_c)
  RETURN 
END SUBROUTINE  deallocate_bands

SUBROUTINE deallocate_exc(a)
  IMPLICIT NONE 
  TYPE(exc) :: a
  IF(associated(a%avc)) DEALLOCATE(a%avc)
  nullify(a%avc)
  RETURN 
END SUBROUTINE  deallocate_exc
      
SUBROUTINE deallocate_product(pd)
  IMPLICIT NONE 
  TYPE(product) :: pd
  IF(associated(pd%fij)) DEALLOCATE(pd%fij)
  nullify(pd%fij)
  RETURN 
END SUBROUTINE  deallocate_product

SUBROUTINE initialize_product(pd)
  IMPLICIT NONE 
  TYPE(product) :: pd
  nullify(pd%fij)
END SUBROUTINE  initialize_product


SUBROUTINE initialize_potential(pt)
  IMPLICIT NONE 
  TYPE(potential) :: pt
  nullify(pt%vpot)
  nullify(pt%vpotq)
  nullify(pt%ijk)
  RETURN 
END SUBROUTINE  initialize_potential


SUBROUTINE deallocate_potential(pt)
  IMPLICIT NONE 
  TYPE(potential) :: pt
  IF(associated(pt%vpot)) DEALLOCATE(pt%vpot)
  nullify(pt%vpot)
  IF(associated(pt%vpotq)) DEALLOCATE(pt%vpotq)
  nullify(pt%vpotq)
  IF(associated(pt%wpotq)) DEALLOCATE(pt%wpotq)
  nullify(pt%wpotq)
  IF(associated(pt%ijk)) DEALLOCATE(pt%ijk)
  nullify(pt%ijk)
  RETURN 
END SUBROUTINE  deallocate_potential


!this SUBROUTINE read in the bands structure from disk
!and distribute it among processors
SUBROUTINE read_bands(simple_in,bd)
  USE input_simple_exc, ONLY : input_options 
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE io_files,  ONLY : tmp_dir
  USE io_global, ONLY : ionode_id, ionode, stdout
  USE mp_world,  ONLY : mpime, nproc
  USE constants, ONLY : rytoev

  IMPLICIT NONE  
  INTEGER, EXTERNAL :: find_free_unit    
      
  TYPE(input_options) :: simple_in
  TYPE(bands) :: bd

  INTEGER :: iun,ik,i
  INTEGER :: l_blk
  REAL(kind=DP) :: xk(3)
  REAL(kind=DP), ALLOCATABLE :: et(:)
  COMPLEX(kind=DP), ALLOCATABLE :: omat(:,:)

  IF(ionode) THEN 
    iun = find_free_unit()
    open( unit=iun, file=trim(tmp_dir)//trim(simple_in%prefix)//'.wfc_basis', status='old',form='unformatted')
    read(iun) bd%nk
    read(iun) bd%numv
    read(iun) bd%numc
    read(iun) bd%ntot_e
  ENDIF  

  CALL mp_bcast(bd%nk, ionode_id, world_comm)
  CALL mp_bcast(bd%numv, ionode_id, world_comm)
  CALL mp_bcast(bd%numc, ionode_id, world_comm)
  CALL mp_bcast(bd%ntot_e, ionode_id, world_comm)

  WRITE(stdout,*) 'NUMBER OF K POINTS : ', bd%nk
  WRITE(stdout,*) 'NUMBER OF VALENCE STATES : ', bd%numv
  WRITE(stdout,*) 'NUMBER OF CONDUCTION STATES : ', bd%numc
  WRITE(stdout,*) 'NUMBER OF GLOBAL STATES : ', bd%ntot_e

  bd%num = bd%numv + bd%numc

  l_blk=bd%nk/nproc
  IF(l_blk*nproc<bd%nk) l_blk=l_blk+1

  IF(l_blk*mpime+1 <= bd%nk) THEN 
    bd%ik_first=l_blk*mpime+1 
    bd%ik_last=bd%ik_first+l_blk-1
    IF(bd%ik_last>bd%nk) bd%ik_last=bd%nk
    bd%nk_loc=bd%ik_last-bd%ik_first+1
  ELSE 
    bd%nk_loc=0
    bd%ik_first=0
    bd%ik_last=-1
  ENDIF 

!      WRITE(stdout,*) 'DEBUG nk_loc :', bd%nk_loc
!nk_loc ik_first ik_last

  IF(bd%nk_loc>0) THEN 
    ALLOCATE(bd%k(3,bd%nk_loc))
    ALLOCATE(bd%omat(bd%ntot_e,bd%num,bd%nk_loc))
    ALLOCATE(bd%en_v(bd%numv,bd%nk_loc))
    ALLOCATE(bd%en_c(bd%numc,bd%nk_loc))
  ELSE 
    nullify(bd%k)
    nullify(bd%omat)
    nullify(bd%en_v)
    nullify(bd%en_c)
  ENDIF 

  ALLOCATE(et(bd%num))
  ALLOCATE(omat(bd%ntot_e,bd%num))

  DO ik=1,bd%nk
     IF(ionode) THEN 
        read(iun) xk(1:3)
        read(iun) et(1:bd%num)
        DO i=1,bd%num
           read(iun) omat(1:bd%ntot_e,i)
        ENDDO 
     ENDIF 
     CALL mp_bcast(xk, ionode_id, world_comm )
     CALL mp_bcast(et, ionode_id, world_comm )
     CALL mp_bcast(omat, ionode_id, world_comm )
     IF(ik>=bd%ik_first .and. ik <= bd%ik_last) THEN 
        bd%k(1:3,ik-bd%ik_first+1)=xk(1:3)
        bd%en_v(1:bd%numv,ik-bd%ik_first+1)=et(1:bd%numv)
        bd%en_c(1:bd%numc,ik-bd%ik_first+1)=et(bd%numv+1:bd%num)
        bd%omat(1:bd%ntot_e,1:bd%num,ik-bd%ik_first+1)=omat(1:bd%ntot_e,1:bd%num)
     ENDIF 
  ENDDO 

  IF(ionode) close(iun)

  DEALLOCATE(et,omat)

  bd%scissor=simple_in%scissor/rytoev

  RETURN 
END SUBROUTINE  read_bands


!read product terms from disk
SUBROUTINE read_product(simple_in,pd)
  USE input_simple_exc, ONLY : input_options
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE io_files,  ONLY : tmp_dir
  USE io_global, ONLY : ionode_id, ionode, stdout
  USE mp_world,  ONLY : mpime, nproc
  IMPLICIT NONE 
  TYPE(input_options) :: simple_in
  TYPE(product) :: pd
  INTEGER, EXTERNAL :: find_free_unit
  INTEGER :: iun
  INTEGER :: ii,jj
  IF(ionode) THEN 
    iun = find_free_unit()
    open( unit=iun, file=trim(tmp_dir)//trim(simple_in%prefix)//'.product_basis', status='old',form='unformatted')
    read(iun) pd%nprod_e
    read(iun) pd%ntot_e
  ENDIF  
  CALL mp_bcast(pd%nprod_e, ionode_id, world_comm)
  CALL mp_bcast(pd%ntot_e, ionode_id, world_comm)
  WRITE(stdout,*) 'NUMBER OF PRODUCTS : ', pd%nprod_e
  WRITE(stdout,*) 'NUMBER OF GLOBAL STATES : ', pd%ntot_e
  ALLOCATE(pd%fij(pd%nprod_e,pd%ntot_e,pd%ntot_e))
  !note the order important for index corresponding to complex conjugate
  IF(ionode) THEN 
    DO ii=1,pd%ntot_e
      DO jj=1,pd%ntot_e
        read(iun) pd%fij(1:pd%nprod_e,ii,jj)
      ENDDO 
    ENDDO 
    close(iun)
  ENDIF 
  CALL mp_bcast(pd%fij, ionode_id, world_comm)
  RETURN 
END SUBROUTINE  read_product


SUBROUTINE read_potential(simple_in,pt)
  USE input_simple_exc, ONLY : input_options
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE io_files,  ONLY : tmp_dir
  USE io_global, ONLY : ionode_id, ionode, stdout
  USE mp_world,  ONLY : mpime, nproc
  IMPLICIT NONE 
  INTEGER, EXTERNAL :: find_free_unit
  TYPE(input_options) :: simple_in
  TYPE(potential) :: pt
  INTEGER :: iun,ii,ik,jk,kk
  INTEGER :: nktot
  IF(ionode) THEN 
    iun = find_free_unit()
    open( unit=iun, file=trim(tmp_dir)//trim(simple_in%prefix)//'.v_mat0', status='old',form='unformatted')
    read(iun) pt%nprod_e
  ENDIF  
  CALL mp_bcast(pt%nprod_e, ionode_id, world_comm)
  ALLOCATE(pt%vpot(pt%nprod_e,pt%nprod_e))
  IF(ionode) THEN 
    DO ii=1,pt%nprod_e
      read(iun) pt%vpot(1:pt%nprod_e,ii)
    ENDDO 
  ENDIF 
  CALL mp_bcast(pt%vpot, ionode_id, world_comm)
  IF(ionode) read(iun) pt%nkp(1:3)
  CALL mp_bcast(pt%nkp, ionode_id, world_comm)
  pt%nk=pt%nkp(1)*pt%nkp(2)*pt%nkp(3)
  ALLOCATE(pt%ijk(3,pt%nk,pt%nk),pt%vpotq(pt%nprod_e,pt%nprod_e,pt%nkp(1),pt%nkp(2),pt%nkp(3)))
  ALLOCATE(pt%wpotq(pt%nprod_e,pt%nprod_e,pt%nkp(1),pt%nkp(2),pt%nkp(3)))
  IF(ionode) THEN 
    DO ik=1,pt%nk!on k'
      DO jk=1,pt%nk! on k
        read(iun) pt%ijk(1:3,jk,ik)
      ENDDO 
    ENDDO 
  ENDIF 
  CALL mp_bcast(pt%ijk, ionode_id, world_comm)
  DO ik=1,pt%nkp(1)
    DO jk=1,pt%nkp(2)
      DO kk=1,pt%nkp(3)
        IF(ionode) THEN 
          DO ii=1,pt%nprod_e
            read(iun) pt%vpotq(1:pt%nprod_e,ii,ik,jk,kk)
          ENDDO 
        ENDIF 
        CALL mp_bcast(pt%vpotq(:,:,ik,jk,kk), ionode_id, world_comm)
      ENDDO 
    ENDDO 
  ENDDO 
  DO ik=1,pt%nkp(1)
    DO jk=1,pt%nkp(2)
      DO kk=1,pt%nkp(3)
        IF(ionode) THEN 
          DO ii=1,pt%nprod_e
            read(iun) pt%wpotq(1:pt%nprod_e,ii,ik,jk,kk)
          ENDDO 
        ENDIF 
        CALL mp_bcast(pt%wpotq(:,:,ik,jk,kk), ionode_id, world_comm)
      ENDDO 
    ENDDO 
  ENDDO 

  IF(ionode) close(iun)
  RETURN 
END SUBROUTINE read_potential


!from bands object sets up exciton a
!-------------------------
SUBROUTINE setup_exc(bd,a)
!-------------------------
  IMPLICIT NONE 
  TYPE(bands), INTENT(in) :: bd
  TYPE(exc), INTENT(out) :: a
  a%numv=bd%numv
  a%numc=bd%numc
  a%num=bd%num
  a%nk=bd%nk
  a%nk_loc=bd%nk_loc
  a%ik_first=bd%ik_first
  a%ik_last=bd%ik_last
  IF(a%nk_loc>0) THEN 
    ALLOCATE(a%avc(a%numv,a%numc,a%nk_loc))
  ELSE 
    nullify(a%avc)
  ENDIF 
  RETURN 
END SUBROUTINE  setup_exc

    FUNCTION sum_exc(a,b)
      USE io_global, ONLY : stdout
      ! a+b
      IMPLICIT NONE 
      TYPE(exc), INTENT(in) :: a,b
      TYPE(exc) :: sum_exc
      !check for consistency 
      IF(a%numv/=b%numv .or. a%numc/=b%numc.or. a%num/=b%num .or. a%nk/=b%nk &
           &.or.a%nk_loc/=b%nk_loc .or. a%ik_first/=b%ik_first .or. a%ik_last/=b%ik_last) THEN 
         WRITE(stdout,*) 'Problem with sum_exc: inconsistency'
         stop
      ENDIF 
      sum_exc%numv=b%numv
      sum_exc%numc=b%numc
      sum_exc%num=b%num
      sum_exc%nk=b%nk
      sum_exc%nk_loc=b%nk_loc
      sum_exc%ik_first=b%ik_first
      sum_exc%ik_last=b%ik_last
      IF(associated(sum_exc%avc)) DEALLOCATE(sum_exc%avc)
      IF(sum_exc%nk_loc>0) THEN 
         ALLOCATE(sum_exc%avc(sum_exc%numv,sum_exc%numc,sum_exc%nk_loc))
      ELSE
         nullify(sum_exc%avc)
      ENDIF 
      IF(a%nk_loc>0) THEN 
         sum_exc%avc(1:a%numv,1:a%numc,1:a%nk_loc)=a%avc(1:a%numv,1:a%numc,1:a%nk_loc)+&
              &b%avc(1:a%numv,1:a%numc,1:a%nk_loc)
      ENDIF 
      
      RETURN 
    END FUNCTION sum_exc

! a+b    
SUBROUTINE  sum_exc_sub(sum_exc, a,b)
  USE io_global, ONLY : stdout                                                          
  IMPLICIT NONE 
  TYPE(exc), INTENT(in) :: a,b
  TYPE(exc) :: sum_exc
  !check for consistency                                                      
  IF(a%numv/=b%numv .or. a%numc/=b%numc.or. a%num/=b%num .or. a%nk/=b%nk &
       &.or.a%nk_loc/=b%nk_loc .or. a%ik_first/=b%ik_first .or. a%ik_last/=b%ik_last) THEN 
    WRITE(stdout,*) 'Problem with sum_exc: inconsistency'
    stop
  ENDIF 

  sum_exc%numv=b%numv
  sum_exc%numc=b%numc
  sum_exc%num=b%num
  sum_exc%nk=b%nk
  sum_exc%nk_loc=b%nk_loc
  sum_exc%ik_first=b%ik_first
  sum_exc%ik_last=b%ik_last
  
  IF(associated(sum_exc%avc)) DEALLOCATE(sum_exc%avc)
  IF(sum_exc%nk_loc>0) THEN 
     ALLOCATE(sum_exc%avc(sum_exc%numv,sum_exc%numc,sum_exc%nk_loc))
  else
     nullify(sum_exc%avc)
  ENDIF 
  
  IF(a%nk_loc>0) THEN 
    sum_exc%avc(1:a%numv,1:a%numc,1:a%nk_loc) = &
        a%avc(1:a%numv,1:a%numc,1:a%nk_loc) + b%avc(1:a%numv,1:a%numc,1:a%nk_loc)
  ENDIF 

  RETURN 
END SUBROUTINE sum_exc_sub


FUNCTION prod_exc(a,b)
  USE io_global, ONLY : stdout
  USE mp,                   ONLY : mp_sum
  USE mp_world,             ONLY : world_comm
  ! scalar produc
  IMPLICIT NONE 
  TYPE(exc), INTENT(in) :: a,b
  COMPLEX(kind=DP) :: prod_exc
  COMPLEX(kind=DP), EXTERNAL :: zdotc
  !check for consistency
  IF(a%numv/=b%numv .or. a%numc/=b%numc.or. a%num/=b%num .or. a%nk/=b%nk &
       &.or.a%nk_loc/=b%nk_loc .or. a%ik_first/=b%ik_first .or. a%ik_last/=b%ik_last) THEN 
    WRITE(*,*) 'Problem with prod_exc: inconsistency'
    stop
  ENDIF 

      IF(a%nk_loc>0) THEN 
         prod_exc= zdotc(a%numv*a%numc*a%nk_loc,a%avc,1,b%avc,1)
      else
         prod_exc=(0.d0,0.d0)
      ENDIF 
      CALL mp_sum(prod_exc,world_comm)


      RETURN 
    END FUNCTION prod_exc
    
    FUNCTION prod_c_exc(z,a)
!complex * exciton
      USE io_global, ONLY : stdout
      USE mp,                   ONLY : mp_sum
      USE mp_world,             ONLY : world_comm

      IMPLICIT NONE 
      TYPE(exc), INTENT(in) :: a
      COMPLEX(kind=DP), INTENT(in) :: z
      TYPE(exc) :: prod_c_exc


      prod_c_exc%numv=a%numv
      prod_c_exc%numc=a%numc
      prod_c_exc%num=a%num
      prod_c_exc%nk=a%nk
      prod_c_exc%nk_loc=a%nk_loc
      prod_c_exc%ik_first=a%ik_first
      prod_c_exc%ik_last=a%ik_last




      !IF(associated(prod_c_exc%avc)) DEALLOCATE(prod_c_exc%avc)
      IF(prod_c_exc%nk_loc>0) THEN 
         ALLOCATE(prod_c_exc%avc(prod_c_exc%numv,prod_c_exc%numc,prod_c_exc%nk_loc))
      else
         nullify(prod_c_exc%avc)
      ENDIF 
      IF(prod_c_exc%nk_loc>0) THEN 
         prod_c_exc%avc(1:a%numv,1:a%numc,1:a%nk_loc)=z*a%avc(1:a%numv,1:a%numc,1:a%nk_loc)

      ENDIF 


      RETURN 
    END FUNCTION prod_c_exc

    SUBROUTINE assign_exc(a,b)
!x=a
      IMPLICIT NONE 

      TYPE(exc), INTENT(out) :: a
      TYPE(exc),INTENT(in) ::b


      a%numv=b%numv
      a%numc=b%numc
      a%num=b%num
      a%nk=b%nk
      a%nk_loc=b%nk_loc
      a%ik_first=b%ik_first
      a%ik_last=b%ik_last




      IF(associated(a%avc)) DEALLOCATE(a%avc)
      IF(a%nk_loc>0) THEN 
         ALLOCATE(a%avc(a%numv,a%numc,a%nk_loc))
      else
         nullify(a%avc)
      ENDIF 
      IF(a%nk_loc>0) THEN 
         a%avc(1:a%numv,1:a%numc,1:a%nk_loc)=b%avc(1:a%numv,1:a%numc,1:a%nk_loc)

      ENDIF 

      RETURN 
    END SUBROUTINE assign_exc

    SUBROUTINE randomize_exc(a)
!this SUBROUTINE set an exc object randomly

      USE random_numbers, ONLY : randy

      IMPLICIT NONE 
      TYPE(exc), INTENT(inout) :: a
      INTEGER :: i,j,k

      IF(a%nk_loc > 0) THEN 
         IF(a%nk_loc >  0 ) THEN 
            DO i=1,a%numv
               DO j=1,a%numc
                  DO k=1,a%nk_loc
                     a%avc(i,j,k)=dcmplx(randy(),randy())
                  ENDDO 
               ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
      RETURN 
    END SUBROUTINE randomize_exc

!normalize the avc vector
SUBROUTINE normalize_exc(a)
  IMPLICIT NONE 
  TYPE(exc), INTENT(inout) :: a
  COMPLEX(kind=DP) :: csca
  REAL(kind=DP) :: sca
  csca = a*a
  sca = dble(csca)
  csca = cmplx(1.d0/sqrt(sca),0.d0)
  IF(a%nk_loc > 0) a%avc=csca*a%avc
  RETURN 
END SUBROUTINE normalize_exc

!this function applies the diagonal part of the excitonic Hamiltonian
!TO DO SCISSOR OR SIMILA
FUNCTION h_diag(bd,a)
  USE io_global, ONLY : stdout,ionode
      IMPLICIT NONE 
      TYPE(exc) :: h_diag
      TYPE(exc), INTENT(in) :: a
      TYPE(bands), INTENT(in) :: bd
      INTEGER :: iv,ic,ik

      h_diag%numv=a%numv
      h_diag%numc=a%numc
      h_diag%num=a%num
      h_diag%nk=a%nk
      h_diag%nk_loc=a%nk_loc
      h_diag%ik_first=a%ik_first
      h_diag%ik_last=a%ik_last



      
      IF(h_diag%nk_loc>0) THEN 
         !IF(associated(h_diag%avc)) DEALLOCATE(h_diag%avc)
         ALLOCATE(h_diag%avc(h_diag%numv,h_diag%numc,h_diag%nk_loc))
      else
         nullify(h_diag%avc)
      ENDIF 
      IF(ionode) THEN 
 !        WRITE(stdout,*) 'DEBUG:',h_diag%nk_loc,h_diag%numc,h_diag%numv,bd%en_c(1,1)
 !        WRITE(stdout,*) 'DEBUG:',bd%nk_loc,bd%numc,bd%numv,bd%en_v(1,1)
      ENDIF 
      
      IF(h_diag%nk_loc>0) THEN 
         DO ik=1,h_diag%nk_loc
            DO ic=1,h_diag%numc
                DO iv=1,h_diag%numv
                   h_diag%avc(iv,ic,ik)=(bd%en_c(ic,ik)-bd%en_v(iv,ik)+bd%scissor)*a%avc(iv,ic,ik)  
                ENDDO 
             ENDDO 
          ENDDO 

      ENDIF 



      RETURN 
    END FUNCTION h_diag




!==============================================================================
END MODULE simple_objects
!==============================================================================