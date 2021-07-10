include 'prepare_all.f90'

PROGRAM main
  USE atom, ONLY : rgrid
  USE uspp_param, ONLY : upf, nbetam
  USE ions_base, ONLY : ntyp => nsp !, ityp, nat
  USE cell_base, ONLY : omega
  USE gvecw, ONLY : ecutwfc
  USE constants, ONLY : fpi
  IMPLICIT NONE 
  INTEGER :: ndm
  INTEGER :: nqx ! size of interpolation table
  REAL(8), PARAMETER:: dq = 0.01d0
  !! space between points in the pseudopotential tab.
  REAL(8) :: cell_factor
  REAL(8), allocatable :: aux(:), besr(:) ! various work space
  REAL(8) :: pref
  !
  REAL(8), ALLOCATABLE :: tab(:,:,:) !interpolation table for PPs
  REAL(8) :: vqint, qi
  INTEGER :: nt, nb, l, iq, ir
  INTEGER :: i0, i1, i2, i3
  REAL(8) :: Gm, ux, vx, wx, Vq, px


  CALL prepare_all()

  cell_factor = 1.d0

  WRITE(*,*) 'ecutwfc = ', ecutwfc
  WRITE(*,*) 'dq = ', dq

  nqx = INT( (SQRT(ecutwfc) / dq + 4) * cell_factor )
  WRITE(*,*) 'nqx = ', nqx

  WRITE(*,*) 'upf(:)%nbeta = ', upf(1:ntyp)%nbeta
  WRITE(*,*) 'nbetam = ', nbetam

  ALLOCATE( tab(nqx,nbetam,ntyp) )

  !
  ! fill the interpolation table tab
  !
  ndm = MAXVAL( upf(:)%kkbeta )

  WRITE(*,*) 'upf(:)%kkbeta = ', upf(1:ntyp)%kkbeta
  WRITE(*,*) 'ndm = ', ndm

  ALLOCATE( aux(ndm) )
  ALLOCATE( besr(ndm) )

  pref = fpi/sqrt(omega)

  WRITE(*,*) 'omega = ', omega
  WRITE(*,*) 'pref = ', pref

  tab(:,:,:) = 0.d0
  DO nt = 1, ntyp
    IF( upf(nt)%is_gth ) cycle
    WRITE(*,*)
    WRITE(*,*) 'nt = ', nt
    WRITE(*,*) 'shape rgrid%r = ', shape(rgrid(nt)%r)
    DO nb = 1, upf(nt)%nbeta
      l = upf(nt)%lll(nb)
      DO iq = 1,nqx
        qi = (iq - 1) * dq
        call sph_bes( upf(nt)%kkbeta, rgrid(nt)%r, qi, l, besr ) ! calculate j_l(1:kkbeta)
        DO ir = 1,upf(nt)%kkbeta
          aux(ir) = upf(nt)%beta(ir, nb) * besr(ir) * rgrid(nt)%r(ir)
        ENDDO
        CALL simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
        tab(iq, nb, nt) = vqint * pref
      ENDDO 
    ENDDO 
    WRITE(*,*) 'kkbeta = ', upf(nt)%kkbeta
  ENDDO 

  write(*,*) 'rgrid(1)%r = ', rgrid(1)%r(1:4)

  ! For testing
  nb = 1
  nt = 2

  ! Write to file
  DO iq = 1,nqx
    qi = (iq-1)*dq
    WRITE(1000,'(1x,2F18.10)') qi, tab(iq,nb,nt) 
  ENDDO 

  DO iq = 1,5
    qi = (iq-1)*dq
    WRITE(*,'(1x,2F18.10)') qi, tab(iq,nb,nt) 
  ENDDO 
  
  DEALLOCATE(besr)
  DEALLOCATE(aux)

  Gm = 2.0d0 ! already multiplied by tpiba in QE
  ! Interpolation procedure
  px = Gm/dq - INT(Gm/dq)
  write(*,*) 'px = ', px
  ux = 1.d0 - px
  vx = 2.d0 - px
  wx = 3.d0 - px
  i0 = int(Gm/dq) + 1
  write(*,*) 'i0 = ', i0
  i1 = i0 + 1
  i2 = i0 + 2
  i3 = i0 + 3
  Vq = tab(i0,nb,nt) * ux * vx * wx / 6.d0 + &
       tab(i1,nb,nt) * px * vx * wx / 2.d0 - &
       tab(i2,nb,nt) * px * ux * wx / 2.d0 + &
       tab(i3,nb,nt) * px * ux * vx / 6.d0

  WRITE(*,'(1x,A,F18.10)') 'Vq = ', Vq

END PROGRAM 
