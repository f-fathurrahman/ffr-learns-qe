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
          write(*,'(I8,4F18.10)') ir, qi, rgrid(nt)%r(ir), besr(ir), upf(nt)%beta(ir,nb)
        ENDDO
        write(*,*) 'sum(aux(1:kkbeta) = ', sum(aux(1:upf(nt)%kkbeta))
        stop 'ffr 69'
        CALL simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
        tab(iq, nb, nt) = vqint * pref
      ENDDO 
    ENDDO 
    WRITE(*,*) 'kkbeta = ', upf(nt)%kkbeta
  ENDDO 

  write(*,*) 'rgrid(1)%r = ', rgrid(1)%r(1:4)

  ! Write to file
  nb = 1
  nt = 2
  DO iq = 1,nqx
    qi = (iq-1)*dq
    WRITE(1000,'(1x,2F18.10)') qi, tab(iq,nb,nt) 
  ENDDO 

  DO iq = 1,5
    qi = (iq-1)*dq
    WRITE(*,'(1x,2F18.10)') qi, tab(iq,nb,nt) 
  ENDDO 

  !write(*,*) 'startq = ', startq
  !write(*,*) 'lastq  = ', lastq
  !write(*,*) 'nqx    = ', nqx
  !write(*,*) 'shape tab: ', shape(tab)
  
  DEALLOCATE(besr)
  DEALLOCATE(aux)


END PROGRAM 
