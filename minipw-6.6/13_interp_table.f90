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

  call prepare_all()

  cell_factor = 1.d0

  write(*,*) 'ecutwfc = ', ecutwfc

  nqx = INT( (SQRT(ecutwfc) / dq + 4) * cell_factor )
  write(*,*) 'nqx = ', nqx

  write(*,*) 'upf(:)%nbeta = ', upf(1:ntyp)%nbeta
  write(*,*) 'nbetam = ', nbetam

  ALLOCATE( tab(nqx,nbetam,ntyp) )

  !
  ! fill the interpolation table tab
  !
  ndm = MAXVAL( upf(:)%kkbeta )

  write(*,*) 'upf(:)%kkbeta = ', upf(1:ntyp)%kkbeta
  write(*,*) 'ndm = ', ndm
  
  allocate( aux(ndm) )
  allocate( besr(ndm) )

  pref = fpi/sqrt(omega)


  tab(:,:,:) = 0.d0
  do nt = 1, ntyp
    if( upf(nt)%is_gth ) cycle
    write(*,*)
    write(*,*) 'nt = ', nt
    write(*,*) 'shape rgrid%r = ', shape(rgrid(nt)%r)
    do nb = 1, upf(nt)%nbeta
      l = upf(nt)%lll(nb)
      do iq = 1,nqx
        qi = (iq - 1) * dq
        call sph_bes( upf(nt)%kkbeta, rgrid(nt)%r, qi, l, besr ) ! calculate j_l(1:nr)
        do ir = 1,upf(nt)%kkbeta
           aux(ir) = upf(nt)%beta (ir, nb) * besr(ir) * rgrid(nt)%r(ir)
        enddo
        call simpson(upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
        tab(iq, nb, nt) = vqint * pref
      enddo
    enddo
    write(*,*) 'kkbeta = ', upf(nt)%kkbeta
  enddo

  ! Write to file
  nb = 2
  nt = 1
  do iq = 1,nqx
    qi = (iq-1)*dq
    write(1000,'(1x,2F18.10)') qi, tab(iq,nb,nt) 
  enddo
  !write(*,*) 'startq = ', startq
  !write(*,*) 'lastq  = ', lastq
  !write(*,*) 'nqx    = ', nqx
  !write(*,*) 'shape tab: ', shape(tab)
  !deallocate (besr)
  !deallocate (aux)


END PROGRAM 
