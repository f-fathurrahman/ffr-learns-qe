!===========================================================================
!
! Module mtxelmultiply_m
!
! (1) mtxelMultiply() Refactored from main.f90 By (PWD)
! Last Modified 10/21/2010 (PWD)
!
! Combine the <c,k|e^(-i(q+G).r)|v,k+q><v,k+q|e^(i(q+G`).r)|c,k>
! -------------------------------------------------
! energy denominator
! Input:
! pol%eden(band,spin) = 1/(e_val-e_cond) = energy denominators
! pol%gme(band,g-vector,spin) = plane wave matrix elements
! pol%isrtx orders the |G(i)|^2 i=1,pol%nmtx
! vwfn%isort orders |qk+g|^2 (in vwfn type)
!
! energies are apparently assumed in Rydbergs.
!
!===========================================================================

module mtxelmultiply_m
  use algos_epsilon_m
  use acc_linalg_m
  use blas_m
  use global_m
  use scalapack_m
  use timing_m, only: timing => epsilon_timing
  implicit none
  private
  public :: &
    mtxelMultiply
contains
  subroutine mtxelmultiply(scal,ntot,nrk,nst,fact,vwfn,gmetempr,gmetempc,chilocal, &
    polgme,pol,indt,pht,ipe,ispin)
    type (scalapack), intent(in) :: scal
    integer, intent(in) :: ntot
    integer, intent(in) :: nrk
    integer, intent(in) :: nst(:) !< (nrk)
    real(DP), intent(in) :: fact
    type (valence_wfns), intent(in) :: vwfn
    real(DP), dimension(:,:), intent(inout) :: gmetempr(:,:),gmetempc(:,:)
    real(DP), dimension(:,:), intent(inout) :: chilocal(:,:)
    real(DP), intent(in) :: polgme(:,:,:,:,:,:)
    type (polarizability), intent(inout) :: pol
    integer, dimension(:,:,:), intent(in) :: indt
    real(DP), dimension(:,:,:), intent(in) :: pht
    integer, intent(in) :: ipe
    integer, intent(in) :: ispin
    real(DP) :: negfact
    integer, allocatable :: tmprowindex(:),tmpcolindex(:)
    real(DP), allocatable :: tmprowph(:),tmpcolph(:)
    integer :: irk, iv, ilimit, j, it, icurr, itot, mytot
   
    itot=0
    negfact = -1D0*fact
    call timing%start(timing%chi_sum_prep)
    allocate(tmprowindex (scal%nprd(ipe)))
    allocate(tmpcolindex (scal%npcd(ipe)))
    allocate(tmprowph (scal%nprd(ipe)))
    allocate(tmpcolph (scal%npcd(ipe)))
    do irk = 1, nrk
      do it = 1, nst(irk)
        do icurr=1,scal%nprd(ipe)
          tmprowindex(icurr) = indt(scal%imyrowd(icurr,ipe),it,irk)
          tmprowph(icurr) = pht(scal%imyrowd(icurr,ipe),it,irk)
        enddo
        do icurr=1,scal%npcd(ipe)
          tmpcolindex(icurr) = indt(scal%imycold(icurr,ipe),it,irk)
          tmpcolph(icurr) = pht(scal%imycold(icurr,ipe),it,irk)
        enddo
!disabled PARALLEL DO collapse(2) private (mytot,iv,j,icurr)
        do iv = 1,peinf%nvownactual
          do j = 1, peinf%ncownactual
            mytot = itot + (iv-1)*peinf%ncownactual + j
! JRD XXX the index here probably generates gather instructions
! May want to also consider streaming stores
            do icurr=1,scal%nprd(ipe)
              gmetempr(icurr,mytot)=polgme( &
                tmprowindex(icurr),j,iv, &
                ispin,irk,pol%nfreq_group)* &
                tmprowph(icurr)
            enddo
            do icurr=1,scal%npcd(ipe)
              gmetempc(icurr,mytot)= &
                (polgme(tmpcolindex(icurr),j,iv,ispin,irk,pol%nfreq_group)*tmpcolph(icurr))
            enddo
          enddo ! j
        enddo ! iv
!disabled END PARALLEL DO
        itot = itot + peinf%nvownactual*peinf%ncownactual
      enddo ! it
    enddo ! irk
    call timing%stop(timing%chi_sum_prep)
    call timing%start(timing%chi_sum_gemm)
    if (ntot.ne.0) then
      !call dgemm('n','n',scal%nprd(ipe),scal%npcd(ipe),ntot, &
      ! negfact,gmetempr,scal%nprd(ipe),gmetempc,ntot,0.0d0,chilocal,scal%nprd(ipe))
      if (chi_summation_algo .eq. OPENACC_ALGO) then
        !$acc enter data copyin(gmetempr, gmetempc) create(chilocal)
      end if
      call acc_xgemm('n', 't', &
                     scal%nprd(ipe), scal%npcd(ipe), ntot, &
                     negfact, &
                     gmetempr, scal%nprd(ipe), &
                     gmetempc, scal%npcd(ipe), &
                     0.0d0, &
                     chilocal, scal%nprd(ipe), &
                     chi_summation_algo)
      if (chi_summation_algo .eq. OPENACC_ALGO) then
        !$acc exit data delete(gmetempr, gmetempc) copyout(chilocal)
      end if
    end if
    call timing%stop(timing%chi_sum_gemm)
    if(allocated(tmprowindex))then;deallocate(tmprowindex);endif
    if(allocated(tmpcolindex))then;deallocate(tmpcolindex);endif
    if(allocated(tmprowph))then;deallocate(tmprowph);endif
    if(allocated(tmpcolph))then;deallocate(tmpcolph);endif
   
    return
  end subroutine mtxelMultiply
end module mtxelmultiply_m
