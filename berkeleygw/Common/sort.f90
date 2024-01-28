!>==========================================================================
!!
!! Module sort_m:
!!
!! (1) gcutoff
!!
!! Given G-vectors sorted by kinetic energy and an energy cutoff,
!! find the corresponding G-vector cutoff.
!!
!! (2,3) sortrx, sortix
!!
!! Sorts an array by the quicksort method. real(DP) and integer versions.
!! See included sort_inc.f90.
!!
!! (4) make_identity_symmetry_first
!!
!! The identity must always be the first symmetry, as assumed in various places
!! in the code. We enforce this by swapping it with op #1 if it is not first.
!!
!! (5) sort_symmetries
!!
!! Bring symmetries into a standardized order. We are not currently using this.
!!
!!==========================================================================

module sort_m
  use global_m
  implicit none
  private
  public :: &
    gcutoff, &
    sortrx, &
    sortix, &
    make_identity_symmetry_first, &
    sort_symmetries, &
    sort_sequential, &
    sort_threaded, &
    get_GK_array_from_gvecs
  interface sortrx
    module procedure sortrx_gvec, sortrx_no_gvec
  end interface
  interface sortix
    module procedure sortix_gvec, sortix_no_gvec
  end interface
  interface sort_sequential
    module procedure &
      sortrx_no_gvec, sortrx_gvec_GK, sortrx_gvec, &
      sortix_no_gvec, sortix_gvec_GK, sortix_gvec
  end interface sort_sequential
  interface sort_threaded
    module procedure &
      sortrx_no_gvec_threaded, &
      sortrx_gvec_threaded_GK, &
      sortrx_gvec_threaded, &
      sortix_no_gvec_threaded, &
      sortix_gvec_threaded_GK, &
      sortix_gvec_threaded
  end interface sort_threaded
contains
  !> Given G-vectors sorted by kinetic energy and an energy cutoff, find the corresponding G-vector cutoff
  !! such that all(ekin(isrtrq(ig)) <= ecutoff) for ig <= gcutoff.
  integer function gcutoff(ng, ekin, isrtrq, ecutoff)
    integer, intent(in) :: ng !< number of G-vectors
    real(DP), intent(in) :: ekin(:) !< (ng) kinetic energies, should be sorted already
    integer, intent(in) :: isrtrq(:) !< (ng) this is the index array returned by sorting ekin
    real(DP), intent(in) :: ecutoff !< energy cutoff, in same units as ekin (Ry typically)
    integer :: gup, gdn, gmid, ig
   
    ! perhaps all G-vectors fall within the cutoff
    if(ekin(isrtrq(ng)) < ecutoff) then
      gcutoff = ng
     
      return
    endif
    ! otherwise, use bisection
    gup = ng
    gdn = 1
    do ig = 1, ng
      gmid = (gup + gdn) / 2
      if(gmid == gdn) exit
      if(ekin(isrtrq(gmid)) > ecutoff) then
        gup = gmid
      else
        gdn = gmid
      endif
    enddo
    gcutoff = gdn
   
    return
  end function gcutoff
!=====================================================================
!> The identity must always be the first symmetry, as assumed in various places in the code.
  subroutine make_identity_symmetry_first(nsyms, mtrx, tnp)
    integer, intent(in) :: nsyms
    integer, intent(inout) :: mtrx(3, 3, 48)
    real(DP), intent(inout) :: tnp(3, 48)
    integer :: isym, mtrx_temp(3, 3), identity(3,3)
    real(DP) :: tnp_temp(3)
    logical :: found
   
    identity = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), shape(identity))
    found = all(mtrx(1:3, 1:3, 1) == identity(1:3, 1:3))
    do isym = 2, nsyms
      if(all(mtrx(1:3, 1:3, isym) == identity(1:3, 1:3))) then
        if(.not. found) then
          ! if identity is not first, swap
          mtrx_temp(1:3, 1:3) = mtrx(1:3, 1:3, 1)
          mtrx(1:3, 1:3, 1) = mtrx(1:3, 1:3, isym)
          mtrx(1:3, 1:3, isym) = mtrx_temp(1:3, 1:3)
          tnp_temp(1:3) = tnp(1:3, 1)
          tnp(1:3, 1) = tnp(1:3, isym)
          tnp(1:3, isym) = tnp_temp(1:3)
          found = .true.
          write(0,'(a,i2)') 'WARNING: making identity first by swapping with symmetry op #', isym
        else
          call die("There is a duplicate identity in the symmetry operations.")
        endif
      endif
    enddo
    if(.not. found) then
      call die("Identity is not present in the list of symmetries.")
    endif
   
    return
  end subroutine make_identity_symmetry_first
!=====================================================================
!> Bring symmetries into a standardized order.
!! The identity is always the first one.
  subroutine sort_symmetries(nsyms, mtrx, tnp)
    integer, intent(in) :: nsyms
    integer, intent(inout) :: mtrx(3, 3, 48)
    real(DP), intent(inout) :: tnp(3, 48)
    integer :: isym, ii, jj, factor, hash(48), order(48), mtrx_temp(3, 3, 48), identity(3,3)
    real(DP) :: tnp_temp(3, 48)
   
    identity = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), shape(identity))
    do isym = 1, nsyms
      ! make sure the identity comes first
      if(all(mtrx(1:3, 1:3, isym) == identity(1:3, 1:3))) then
        hash(isym) = -1d9
        cycle
      endif
      hash(isym) = 0
      factor = 1
      do jj = 1, 3
        if(jj > 1) factor = factor * 3
        do ii = 1, 3
          if(ii > 1) factor = factor * 3
          hash(isym) = hash(isym) + mtrx(4 - ii, 4 - jj, isym) * factor
        enddo
      enddo
    enddo
    call sortix(nsyms, hash, order)
    do isym = 1, nsyms
      mtrx_temp(1:3, 1:3, isym) = mtrx(1:3, 1:3, order(isym))
      tnp_temp(1:3, isym) = tnp(1:3, order(isym))
    enddo
    mtrx(1:3, 1:3, 1:nsyms) = mtrx_temp(1:3, 1:3, 1:nsyms)
    tnp(1:3, 1:nsyms) = tnp_temp(1:3, 1:nsyms)
   
    return
  end subroutine sort_symmetries
!=====================================================================
!> Internal subroutine to get the ranking array GK from the G-vectors
!! array gvec. This is used internally by the sorting routines, and should
!! be transparent to most developers.
subroutine get_GK_array_from_gvecs(NVAL, gvec, GK)
  integer, intent(in) :: NVAL
  integer, intent(in) :: gvec(3,NVAL)
  integer, intent(out) :: GK(NVAL)
  integer :: eff_grid(3), ii
 
  eff_grid(1:3) = maxval(gvec(1:3,1:NVAL), 2) - minval(gvec(1:3,1:NVAL), 2) + 1
  !!disabled PARALLEL DO DEFAULT(SHARED) PRIVATE(ii)
  do ii = 1, NVAL
    GK(ii) = gvec(3,ii) + eff_grid(3)*(gvec(2,ii) + eff_grid(2)*gvec(1,ii))
  enddo
  !!disabled END PARALLEL DO
 
end subroutine get_GK_array_from_gvecs
! FHJ: Use the preprocessor to create the following routines:
! sortix_gvec, sortix_no_gvec, sortrx_gvec, sortrx_no_gvec
! This file is based on the work of Michel Olagnon.
! The original code for the MRGRNK subroutine is available at:
! http://fortran-2000.com/rank/mrgrnk.f90
! MRGRNK - Copyright (c) Michel Olagnon
! Copying and distribution of this file, with or without modification,
! are permitted in any medium without royalty provided the copyright
! notice and this notice are preserved. This file is offered as-is,
! without any warranty.
! FHJ: WARNING - make sure you don`t change the following lines too much,
! otherwise they will be longer than 120 characters after the preprocessors kicks in.
! Note that, if there the extra "gvec" argument, we use a tolerance to figure
! out if the two items AA(ii) and AA(jj) are degenerate.
! FHJ: Note - we need to nest the "JOIN" macro otherwise the symbol sortix_gvec
! doesn`t get expanded by the C preprocessor.
! LABEL_GK is the sorting function without automatically generating the GK
! array out of the gvec array. This is useful for threaded sorting, when we
! first generate the GK array, but then call the sorting routine many times.
! LABEL_GK is the same as sortix_gvec if there is no gvec input.
! This is the kernel of the sourting routine
! These functions are threaded wrappers for the sorting routines
!> Sorts (actually, ranks) a small array AA using the insert sort algorithm.
subroutine sortix_gvec_insertsort(NVAL, AA, ord&
, GK&
)
  integer, intent(in) :: NVAL
  integer, intent(in) :: AA(NVAL)
  integer, intent(inout) :: ord(NVAL)
  integer, intent(in) :: GK(NVAL)
  integer, parameter :: TOL=TOL_ZERO
  integer :: ii, jj, tord
 
  do ii = 2, NVAL
    tord = ord(ii)
    jj = ii - 1
    do while (jj>0)
      if (.not.(AA(ord(jj))-AA(tord)>TOL.or.(AA(ord(jj))-AA(tord)>-TOL.and.GK(ord(jj))>GK(tord)))) exit
      ord(jj+1) = ord(jj)
      jj = jj - 1
    enddo
    ord(jj+1) = tord
  enddo
 
end subroutine sortix_gvec_insertsort
subroutine sortix_gvec(NVAL, AA, ord, gvec)
  integer, intent(in) :: NVAL
  integer, intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer, intent(in) :: gvec(3,NVAL) !< (3, N) G-vectors, used to break tie for equal data
  integer :: GK(NVAL)
 
  call get_GK_array_from_gvecs(NVAL, gvec, GK)
  call sortix_gvec_GK(NVAL, AA, ord, GK)
 
end subroutine sortix_gvec
!> Sorts (actually, ranks) a real/integer array AA.
!! The rank is written to the output array ord.
!! This subroutine is based on the routine MRGRNK by Michel Olagnon, which
!! uses the merge sort algorithm.
subroutine sortix_gvec_GK(NVAL, AA, ord&
, GK&
)
  integer, intent(in) :: NVAL
  integer, intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer, intent(in) :: GK(NVAL)
  integer, parameter :: TOL=TOL_ZERO
  integer :: JT(NVAL)
  integer :: LMTNA, LMTNC, IRNG1, IRNG2
  integer :: IIND, ID, IWRK, IWRKF, JINDA, IA, IB
 
!
! Fill-in the index array, creating ordered couples
!
  Do IIND = 2, NVAL, 2
    If (&
(AA(IIND)-AA(IIND-1)>TOL.or.(AA(IIND)-AA(IIND-1)>-TOL.and.GK(IIND)>GK(IIND-1)))&
    ) Then
      ord (IIND-1) = IIND - 1
      ord (IIND) = IIND
    Else
      ord (IIND-1) = IIND
      ord (IIND) = IIND - 1
    End If
  End Do
  If (Modulo(NVAL, 2) /= 0) Then
    ord (NVAL) = NVAL
  End If
  ! FHJ - shortcut if the array is small enough
  if (NVAL<16) then
    call sortix_gvec_insertsort(NVAL, AA, ord&
, GK&
)
   
    return
  endif
!
! We will now have ordered subsets A - B - A - B - ...
! and merge A and B couples into C - C - ...
!
  LMTNA = 2
  LMTNC = 4
!
! First iteration. The length of the ordered subsets goes from 2 to 4
!
  Do
    If (NVAL <= 2) Exit
!
! Loop on merges of A and B into C
!
    Do ID = 0, NVAL - 1, 4
      If ((ID+4) > NVAL) Then
        If ((ID+2) >= NVAL) Exit
!
! 1 2 3
!
        If (&
(AA(ord(ID+3))-AA(ord(ID+2))>TOL.or.(AA(ord(ID+3))-AA(ord(ID+2))>-TOL.and.GK(ord(ID+3))>GK(ord(ID+2))))&
        ) Exit
!
! 1 3 2
!
        If (&
(AA(ord(ID+3))-AA(ord(ID+1))>TOL.or.(AA(ord(ID+3))-AA(ord(ID+1))>-TOL.and.GK(ord(ID+3))>GK(ord(ID+1))))&
        ) Then
          IRNG2 = ord (ID+2)
          ord (ID+2) = ord (ID+3)
          ord (ID+3) = IRNG2
!
! 3 1 2
!
        Else
          IRNG1 = ord (ID+1)
          ord (ID+1) = ord (ID+3)
          ord (ID+3) = ord (ID+2)
          ord (ID+2) = IRNG1
        End If
        Exit
      End If
!
! 1 2 3 4
!
      If (&
(AA(ord(ID+3))-AA(ord(ID+2))>TOL.or.(AA(ord(ID+3))-AA(ord(ID+2))>-TOL.and.GK(ord(ID+3))>GK(ord(ID+2))))&
      ) Cycle
!
! 1 3 x x
!
      If (&
(AA(ord(ID+3))-AA(ord(ID+1))>TOL.or.(AA(ord(ID+3))-AA(ord(ID+1))>-TOL.and.GK(ord(ID+3))>GK(ord(ID+1))))&
      ) Then
        IRNG2 = ord (ID+2)
        ord (ID+2) = ord (ID+3)
        If (&
(AA(ord(ID+4))-AA(IRNG2)>TOL.or.(AA(ord(ID+4))-AA(IRNG2)>-TOL.and.GK(ord(ID+4))>GK(IRNG2)))&
        ) Then
! 1 3 2 4
          ord (ID+3) = IRNG2
        Else
! 1 3 4 2
          ord (ID+3) = ord (ID+4)
          ord (ID+4) = IRNG2
        End If
!
! 3 x x x
!
      Else
        IRNG1 = ord (ID+1)
        IRNG2 = ord (ID+2)
        ord (ID+1) = ord (ID+3)
        If (&
(AA(ord(ID+4))-AA(IRNG1)>TOL.or.(AA(ord(ID+4))-AA(IRNG1)>-TOL.and.GK(ord(ID+4))>GK(IRNG1)))&
        ) Then
          ord (ID+2) = IRNG1
          If (&
(AA(ord(ID+4))-AA(IRNG2)>TOL.or.(AA(ord(ID+4))-AA(IRNG2)>-TOL.and.GK(ord(ID+4))>GK(IRNG2)))&
          ) Then
! 3 1 2 4
            ord (ID+3) = IRNG2
          Else
! 3 1 4 2
            ord (ID+3) = ord (ID+4)
            ord (ID+4) = IRNG2
          End If
        Else
! 3 4 1 2
          ord (ID+2) = ord (ID+4)
          ord (ID+3) = IRNG1
          ord (ID+4) = IRNG2
        End If
      End If
    End Do
!
! The Cs become As and Bs
!
    LMTNA = 4
    Exit
  End Do
!
! Iteration loop. Each time, the length of the ordered subsets
! is doubled.
!
  Do
    If (LMTNA >= NVAL) Exit
    IWRKF = 0
    LMTNC = 2 * LMTNC
!
! Loop on merges of A and B into C
!
    Do
      IWRK = IWRKF
      ID = IWRKF + 1
      JINDA = IWRKF + LMTNA
      IWRKF = IWRKF + LMTNC
      If (IWRKF >= NVAL) Then
        If (JINDA >= NVAL) Exit
        IWRKF = NVAL
      End If
      IA = 1
      IB = JINDA + 1
!
! Shortcut for the case when the max of A is smaller
! than the min of B. This line may be activated when the
! initial set is already close to sorted.
!
      IF (&
(AA(ord(IB))-AA(ord(JINDA))>TOL.or.(AA(ord(IB))-AA(ord(JINDA))>-TOL.and.GK(ord(IB))>GK(ord(JINDA))))&
      ) CYCLE
!
! One steps in the C subset, that we build in the final rank array
!
! Make a copy of the rank array for the merge iteration
!
      JT (1:LMTNA) = ord (ID:JINDA)
!
      Do
        IWRK = IWRK + 1
!
! We still have unprocessed values in both A and B
!
        If (&
(AA(JT(IA))-AA(ord(IB))>TOL.or.(AA(JT(IA))-AA(ord(IB))>-TOL.and.GK(JT(IA))>GK(ord(IB))))&
        ) Then
          ord (IWRK) = ord (IB)
          IB = IB + 1
          If (IB > IWRKF) Then
! Only A still with unprocessed values
            ord (IWRK+1:IWRKF) = JT (IA:LMTNA)
            Exit
          End If
        Else
          ord (IWRK) = JT (IA)
          IA = IA + 1
          If (IA > LMTNA) Exit! Only B still with unprocessed values
        End If
!
      End Do
    End Do
!
! The Cs become As and Bs
!
    LMTNA = 2 * LMTNA
  End Do
!
 
!
End Subroutine sortix_gvec_GK
subroutine sortix_gvec_threaded(NVAL, AA, ord, gvec)
  integer, intent(in) :: NVAL
  integer, intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer, intent(in) :: gvec(3,NVAL) !< (3, N) G-vectors, used to break tie for equal data
  integer :: GK(NVAL)
 
  call get_GK_array_from_gvecs(NVAL, gvec, GK)
  call sortix_gvec_threaded_GK(NVAL, AA, ord, GK)
 
end subroutine sortix_gvec_threaded
!> Threaded sort (actually, ranking) of a real/integer array AA.
!! Public code adapted from:
!! https://github.com/cphyc/Fortran-parallel-sort/blob/master/mod_sort.f90
subroutine sortix_gvec_threaded_GK(NVAL, AA, ord&
, GK&
)
  integer, intent(in) :: NVAL
  integer, intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer, intent(in) :: GK(NVAL)
  integer :: nthreads
 
  ! FHJ: TODO - nthreads should be different that peinf%nthreads_sort if we
  ! have nested OMPs. However, starting a disabled PARALLEL region to figure out
  ! th number of threads is just too much overhead.
  nthreads = peinf%nthreads_sort
  if (nthreads==1.or.NVAL<nthreads*100) then
    call sortix_gvec_GK(NVAL, AA, ord&
, GK&
)
   
    return
  endif
 
end subroutine sortix_gvec_threaded_GK
! Merge two parts of A, ordered by order from left to right around middle.
subroutine sortix_gvec_threaded_merge(NVAL, AA, ord, &
GK, &
left, middle, right)
  integer, intent(in) :: NVAL
  integer, intent(in) :: AA(NVAL)
  integer, intent(inout) :: ord(NVAL)
  integer, intent(in) :: GK(NVAL)
  integer, intent(in) :: left, middle, right
  integer :: leftA, rightA, leftB, rightB
  integer :: iA, iB, i
  integer :: lenA, lenB
  integer :: idxA, idxB
  integer :: orderA(left:middle)
  integer :: orderB(middle+1:right)
  integer, parameter :: TOL=TOL_ZERO
 
  ! copy order
  orderA = ord(left:middle)
  orderB = ord(middle+1:right)
  ! more explicit variables
  leftA = left
  rightA = middle
  leftB = middle+1
  rightB = right
  ! initialize iA, iB to their leftmost position
  iA = leftA
  iB = leftB
  i = leftA
  do while ((iA <= rightA) .and. (iB <= rightB))
    idxA = orderA(iA)
    idxB = orderB(iB)
    if ((AA(idxA)-AA(idxB)>TOL.or.(AA(idxA)-AA(idxB)>-TOL.and.GK(idxA)>GK(idxB)))) then
       ord(i) = idxB
       iB = iB + 1
    else
       ord(i) = idxA
       iA = iA + 1
    end if
    i = i + 1
  end do
  ! either A or B still have elements, append them to the new order
  do while (iA <= rightA)
    ord(i) = orderA(iA)
    iA = iA + 1
    i = i + 1
  end do
  do while (iB <= rightB)
    ord(i) = orderB(iB)
    iB = iB + 1
    i = i + 1
  end do
 
end subroutine sortix_gvec_threaded_merge
! This file is based on the work of Michel Olagnon.
! The original code for the MRGRNK subroutine is available at:
! http://fortran-2000.com/rank/mrgrnk.f90
! MRGRNK - Copyright (c) Michel Olagnon
! Copying and distribution of this file, with or without modification,
! are permitted in any medium without royalty provided the copyright
! notice and this notice are preserved. This file is offered as-is,
! without any warranty.
! FHJ: WARNING - make sure you don`t change the following lines too much,
! otherwise they will be longer than 120 characters after the preprocessors kicks in.
! Note that, if there the extra "gvec" argument, we use a tolerance to figure
! out if the two items AA(ii) and AA(jj) are degenerate.
! FHJ: Note - we need to nest the "JOIN" macro otherwise the symbol sortix_no_gvec
! doesn`t get expanded by the C preprocessor.
! LABEL_GK is the sorting function without automatically generating the GK
! array out of the gvec array. This is useful for threaded sorting, when we
! first generate the GK array, but then call the sorting routine many times.
! LABEL_GK is the same as sortix_no_gvec if there is no gvec input.
! This is the kernel of the sourting routine
! These functions are threaded wrappers for the sorting routines
!> Sorts (actually, ranks) a small array AA using the insert sort algorithm.
subroutine sortix_no_gvec_insertsort(NVAL, AA, ord&
)
  integer, intent(in) :: NVAL
  integer, intent(in) :: AA(NVAL)
  integer, intent(inout) :: ord(NVAL)
  integer :: ii, jj, tord
 
  do ii = 2, NVAL
    tord = ord(ii)
    jj = ii - 1
    do while (jj>0)
      if (.not.(AA(ord(jj))>AA(tord))) exit
      ord(jj+1) = ord(jj)
      jj = jj - 1
    enddo
    ord(jj+1) = tord
  enddo
 
end subroutine sortix_no_gvec_insertsort
!> Sorts (actually, ranks) a real/integer array AA.
!! The rank is written to the output array ord.
!! This subroutine is based on the routine MRGRNK by Michel Olagnon, which
!! uses the merge sort algorithm.
subroutine sortix_no_gvec(NVAL, AA, ord&
)
  integer, intent(in) :: NVAL
  integer, intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer :: JT(NVAL)
  integer :: LMTNA, LMTNC, IRNG1, IRNG2
  integer :: IIND, ID, IWRK, IWRKF, JINDA, IA, IB
 
!
! Fill-in the index array, creating ordered couples
!
  Do IIND = 2, NVAL, 2
    If (&
(AA(IIND)>AA(IIND-1))&
    ) Then
      ord (IIND-1) = IIND - 1
      ord (IIND) = IIND
    Else
      ord (IIND-1) = IIND
      ord (IIND) = IIND - 1
    End If
  End Do
  If (Modulo(NVAL, 2) /= 0) Then
    ord (NVAL) = NVAL
  End If
  ! FHJ - shortcut if the array is small enough
  if (NVAL<16) then
    call sortix_no_gvec_insertsort(NVAL, AA, ord&
)
   
    return
  endif
!
! We will now have ordered subsets A - B - A - B - ...
! and merge A and B couples into C - C - ...
!
  LMTNA = 2
  LMTNC = 4
!
! First iteration. The length of the ordered subsets goes from 2 to 4
!
  Do
    If (NVAL <= 2) Exit
!
! Loop on merges of A and B into C
!
    Do ID = 0, NVAL - 1, 4
      If ((ID+4) > NVAL) Then
        If ((ID+2) >= NVAL) Exit
!
! 1 2 3
!
        If (&
(AA(ord(ID+3))>AA(ord(ID+2)))&
        ) Exit
!
! 1 3 2
!
        If (&
(AA(ord(ID+3))>AA(ord(ID+1)))&
        ) Then
          IRNG2 = ord (ID+2)
          ord (ID+2) = ord (ID+3)
          ord (ID+3) = IRNG2
!
! 3 1 2
!
        Else
          IRNG1 = ord (ID+1)
          ord (ID+1) = ord (ID+3)
          ord (ID+3) = ord (ID+2)
          ord (ID+2) = IRNG1
        End If
        Exit
      End If
!
! 1 2 3 4
!
      If (&
(AA(ord(ID+3))>AA(ord(ID+2)))&
      ) Cycle
!
! 1 3 x x
!
      If (&
(AA(ord(ID+3))>AA(ord(ID+1)))&
      ) Then
        IRNG2 = ord (ID+2)
        ord (ID+2) = ord (ID+3)
        If (&
(AA(ord(ID+4))>AA(IRNG2))&
        ) Then
! 1 3 2 4
          ord (ID+3) = IRNG2
        Else
! 1 3 4 2
          ord (ID+3) = ord (ID+4)
          ord (ID+4) = IRNG2
        End If
!
! 3 x x x
!
      Else
        IRNG1 = ord (ID+1)
        IRNG2 = ord (ID+2)
        ord (ID+1) = ord (ID+3)
        If (&
(AA(ord(ID+4))>AA(IRNG1))&
        ) Then
          ord (ID+2) = IRNG1
          If (&
(AA(ord(ID+4))>AA(IRNG2))&
          ) Then
! 3 1 2 4
            ord (ID+3) = IRNG2
          Else
! 3 1 4 2
            ord (ID+3) = ord (ID+4)
            ord (ID+4) = IRNG2
          End If
        Else
! 3 4 1 2
          ord (ID+2) = ord (ID+4)
          ord (ID+3) = IRNG1
          ord (ID+4) = IRNG2
        End If
      End If
    End Do
!
! The Cs become As and Bs
!
    LMTNA = 4
    Exit
  End Do
!
! Iteration loop. Each time, the length of the ordered subsets
! is doubled.
!
  Do
    If (LMTNA >= NVAL) Exit
    IWRKF = 0
    LMTNC = 2 * LMTNC
!
! Loop on merges of A and B into C
!
    Do
      IWRK = IWRKF
      ID = IWRKF + 1
      JINDA = IWRKF + LMTNA
      IWRKF = IWRKF + LMTNC
      If (IWRKF >= NVAL) Then
        If (JINDA >= NVAL) Exit
        IWRKF = NVAL
      End If
      IA = 1
      IB = JINDA + 1
!
! Shortcut for the case when the max of A is smaller
! than the min of B. This line may be activated when the
! initial set is already close to sorted.
!
      IF (&
(AA(ord(IB))>AA(ord(JINDA)))&
      ) CYCLE
!
! One steps in the C subset, that we build in the final rank array
!
! Make a copy of the rank array for the merge iteration
!
      JT (1:LMTNA) = ord (ID:JINDA)
!
      Do
        IWRK = IWRK + 1
!
! We still have unprocessed values in both A and B
!
        If (&
(AA(JT(IA))>AA(ord(IB)))&
        ) Then
          ord (IWRK) = ord (IB)
          IB = IB + 1
          If (IB > IWRKF) Then
! Only A still with unprocessed values
            ord (IWRK+1:IWRKF) = JT (IA:LMTNA)
            Exit
          End If
        Else
          ord (IWRK) = JT (IA)
          IA = IA + 1
          If (IA > LMTNA) Exit! Only B still with unprocessed values
        End If
!
      End Do
    End Do
!
! The Cs become As and Bs
!
    LMTNA = 2 * LMTNA
  End Do
!
 
!
End Subroutine sortix_no_gvec
!> Threaded sort (actually, ranking) of a real/integer array AA.
!! Public code adapted from:
!! https://github.com/cphyc/Fortran-parallel-sort/blob/master/mod_sort.f90
subroutine sortix_no_gvec_threaded(NVAL, AA, ord&
)
  integer, intent(in) :: NVAL
  integer, intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer :: nthreads
 
  ! FHJ: TODO - nthreads should be different that peinf%nthreads_sort if we
  ! have nested OMPs. However, starting a disabled PARALLEL region to figure out
  ! th number of threads is just too much overhead.
  nthreads = peinf%nthreads_sort
  if (nthreads==1.or.NVAL<nthreads*100) then
    call sortix_no_gvec(NVAL, AA, ord&
)
   
    return
  endif
 
end subroutine sortix_no_gvec_threaded
! Merge two parts of A, ordered by order from left to right around middle.
subroutine sortix_no_gvec_threaded_merge(NVAL, AA, ord, &
left, middle, right)
  integer, intent(in) :: NVAL
  integer, intent(in) :: AA(NVAL)
  integer, intent(inout) :: ord(NVAL)
  integer, intent(in) :: left, middle, right
  integer :: leftA, rightA, leftB, rightB
  integer :: iA, iB, i
  integer :: lenA, lenB
  integer :: idxA, idxB
  integer :: orderA(left:middle)
  integer :: orderB(middle+1:right)
 
  ! copy order
  orderA = ord(left:middle)
  orderB = ord(middle+1:right)
  ! more explicit variables
  leftA = left
  rightA = middle
  leftB = middle+1
  rightB = right
  ! initialize iA, iB to their leftmost position
  iA = leftA
  iB = leftB
  i = leftA
  do while ((iA <= rightA) .and. (iB <= rightB))
    idxA = orderA(iA)
    idxB = orderB(iB)
    if ((AA(idxA)>AA(idxB))) then
       ord(i) = idxB
       iB = iB + 1
    else
       ord(i) = idxA
       iA = iA + 1
    end if
    i = i + 1
  end do
  ! either A or B still have elements, append them to the new order
  do while (iA <= rightA)
    ord(i) = orderA(iA)
    iA = iA + 1
    i = i + 1
  end do
  do while (iB <= rightB)
    ord(i) = orderB(iB)
    iB = iB + 1
    i = i + 1
  end do
 
end subroutine sortix_no_gvec_threaded_merge
! This file is based on the work of Michel Olagnon.
! The original code for the MRGRNK subroutine is available at:
! http://fortran-2000.com/rank/mrgrnk.f90
! MRGRNK - Copyright (c) Michel Olagnon
! Copying and distribution of this file, with or without modification,
! are permitted in any medium without royalty provided the copyright
! notice and this notice are preserved. This file is offered as-is,
! without any warranty.
! FHJ: WARNING - make sure you don`t change the following lines too much,
! otherwise they will be longer than 120 characters after the preprocessors kicks in.
! Note that, if there the extra "gvec" argument, we use a tolerance to figure
! out if the two items AA(ii) and AA(jj) are degenerate.
! FHJ: Note - we need to nest the "JOIN" macro otherwise the symbol sortrx_gvec
! doesn`t get expanded by the C preprocessor.
! LABEL_GK is the sorting function without automatically generating the GK
! array out of the gvec array. This is useful for threaded sorting, when we
! first generate the GK array, but then call the sorting routine many times.
! LABEL_GK is the same as sortrx_gvec if there is no gvec input.
! This is the kernel of the sourting routine
! These functions are threaded wrappers for the sorting routines
!> Sorts (actually, ranks) a small array AA using the insert sort algorithm.
subroutine sortrx_gvec_insertsort(NVAL, AA, ord&
, GK&
)
  integer, intent(in) :: NVAL
  real(DP), intent(in) :: AA(NVAL)
  integer, intent(inout) :: ord(NVAL)
  integer, intent(in) :: GK(NVAL)
  real(DP), parameter :: TOL=TOL_ZERO
  integer :: ii, jj, tord
 
  do ii = 2, NVAL
    tord = ord(ii)
    jj = ii - 1
    do while (jj>0)
      if (.not.(AA(ord(jj))-AA(tord)>TOL.or.(AA(ord(jj))-AA(tord)>-TOL.and.GK(ord(jj))>GK(tord)))) exit
      ord(jj+1) = ord(jj)
      jj = jj - 1
    enddo
    ord(jj+1) = tord
  enddo
 
end subroutine sortrx_gvec_insertsort
subroutine sortrx_gvec(NVAL, AA, ord, gvec)
  integer, intent(in) :: NVAL
  real(DP), intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer, intent(in) :: gvec(3,NVAL) !< (3, N) G-vectors, used to break tie for equal data
  integer :: GK(NVAL)
 
  call get_GK_array_from_gvecs(NVAL, gvec, GK)
  call sortrx_gvec_GK(NVAL, AA, ord, GK)
 
end subroutine sortrx_gvec
!> Sorts (actually, ranks) a real/integer array AA.
!! The rank is written to the output array ord.
!! This subroutine is based on the routine MRGRNK by Michel Olagnon, which
!! uses the merge sort algorithm.
subroutine sortrx_gvec_GK(NVAL, AA, ord&
, GK&
)
  integer, intent(in) :: NVAL
  real(DP), intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer, intent(in) :: GK(NVAL)
  real(DP), parameter :: TOL=TOL_ZERO
  integer :: JT(NVAL)
  integer :: LMTNA, LMTNC, IRNG1, IRNG2
  integer :: IIND, ID, IWRK, IWRKF, JINDA, IA, IB
 
!
! Fill-in the index array, creating ordered couples
!
  Do IIND = 2, NVAL, 2
    If (&
(AA(IIND)-AA(IIND-1)>TOL.or.(AA(IIND)-AA(IIND-1)>-TOL.and.GK(IIND)>GK(IIND-1)))&
    ) Then
      ord (IIND-1) = IIND - 1
      ord (IIND) = IIND
    Else
      ord (IIND-1) = IIND
      ord (IIND) = IIND - 1
    End If
  End Do
  If (Modulo(NVAL, 2) /= 0) Then
    ord (NVAL) = NVAL
  End If
  ! FHJ - shortcut if the array is small enough
  if (NVAL<16) then
    call sortrx_gvec_insertsort(NVAL, AA, ord&
, GK&
)
   
    return
  endif
!
! We will now have ordered subsets A - B - A - B - ...
! and merge A and B couples into C - C - ...
!
  LMTNA = 2
  LMTNC = 4
!
! First iteration. The length of the ordered subsets goes from 2 to 4
!
  Do
    If (NVAL <= 2) Exit
!
! Loop on merges of A and B into C
!
    Do ID = 0, NVAL - 1, 4
      If ((ID+4) > NVAL) Then
        If ((ID+2) >= NVAL) Exit
!
! 1 2 3
!
        If (&
(AA(ord(ID+3))-AA(ord(ID+2))>TOL.or.(AA(ord(ID+3))-AA(ord(ID+2))>-TOL.and.GK(ord(ID+3))>GK(ord(ID+2))))&
        ) Exit
!
! 1 3 2
!
        If (&
(AA(ord(ID+3))-AA(ord(ID+1))>TOL.or.(AA(ord(ID+3))-AA(ord(ID+1))>-TOL.and.GK(ord(ID+3))>GK(ord(ID+1))))&
        ) Then
          IRNG2 = ord (ID+2)
          ord (ID+2) = ord (ID+3)
          ord (ID+3) = IRNG2
!
! 3 1 2
!
        Else
          IRNG1 = ord (ID+1)
          ord (ID+1) = ord (ID+3)
          ord (ID+3) = ord (ID+2)
          ord (ID+2) = IRNG1
        End If
        Exit
      End If
!
! 1 2 3 4
!
      If (&
(AA(ord(ID+3))-AA(ord(ID+2))>TOL.or.(AA(ord(ID+3))-AA(ord(ID+2))>-TOL.and.GK(ord(ID+3))>GK(ord(ID+2))))&
      ) Cycle
!
! 1 3 x x
!
      If (&
(AA(ord(ID+3))-AA(ord(ID+1))>TOL.or.(AA(ord(ID+3))-AA(ord(ID+1))>-TOL.and.GK(ord(ID+3))>GK(ord(ID+1))))&
      ) Then
        IRNG2 = ord (ID+2)
        ord (ID+2) = ord (ID+3)
        If (&
(AA(ord(ID+4))-AA(IRNG2)>TOL.or.(AA(ord(ID+4))-AA(IRNG2)>-TOL.and.GK(ord(ID+4))>GK(IRNG2)))&
        ) Then
! 1 3 2 4
          ord (ID+3) = IRNG2
        Else
! 1 3 4 2
          ord (ID+3) = ord (ID+4)
          ord (ID+4) = IRNG2
        End If
!
! 3 x x x
!
      Else
        IRNG1 = ord (ID+1)
        IRNG2 = ord (ID+2)
        ord (ID+1) = ord (ID+3)
        If (&
(AA(ord(ID+4))-AA(IRNG1)>TOL.or.(AA(ord(ID+4))-AA(IRNG1)>-TOL.and.GK(ord(ID+4))>GK(IRNG1)))&
        ) Then
          ord (ID+2) = IRNG1
          If (&
(AA(ord(ID+4))-AA(IRNG2)>TOL.or.(AA(ord(ID+4))-AA(IRNG2)>-TOL.and.GK(ord(ID+4))>GK(IRNG2)))&
          ) Then
! 3 1 2 4
            ord (ID+3) = IRNG2
          Else
! 3 1 4 2
            ord (ID+3) = ord (ID+4)
            ord (ID+4) = IRNG2
          End If
        Else
! 3 4 1 2
          ord (ID+2) = ord (ID+4)
          ord (ID+3) = IRNG1
          ord (ID+4) = IRNG2
        End If
      End If
    End Do
!
! The Cs become As and Bs
!
    LMTNA = 4
    Exit
  End Do
!
! Iteration loop. Each time, the length of the ordered subsets
! is doubled.
!
  Do
    If (LMTNA >= NVAL) Exit
    IWRKF = 0
    LMTNC = 2 * LMTNC
!
! Loop on merges of A and B into C
!
    Do
      IWRK = IWRKF
      ID = IWRKF + 1
      JINDA = IWRKF + LMTNA
      IWRKF = IWRKF + LMTNC
      If (IWRKF >= NVAL) Then
        If (JINDA >= NVAL) Exit
        IWRKF = NVAL
      End If
      IA = 1
      IB = JINDA + 1
!
! Shortcut for the case when the max of A is smaller
! than the min of B. This line may be activated when the
! initial set is already close to sorted.
!
      IF (&
(AA(ord(IB))-AA(ord(JINDA))>TOL.or.(AA(ord(IB))-AA(ord(JINDA))>-TOL.and.GK(ord(IB))>GK(ord(JINDA))))&
      ) CYCLE
!
! One steps in the C subset, that we build in the final rank array
!
! Make a copy of the rank array for the merge iteration
!
      JT (1:LMTNA) = ord (ID:JINDA)
!
      Do
        IWRK = IWRK + 1
!
! We still have unprocessed values in both A and B
!
        If (&
(AA(JT(IA))-AA(ord(IB))>TOL.or.(AA(JT(IA))-AA(ord(IB))>-TOL.and.GK(JT(IA))>GK(ord(IB))))&
        ) Then
          ord (IWRK) = ord (IB)
          IB = IB + 1
          If (IB > IWRKF) Then
! Only A still with unprocessed values
            ord (IWRK+1:IWRKF) = JT (IA:LMTNA)
            Exit
          End If
        Else
          ord (IWRK) = JT (IA)
          IA = IA + 1
          If (IA > LMTNA) Exit! Only B still with unprocessed values
        End If
!
      End Do
    End Do
!
! The Cs become As and Bs
!
    LMTNA = 2 * LMTNA
  End Do
!
 
!
End Subroutine sortrx_gvec_GK
subroutine sortrx_gvec_threaded(NVAL, AA, ord, gvec)
  integer, intent(in) :: NVAL
  real(DP), intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer, intent(in) :: gvec(3,NVAL) !< (3, N) G-vectors, used to break tie for equal data
  integer :: GK(NVAL)
 
  call get_GK_array_from_gvecs(NVAL, gvec, GK)
  call sortrx_gvec_threaded_GK(NVAL, AA, ord, GK)
 
end subroutine sortrx_gvec_threaded
!> Threaded sort (actually, ranking) of a real/integer array AA.
!! Public code adapted from:
!! https://github.com/cphyc/Fortran-parallel-sort/blob/master/mod_sort.f90
subroutine sortrx_gvec_threaded_GK(NVAL, AA, ord&
, GK&
)
  integer, intent(in) :: NVAL
  real(DP), intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer, intent(in) :: GK(NVAL)
  integer :: nthreads
 
  ! FHJ: TODO - nthreads should be different that peinf%nthreads_sort if we
  ! have nested OMPs. However, starting a disabled PARALLEL region to figure out
  ! th number of threads is just too much overhead.
  nthreads = peinf%nthreads_sort
  if (nthreads==1.or.NVAL<nthreads*100) then
    call sortrx_gvec_GK(NVAL, AA, ord&
, GK&
)
   
    return
  endif
 
end subroutine sortrx_gvec_threaded_GK
! Merge two parts of A, ordered by order from left to right around middle.
subroutine sortrx_gvec_threaded_merge(NVAL, AA, ord, &
GK, &
left, middle, right)
  integer, intent(in) :: NVAL
  real(DP), intent(in) :: AA(NVAL)
  integer, intent(inout) :: ord(NVAL)
  integer, intent(in) :: GK(NVAL)
  integer, intent(in) :: left, middle, right
  integer :: leftA, rightA, leftB, rightB
  integer :: iA, iB, i
  integer :: lenA, lenB
  integer :: idxA, idxB
  integer :: orderA(left:middle)
  integer :: orderB(middle+1:right)
  real(DP), parameter :: TOL=TOL_ZERO
 
  ! copy order
  orderA = ord(left:middle)
  orderB = ord(middle+1:right)
  ! more explicit variables
  leftA = left
  rightA = middle
  leftB = middle+1
  rightB = right
  ! initialize iA, iB to their leftmost position
  iA = leftA
  iB = leftB
  i = leftA
  do while ((iA <= rightA) .and. (iB <= rightB))
    idxA = orderA(iA)
    idxB = orderB(iB)
    if ((AA(idxA)-AA(idxB)>TOL.or.(AA(idxA)-AA(idxB)>-TOL.and.GK(idxA)>GK(idxB)))) then
       ord(i) = idxB
       iB = iB + 1
    else
       ord(i) = idxA
       iA = iA + 1
    end if
    i = i + 1
  end do
  ! either A or B still have elements, append them to the new order
  do while (iA <= rightA)
    ord(i) = orderA(iA)
    iA = iA + 1
    i = i + 1
  end do
  do while (iB <= rightB)
    ord(i) = orderB(iB)
    iB = iB + 1
    i = i + 1
  end do
 
end subroutine sortrx_gvec_threaded_merge
! This file is based on the work of Michel Olagnon.
! The original code for the MRGRNK subroutine is available at:
! http://fortran-2000.com/rank/mrgrnk.f90
! MRGRNK - Copyright (c) Michel Olagnon
! Copying and distribution of this file, with or without modification,
! are permitted in any medium without royalty provided the copyright
! notice and this notice are preserved. This file is offered as-is,
! without any warranty.
! FHJ: WARNING - make sure you don`t change the following lines too much,
! otherwise they will be longer than 120 characters after the preprocessors kicks in.
! Note that, if there the extra "gvec" argument, we use a tolerance to figure
! out if the two items AA(ii) and AA(jj) are degenerate.
! FHJ: Note - we need to nest the "JOIN" macro otherwise the symbol sortrx_no_gvec
! doesn`t get expanded by the C preprocessor.
! LABEL_GK is the sorting function without automatically generating the GK
! array out of the gvec array. This is useful for threaded sorting, when we
! first generate the GK array, but then call the sorting routine many times.
! LABEL_GK is the same as sortrx_no_gvec if there is no gvec input.
! This is the kernel of the sourting routine
! These functions are threaded wrappers for the sorting routines
!> Sorts (actually, ranks) a small array AA using the insert sort algorithm.
subroutine sortrx_no_gvec_insertsort(NVAL, AA, ord&
)
  integer, intent(in) :: NVAL
  real(DP), intent(in) :: AA(NVAL)
  integer, intent(inout) :: ord(NVAL)
  integer :: ii, jj, tord
 
  do ii = 2, NVAL
    tord = ord(ii)
    jj = ii - 1
    do while (jj>0)
      if (.not.(AA(ord(jj))>AA(tord))) exit
      ord(jj+1) = ord(jj)
      jj = jj - 1
    enddo
    ord(jj+1) = tord
  enddo
 
end subroutine sortrx_no_gvec_insertsort
!> Sorts (actually, ranks) a real/integer array AA.
!! The rank is written to the output array ord.
!! This subroutine is based on the routine MRGRNK by Michel Olagnon, which
!! uses the merge sort algorithm.
subroutine sortrx_no_gvec(NVAL, AA, ord&
)
  integer, intent(in) :: NVAL
  real(DP), intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer :: JT(NVAL)
  integer :: LMTNA, LMTNC, IRNG1, IRNG2
  integer :: IIND, ID, IWRK, IWRKF, JINDA, IA, IB
 
!
! Fill-in the index array, creating ordered couples
!
  Do IIND = 2, NVAL, 2
    If (&
(AA(IIND)>AA(IIND-1))&
    ) Then
      ord (IIND-1) = IIND - 1
      ord (IIND) = IIND
    Else
      ord (IIND-1) = IIND
      ord (IIND) = IIND - 1
    End If
  End Do
  If (Modulo(NVAL, 2) /= 0) Then
    ord (NVAL) = NVAL
  End If
  ! FHJ - shortcut if the array is small enough
  if (NVAL<16) then
    call sortrx_no_gvec_insertsort(NVAL, AA, ord&
)
   
    return
  endif
!
! We will now have ordered subsets A - B - A - B - ...
! and merge A and B couples into C - C - ...
!
  LMTNA = 2
  LMTNC = 4
!
! First iteration. The length of the ordered subsets goes from 2 to 4
!
  Do
    If (NVAL <= 2) Exit
!
! Loop on merges of A and B into C
!
    Do ID = 0, NVAL - 1, 4
      If ((ID+4) > NVAL) Then
        If ((ID+2) >= NVAL) Exit
!
! 1 2 3
!
        If (&
(AA(ord(ID+3))>AA(ord(ID+2)))&
        ) Exit
!
! 1 3 2
!
        If (&
(AA(ord(ID+3))>AA(ord(ID+1)))&
        ) Then
          IRNG2 = ord (ID+2)
          ord (ID+2) = ord (ID+3)
          ord (ID+3) = IRNG2
!
! 3 1 2
!
        Else
          IRNG1 = ord (ID+1)
          ord (ID+1) = ord (ID+3)
          ord (ID+3) = ord (ID+2)
          ord (ID+2) = IRNG1
        End If
        Exit
      End If
!
! 1 2 3 4
!
      If (&
(AA(ord(ID+3))>AA(ord(ID+2)))&
      ) Cycle
!
! 1 3 x x
!
      If (&
(AA(ord(ID+3))>AA(ord(ID+1)))&
      ) Then
        IRNG2 = ord (ID+2)
        ord (ID+2) = ord (ID+3)
        If (&
(AA(ord(ID+4))>AA(IRNG2))&
        ) Then
! 1 3 2 4
          ord (ID+3) = IRNG2
        Else
! 1 3 4 2
          ord (ID+3) = ord (ID+4)
          ord (ID+4) = IRNG2
        End If
!
! 3 x x x
!
      Else
        IRNG1 = ord (ID+1)
        IRNG2 = ord (ID+2)
        ord (ID+1) = ord (ID+3)
        If (&
(AA(ord(ID+4))>AA(IRNG1))&
        ) Then
          ord (ID+2) = IRNG1
          If (&
(AA(ord(ID+4))>AA(IRNG2))&
          ) Then
! 3 1 2 4
            ord (ID+3) = IRNG2
          Else
! 3 1 4 2
            ord (ID+3) = ord (ID+4)
            ord (ID+4) = IRNG2
          End If
        Else
! 3 4 1 2
          ord (ID+2) = ord (ID+4)
          ord (ID+3) = IRNG1
          ord (ID+4) = IRNG2
        End If
      End If
    End Do
!
! The Cs become As and Bs
!
    LMTNA = 4
    Exit
  End Do
!
! Iteration loop. Each time, the length of the ordered subsets
! is doubled.
!
  Do
    If (LMTNA >= NVAL) Exit
    IWRKF = 0
    LMTNC = 2 * LMTNC
!
! Loop on merges of A and B into C
!
    Do
      IWRK = IWRKF
      ID = IWRKF + 1
      JINDA = IWRKF + LMTNA
      IWRKF = IWRKF + LMTNC
      If (IWRKF >= NVAL) Then
        If (JINDA >= NVAL) Exit
        IWRKF = NVAL
      End If
      IA = 1
      IB = JINDA + 1
!
! Shortcut for the case when the max of A is smaller
! than the min of B. This line may be activated when the
! initial set is already close to sorted.
!
      IF (&
(AA(ord(IB))>AA(ord(JINDA)))&
      ) CYCLE
!
! One steps in the C subset, that we build in the final rank array
!
! Make a copy of the rank array for the merge iteration
!
      JT (1:LMTNA) = ord (ID:JINDA)
!
      Do
        IWRK = IWRK + 1
!
! We still have unprocessed values in both A and B
!
        If (&
(AA(JT(IA))>AA(ord(IB)))&
        ) Then
          ord (IWRK) = ord (IB)
          IB = IB + 1
          If (IB > IWRKF) Then
! Only A still with unprocessed values
            ord (IWRK+1:IWRKF) = JT (IA:LMTNA)
            Exit
          End If
        Else
          ord (IWRK) = JT (IA)
          IA = IA + 1
          If (IA > LMTNA) Exit! Only B still with unprocessed values
        End If
!
      End Do
    End Do
!
! The Cs become As and Bs
!
    LMTNA = 2 * LMTNA
  End Do
!
 
!
End Subroutine sortrx_no_gvec
!> Threaded sort (actually, ranking) of a real/integer array AA.
!! Public code adapted from:
!! https://github.com/cphyc/Fortran-parallel-sort/blob/master/mod_sort.f90
subroutine sortrx_no_gvec_threaded(NVAL, AA, ord&
)
  integer, intent(in) :: NVAL
  real(DP), intent(in) :: AA(NVAL)
  integer, intent(out) :: ord(NVAL)
  integer :: nthreads
 
  ! FHJ: TODO - nthreads should be different that peinf%nthreads_sort if we
  ! have nested OMPs. However, starting a disabled PARALLEL region to figure out
  ! th number of threads is just too much overhead.
  nthreads = peinf%nthreads_sort
  if (nthreads==1.or.NVAL<nthreads*100) then
    call sortrx_no_gvec(NVAL, AA, ord&
)
   
    return
  endif
 
end subroutine sortrx_no_gvec_threaded
! Merge two parts of A, ordered by order from left to right around middle.
subroutine sortrx_no_gvec_threaded_merge(NVAL, AA, ord, &
left, middle, right)
  integer, intent(in) :: NVAL
  real(DP), intent(in) :: AA(NVAL)
  integer, intent(inout) :: ord(NVAL)
  integer, intent(in) :: left, middle, right
  integer :: leftA, rightA, leftB, rightB
  integer :: iA, iB, i
  integer :: lenA, lenB
  integer :: idxA, idxB
  integer :: orderA(left:middle)
  integer :: orderB(middle+1:right)
 
  ! copy order
  orderA = ord(left:middle)
  orderB = ord(middle+1:right)
  ! more explicit variables
  leftA = left
  rightA = middle
  leftB = middle+1
  rightB = right
  ! initialize iA, iB to their leftmost position
  iA = leftA
  iB = leftB
  i = leftA
  do while ((iA <= rightA) .and. (iB <= rightB))
    idxA = orderA(iA)
    idxB = orderB(iB)
    if ((AA(idxA)>AA(idxB))) then
       ord(i) = idxB
       iB = iB + 1
    else
       ord(i) = idxA
       iA = iA + 1
    end if
    i = i + 1
  end do
  ! either A or B still have elements, append them to the new order
  do while (iA <= rightA)
    ord(i) = orderA(iA)
    iA = iA + 1
    i = i + 1
  end do
  do while (iB <= rightB)
    ord(i) = orderB(iB)
    iB = iB + 1
    i = i + 1
  end do
 
end subroutine sortrx_no_gvec_threaded_merge
end module sort_m
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
