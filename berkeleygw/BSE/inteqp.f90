!===============================================================================
!
! Routines:
!
! (1) inteqp Originally by JRD Last Edited: 10/11/2010 (JRD)
!
! Extrapolates Eqp corrections from the coarse grid to the fine grid by
! a wavefunction-based plus linear interpolation scheme that preserves band crossings/character.
! The code reads eqp_co.dat and writes eqp.dat and
! eqp_q.dat. The DFT wavefunctions for the interpolation are read
! from WFN_co, WFN_fi, and WFNq_fi files.
!
! See absorption.inp for options and keywords.
!
!================================================================================
!The following macro puts any point/array in the [-0.5, 0.5) range:
!The following macro puts any point/array in the [0, 1) range:
!Integer division of a/b rounded up*/
!Rounds a up to the smallest multiple of b*/
! disable Fortran OMP pragmas if not -DOMP*/
! note: C standard does not permit $ in identifiers, however this seems acceptable
! as an extension, for all versions of cpp I tried. --DAS
! truncate spaces in string
!#!define TRUNC(s) trim(adjustl(s))
! Oracle compiler has a length limit of 132 characters and won`t support these macros
! No checking for faster performance, if not in debug mode
! Use this instead of the intrinsic 'deallocate' for pointers
! Use this instead of the intrinsic 'deallocate' for arrays
!the TOSTRING macro converts a macro into a string
! deprecated identifiers
! Created Sept 2011 by DAS.
! Define characteristics of various compilers, via compiler symbols (e.g. -DGNU)
! to be used directly from the arch.mk files, and then defining what we need to do
! for that compiler via the symbols for various properties (e.g. NOSIZEOF).
! Ideally, to support a new compiler, one need only change this file, adding a
! new block to define what -DNEWCOMPILER would mean.
! NOTE: of course, Makefile-level issues still need to be handled in common-rules.mk
! very ancient version may require NOSIZEOF
! FHJ: Support for Open64 will be removed shortly in favor of OpenUH
! open64 is very similar to path, it is an open-sourced version of it
! omp_lib.f90 needed to do OpenMP, see common-rules.mk.
! cce 7.4.4 and before support sizeof for intrinsic types, but need NOSIZEOF_TYPE
! cce 8.0.0 and later do not allow sizeof for multidimensional arrays, requiring us
! to turn sizeof off everywhere. Why would Cray do this?
! It is considered a bug in OPEN64 that sizeof will not work in our code.
! on some platforms there is a different return value for sizeof if build is 64-bit
! Intrinsic module for OpenMP. Almost all compilers that support OpenMP provide
! a "omp_lib.mod" module, though the OpenMP standard allow them to only ship a
! "omp_lib.h" Fortran header.
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
program inteqp
  use bse_init_m
  use fullbz_m
  use global_m
  use input_fi_m
  use input_q_m
  use inread_m
  use intwfn_m
  use timing_m, only: timing => bse_timing
  use write_program_header_m
  implicit none
  type (crystal) :: crys
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (eqpinfo) :: eqp
  type (xctinfo) :: xct
  type (flags) :: flag
  type (grid) :: kg_fi,kgq_fi,kg_co,kgq_co
  type (kpoints) :: kp_fi,kpq_fi
  type (int_wavefunction) :: intwfnc
  type (int_wavefunction) :: intwfnv
  integer :: ii,ncount
  integer :: iunit_c,iunit_v
  real(DP) :: vol,omega_plasma
  real(DP) :: tsec(2),tmin(2),tmax(2)
  character*16, allocatable :: routnam(:)
  integer, allocatable :: routsrt(:)
  integer, allocatable :: fi2co_wfn(:,:),indexq_fi(:)
  real(DP), allocatable :: kco(:,:)
  character :: filename*20
  real(DP), allocatable :: dcc(:,:,:,:,:),dvv(:,:,:,:,:)
  real(DP), allocatable :: intp_coefs(:,:)
  call peinfo_init()
!----------------------
! Initialize timer
  call timing%init()
  call timing%start(timing%total)
!---------------------------
! Write header
  call write_program_header('BSE/IntEqp', .false.)
!---------------------------
! Read inteqp.inp
  call logit('Calling inread_inteqp')
  call open_file(8,file='inteqp.inp',form='formatted',status='old')
  call inread(eqp,xct,flag)
  call close_file(8)
! FHJ: Initialize xct%nkpt_co and dimentionality of the problem
  call bse_init(xct,flag)
!--------------------------
! Read wavefunctions on the fine grid
  call logit('Calling input')
  call timing%start(timing%input)
  call input_fi(crys,gvec,kg_fi,kp_fi,syms,eqp,xct,flag,omega_plasma,.false.,intwfnc)
  vol = xct%nktotal*crys%celvol
  if (peinf%inode.eq.0) then
    write(6,*) ' '
    write(6,*) 'More Job Parameters: '
    write(6,'(a,f32.14,a)') ' Crystal volume = ',vol,' a.u.'
    write(6,*) 'Number of valence bands = ',xct%nvb_fi
    write(6,*) 'Number of cond. bands   = ',xct%ncb_fi
    write(6,*) 'Number of spins   = ',xct%nspin
    write(6,*) ' '
  endif
  call timing%stop(timing%input)
  allocate(indexq_fi (xct%nkpt_fi))
  allocate(xct%indexq_fi (xct%nkpt_fi))
! if (flag%vm.ne.1.or.flag%dtm.ne.1) then ! both are always 0 in this code --DAS
  call timing%start(timing%input_q)
  call logit('Calling input_q')
  call input_q(kp_fi,crys,gvec,kg_fi,kgq_fi,kpq_fi,syms,xct,indexq_fi,eqp,flag,intwfnv)
  call timing%stop(timing%input_q)
!------------------------------
! Calculate the transformation matrices from coarse grid wavefunctions
! FHJ: These are the final transformation coefs that will be used to interpolate
! the kernel. However, we might use an unrestricted version of dvv/dcc to
! interpolate eqp if xct%unrestricted_transf==.true..
  allocate(dvv (xct%nvb_fi,xct%n1b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
  allocate(dcc (xct%ncb_fi,xct%n2b_co,xct%nspin,xct%nkpt_fi,xct%npts_intp_kernel))
  allocate(kco (3,xct%nkpt_co))
  allocate(fi2co_wfn (xct%npts_intp_kernel,xct%nkpt_fi))
  allocate(intp_coefs (xct%npts_intp_kernel, xct%nkpt_fi))
  call logit('Calling intwfn')
  call timing%start(timing%intwfn)
  call intwfn(kp_fi,crys,syms,xct,flag,gvec,kg_fi,kgq_fi,kg_co,kgq_co,dcc,dvv,&
    kco,fi2co_wfn,indexq_fi,eqp,intwfnv,intwfnc,intp_coefs)
  call timing%stop(timing%intwfn)
  call kp_fi%free()
  call kpq_fi%free()
  call dealloc_grid(kg_fi)
  call dealloc_grid(kgq_fi)
  if(associated(xct%ifmax))then;deallocate(xct%ifmax);nullify(xct%ifmax);endif
  if(associated(xct%ifmaxq))then;deallocate(xct%ifmaxq);nullify(xct%ifmaxq);endif
  if(associated(intwfnc%cgk))then;deallocate(intwfnc%cgk);nullify(intwfnc%cgk);endif
  if(associated(intwfnv%cgk))then;deallocate(intwfnv%cgk);nullify(intwfnv%cgk);endif
  if(associated(intwfnc%isort))then;deallocate(intwfnc%isort);nullify(intwfnc%isort);endif
  if(associated(intwfnv%isort))then;deallocate(intwfnv%isort);nullify(intwfnv%isort);endif
  if(associated(eqp%ecqp))then;deallocate(eqp%ecqp);nullify(eqp%ecqp);endif
  if(associated(eqp%evqp))then;deallocate(eqp%evqp);nullify(eqp%evqp);endif
!--------------------------------
! Time accounting
  call timing%stop(timing%total)
  call timing%print()
  call write_memory_usage()
end program inteqp
