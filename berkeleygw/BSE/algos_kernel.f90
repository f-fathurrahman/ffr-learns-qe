! Module to control selection of various accelerated algos in Kernel
! Any build-specific preprocessor directives should be as high-level
! as possible, i.e. no explicit framework-based calls
module algos_kernel_m
  use ISO_C_BINDING
  use algos_common_m
  use message_m
  use nrtype_m, only: DP, DPC
  use peinfo_m, only: peinf
  use push_pop_m
  implicit none
  private
  ! Publicly exposed subroutines
  public :: CPU_ALGO, OPENACC_ALGO, OMP_TARGET_ALGO
  public :: output_algos
  public :: set_algos_to_cpu
  public :: set_algos_to_best_available_gpu
  public :: verify_gpu_settings
  public :: initialize_gpu
  public :: algos_inread
  ! General control variables
  integer(kind(CPU_ALGO)), public :: calc_charge_matrix_algo
  integer(kind(CPU_ALGO)), public :: w_sum_algo
  integer(kind(CPU_ALGO)), public :: g_sum_algo
  ! Variables for OpenACC acceleration
  ! WPH: We put the OpenACC variables here, rather than in their own separate
  ! OpenACC module, because OpenACC is valid CPU code.
  real(DP), target, allocatable, public :: tempb_acc(:,:,:,:), temph_acc(:,:,:)
  real(DP), target, allocatable, public :: tempw_acc(:,:,:,:)
  real(DP), target, allocatable, public :: bsedbody_acc(:,:,:), bsedwing_acc(:,:,:), bsedhead_acc(:,:,:)
  !$ACC DECLARE CREATE(tempb_acc, temph_acc, tempw_acc, bsedbody_acc, bsedwing_acc, bsedhead_acc)
contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! __ ____ _ ___ __ _ !
  ! / / / / /_(_) (_) /_(_)__ _____ !
  ! / / / / __/ / / / __/ / _ \/ ___/ !
  ! / /_/ / /_/ / / / /_/ / __(__ ) !
  ! \____/\__/_/_/_/\__/_/\___/____/ !
  ! !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine output_algos()
   
    if (peinf%inode.eq.0) then
      call output_algos_common()
      write(6,*)
      write(6,'(1x,a)') 'Algorithms used:'
      write(6,'(1x,a,a)') '- calc_charge_matrix: ', trim(get_algo_str(calc_charge_matrix_algo))
      write(6,'(1x,a,a)') '- g_sum             : ', trim(get_algo_str(g_sum_algo))
      write(6,'(1x,a,a)') '- w_sum             : ', trim(get_algo_str(w_sum_algo))
    endif
   
  end subroutine
  subroutine set_algos_to_cpu()
   
    call set_algos_to_cpu_common()
    calc_charge_matrix_algo = CPU_ALGO
    g_sum_algo = CPU_ALGO
    w_sum_algo = CPU_ALGO
   
  end subroutine
  subroutine set_algos_to_best_available_gpu()
   
    call set_algos_to_best_available_gpu_common()
   
  end subroutine
  subroutine verify_gpu_settings()
   
    call verify_gpu_settings_common()
    ! First check if the desired algorithm is compiled in
    if (calc_charge_matrix_algo == OPENACC_ALGO .or. &
        w_sum_algo == OPENACC_ALGO .or. &
        g_sum_algo == OPENACC_ALGO) then
      call die_algos("OpenACC")
    end if
    ! Then check if we`re using the GPU
    if (.not. use_gpu .and. (calc_charge_matrix_algo /= CPU_ALGO .or. &
                             w_sum_algo /= CPU_ALGO .or. &
                             g_sum_algo /= CPU_ALGO)) then
      call die("You have specified one or more GPU-accelerated algorithms, but&
               & the GPU is not enabled.", only_root_writes=.true.)
    end if
    ! Then check for unimplemented algorithms
    if (calc_charge_matrix_algo == OPENACC_ALGO) then
      call die("You have specified the OpenACC version of calc_charge_matrix(),&
               & but this has not been implemented.", only_root_writes=.true.)
    end if
    ! Finally, check for algorithms that don`t play nicely with one another
    if (w_sum_algo == OPENACC_ALGO .and. g_sum_algo /= OPENACC_ALGO) then
      call die("You have specified the OpenACC version of w_sum() and some&
               & other version of g_sum(). The OpenACC version of g_sum() must&
               & be selected as well.", only_root_writes=.true.)
    end if
    if (w_sum_algo /= OPENACC_ALGO .and. g_sum_algo == OPENACC_ALGO) then
      call die("You have specified the OpenACC version of g_sum() and some&
               & other version of w_sum(). The OpenACC version of w_sum() must&
               & be selected as well.", only_root_writes=.true.)
    end if
   
  end subroutine
  subroutine initialize_gpu(algo_type)
    integer(kind(CPU_ALGO)), intent(in) :: algo_type
   
    call initialize_gpu_common(algo_type)
   
  end subroutine initialize_gpu
  subroutine algos_inread(keyword, line, found)
    character(len=*), intent(in) :: keyword, line
    logical, intent(out) :: found
    character*256 :: errmsg
   
    call algos_inread_common(keyword, line, found)
    if (.not. found) then
      found = .true.
      if (trim(keyword) .eq. 'use_gpu') then
        if (trim(line) .eq. '.true.') then
          write (errmsg,'(a)') 'GPU acceleration not compiled into this executable, use_gpu keyword is invalid.'
          call die(errmsg, only_root_writes = .true.)
        else if(trim(line) .eq. '.false.') then
          call set_algos_to_cpu()
        else
          write(errmsg,'(3a)') 'Unexpected parameter "', trim(line), '" for keyword "use_gpu" was found in kernel.inp.'
          call die(errmsg, only_root_writes = .true.)
        end if
      else if(trim(keyword) .eq. 'calc_charge_matrix_algo') then
        calc_charge_matrix_algo = get_algo(trim(line))
      else if(trim(keyword) .eq. 'g_sum_algo') then
        g_sum_algo = get_algo(trim(line))
      else if(trim(keyword) .eq. 'w_sum_algo') then
        w_sum_algo = get_algo(trim(line))
      else
        found = .false.
      end if
    end if
   
  end subroutine algos_inread
end module algos_kernel_m
