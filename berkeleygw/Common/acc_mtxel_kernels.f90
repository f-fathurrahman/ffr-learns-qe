
module acc_mtxel_kernels_m
  use ISO_C_BINDING
  use algos_common_m
  use message_m
  use nrtype_m
  use peinfo_m, only: peinf
  use push_pop_m
  implicit none
  private
  integer, private :: mtxel_algo = -1
  public :: allocate_acc_mtxel_sig, deallocate_acc_mtxel_sig
  public :: acc_zero_box, acc_put_into_fftbox, acc_box_multiply, acc_get_from_fftbox
  public :: acc_run_fft
  ! Type definition
  type :: fft_box_type
    complex(DPC), allocatable :: box(:,:,:)
  end type
  type :: fft_vec_type
    complex(DPC), allocatable :: vec(:)
  end type
  type :: acc_mtxel_sig_type
    ! gvec info
    integer, allocatable :: gvec_comp(:,:)
    integer, allocatable :: gvec_isrtx_k(:), gvec_isrtx_q(:), gvec_isrtx_mtxel(:)
    ! eqp aux function (outer WFN)
    complex(DPC), allocatable :: eqp_box(:,:,:), eqp_vec(:)
    ! bands aux function (inner WFN)
    type(fft_vec_type), allocatable :: bands_vec(:)
    ! output (MTXEL)
    type(fft_box_type), allocatable :: mtxel_box(:)
    type(fft_vec_type), allocatable :: mtxel_vec(:)
    ! fft plan
    integer :: eqp_plan
    integer, allocatable :: mtxel_plan(:) ! do we need many plans??
    ! host stuff
    complex(DPC), allocatable :: fftbox_aux(:,:,:)
    real(DP), allocatable :: vec_aux(:)
    ! control size
    integer :: Nfft(3)
    integer :: mtxel_band_block_size
  end type acc_mtxel_sig_type
  type(acc_mtxel_sig_type), public :: acc_mtxel_sig
contains
  subroutine set_mtxel_algo( which_mtxel_algo )
    integer :: which_mtxel_algo
   
    mtxel_algo = which_mtxel_algo
    if ( mtxel_algo /= OPENACC_ALGO .and. mtxel_algo /= OMP_TARGET_ALGO ) then
      call die("MTXEL Algorithm not implemented", only_root_writes = .true.)
    end if
   
  end subroutine set_mtxel_algo
  subroutine acc_run_fft( box, plan, fsign )
    complex(DPC), allocatable :: box(:,:,:)
    integer :: plan, fsign
    integer :: ierr
   
    select case( mtxel_algo )
      case (OPENACC_ALGO)
        call die_algos("OpenACC")
      case(OMP_TARGET_ALGO)
        call die_algos("OpenMP Target")
    end select
   
  end subroutine acc_run_fft
  subroutine acc_zero_box(box, Nfft, queue )
    complex(DPC), allocatable :: box(:,:,:)
    integer :: Nfft(3), queue
    integer :: ix, iy, iz
   
    select case( mtxel_algo )
      case (OPENACC_ALGO)
        call die_algos("OpenACC")
      case (OMP_TARGET_ALGO)
        call die_algos("OpenMP Target")
    end select
   
  end subroutine acc_zero_box
  subroutine acc_box_multiply(box_out, box_in, Nfft, conj, queue )
    complex(DPC), allocatable :: box_out(:,:,:), box_in(:,:,:)
    integer :: Nfft(3), conj, queue
    integer :: ix, iy, iz
   
    select case( mtxel_algo )
      case (OPENACC_ALGO)
        call die_algos("OpenACC")
      case (OMP_TARGET_ALGO)
        call die_algos("OpenMP Target")
    end select
   
  end subroutine acc_box_multiply
  subroutine acc_put_into_fftbox( Ng, vec, g_comp, g_index, box, Nfft, alpha, queue )
    integer :: Ng, Nfft(3)
    complex(DPC), allocatable :: vec(:)
    integer, allocatable :: g_comp(:,:), g_index(:)
    complex(DPC), allocatable :: box(:,:,:)
    real(DP) :: alpha
    integer :: queue
    integer :: iii, bidx(3)
   
    select case( mtxel_algo )
      case (OPENACC_ALGO)
        call die_algos("OpenACC")
      case(OMP_TARGET_ALGO)
        call die_algos("OpenMP Target")
    end select
   
  end subroutine acc_put_into_fftbox
  subroutine acc_get_from_fftbox( Ng, vec, g_comp, g_index, box, Nfft, alpha, queue )
    integer :: Ng, Nfft(3)
    complex(DPC), allocatable :: vec(:)
    integer, allocatable :: g_comp(:,:), g_index(:)
    complex(DPC), allocatable :: box(:,:,:)
    real(DP) :: alpha
    integer :: queue
    integer :: iii, bidx(3)
   
    select case( mtxel_algo )
      case (OPENACC_ALGO)
        call die_algos("OpenACC")
      case (OMP_TARGET_ALGO)
        call die_algos("OpenMP Target")
    end select
   
  end subroutine acc_get_from_fftbox
  subroutine allocate_acc_mtxel_sig ( Nfft, NG_q, NG_k, ng, ncoul, which_mtxel_algo )
    integer, intent(in) :: Nfft(3), NG_q, NG_k, ng, ncoul
    integer, intent(in) :: which_mtxel_algo
    integer :: nbands, n1, n2, n3
    integer :: nloc, ierr
   
    call set_mtxel_algo( which_mtxel_algo )
    nbands = acc_mtxel_sig%mtxel_band_block_size
    n1 = Nfft(1)
    n2 = Nfft(2)
    n3 = Nfft(3)
    acc_mtxel_sig%Nfft = Nfft
    allocate(acc_mtxel_sig%gvec_comp (3,ng))
    allocate(acc_mtxel_sig%gvec_isrtx_k (ng))
    allocate(acc_mtxel_sig%gvec_isrtx_q (ng))
    allocate(acc_mtxel_sig%gvec_isrtx_mtxel (ng))
    acc_mtxel_sig%gvec_comp = 0
    acc_mtxel_sig%gvec_isrtx_k = 0
    acc_mtxel_sig%gvec_isrtx_q = 0
    acc_mtxel_sig%gvec_isrtx_mtxel = 0
    allocate(acc_mtxel_sig%bands_vec (nbands))
    allocate(acc_mtxel_sig%mtxel_vec (nbands))
    allocate(acc_mtxel_sig%mtxel_box (nbands))
    allocate(acc_mtxel_sig%mtxel_plan (nbands))
    allocate(acc_mtxel_sig%eqp_vec (NG_q))
    acc_mtxel_sig%eqp_vec = cmplx(0.0D+00,0.0D+00,kind=DPC)
    allocate(acc_mtxel_sig%eqp_box (n1,n2,n3))
    acc_mtxel_sig%eqp_box = cmplx(0.0D+00,0.0D+00,kind=DPC)
    select case( mtxel_algo )
      case (OPENACC_ALGO)
        call die_algos("OpenACC")
      case (OMP_TARGET_ALGO)
        call die_algos("OpenMP Target")
    end select
    ! is this reverse order?
    do nloc = 1, nbands
      allocate(acc_mtxel_sig%bands_vec(nloc)%vec (NG_k))
      allocate(acc_mtxel_sig%mtxel_vec(nloc)%vec (ncoul))
      allocate(acc_mtxel_sig%mtxel_box(nloc)%box (n1,n2,n3))
      acc_mtxel_sig%bands_vec(nloc)%vec = cmplx(0.0D+00,0.0D+00,kind=DPC)
      acc_mtxel_sig%mtxel_vec(nloc)%vec = cmplx(0.0D+00,0.0D+00,kind=DPC)
      acc_mtxel_sig%mtxel_box(nloc)%box = cmplx(0.0D+00,0.0D+00,kind=DPC)
      select case( mtxel_algo )
        case (OPENACC_ALGO)
        call die_algos("OpenACC")
        case (OMP_TARGET_ALGO)
        call die_algos("OpenMP Target")
      end select
    end do
   
  end subroutine allocate_acc_mtxel_sig
  subroutine deallocate_acc_mtxel_sig ( )
    integer :: nloc, nbands, ierr
   
    nbands = acc_mtxel_sig%mtxel_band_block_size
    do nloc = 1, nbands
      select case( mtxel_algo )
        case (OPENACC_ALGO)
        call die_algos("OpenACC")
        case (OMP_TARGET_ALGO)
        call die_algos("OpenMP Target")
      end select
      if(allocated(acc_mtxel_sig%bands_vec(nloc)%vec))then;deallocate(acc_mtxel_sig%bands_vec(nloc)%vec);endif
      if(allocated(acc_mtxel_sig%mtxel_vec(nloc)%vec))then;deallocate(acc_mtxel_sig%mtxel_vec(nloc)%vec);endif
      if(allocated(acc_mtxel_sig%mtxel_box(nloc)%box))then;deallocate(acc_mtxel_sig%mtxel_box(nloc)%box);endif
    end do
    select case( mtxel_algo )
      case (OPENACC_ALGO)
        call die_algos("OpenACC")
      case (OMP_TARGET_ALGO)
        call die_algos("OpenMP Target")
    end select
    if(allocated(acc_mtxel_sig%gvec_comp))then;deallocate(acc_mtxel_sig%gvec_comp);endif
    if(allocated(acc_mtxel_sig%gvec_isrtx_k))then;deallocate(acc_mtxel_sig%gvec_isrtx_k);endif
    if(allocated(acc_mtxel_sig%gvec_isrtx_q))then;deallocate(acc_mtxel_sig%gvec_isrtx_q);endif
    if(allocated(acc_mtxel_sig%gvec_isrtx_mtxel))then;deallocate(acc_mtxel_sig%gvec_isrtx_mtxel);endif
    if(allocated(acc_mtxel_sig%bands_vec))then;deallocate(acc_mtxel_sig%bands_vec);endif
    if(allocated(acc_mtxel_sig%mtxel_vec))then;deallocate(acc_mtxel_sig%mtxel_vec);endif
    if(allocated(acc_mtxel_sig%mtxel_box))then;deallocate(acc_mtxel_sig%mtxel_box);endif
    if(allocated(acc_mtxel_sig%mtxel_plan))then;deallocate(acc_mtxel_sig%mtxel_plan);endif
    if(allocated(acc_mtxel_sig%eqp_vec))then;deallocate(acc_mtxel_sig%eqp_vec);endif
    if(allocated(acc_mtxel_sig%eqp_box))then;deallocate(acc_mtxel_sig%eqp_box);endif
   
  end subroutine deallocate_acc_mtxel_sig
end module acc_mtxel_kernels_m
