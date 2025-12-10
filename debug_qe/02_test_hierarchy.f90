program test_mp
  use io_global, only: ionode
  use mp_global, only: mp_startup, mp_global_end
  !USE environment, only : environment_start, environment_end
  implicit none

  !call mp_startup(start_images=.true., images_only=.true.)
  call mp_startup()

  call print_mpi_hierarchy()

  call mp_global_end()
  stop

end program


!===============================================================
!  Minimal MPI hierarchy printer for Quantum ESPRESSO
!===============================================================
SUBROUTINE print_mpi_hierarchy()
    USE kinds, only: DP
    USE mp, only: mp_rank, mp_barrier
    USE mp_world, only: world_comm
    use mp_pools, only: inter_pool_comm, intra_pool_comm
    USE mp_bands, only: intra_bgrp_comm
    USE io_global, only: stdout

    IMPLICIT NONE
    INTEGER :: me, np
    integer :: ierr

    CALL mp_barrier(world_comm)

    !=======================
    ! world_comm
    !=======================
    CALL mpi_comm_size(world_comm, np, ierr)
    me = mp_rank(world_comm)
    write(*,*) 'me = ', me
    IF (me == 0) THEN
        WRITE(stdout,'(/,a,i6)') 'world_comm size: ', np
        WRITE(stdout,'(a)') '-----------------------------------------'
    END IF
    CALL mp_barrier(world_comm)

    !=======================
    ! inter_pool_comm
    !=======================
    !CALL mpi_comm_size(inter_pool_comm, np)
    !!me = mp_rank(inter_pool_comm)
    !IF (me == 0) THEN
    !    WRITE(stdout,'(a,i6)') 'inter_pool_comm size: ', np
    !END IF
    !CALL mp_barrier(world_comm)

    !=======================
    ! intra_pool_comm
    !=======================
    !CALL mpi_comm_size(intra_pool_comm, np)
    !me = mp_rank(intra_pool_comm)
    !IF (me == 0) THEN
    !    WRITE(stdout,'(a,i6)') 'intra_pool_comm size: ', np
    !END IF
    !CALL mp_barrier(world_comm)

    !=======================
    ! intra_bgrp_comm
    !=======================
    !CALL mpi_comm_size(intra_bgrp_comm, np)
    !me = mp_rank(intra_bgrp_comm)
    !IF (me == 0) THEN
    !    WRITE(stdout,'(a,i6)') 'intra_bgrp_comm size: ', np
    !END IF
    !CALL mp_barrier(world_comm)

    !=======================
    ! intra_pgrp_comm
    !=======================
    !CALL mp_comm_size(intra_pgrp_comm, np)
    !me = CALL mp_rank(intra_pgrp_comm)
    !IF (me == 0) THEN
    !    WRITE(stdout,'(a,i6)') 'intra_pgrp_comm size: ', np
    !END IF
    !CALL mp_barrier(world_comm)

    !=======================
    ! intra_tg_comm
    !=======================
    !CALL mpi_comm_size(intra_tg_comm, np)
    !me = mp_rank(intra_tg_comm)
    !IF (me == 0) THEN
    !    WRITE(stdout,'(a,i6)') 'intra_tg_comm size: ', np
    !END IF
    !CALL mp_barrier(world_comm)

    !=======================
    ! Finished
    !=======================
    !me = mp_rank(world_comm)
    !IF (me == 0) WRITE(stdout,'(a,/)') '=== Finished printing MPI hierarchy ==='

    CALL mp_barrier(world_comm)
END SUBROUTINE print_mpi_hierarchy
