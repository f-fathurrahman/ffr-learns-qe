!===============================================================================
!
! MODULE:
!
! tile_m Last Modified: Dec/2012 (FHJ)
!
!> Tiling and tessellation routines for the R^2/R^3. Currently, only Delaunay
!! triangulation is available, but Voronoi tessellation should be added soon.
!
! DESCRIPTION:
!
! This module is an abstraction layer that provides tessellation methods in
! arbitrary dimensions. The "API" declared here (in the `contains` section)
! should be generic, and the methods may be implemented by in external
! libraries. In this case, please write appropriate wrappers so that the
! routines declared here do not depend on the technical details of the
! external back end.
!
! We currently use Qhull as the back end for the Delaunay tessellation.
!
!===============================================================================

module tile_m
  use global_m
  implicit none
  public :: &
    init_delaunay, &
    find_delaunay_simplex, &
    get_num_simplices, &
    get_simplices, &
    get_neighbors, &
    init_second_neighbors, &
    get_max_num_second_neighbors, &
    find_delaunay_simplex_with_second_neighbors, &
    free_delaunay
  private
  interface
    !---------------------------------------------------------------------------
    ! FHJ: these tessellation functions are implemented in qhull/libtile_qhull.a
    integer function qhull_init_delaunay(points, num_points, dimensions, inode)
      use global_m
      implicit none
      real(DP), intent(in) :: points
      integer, intent(in) :: num_points
      integer, intent(in) :: dimensions
      integer, intent(in) :: inode
    end function qhull_init_delaunay
    !> Finds the Delaunay triangle/tetrahedron that encloses `point`.
    !!
    !! @param point [in] array (ndim).
    !! @param indices [out] array (ndim+1) of the (Fortran) indices of the vertices.
    !! @param coefs [out] coefficients of point in barycentric coordinates.
    integer function qhull_find_delaunay_simplex(point, indices, coefs)
      use global_m
      implicit none
      real(DP), intent(in) :: point
      integer, intent(out) :: indices
      real(DP), intent(out) :: coefs
    end function qhull_find_delaunay_simplex
    integer function qhull_get_num_simplices(num_simplices)
      implicit none
      integer, intent(out) :: num_simplices
    end function qhull_get_num_simplices
    integer function qhull_get_simplices(indices)
      implicit none
      integer, intent(out) :: indices
    end function qhull_get_simplices
    integer function qhull_get_neighbors(neighbors)
      implicit none
      integer, intent(out) :: neighbors
    end function qhull_get_neighbors
    integer function qhull_init_second_neighbors()
      implicit none
    end function qhull_init_second_neighbors
    integer function qhull_get_max_num_second_neighbors(n_nei)
      implicit none
      integer, intent(out) :: n_nei
    end function qhull_get_max_num_second_neighbors
    !> Computes second-nearest neighbors and barycentric coefficients for point `point`.
    !!
    !! The arrrays `indices` and `coefs` must be ndim+1 + `max_num_second_neighbors`
    !! big. The first ndim+1 points will correspond to the first neighbors.
    !! `nnei` is set to the total number of neighbors found, with
    !! `nnei` <= ndim + 1 + max_num_second_neighbors.
    integer function qhull_find_delaunay_simplex_with_second_neighbors(point, indices, coefs, nnei)
      use global_m
      implicit none
      real(DP), intent(in) :: point
      integer, intent(out) :: indices
      real(DP), intent(out) :: coefs
      integer, intent(out) :: nnei
    end function qhull_find_delaunay_simplex_with_second_neighbors
    integer function qhull_free_delaunay()
      implicit none
    end function qhull_free_delaunay
  end interface
contains
  !> Calculates the Delaunay triangulation/tetrahedralization of a set of points.
  !!
  !! For performance issues, only one Delaunay triangulation can be active at
  !! a given time.
  !!
  !! \param points [inout] coarse points (dim, npts)
  !! \param num_points [in] number of points.
  !! \param dim [in] number of dimensions
  !! \return 0 for success
  integer function init_delaunay(points, num_points, dimensions)
    real(DP), intent(in) :: points(:,:)
    integer, intent(in) :: num_points
    integer, intent(in) :: dimensions
   
    init_delaunay = qhull_init_delaunay(points(1,1), num_points, dimensions, peinf%inode)
   
    return
  end function init_delaunay
  !> Finds the Delaunay triangle/tetrahedron that encloses point.
  !!
  !! \param point [in] real array (dim) containing the coordinates.
  !! \param indices [out] integer array (dim+1) with the indices of the vertices.
  !! \param coefs [out] coefficients of point in barycentric coordinates.
  integer function find_delaunay_simplex(point, indices, coefs)
    real(DP), intent(in) :: point(:)
    integer, intent(out) :: indices(:)
    real(DP), intent(out) :: coefs(:)
   
    find_delaunay_simplex = qhull_find_delaunay_simplex(point(1), indices(1), coefs(1))
   
    return
  end function find_delaunay_simplex
  !> Returns the total number of simplices obtained from the Delaunay
  !! triangulation.
  !!
  !! \param num_simplices [out] number of simplices.
  integer function get_num_simplices(num_simplices)
    integer, intent(out) :: num_simplices
   
    get_num_simplices = qhull_get_num_simplices(num_simplices)
   
    return
  end function get_num_simplices
  !> Returns the indices of the points that define each simplex.
  !!
  !! \param simplices [out] (ndims+1, num_simplices) Vertex ivert of
  !! simplex isimp corresponds to the point simplices(ivert,isimp) of
  !! the original list of points.
  integer function get_simplices(simplices)
    integer, intent(out) :: simplices(:,:)
   
    get_simplices = qhull_get_simplices(simplices(1,1))
   
    return
  end function get_simplices
  !> Returns the indices of the neighbors for each simplex.
  !!
  !! \param neighbors [out] (ndims+1, num_simplices) The neighbor ivert
  !! of simplex isimp is neighbors(ivert,isimp), and corresponds to the
  !! simplex reached by following the ridge opposite to vertex ivert.
  integer function get_neighbors(neighbors)
    integer, intent(out) :: neighbors(:,:)
   
    get_neighbors = qhull_get_neighbors(neighbors(1,1))
   
    return
  end function get_neighbors
  !> Initializes extra buffers required to compute second nearest neighbors
  !! You don`t need to do any extra procedure to free this functin, just
  !! call `free_delaunay` at the end as usual.
  integer function init_second_neighbors()
   
    init_second_neighbors = qhull_init_second_neighbors()
   
    return
  end function init_second_neighbors
  !> Get maximum number of second nearest neighbors
  integer function get_max_num_second_neighbors(n_nei)
    integer, intent(out) :: n_nei
   
    get_max_num_second_neighbors = qhull_get_max_num_second_neighbors(n_nei)
   
    return
  end function get_max_num_second_neighbors
  !> Finds the Delaunay simplex that encloses a point, plus second neighbors.
  !!
  !! \param point [in] real array (dim) containing the coordinates.
  !! \param indices [out] integer array (dim+1) with the indices of the vertices.
  !! The first ndim+1 points are the first neighbors, and remaining are
  !! second neighbors. The size of the array indices MUST BE that given
  !! by get_max_num_second_neighbors().
  !! \param coefs [out] coefficients of point in barycentric coordinates.
  !! The first ndim+1 points are the first neighbors and will have the
  !! regular coefficients for barycentric interpolation. The remaining
  !! coefficients are zero.
  !! \param nnei [out] total number of neighbors (first+second) around this
  !! particular input point. The value of the array indices(:) will
  !! be zero for points outside nnei.
  integer function find_delaunay_simplex_with_second_neighbors(point, indices, coefs, nnei)
    real(DP), intent(in) :: point(:)
    integer, intent(out) :: indices(:)
    real(DP), intent(out) :: coefs(:)
    integer, intent(out) :: nnei
   
    find_delaunay_simplex_with_second_neighbors = &
      qhull_find_delaunay_simplex_with_second_neighbors(point(1), indices(1), coefs(1), nnei)
   
    return
  end function find_delaunay_simplex_with_second_neighbors
  !> Frees buffers associated to the Delaunay triangulation.
  integer function free_delaunay()
   
    free_delaunay = qhull_free_delaunay()
   
    return
  end function free_delaunay
end module tile_m
