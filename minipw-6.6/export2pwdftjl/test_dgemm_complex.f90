program test

implicit none

integer :: nij, nab, ngm
real(8) :: fact
complex(8), allocatable :: qgm(:,:), aux(:,:)
real(8), allocatable :: deeaux(:,:)


fact = 1.d0
nij = 4
nab = 3
ngm = 5

allocate( qgm(ngm,nij) )
allocate( deeaux(nij,nab) )
allocate( aux(ngm,nab) )

qgm(:,:) = cmplx(1.1d0, 2.3d0, kind=8)
qgm(2,2) = cmplx(7.d0, 0.d0, kind=8)

aux(:,:) = cmplx(5.d0, 3.d0, kind=8)
aux(2,3) = 9.7

deeaux(:,:) = 0.d0

CALL DGEMM('C', 'N', nij, nab, 2*ngm, fact, &
    & qgm, 2*ngm, &
    & aux, 2*ngm, &
    & 0.d0, deeaux, nij )

write(100,*) deeaux

deallocate( qgm )
deallocate( deeaux )
deallocate( aux )



end program