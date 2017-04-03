# GTH pseudopotential

GTH pseudopotential is defined in file `Modules/gth.f90`.

There is one user-defined type in this module:

```fortran
type gth_parameters
   integer  :: itype, lloc, lmax
   real(dp) :: rloc, cc(4)
   integer,  pointer :: lll(:), ipr(:)
   real(dp), pointer :: rrl(:)
end type gth_parameters
```

One global pointer is defined.

```fortran
type (gth_parameters), pointer, dimension(:), save :: gth_p
```

Note that I have removed the attribute `private` of this global
pointer.


# Module `uspp_param`

File: `Modules/uspp.f90`

Despite its name, this module contains parameters for ultrasoft and
norm-converving pseudopotential.

Probably the most important is global array `upf` which is of `pseudo_upf`
type.

```
TYPE (pseudo_upf),  ALLOCATABLE, TARGET :: upf(:)
```

This module also contains the following variables:

```fortran
! number of beta functions per atomic type
INTEGER :: nh(npsx)

! max number of different beta functions per atom
INTEGER :: nhm
```

# Module `uspp`

In file `Modules/uspp.f90`

Variables

```fortran
! Total number of beta functions, with structure factor
INTEGER :: nkb

! correspondence n <-> angular momentum l
INTEGER, ALLOCATABLE :: nhtol(:,:)

! correspondence n <-> combined lm index for (l,m)
INTEGER, ALLOCATABLE :: nhtolm(:,:)

! indes linking  atomic beta's to beta's in the solid
INTEGER, ALLOCATABLE :: indv(:,:)
```



# Type `pseudo_upf`

The type `pseudo_upf` is defined in module `pseudo_types`, in file
`Modules/pseudo_types.f90`.

This type has many fields.
