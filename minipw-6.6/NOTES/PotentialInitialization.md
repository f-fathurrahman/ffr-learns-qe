Subroutine: potinit

potinit is called in init_run

In the case of fresh calculation, electron density is initialized
from superposition of free atomic electron density.
Relevant subroutine: atomic_rho_g

```fortran
CALL atomic_rho_g( rho%of_g, nspin )
```
