program DAMASK

#include <petsc/finclude/petscsys.h>
  use PETScSys
!#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
!  use MPI_f08
!#endif

  use grid_mechanical_FEM
  use grid_thermal_spectral

  call grid_mechanical_FEM_init()
  call grid_thermal_spectral_init()

end program
