#define PETSC_MAJOR 3
#define PETSC_MINOR_MIN 12
#define PETSC_MINOR_MAX 16


module prec
  use, intrinsic :: IEEE_arithmetic
  use, intrinsic :: ISO_C_binding

#ifdef PETSC
#include <petsc/finclude/petscsys.h>
  use PETScSys
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif


  implicit none
  public

  ! https://stevelionel.com/drfortran/2017/03/27/doctor-fortran-in-it-takes-all-kinds
  integer,     parameter :: pReal      = IEEE_selected_real_kind(15,307)                            !< number with 15 significant digits, up to 1e+-307 (typically 64 bit)
  integer,     parameter :: pI32       = selected_int_kind(9)                                       !< number with at least up to +-1e9 (typically 32 bit)
  integer,     parameter :: pI64       = selected_int_kind(18)                                      !< number with at least up to +-1e18 (typically 64 bit)
#ifdef PETSC
  PetscInt,    public   :: dummy
  integer,     parameter :: pPETSCINT  = kind(dummy)
#endif

  integer, parameter, public :: &
    worldrank = 0, &                                                               !< MPI dummy worldrank
    worldsize = 1                                                                  !< MPI dummy worldsize

  integer,     dimension(3), public, protected :: &
    grid = (/16,16,16/)                                                                                            !< (global) grid
  integer,      public, protected :: &
    grid3 = 16 , &
    grid3Offset = 0                                                                                     !< (local) grid offset in 3rd direction
  real(pReal), dimension(3), public, protected :: &
    geomSize = (/1.0,1.0,1.0/)                                                                                       !< (global) physical size


end module
