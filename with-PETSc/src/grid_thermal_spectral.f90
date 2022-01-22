!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Shaokang Zhang, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Spectral solver for thermal conduction
!--------------------------------------------------------------------------------------------------
module grid_thermal_spectral
#include <petsc/finclude/petscsnes.h>
#include <petsc/finclude/petscdmda.h>
  use PETScDMDA
  use PETScSNES
!#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
!  use MPI_f08
!#endif

  use prec

  implicit none
  private

!--------------------------------------------------------------------------------------------------
! PETSc data
  SNES :: SNES_thermal
  Vec :: solution_vec
  real(pReal), dimension(:,:,:), allocatable :: &
    T_current                                                                                      !< field of current temperature

  public :: &
    grid_thermal_spectral_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all neccessary fields and fills them with data
! ToDo: Restart not implemented
!--------------------------------------------------------------------------------------------------
subroutine grid_thermal_spectral_init()

  PetscInt, dimension(0:worldsize-1) :: localK
  integer :: i, j, k, ce
  DM :: thermal_grid
  PetscScalar, dimension(:,:,:), pointer :: T_PETSc
  PetscErrorCode :: err_PETSc

  print'(/,1x,a)', '<<<+-  grid_thermal_spectral init  -+>>>'

!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
 call PetscOptionsInsertString(PETSC_NULL_OPTIONS,'-thermal_snes_type newtonls -thermal_snes_mf &
                               &-thermal_snes_ksp_ew -thermal_ksp_type fgmres',err_PETSc)
 CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! init fields
  allocate(T_current(grid(1),grid(2),grid3), source=0.0_pReal)

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,SNES_thermal,err_PETSc); CHKERRQ(err_PETSc)
  call SNESSetOptionsPrefix(SNES_thermal,'thermal_',err_PETSc);CHKERRQ(err_PETSc)
  localK            = 0
  localK(worldrank) = grid3
  call DMDACreate3D(PETSC_COMM_WORLD, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                    ! cut off stencil at boundary
         DMDA_STENCIL_BOX, &                                                                        ! Moore (26) neighborhood around central point
         grid(1),grid(2),grid(3), &                    ! global grid
         1, 1, worldsize, &
         1, 0, &                                                                ! #dof (T, scalar), ghost boundary width (domain overlap)
         [grid(1)],[grid(2)],localK, &                                ! local grid
         thermal_grid,err_PETSc)                                                                    ! handle, error
  CHKERRQ(err_PETSc)
  call DMsetFromOptions(thermal_grid,err_PETSc); CHKERRQ(err_PETSc)
  call DMsetUp(thermal_grid,err_PETSc); CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(thermal_grid,solution_vec,err_PETSc)                                    ! global solution vector (grid x 1, i.e. every def grad tensor)
  CHKERRQ(err_PETSc)
  call DMDASNESSetFunctionLocal(thermal_grid,INSERT_VALUES,formResidual,PETSC_NULL_SNES,err_PETSc)  ! residual vector of same shape as solution vector
  CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_thermal,thermal_grid,err_PETSc); CHKERRQ(err_PETSc)
  call SNESSetFromOptions(SNES_thermal,err_PETSc); CHKERRQ(err_PETSc)                               ! pull it all together with additional CLI arguments
  call DMDAVecGetArrayF90(thermal_grid,solution_vec,T_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)
  T_PETSc = T_current
  call DMDAVecRestoreArrayF90(thermal_grid,solution_vec,T_PETSc,err_PETSc)
  CHKERRQ(err_PETSc)


end subroutine grid_thermal_spectral_init


!--------------------------------------------------------------------------------------------------
!> @brief forms the spectral thermal residual vector
!--------------------------------------------------------------------------------------------------
subroutine formResidual(in,x_scal,r,dummy,err_PETSc)

  DMDALocalInfo, dimension(DMDA_LOCAL_INFO_SIZE) :: &
    in

  PetscScalar, dimension( &
    XG_RANGE,YG_RANGE,ZG_RANGE), intent(in) :: &
    x_scal
  PetscScalar, dimension( &
    X_RANGE,Y_RANGE,Z_RANGE), intent(out) :: &
    r
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc

  T_current = x_scal

!--------------------------------------------------------------------------------------------------
! constructing residual
  r = T_current
  err_PETSc = 0

end subroutine formResidual


end module grid_thermal_spectral
