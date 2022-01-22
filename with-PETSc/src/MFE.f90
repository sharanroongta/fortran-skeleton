program minimum_failing_example
  
  use, intrinsic :: IEEE_arithmetic
#include <petsc/finclude/petsc.h>
  use petsc
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>14) && !defined(PETSC_HAVE_MPI_F90MODULE_VISIBILITY)
  use MPI_f08
#endif

  implicit none
  integer, parameter :: pReal = IEEE_selected_real_kind(15,307)                                     !< number with 15 significant digits, up to 1e+-307 (typically 64 bit)
  
  integer(MPI_INTEGER_KIND) :: &
    worldrank = 0_MPI_INTEGER_KIND, &                                                               !< MPI worldrank (/=0 for MPI simulations only)
    worldsize = 1_MPI_INTEGER_KIND                                                                  !< MPI worldsize (/=1 for MPI simulations only)
  
  DM :: mechanical_grid, damage_grid
  SNES :: SNES_mechanical, SNES_damage
  PetscErrorCode :: err_PETSc
  
  integer, dimension(3), parameter :: grid = 4
  integer, parameter :: grid3= 4, grid3offset = 0 

  real(pReal), parameter, dimension(3) :: geomSize = 1.0_pReal
  

  call PetscInitialize(PETSC_NULL_CHARACTER,err_PETSc)
  CHKERRA(err_PETSc)
 
! if (.false.) then
  ! -----------------------------------------------------------------------------------------------
  ! BEGIN BLOCK 1 (works by itself)
  call SNESCreate(PETSC_COMM_WORLD,SNES_mechanical,err_PETSc)
  CHKERRA(err_PETSc)
  call DMDACreate3d(PETSC_COMM_WORLD, &
         DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &
         DMDA_STENCIL_BOX, &
         grid(1),grid(2),grid(3), &
         1, 1, worldsize, &
         3, 1, &
         [grid(1)],[grid(2)],[grid3], &
         mechanical_grid,err_PETSc)
  CHKERRA(err_PETSc)
  call DMsetUp(mechanical_grid,err_PETSc)
  CHKERRA(err_PETSc)
  call DMDASetUniformCoordinates(mechanical_grid,0.0_pReal,geomSize(1),0.0_pReal,geomSize(2),0.0_pReal,geomSize(3),err_PETSc)
  CHKERRA(err_PETSc)
  call DMSNESSetFunctionLocal(mechanical_grid,residual_mech,PETSC_NULL_SNES,err_PETSc)
  CHKERRA(err_PETSc)
  call DMSNESSetJacobianLocal(mechanical_grid,jacobian_mech,PETSC_NULL_SNES,err_PETSc)
  CHKERRA(err_PETSc)
  ! END BLOCK 1
  ! -----------------------------------------------------------------------------------------------
! end if 

! if (.false.) then
  ! -----------------------------------------------------------------------------------------------
  ! BEGIN BLOCK 2 (works by itself)
  call SNESCreate(PETSC_COMM_WORLD,SNES_damage,err_PETSc)
  CHKERRA(err_PETSc)
  call DMDACreate3D(PETSC_COMM_WORLD, &
         DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, &                                    ! cut off stencil at boundary
         DMDA_STENCIL_BOX, &                                                                        ! Moore (26) neighborhood around central point
         grid(1),grid(2),grid(3), &                                                                 ! global grid
         1, 1, worldsize, &
         1, 0, &                                                                                    ! #dof (phi, scalar), ghost boundary width (domain overlap)
         [grid(1)],[grid(2)],[grid3], &
         damage_grid,err_PETSc)                                                                     ! handle, error
  CHKERRA(err_PETSc)
  call DMsetUp(damage_grid,err_PETSc)
  CHKERRA(err_PETSc)
  call DMDASNESSetFunctionLocal(damage_grid,INSERT_VALUES,residual_damage,PETSC_NULL_SNES,err_PETSc)   ! residual vector of same shape as solution vector
  CHKERRA(err_PETSc)
  ! END BLOCK 2
  ! -----------------------------------------------------------------------------------------------
! end if

contains

subroutine residual_mech(da_local,x_local, f_local,dummy,err_PETSc)
#include <petsc/finclude/petsc.h>

  DM  :: da_local
  Vec :: x_local, f_local
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc
  
  print*, 'residual_mech'; flush(6)

end subroutine residual_mech

subroutine residual_damage(in,x_scal,r,dummy,err_PETSc)
#include <petsc/finclude/petsc.h>

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
  
  print*, 'residual_damage'; flush(6)

end subroutine residual_damage

subroutine jacobian_mech(da_local,x_local,Jac_pre,Jac,dummy,err_PETSc)
#include <petsc/finclude/petsc.h>

  DM             :: da_local
  Vec            :: x_local
  Mat            :: Jac_pre, Jac
  PetscObject    :: dummy
  PetscErrorCode :: err_PETSc
  
  print*, 'jacobian_mech'; flush(6)

end subroutine jacobian_mech

end program minimum_failing_example
