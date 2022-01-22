!--------------------------------------------------------------------------------------------------
!> @author Arko Jyoti Bhattacharjee, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Grid solver for mechanics: FEM
!--------------------------------------------------------------------------------------------------
module grid_mechanical_FEM
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
  DM   :: mechanical_grid
  SNES :: SNES_mechanical
  Vec  :: solution_current, solution_lastInc, solution_rate

!--------------------------------------------------------------------------------------------------
! common pointwise data
  real(pReal), dimension(:,:,:,:,:), allocatable :: F, P_current, F_lastInc
  real(pReal) :: detJ
  real(pReal), dimension(3)   :: delta
  real(pReal), dimension(3,8) :: BMat
  real(pReal), dimension(8,8) :: HGMat

!--------------------------------------------------------------------------------------------------
! stress, stiffness and compliance average etc.
  real(pReal), dimension(3,3,3,3) :: &
    S = 0.0_pReal                                                                                   !< current compliance (filled up with zeros)

  public :: &
    grid_mechanical_FEM_init

contains

!--------------------------------------------------------------------------------------------------
!> @brief allocates all necessary fields and fills them with data, potentially from restart info
!--------------------------------------------------------------------------------------------------
subroutine grid_mechanical_FEM_init

  real(pReal), parameter :: HGCoeff = 0.0e-2_pReal
  real(pReal), parameter, dimension(4,8) :: &
    HGcomp = reshape([ 1.0_pReal, 1.0_pReal, 1.0_pReal,-1.0_pReal, &
                       1.0_pReal,-1.0_pReal,-1.0_pReal, 1.0_pReal, &
                      -1.0_pReal, 1.0_pReal,-1.0_pReal, 1.0_pReal, &
                      -1.0_pReal,-1.0_pReal, 1.0_pReal,-1.0_pReal, &
                      -1.0_pReal,-1.0_pReal, 1.0_pReal, 1.0_pReal, &
                      -1.0_pReal, 1.0_pReal,-1.0_pReal,-1.0_pReal, &
                       1.0_pReal,-1.0_pReal,-1.0_pReal,-1.0_pReal, &
                       1.0_pReal, 1.0_pReal, 1.0_pReal, 1.0_pReal], [4,8])
  real(pReal), dimension(3,3,3,3) :: devNull
  PetscErrorCode :: err_PETSc
  PetscScalar, pointer, dimension(:,:,:,:) :: &
    u_current,u_lastInc
  PetscInt, dimension(0:worldsize-1) :: localK


!--------------------------------------------------------------------------------------------------
! set default and user defined options for PETSc
  call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
                                '-mechanical_snes_type newtonls -mechanical_ksp_type fgmres &
                                &-mechanical_ksp_max_it 25', &
                                err_PETSc)
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! allocate global fields
  allocate(F (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
  allocate(P_current (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)
  allocate(F_lastInc (3,3,grid(1),grid(2),grid3),source = 0.0_pReal)

!--------------------------------------------------------------------------------------------------
! initialize solver specific parts of PETSc
  call SNESCreate(PETSC_COMM_WORLD,SNES_mechanical,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetOptionsPrefix(SNES_mechanical,'mechanical_',err_PETSc)
  CHKERRQ(err_PETSc)
  localK            = 0
  localK(worldrank) = grid3
  call DMDACreate3d(PETSC_COMM_WORLD, &
         DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &
         DMDA_STENCIL_BOX, &
         grid(1),grid(2),grid(3), &                    ! global grid
         1, 1, worldsize, &
         3, 1, &                                                                ! #dof (u, vector), ghost boundary width (domain overlap)
         [grid(1)],[grid(2)],localK, &                                ! local grid
         mechanical_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMsetFromOptions(mechanical_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMsetUp(mechanical_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDASetUniformCoordinates(mechanical_grid,0.0_pReal,geomSize(1),0.0_pReal,geomSize(2),0.0_pReal,geomSize(3),err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(mechanical_grid,solution_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(mechanical_grid,solution_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMCreateGlobalVector(mechanical_grid,solution_rate   ,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSNESSetFunctionLocal(mechanical_grid,formResidual,PETSC_NULL_SNES,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMSNESSetJacobianLocal(mechanical_grid,formJacobian,PETSC_NULL_SNES,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetMaxLinearSolveFailures(SNES_mechanical, huge(1), err_PETSc)                 ! ignore linear solve failures
  CHKERRQ(err_PETSc)
  call SNESSetDM(SNES_mechanical,mechanical_grid,err_PETSc)
  CHKERRQ(err_PETSc)
  call SNESSetFromOptions(SNES_mechanical,err_PETSc)                                                ! pull it all together with additional cli arguments
  CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! init fields
  call VecSet(solution_current,0.0_pReal,err_PETSc);CHKERRQ(err_PETSc)
  call VecSet(solution_lastInc,0.0_pReal,err_PETSc);CHKERRQ(err_PETSc)
  call VecSet(solution_rate   ,0.0_pReal,err_PETSc);CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(mechanical_grid,solution_current,u_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(mechanical_grid,solution_lastInc,u_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

  delta = geomSize/real(grid,pReal)                                                                 ! grid spacing
  detJ = product(delta)                                                                             ! cell volume

  BMat = reshape(real([-1.0_pReal/delta(1),-1.0_pReal/delta(2),-1.0_pReal/delta(3), &
                        1.0_pReal/delta(1),-1.0_pReal/delta(2),-1.0_pReal/delta(3), &
                       -1.0_pReal/delta(1), 1.0_pReal/delta(2),-1.0_pReal/delta(3), &
                        1.0_pReal/delta(1), 1.0_pReal/delta(2),-1.0_pReal/delta(3), &
                       -1.0_pReal/delta(1),-1.0_pReal/delta(2), 1.0_pReal/delta(3), &
                        1.0_pReal/delta(1),-1.0_pReal/delta(2), 1.0_pReal/delta(3), &
                       -1.0_pReal/delta(1), 1.0_pReal/delta(2), 1.0_pReal/delta(3), &
                        1.0_pReal/delta(1), 1.0_pReal/delta(2), 1.0_pReal/delta(3)],pReal), [3,8])/4.0_pReal ! shape function derivative matrix

  HGMat = matmul(transpose(HGcomp),HGcomp) &
        * HGCoeff*(delta(1)*delta(2) + delta(2)*delta(3) + delta(3)*delta(1))/16.0_pReal            ! hourglass stabilization matrix

!--------------------------------------------------------------------------------------------------
! init fields

  call DMDAVecRestoreArrayF90(mechanical_grid,solution_current,u_current,err_PETSc)
  CHKERRQ(err_PETSc)
  call DMDAVecRestoreArrayF90(mechanical_grid,solution_lastInc,u_lastInc,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine grid_mechanical_FEM_init


!--------------------------------------------------------------------------------------------------
!> @brief forms the residual vector
!--------------------------------------------------------------------------------------------------
subroutine formResidual(da_local,x_local, &
                        f_local,dummy,err_PETSc)

  DM                   :: da_local
  Vec                  :: x_local, f_local
  PetscScalar, pointer,dimension(:,:,:,:) :: x_scal, r
  PetscScalar, dimension(8,3) :: x_elem,  f_elem
  PetscInt :: &
    PETScIter
  PetscObject :: dummy
  PetscErrorCode :: err_PETSc

  call DMDAVecGetArrayF90(da_local,x_local,x_scal,err_PETSc);CHKERRQ(err_PETSc)
  call DMDAVecRestoreArrayF90(da_local,x_local,x_scal,err_PETSc);CHKERRQ(err_PETSc)

!--------------------------------------------------------------------------------------------------
! constructing residual
  call VecSet(f_local,0.0_pReal,err_PETSc);CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(da_local,f_local,r,err_PETSc);CHKERRQ(err_PETSc)
  call DMDAVecGetArrayF90(da_local,x_local,x_scal,err_PETSc);CHKERRQ(err_PETSc)
  r = 0.0
  call DMDAVecRestoreArrayF90(da_local,x_local,x_scal,err_PETSc);CHKERRQ(err_PETSc)
  call DMDAVecRestoreArrayF90(da_local,f_local,r,err_PETSc);CHKERRQ(err_PETSc)

end subroutine formResidual


!--------------------------------------------------------------------------------------------------
!> @brief forms the FEM stiffness matrix
!--------------------------------------------------------------------------------------------------
subroutine formJacobian(da_local,x_local,Jac_pre,Jac,dummy,err_PETSc)

  DM                                   :: da_local
  Vec                                  :: x_local
  Mat                                  :: Jac_pre, Jac
  MatStencil,dimension(4,24)           :: row, col
  PetscScalar,pointer,dimension(:,:,:,:) :: x_scal
  PetscScalar,dimension(24,24)         :: K_ele
  PetscScalar,dimension(9,24)          :: BMatFull
  PetscInt                             :: i, ii, j, jj, k, kk, ctr, ce
  PetscObject                          :: dummy
  PetscErrorCode                       :: err_PETSc

  BMatFull = 0.0
  BMatFull(1:3,1 :8 ) = BMat
  BMatFull(4:6,9 :16) = BMat
  BMatFull(7:9,17:24) = BMat
  call MatSetOption(Jac,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatSetOption(Jac,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatZeroEntries(Jac,err_PETSc); CHKERRQ(err_PETSc)
  ctr = 1
  col(MatStencil_i,ctr   ) = 1
  col(MatStencil_j,ctr   ) = 1
  col(MatStencil_k,ctr   ) = 1
  col(MatStencil_c,ctr   ) = 0
  col(MatStencil_i,ctr+8 ) = 1
  col(MatStencil_j,ctr+8 ) = 1
  col(MatStencil_k,ctr+8 ) = 1
  col(MatStencil_c,ctr+8 ) = 1
  col(MatStencil_i,ctr+16) = 1
  col(MatStencil_j,ctr+16) = 1
  col(MatStencil_k,ctr+16) = 1
  col(MatStencil_c,ctr+16) = 2
  row = col
  K_ele = 1.0
  call MatSetValuesStencil(Jac,24,row,24,col,K_ele,ADD_VALUES,err_PETSc)
  CHKERRQ(err_PETSc)
  call MatAssemblyBegin(Jac,MAT_FINAL_ASSEMBLY,err_PETSc); CHKERRQ(err_PETSc)
  call MatAssemblyEnd(Jac,MAT_FINAL_ASSEMBLY,err_PETSc); CHKERRQ(err_PETSc)
  call MatAssemblyBegin(Jac_pre,MAT_FINAL_ASSEMBLY,err_PETSc); CHKERRQ(err_PETSc)
  call MatAssemblyEnd(Jac_pre,MAT_FINAL_ASSEMBLY,err_PETSc); CHKERRQ(err_PETSc)

end subroutine formJacobian

end module grid_mechanical_FEM
