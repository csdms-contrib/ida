/* Compare calculated drainage area with provided true drainage area */
#include <petscsys.h>
#include <petscvec.h>
#include "io.h"

PetscErrorCode
compare_true (Vec x, const char * const truefilename, 
    const PetscInt m, const PetscInt n, const PetscInt M, const PetscInt N,
    const MPI_Datatype doublefiletype, const PetscInt its, const int rank)
{

  Vec xtrue;
  PetscErrorCode ierr;
  PetscScalar *xptr;
  PetscReal norm;
  PetscReal max_error;
  PetscReal num_less_than_half;
  PetscInt i, j;

  /* Allocate a vector of the same size as x */
  ierr = VecDuplicate (x, &xtrue);
  CHKERRQ (ierr);
  ierr = VecSet (xtrue, 0.0);
  CHKERRQ (ierr);

  /* Get the pointer to the newly created vector */
  ierr = VecGetArray (xtrue, &xptr);
  CHKERRQ (ierr);

  /* Load the provided true drainage area file into xtrue */
  ierr = LoadFile (truefilename, "true", MPI_DOUBLE, doublefiletype,
      m * n, sizeof (double), xptr, rank);
  CHKERRQ (ierr);

  /* Convert the pointer back into a vector */
  ierr = VecRestoreArray (xtrue, &xptr);
  CHKERRQ (ierr);

  /* xtrue = xtrue - x */
  ierr = VecAXPY (xtrue, -1.0, x);
  CHKERRQ (ierr);
  /* norm = Norm of error */
  ierr = VecNorm (xtrue, NORM_2, &norm);
  CHKERRQ (ierr);
  /* max_error = Maximum absolute value of error */
  ierr = VecAbs (xtrue);
  CHKERRQ (ierr);
  ierr = VecMax (xtrue, PETSC_NULL, &max_error);

  /* Convert xtrue back to pointer to determine the percentage of cells
   * where the error is greater than 0.5 */
  ierr = VecGetArray (xtrue, &xptr);
  CHKERRQ (ierr);
  num_less_than_half = 0.0;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      if (xptr[i * m + j] < 0.5)
      {
        num_less_than_half += 1.0;
      }
    }
  }
  /* num_less_than_half contains local sum; do global reduce to sum globally */
  if (rank == 0)
    MPI_Reduce (MPI_IN_PLACE, &num_less_than_half, 1, MPI_DOUBLE, MPI_SUM,
        0, PETSC_COMM_WORLD);
  else
    MPI_Reduce (&num_less_than_half, NULL, 1, MPI_DOUBLE, MPI_SUM, 0,
        PETSC_COMM_WORLD);
  /* Divide by global number of cells to get percentage */
  num_less_than_half /= (M * N);
  /* Return pointer to vector */
  ierr = VecRestoreArray (xtrue, &xptr);
  CHKERRQ (ierr);
  /* Free xtrue */
  ierr = VecDestroy (&xtrue);
  CHKERRQ (ierr);

  /* Print result */
  ierr = PetscPrintf (PETSC_COMM_WORLD,
      "Error norm %g, Iterations %D, Max error %g, "
      "Fraction with error < 0.5 %g\n",
      norm, its, max_error, num_less_than_half);
  CHKERRQ (ierr);

  return 0;
}
