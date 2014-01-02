/* Regular IDA: Solve the whole system, even though many of the cells already
 * have the correct value */

#include <stdio.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscao.h>

#include "form_A_matrix.h"
#include "initial_guess.h"

PetscErrorCode
regular_ida (Mat A, const PetscInt M, const PetscInt N,
    const PetscInt m, const PetscInt n, const PetscInt mg,
    const PetscInt localfirstcol_l, const PetscInt localfirstrow_l,
    const PetscInt localfirstcol_g, const PetscInt localfirstrow_g,
    const PetscBool init_guess_flg,
    const PetscBool init_guess_file_flg, const PetscBool initguesscalc,
    const char * restrict const initguessfilename,
    const MPI_Datatype doublefiletype, const unsigned char * restrict const W,
    Vec x, const int rank, KSP * ksp)
{

  Vec b;
  PetscErrorCode ierr;

  /* Form the A matrix */
  ierr = formA (A, M, N, m, n, mg, localfirstcol_l, localfirstrow_l,
      localfirstcol_g, localfirstrow_g, W);
  CHKERRQ (ierr);

  /* Allocate corresponding b vector */
  ierr = MatGetVecs (A, PETSC_NULL, &b);
  CHKERRQ (ierr);

  /* Set the b vector (RHS) to be all ones */
  ierr = VecSet (b, 1.0);
  CHKERRQ (ierr);

  /* If using an initial guess (input file, or local solution), set x vector
   * to that, otherwise set it to 0 */
  if (init_guess_flg)
  {
    ierr = setInitialGuess (m, n, init_guess_file_flg, initguesscalc, 
        initguessfilename, doublefiletype, W, localfirstrow_l, mg, 
        localfirstcol_l, x, rank);
  }
  else
  {
    ierr = VecSet (x, 0.0);
  }
  CHKERRQ (ierr);

  /* Setup the solver */
  ierr = KSPCreate (PETSC_COMM_WORLD, ksp);
  CHKERRQ (ierr);
  ierr = KSPSetOperators (*ksp, A, A, DIFFERENT_NONZERO_PATTERN);
  CHKERRQ (ierr);
  ierr = KSPSetInitialGuessNonzero (*ksp, init_guess_flg);
  ierr = KSPSetFromOptions (*ksp);
  CHKERRQ (ierr);

  /* Run the solver */
  ierr = KSPSolve (*ksp, b, x);
  CHKERRQ (ierr);

  /* Free memory */
  ierr = VecDestroy (&b);
  CHKERRQ (ierr);

  return 0;

}
