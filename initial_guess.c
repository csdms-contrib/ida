/* Provide an initial guess of the drainage area to the solver. This is
 * either done by loading an initial guess from a file, or by using the
 * serial algorithm to solve for drainage area locally on each processor
 * (ignoring flow from outside the processor's domain).*/
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <petscsys.h>
#include <petscvec.h>
#include "exact_solve.h"
#include "io.h"

/* Forward function declarations */
static PetscErrorCode
set_init_guess (const PetscBool init_guess_file_flg,
    const char * const initguessfilename, const MPI_Datatype doublefiletype,
    const PetscInt m, const PetscInt n, double * restrict const tmp,
    const int rank, const PetscBool initguesscalc,
    const unsigned char * restrict const W);

static PetscErrorCode
copy_init_guess_to_x (const PetscInt m, const PetscInt n, const PetscInt mg,
    const PetscInt localfirstrow_l, const PetscInt localfirstcol_l,
    Vec x, const double * restrict const tmp);

/* Set the x vector to an initial guess, either a provided drainage area
 * file, or a local (serial) computation on each processor. This is the
 * regular IDA version. */
  PetscErrorCode
setInitialGuess (const PetscInt m, const PetscInt n, 
    const PetscBool init_guess_file_flg, const PetscBool initguesscalc, 
    const char * const initguessfilename, const MPI_Datatype doublefiletype,
    const unsigned char * restrict const W, const PetscInt localfirstrow_l, 
    const PetscInt mg, const PetscInt localfirstcol_l, Vec x, const int rank)
{
  PetscErrorCode ierr;
  PetscScalar *tmp;

  /* Allocate a temporary array to store initial guess */
  tmp = malloc (sizeof (PetscScalar) * m * n);
  if (tmp == NULL)
    SETERRQ (PETSC_COMM_WORLD, 55, "Error allocating tmp");
  memset (tmp, 0, sizeof (PetscScalar) * m * n);

  /* Get initial guess */
  ierr = set_init_guess (init_guess_file_flg, initguessfilename, doublefiletype,
      m, n, tmp, rank, initguesscalc, W);
  CHKERRQ (ierr);

  /* Copy initial guess to x vector */
  ierr = copy_init_guess_to_x (m, n, mg, localfirstrow_l, localfirstcol_l, x,
      tmp);
  CHKERRQ (ierr);

  /* Free temporary arrays */
  free (tmp);

  return 0;
}


/* Set the x vector to an initial guess, either a provided drainage area
 * file, or a local (serial) computation on each processor. This is the
 * hybrid IDA version. */
  PetscErrorCode
setInitialGuess_crossborder (const PetscInt m, const PetscInt n, 
    const PetscBool init_guess_file_flg, const PetscBool initguesscalc, 
    const char * const initguessfilename, const MPI_Datatype doublefiletype,
    const unsigned char * restrict const W, const PetscInt localfirstrow_l, 
    const PetscInt mg, const PetscInt localfirstcol_l, Vec x, 
    const int rank, Vec rcv_across_border,
    PetscInt * restrict const num_only_send_across_border, 
    PetscInt * restrict const num_send_rcv_across_border,
    const PetscInt num_send_across_border,
    const PetscInt * restrict const send_across_border_idx,
    PetscInt * restrict const only_send_across_border_idx,
    const PetscInt M, const PetscInt localfirstrow_g,
    const PetscInt localfirstcol_g,
    PetscScalar * restrict const only_send_across_border,
    PetscInt * restrict const num_rcv_across_border,
    PetscInt * restrict const rcv_across_border_idx)
{
  PetscErrorCode ierr;
  PetscScalar *tmp;
  PetscScalar *rcv_across_border_ptr;
  PetscInt global_idx;
  PetscInt local_row_l;
  PetscInt local_col_l;
  PetscInt local_idx;
  int i, j;

  /* Allocate a temporary array to store initial guess */
  tmp = malloc (sizeof (PetscScalar) * m * n);
  if (tmp == NULL)
    SETERRQ (PETSC_COMM_WORLD, 55, "Error allocating tmp");
  memset (tmp, 0, sizeof (PetscScalar) * m * n);

  /* Set initial guess for cells that receive flow from across a processor
   * boundary to be NAN, so the cells downstream from these cells will
   * also be set to NAN, allowing cells without the correct initial
   * guess to be identified */
  ierr = VecGetArray (rcv_across_border, &rcv_across_border_ptr);
  CHKERRQ (ierr);
  /* Loop over borders, since those are the only locations that could receive
   * cross border flow */
  for (i = 0; i < m; i++)
  {
    /* top */
    if (rcv_across_border_ptr [i] > 0.5) /* >= 1.0 */
    {
      tmp [i] = NAN;
    }
    /* bottom */
    if (rcv_across_border_ptr [(n - 1) * m + i] > 0.5) /* >= 1.0 */
    {
      tmp [(n - 1) * m + i] = NAN;
    }
  }
  for (i = 0; i < n; i++)
  {
    /* left */
    if (rcv_across_border_ptr [i * m] > 0.5) /* >= 1.0 */
    {
      tmp [i * m] = NAN;
    }
    /* right */
    if (rcv_across_border_ptr [(m - 1) + i * m] > 0.5) /* >= 1.0 */
    {
      tmp [(m - 1) + i * m] = NAN;
    }
  }
  ierr = VecRestoreArray (rcv_across_border, &rcv_across_border_ptr);
  CHKERRQ (ierr);

  /* Get initial guess */
  ierr = set_init_guess (init_guess_file_flg, initguessfilename, doublefiletype,
      m, n, tmp, rank, initguesscalc, W);
  CHKERRQ (ierr);

  /* For the cells that send flow across processor boundaries, check
   * whether their initial guess is NAN (i.e. they also receive flow that
   * crossed processor boundaries). If not, the cell only sends and does
   * not receive flow from across processor boundaries, so it has the
   * correct value. */
  *num_only_send_across_border = 0;
  *num_send_rcv_across_border = 0;
  for (i = 0; i < num_send_across_border; i++)
  {
    global_idx = send_across_border_idx [i];
    local_row_l = global_idx / M - localfirstrow_g;
    local_col_l = global_idx % M - localfirstcol_g;
    local_idx = local_row_l * m + local_col_l;
    if (! isnan (tmp [local_idx]))
    {
      only_send_across_border_idx [*num_only_send_across_border] = 
        send_across_border_idx [i];
      only_send_across_border [*num_only_send_across_border] =
        tmp [local_idx];
      (*num_only_send_across_border)++;
    }
    else
    {
      (*num_send_rcv_across_border)++;
    }
  }

  /* Find all cells that are NAN, and add them to the list of cells that
   * receive flow from across a border */
  *num_rcv_across_border = 0;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      if (isnan (tmp [i * m + j]))
      {
        rcv_across_border_idx [*num_rcv_across_border] =
          (i + localfirstrow_g) * M + (j + localfirstcol_g);
        (*num_rcv_across_border)++;
        /* Set the cells that were NAN to be 0 so that when we run the
         * serial solver again (below) these cells are updated with a guess */
        tmp [i * m + j] = 0.0;
      }
    }
  }

  /* Run serial method again to get guess for cells that receive flow from
   * across a processor domain boundary. The values calculated will be
   * incorrect, as the cross border flow is neglected. It may be faster to
   * skip this step, as the reduction in runtime of the solver may not
   * outweigh its cost. */
  exact_solve (tmp, W, m, n);

  /* Copy initial guess to x vector */
  ierr = copy_init_guess_to_x (m, n, mg, localfirstrow_l, localfirstcol_l, x,
      tmp);
  CHKERRQ (ierr);

  /* Free temporary arrays */
  free (tmp);

  return 0;
}

/* Load the initial guess from a file, or run the serial algorithm, depending
 * on what the user has requested */
  static PetscErrorCode
set_init_guess (const PetscBool init_guess_file_flg,
    const char * const initguessfilename, const MPI_Datatype doublefiletype,
    const PetscInt m, const PetscInt n, double * restrict const tmp,
    const int rank, const PetscBool initguesscalc,
    const unsigned char * restrict const W)
{

  PetscErrorCode ierr;

  if (init_guess_file_flg)
  {
    /* Use provided file */
    ierr = LoadFile (initguessfilename, "init guess", MPI_DOUBLE, 
        doublefiletype, m * n, sizeof (double), tmp, rank);
    CHKERRQ (ierr);
  }
  else if (initguesscalc)
  {
    /* Solve drainage area locally on each processor (drainage across
     * processor boundaries will be incorrect) */

    exact_solve (tmp, W, m, n);
  }
  else
  {
    /* Should never get here */
    SETERRQ (PETSC_COMM_WORLD, 56, "init_guess_flg set, but neither"
        "init_guess_file_flg nor initguesscalc set");
  }

  return 0;

}

/* Copy the initial guess into the x vector */
  static PetscErrorCode
copy_init_guess_to_x (const PetscInt m, const PetscInt n, const PetscInt mg,
    const PetscInt localfirstrow_l, const PetscInt localfirstcol_l,
    Vec x, const double * restrict const tmp)
{

  PetscInt *indices;
  PetscErrorCode ierr;
  int i, j, k;

  /* Create array of local coordinate indices for initial guess values */
  indices = malloc (sizeof (PetscInt) * m * n);
  if (indices == NULL)
    SETERRQ (PETSC_COMM_WORLD, 55, "Error allocating indices");

  k = 0;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      indices[k] = (localfirstrow_l + i) * mg + (localfirstcol_l + j);
      k++;
    }
  }

  /* Insert initial guess values into x vector */
  ierr = VecSetValuesLocal (x, m * n, indices, tmp, INSERT_VALUES);
  CHKERRQ (ierr);

  /* Free temporary arrays */
  free (indices);

  /* Assemble x vector */
  ierr = VecAssemblyBegin (x);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (x);
  CHKERRQ (ierr);

  return 0;

}
