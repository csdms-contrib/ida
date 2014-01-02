/* Hybrid IDA: We are only solving for the values of cells that could not be
 * determined with the local solver because they are downstream from
 * flow between processors; create the appropriate vectors/matrices
 * and solve for this reduced system
 */

#include <stdio.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscao.h>

#include "form_A_matrix.h"
#include "initial_guess.h"

/* Forward function declarations */
static PetscErrorCode
determine_cross_border_flow (const PetscInt M, const PetscInt N,
    const PetscInt m, const PetscInt n,
    const PetscInt localfirstcol_g, const PetscInt localfirstrow_g,
    const unsigned char * restrict const W,
    PetscInt * num_send_across_border, PetscInt * send_across_border_idx,
    Vec rcv_across_border, const AO dmao);


  PetscErrorCode
hybrid_ida (const PetscInt M, const PetscInt N,
    const PetscInt m, const PetscInt n, const PetscInt mg,
    const PetscInt localfirstcol_l, const PetscInt localfirstrow_l,
    const PetscInt localfirstcol_g, const PetscInt localfirstrow_g,
    const PetscBool init_guess_file_flg, const PetscBool initguesscalc,
    const char * restrict const initguessfilename,
    const MPI_Datatype doublefiletype, const unsigned char * restrict const W,
    Vec x, const int rank, KSP * ksp, const DM da)
{

  Mat Ar;
  Vec br;
  Vec xr;
  Vec rcv_across_border;
  PetscErrorCode ierr;
  PetscInt num_send_across_border;
  PetscInt num_only_send_across_border;
  PetscInt num_send_rcv_across_border;
  PetscInt num_rcv_across_border;
  PetscInt num_ao;
  PetscInt num_br_vals;
  const PetscInt num_border_cells = 2 * m + 2 * (n - 2);
  PetscInt i;
  PetscInt * send_across_border_idx;
  PetscInt * only_send_across_border_idx;
  PetscInt * rcv_across_border_idx;
  PetscInt * br_vals_idx;
  PetscInt * aoordering;
  PetscInt * dmordering;
  PetscScalar * only_send_across_border;
  PetscScalar * br_vals;
  PetscScalar *xptr;
  AO ao;
  AO dmao;

  /* Allocate memory to store bookkeeping vectors related to only using
   * the linear solver to set cells that cross borders */

  /* rcv_across_border: 1.0 if directly receive flow from across process
   *  boundary, 0.0 otherwise */
  /* TODO: This vector could be significantly smaller: it is currently
   * as large as the local domain (m * n), but only the elements on the
   * borders (2 * m + 2 * (n - 2)) have the possibility of receiving flow
   * from across a boundary, and so those are the only ones that actually
   * need to be stored. Changing this would require also changing the section
   * in determine_cross_border_flow where the values are set and
   * setInitialGuess_crossborder, where they are used. */
  ierr = VecCreate (PETSC_COMM_WORLD, &rcv_across_border);
  CHKERRQ (ierr);
  ierr = VecSetSizes (rcv_across_border, m * n, PETSC_DECIDE);
  CHKERRQ (ierr);
  ierr = VecSetFromOptions (rcv_across_border);
  CHKERRQ (ierr);
  ierr = VecSet (rcv_across_border, 0.0);
  CHKERRQ (ierr);

  /* send_across_border_idx: global indices of cells directly sending flow
   *  across process border */
  send_across_border_idx = malloc (sizeof (PetscInt) * num_border_cells);
  if (send_across_border_idx == NULL)
    SETERRQ (PETSC_COMM_WORLD, 55, "Error allocating send_across_border_idx");

  /* only_send_across_border_idx: global indices of cells that directly
   *  send flow across process borders, but do not also receive any flow
   *  from across process borders (so will be set correctly by local
   *  serial calculation */
  only_send_across_border_idx = malloc (sizeof (PetscInt) * num_border_cells);
  if (only_send_across_border_idx == NULL)
    SETERRQ (PETSC_COMM_WORLD, 55, "Error allocating "
        "only_send_across_border_idx");

  /* rcv_across_border_idx: */
  rcv_across_border_idx = malloc (sizeof (PetscInt) * m * n);
  if (rcv_across_border_idx == NULL)
    SETERRQ (PETSC_COMM_WORLD, 55, "Error allocating rcv_across_border_idx");

  /* only_send_across_border: */
  only_send_across_border = malloc (sizeof (PetscScalar) * num_border_cells);
  if (only_send_across_border == NULL)
    SETERRQ (PETSC_COMM_WORLD, 55,
        "Error allocating only_send_across_border");


  /* Get the mapping between the DMDA indices and PETSCs's indices */
  ierr = DMDAGetAO (da, &dmao);

  /* Determine cells that send and receive across process borders */
  ierr = determine_cross_border_flow (M, N, m, n, 
      localfirstcol_g, localfirstrow_g, W,
      &num_send_across_border, send_across_border_idx,
      rcv_across_border, dmao);
  CHKERRQ (ierr);

  /* If using an initial guess (input file, or local solution), set x vector
   * to that, otherwise set it to 0 */
  ierr = setInitialGuess_crossborder (m, n, init_guess_file_flg, initguesscalc,
      initguessfilename, doublefiletype, W, localfirstrow_l, mg, 
      localfirstcol_l, x, rank, rcv_across_border,
      &num_only_send_across_border, &num_send_rcv_across_border,
      num_send_across_border, send_across_border_idx,
      only_send_across_border_idx, M, localfirstrow_g,
      localfirstcol_g, only_send_across_border, &num_rcv_across_border,
      rcv_across_border_idx);
  CHKERRQ (ierr);

  /* Allocate a vector large enough to store the global index of each cell */
  if (num_only_send_across_border + num_send_rcv_across_border !=
      num_send_across_border)
    SETERRQ (PETSC_COMM_WORLD, 85, "only_send + send_rcv != send");
  num_ao = num_rcv_across_border + num_only_send_across_border;
  aoordering = malloc (sizeof (PetscInt) * num_ao);

  /* Determine the number of cells that could not be solved for locally,
   * and copy their global index into the aoordering vector */
  for (i = 0; i < num_only_send_across_border; i++)
  {
    aoordering [i] = only_send_across_border_idx [i];
  }
  for (i = 0; i < num_rcv_across_border; i++)
  {
    aoordering [num_only_send_across_border + i] = rcv_across_border_idx [i];
  }

  /* Create AO ordering */
  ierr = AOCreateMapping (PETSC_COMM_WORLD, num_ao, aoordering, NULL,
      &ao);
  CHKERRQ (ierr);

  /* Allocate reduced A: Ar */
  /* TODO: could use smaller, more precise preallocation */
  ierr = MatCreate (PETSC_COMM_WORLD, &Ar);
  CHKERRQ (ierr);
  ierr = MatSetSizes (Ar, num_ao, num_ao, PETSC_DECIDE, PETSC_DECIDE);
  CHKERRQ (ierr);
  ierr = MatSetFromOptions (Ar);
  CHKERRQ (ierr);
  ierr = MatMPIAIJSetPreallocation (Ar, 8, PETSC_NULL, 8, PETSC_NULL);
  CHKERRQ (ierr);
  ierr = MatSeqAIJSetPreallocation (Ar, 8, PETSC_NULL);
  CHKERRQ (ierr);

  /* br_vals: */
  br_vals = malloc (sizeof (PetscScalar) * num_ao * 2);
  if (br_vals == NULL)
    SETERRQ (PETSC_COMM_WORLD, 55, "Error allocating br_vals");

  /* br_vals_idx: */
  br_vals_idx = malloc (sizeof (PetscInt) * num_ao * 2);
  if (br_vals_idx == NULL)
    SETERRQ (PETSC_COMM_WORLD, 55, "Error allocating br_vals_idx");

  /* Form reduced A */
  ierr = formAr (Ar, M, N, m, n,
      localfirstcol_g, localfirstrow_g, W, ao, aoordering, num_ao,
      br_vals, br_vals_idx, &num_br_vals, x);
  CHKERRQ (ierr);

  /* Allocate reduced x and b */
  ierr = VecCreate (PETSC_COMM_WORLD, &xr);
  CHKERRQ (ierr);
  ierr = VecSetSizes (xr, num_ao, PETSC_DECIDE);
  CHKERRQ (ierr);
  ierr = VecSetFromOptions (xr);
  CHKERRQ (ierr);
  ierr = VecSet (xr, 0.0);
  CHKERRQ (ierr);
  ierr = VecCreate (PETSC_COMM_WORLD, &br);
  CHKERRQ (ierr);
  ierr = VecSetSizes (br, num_ao, PETSC_DECIDE);
  CHKERRQ (ierr);
  ierr = VecSetFromOptions (br);
  CHKERRQ (ierr);
  ierr = VecSet (br, 0.0);
  CHKERRQ (ierr);
  ierr = VecSetValues (br, num_br_vals, br_vals_idx, br_vals, ADD_VALUES);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (br);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (br);
  CHKERRQ (ierr);

  /* Setup the solver */
  ierr = KSPCreate (PETSC_COMM_WORLD, ksp);
  CHKERRQ (ierr);
  ierr = KSPSetOperators (*ksp, Ar, Ar, DIFFERENT_NONZERO_PATTERN);
  CHKERRQ (ierr);
  ierr = KSPSetInitialGuessNonzero (*ksp, PETSC_FALSE);
  ierr = KSPSetFromOptions (*ksp);
  CHKERRQ (ierr);

  /* Run the solver */
  ierr = KSPSolve (*ksp, br, xr);
  CHKERRQ (ierr);

  /* Copy the values from the reduced system into the full system */
  if (num_rcv_across_border > 0)
  {
    ierr = VecGetArray (xr, &xptr);
    CHKERRQ (ierr);
    dmordering = malloc (sizeof (PetscInt) * num_rcv_across_border);
    memcpy (dmordering, aoordering + num_only_send_across_border,
        num_rcv_across_border * sizeof (PetscInt));
    ierr = AOApplicationToPetsc (dmao, num_rcv_across_border, dmordering);
    CHKERRQ (ierr);
    ierr = VecSetValues (x, num_rcv_across_border, 
        dmordering, xptr+num_only_send_across_border, INSERT_VALUES);
    CHKERRQ (ierr);
    ierr = VecRestoreArray (xr, &xptr);
    CHKERRQ (ierr);
    (void) free (dmordering);
  }

  /* Free memory */
  (void) free (aoordering);
  (void) free (send_across_border_idx);
  (void) free (only_send_across_border_idx);
  (void) free (only_send_across_border);
  (void) free (rcv_across_border_idx);
  (void) free (br_vals);
  (void) free (br_vals_idx);
  ierr = MatDestroy (&Ar);
  CHKERRQ (ierr);
  ierr = VecDestroy (&br);
  CHKERRQ (ierr);
  ierr = VecDestroy (&xr);
  CHKERRQ (ierr);
  ierr = VecDestroy (&rcv_across_border);
  CHKERRQ (ierr);

  return 0;
}

/* Determine which cells directly send and receive flow across a
 * processor domain border */
  static PetscErrorCode
determine_cross_border_flow (const PetscInt M, const PetscInt N,
    const PetscInt m, const PetscInt n,
    const PetscInt localfirstcol_g, const PetscInt localfirstrow_g,
    const unsigned char * restrict const W,
    PetscInt * num_send_across_border, PetscInt * send_across_border_idx,
    Vec rcv_across_border, const AO dmao)
{
  PetscInt w_idx, row, col;
  PetscInt dest_row, dest_col;
  PetscInt global_dest_row, global_dest_col;
  PetscErrorCode ierr;
  const int numOnProc = m * n;
  PetscInt i;
  PetscInt * dest_across_border_idx;
  PetscScalar one = 1.0;
  PetscScalar * ones;

  /* Allocate space to store column and row numbers and the associated value */
  dest_across_border_idx = malloc (sizeof (PetscInt) * numOnProc);

  *num_send_across_border = 0;
  for (row = 0; row < n; row++)
  {
    for (col = 0; col < m; col++)
    {

      /* Ones along diagonal: localfirst*_l are the first coordinates of the
       * local domain in local coordinates - 0 or 1 depending on whether
       * there are ghosts cells to the bottom and left of the local domain. The
       * local domain is of width mg (including ghost cells). */

      /* Determine coordinates of cell that I drain to */
      w_idx = row * m + col;
      switch ((int) W[w_idx])
      {

        case 1:
          dest_col = col + 1;
          dest_row = row;
          break;

        case 2:
          dest_col = col + 1;
          dest_row = row + 1;
          break;

        case 4:
          dest_col = col;
          dest_row = row + 1;
          break;

        case 8:
          dest_col = col - 1;
          dest_row = row + 1;
          break;

        case 16:
          dest_col = col - 1;
          dest_row = row;
          break;

        case 32:
          dest_col = col - 1;
          dest_row = row - 1;
          break;

        case 64:
          dest_col = col;
          dest_row = row - 1;
          break;

        case 128:
          dest_col = col + 1;
          dest_row = row - 1;
          break;

        default:
          continue;
          break;

      }
      global_dest_col = dest_col + localfirstcol_g;
      global_dest_row = dest_row + localfirstrow_g;

      /* Check that we are not draining to a cell outside of domain */
      if (global_dest_row >= 0 && global_dest_row < N &&
          global_dest_col >= 0 && global_dest_col < M)
      {
        if ((dest_col < 0) || (dest_col >= m) || 
            (dest_row < 0) || (dest_row >= n)) /* outside local domain */
        {
          send_across_border_idx[*num_send_across_border]=
            (row + localfirstrow_g) * M + (col + localfirstcol_g);
          dest_across_border_idx[*num_send_across_border]=
            global_dest_row * M + global_dest_col;
          (*num_send_across_border)++;
        }
      }

    }
  }


  /* Insert into rcv_across_border vector. This is used in
   * setInitialGuess_crossborder to determine which cells receive flow
   * from across borders. */
  ierr = AOApplicationToPetsc (dmao, *num_send_across_border,
      dest_across_border_idx);
  CHKERRQ (ierr);

  ones = malloc (sizeof (PetscScalar) * *num_send_across_border);
  if (ones == NULL)
    SETERRQ (PETSC_COMM_WORLD, 55, "Error allocating ones");
  for (i = 0; i < *num_send_across_border; i++)
  {
    ones [i] = one;
  }

  VecSetValues (rcv_across_border, *num_send_across_border, 
      dest_across_border_idx, ones, ADD_VALUES);
  ierr = VecAssemblyBegin (rcv_across_border);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (rcv_across_border);
  CHKERRQ (ierr);

  free (dest_across_border_idx);
  free (ones);

  return 0;

}
