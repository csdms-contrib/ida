/* Form the A matrix (regular IDA) or Ar matrix (hybrid IDA), such that
 * A * x = b or Ar * xr = br. A(r) therefore describes how a cell's drainage
 * area is related to that of other cells, which is derived from the flow
 * directions. */
#include <petscsys.h>
#include <petscmat.h>
#include <petscao.h>

/* Convert matrix of drainage directions into the A matrix */
  PetscErrorCode
formA (Mat A, const PetscInt M, const PetscInt N,
    const PetscInt m, const PetscInt n, const PetscInt mg,
    const PetscInt localfirstcol_l, const PetscInt localfirstrow_l,
    const PetscInt localfirstcol_g, const PetscInt localfirstrow_g,
    const unsigned char * restrict const W)
{
  PetscInt w_idx, idx, row, col;
  PetscInt dest_row, dest_col;
  PetscInt global_dest_row, global_dest_col;
  PetscScalar one = 1.0, neg_one = -1.0;
  PetscErrorCode ierr;
  PetscInt *cols, *rows;
  PetscScalar *values;
  const int numOnProc = m * n;
  PetscInt i;

  /* Allocate space to store column and row numbers and the associated value */
  cols = malloc (sizeof (PetscInt) * numOnProc * 2);
  rows = malloc (sizeof (PetscInt) * numOnProc * 2);
  values = malloc (sizeof (PetscScalar) * numOnProc * 2);

  idx = 0;
  for (row = 0; row < n; row++)
  {
    for (col = 0; col < m; col++)
    {

      /* Ones along diagonal: localfirst*_l are the first coordinates of the
       * local domain in local coordinates - 0 or 1 depending on whether
       * there are ghosts cells to the bottom and left of the local domain. The
       * local domain is of width mg (including ghost cells). */
      cols[idx] = (localfirstrow_l + row) * mg + (localfirstcol_l + col);
      rows[idx] = cols[idx];
      values[idx] = one;
      idx++;

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
        /* Column is my cell index, row is index of cell that I drain to */
        cols[idx] =
          (localfirstrow_l + row) * mg + (localfirstcol_l + col);
        rows[idx] =
          (localfirstrow_l + dest_row) * mg + (localfirstcol_l + dest_col);
        /* Put -1 at this location in the A matrix */
        values[idx] = neg_one;
        idx++;
      }

    }
  }

  /* Insert all of these values into the global A matrix, using local
   * coordinates */
  for (i = 0; i < idx; i++)
  {
    ierr =
      MatSetValuesLocal (A, 1, rows + i, 1, cols + i, values + i,
          INSERT_VALUES);
    CHKERRQ (ierr);
  }

  free (cols);
  free (rows);
  free (values);

  /* Assemble the A matrix */
  ierr = MatAssemblyBegin (A, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);
  ierr = MatAssemblyEnd (A, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);

  return 0;

}


/* Convert matrix of drainage directions into the Ar matrix, and also
 * determine the values of the br vector. */
  PetscErrorCode
formAr (Mat Ar, const PetscInt M, const PetscInt N,
    const PetscInt m, const PetscInt n, 
    const PetscInt localfirstcol_g, const PetscInt localfirstrow_g,
    const unsigned char * restrict const W, const AO ao,
    const PetscInt * restrict const aoordering, const PetscInt num_ao,
    PetscScalar * restrict const br_vals, 
    PetscInt * restrict const br_vals_idx,
    PetscInt * restrict const num_br_vals, Vec x)
{
  PetscInt idx_full_l, idx, i;
  PetscInt dest_row, dest_col;
  PetscInt dest_row_l, dest_col_l;
  PetscScalar one = 1.0, neg_one = -1.0;
  PetscErrorCode ierr;
  PetscInt *cols, *rows;
  PetscScalar *values;
  PetscInt rowidx;
  PetscInt idx_full_g;
  PetscInt idx_r_g;
  PetscInt row_full_g;
  PetscInt col_full_g;
  PetscInt row_full_l;
  PetscInt col_full_l;
  PetscInt dest_idx_full_g;
  PetscInt dest_idx_r_g;
  PetscScalar *xptr;

  /* Allocate space to store column and row numbers and the associated value */
  cols = malloc (sizeof (PetscInt) * m * n * 2);
  rows = malloc (sizeof (PetscInt) * m * n * 2);
  values = malloc (sizeof (PetscScalar) * m * n * 2);

  ierr = VecGetArray (x, &xptr);
  CHKERRQ (ierr);

  idx = 0;
  *num_br_vals = 0;
  for (rowidx = 0; rowidx < num_ao; rowidx++)
  {
    idx_full_g = aoordering [rowidx];
    idx_r_g = idx_full_g;
    ierr = AOApplicationToPetsc (ao, 1, &idx_r_g);
    CHKERRQ (ierr);
    row_full_g = idx_full_g / M;
    col_full_g = idx_full_g % M;
    row_full_l = row_full_g - localfirstrow_g;
    col_full_l = col_full_g - localfirstcol_g;
    idx_full_l = row_full_l * m + col_full_l;

    /* Ones along diagonal: localfirst*_l are the first coordinates of the
     * local domain in local coordinates - 0 or 1 depending on whether
     * there are ghosts cells to the bottom and left of the local domain. The
     * local domain is of width mg (including ghost cells). */
    cols[idx] = idx_r_g;
    rows[idx] = cols[idx];
    values[idx] = one;
    idx++;
    br_vals [*num_br_vals] = xptr [idx_full_l]; 
    br_vals_idx [*num_br_vals] = idx_r_g;
    (*num_br_vals)++;

    /* Determine coordinates of cell that I drain to */

    switch ((int) W[idx_full_l])
    {

      case 1:
        dest_col = col_full_g + 1;
        dest_row = row_full_g;
        break;

      case 2:
        dest_col = col_full_g + 1;
        dest_row = row_full_g + 1;
        break;

      case 4:
        dest_col = col_full_g;
        dest_row = row_full_g + 1;
        break;

      case 8:
        dest_col = col_full_g - 1;
        dest_row = row_full_g + 1;
        break;

      case 16:
        dest_col = col_full_g - 1;
        dest_row = row_full_g;
        break;

      case 32:
        dest_col = col_full_g - 1;
        dest_row = row_full_g - 1;
        break;

      case 64:
        dest_col = col_full_g;
        dest_row = row_full_g - 1;
        break;

      case 128:
        dest_col = col_full_g + 1;
        dest_row = row_full_g - 1;
        break;

      default:
        continue;
        break;

    }
    /* Check that we are not draining to a cell outside of domain */
    if (dest_row >= 0 && dest_row < N &&
        dest_col >= 0 && dest_col < M)
    {
      /* Column is my cell index, row is index of cell that I drain to */
      dest_idx_full_g = dest_row * M + dest_col;
      dest_idx_r_g = dest_idx_full_g;
      ierr = AOApplicationToPetsc (ao, 1, &dest_idx_r_g);
      CHKERRQ (ierr);
      cols[idx] = idx_r_g;
      rows[idx] = dest_idx_r_g;
      /* Put -1 at this location in the A matrix */
      values[idx] = neg_one;
      /* br_vals */
      dest_row_l = dest_row - localfirstrow_g;
      dest_col_l = dest_col - localfirstcol_g;
      if (dest_row_l >= 0 && dest_row_l < n &&
          dest_col_l >= 0 && dest_col_l < m)
      {
        br_vals [*num_br_vals] = -xptr [idx_full_l];
        br_vals_idx [*num_br_vals] = dest_idx_r_g;
        (*num_br_vals)++;
      }
      idx++;
    }
  }

  ierr = VecRestoreArray (x, &xptr);
  CHKERRQ (ierr);


  /* Insert all of these values into the global Ar matrix */
  for (i = 0; i < idx; i++)
  {
    ierr =
      MatSetValues (Ar, 1, rows + i, 1, cols + i, values + i,
          INSERT_VALUES);
    CHKERRQ (ierr);
  }

  free (cols);
  free (rows);
  free (values);

  /* Assemble the Ar matrix */
  ierr = MatAssemblyBegin (Ar, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);
  ierr = MatAssemblyEnd (Ar, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);

  return 0;

}
