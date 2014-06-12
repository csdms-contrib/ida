#include <unistd.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscao.h>

#include "io.h"
#include "compare_true.h"
#include "check.h"
#include "input_params.h"
#include "hybrid_ida.h"
#include "regular_ida.h"

/* IDA: Implicit Drainage Area
 * Alan Richardson, Chris Hill, Taylor Perron
 * Version 2013/12/19
 *
 * Compute drainage area for an input flow direction file
 *
 * fullN = height of matrix on disk
 * fullM = width of matrix on disk
 * globalfirstcol_g = first column of the matrix on disk that is to be computed
 * globalfirstrow_g = first row of the matrix on disk that is to be computed
 * M = width (number of columns) of matrix to be computed
 * N = height (number of rows) of matrix to be computed
 * m = width of my portion
 * n = height of my portion
 * localfullfirstcol_g = first column of my portion of matrix on disk in global
 *  coords
 * localfullfirstrow_g = first row of my portion of matrix on disk in global 
 *  coords
 * localfirstcol_g = first column of my portion of computed grid in global 
 *  coords
 * localfirstrow_g = first row of my portion of computed grid in global coords
 * localfirstcolghost_g = first column of my portion of computed grid, with 
 *  halo, in global coords
 * localfirstrowghost_g = first row of my portion of computed grid, with halo, 
 *  in global coords
 * localfirstcol_l = first column of my portion in local coords (may be 1 if 
 *  halo)
 * localfirstrow_l = first row of my portion in local coords (may be 1 if halo)
 */

static char help[] = "Solve for D8 area\n\
\tRequired input parameters:\n\
\t\t-in: Path to input flow directions file (8-bit int)\n\
\t\t-N: Number of landscape rows to be computed\n\
\t\t-M: Number of landscape columns to be computed\n\
\n\
\tOptional input parameters:\n\
\t\t-out: Path to output drainage area file (64-bit float)\n\
\t\t-true: Path to true drainage area file for comparison (64-bit float)\n\
\t\t-initguess: Path to initial guess of drainage area file (64-bit float)\n\
\t\t-npy: Number of processes to divide rows over\n\
\t\t-npx: Number of processes to divide columns over\n\
\t\t-fullN: Number of rows in input landscape\n\
\t\t-fullM: Number of columns in input landscape\n\
\t\t-firstrow: First row of input landscape to compute\n\
\t\t-firstcol: First column of input landscape to compute\n\
\t\t-initguesscalc: Locally calculate drainage area on each process with\n\
\t\t\tserial algorithm and use as an initial guess\n\
\t\t-onlycrossborder: Use hybrid IDA algorithm (recommended)\n\
\n\
\tSelected optional PETSc parameters:\n\
\t\t-ksp_type: solver (e.g. richardson)\n\
\t\t-pc_type: preconditioner (e.g. asm)\n\
\t\t-log_summary: Information (including performance data) about run\n";

static PetscErrorCode get_my_rows_and_cols (const PetscInt globalfirstrow_g,
    const PetscInt globalfirstcol_g,
    const PetscInt M,
    const PetscInt N,
    const PetscInt localfirstrow_g,
    const PetscInt localfirstcol_g,
    const PetscInt
    localfirstrowghost_g,
    const PetscInt
    localfirstcolghost_g,
    const PetscInt m,
    const PetscInt n,
    PetscInt * localfullfirstrow_g,
    PetscInt * localfullfirstcol_g,
    PetscInt * localfirstrow_l,
    PetscInt * localfirstcol_l);

  int
main (int argc, char **argv)
{

  Mat A;
  Vec x;
  KSP ksp;
  PetscErrorCode ierr;
  PetscInt its;
  PetscScalar *xptr;
  PetscBool compare_true_flg;
  PetscBool write_outfile_flg;
  PetscBool init_guess_flg;
  PetscBool init_guess_file_flg;
  PetscBool init_guess_calc_flg;
  PetscBool initguesscalc;
  PetscBool crossborder;
  PetscBool npxflg;
  PetscBool npyflg;
  PetscInt npx, npy;
  PetscInt fullN, fullM;
  PetscInt N, M;
  PetscInt n, m;
  PetscInt ng, mg;
  PetscInt globalfirstcol_g, globalfirstrow_g;
  PetscInt localfullfirstrow_g, localfullfirstcol_g;
  PetscInt localfirstcol_g, localfirstrow_g;
  PetscInt localfirstcolghost_g, localfirstrowghost_g;
  PetscInt localfirstcol_l, localfirstrow_l;
  PetscReal norm;
  DM da;
  MPI_Datatype int8filetype, doublefiletype, outdoublefiletype;
  int rank, size;
  unsigned char *W;
  char infilename[1024];
  char outfilename[1024];
  char truefilename[1024];
  char initguessfilename[1024];
  long long int t1, t2;
  struct timeval time;

  /* Initialize */
  PetscInitialize (&argc, &argv, (char *) 0, help);
  MPI_Comm_rank (PETSC_COMM_WORLD, &rank);
  MPI_Comm_size (PETSC_COMM_WORLD, &size);

  /* Read command line options */
  ierr = input_params (&npx, &npxflg, &npy, &npyflg, &N, &M, &fullN, &fullM,
    &globalfirstrow_g, &globalfirstcol_g, infilename, outfilename, truefilename,
    &write_outfile_flg, &compare_true_flg, &crossborder,
    initguessfilename, &init_guess_file_flg, &initguesscalc, &init_guess_flg,
    &init_guess_calc_flg);
  CHKERRQ (ierr);

  /* Initial checks such as that datatypes are of expected size */
  ierr = initial_check (write_outfile_flg, outfilename);
  CHKERRQ (ierr);

  /* Setup distributed arrays (DMDA) to store data */
  if (npxflg && npyflg)
  {
    /* User has specified x-y processor grid */
    ierr = DMDACreate2d (PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE,
        DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, M, N,
        npx, npy, 1, 1, PETSC_NULL, PETSC_NULL, &da);
  }
  else
  {
    /* Let PETSc decide x-y processor grid */
    ierr = DMDACreate2d (PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE,
        DMDA_BOUNDARY_NONE, DMDA_STENCIL_BOX, M, N,
        PETSC_DECIDE, PETSC_DECIDE, 1, 1, PETSC_NULL,
        PETSC_NULL, &da);
  }
  CHKERRQ (ierr);

  /* Create A matrix */
  ierr = DMCreateMatrix (da, MATAIJ, &A);
  CHKERRQ (ierr);

  /* Allocate corresponding x vector */
  ierr = MatGetVecs (A, &x, PETSC_NULL);
  CHKERRQ (ierr);

  /* Get global indices of lower left corner of my portion of domain,
   * including ghost cells (if present) */
  ierr = DMDAGetGhostCorners (da, &localfirstcolghost_g,
      &localfirstrowghost_g, PETSC_NULL, &mg, &ng, PETSC_NULL);
  CHKERRQ (ierr);

  /* As above, but exluding ghost cells */
  ierr = DMDAGetCorners (da, &localfirstcol_g, 
      &localfirstrow_g, PETSC_NULL, &m, &n, PETSC_NULL);
  CHKERRQ (ierr);

  /* Calculate information about the domain and my place in it */
  ierr = get_my_rows_and_cols (globalfirstrow_g, globalfirstcol_g,
      M, N, localfirstrow_g, localfirstcol_g,
      localfirstrowghost_g, localfirstcolghost_g, m,
      n, &localfullfirstrow_g, &localfullfirstcol_g,
      &localfirstrow_l, &localfirstcol_l);
  CHKERRQ (ierr);

  /* Create MPI datatypes for I/O */
  ierr = create_filetypes (fullM, fullN, localfullfirstrow_g,
      localfullfirstcol_g, localfirstrow_g,
      localfirstcol_g, M, N, m, n, &int8filetype,
      &doublefiletype, &outdoublefiletype);
  CHKERRQ (ierr);

  /* Load input drainage directions */
  W = malloc (sizeof (unsigned char) * m * n);
  if (W == NULL)
    SETERRQ (PETSC_COMM_WORLD, 55, "Error allocating W");
  ierr = LoadFile (infilename, "in", MPI_BYTE, int8filetype, m * n, 1, W, rank);
  CHKERRQ (ierr);

  /* Run IDA to solve for drainage area, either regular IDA or
   * the hybrid version */
  gettimeofday (&time, NULL);
  t1 = time.tv_sec * 1e6 + time.tv_usec;
  if (crossborder)
  {
    ierr = hybrid_ida (M, N, m, n, mg, localfirstcol_l, localfirstrow_l,
      localfirstcol_g, localfirstrow_g,
      init_guess_file_flg, initguesscalc, initguessfilename, doublefiletype, W,
      x, rank, &ksp);
  }
  else
  {
    ierr = regular_ida (A, M, N, m, n, mg, localfirstcol_l, localfirstrow_l,
      localfirstcol_g, localfirstrow_g, init_guess_flg,
      init_guess_file_flg, initguesscalc, initguessfilename, doublefiletype, W,
      x, rank, &ksp);
  }
  CHKERRQ (ierr);

  gettimeofday (&time, NULL);
  t2 = time.tv_sec * 1e6 + time.tv_usec;
  ierr = PetscPrintf (PETSC_COMM_WORLD, "IDA runtime: %lld\n", t2 - t1);
  CHKERRQ (ierr);

  /* Get the number of iterations that the solver ran for */
  ierr = KSPGetIterationNumber (ksp, &its);
  CHKERRQ (ierr);

  if (compare_true_flg == PETSC_TRUE)
  {
    /* Compare calculated drainage area with true solution */
    ierr = compare_true (x, truefilename, m, n, M, N, doublefiletype, its, 
        rank);
  }
  else
  {
    /* Just print final residual and number of iterations */
    ierr = KSPGetResidualNorm (ksp, &norm);
    CHKERRQ (ierr);
    ierr = PetscPrintf (PETSC_COMM_WORLD,
        "Error norm %g, Iterations %D\n", norm, its);
  }
  CHKERRQ (ierr);

  /* Free memory that is no longer needed */
  (void) free (W);
  ierr = MatDestroy (&A);
  CHKERRQ (ierr);
  ierr = KSPDestroy (&ksp);
  CHKERRQ (ierr);
  ierr = DMDestroy (&da);
  CHKERRQ (ierr);

  if (write_outfile_flg == PETSC_TRUE)
  {
    /* Write output drainage area to file */
    ierr = VecGetArray (x, &xptr);
    CHKERRQ (ierr);

    ierr = WriteFile (outfilename, "out", MPI_DOUBLE, outdoublefiletype,
        m * n, xptr, rank);
    CHKERRQ (ierr);

    ierr = VecRestoreArray (x, &xptr);
    CHKERRQ (ierr);
  }

  /* Free the x vector */
  ierr = VecDestroy (&x);
  CHKERRQ (ierr);

  /* Free MPI datatypes */
  ierr = destroy_filetypes (&int8filetype, &doublefiletype, &outdoublefiletype);
  CHKERRQ (ierr);

  /* Finalize */
  PetscFinalize ();

  return EXIT_SUCCESS;
}

  static PetscErrorCode
get_my_rows_and_cols (const PetscInt globalfirstrow_g,
    const PetscInt globalfirstcol_g, const PetscInt M,
    const PetscInt N, const PetscInt localfirstrow_g,
    const PetscInt localfirstcol_g,
    const PetscInt localfirstrowghost_g,
    const PetscInt localfirstcolghost_g, const PetscInt m,
    const PetscInt n, PetscInt * localfullfirstrow_g,
    PetscInt * localfullfirstcol_g,
    PetscInt * localfirstrow_l, PetscInt * localfirstcol_l)
{

  /* Coordinates of first local row and column in the full file */
  *localfullfirstrow_g = localfirstrow_g + globalfirstrow_g;
  *localfullfirstcol_g = localfirstcol_g + globalfirstcol_g;

  /* First coordinates of local domain in local coordinates: 0 or 1 depending
   * on presence of boundary of ghost cells */
  *localfirstrow_l = localfirstrow_g - localfirstrowghost_g;
  *localfirstcol_l = localfirstcol_g - localfirstcolghost_g;

  /* Ensure that my local domain does not extend beyond global domain */
  if (localfirstrow_g + n > N)
    return 1;
  if (localfirstcol_g + m > M)
    return 1;

  return 0;

}
