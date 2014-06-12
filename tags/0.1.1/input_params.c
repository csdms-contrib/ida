/* Read and process the input parameters specified by the user */
#include <petscsys.h>

PetscErrorCode
input_params (PetscInt * npx, PetscBool * npxflg,
    PetscInt * npy, PetscBool * npyflg, PetscInt * N, PetscInt * M,
    PetscInt * fullN, PetscInt * fullM,
    PetscInt * globalfirstrow_g, PetscInt * globalfirstcol_g,
    char * infilename, char * outfilename, char * truefilename,
    PetscBool * write_outfile_flg, PetscBool * compare_true_flg,
    PetscBool * crossborder, char * initguessfilename,
    PetscBool * init_guess_file_flg, PetscBool * initguesscalc,
    PetscBool * init_guess_flg, PetscBool * init_guess_calc_flg)
{

  PetscBool flg;
  PetscErrorCode ierr;

  ierr = PetscOptionsGetInt (PETSC_NULL, "-npx", npx, npxflg);
  CHKERRQ (ierr);
  ierr = PetscOptionsGetInt (PETSC_NULL, "-npy", npy, npyflg);
  CHKERRQ (ierr);
  ierr = PetscOptionsGetInt (PETSC_NULL, "-N", N, &flg);
  if (flg != PETSC_TRUE)
    SETERRQ (PETSC_COMM_WORLD, 85, "Need to set N");
  ierr = PetscOptionsGetInt (PETSC_NULL, "-M", M, &flg);
  if (flg != PETSC_TRUE)
    SETERRQ (PETSC_COMM_WORLD, 85, "Need to set M");
  ierr = PetscOptionsGetInt (PETSC_NULL, "-fullN", fullN, &flg);
  if (flg != PETSC_TRUE)
    *fullN = *N;
  ierr = PetscOptionsGetInt (PETSC_NULL, "-fullM", fullM, &flg);
  if (flg != PETSC_TRUE)
    *fullM = *M;
  ierr =
    PetscOptionsGetInt (PETSC_NULL, "-firstrow", globalfirstrow_g, &flg);
  if (flg != PETSC_TRUE)
    *globalfirstrow_g = 0;
  ierr =
    PetscOptionsGetInt (PETSC_NULL, "-firstcol", globalfirstcol_g, &flg);
  if (flg != PETSC_TRUE)
    *globalfirstcol_g = 0;
  ierr = PetscOptionsGetString (PETSC_NULL, "-in", infilename, 1024, &flg);
  if (flg != PETSC_TRUE)
    SETERRQ (PETSC_COMM_WORLD, 85, "Need to set -in");
  ierr = PetscOptionsGetString (PETSC_NULL, "-out", outfilename, 1024, 
      write_outfile_flg);
  ierr = PetscOptionsGetString (PETSC_NULL, "-true", truefilename, 1024, 
      compare_true_flg);
  *crossborder = PETSC_FALSE;
  ierr = PetscOptionsGetBool (PETSC_NULL, "-onlycrossborder", crossborder, 
      &flg);
  ierr = PetscOptionsGetString (PETSC_NULL, "-initguess", initguessfilename,
      1024, init_guess_file_flg);
  if (*crossborder)
  {
    *initguesscalc = PETSC_TRUE;
    *init_guess_flg = PETSC_TRUE;
  }
  else
  {
    *initguesscalc = PETSC_FALSE;
    ierr = PetscOptionsGetBool (PETSC_NULL, "-initguesscalc", initguesscalc,
        init_guess_calc_flg);
    *init_guess_flg = *init_guess_file_flg || *initguesscalc;
  }
  if (*init_guess_file_flg && *init_guess_calc_flg)
  {
    SETERRQ (PETSC_COMM_WORLD, 85,
        "Cannot use both -initguess and -initguesscalc");
  }

  return 0;

}
