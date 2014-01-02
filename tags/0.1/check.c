/* Perform some simple checks to detect potential problems */

#include <stdio.h>
#include <petscsys.h>

/* Forward function declarations */
static PetscErrorCode
check_datatypes (void);
static PetscErrorCode
check_output_openable (const char * const outfilename);

PetscErrorCode initial_check (const PetscBool write_outfile_flg,
    const char * const outfilename)
{

  PetscErrorCode ierr;

  /* Check data types */
  ierr = check_datatypes ();
  CHKERRQ (ierr);

  /* Check can open output file, if any */
  if (write_outfile_flg == PETSC_TRUE)
  {
    ierr = check_output_openable (outfilename);
    CHKERRQ (ierr);
  }

  return 0;
}

static PetscErrorCode
check_datatypes (void)
{

  /* data type used for input flow directions */
  if (sizeof (unsigned char) != 1)
  {
    SETERRQ1 (PETSC_COMM_WORLD, 85, "Code assumes that"
        " sizeof (unsigned char) == 1,"
        " but, it is %zu on this system. Consider altering the code or"
        " running on a different system\n", sizeof (unsigned char));
  }

  /* data type used for output, and input true answer and initial guess */
  if (sizeof (PetscScalar) != 8)
  {
    SETERRQ1 (PETSC_COMM_WORLD, 85, "Code assumes that"
        " sizeof (PetscScalar) == 8,"
        " but, it is %zu on this system. Consider altering the code or"
        " re-configuring PETSc\n", sizeof (PetscScalar));
  }

  return 0;
}

static PetscErrorCode
check_output_openable (const char * const outfilename)
{

  MPI_File fh;

  /* Open the output file */
  MPI_File_open (PETSC_COMM_WORLD, outfilename,
      MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  if (fh == NULL)
    SETERRQ1 (PETSC_COMM_WORLD, 65, "Error opening output file %s\n",
        outfilename);
  /* Close it again */
  MPI_File_close (&fh);

  return 0;
}
