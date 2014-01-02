/* File I/O: reading input files (flow directions, initial guess, true answer),
 * and write output files (drainage area computed with IDA) */
#include <sys/stat.h>
#include <petscsys.h>

static PetscInt fullM;
static PetscInt fullN;

/* Create MPI datatypes to allow each processor to read and write to the
 * locations in the files corresponding to its local domain */
  PetscErrorCode
create_filetypes (const PetscInt fullM_in, const PetscInt fullN_in,
    const PetscInt localfullfirstrow_g,
    const PetscInt localfullfirstcol_g,
    const PetscInt localfirstrow_g,
    const PetscInt localfirstcol_g, const PetscInt M,
    const PetscInt N, const PetscInt m, const PetscInt n,
    MPI_Datatype * int8filetype, MPI_Datatype * doublefiletype,
    MPI_Datatype * outdoublefiletype)
{

  int sizes[2] = { fullN_in, fullM_in };
  int subsizes[2] = { n, m };
  int starts[2] = { localfullfirstrow_g, localfullfirstcol_g };

  /* Put full file size values in file scope variables for use later in
   * LoadFile */
  fullM = fullM_in;
  fullN = fullN_in;

  /* Input flow directions, stored in int8 datatype */
  MPI_Type_create_subarray (2, sizes, subsizes, starts,
      MPI_ORDER_C, MPI_BYTE, int8filetype);

  /* Input drainage area, either an initial guess, or the true values, in
   * double precision (64 bit) floating point datatype */
  MPI_Type_create_subarray (2, sizes, subsizes, starts,
      MPI_ORDER_C, MPI_DOUBLE, doublefiletype);

  /* Output file is only the size of the global domain, which might be
   * different from the size of the full input files */
  sizes[0] = N;
  sizes[1] = M;
  subsizes[0] = n;
  subsizes[1] = m;
  starts[0] = localfirstrow_g;
  starts[1] = localfirstcol_g;

  /* Output drainage area in double precision datatype */
  MPI_Type_create_subarray (2, sizes, subsizes, starts,
      MPI_ORDER_C, MPI_DOUBLE, outdoublefiletype);

  /* Commit MPI datatypes */
  MPI_Type_commit (int8filetype);
  MPI_Type_commit (doublefiletype);
  MPI_Type_commit (outdoublefiletype);

  return 0;

}

/* Free MPI datatypes */
  PetscErrorCode
destroy_filetypes (MPI_Datatype * int8filetype, MPI_Datatype * doublefiletype,
    MPI_Datatype * outdoublefiletype)
{
  MPI_Type_free (int8filetype);
  MPI_Type_free (doublefiletype);
  MPI_Type_free (outdoublefiletype);

  return 0;
}

/* Load an input file: flow directions, initial guess, true drainage area.
 * MPI_IO is used to enable each process to load its local portion in
 * parallel */
  PetscErrorCode
LoadFile (const char * const filename, const char * const fileDescription, 
    const MPI_Datatype datatype, const MPI_Datatype filetype, 
    const PetscInt numToRead, const int bytesPerElement,
    void * restrict const data, const int rank)
{

  MPI_File fh;
  MPI_Status status;
  int count;
  struct stat file_stat;

  /* Check that the file is the right size */
  if (rank == 0)
  {
    stat (filename, &file_stat);
    if (file_stat.st_size != (off_t) (fullM * fullN * bytesPerElement))
      SETERRQ4 (PETSC_COMM_WORLD, 65, "Expected %s file %s to contain "
          "%D bytes, but contains %jd bytes", fileDescription, filename,
          fullM * fullN * bytesPerElement, (intmax_t) file_stat.st_size);
  }

  /* Open the file */
  MPI_File_open (PETSC_COMM_WORLD, filename,
      MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  if (fh == NULL)
    SETERRQ2 (PETSC_COMM_WORLD, 65, "Error opening %s file %s\n",
        fileDescription, filename);

  /* Specify the part of the file to load onto this process */
  MPI_File_set_view (fh, 0, datatype, filetype, "native",
      MPI_INFO_NULL);

  /* Read the relevant portion of the file and close it */
  MPI_File_read_all (fh, data, numToRead, datatype, &status);
  MPI_File_close (&fh);

  /* Ensure that the correct number of elements were read */
  MPI_Get_count (&status, datatype, &count);
  if (count != (int) numToRead)
  {
    SETERRQ5 (PETSC_COMM_WORLD, 66, "Rank %d: Read %d from %s file %s, "
        "expected %D\n", rank, count, fileDescription, filename, numToRead);
  }

  return 0;
}

/* Write the output drainage area */
  PetscErrorCode
WriteFile (const char * const filename, const char * const fileDescription, 
    const MPI_Datatype datatype, const MPI_Datatype filetype, 
    const PetscInt numToWrite, const PetscScalar * restrict const data,
    const int rank)
{

  MPI_File fh;
  MPI_Status status;
  int count;

  /* Open the output file */
  MPI_File_open (PETSC_COMM_WORLD, filename,
      MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  if (fh == NULL)
    SETERRQ2 (PETSC_COMM_WORLD, 65, "Error opening %s file %s\n",
        fileDescription, filename);

  /* Specify the portion this process should write to */
  MPI_File_set_view (fh, 0, datatype, filetype,
      "native", MPI_INFO_NULL);

  /* Write data and close file */
  MPI_File_write_all (fh, data, numToWrite, datatype, &status);
  MPI_File_close (&fh);

  /* Check that the correct number of elements were written */
  MPI_Get_count (&status, datatype, &count);
  if (count != (int) numToWrite)
  {
    SETERRQ5 (PETSC_COMM_WORLD, 67, "Rank %d: Wrote %d to %s file %s, "
        "expected %D\n", rank, count, fileDescription, filename, numToWrite);
  }

  return 0;

}
