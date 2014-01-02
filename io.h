PetscErrorCode create_filetypes (const PetscInt fullM,
    const PetscInt fullN,
    const PetscInt localfullfirstrow_g,
    const PetscInt localfullfirstcol_g,
    const PetscInt localfirstrow_g,
    const PetscInt localfirstcol_g,
    const PetscInt M, const PetscInt N,
    const PetscInt m, const PetscInt n,
    MPI_Datatype * int8filetype,
    MPI_Datatype * doublefiletype,
    MPI_Datatype * outdoublefiletype);
PetscErrorCode destroy_filetypes (MPI_Datatype * int8filetype,
    MPI_Datatype * doublefiletype, MPI_Datatype * outdoublefiletype);
PetscErrorCode
LoadFile (const char * const filename, const char * const fileDescription, 
    const MPI_Datatype datatype, const MPI_Datatype filetype, 
    const PetscInt numToRead, const int bytesPerElement,
    void * restrict const data, const int rank);
PetscErrorCode
WriteFile (const char * const filename, const char * const fileDescription, 
    const MPI_Datatype datatype, const MPI_Datatype filetype, 
    const PetscInt numToWrite, const PetscScalar * restrict const data,
    const int rank);
