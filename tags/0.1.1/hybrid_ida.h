  PetscErrorCode
hybrid_ida (const PetscInt M, const PetscInt N,
    const PetscInt m, const PetscInt n, const PetscInt mg,
    const PetscInt localfirstcol_l, const PetscInt localfirstrow_l,
    const PetscInt localfirstcol_g, const PetscInt localfirstrow_g,
    const PetscBool init_guess_file_flg, const PetscBool initguesscalc,
    const char * restrict const initguessfilename,
    const MPI_Datatype doublefiletype, const unsigned char * restrict const W,
    Vec x, const int rank, KSP * ksp);
