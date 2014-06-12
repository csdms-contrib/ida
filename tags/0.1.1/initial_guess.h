PetscErrorCode
setInitialGuess (const PetscInt m, const PetscInt n, 
    const PetscBool init_guess_file_flg, const PetscBool initguesscalc, 
    const char * const initguessfilename, const MPI_Datatype doublefiletype,
    const unsigned char * restrict const W, const PetscInt localfirstrow_l, 
    const PetscInt mg, const PetscInt localfirstcol_l, Vec x, const int rank);

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
    PetscInt * restrict const rcv_across_border_idx);
