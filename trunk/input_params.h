PetscErrorCode
input_params (PetscInt * npx, PetscBool * npxflg,
    PetscInt * npy, PetscBool * npyflg, PetscInt * N, PetscInt * M,
    PetscInt * fullN, PetscInt * fullM,
    PetscInt * globalfirstrow_g, PetscInt * globalfirstcol_g,
    char * infilename, char * outfilename, char * truefilename,
    PetscBool * write_outfile_flg, PetscBool * compare_true_flg,
    PetscBool * crossborder, char * initguessfilename,
    PetscBool * init_guess_file_flg, PetscBool * initguesscalc,
    PetscBool * init_guess_flg, PetscBool * init_guess_calc_flg);
