PetscErrorCode formA (Mat A, const PetscInt M, const PetscInt N,
    const PetscInt m, const PetscInt n,
    const PetscInt mg,
    const PetscInt localfirstcol_l,
    const PetscInt localfirstrow_l,
    const PetscInt localfirstcol_g,
    const PetscInt localfirstrow_g,
    const unsigned char * restrict const W);
PetscErrorCode
formAr (Mat Ar, const PetscInt M, const PetscInt N,
    const PetscInt m, const PetscInt n,
    const PetscInt localfirstcol_g, const PetscInt localfirstrow_g,
    const unsigned char * restrict const W, const AO ao,
    const PetscInt * restrict const aoordering, const PetscInt num_ao,
    PetscScalar * restrict const br_vals, 
    PetscInt * restrict const br_vals_idx,
    PetscInt * restrict const num_br_vals, Vec x);
