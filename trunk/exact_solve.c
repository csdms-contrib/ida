/* Solve for drainage locally using the serial algorithm. Flow from outside
 * the local domain is not included, and so cells that receive such flow
 * will have the wrong value */
#include <petscsys.h>

/* Relative i and j coordinates of eight immediate neighbors */
static const int di[] = { 0, -1, -1, -1, 0, 1, 1, 1 };
static const int dj[] = { 1, 1, 0, -1, -1, -1, 0, 1 };
/* Value neighbor should have if it drains to me */
static const unsigned char w[] = { 16, 8, 4, 2, 1, 128, 64, 32 };

static void getArea (const PetscInt i, const PetscInt j,
    PetscScalar * restrict const A,
    const unsigned char *restrict const W, const PetscInt n,
    const PetscInt m);

/* Recursively determine drainage area (conventional serial calculation) */
  static void
getArea (const PetscInt i, const PetscInt j,
    PetscScalar * restrict const A,
    const unsigned char *restrict const W, const PetscInt n,
    const PetscInt m)
{

  int k;
  PetscInt p, q;

  if (A[i * m + j] != 0)
    /* Area for cell has already been calculated */
    return;

  /* Start with drainage area = 1 (self) */
  A[i * m + j] = 1.0;

  /* Loop over eight neighbors to find ones that drain to me */
  for (k = 0; k < 8; k++)
  {
    /* Calculate neighbor's coordinates */
    p = i + di[k];
    q = j + dj[k];

    /* If neighbor is within domain, and it drains to me, find its
     * drainage area and then add it to my own drainage area value */
    if (p < n && p >= 0 && q < m && q >= 0)
    {
      if (W[p * m + q] == w[k])
      {
        getArea (p, q, A, W, n, m);
        A[i * m + j] += A[p * m + q];
      }
    }
  }

}


/* Loop over all cells in domain and calculate their drainage area */
  void
exact_solve (PetscScalar * restrict const A, 
    const unsigned char * restrict const W, const PetscInt m, const PetscInt n)
{

  PetscInt i, j;

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
    {
      getArea (i, j, A, W, n, m);
    }
  }

}
