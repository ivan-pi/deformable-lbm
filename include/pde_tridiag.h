// pde_tridiag.h
/* Helper routines for tridiagonal matrix systems */

// Account for compiler name mangling
#if !defined(_WIN32)
#define dgtsv  dgtsv_
#define dgttrf dgttrf_
#define dgttrs dgttrs_
#endif


static const int NRHS = 1;

/* Solve tridiagonal system of equations */
static int factor_and_solve_tridiag(
        mwSignedIndex N, 
        double *A, 
        double *b) 
{
    int ierr;

    dgtsv(&N,&NRHS,&A[0],&A[N],&A[2*N],&A[3*N],b,&N,&ierr);
    mxAssert(ierr >= 0, "Illegal argument in call to dgtsv!");

    return ierr; // 0 if success
}

/* Factorize tridiagonal system using Gaussian elimination */
static int factor_tridiag(
        mwSignedIndex N,
        const double *A,
        const mwSignedIndex *ipiv)
{
    int ierr;

    // Solve system of equations for edge values
    dgttrf(
        &N, &NRHS,
        &A[  0], /* dl  */
        &A[  N], /* d   */
        &A[2*N], /* du  */
        &A[3*N], /* du2 */
        ipiv,
        b, &N, 
        &ierr);

    mxAssert(ierr >= 0, "Illegal argument in call to dgttrs!");

    return ierr; // 0 if success
}

/* Solve system of tridiagonal equations using factorized matrix */
static void solve_tridiag(
        mwSignedIndex N,
        const double *A,
        const mwSignedIndex *ipiv,
        double *b)
{
    const char *trans = "N";
    int ierr;

    // Solve system of equations for edge values
    dgttrs(trans,
        &N, &NRHS,
        &A[  0], /* dl  */
        &A[  N], /* d   */
        &A[2*N], /* du  */
        &A[3*N], /* du2 */
        ipiv,
        b, &N, 
        &ierr);

    mxAssert(ierr == 0, "Illegal argument in call to dgttrs!");
}