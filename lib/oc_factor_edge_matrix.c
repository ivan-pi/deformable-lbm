#include "mex.h"
#include "matrix.h"

// Details about calling LAPACK and BLAS functions in MATLAB are found in
// https://www.mathworks.com/help/matlab/matlab_external/calling-lapack-and-blas-functions-from-mex-files.html

// Helpers for edge matrix
#include "oc_helpers.h"

// Helpers for tridiagonal matrix factorization (requires LAPACK)
#include "pde_tridiag.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Value of independent variable
    const double t     = mxGetScalar(prhs[0]);

    // Get pointer to dependent variables
    const double *y    = mxGetDoubles(prhs[1]);
    const int m        = mxGetM(prhs[1]); /* number of rows in y */

    // X .... 3*nel + 1
    // r ....   nel + 1
    // ----------------
    //        4*nel + 2

    mxAssert((m - 2) % 4 == 0, 
        "Number of equations m does not match expectations of the routine!");

    // Number of edges
    const int N = nel + 1;

    // Unpack other input Parameters
    const double h      = mxGetScalar(prhs[2]);
    const double *zeta  = mxGetDoubles(prhs[3]);
    const double D      = mxGetScalar(prhs[4]);
    const double alpha  = mxGetScalar(prhs[5]);

    // Create output arrays 
    plhs[0] = mxCreateDoubleMatrix(m, 4, mxREAL);
    plhs[1] = mxCreateNumericMatrix(m, 1, mxINT32_CLASS, mxREAL);

    // Unpack pointers to buffers
    double *A = mxGetDoubles(plhs[0]);
    int *ipiv = mxGetInt32s(plhs[1]);

    //
    // ---- Function logic ----
    //

    // Prepare tridiagonal matrix
    {
        double *Xpdum, *bdum; // unused
        const int mode = 1; // rhs only
        prepare_edge_matrix_and_rhs(mode,N,Xpdum,A,bdum,D,h,zeta);
    }

    // The matrix A has been filled. We can proceed to factorize it.
    const int ierr = factor_tridiag(N,A,ipiv);
    mxAssert(ierr == 0, "The factorization of the tridiagonal matrix failed!");

}