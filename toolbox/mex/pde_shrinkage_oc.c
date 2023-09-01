// pde_shrinkage.c

#include "mex.h"
#include "matrix.h"

// Details about calling LAPACK and BLAS functions in MATLAB are found in
// https://www.mathworks.com/help/matlab/matlab_external/calling-lapack-and-blas-functions-from-mex-files.html
//
// * Use -lmwlapack -lmwblas for linking
// * Only 64-bit integers are supported for matrix dimensions

// Helpers for tridiagonal matrix factorization
#include "pde_tridiag.h"

// Helpers for edge matrix
#include "oc_helpers.h"


/* Evaluates the semi-discrete PDE for diffusive drying with shrinkage
 *
 *  X[3*nel + 1]
 *  r[  nel + 1]
 * Xp[3*nel + 1]
 * rp[  nel + 1]
 * zeta[nel + 1]
 */
static void pde_shrinkage_oc(
        int nel, 
        double t,
        const double *X, 
        const double *r 
        double *Xp, 
        double *rp,
        double h,
        const double *zeta,
        double D,
        double alpha
        double *A /* tridiagonal matrix */
        int *ipiv /* pivot vector */) 
{

// Generally, D will be a function of X and T, D := D(X,T)
#define Dz(X) (D/(1.0 + alpha*X))
#define cntr(x1,x2) (0.5*(x1 + x2))

    //
    // Evaluate moisture change rate
    //

    // TODO: VERIFY THIS EQUATION!
    // Concentration at center node, r = xi = 0
    Xp[0] = 4 * Dz(X[0]) * (X[1] - X[0]) / (h*h);
    
    // Bulk residuals
    for (int i = 1; i < N-1; i++) {
        Xp[i] = 
    }

    // Dirichlet boundary at the surface
    Xp[N] = 0; 

    // Solve A x = b for the edge values
#if VARDIFF
    // For variable diffusivity, we need to assemble and factorize the
    // matrix at each time step:

    const int mode = 2; // matrix and rhs
    prepare_edge_matrix_and_rhs(mode,N,Xp,A,Xrhs,D,h,zeta);
    const int ierr = factor_and_solve_tridiag(N,A,Xrhs);
    mxAssert(ierr == 0, "The factorization of the tridiagonal matrix failed!")
#else
    // For constant diffusivity, we only need to perform back substitution.
    // The matrix A should be factorized beforehand. 
    // We also need the pivot vector ipiv.

    const int mode = 0; // rhs only
    prepare_edge_matrix_and_rhs(mode,N,Xp,A,Xrhs,D,h,zeta);
    solve_tridiag(N,A,ipiv,Xrhs);
#endif

    // Copy results into output array
#if WITH_BLAS
    const mwSignedIndex incrhs = 1, incp = 3;
    dcopy_(&N,Xrhs,&incrhs,Xp,&incp);
#else
    for (int i = 0; i < N; i++) {
        Xp[3*i] = Xrhs[i];
    }
#endif
    
    //
    // Displacement equations (only needed for plotting)
    //

    // Radial displacement at the center
    rp[0] = 0;

    // Displacement in the bulk (we simply use a central difference)
    for (int i = 1; i < N-1; i++) {
        const double dXdzeta = (X[i+1] - X[i-1])/(2*h);
        rp[i] = alpha * Dz(X[i]) * zeta[i] * dXdzeta / r[i];
    }

    // Surface displacement (one-sided difference)
    const double dXdzeta = (0.5*X[N-3] - 2.0*X[N-2] + 1.5*X[N-1])/h;
    rp[N-1] = alpha * Dz(X[N-1]) * zeta[N-1] * dXdzeta / r[N-1];

#undef cntr
#undef Dz

} // pde_shrinkage

// TODO: write the Jacobian sparsity pattern. This will be a nightmare!
//       The edge values are all linked.
//       Each edge value depends on the closest two residuals

//      dXdt/dX = ... nightmare ...
//      dXdt/dr = 0
//
//      drdt/dr = diag(ones)
//
//      drdt/dX = diag(ones,-1:1)
//      drdt/dX|N = ...

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

    // Number of elements
    const int nel = (m - 2) / 4;

    // Number of DOF's in moisture vector
    const int ndof = 3*nel + 1;

    // Unpack other input Parameters
    const double h      = mxGetScalar(prhs[2]);
    const double *zeta  = mxGetDoubles(prhs[3]);
    const double D      = mxGetScalar(prhs[4]);
    const double alpha  = mxGetScalar(prhs[5]);

    // Unpack workspace for tridiagonal matrix system
    const double *w = mxGetDoubles(prhs[6]);
    const int *ipiv = mxGetInt32(prhs[7]);

    // Create output buffer
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);

    // Get pointer to storage for right-hand side, 
    double *dydt = mxGetDoubles(plhs[0]);

    // Evaluate dy/dt = F(t,y)
    pde_shrinkage_oc(nel, t, 
           &y[0],    &y[ndof],   /* inputs */
        &dydt[0], &dydt[ndof],  /* outputs */
        h, zeta, D, alpha,
        wl,wc,wu,wu2,ipiv);
}