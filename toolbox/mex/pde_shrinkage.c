// pde_shrinkage.c

#include "mex.h"
#include "matrix.h"

/* Evaluates the semi-discrete PDE for diffusive drying with shrinkage
 */
static void pde_shrinkage(
        int N, 
        double t,
        const double X[N], 
        const double r[N], 
        double Xp[N], 
        double rp[N],
        double h,
        const double zeta[N],
        double D,
        double alpha) {

#define Dzeta(X) (D/(1.0 + alpha*X))
#define cntr(x1,x2) (0.5*(x1 + x2))

    //
    // Evaluate moisture change rate
    //

    // Concentration at center node, r = xi = 0
    Xp[0] = 4*Dzeta(X[0])*(X[1] - X[0])/(h*h);
    
    // Radial displacement at the center
    rp[0] = 0;

    // Bulk nodes, form differences
    for (int i = 1; i < N-1; i++) {
        
        const double zr = cntr(zeta[i+1],zeta[i]);
        const double zl = cntr(zeta[i],zeta[i-1]);

        const double fr = Dzeta(cntr(X[i+1],X[i])) * zr * (X[i+1] - X[i])/h;
        const double fl = Dzeta(cntr(X[i],X[i-1])) * zl * (X[i] - X[i-1])/h;

        // Moisture change rate
        const double gradient_of_flux = (fr - fl)/h;
        Xp[i] = gradient_of_flux/zeta[i];

        // Displacement in the bulk
        const double dXdzeta = (X[i+1] - X[i-1])/(2*h);
        rp[i] = alpha * Dzeta(X[i]) * zeta[i] * dXdzeta / r[i];
    }

    // Dirichlet boundary at surface
    Xp[N-1] = 0;

    // Surface displacement
    const double dXdzeta = (0.5*X[N-3] - 2.0*X[N-2] + 1.5*X[N-1])/h;
    rp[N-1] = alpha * Dzeta(X[N-1]) * zeta[N-1] * dXdzeta / r[N-1];

#undef cntr
#undef Dzeta

} // pde_shrinkage

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Value of independent variable
    const double t     = mxGetScalar(prhs[0]);

    // Get pointer to dependent variables
    const double *y    = mxGetDoubles(prhs[1]);
    const int m        = mxGetM(prhs[1]); /* number of rows in y */

    mxAssert(m % 2 == 0, 
        "Number of equations m is not divisible by 2 as expected!");

    // Unpack other input Parameters
    const double h     = mxGetScalar(prhs[2]);
    const double *zeta  = mxGetDoubles(prhs[3]);
    const double D     = mxGetScalar(prhs[4]);
    const double alpha = mxGetScalar(prhs[5]);

    // Create output buffer
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);

    // Get pointer to storage for right-hand side, 
    double *dydt = mxGetDoubles(plhs[0]);

    // Number of nodes or variables per discrete field
    const int N = m / 2;

    // Evaluate dy/dt = F(t,y)
    pde_shrinkage(N, t, 
           &y[0],    &y[N],      /* inputs */
        &dydt[0], &dydt[N],  /* outputs */
        h, zeta, D, alpha);
}