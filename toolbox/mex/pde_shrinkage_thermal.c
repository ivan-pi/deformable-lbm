// pde_shrinkage.c

#include "mex.h"
#include "matrix.h"

/* Evaluates the semi-discrete PDE for diffusive drying with shrinkage
 */
static void pde_shrinkage_thermal(
        int N, 
        double time,
        double T,               /* temperature */
        const double X[N],      /* moisture content */
        const double r[N],      /* radial position */
        double *Tp,             /* temperature change rate */
        double Xp[N],           /* moisture change rate */
        double rp[N],           /* radial change rate */
        double h,               /* grid spacing */
        const double zeta[N],   /* grid in solid coordinates */
        double D,               /* diffusivity */
        double alpha)           /* shrinkage coefficent */
{

#define Dzeta(X) (D/(1 + alpha*X))
#define cntr(x1,x2) (0.5*(x1 + x2))

    //
    // Evaluate moisture change rate
    //

    // Concentration at center node, r = xi = 0
    Xp[0] = 4*Dzeta(X[0])*(X[1] - X[0])/(h*h);

    // Bulk nodes, form differences
    for (int i = 1; i < N-1; i++) {
        
        const double zr = cntr(zeta[i+1],zeta[i]);
        const double zl = cntr(zeta[i],zeta[i-1]);

        const double fr = Dzeta(cntr(X[i+1],X[i])) * zr * (X[i+1] - X[i])/h;
        const double fl = Dzeta(cntr(X[i],X[i-1])) * zl * (X[i] - X[i-1])/h;

        const double gradient_of_flux = (fr - fl)/h;

        Xp[i] = gradient_of_flux/zeta[i];
    }

    // Robyn boundary at the surface
    const double surface_grad = (0.5*X[N-3] - 2*X[N-2] + 1.5*X[N-1])/(h*h);
    Xp[N-1] = ;

    //
    // Evaluate temperature change
    //
    *Tp = 

    //
    // Drying evolution is done, now we evaluate the displacements ...
    //

    // Radial displacement at the center
    rp[0] = 0;

    // Displacements in the bulk
    for (int i = 1; i < N-1; i++) {
        const double dXdzeta = (X[i+1] - X[i-1])/(2*h);
        rp[i] = alpha * Dzeta(X[i]) * zeta[i] * dXdzeta / r[i];
    }

    // Surface displacement
    const double dXdzeta = (0.5*X[N-3] - 2.0*X[N-2] + 1.5*X[N-1])/h;
    rp[N-1] = alpha * Dzeta(X[N-1]) * zeta[N-1] * dXdzeta / r[N-1];

#undef cntr
#undef Dzeta

} // pde_shrinkage

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Value of independent variable
    const double *t     = mxGetDoubles(prhs[0]);


    // Get pointer to dependent variables
    const double *y     = mxGetDoubles(prhs[1]);

    // Get number of rows in column vector y
    const int m = mxGetM(prhs[1]);

    // Unpack other input Parameters
    const double *h     = mxGetDoubles(prhs[2]);
    const double *zeta  = mxGetDoubles(prhs[3]);
    const double *D     = mxGetDoubles(prhs[4]);
    const double *alpha = mxGetDoubles(prhs[5]);

    // Create output buffer
    plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL);

    // Get pointer to storage for right-hand side, dy/dt = F(t,y)
    double *dydt = mxGetDoubles(plhs[0]);

    mxAssert(0 == m % 2, "Number of equations m is not divisible by 2!");
    const int N = m/2;

    pde_shrinkage(N,
        *t,&y[0],&y[N],&y[N+1],
        &dydt[0],&dydt[N],&dydt[N+1],
        *h, zeta, *D, *alpha);

}