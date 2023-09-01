#include "mex.h"
#include "matrix.h"

#include "dlbm.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    // TODO: figure out we are expecting row or column matrix

    const double *rw = mxGetDoubles(prhs[0]);
    const int m     = mxGetM(prhs[0]);

    const double *X = mxGetDoubles(prhs[1]);

    mxAssert(m == mxGetM(prhs[1]), "Inputs must be of same length!");

    // Perform integration
    double length, res;
    res = dlbm_integrate(&m,rw,X,&length);

    plhs[0] = mxCreateDoubleScalar(length);
    plhs[1] = mxCreateDoubleScalar(res);

}

