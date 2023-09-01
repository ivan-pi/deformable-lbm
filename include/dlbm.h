
// Initialize lattice
void dlbm_init(
    int n,
    double csqr, double dt, double X0, double L0,
    double *X, double *us, double *r, double *Npdf);

// Main time-stepping
void dlbm_step(
    int nsteps,
    double D,
    double csqr,
    int n,
    double *X,
    double *us,
    double *r,
    double *N);

// Midpoint rule for integration
double dlbm_integrate(
    int *n,
    const double *r
    const double *X
    double *radius);