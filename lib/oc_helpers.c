// oc_helpers.c

#include "oc_helpers.h"

// Only for Gauss polynomial of order P = 2
void prepare_edge_matrix_and_rhs(
        int mode, int N, 
        const double *Xp, double *A, double Xrhs[N], 
        double D, double h, const double *zeta)
{
    // mode = 0 : fill rhs only
    // mode = 1 : fill matrix only
    // mode = 2 : fill matrix and rhs

    switch(mode) {
    case 0:
    case 2:
        //
        // Fill right-hand side of tridiagonal system for edge values
        //
        //   The purpose of this step is to guaranteee continuity of the
        //   flux (or concentration) at the edges.

        // Symmetry boundary condition at the center
        Xrhs[i] = 0;

        for (int i = 1; i < N-1; i++) {
            Xrhs[i] = 
        }
        
        // Constant Dirichlet boundary
        // (for time-varying concentration, we need the derivative!)
        Xrhs[N] = 0;

        if (mode == 0) break;
        // [[fallthrough]]
    case 1:
        //
        // Fill diagonals of the tridiagonal matrix system
        //

// Helper macros
#define AL(i)  A[i + N]
#define AC(i)  A[i + N*2]
#define AU(i)  A[i + N*3]

        AC(0) = ; 
        AU(0) = ;
        
        for (int i = 1; i < nel; i++) {
            AC(i) =
            AU(i) =
            AL(i-1) =
        }
        
        AC(N) = ; 
        AL(N-1) = ;

#undef AL
#undef AC
#undef AU

    } // switch

}