#include "Solver.h"
#include "Solver_HLLC.h"
#include "Parallel.h"

int main() {
    // Parameters for initializing the simulation:
    //   1. Initial time
    //   2. End time
    //   3. Number of grid points in the x-direction
    //   4. Number of grid points in the y-direction
    //   5. CFL number (Courant–Friedrichs–Lewy condition number)
    //   6. Starting point in the x-direction
    //   7. Ending point in the x-direction
    //   8. Starting point in the y-direction
    //   9. Ending point in the y-direction
    //  10. Boundary condition code
    //  11. Enable non-physical point detection (0: No, 1: Yes)
    //  12. Enable timing execution (0: No, 1: Yes)
    //
    // Boundary condition codes:
    //   0: Default boundary condition
    //   1: Periodic boundary condition
    //   2: Reflective boundary condition
    //   3: Solid wall boundary condition

    InitState b(
        0,      // Initial time
        0.25,   // End time
        100,    // Number of grid points in x-direction
        100,    // Number of grid points in y-direction
        0.8,    // CFL number
        0,      // Starting point in x-direction
        2,      // Ending point in x-direction
        0,      // Starting point in y-direction
        2,      // Ending point in y-direction
        0,      // Boundary condition code
        0,      // Non-physical point detection disabled
        1       // Timing execution enabled
    );

    // Parallel computing parameters:
    //   First parameter: Whether parallel is required
    //   Second parameter: Number of processes
    Parallel parallel(
        1,  // Whether parallel is required (0: No, 1: Yes)
        4   // Number of processes
    );

    // To use the SLIC Solver, uncomment the following lines:
    
    Solver2D a(b, parallel);
    a.solve(1); // 1:The progress bar is displayed 0:No progress bar is displayed
    

    // Using the HLLC Solver
    // Solver_HLLC_2D a(b, parallel);
    // a.solve(1); // 1:The progress bar is displayed 0:No progress bar is displayed
    
    return 0;
}