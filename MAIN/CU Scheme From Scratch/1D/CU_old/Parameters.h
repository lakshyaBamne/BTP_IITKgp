#pragma once

#include<iostream>
#include<string>

using namespace std;

namespace Parameters{
    /*
        Select the problem to be solved
    */
    // const int PROBLEM = 1; // Moving Contact Wave (MCW)
    // const int PROBLEM = 2; // Stationary Contact Wave (SCW)
    // const int PROBLEM = 3; // Lax Problem (LAX)
    // const int PROBLEM = 4; // Shu Osher Problem (SOS)
    const int PROBLEM = 5; // Blast Wave Problem (BLW)

    /*
        Select appropriate boundary conditons
    */
    // const int BC = 1; // FREE
    const int BC = 2; // SOLIDWALL

    /*
        Select the appropriate simulation type
    */
    const int SIMTYPE = 1; // PLOT
    // const int SIMTYPE = 2; // ANIMATE

    /*
        Select the appropriate mode to run the simulation
    */
    // const int MODE = 1; // CU
    const int MODE = 2; // CU-REF
    // const int MODE = 3; // CURH
    // const int MODE = 4; // CURH-REF

}



