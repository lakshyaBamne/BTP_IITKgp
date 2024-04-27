#pragma once

#include<iostream>

using namespace std;

// Namespace definition
namespace Constants{
    const double GAMMA = 1.4;
    const double THETA = 1.3; // Tuning parameter to calculate the Piecewise Linear Reconstruction
    const double EPSILON = 1.0E-12;
    const double PI = 3.14159;

    // CFL number used in calculation of the time steps
    const double CFL = 0.475;
    
    //! replace by 0.1 when using Euler Forward Difference for Time Discretization
    // const double CFL = 0.1;
}

