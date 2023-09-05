/*
    Namespace to store constant values related to the program
*/
#pragma once

#include<iostream>

using namespace std;

// Namespace definition
namespace Constants{
    const double GAMMA = 1.4;
    const double THETA = 1.3; // Tuning parameter to calculate the Piecewise Linear Reconstruction
    const double EPSILON = 1.0E-12;

    // CFL number used in calculation of the time steps
    //! replace by 0.1 when using Euler Forward Difference for Time Discretization
    const double CFL = 0.475;
    // const double CFL = 0.1;  

    // function to print all the constants being used in the program
    void print_constants();
}

/*
    Function implementations for the associated functions
*/

void Constants::print_constants(){
    cout << "-----------------------CONSTANTS------------------------" << endl;
    cout << "GAMMA : " << GAMMA << endl;
    cout << "THETA : " << THETA << endl;
    cout << "EPSILON : " << EPSILON << endl;
    cout << "CFL NUMBER : " << CFL << endl;
    cout << "--------------------------------------------------------" << endl;
}
