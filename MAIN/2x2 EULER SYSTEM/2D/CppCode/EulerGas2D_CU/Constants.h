/*
    * @author lakshya Bamne (20MA20029)
    * @supervisor Prof. Naveen Kumar Garg

    ! Implementation of CU scheme for modified Flux functions in Euler's Equations
*/

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
    const double PI = 3.14159;
    const double CFL = 0.475;
    const double SMOOTHP = 0.00625;

    // Constants for Raleigh Taylor Instability 
    const double GAMMA_RT = 5/3;
    const double THETA_RT = 1.3;
    const double EPSILON_RT = 1.0E-12;
    const double CFL_RT = 0.475;

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
    cout << "GAMMA RT : " << GAMMA_RT << endl;
    cout << "THETA RT : " << THETA_RT << endl;
    cout << "EPSILON RT : " << EPSILON_RT << endl;
    cout << "CFL NUMBER RT : " << CFL_RT << endl;
    cout << "--------------------------------------------------------" << endl;
}
