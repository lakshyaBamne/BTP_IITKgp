/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

/*
    Namespace to store constant values related to the program
*/
#pragma once

#include<iostream>

using namespace std;

namespace Constants{
    const double GAMMA = 1.4;
    const double THETA = 1.3; // Tuning parameter to calculate the Piecewise Linear Reconstruction
    const double EPSILON = 1.0E-12;
    const double PI = 3.14159;
    const double CFL = 0.475;
    const double SMOOTHP = 0.00625;

    // function to print all the constants being used in the program
    void print_constants();
}

// Constants for Raleigh Taylor Instability 
namespace ConstantsRT{
    const double GAMMA = 5.0/3.0;
    const double THETA = 1.3;
    const double EPSILON = 1.0E-12;
    const double PI = 3.14159;
    const double CFL = 0.475;
    const double SMOOTHP = 0.00625;

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
    cout << "PI : " << PI << endl;
    cout << "SMOOTHING PARAMETER : " << SMOOTHP << endl;
}

void ConstantsRT::print_constants(){
    cout << "----------------------CONSTANTS-RT----------------------" << endl;
    cout << "GAMMA RT : " << GAMMA << endl;
    cout << "THETA RT : " << THETA << endl;
    cout << "EPSILON RT : " << EPSILON << endl;
    cout << "CFL NUMBER RT : " << CFL << endl;
    cout << "PI RT : " << CFL << endl;
    cout << "SMOOTHING PARAMETER RT : " << CFL << endl;
    cout << "--------------------------------------------------------" << endl;
}


