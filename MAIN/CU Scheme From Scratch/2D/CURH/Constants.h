/*
    ! @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-8 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme - with RH optimization for Euler Equations of Gas Dynamics (CURH Scheme)
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

/*
    Header file containing the constant values related to the program
*/
#pragma once

#include<iostream>
#include<string>

using namespace std;

namespace Constants{
    //! Constants for all the other problems
    const double GAMMA = 1.4;
    const double THETA = 1.3;
    const double EPSILON = 1.0E-12;
    const double PI = 3.14159;
    const double CFL = 0.45;
    const double SMOOTHP = 0.00625;

    //! Constants for Raleigh Taylor Instability 
    // const double GAMMA = 5/3;
    // const double THETA = 1.3;
    // const double EPSILON = 1.0E-12;
    // const double CFL = 0.475;
    
    //! Switch for CURH Optimization block
    // const bool ShrinkSw = true;
    // const bool ScaleSw = true;

    const bool ShrinkSw = false;
    const bool ScaleSw = false;
}

