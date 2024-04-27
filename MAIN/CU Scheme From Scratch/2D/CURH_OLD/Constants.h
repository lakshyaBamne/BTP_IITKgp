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
#include<string>

using namespace std;

// Namespace definition
namespace Constants{
    const double GAMMA = 1.4;
    const double THETA = 1.3; // Tuning parameter to calculate the Piecewise Linear Reconstruction
    const double EPSILON = 1.0E-12;
    const double PI = 3.14159;
    const double CFL = 0.475;
    const double SMOOTHP = 0.00625;
}