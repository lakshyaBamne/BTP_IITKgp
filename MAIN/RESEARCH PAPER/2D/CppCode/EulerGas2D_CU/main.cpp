/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#include<iostream>
#include<vector>

#include "InitializeRiemannProblem.h"

using namespace std;

int main(){

    // class RP (RiemannProblem) is responsible for initializing variables and running Numerical Schemes
    RP rp1("CU", "MCW"); // CU
    RP rp2("REF", "MCW"); // Reference


    return 0;
}

