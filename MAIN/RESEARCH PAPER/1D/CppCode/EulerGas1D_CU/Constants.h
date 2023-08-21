/*
    Namespace to store constant values related to the program
*/
#pragma once

#include<iostream>

using namespace std;

// Namespace definition
namespace Constants{
    const double GAMMA = 1.4;
    const double THETA = 1.3;
    const double EPSILON = 1.0E-12;

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
    cout << "--------------------------------------------------------" << endl;
}
