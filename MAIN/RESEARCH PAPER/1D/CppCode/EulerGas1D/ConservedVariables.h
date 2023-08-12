#pragma once

#include<iostream>
#include<vector>

#define ll long long int

#include "Domain.h"

using namespace std;

/*
    Class to initialize and encapsulate the vectors of Conserved Variables
    -> Density
    -> Momentum
    -> Energy
*/
class ConservedVariables{
private:
    vector<vector<double>> conserved_variables;

    vector<double> density;
    vector<double> momentum;
    vector<double> energy;

public:
    // constructor takes the domain object and initializes conserved variables using that
    ConservedVariables(Domain* domain);

    // function to show the vectors storing the Conserved Variables
    void show_conserved_variables();

};
