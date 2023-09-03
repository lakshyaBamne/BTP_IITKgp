/*
    Namespace to calculate the next iteration using simple Euler Forward Difference
*/
#pragma once

#include<iostream>
#include<vector>

#include "ExtendCells.h"
#include "CUNumericalFlux.h"

using namespace std;

namespace CUF = CUNumericalFlux;
namespace EXC = ExtendCells;

namespace EulerForwardDifference{
    vector<vector<double>> get_next_cons_vars(vector<vector<double>>& cons_vars, vector<vector<double>>& cu_flux, double LMBDA, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time);
}

vector<vector<double>> EulerForwardDifference::get_next_cons_vars(vector<vector<double>>& cons_vars, vector<vector<double>>& cu_flux, double LMBDA, pair<string,string>& initial_conditions, double& dt, double& dx, double& t, pair<double,double>& time){

    // find the flux for the present conservative vector
    // value of dt should not change in this step

    int size = cons_vars[0].size()-2;

    vector<vector<double>> new_cons;

    vector<double> u1(size);
    vector<double> u2(size);
    vector<double> u3(size);

    for(int i=1 ; i<=size ; i++){
        u1[i-1] = cons_vars[0][i] - LMBDA*( cu_flux[0][i] - cu_flux[0][i-1] );
        u2[i-1] = cons_vars[1][i] - LMBDA*( cu_flux[1][i] - cu_flux[1][i-1] );
        u3[i-1] = cons_vars[2][i] - LMBDA*( cu_flux[2][i] - cu_flux[2][i-1] );
    }

    new_cons.push_back(u1);
    new_cons.push_back(u2);
    new_cons.push_back(u3);

    EXC::extend_cells(initial_conditions.second , new_cons);

    return new_cons;
}



