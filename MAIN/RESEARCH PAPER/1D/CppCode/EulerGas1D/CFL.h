/*
    CFL - Courant, Friedrics and Levy conditions

    Namespace to calculate dt for the next iteration given the local speeds of propagation
*/
#pragma once

#include<iostream>
#include<vector>

using namespace std;

/*
    Namespace definition
*/
namespace CFL{
    double get_new_dt(vector<vector<double>>& lsp, double T_final, double dx);
}   

/*
    Function implementations
*/
double CFL::get_new_dt(vector<vector<double>>& lsp, double T_final, double dx){
    int size = lsp[0].size();
    
    double dT;

    double a_max=0.0;
    for(int i=0 ; i<size ; i++){
        a_max = max( a_max , max( lsp[0][i], -1*lsp[1][i] ) );
    }

    // Taken from the fortran code directly
    dT = 0.475*dx / a_max;

    return dT;
}
