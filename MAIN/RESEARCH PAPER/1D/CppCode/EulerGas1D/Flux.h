/*
    Namespace contains function to calculate the Flux vector
    given the Conserved variables
*/
#pragma once

#include<iostream>
#include<vector>

#include "PiecewiseConstruction.h"

using namespace std;

namespace PLR=PiecewiseLinearReconstruction;

// namespace definition
namespace Flux{
    // Function to calculate Flux vector elements
    vector<vector<double>> get_flux(int num, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr);
}

// Function implementations for the namespace
vector<vector<double>> Flux::get_flux(int num, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr){
    vector<vector<double>> flux;
    
    vector<vector<double>> u_plr = PLR::primitive_plr(rho_plr, m_plr);
    vector<vector<double>> p_plr = PLR::primitive_plr(u_plr, rho_plr, E_plr);

    int size = u_plr[0].size();

    vector<double> fE(size);    
    vector<double> fW(size);    

    switch (num){
        //! f1 = rho*u
        case 1:
            for(int i=0 ; i<size ; i++){
                fE[i] = rho_plr[0][i]*u_plr[0][i];
                fW[i] = rho_plr[1][i+1]*u_plr[0][i];
            }

            break;
        //! f2 = m*u + p
        case 2:
            for(int i=0 ; i<size ; i++){
                fE[i] = m_plr[0][i]*u_plr[0][i] + p_plr[0][i];
                fW[i] = m_plr[1][i+1]*u_plr[1][i] + p_plr[1][i];
            }

            break;
        //! f3 = u*(E + p)
        case 3:
            for(int i=0 ; i<size ; i++){
                fE[i] = u_plr[0][i]*( E_plr[0][i] + p_plr[0][i] );
                fW[i] = u_plr[1][i]*( E_plr[1][i+1] + p_plr[1][i] );
            }

            break;
    }

    flux.push_back(fE);
    flux.push_back(fW);

    return flux;

}


