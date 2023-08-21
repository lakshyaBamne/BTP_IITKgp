/*
    Namespace to calculate the CU Numerical Flux
*/
#pragma once

#include<iostream>
#include<vector>

#include "Constants.h"
#include "AntiDiffusionCalculator.h"

using namespace std;

namespace ADC=AntiDuffusionCalculator;
namespace CTS=Constants;

/*
    Namespace definition
*/
namespace CUNumericalFlux{
    vector<vector<double>> get_cu_flux(vector<vector<double>>& lsp, vector<vector<double>>& ADTerm, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr, vector<vector<double>>& f1, vector<vector<double>>& f2, vector<vector<double>>& f3);
}

/*
    Function implementation
*/
vector<vector<double>> CUNumericalFlux::get_cu_flux(vector<vector<double>>& lsp, vector<vector<double>>& ADTerm, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr, vector<vector<double>>& f1, vector<vector<double>>& f2, vector<vector<double>>& f3){
    // value of the flux depends on the difference between the Local Speeds of Propagation
    vector<vector<double>> cu_flux;

    int size = lsp[0].size();

    vector<double> cu_f1(size);
    vector<double> cu_f2(size);
    vector<double> cu_f3(size);

    for(int i=0 ; i<size ; i++){
        if( lsp[0][i]-lsp[1][i] > CTS::EPSILON ){
            // a+ - a- appears in the denominator so it cant be too small

            cu_f1[i] = ( lsp[0][i]*f1[0][i] - lsp[1][i]*f1[1][i] + lsp[0][i]*lsp[1][i]*( rho_plr[1][i+1] - rho_plr[0][i] - ADTerm[0][i] ) ) / ( lsp[0][i]-lsp[1][i] );
            cu_f2[i] = ( lsp[0][i]*f2[0][i] - lsp[1][i]*f2[1][i] + lsp[0][i]*lsp[1][i]*( m_plr[1][i+1] - m_plr[0][i] - ADTerm[1][i] ) ) / ( lsp[0][i]-lsp[1][i] );
            cu_f3[i] = ( lsp[0][i]*f3[0][i] - lsp[1][i]*f3[1][i] + lsp[0][i]*lsp[1][i]*( E_plr[1][i+1] - E_plr[0][i] - ADTerm[2][i] ) ) / ( lsp[0][i]-lsp[1][i] );
        }
        else{
            cu_f1[i] = 0.5 * ( f1[0][i] + f1[1][i] );
            cu_f2[i] = 0.5 * ( f2[0][i] + f2[1][i] );
            cu_f3[i] = 0.5 * ( f3[0][i] + f3[1][i] );
        }
    }

    cu_flux.push_back(cu_f1);
    cu_flux.push_back(cu_f2);
    cu_flux.push_back(cu_f3);

    return cu_flux;
}



