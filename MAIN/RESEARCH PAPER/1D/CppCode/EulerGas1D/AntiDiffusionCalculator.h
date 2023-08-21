/*
    Namespace encapsulating functions to Calculate the built in anti diffusion term
*/
#pragma once

#include<iostream>
#include<vector>

#include "Utility.h" // Namespace containing useful utility fucntions
#include "Constants.h" // Namespace containing the constants used in the program

using namespace std;

namespace UTL=Utility;
namespace CTS=Constants;

namespace AntiDuffusionCalculator{
    // Functions to calculate the U* values used in calculation of the Anti-diffusion term
    vector<vector<double>> CalculateStar(vector<vector<double>>& lsp, vector<vector<double>>& f1, vector<vector<double>>& f2, vector<vector<double>>& f3, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr);

    // Function to calculate the Anti-diffusion term
    vector<vector<double>> CalculateADT(vector<vector<double>>& Star, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr);
}

// Function implementations

// Function to Calculate the U* values used for the Anti-diffusion term
vector<vector<double>> AntiDuffusionCalculator::CalculateStar(vector<vector<double>>& lsp, vector<vector<double>>& f1, vector<vector<double>>& f2, vector<vector<double>>& f3, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr ){
    vector<vector<double>> Star;

    int size = lsp[0].size();

    vector<double> rho_star(size);
    vector<double> m_star(size);
    vector<double> E_star(size);

    for(int i=0 ; i<size ; i++){
        rho_star[i] = ( lsp[0][i]*rho_plr[1][i+1] - lsp[1][i]*rho_plr[0][i] - (f1[1][i]-f1[0][i]) ) / ( lsp[1][i]-lsp[0][i] );
        m_star[i] = ( lsp[0][i]*m_plr[1][i+1] - lsp[1][i]*m_plr[0][i] - (f1[1][i]-f1[0][i]) ) / ( lsp[1][i]-lsp[0][i] );
        E_star[i] = ( lsp[0][i]*E_plr[1][i+1] - lsp[1][i]*E_plr[0][i] - (f1[1][i]-f1[0][i]) ) / ( lsp[1][i]-lsp[0][i] );
    }

    Star.push_back(rho_star);
    Star.push_back(m_star);
    Star.push_back(E_star);

    return Star;
}

// Function to calculate the Anti-diffusion term
vector<vector<double>> AntiDuffusionCalculator::CalculateADT(vector<vector<double>>& Star, vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr){
    int size = Star[0].size();

    vector<vector<double>> ADTerm;

    vector<double> d_rho(size);
    vector<double> d_m(size);
    vector<double> d_E(size);

    for(int i=0 ; i<size ; i++){
        d_rho[i] = UTL::minmod( rho_plr[1][i+1]-Star[0][i] , Star[0][i] - rho_plr[0][i] );
        d_m[i] = UTL::minmod( m_plr[1][i+1]-Star[1][i] , Star[1][i] - m_plr[0][i] );
        d_E[i] = UTL::minmod( E_plr[1][i+1]-Star[2][i] , Star[2][i] - E_plr[0][i] );
    }

    ADTerm.push_back(d_rho);
    ADTerm.push_back(d_m);
    ADTerm.push_back(d_E);

    return ADTerm;
}



