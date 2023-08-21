/*
    Header file to implement calculation of the 
    One Sided Local Speeds of Propagation

    -> CU Scheme : Using highest and lowest Eigen values of the Flux Jacobian
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<cmath>

#include "Constants.h"
#include "PiecewiseConstruction.h"

using namespace std;

namespace C=Constants;
namespace PLR=PiecewiseLinearReconstruction;

// Namespace definition
namespace LocalSpeedsOfPropagation{
    // 1) CU Scheme
    // function to find the Local speeds of propagation
    // using the highest and lowest Eigen values for the Flux Jacobian
    vector<vector<double>> cu_lsp(vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr);

    // Function to print the One sided local speeds of propagation
    void print_lsp(vector<vector<double>>& lsp);
}

/*
    Function implementations
*/

// CU Local speeds of propagation a+ and a-
//! a+ and a- are calculated for all the cells other than the right ghost value cell
vector<vector<double>> LocalSpeedsOfPropagation::cu_lsp(vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr){
    vector<vector<double>> lsp;

    // To apply the formula we need u,p,rho which need to be calculated using vector cv
    vector<vector<double>> u_plr = PLR::primitive_plr(rho_plr, m_plr);
    vector<vector<double>> p_plr = PLR::primitive_plr(u_plr, rho_plr, E_plr);

    int size = rho_plr[0].size()-1;

    vector<double> ap(size); // a+
    vector<double> am(size); // a-

    for(int i=0 ; i<size ; i++){
        double speed_sound_E =  sqrtf( (C::GAMMA*p_plr[0][i])/rho_plr[0][i] );
        double speed_sound_W = sqrtf( (C::GAMMA*p_plr[1][i])/rho_plr[1][i+1] );

        // a+
        ap[i] = max( 0.0 , max( (u_plr[0][i] + speed_sound_E) , (u_plr[1][i] + speed_sound_W) ) );

        // a-
        am[i] = min( 0.0 , min( (u_plr[0][i] - speed_sound_E) , (u_plr[1][i] - speed_sound_W) ) );
    }

    // output
    lsp.push_back(ap);
    lsp.push_back(am);

    return lsp;
}

void LocalSpeedsOfPropagation::print_lsp(vector<vector<double>>& lsp){
    cout << "One sided local speeds of propagation" << endl;
    cout << "a+ : ";
    for(auto i : lsp[0]){
        cout << i << " ";
    }
    cout << endl;

    cout << "a- : ";
    for(auto i : lsp[1]){
        cout << i << " ";
    }
    cout << endl;
}

