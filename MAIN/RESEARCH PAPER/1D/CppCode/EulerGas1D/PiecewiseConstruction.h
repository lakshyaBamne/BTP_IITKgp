/*
    Namespace to encapsulate the functions to find the Piecewise Linear Reconstruction
    
    1) Calculate the slopes using the minmod function
    2) also extend the cells for ghost values
    3) Now calculate the Respective piecewise linear reconstruction inside each cell
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>

#include "Constants.h"
#include "Utility.h"

using namespace std;

namespace CTS=Constants;
namespace UTL=Utility;

namespace PiecewiseLinearReconstruction{
    // function to construct the Piecewise Linear Reconstruction of a Conserved variable
    // in 1D i.e., returns the East and West Extended vectors
    vector<vector<double>> construct_plr(vector<double>& var);

    // function which finds the PLR of primitive variables like velocity(u) and pressure(p)
    // given the PLR for the conserved variables (using polymorphism for different primitive variables)
   
    // -> for calculating velocity plr
    vector<vector<double>> primitive_plr(vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr);

    // -> for calculating pressure plr
    vector<vector<double>> primitive_plr(vector<vector<double>>& u_plr ,vector<vector<double>>& rho_plr, vector<vector<double>>& E_plr);

    // Function to print the Piecewise Linear Reconstruction for all the primitive variables
    void print_plr_conserved(vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr);
}

/*
    Function implementations
*/

// Function to find the Piecewise Linear Reconstruction for the Conserved variables
vector<vector<double>> PiecewiseLinearReconstruction::construct_plr(vector<double>& var){

    vector<vector<double>> plr; // 2D vector contains the East and West PLR

    vector<double> East;
    vector<double> West;

    for(int i=1 ; i<var.size()-1 ; i++){
        auto v1 = CTS::THETA*( var[i+1] - var[i] );
        auto v2 = ( var[i+1] - var[i-1] )/2;
        auto v3 = CTS::THETA*( var[i] - var[i-1] );

        auto slope_val = UTL::minmod(v1,v2,v3);

        East.push_back(var[i] + slope_val/2);        
        West.push_back(var[i] - slope_val/2);
    }

    // Extend the cells for ghost values for the plr vectors (East and West)
    //! here we are extending using the FREE boundary conditions which will be changed later
    //! and should be given an option for the user to choose the Boundary Conditions for a problem

    East.insert(East.begin(), East[0]);
    East.push_back(East.back());

    West.insert(West.begin(), West[0]);
    West.push_back(West.back());

    // output is returned as a 2D vector containing East, West in indices 0,1 respectively
    plr.push_back(East);
    plr.push_back(West);

    return plr;

}

// Function to find the Piecewise Linear Reconstruction for velocity
vector<vector<double>> PiecewiseLinearReconstruction::primitive_plr(vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr){
    vector<vector<double>> u_plr;
    
    int size = rho_plr[0].size()-1;

    vector<double> u_E(size);
    vector<double> u_W(size);

    for(int i=0 ; i<size ; i++){
        u_E[i] = m_plr[0][i] / rho_plr[0][i];
        u_W[i] = m_plr[1][i+1] / rho_plr[1][i];
    }

    // Ghost values are not required for reconstruction of Primitive variables

    u_plr.push_back(u_E);
    u_plr.push_back(u_W);

    return u_plr;
}


// Function to find the Piecewise Linear Reconstruction for Pressure
vector<vector<double>> PiecewiseLinearReconstruction::primitive_plr(vector<vector<double>>& u_plr, vector<vector<double>>& rho_plr, vector<vector<double>>& E_plr){
    vector<vector<double>> p_plr;

    int size = u_plr[0].size()-1;

    vector<double> p_E(size);
    vector<double> p_W(size);

    for(int i=0 ; i<size ; i++){
        p_E[i] = (CTS::GAMMA-1)*( E_plr[0][i] - ( (rho_plr[0][i]*u_plr[0][i]*u_plr[0][i] )/2 ) );
        p_W[i] = (CTS::GAMMA-1)*( E_plr[1][i+1] - ( (rho_plr[1][i+1]*u_plr[1][i]*u_plr[1][i] )/2 ) );
    }

    // Ghost values are not required for reconstruction of Primitive variables

    p_plr.push_back(p_E);
    p_plr.push_back(p_W);

    return p_plr;
}

// Function to print the PLR for all the conserved variables
void PiecewiseLinearReconstruction::print_plr_conserved(vector<vector<double>>& rho_plr, vector<vector<double>>& m_plr, vector<vector<double>>& E_plr){
    cout << "Piecewise Linear Reconstruction for conserved variables (len=" << rho_plr[0].size() << ")" << endl;
    
    cout << "rho_East : ";
    for(auto i : rho_plr[0]){
        cout << i << " ";
    }
    cout << endl;

    cout << "rho_West : ";
    for(auto i : rho_plr[1]){
        cout << i << " ";
    }
    cout << endl;

    cout << "m_East : ";
    for(auto i : m_plr[0]){
        cout << i << " ";
    }
    cout << endl;

    cout << "m_West : ";
    for(auto i : m_plr[1]){
        cout << i << " ";
    }
    cout << endl;

    cout << "E_East : ";
    for(auto i : E_plr[0]){
        cout << i << " ";
    }
    cout << endl;

    cout << "E_West : ";
    for(auto i : E_plr[1]){
        cout << i << " ";
    }
    cout << endl;
}

