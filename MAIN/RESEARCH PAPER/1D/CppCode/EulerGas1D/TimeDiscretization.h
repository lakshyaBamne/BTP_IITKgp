/*
    Namespace to calculate the new values of the Conserved variables
    using the 3-stage SSP-RK scheme
    ! SSP-RK : Strong Stability preserving Runge Kutta Methods
*/
#pragma once

#include<iostream>
#include<vector>

using namespace std;

namespace ThreeStageSSPRK{
    // Function to get the new values for the Conserved variables using 3-stage SSPRK
    vector<vector<double>> cu_update(vector<vector<double>>& old, vector<vector<double>>& cu_flux, double dx, double dt);
}

// Function implementations
vector<vector<double>> ThreeStageSSPRK::cu_update(vector<vector<double>>& old, vector<vector<double>>& cu_flux, double dx, double dt){
    int size = old[0].size()-2;
    
    double LAMBDA = dt/dx;

    vector<vector<double>> step1; // U1
    vector<vector<double>> step2; // U2
    vector<vector<double>> step_new; // Unew

    // STEP-1 U1 calculation
    vector<double> step11(size);
    vector<double> step12(size);
    vector<double> step13(size);

    for(int i=0 ; i<size ; i++){
        step11[i] = old[0][i+1] - LAMBDA*( cu_flux[0][i+1] - cu_flux[0][i] );
        step12[i] = old[1][i+1] - LAMBDA*( cu_flux[1][i+1] - cu_flux[1][i] );
        step13[i] = old[2][i+1] - LAMBDA*( cu_flux[2][i+1] - cu_flux[2][i] );
    }

    // extend these cells using the given conditions (here FREE)
    step11.insert(step11.begin(), step11[0]);
    step11.push_back(step11.back());

    step12.insert(step12.begin(), step12[0]);
    step12.push_back(step12.back());

    step13.insert(step13.begin(), step13[0]);
    step13.push_back(step13.back());

    // STEP-2 U2 calculation -> to calculate this we need the CU Numerical Flux for the U1 vector
    vector<double> step21(size);
    vector<double> step21(size);
    vector<double> step21(size);
    
    // STEP-3 U(n+1) calculation -> to calculate this we need the CU Numerical Flux for the U2 vector

}
