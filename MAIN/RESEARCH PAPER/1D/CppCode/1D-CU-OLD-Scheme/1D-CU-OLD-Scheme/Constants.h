#pragma once

#include<iostream>
#include<string>
#include<vector>
#include<unordered_map>
#include<utility>

using namespace std;

#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define pdd pair<double,double>
#define pss pair<string,string>
#define umap_ss unordered_map<string,string>

// Namespace definition
namespace Constants{
    // const string FLUX = "LLF";
    const string FLUX = "HLL";
    // const string FLUX = "CU-OLD";
    // const string FLUX = "CU";
    // const string FLUX = "CURH";

    const double GAMMA = 1.4;
    const double THETA = 1.3; // Tuning parameter to calculate the Piecewise Linear Reconstruction
    const double EPSILON = 1.0E-12;
    const double PI = 3.14159;

    // CFL number used in calculation of the time steps
    // const double CFL = 0.475;
    // const double CFL = 0.45;
    const double CFL = 0.9;
    
    //! replace by 0.1 when using Euler Forward Difference for Time Discretization
    // const double CFL = 0.1;
}

