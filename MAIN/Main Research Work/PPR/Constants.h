/*
    ! Namespace contains parameters used in the program for simulations
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>

using namespace std;

#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define pdd pair<double,double>
#define pss pair<string,string>
#define umap_ss unordered_map<string,string>

namespace Constants{
    //! Type of approximation for the simulations
    // string APPROX = "FIRST";
    string APPROX = "SECOND";

    //! Flux to be used
    string FLUX = "LLF";
    // string FLUX = "HLL";
    // string FLUX = "CU-OLD";
    // string FLUX = "CU";
    // string FLUX = "CURH";
    
    //! problem to solve
    string PROBLEM = "MCW";
    // string PROBLEM = "SCW";
    // string PROBLEM = "BLW";
    // string PROBLEM = "TORO-1";
    // string PROBLEM = "TORO-2";
    // string PROBLEM = "TORO-3";
    // string PROBLEM = "TORO-4";

    //! boundary conditions
    string BC = "FREE";
    // string BC = "SOLIDWALL";

    //! Time Integration Scheme
    // string TI = "EFD";
    string TI = "SSPRK";

    //! simulation parameters
    const double GAMMA = 1.4;
    const double THETA = 1.3;
    // const double CFL = 0.9;
    const double CFL = 0.45;

    //! Simulation related parameters
    
    //! Moving Contact Wave (MCW)
    int n = 100;
    pdd domx = {0.0, 1.0};
    pdd time = {0.0, 2.0};

    //! Stationary Contact Wave (SCW)
    // int n = 100;
    // pdd domx = {0.0, 1.0};
    // pdd time = {0.0, 0.012};

    //! Blastwave Problem (BLW)
    // int n = 500;
    // pdd domx = {0.0, 1.0};
    // pdd time = {0.0, 0.038};

    //! Shock Density Wave Problem (SDW)
    // int n = 500;
    // pdd domx = {0.0, 1.0};
    // pdd time = {0.0, 0.038};

    //! Shock Entropy Wave Problem (SEW)
    // int n = 500;
    // pdd domx = {0.0, 1.0};
    // pdd time = {0.0, 0.038};

    //! Toro Test-1
    // int n = 100;
    // pdd domx = {0.0, 1.0};
    // pdd time = {0.0, 0.2};

    //! Toro Test-2
    // int n = 100;
    // pdd domx = {0.0, 1.0};
    // pdd time = {0.0, 0.15};

    //! Toro Test-3
    // int n = 100;
    // pdd domx = {0.0, 1.0};
    // pdd time = {0.0, 0.012};

    //! Toro Test-4
    // int n = 100;
    // pdd domx = {0.0, 1.0};
    // pdd time = {0.0, 0.035};

    int nref = 4*n;
    double dx = ( domx.second - domx.first ) / n;

}

