#pragma once

#include<iostream>
#include<string>
#include<utility>
#include<unordered_map>
#include<cmath>
#include<vector>

using namespace std;

#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define vvvvd vector< vector< vector< vector<double> > > >
#define umap_ss unordered_map<string,string>
#define pdd pair<double,double>
#define pii pair<int,int>

namespace Constants{

    //! Moving Contact Wave
    int n = 100;
    pdd time = {0,2.0};
    pdd dom = {0,1.0};

    //! Stationary Contact Wave
    // int n = 100;
    // pdd time = {0,0.038};
    // pdd dom = {0,1.0};

    //! Blast Wave Problem
    // int n = 100;
    // pdd time = {0,0.038};
    // pdd dom = {0,1.0};

    //! Shu Osher Problem
    // int n = 100;
    // pdd time = {0,0.038};
    // pdd dom = {0,1.0};

    //! Lax Problem
    // int n = 100;
    // pdd time = {0,0.038};
    // pdd dom = {0,1.0};

    double dx = ( dom.second - dom.first ) / n;

    //! Constants for all the other problems
    const double GAMMA = 1.4;
    const double THETA = 1.3;
    const double EPSILON = 1.0E-12;
    const double CFL = 0.475;

    //! Information on the problem to be solved using the simulations
    const string PROBLEM = "MCW";
    // const string PROBLEM = "SCW";
    // const string PROBLEM = "BLW";
    // const string PROBLEM = "SOS";
    // const string PROBLEM = "LAX";

    const int MODE = 1; // CU
    // const int MODE = 2; // REF

    const string BC = "FREE";
    // const string BC = "SOLIDWALL";
}

