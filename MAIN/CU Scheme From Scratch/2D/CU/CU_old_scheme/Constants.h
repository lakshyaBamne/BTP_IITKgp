#pragma once

#include<iostream>
#include<string>

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
    int nx = 100;
    int ny = 100;
    pdd time = {0, 2.0};
    pdd domx = {-0.2, 0.2};
    pdd domy = {0, 0.8};

    //! 2D Riemann Problem (Configuration-3)
    // int nx = 100;
    // int ny = 100;
    // pdd time = {0, 1.0};
    // pdd domx = {0, 1.2};
    // pdd domy = {0, 1.2};

    //! Explosion Problem
    // int nx = 200;
    // int ny = 200;
    // pdd time = {0, 3.2};
    // pdd domx = {0, 1.5};
    // pdd domy = {0, 1.5};

    //! Implosion Problem - 1
    // int nx = 200;
    // int ny = 200;
    // pdd time = {0, 1.0};
    // pdd domx = {0, 1.2};
    // pdd domy = {0, 1.2};

    //! Implosion Problem - 2
    // int nx = 200;
    // int ny = 200;
    // pdd time = {0, 1.0};
    // pdd domx = {0, 1.2};
    // pdd domy = {0, 1.2};

    //! Kelvin Helmholtz Instability - 1
    // int nx = 200;
    // int ny = 200;
    // pdd time = {0, 1.0};
    // pdd domx = {0, 1.2};
    // pdd domy = {0, 1.2};

    //! Kelvin Helmholtz Instability - 2
    // int nx = 200;
    // int ny = 200;
    // pdd time = {0, 1.0};
    // pdd domx = {0, 1.2};
    // pdd domy = {0, 1.2};

    //! Kelvin Helmholtz Instability - 3
    // int nx = 200;
    // int ny = 200;
    // pdd time = {0, 1.0};
    // pdd domx = {0, 1.2};
    // pdd domy = {0, 1.2};

    //! Raleigh Taylor Instability
    // int nx = 200;
    // int ny = 200;
    // pdd time = {0, 1.0};
    // pdd domx = {0, 1.2};
    // pdd domy = {0, 1.2};


    double dx = (domx.second - domx.first) / nx;
    double dy = (domy.second - domy.first) / ny;

    //! Constants for all the other problems
    const double GAMMA = 1.4;
    const double THETA = 1.3;
    const double EPSILON = 1.0E-12;
    const double PI = 3.14159;
    const double CFL = 0.475;
    const double SMOOTHP = 0.00625;

    //! Constants for Raleigh Taylor Instability 
    // const double GAMMA = 5/3;
    // const double THETA = 1.3;
    // const double EPSILON = 1.0E-12;
    // const double CFL = 0.475;
    
    //! Switch for CURH Optimization block
    // const bool ShrinkSw = true;
    // const bool ScaleSw = true;

    const bool ShrinkSw = false;
    const bool ScaleSw = false;

    //! Information on the problem to be solved using the simulations
    const string PROBLEM = "MCW";
    // const string PROBLEM = "2DR";
    // const string PROBLEM = "EXP";

    const int MODE = 1; // CU
    // const int MODE = 2; // REF

    const string BC = "FREE";
    // const string BC = "REFLECTFREE";
    // const string BC = "SOLIDWALL";
    // const string BC = "PERIODIC";
}

