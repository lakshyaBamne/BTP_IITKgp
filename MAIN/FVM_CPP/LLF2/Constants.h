#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<cmath>
#include<utility>
#include<fstream>
#include<sstream>

using namespace std;

#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define pdd pair<double,double>
#define pss pair<string,string>

namespace Constants{
    //! Global constants
    double GAMMA = 1.4;
    double THETA = 1.3;
    double EPS = 1.0E-12;

    //! CFL number
    // double CFL = 0.9; // 1st order schemes
    double CFL = 0.45; // 2nd order schemes

    //! domain endpoints
    pdd domainX = {0, 1};
    // pdd domainX = {-5, 5}; // LAX
    // pdd domainX = {-5, 15}; // SDW
    // pdd domainX = {-5,5}; // SEW

    //! time endpoints
    pdd time = {0, 2}; // MCW
    // pdd time = {0, 0.012}; // SCW
    // pdd time = {0, 0.038}; // BLW
    // pdd time = {0, 1.3}; // LAX
    // pdd time = {0, 0.2}; // TORO-1
    // pdd time = {0, 0.15}; // TORO-2
    // pdd time = {0, 0.012}; // TORO-3
    // pdd time = {0, 0.035}; // TORO-4
    // pdd time = {0, 5}; // SDW
    // pdd time = {0, 5}; // SEW

    //! problem 
    string PROBLEM = "MCW";
    // string PROBLEM = "SCW";
    // string PROBLEM = "BLW";
    // string PROBLEM = "LAX";
    // string PROBLEM = "TORO-1";
    // string PROBLEM = "TORO-2";
    // string PROBLEM = "TORO-3";
    // string PROBLEM = "TORO-4";
    // string PROBLEM = "SDW";
    // string PROBLEM = "SEW";

    string BC = "FREE";
    // string BC = "REFLECTIVE";

    //! number of finite volume cells
    int N = 1000;
    int Nref = 4000;

    // length of an individual finite volume cell
    double dx = (domainX.second-domainX.first)/N;
}
