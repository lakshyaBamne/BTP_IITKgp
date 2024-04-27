#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>
#include<cmath>

using namespace std;

#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define umap_ss unordered_map<string,string>
#define pss pair<string,string>
#define pdd pair<double,double>
#define pii pair<int,int>

namespace Constants{
    //! Parameters for the problem
    const double CFL = 0.9; // for 1st order approximations
    // const double CFL = 0.45; // for 2nd order approximations
    const double GAMMA = 1.4;
    const double THETA = 1.3;

    //! problem specific variables
    const int n = 100;
    const pdd domx = {0.0, 1.0};
    const pdd time = {0.0, 0.012};
    const string problem = "SCW";
    const string bc = "FREE";
    
    const double dx = (domx.second - domx.first) / n;
}

