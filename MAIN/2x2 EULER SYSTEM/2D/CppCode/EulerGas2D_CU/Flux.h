/*
    * @author - lakshya Bamne (20MA20029)
    * @supervisor - Prof. Naveen Kumar Garg

    ! Implementation of CU scheme for modified Flux functions in Euler's Equations
*/

#pragma once

#include<iostream>
#include<vector>
#include<utility>
#include<unordered_map>
#include<string>
#include<cmath>

#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >

using namespace std;

/*
    Namespace containing functions to calculate flux vectors
*/
namespace Flux{
    vvvd fflux(vvd u, vvd rho);
    vvvd gflux(vvd u, vvd rho);
}

/*
    Function implementations
*/
vvvd Flux::fflux(vvd u, vvd rho){
    int row = rho.size();
    int col = rho[0].size();

    vvvd flux;

    vvd f1(row, vd(col));
    vvd f2(row, vd(col));

    for(int i=0 ; i<row ; i++){
        for(int j=0 ; j<col ; j++){
            f1[i][j] = pow(u[i][j],2);
            f2[i][j] = u[i][j]*rho[i][j];
        }
    }

    flux.push_back(f1);
    flux.push_back(f2);

    return flux;
}

vvvd Flux::gflux(vvd u, vvd rho){
    int row = rho.size();
    int col = rho[0].size();

    vvvd flux;

    vvd f1(row, vd(col));
    vvd f2(row, vd(col));

    for(int i=0 ; i<row ; i++){
        for(int j=0 ; j<col ; j++){
            f1[i][j] = pow(u[i][j],2);
            f2[i][j] = u[i][j]*rho[i][j];
        }
    }

    flux.push_back(f1);
    flux.push_back(f2);

    return flux;
}