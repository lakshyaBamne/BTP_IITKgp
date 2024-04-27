/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/
#pragma once

#include<iostream>
#include<vector>
#include<utility>
#include<unordered_map>
#include<string>

#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >

#include "PrimitiveVariables.h"

using namespace std;

namespace PRV = PrimitiveVariables;

/*
    Namespace containing functions to calculate flux vectors
*/
namespace Flux{
    vvvd fflux(string problem, vvd rho, vvd mx, vvd my, vvd E);
    vvvd gflux(string problem, vvd rho, vvd mx, vvd my, vvd E);
}

/*
    Function implementations
*/
vvvd Flux::fflux(string problem, vvd rho, vvd mx, vvd my, vvd E){
    int row = rho.size();
    int col = rho[0].size();

    vvvd flux;
    
    // get primitive variables for calculating flux easily
    vvd u = PRV::get_u(problem, rho, mx);
    vvd v = PRV::get_v(problem, rho, my);
    vvd p = PRV::get_p(problem, rho, u, v, E);

    vvd f1(row, vd(col));
    vvd f2(row, vd(col));
    vvd f3(row, vd(col));
    vvd f4(row, vd(col));

    for(int i=0 ; i<row ; i++){
        for(int j=0 ; j<col ; j++){
            f1[i][j] = mx[i][j];
            f2[i][j] = mx[i][j]*u[i][j] + p[i][j];
            f3[i][j] = mx[i][j]*v[i][j];
            f4[i][j] = u[i][j]*( E[i][j] + p[i][j] );
        }
    }

    flux.push_back(f1);
    flux.push_back(f2);
    flux.push_back(f3);
    flux.push_back(f4);

    return flux;
}

vvvd Flux::gflux(string problem, vvd rho, vvd mx, vvd my, vvd E){
    int row = rho.size();
    int col = rho[0].size();

    vvvd flux;
    
    // get primitive variables for calculating flux easily
    vvd u = PRV::get_u(problem, rho, mx);
    vvd v = PRV::get_v(problem, rho, my);
    vvd p = PRV::get_p(problem, rho, u, v, E);

    vvd f1(row, vd(col));
    vvd f2(row, vd(col));
    vvd f3(row, vd(col));
    vvd f4(row, vd(col));

    for(int i=0 ; i<row ; i++){
        for(int j=0 ; j<col ; j++){
            f1[i][j] = my[i][j];
            f2[i][j] = my[i][j]*u[i][j];
            f3[i][j] = my[i][j]*v[i][j] + p[i][j];
            f4[i][j] = v[i][j]*( E[i][j] + p[i][j] );
        }
    }

    flux.push_back(f1);
    flux.push_back(f2);
    flux.push_back(f3);
    flux.push_back(f4);

    return flux;
}
