#pragma once

#include<iostream>
#include<vector>
#include<utility>
#include<string>
#include<cmath>

#include "Constants.h"
#include "ExtendCells.h"
#include "OutputResult.h"

#define ll long long int

using namespace std;

namespace CTS = Constants;
namespace EXC = ExtendCells;
namespace OPF = OutputResult;

namespace InitRiemannProblem{
    //! 2D code
    // Function to initialize the computational grids in X and Y directions
    vector< vector<double> > make_grids(pair<double,double>& domainX, pair<double,double>& domainY, ll Nx, ll Ny, string result);
    // Function to initialize the conserved variables vector based on the Riemann Problem
    vector< vector< vector<double> > > get_conserved_variables(pair<string,string>& initial_conditions, pair<double,double>& domainX, pair<double,double>& domainY, ll Nx, ll Ny, string result);
}

vector< vector<double> > InitRiemannProblem::make_grids(pair<double,double>& domainX, pair<double,double>& domainY, ll Nx, ll Ny, string result){
    double dx = ( domainX.second - domainX.first ) / Nx;
    double dy = ( domainY.second - domainY.first ) / Ny;

    string file_grid = result + "/ComputationalDomain.txt";

    vector<double> gridX(Nx);
    vector<double> gridY(Ny);

    for(int i=0 ; i<Nx ; i++){
        gridX[i] = domainX.first + (i + 0.5)*dx;
    }

    for(int i=0 ; i<Ny ; i++){
        gridY[i] = domainY.first + (i + 0.5)*dy;
    }

    // write the grids into output file
    OPF::write_vector(gridX, file_grid);
    OPF::write_vector(gridY, file_grid);

    return {gridX, gridY};

}

vector< vector< vector<double> > > InitRiemannProblem::get_conserved_variables(pair<string,string>& initial_conditions, pair<double,double>& domainX, pair<double,double>& domainY, ll Nx, ll Ny, string result){
    vector< vector<double> > grid = InitRiemannProblem::make_grids(domainX, domainY, Nx, Ny, result);

    vector< vector<double> > rho(Nx, vector<double>(Ny, 0));
    vector< vector<double> > mx(Nx, vector<double>(Ny, 0));
    vector< vector<double> > my(Nx, vector<double>(Ny, 0));
    vector< vector<double> > E(Nx, vector<double>(Ny, 0));

    if( initial_conditions.first == "MCW" ){
        for(int i=0 ; i<Nx ; i++){
            for(int j=0 ; j<Ny ; j++){

                bool c1 = grid[0][i]<=0.1 && grid[0][i]>=-0.1 && grid[1][j]<=0.02 && grid[1][j]>=0;
                bool c2 = grid[0][i]<=0.02 && grid[0][i]>=-0.02 && grid[1][j]<=0.1 && grid[1][j]>=0.02;
                bool c3 = pow(grid[0][i]+0.02,2) + pow(grid[1][j]-0.02,2) <= pow(0.08,2);
                bool c4 = pow(grid[0][i]-0.02,2) + pow(grid[1][j]-0.02,2) <= pow(0.08,2);

                if( c1 || c2 || c3 || c4 ){
                    rho[i][j] = 1.4;
                    mx[i][j] = 0.0;
                    my[i][j] = 0.28;
                    E[i][j] = 1.0/(CTS::GAMMA-1) + 0.5*( pow(mx[i][j],2) + pow(my[i][j],2) )/rho[i][j];
                }
                else{
                    rho[i][j] = 1.0;
                    mx[i][j] = 0.0;
                    my[i][j] = 0.2;
                    E[i][j] = 1.0/(CTS::GAMMA-1) + 0.5*( pow(mx[i][j],2) + pow(my[i][j],2) )/rho[i][j];
                }

            }
        }
    }
    else if( initial_conditions.first == "2DR-CFG3" ){
        for(int i=0 ; i<Nx ; i++){
            for(int j=0 ; j<Ny ; j++){
                double u, v, p;

                if( grid[0][i]>=1.0 && grid[1][j]>=1.0 ){
                    rho[i][j] = 1.5; 
                    u = 0;
                    v = 0;
                    p = 1.5;
                }
                else if( grid[0][i]<1.0 && grid[1][j]>1.0 ){
                    rho[i][j] = 0.5323;
                    u = 1.206;
                    v = 0;
                    p = 0.3;
                }
                else if( grid[0][i]<1.0 && grid[1][j]<1.0 ){
                    rho[i][j] = 0.138;
                    u = 1.206;
                    v = 1.206;
                    p = 0.029;
                }
                else{
                    rho[i][j] = 0.5323;
                    u = 0;
                    v = 1.206;
                    p = 0.3;
                }

                mx[i][j] = rho[i][j] * u;
                my[i][j] = rho[i][j] * v;
                E[i][j] = p/(CTS::GAMMA-1) + 0.5*( pow(mx[i][j],2) + pow(my[i][j],2) )/rho[i][j]; 
            }
        } 
    }
    else if( initial_conditions.first == "EXP" ){
        for(int i=0 ; i<Nx ; i++){
            for(int j=0 ; j<Ny ; j++){
                
                if( pow(grid[0][i],2) + pow(grid[1][j],2) < 0.16 ){
                    rho[i][j] = 1.0;
                    mx[i][j] = 0;
                    my[i][j] = 0;
                    E[i][j] = 1.0/(CTS::GAMMA-1);
                }
                else{
                    rho[i][j] = 0.125;
                    mx[i][j] = 0;
                    my[i][j] = 0;
                    E[i][j] = 0.1/(CTS::GAMMA-1);
                }

            }
        }
    }
    else{
        cout << "---ERROR--- Please select correct problem ---" << endl;
    }

    return {rho, mx, my, E};


}

