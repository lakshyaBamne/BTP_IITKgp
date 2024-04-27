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
    // Function to initialize the grid based on initial conditions entered by the user
    vector<double> make_grid(pair<double,double>& domain, ll Nx, string result);
    // Function to Initialize the Conserved variables vectors based on the Riemann Problem
    vector<vector<double>> get_conserved_variables(pair<string,string>& initial_conditions, pair<double,double>& domain, ll Nx, string result);
}

// Function to initialize the computational grid
vector<double> InitRiemannProblem::make_grid(pair<double,double>& domain, ll Nx, string result){
    double xlen = domain.second - domain.first;
    double dx = xlen/Nx;

    string file_grid = result + "/ComputationalDomain.txt";

    vector<double> grid(Nx);

    for(int i=0 ; i<Nx ; i++){
        grid[i] = domain.first + (i+0.5)*dx;
    }

    // write the grid in a file
    OPF::write_vector(grid, file_grid);

    return grid;
}

// Function to get the Conserved variable vectors according to the Given Riemann problem
vector<vector<double>> InitRiemannProblem::get_conserved_variables(pair<string,string>& initial_conditions, pair<double,double>& domain, ll Nx, string result){
    vector<double> grid = InitRiemannProblem::make_grid(domain, Nx, result); // initialize the grid

    vector<vector<double>> cons_vars;

    vector<double> rho(grid.size()); // density
    vector<double> u(grid.size()); // velocity
    vector<double> p(grid.size()); // pressure
    vector<double> m(grid.size()); // momentum
    vector<double> E(grid.size()); // Energy

    if( initial_conditions.first == "MCW" ){
        // Moving Contact Wave (MCW)
        // x[0,1] -> 0.3
        // t[0,2]
        for(int i=0 ; i<grid.size() ; i++){
            if( grid[i] < 0.3 ){
                rho[i] = 1.4;
            }
            else{
                rho[i] = 1.0;
            }

            u[i] = 0.1;
            p[i] = 1.0;

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }
    }
    else if( initial_conditions.first == "SCW" ){
        // Stationary contact wave and travelling shock and rarefaction (SCW)
        // x[0,1] -> 0.8
        // t[0,0.012]
        for(int i=0 ; i<grid.size() ; i++){
            rho[i] = 1.0;
            u[i] = -19.59745;
            
            if( grid[i] < 0.8 ){
                p[i] = 1000.0;
            }
            else{
                p[i] = 0.01;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;    
        }

    }
    else if( initial_conditions.first == "BLW" ){
        // Blast wave problem
        // x[0,1] -> 0.1 and 0.9
        // t[0, 0.038]
        // Reflective Boundary Conditions
        for(int i=0 ; i<grid.size() ; i++){
            rho[i] = 1.0;
            u[i] = 0.0;
            
            if( grid[i] < 0.1 ){
                p[i] = 1000.0;
            }
            else if( grid[i]>=0.1 && grid[i]<=0.9 ){
                p[i] = 0.01;
            }
            else{
                p[i] = 100.0;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }
    }
    else if( initial_conditions.first == "LAX" ){
        // Lax problem
        // x[-5,5] -> 0.0
        // t[0,1.3]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < 0.0 ){
                rho[i] = 0.445;
                u[i] = 0.698;
                p[i] = 3.528;
            }
            else{
                rho[i] = 0.5;
                u[i] = 0.0;
                p[i] = 0.571;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }
    }
    else if( initial_conditions.first == "SOD" ){
        // Sod's shock tube problem
        // x[-10,10] -> 0
        // t[0,0.01]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < 0.0 ){
                rho[i] = 1.0;
                u[i] = 0.0;
                p[i] = 100000.0;
            }
            else{
                rho[i] = 0.125;
                u[i] = 0.0;
                p[i] = 10000.0;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }

    }
    else if( initial_conditions.first == "SPP" ){
        // Sonic Point Problem
        // x[0,1] -> 0.3
        // t[0,0.2]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < 0.3 ){
                rho[i] = 1.0;
                u[i] = 0.75;
                p[i] = 1.0;
            }
            else{
                rho[i] = 0.125;
                u[i] = 0.0;
                p[i] = 0.1;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }

    }
    else if( initial_conditions.first == "SSW" ){
        // Strong Shock wave
        // x[0,1] -> 0.5
        // t[0,0.012]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < 0.5 ){
                rho[i] = 1.0;
                u[i] = 0.0;
                p[i] = 1000.0;
            }
            else{
                rho[i] = 1.0;
                u[i] = 0.0;
                p[i] = 0.01;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }

    }
    else if( initial_conditions.first == "MA3" ){
        // Mach-3 Problem
        // x[0,1] -> 0.4
        // t[0,0.1]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < 0.4 ){
                rho[i] = 3.857;
                u[i] = 0.92;
                p[i] = 10.333;
            }
            else{
                rho[i] = 1.0;
                u[i] = 3.55;
                p[i] = 1.0;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }

    }
    else if( initial_conditions.first == "SDW" ){
        // Shock density wave problem
        // x[-5,15] -> -4
        // t[0,5]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < -4 ){
                rho[i] = 27.0/7.0;
                u[i] = ( 4*sqrtf(35) )/9.0;
                p[i] = 31.0/3.0;
            }
            else{
                rho[i] = 1.0 + 0.2*( sin( 5*grid[i] ) );
                u[i] = 0.0;
                p[i] = 1.0;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }

    }
    else if( initial_conditions.first == "SEW" ){
        // Shock entropy wave problem
        // x[-5,5] -> -4.5
        // t[0,5]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < -4.5 ){
                rho[i] = 1.51695;
                u[i] = 0.523346;
                p[i] = 1.805;
            }
            else{
                rho[i] = 1 + 0.1*( sin( 20*grid[i] ) );
                u[i] = 0.0;
                p[i] = 1.0;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }

    }
    else if( initial_conditions.first == "TORO-1" ){
        // Toro test-1
        // x[0,1] -> 0.3
        // t[0,0.2]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < 0.3 ){
                rho[i] = 1.0;
                u[i] = 0.75;
                p[i] = 1.0;
            }
            else{
                rho[i] = 0.125;
                u[i] = 0.0;
                p[i] = 0.1;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }
    }
    else if( initial_conditions.first == "TORO-2" ){
        // Toro test-2
        // x[0,1] -> 0.5
        // t[0,0.15]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < 0.5 ){
                rho[i] = 1.0;
                u[i] = -2.0;
                p[i] = 0.4;
            }
            else{
                rho[i] = 1.0;
                u[i] = 2.0;
                p[i] = 0.4;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }
    }
    else if( initial_conditions.first == "TORO-3" ){
        // Toro test-3
        // x[0,1] -> 0.5
        // t[0,0.012]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < 0.5 ){
                rho[i] = 1.0;
                u[i] = 0.0;
                p[i] = 1000.0;
            }
            else{
                rho[i] = 1.0;
                u[i] = 0.0;
                p[i] = 0.01;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }
    }
    else if( initial_conditions.first == "TORO-4" ){
        // Toro test-4
        // x[0,1] -> 0.4
        // t[0,0.035]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < 0.4 ){
                rho[i] = 5.99924;
                u[i] = 19.5975;
                p[i] = 460.894;
            }
            else{
                rho[i] = 5.99242;
                u[i] = -6.19633;
                p[i] = 46.0950;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }
    }
    else if( initial_conditions.first == "TORO-5" ){
        // Toro test-5
        // x[0,1] -> 0.4
        // t[0,0.035]
        for(int i=0 ; i<grid.size() ; i++){
            
            if( grid[i] < 0.4 ){
                rho[i] = 5.99924;
                u[i] = 19.5975;
                p[i] = 460.894;
            }
            else{
                rho[i] = 5.99242;
                u[i] = -6.19633;
                p[i] = 46.0950;
            }

            m[i] = rho[i]*u[i];
            E[i] = p[i]/(CTS::GAMMA-1) + m[i]*u[i]*0.5;
        }
    }
    else{
        cout << "---ERROR--- Please enter correct Riemann Problem ---" << endl;
    }

    cons_vars.push_back(rho);
    cons_vars.push_back(m);
    cons_vars.push_back(E);

    // We need to extend the cells for the Ghost values as well
    EXC::extend_cells(initial_conditions.second, cons_vars);

    return cons_vars;
}

