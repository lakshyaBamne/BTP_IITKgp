/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>
#include<cmath>
#include<sstream>
#include<fstream>

#include "ExtendCells.h"
#include "OutputResult.h"
#include "Constants.h"
#include "PrimitiveVariables.h"
#include "Utility.h"
#include "Flux.h"

#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define vvvvd vector< vector< vector< vector<double> > > >
#define umap_ss unordered_map<string,string>

using namespace std;

namespace EXC = ExtendCells;
namespace OPR = OutputResult;
namespace CTS = Constants;
namespace PRV = PrimitiveVariables;
namespace UTL = Utility;
namespace FLX = Flux;

// class RiemannProblem(RP) to initialize the 2D riemann problem for a given problem
class RP{
public: // everything is public

    // required later for various purposes
    string mode;
    string problem;

    pair<double,double> domX;
    pair<double,double> domY;

    pair<double,double> time;

    int Nx;
    int Ny;

    double dx;
    double dy;

    double t;
    double dt; // time discretization step updated using CFL conditions

    // grid vectors are initialized after initializing the variables
    vector< double > gridX;
    vector< double > gridY;

    // Variable matrices to store information
    vector< vector<double> > rho;   
    vector< vector<double> > mx;
    vector< vector<double> > my;
    vector< vector<double> > E;

    unordered_map<string,string> bc; // boundary conditions for N,S,E,W in order

    /*
        Constructor takes the problem string as input to initialize the variables
    
        -> [MCW] Moving Contact Wave
        -> [2DR] 2-Dimensional Riemann Problem
        -> [EXP] Explosion Problem
        -> [IMP] Implosion Problem
        -> [KHI] Kelvin Helmholtz Instability
        -> [RTI] Raleigh Taylor Instability
    */
    RP(string mode, string problem) : mode(mode), problem(problem){

        cout << "---LOG--- Constructor running ---" << endl;

        // initialize the environment and set up initial variables and matrices
        initialize_environment();        

        cout << "---LOG--- environment write successful ---" << endl;
    }

    // Function to initialize environment variables and set up initial matrices
    void initialize_environment(){
        if( mode == "CU" ){
            // CU plot is the plot for the Numerical Scheme itself
            if( problem == "MCW" ){
                domX.first = -0.2;
                domX.second = 0.2;

                domY.first = 0;
                domY.second = 0.8;

                dx = (double)1/250;
                dy = (double)1/250;

                time.first = 0;
                time.second = 2;

                t = time.first;

                bc["N"] = "FREE";
                bc["S"] = "FREE";
                bc["E"] = "FREE";
                bc["W"] = "FREE";
            }
            else if( problem == "2DR" ){
                domX.first = 0;
                domX.second = 1.2;

                domY.first = 0;
                domY.second = 1.2;

                dx = (double)3/2500;
                dy = (double)3/2500;

                time.first = 0;
                time.second = 1;

                t = time.first;

                bc["N"] = "FREE";
                bc["S"] = "FREE";
                bc["E"] = "FREE";
                bc["W"] = "FREE";
            }
            else if( problem == "EXP" ){
                domX.first = 0;
                domX.second = 1.5;

                domY.first = 0;
                domY.second = 1.5;

                dx = (double)3/800;
                dy = (double)3/800;

                time.first = 0;
                time.second = 3.2;

                t = time.first;

                bc["N"] = "FREE";
                bc["S"] = "REF";
                bc["E"] = "FREE";
                bc["W"] = "REF";
            }
            else if( problem == "IMP" ){
                domX.first = 0;
                domX.second = 0.3;

                domY.first = 0;
                domY.second = 0.3;

                dx = (double)3/6000;
                dy = (double)3/6000;

                time.first = 0;
                time.second = 2.5;

                t = time.first;

                bc["N"] = "REF";
                bc["S"] = "REF";
                bc["E"] = "REF";
                bc["W"] = "REF";
            }
            else if( problem == "KHI" ){
                domX.first = -0.5;
                domX.second = 0.5;

                domY.first = -0.5;
                domY.second = 0.5;

                dx = (double)3/1024;
                dy = (double)3/1024;

                time.first = 0;
                time.second = 2.5; // evaluate results at t=1, t=2.5 and t=4

                t = time.first;

                bc["N"] = "PER";
                bc["S"] = "PER";
                bc["E"] = "PER";
                bc["W"] = "PER";
            }
            else if( problem == "RTI" ){
                domX.first = 0;
                domX.second = 0.25;

                domY.first = 0;
                domY.second = 1.0;

                dx = (double)1/1024;
                dy = (double)1/1024;

                time.first = 0;
                time.second = 2.95; // also evaluate results at t=2.5 and t=4

                t = time.first;

                bc["N"] = "DIR";
                bc["S"] = "DIR";
                bc["E"] = "REF";
                bc["W"] = "REF";
            }
            else{
                cout << "---ERROR--- Please enter correct problem ---" << endl;
            }
        }
        else if( mode == "REF" ){
            // reference is taken as 4 times finer grid
            if( problem == "MCW" ){
                domX.first = -0.2;
                domX.second = 0.2;

                domY.first = 0;
                domY.second = 0.8;

                dx = (double)1/1000;
                dy = (double)1/1000;

                time.first = 0;
                time.second = 2;

                t = time.first;

                bc["N"] = "FREE";
                bc["S"] = "FREE";
                bc["E"] = "FREE";
                bc["W"] = "FREE";
            }
            else if( problem == "2DR" ){
                domX.first = 0;
                domX.second = 1.2;

                domY.first = 0;
                domY.second = 1.2;

                dx = (double)3/10000;
                dy = (double)3/10000;

                time.first = 0;
                time.second = 1;

                t = time.first;

                bc["N"] = "FREE";
                bc["S"] = "FREE";
                bc["E"] = "FREE";
                bc["W"] = "FREE";
            }
            else if( problem == "EXP" ){
                domX.first = 0;
                domX.second = 1.5;

                domY.first = 0;
                domY.second = 1.5;

                dx = (double)3/3200;
                dy = (double)3/3200;

                time.first = 0;
                time.second = 3.2;

                t = time.first;

                bc["N"] = "FREE";
                bc["S"] = "REF";
                bc["E"] = "FREE";
                bc["W"] = "REF";
            }
            else if( problem == "IMP" ){
                domX.first = 0;
                domX.second = 0.3;

                domY.first = 0;
                domY.second = 0.3;

                dx = (double)3/24000;
                dy = (double)3/24000;

                time.first = 0;
                time.second = 2.5;

                t = time.first;

                bc["N"] = "REF";
                bc["S"] = "REF";
                bc["E"] = "REF";
                bc["W"] = "REF";
            }
            else if( problem == "KHI" ){
                domX.first = -0.5;
                domX.second = 0.5;

                domY.first = -0.5;
                domY.second = 0.5;

                dx = (double)3/4096;
                dy = (double)3/4096;

                time.first = 0;
                time.second = 1; // also evaluate results at t=2.5 and t=4

                t = time.first;

                bc["N"] = "PER";
                bc["S"] = "PER";
                bc["E"] = "PER";
                bc["W"] = "PER";
            }
            else if( problem == "RTI" ){
                domX.first = 0;
                domX.second = 0.25;

                domY.first = 0;
                domY.second = 1.0;

                dx = (double)1/4096;
                dy = (double)1/4096;

                time.first = 0;
                time.second = 2.95; // also evaluate results at t=2.5 and t=4

                t = time.first;

                bc["N"] = "DIR";
                bc["S"] = "DIR";
                bc["E"] = "REF";
                bc["W"] = "REF";
            }
            else{
                cout << "---ERROR--- Please enter correct problem ---" << endl;
            }
        }
        else{
            cout << "---ERROR--- Please enter correct mode for initialization ---" << endl;
        }

        // calculating the number of points in the grid using the cell volumes
        Nx = (int)( (domX.second-domX.first) / dx );
        Ny = (int)( (domY.second-domY.first) / dy );

        // LOG
        show_var();

        cout << "---LOG--- Initialized Riemann problem variables ---" << endl;

        // initialize grid vectors
        make_grids();

        cout << "---LOG--- Initialized computational grid ---" << endl;

        OPR::write_grids(mode, gridX, gridY);

        // initialize variable matrices
        init_vars();

        cout << "---LOG--- Initialized Variable matrices ---" << endl;

        // export environment
        export_env(mode);
    }

    // function to initialize the computational grid
    void make_grids(){
        // gridX
        for(int i=0 ; i<Nx ; i++){
            gridX.push_back( domX.first + dx*( i + 0.5 ) );
        }

        // gridY
        for(int i=0 ; i<Ny ; i++){
            gridY.push_back( domY.first + dy*( i + 0.5 ) );
        }
    }

    // function to initialize the variable grids (rho, mx, my, E)
    void init_vars(){

        if( problem == "MCW" ){
            
            for(auto j : gridY){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : gridX){
                    bool c1 = i>-0.1 && i<0.1 && j>0 && j<0.02;
                    bool c2 = i>-0.02 && i<0.02 && j>0.02 && j<0.1;
                    bool c3 = pow(i+0.02,2) + pow(j-0.02,2) < pow(0.08,2);
                    bool c4 = pow(i-0.02,2) + pow(j-0.02,2) < pow(0.08,2);

                    if( c1 || c2 || c3 || c4 ){
                        rho1.push_back( 1.4 );
                    }
                    else{
                        rho1.push_back( 1.0 );
                    }

                    u = 0.0;
                    v = 0.2;
                    p = 1;

                    mx1.push_back( rho1.back() * u );
                    my1.push_back( rho1.back() * v );

                    E1.push_back( (p/CTS::GAMMA) + 0.5*(mx1.back()*u + my1.back()*v) );

                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }
        }
        else if( problem == "2DR" ){

            for(auto j : gridY){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : gridX){
                    
                    if( i>1 && j>1 ){
                        rho1.push_back( 1.5 );
                        u = 0.0;
                        v = 0.0;
                        p = 1.5;
                    }
                    else if( i<1 && j>1 ){
                        rho1.push_back( 0.5323 );
                        u = 1.206;
                        v = 0.0;
                        p = 0.3;
                    }
                    else if( i<1 && j<1 ){
                        rho1.push_back( 0.138 );
                        u = 1.206;
                        v = 1.206;
                        p = 0.029;
                    }
                    else if( i>1 && j<1 ){
                        rho1.push_back( 0.5323 );
                        u = 0.0;
                        v = 1.206;
                        p = 0.3;
                    }

                    mx1.push_back( rho1.back() * u );
                    my1.push_back( rho1.back() * v );

                    E1.push_back( (p/CTS::GAMMA) + 0.5*(mx1.back()*u + my1.back()*v) );

                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }

        }
        else if( problem == "EXP" ){

            for(auto j : gridY){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : gridX){
                    
                    if( pow(i,2)+pow(j,2) < 0.16 ){
                        rho1.push_back(1);
                        u = 0;
                        v = 0;
                        p = 1;
                    }
                    else{
                        rho1.push_back(0.125);
                        u = 0;
                        v = 0;
                        p = 0.1;
                    }

                    mx1.push_back( rho1.back() * u );
                    my1.push_back( rho1.back() * v );

                    E1.push_back( (p/CTS::GAMMA) + 0.5*(mx1.back()*u + my1.back()*v) );

                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }

        }
        else if( problem == "IMP" ){

            for(auto j : gridY){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : gridX){
                    
                    if( abs(i) + abs(j) < 0.15 ){
                        rho1.push_back(0.125);
                        u = 0;
                        v = 0;
                        p = 0.14;
                    }
                    else{
                        rho1.push_back(1);
                        u = 0;
                        v = 0;
                        p = 1;
                    }

                    mx1.push_back( rho1.back() * u );
                    my1.push_back( rho1.back() * v );

                    E1.push_back( (p/CTS::GAMMA) + 0.5*(mx1.back()*u + my1.back()*v) );

                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }

        }
        else if( problem == "KHI" ){

            for(auto j : gridY){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : gridX){
                    
                    if( j>=-0.5 && j<-0.25 ){
                        rho1.push_back(1);
                        u = -0.5 + 0.5*expf( (j+0.25)/CTS::SMOOTHP );
                    }
                    else if( j>=-0.25 && j<0 ){
                        rho1.push_back(2);
                        u = 0.5 - 0.5*expf( (-j-0.25)/CTS::SMOOTHP );
                    }
                    else if( j>=0 && j<0.25 ){
                        rho1.push_back(2);
                        u = 0.5 - 0.5*expf( (j-0.25)/CTS::SMOOTHP );
                    }
                    else if( j>=0.25 && j<=0.5 ){
                        rho1.push_back(1);
                        u = -0.5 + 0.5*expf( (0.25-j)/CTS::SMOOTHP );
                    }
                    
                    v = 0.01*sin(4*CTS::PI*i);
                    p = 1.5;

                    mx1.push_back( rho1.back() * u );
                    my1.push_back( rho1.back() * v );

                    E1.push_back( (p/CTS::GAMMA) + 0.5*(mx1.back()*u + my1.back()*v) );

                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }

        }
        else if( problem == "RTI" ){

            for(auto j : gridY){
                vector< double > rho1;
                vector< double > mx1;
                vector< double > my1;
                vector< double > E1;

                double u,v,p;

                for(auto i : gridX){
                    
                    if( j<0.5 ){
                        rho1.push_back( 2 );
                        u = 0;
                        p = 2*j + 1;
                    }
                    else{
                        rho1.push_back( 1 );
                        u = 0;
                        p = j + 1.5;
                    }

                    double c = sqrtf( (CTS::GAMMA_RT*p)/rho1.back() );

                    v = -0.025*cos(8*CTS::PI*i)*c;

                    mx1.push_back( rho1.back() * u );
                    my1.push_back( rho1.back() * v );

                    E1.push_back( (p/CTS::GAMMA) + 0.5*(mx1.back()*u + my1.back()*v) );

                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }

        }   
        else{
            cout << "---ERROR--- Please enter correct problem ---" << endl;
        }

        cout << "---LOG--- Initialized the variables ---" << endl;

        // extend cells
        EXC::extend_matrix(rho, bc);
        EXC::extend_matrix(mx, bc);
        EXC::extend_matrix(my, bc);
        EXC::extend_matrix(E, bc);
        
        cout << "---LOG--- Extended cells successfully ---" << endl;

    }

    // function to export the environment variables in a file -> Environment.txt
    void export_env(string mode){
        string file_name;

        if( mode == "CU" ){
            file_name = "env1/Environment.txt";
        }
        else{
            file_name = "env2/Environment.txt";
        }
        
        stringstream ss;
        
        ofstream fout;
        fout.open(file_name , ios::app);

        ss << mode << endl;
        ss << problem << endl;
        ss << domX.first << " " << domX.second << endl;
        ss << domY.first << " " << domY.second << endl;
        ss << time.first << " " << time.second << endl;
        ss << Nx << endl;
        ss << Ny << endl;

        fout << ss.str() << endl;

        fout.close();
    }

    // function to print the environment for the class
    void show_var(){
        cout << "---------------------------------Class Variables--------------------------------------------" << endl;

        cout << "MODE : " << mode << endl;
        cout << "PROBLEM : " << problem << endl;
        cout << "X -> [ " << domX.first << " , " << domX.second << " ]" << endl;
        cout << "Y -> [ " << domY.first << " , " << domY.second << " ]" << endl;
        cout << "T -> [ " << time.first << " , " << time.second << " ]" << endl;
        cout << "Nx : " << Nx << endl;
        cout << "Ny : " << Ny << endl;
        cout << "dx : " << dx << endl;
        cout << "dy : " << dy << endl;
        cout << "Boundary Conditions :-" << endl;
        cout << "N : " << bc["N"] << endl;
        cout << "S : " << bc["S"] << endl;
        cout << "E : " << bc["E"] << endl;
        cout << "W : " << bc["W"] << endl;

        cout << "--------------------------------------------------------------------------------------------" << endl;
    }

    // function to run the CU Scheme
    // and store data for every iteration
    void RunCU_complete(){

        vector< vector<double> > u;
        vector< vector<double> > v;
        vector< vector<double> > p;

        while( t < time.second ){
            cout << "t = " << t << " | dt = " << dt << endl;
            
            //! write variable values at each step for making animations
            OPR::write_matrix(mode, "Density", rho);
            OPR::write_matrix(mode, "MomentumX", mx);
            OPR::write_matrix(mode, "MomentumY", my);
            OPR::write_matrix(mode, "Energy", E);

            // get the primitive variables
            u = PRV::get_u(problem, rho, mx);
            v = PRV::get_v(problem, rho, my);
            p = PRV::get_p(problem, rho, u, v, E);

            OPR::write_matrix(mode, "VelocityX", u);
            OPR::write_matrix(mode, "VelocityY", v);
            OPR::write_matrix(mode, "Pressure", p);

            // function to update the conservative variables internally
            // also update dt for next iteration
            get_next_vars();

            // update new time step calculated using CFL conditions
            t += dt;

        }
    }

    // function to run the CU Scheme
    // and store data only for initial and final iterations
    void RunCU_partial(){

        //! write initial values for plotting
        OPR::write_matrix(mode, "Density", rho);
        OPR::write_matrix(mode, "MomentumX", mx);
        OPR::write_matrix(mode, "MomentumY", my);
        OPR::write_matrix(mode, "Energy", E);

        // get the primitive variables
        vector< vector<double> > u = PRV::get_u(problem, rho, mx);
        vector< vector<double> > v = PRV::get_v(problem, rho, my);
        vector< vector<double> > p = PRV::get_p(problem, rho, u, v, E);

        OPR::write_matrix(mode, "VelocityX", u);
        OPR::write_matrix(mode, "VelocityY", v);
        OPR::write_matrix(mode, "Pressure", p);        

        while( t < time.second ){
            cout << "t = " << t << " | dt = " << dt << endl;
            
            // function to update the conservative variables internally
            // also update dt for next iteration
            get_next_vars();

            // update new time step calculated using CFL conditions
            t += dt;

        }

        //! write final values for plotting
        OPR::write_matrix(mode, "Density", rho);
        OPR::write_matrix(mode, "MomentumX", mx);
        OPR::write_matrix(mode, "MomentumY", my);
        OPR::write_matrix(mode, "Energy", E);

        // get the primitive variables
        u = PRV::get_u(problem, rho, mx);
        v = PRV::get_v(problem, rho, my);
        p = PRV::get_p(problem, rho, u, v, E);

        OPR::write_matrix(mode, "VelocityX", u);
        OPR::write_matrix(mode, "VelocityY", v);
        OPR::write_matrix(mode, "Pressure", p);       
    }

    /*
        Functions to Run the CU Scheme
    */

    // Function to run the three stage SSPRK scheme to get the next iteration values for the conserved variables
    // vars -> already present values for the Conservative variables
    // next_vars -> new values for the Conservative variables found using SSPRK and CU Flux
    void get_next_vars(){
        // using class variables as input
        vvvvd cu_flux = get_cu_flux(rho, mx, my, E);
        
        vvvd F = cu_flux[0];
        vvvd G = cu_flux[1];

        const double old_dt = dt; // used because dt is changed later which is undesirable

        double LAMBDA = dt/dx;
        double MU = dt/dy;

        int row = rho.size();
        int col = rho[0].size();

        // STAGE-1
        vvd rho1(row-2, vd(col-2));
        vvd mx1(row-2, vd(col-2));
        vvd my1(row-2, vd(col-2));
        vvd E1(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1; j<col-1 ; j++){
                rho1[i-1][j-1] = rho[i][j] - ( LAMBDA*( F[0][i][j] - F[0][i][j-1] ) + MU*( G[0][i][j] - G[0][i-1][j] ) );
                mx1[i-1][j-1] = mx[i][j] - ( LAMBDA*( F[1][i][j] - F[1][i][j-1] ) + MU*( G[1][i][j] - G[1][i-1][j] ) );
                my1[i-1][j-1] = my[i][j] - ( LAMBDA*( F[2][i][j] - F[2][i][j-1] ) + MU*( G[2][i][j] - G[2][i-1][j] ) );
                E1[i-1][j-1] = E[i][j] - ( LAMBDA*( F[3][i][j] - F[3][i][j-1] ) + MU*( G[3][i][j] - G[3][i-1][j] ) );
            }
        }

        EXC::extend_matrix(rho1, bc);
        EXC::extend_matrix(mx1, bc);
        EXC::extend_matrix(my1, bc);
        EXC::extend_matrix(E1, bc);

        // STAGE-2
        vvvvd cu_flux_1 = get_cu_flux(rho1, mx1, my1, E1);

        dt = old_dt; // restore value of dt to old value

        vvvd F1 = cu_flux_1[0];
        vvvd G1 = cu_flux_1[1];

        vvd rho2(row-2, vd(col-2));
        vvd mx2(row-2, vd(col-2));
        vvd my2(row-2, vd(col-2));
        vvd E2(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1; j<col-1 ; j++){
                rho2[i-1][j-1] = ( 3*rho[i][j] + rho1[i][j] - ( LAMBDA*( F1[0][i][j] - F1[0][i][j-1] ) + MU*( G1[0][i][j] - G1[0][i-1][j] ) ) ) / 4;
                mx2[i-1][j-1] = ( 3*mx[i][j] + mx1[i][j] - ( LAMBDA*( F1[1][i][j] - F1[1][i][j-1] ) + MU*( G1[1][i][j] - G1[1][i-1][j] ) ) ) / 4;
                my2[i-1][j-1] = ( 3*my[i][j] + my1[i][j] - ( LAMBDA*( F1[2][i][j] - F1[2][i][j-1] ) + MU*( G1[2][i][j] - G1[2][i-1][j] ) ) ) / 4;
                E2[i-1][j-1] = ( 3*E[i][j] + E1[i][j] - ( LAMBDA*( F1[3][i][j] - F1[3][i][j-1] ) + MU*( G1[3][i][j] - G1[3][i-1][j] ) ) ) / 4;
            }
        }

        EXC::extend_matrix(rho2, bc);
        EXC::extend_matrix(mx2, bc);
        EXC::extend_matrix(my2, bc);
        EXC::extend_matrix(E2, bc);

        // STAGE-3
        vvvvd cu_flux_new = get_cu_flux(rho2, mx2, my2, E2);

        dt = old_dt; // restore value of dt to old value

        vvvd FN = cu_flux_new[0];
        vvvd GN = cu_flux_new[1];

        vvd rhon(row-2, vd(col-2));
        vvd mxn(row-2, vd(col-2));
        vvd myn(row-2, vd(col-2));
        vvd En(row-2, vd(col-2));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1; j<col-1 ; j++){
                rhon[i-1][j-1] = ( rho[i][j] + 2*rho2[i][j] - 2*( LAMBDA*( FN[0][i][j] - FN[0][i][j-1] ) + MU*( GN[0][i][j] - GN[0][i-1][j] ) ) ) / 3;
                mxn[i-1][j-1] = ( mx[i][j] + 2*mx2[i][j] - 2*( LAMBDA*( FN[1][i][j] - FN[1][i][j-1] ) + MU*( GN[1][i][j] - GN[1][i-1][j] ) ) ) / 3;
                myn[i-1][j-1] = ( my[i][j] + 2*my2[i][j] - 2*( LAMBDA*( FN[2][i][j] - FN[2][i][j-1] ) + MU*( GN[2][i][j] - GN[2][i-1][j] ) ) ) / 3;
                En[i-1][j-1] = ( E[i][j] + 2*E2[i][j] - 2*( LAMBDA*( FN[3][i][j] - FN[3][i][j-1] ) + MU*( GN[3][i][j] - GN[3][i-1][j] ) ) ) / 3;
            }
        }

        EXC::extend_matrix(rhon, bc);
        EXC::extend_matrix(mxn, bc);
        EXC::extend_matrix(myn, bc);
        EXC::extend_matrix(En, bc);

        // Updation of conserved variables to new values
        rho = rhon;
        mx = mxn;
        my = myn;
        E = En;


    }

    // function to get the CU flux final for any given 2D grid of Conserved variables
    // no arguments required because the grids are already class variables accessible to all functions
    //! returns the flux matrices F and G which are then used in SSPRK scheme to construct next iteration
    vvvvd get_cu_flux(vvd& rho, vvd& mx, vvd& my, vvd& E){
        /*
            STEP-1 : Piecewise Linear Reconstruction [N, S, E, W, NE, NW, SE, SW]
        */
        vvvd rho_plr = get_plr(rho);
        vvvd mx_plr = get_plr(mx);
        vvvd my_plr = get_plr(my);
        vvvd E_plr = get_plr(E);

        /*
            STEP-2 : Piecewise Linear Reconstruction [a+, a-, b+, b-]
        */
        vvvd lsp_x = get_lsp_a(rho_plr, mx_plr, my_plr, E_plr);
        vvvd lsp_y = get_lsp_b(rho_plr, mx_plr, my_plr, E_plr);

        /*
            STEP-3 : Update time step using CFL Conditions
            -> since dt is a class variable it is updated internally
        */
        new_dt(lsp_x, lsp_y);

        /*
            STEP-4 : Flux vectors
        */
        vvvd fflux_E = FLX::fflux(problem, rho_plr[2], mx_plr[2], my_plr[2], E_plr[2]);
        vvvd fflux_W = FLX::fflux(problem, rho_plr[3], mx_plr[3], my_plr[3], E_plr[3]);        
        vvvd gflux_N = FLX::gflux(problem, rho_plr[0], mx_plr[0], my_plr[0], E_plr[0]);
        vvvd gflux_S = FLX::gflux(problem, rho_plr[1], mx_plr[1], my_plr[1], E_plr[1]);

        /*
            STEP-5 : Anti-diffussion Term
        */
        vvvd adt_x = get_adt_x(rho_plr, mx_plr, my_plr, E_plr, lsp_x, fflux_E, fflux_W);
        vvvd adt_y = get_adt_y(rho_plr, mx_plr, my_plr, E_plr, lsp_y, gflux_N, gflux_S);

        /*
            STEP-6 : Calculate CU Flux
            -> F(U)
            -> G(U)
        */
        vvvd F = get_cu_x(rho_plr, mx_plr, my_plr, E_plr, lsp_x, fflux_E, fflux_W, adt_x);
        vvvd G = get_cu_y(rho_plr, mx_plr, my_plr, E_plr, lsp_y, gflux_N, gflux_S, adt_y);
        
        return {F,G};
    }

    // Function to calculate the CU fluxes F(U) and G(U)
    //! checked manually
    vvvd get_cu_x(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr, vvvd& lsp, vvvd& flux_E, vvvd& flux_W, vvvd& adt){
        // Calculate the CU Numerical flux using all the intermediate terms
        vvvd F;
        
        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();

        vvd F1(row, vd(col));
        vvd F2(row, vd(col));
        vvd F3(row, vd(col));
        vvd F4(row, vd(col));

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col-1 ; j++){

                double diff = lsp[0][i][j] - lsp[1][i][j];
                double prod = lsp[0][i][j]*lsp[1][i][j];
                
                if( diff > CTS::EPSILON ){
                    // use ADT

                    F1[i][j] = ( lsp[0][i][j]*flux_E[0][i][j] - lsp[1][i][j]*flux_W[0][i][j+1] + prod*( rho_plr[3][i][j+1] - rho_plr[2][i][j] - adt[0][i][j] ) ) / (diff);
                    F2[i][j] = ( lsp[0][i][j]*flux_E[1][i][j] - lsp[1][i][j]*flux_W[1][i][j+1] + prod*( mx_plr[3][i][j+1] - mx_plr[2][i][j] - adt[1][i][j] ) ) / (diff);
                    F3[i][j] = ( lsp[0][i][j]*flux_E[2][i][j] - lsp[1][i][j]*flux_W[2][i][j+1] + prod*( my_plr[3][i][j+1] - my_plr[2][i][j] - adt[2][i][j] ) ) / (diff);
                    F4[i][j] = ( lsp[0][i][j]*flux_E[3][i][j] - lsp[1][i][j]*flux_W[3][i][j+1] + prod*( E_plr[3][i][j+1] - E_plr[2][i][j] - adt[3][i][j] ) ) / (diff);

                }
                else{
                    // use normal averaging

                    F1[i][j] = 0.5*(flux_E[0][i][j]+flux_W[0][i][j]);
                    F2[i][j] = 0.5*(flux_E[1][i][j]+flux_W[1][i][j]);
                    F3[i][j] = 0.5*(flux_E[2][i][j]+flux_W[2][i][j]);
                    F4[i][j] = 0.5*(flux_E[3][i][j]+flux_W[3][i][j]);

                }

            }
        }

        F.push_back(F1);
        F.push_back(F2);
        F.push_back(F3);
        F.push_back(F4);

        return F;
    }

    //! checked manually
    vvvd get_cu_y(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr, vvvd& lsp, vvvd& flux_N, vvvd& flux_S, vvvd& adt){
        // Calculate the CU Numerical flux using all the intermediate terms
        vvvd G;
        
        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();

        vvd G1(row, vd(col));
        vvd G2(row, vd(col));
        vvd G3(row, vd(col));
        vvd G4(row, vd(col));

        for(int i=0 ; i<row-1 ; i++){
            for(int j=0 ; j<col ; j++){

                double diff = lsp[0][i][j] - lsp[1][i][j];
                double prod = lsp[0][i][j]*lsp[1][i][j];
                
                if( diff > CTS::EPSILON ){
                    // use ADT

                    G1[i][j] = ( lsp[0][i][j]*flux_N[0][i][j] - lsp[1][i][j]*flux_S[0][i+1][j] + prod*( rho_plr[1][i+1][j] - rho_plr[0][i][j] - adt[0][i][j] ) ) / (diff);
                    G2[i][j] = ( lsp[0][i][j]*flux_N[1][i][j] - lsp[1][i][j]*flux_S[1][i+1][j] + prod*( mx_plr[1][i+1][j] - mx_plr[0][i][j] - adt[1][i][j] ) ) / (diff);
                    G3[i][j] = ( lsp[0][i][j]*flux_N[2][i][j] - lsp[1][i][j]*flux_S[2][i+1][j] + prod*( my_plr[1][i+1][j] - my_plr[0][i][j] - adt[2][i][j] ) ) / (diff);
                    G4[i][j] = ( lsp[0][i][j]*flux_N[3][i][j] - lsp[1][i][j]*flux_S[3][i+1][j] + prod*( E_plr[1][i+1][j] - E_plr[0][i][j] - adt[3][i][j] ) ) / (diff);

                }
                else{
                    // use normal averaging

                    G1[i][j] = 0.5*(flux_N[0][i][j]+flux_S[0][i][j]);
                    G2[i][j] = 0.5*(flux_N[1][i][j]+flux_S[1][i][j]);
                    G3[i][j] = 0.5*(flux_N[2][i][j]+flux_S[2][i][j]);
                    G4[i][j] = 0.5*(flux_N[3][i][j]+flux_S[3][i][j]);

                }

            }
        }

        G.push_back(G1);
        G.push_back(G2);
        G.push_back(G3);
        G.push_back(G4);

        return G;
    }

    // Function to calculate new time step using the CFL conditions
    //! checked manually (confirm later)
    void new_dt(vvvd& lsp_x, vvvd& lsp_y){
        double amax=0; 
        double bmax=0;

        // both lsp_x and lsp_y have the same dimensions but different relevant values
        int row = lsp_x[0].size();
        int col = lsp_x[0][0].size();

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col-1 ; j++){
                amax = max(
                    amax,
                    max(
                        lsp_x[0][i][j], -1*lsp_x[1][j][j]
                    )
                );
            }
        }

        for(int i=0 ; i<row-1 ; i++){
            for(int j=0 ; j<col ; j++){
                bmax = max(
                    bmax,
                    max(
                        lsp_y[0][i][j] , -1*lsp_y[1][i][j]
                    )
                );
            }
        }

        dt = CTS::CFL*min( dx/amax , dy/bmax );

        if( t+dt > time.second ){
            dt = time.second - t;
        }


    }

    // function to obtain the PLR of a given conserved variable
    // used as a subroutine in CU flux calculation
    //! checked manually
    vvvd get_plr(vvd& var){
        vvd slx = get_slx(var);        
        vvd sly = get_sly(var);        

        int row = var.size();
        int col = var[0].size();
        
        vvvd plr; // [N, S, E, W, NE, NW, SE, SW] -> respectively

        vvd n(row-2, vd(col-2));
        vvd s(row-2, vd(col-2));
        vvd e(row-2, vd(col-2));
        vvd w(row-2, vd(col-2));
        vvd ne(row-2, vd(col-2));
        vvd nw(row-2, vd(col-2));
        vvd se(row-2, vd(col-2));
        vvd sw(row-2, vd(col-2));

        for( int i=1 ; i<row-1 ; i++ ){
            for(int j=1 ; j<col-1 ; j++){
                n[i-1][j-1] = var[i][j] + 0.5*sly[i][j];
                s[i-1][j-1] = var[i][j] - 0.5*sly[i][j];
                e[i-1][j-1] = var[i][j] + 0.5*slx[i][j];
                w[i-1][j-1] = var[i][j] - 0.5*slx[i][j];
                ne[i-1][j-1] = var[i][j] + 0.5*slx[i][j] + 0.5*sly[i][j];
                nw[i-1][j-1] = var[i][j] - 0.5*slx[i][j] + 0.5*sly[i][j];
                se[i-1][j-1] = var[i][j] + 0.5*slx[i][j] - 0.5*sly[i][j];
                sw[i-1][j-1] = var[i][j] - 0.5*slx[i][j] - 0.5*sly[i][j];
            }
        }

        // extend cells for the piecewise linear reconstructions
        EXC::extend_matrix(n, bc);
        EXC::extend_matrix(s, bc);
        EXC::extend_matrix(e, bc);
        EXC::extend_matrix(w, bc);
        EXC::extend_matrix(ne, bc);
        EXC::extend_matrix(nw, bc);
        EXC::extend_matrix(se, bc);
        EXC::extend_matrix(sw, bc);

        plr.push_back(n);
        plr.push_back(s);
        plr.push_back(e);
        plr.push_back(w);
        plr.push_back(ne);
        plr.push_back(nw);
        plr.push_back(se);
        plr.push_back(sw);

        return plr;

    }

    // functions to get the slx and sly matrices
    // used as a subroutine in PLR calculation
    //! checked manually
    vvd get_slx(vvd& var){
        int row = var.size();
        int col = var[0].size();

        vvd slx(row, vd(col));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1 ; j<col-1 ; j++){
                double v1 = CTS::THETA * (var[i][j] - var[i][j-1]);
                double v2 = 0.5*(var[i][j+1] - var[i][j-1]);
                double v3 = CTS::THETA * (var[i][j+1] - var[i][j]);

                slx[i][j] = UTL::minmod( v1 , v2 , v3 );
            }
        }

        return slx;
    }

    //! checked manually
    vvd get_sly(vvd& var){
        int row = var.size();
        int col = var[0].size();

        vvd sly(row, vd(col));

        for(int i=1 ; i<row-1 ; i++){
            for(int j=1 ; j<col-1 ; j++){
                double v1 = CTS::THETA * (var[i][j] - var[i-1][j]);
                double v2 = 0.5*(var[i+1][j] - var[i-1][j]);
                double v3 = CTS::THETA * (var[i+1][j] - var[i][j]);

                sly[i][j] = UTL::minmod( v1 , v2 , v3 );
            }
        }

        return sly;
    }

    // function to calculate the One sided local speeds of propagations
    // a+, a-
    //! checked manually
    vvvd get_lsp_a(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr){
        double GAMMA;
        
        if( problem == "RTI" ){
            GAMMA = CTS::GAMMA_RT;
        }
        else{
            GAMMA = CTS::GAMMA;
        }

        // first we get the primitive variable values for East and West boundaries
        vvd uE = PRV::get_u(problem, rho_plr[2], mx_plr[2]);
        vvd uW = PRV::get_u(problem, rho_plr[3], mx_plr[3]);
        vvd vE = PRV::get_v(problem, rho_plr[2], my_plr[2]);
        vvd vW = PRV::get_v(problem, rho_plr[3], my_plr[3]);
        vvd pE = PRV::get_p(problem, rho_plr[2], uE, vE, E_plr[2]);        
        vvd pW = PRV::get_p(problem, rho_plr[3], uW, vW, E_plr[3]);        

        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();

        vvvd lsp;

        vvd ap(row, vd(col));
        vvd am(row, vd(col));

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col-1 ; j++){
                ap[i][j] = max(
                    0.0 ,
                    max(
                        uW[i][j+1] + sqrtf( (GAMMA*pW[i][j+1])/rho_plr[3][i][j+1] ),
                        uE[i][j] + sqrtf( (GAMMA*pE[i][j]) / rho_plr[2][i][j] )
                    )
                );

                am[i][j] = min(
                    0.0,
                    min(
                        uW[i][j+1] - sqrtf( (GAMMA*pW[i][j+1])/rho_plr[3][i][j+1] ),
                        uE[i][j] - sqrtf( (GAMMA*pE[i][j])/rho_plr[2][i][j] )
                    )
                );
            }
        }

        lsp.push_back(ap);
        lsp.push_back(am);

        return lsp;
    }   

    // b+, b-
    //! checked manually
    vvvd get_lsp_b(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr){
        double GAMMA;
        
        if( problem == "RTI" ){
            GAMMA = CTS::GAMMA_RT;
        }
        else{
            GAMMA = CTS::GAMMA;
        }

        // first we get the primitive variable values for East and West boundaries
        vvd uE = PRV::get_u(problem, rho_plr[2], mx_plr[2]);
        vvd uW = PRV::get_u(problem, rho_plr[3], mx_plr[3]);
        vvd vE = PRV::get_v(problem, rho_plr[2], my_plr[2]);
        vvd vW = PRV::get_v(problem, rho_plr[3], my_plr[3]);
        vvd pE = PRV::get_p(problem, rho_plr[2], uE, vE, E_plr[2]);        
        vvd pW = PRV::get_p(problem, rho_plr[3], uW, vW, E_plr[3]);        

        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();

        vvvd lsp;

        vvd bp(row, vd(col));
        vvd bm(row, vd(col));

        for(int i=0 ; i<row-1 ; i++){
            for(int j=0 ; j<col ; j++){
                bp[i][j] = max(
                    0.0 ,
                    max(
                        vW[i+1][j] + sqrtf( (GAMMA*pW[i+1][j])/rho_plr[3][i+1][j] ),
                        vE[i][j] + sqrtf( (GAMMA*pE[i][j]) / rho_plr[2][i][j] )
                    )
                );

                bm[i][j] = min(
                    0.0,
                    min(
                        vW[i+1][j] - sqrtf( (GAMMA*pW[i+1][j])/rho_plr[3][i+1][j] ),
                        vE[i][j] - sqrtf( (GAMMA*pE[i][j])/rho_plr[2][i][j] )
                    )
                );
            }
        }

        lsp.push_back(bp);
        lsp.push_back(bm);

        return lsp;
    }

    // Functions to calculate the Anti-diffusion terms
    // rho_adt, mx_adt, my_adt, E_adt
    //! checked manually
    vvvd get_adt_x(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr, vvvd& lsp, vvvd& flux_E, vvvd& flux_W){
        vvd rho_star = get_star_x(rho_plr, lsp, flux_E[0], flux_W[0]);
        vvd mx_star = get_star_x(mx_plr, lsp, flux_E[1], flux_W[1]);
        vvd my_star = get_star_x(my_plr, lsp, flux_E[2], flux_W[2]);
        vvd E_star = get_star_x(E_plr, lsp, flux_E[3], flux_W[3]);

        // calculating anti diffusion terms
        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();
        
        vvvd adt;

        vvd rho_adt(row, vd(col));
        vvd mx_adt(row, vd(col));
        vvd my_adt(row, vd(col));
        vvd E_adt(row, vd(col));

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col-1 ; j++){
                rho_adt[i][j] = UTL::minmod(
                    rho_plr[7][i][j+1]-rho_star[i][j],
                    rho_star[i][j]-rho_plr[6][i][j],
                    rho_plr[5][i][j+1]-rho_star[i][j],
                    rho_star[i][j] - rho_plr[4][i][j]
                );

                mx_adt[i][j] = UTL::minmod(
                    mx_plr[7][i][j+1]-mx_star[i][j],
                    mx_star[i][j]-mx_plr[6][i][j],
                    mx_plr[5][i][j+1]-mx_star[i][j],
                    mx_star[i][j] - mx_plr[4][i][j]
                );

                my_adt[i][j] = UTL::minmod(
                    my_plr[7][i][j+1]-my_star[i][j],
                    my_star[i][j]-my_plr[6][i][j],
                    my_plr[5][i][j+1]-my_star[i][j],
                    my_star[i][j] - my_plr[4][i][j]
                );

                E_adt[i][j] = UTL::minmod(
                    E_plr[7][i][j+1]-E_star[i][j],
                    E_star[i][j]-E_plr[6][i][j],
                    E_plr[5][i][j+1]-E_star[i][j],
                    E_star[i][j] - E_plr[4][i][j]
                );
            }
        }

        adt.push_back(rho_adt);
        adt.push_back(mx_adt);
        adt.push_back(my_adt);
        adt.push_back(E_adt);

        return adt;
    }

    //! checked manually
    vvvd get_adt_y(vvvd& rho_plr, vvvd& mx_plr, vvvd& my_plr, vvvd& E_plr, vvvd& lsp, vvvd& flux_N, vvvd& flux_S){
        vvd rho_star = get_star_y(rho_plr, lsp, flux_N[0], flux_S[0]);
        vvd mx_star = get_star_y(mx_plr, lsp, flux_N[1], flux_S[1]);
        vvd my_star = get_star_y(my_plr, lsp, flux_N[2], flux_S[2]);
        vvd E_star = get_star_y(E_plr, lsp, flux_N[3], flux_S[3]);

        // calculating anti diffusion terms
        int row = rho_plr[0].size();
        int col = rho_plr[0][0].size();
        
        vvvd adt;

        vvd rho_adt(row, vd(col));
        vvd mx_adt(row, vd(col));
        vvd my_adt(row, vd(col));
        vvd E_adt(row, vd(col));

        for(int i=0 ; i<row-1 ; i++){
            for(int j=0 ; j<col ; j++){
                rho_adt[i][j] = UTL::minmod(
                    rho_plr[7][i+1][j]-rho_star[i][j],
                    rho_star[i][j]-rho_plr[5][i][j],
                    rho_plr[6][i+1][j]-rho_star[i][j],
                    rho_star[i][j] - rho_plr[4][i][j]
                );

                mx_adt[i][j] = UTL::minmod(
                    mx_plr[7][i+1][j]-mx_star[i][j],
                    mx_star[i][j]-mx_plr[5][i][j],
                    mx_plr[6][i+1][j]-mx_star[i][j],
                    mx_star[i][j] - mx_plr[4][i][j]
                );

                my_adt[i][j] = UTL::minmod(
                    my_plr[7][i+1][j]-my_star[i][j],
                    my_star[i][j]-my_plr[5][i][j],
                    my_plr[6][i+1][j]-my_star[i][j],
                    my_star[i][j] - my_plr[4][i][j]
                );

                E_adt[i][j] = UTL::minmod(
                    E_plr[7][i+1][j]-E_star[i][j],
                    E_star[i][j]-E_plr[5][i][j],
                    E_plr[6][i+1][j]-E_star[i][j],
                    E_star[i][j] - E_plr[4][i][j]
                );
            }
        }

        adt.push_back(rho_adt);
        adt.push_back(mx_adt);
        adt.push_back(my_adt);
        adt.push_back(E_adt);

        return adt;
    }

    // Functions to calculate the Intermediate star terms used in calculation for the Anti-diffusion terms
    //! checked manually
    vvd get_star_x(vvvd& plr, vvvd& lsp, vvd& flux_E, vvd& flux_W){
        int row = plr[0].size();
        int col = plr[0][0].size();

        vvd star(row, vd(col));

        for(int i=0 ; i<row ; i++){
            for(int j=0 ; j<col-1 ; j++){
                star[i][j] = (lsp[0][i][j]*plr[3][i][j+1] - lsp[1][i][j]*plr[2][i][j] - flux_W[i][j+1] + flux_E[i][j]) / (lsp[0][i][j]-lsp[1][i][j]);
            }
        }

        return star;
    }

    //! checked manually
    vvd get_star_y(vvvd& plr, vvvd& lsp, vvd& flux_N, vvd& flux_S){
        int row = plr[0].size();
        int col = plr[0][0].size();

        vvd star(row, vd(col));

        for(int i=0 ; i<row-1 ; i++){
            for(int j=0 ; j<col ; j++){
                star[i][j] = (lsp[0][i][j]*plr[1][i+1][j] - lsp[1][i][j]*plr[0][i][j] - flux_S[i+1][j] + flux_N[i][j]) / (lsp[0][i][j]-lsp[1][i][j]);
            }
        }

        return star;
    }

};

