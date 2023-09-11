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

using namespace std;

namespace EXC = ExtendCells;
namespace OPR = OutputResult;

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

        // Run CU Numerical Scheme
        
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

                dx = 3/2500;
                dy = 3/2500;

                time.first = 0;
                time.second = 1;

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

                dx = 3/800;
                dy = 3/800;

                time.first = 0;
                time.second = 3.2;

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

                dx = 3/6000;
                dy = 3/6000;

                time.first = 0;
                time.second = 2.5;

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

                dx = 3/1024;
                dy = 3/1024;

                time.first = 0;
                time.second = 1; // also evaluate results at t=2.5 and t=4

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

                dx = 1/1024;
                dy = 1/1024;

                time.first = 0;
                time.second = 2.95; // also evaluate results at t=2.5 and t=4

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

                dx = 3/10000;
                dy = 3/10000;

                time.first = 0;
                time.second = 1;

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

                dx = 3/3200;
                dy = 3/3200;

                time.first = 0;
                time.second = 3.2;

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

                dx = 3/24000;
                dy = 3/24000;

                time.first = 0;
                time.second = 2.5;

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

                dx = 3/4096;
                dy = 3/4096;

                time.first = 0;
                time.second = 1; // also evaluate results at t=2.5 and t=4

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

                dx = 1/4096;
                dy = 1/4096;

                time.first = 0;
                time.second = 2.95; // also evaluate results at t=2.5 and t=4

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

        cout << "---LOG--- Initialized Riemann problem variables ---" << endl;

        // initialize grid vectors
        make_grids();

        cout << "---LOG--- Initialized computational grid ---" << endl;

        OPR::write_grids(mode, gridX, gridY);

        // initialize variable matrices
        init_vars();

        cout << "---LOG--- Initialized Variable matrices ---" << endl;

        OPR::write_matrix(mode, "Density", rho);
        OPR::write_matrix(mode, "MomentumX", mx);
        OPR::write_matrix(mode, "MomentumY", my);
        OPR::write_matrix(mode, "Energy", E);

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

                for(auto i : gridX){
                    bool c1 = i>-0.1 && i<0.1 && j>0 && j<0.02;
                    bool c2 = i>-0.02 && i<0.02 && j>0.02 && j<0.1;
                    bool c3 = pow(i+0.02,2) + pow(j-0.02,2) < pow(0.08,2);
                    bool c4 = pow(i-0.02,2) + pow(j-0.02,2) < pow(0.08,2);

                    if( c1 || c2 || c3 || c4 ){
                        rho1.push_back( 1.4 );
                        mx1.push_back( 0.0 );
                        my1.push_back( 0.2 );
                        E1.push_back( 1.0 );
                    }
                    else{
                        rho1.push_back( 1.0 );
                        mx1.push_back( 0.0 );
                        my1.push_back( 0.2 );
                        E1.push_back( 1.0 );
                    }
                }

                rho.push_back(rho1);
                mx.push_back(mx1);
                my.push_back(my1);
                E.push_back(E1);
            }
        }
        // else if( mode == "2DR" ){

        // }
        // else if( mode == "EXP" ){

        // }
        // else if( mode == "IMP" ){

        // }
        // else if( mode == "KHI" ){

        // }
        // else if( mode == "RTI" ){

        // }
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

};

