/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

/*
    Functions to extend cells for a matrix
    -> [FREE] FREE
    -> [RFREE] REFLECTFREE
    -> [SOL] SOLIDWALL
    -> [PER] PERIODIC
    -> [DIR] DIRICHLET
*/

#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>

using namespace std;

#define vd vector<double>
#define vvd vector< vector<double> >
#define vvvd vector< vector< vector<double> > >
#define umap_ss unordered_map<string,string>

namespace ExtendCells{
    
    // function to extend the cells in the 2D CU scheme
    // since different boundaries can have different BC
    // we have to make a function which can handle it
    void extend_matrix( vvd& matrix , string bc );
    void extend_matrix( vvd& matrix , string bc , string var );

    // function to extend the 1D grids
    void extend_grid( vd& grid , string bc );

    void extend_plr( vvvd& plr , string bc , string var );
}

/*
    Function implementations
*/

void ExtendCells::extend_grid( vd& grid , string bc ){
    
    if( bc == "FREE" ){
        grid.push_back( grid.back() );
        grid.insert( grid.begin() , grid[0] );
    }
    else if( bc == "REF" ){

    }
    else if( bc == "PER" ){
        
    }
    else if( bc == "DIR" ){

    }
    else{
        cout << "---ERROR---  ---" << endl;
    }

}

void ExtendCells::extend_matrix( vvd& matrix , string bc ){

    if( bc == "FREE" ){
        matrix.push_back( matrix.back() );
        matrix.insert( matrix.begin() , matrix[0] );

        for(auto& i : matrix){
            i.push_back( i.back() );
            i.insert( i.begin(), i[0] );
        }
    }
    else if( bc == "RFREE" ){

    }
    else if( bc == "SOLIDWALL" ){
        
    }
    else if( bc == "DIRICHLET" ){

    }
    else if( bc == "PERIODIC" ){

    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions---" << endl;
    }

}

void ExtendCells::extend_matrix( vvd& matrix , string bc , string var ){

    if( bc == "FREE" ){
        matrix.push_back( matrix.back() );
        matrix.insert( matrix.begin() , matrix[0] );

        for(auto& i : matrix){
            i.push_back( i.back() );
            i.insert( i.begin(), i[0] );
        }
    }
    else if( bc == "RFREE" ){

    }
    else if( bc == "SOLIDWALL" ){
        
    }
    else if( bc == "DIRICHLET" ){

    }
    else if( bc == "PERIODIC" ){

    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions---" << endl;
    }
}

void ExtendCells::extend_plr( vvvd& plr , string bc , string var ){

    int row = plr[0].size();
    int col = plr[0][0].size();

    //! Free boundaries on all sides
    if( bc == "FREE" ){

        // E, W, NE, NW, SE, SW
        for(int i=1 ; i<row-1 ; i++){
            plr[2][i][0] = plr[3][i][1];
            plr[3][i][col-1] = plr[2][i][col-2];
        }

        // N, S, SW, SE, NE, NW
        for(int i=1 ; i<col-1 ; i++){
            plr[0][0][i] = plr[1][1][i];
            plr[1][col-1][i] = plr[0][col-2][i];
        }

    }
    else{
        cout << "[ERROR] Fix PLR extension other than FREE Conditions" << endl;
    }

}
