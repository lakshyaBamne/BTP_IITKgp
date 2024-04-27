/*
    ! @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-8 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme - with RH optimization for Euler Equations of Gas Dynamics (CURH Scheme)
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>
#include<utility>

using namespace std;

namespace ExtendCells{

    // function to extend the cells in the 2D CU scheme
    // since different boundaries can have different BC
    // we have to make a function which can handle it
    void extend_matrix( vector< vector<double> >& matrix , unordered_map<string,string>& bc , string var );

    // Function to extend the PLR according to the boundary conditions
    void extend_plr(vector< vector< vector<double> > >& plr, unordered_map<string,string>& bc, string var);
}

void ExtendCells::extend_matrix( vector< vector<double> >& matrix , unordered_map<string,string>& bc , string var ){

    // 1) Extending the North(top) and South(bottom) boundaries

    //! North boundary
    if( bc["N"] == "FREE" ){
        matrix.push_back( matrix.back() );
    }
    else if( bc["N"] == "REF" ){

        matrix.push_back( matrix.back() );
        
        if( var == "MomentumY" ){
            for(auto& i : matrix.back()){
                i = -i;
            }
        }

    }
    else if( bc["N"] == "PER" ){
        
    }
    else if( bc["N"] == "DIR" ){
        
    }
    else if( bc["N"] == "SOL" ){
        
    }
    else{
        cout << "---ERROR--- Please give correct B.C. for North boundary---" << endl;
    }

    //! South boundary
    if( bc["S"] == "FREE" ){
        matrix.insert(matrix.begin(), matrix[0]);
    }
    else if( bc["S"] == "REF" ){
        
        matrix.insert( matrix.begin() , matrix[0] );
        
        if( var == "MomentumY" ){
            for(auto& i : matrix[0]){
                i = -i;
            }
        }

    }
    else if( bc["S"] == "PER" ){
        
    }
    else if( bc["S"] == "DIR" ){
        
    }
    else if( bc["N"] == "SOL" ){
        
    }
    else{
        cout << "---ERROR--- Please give correct B.C. for North boundary---" << endl;
    }

    // 2) Extending the West(left) and East(right) boundaries

    //! East Boundary
    if( bc["E"] == "FREE" ){
        for(auto& i : matrix){
            i.push_back( i.back() );
        }
    }
    else if( bc["E"] == "REF" ){

        if( var == "MomentumX" ){
            for(auto& i : matrix){
                i.push_back( -1*i.back() );
            }
        }
        else{
            for(auto& i : matrix){
                i.push_back( i.back() );
            }
        }

    }
    else if( bc["E"] == "PER" ){

    }
    else if( bc["E"] == "DIR" ){
        
    }
    else if( bc["N"] == "SOL" ){

    }
    else{
        cout << "---ERROR--- Please give correct B.C. for North boundary---" << endl;
    }

    //! West Boundary
    if( bc["W"] == "FREE" ){
        for(auto& i : matrix){
            i.insert( i.begin() , i[0] );
        }
    }
    else if( bc["W"] == "REF" ){

        if( var == "MomentumX" ){
            for(auto& i : matrix){
                i.insert( i.begin() , -1*i[0] );
            }   
        }
        else{
            for(auto& i : matrix){
                i.insert( i.begin() , i[0] );
            }
        }

    }
    else if( bc["W"] == "PER" ){

    }
    else if( bc["W"] == "DIR" ){
        
    }
    else if( bc["N"] == "SOL" ){

    }
    else{
        cout << "---ERROR--- Please give correct B.C. for North boundary---" << endl;
    }

}

void ExtendCells::extend_plr(vector< vector< vector<double> > >& plr, unordered_map<string,string>& bc, string var){
    int row = plr[0].size();
    int col = plr[0][0].size();

    //! Free boundary conditions on all the boundaries
    // -> MCW, 2DR-CFG3
    if( bc["N"]=="FREE" && bc["S"]=="FREE" && bc["E"]=="FREE" && bc["W"]=="FREE" ){
        
        for(int i=1 ; i<row-1 ; i++){

            plr[2][i][0] = plr[3][i][1]; // E-W
            plr[3][i][col-1] = plr[2][i][col-2]; // W-E

        }

        for(int i=1 ; i<col-1 ; i++){

            plr[0][0][i] = plr[1][1][i]; // N-S
            plr[1][row-1][i] = plr[0][row-2][i]; // S-N

        }

    }
    //! Reflectfree boundary conditions
    // -> EXP
    else if( bc["N"]=="FREE" && bc["S"]=="REF" && bc["E"]=="FREE" && bc["W"]=="REF" ){

        if( var == "MomentumX" ){
            for(int i=1 ; i<row-1 ; i++){

                plr[2][i][0] = -plr[3][i][1]; // E-W
                plr[3][i][col-1] = plr[2][i][col-2]; // W-E

            }

            for(int i=1 ; i<col-1 ; i++){

                plr[0][0][i] = plr[1][1][i]; // N-S
                plr[1][row-1][i] = plr[0][row-2][i]; // S-N

            }
        }
        else if( var == "MomentumY" ){
            for(int i=1 ; i<row-1 ; i++){

                plr[2][i][0] = plr[3][i][1]; // E-W
                plr[3][i][col-1] = plr[2][i][col-2]; // W-E

            }

            for(int i=1 ; i<col-1 ; i++){

                plr[0][0][i] = -plr[1][1][i]; // N-S
                plr[1][row-1][i] = plr[0][row-2][i]; // S-N

            }
        }
        else{ // Density and Energy are same as that in Free conditions
            for(int i=1 ; i<row-1 ; i++){

                plr[2][i][0] = plr[3][i][1]; // E-W
                plr[3][i][col-1] = plr[2][i][col-2]; // W-E

            }

            for(int i=1 ; i<col-1 ; i++){

                plr[0][0][i] = plr[1][1][i]; // N-S
                plr[1][row-1][i] = plr[0][row-2][i]; // S-N

            }
        }


    }   

}


