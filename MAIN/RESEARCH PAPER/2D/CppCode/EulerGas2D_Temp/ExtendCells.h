/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

/*
    Functions to extend cells for a matrix
    -> [FREE] FREE
    -> [REF] REFLECTIVE
    -> [PER] PERIODIC
    -> [DIR] DIRICHLET
*/

#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>

using namespace std;

namespace ExtendCells{
    
    // function to extend the cells in the 2D CU scheme
    // since different boundaries can have different BC
    // we have to make a function which can handle it
    void extend_matrix( vector< vector<double> >& matrix , unordered_map<string,string>& bc );
    void extend_matrix( vector< vector<double> >& matrix , unordered_map<string,string>& bc , string var );

    // function to extend the 1D grids
    void extend_grid( vector<double>& grid , string bc );
}

/*
    Function implementations
*/

void ExtendCells::extend_grid( vector<double>& grid , string bc ){
    
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

void ExtendCells::extend_matrix( vector< vector<double> >& matrix , unordered_map<string,string>& bc ){

    // 1) Extending the North(top) and South(bottom) boundaries

    //! North boundary
    if( bc["N"] == "FREE" ){
        matrix.push_back( matrix.back() );
    }
    else if( bc["N"] == "REF" ){

    }
    else if( bc["N"] == "PER" ){
        matrix.push_back( matrix[0] );
    }
    else if( bc["N"] == "DIR" ){
        
    }
    else{
        cout << "---ERROR--- Please give correct B.C. for North boundary---" << endl;
    }

    //! South boundary
    if( bc["S"] == "FREE" ){
        matrix.insert(matrix.begin(), matrix[0]);
    }
    else if( bc["S"] == "REF" ){

    }
    else if( bc["S"] == "PER" ){
        matrix.insert( matrix.begin(), matrix[ matrix.size()-2 ] );
    }
    else if( bc["S"] == "DIR" ){
        
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

    }
    else if( bc["E"] == "PER" ){
        for(auto& i : matrix){
            i.push_back( i[0] );
        }
    }
    else if( bc["E"] == "DIR" ){
        
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

    }
    else if( bc["W"] == "PER" ){
        for(auto& i : matrix){
            i.insert( i.begin() , i[i.size()-2] );
        }
    }
    else if( bc["W"] == "DIR" ){
        
    }
    else{
        cout << "---ERROR--- Please give correct B.C. for North boundary---" << endl;
    }

}

void ExtendCells::extend_matrix( vector< vector<double> >& matrix , unordered_map<string,string>& bc , string var ){

    // 1) Extending the North(top) and South(bottom) boundaries

    //! North boundary
    if( bc["N"] == "FREE" ){
        matrix.push_back( matrix.back() );
    }
    else if( bc["N"] == "REF" ){

        if( var=="MomentumY" ){
            auto temp = matrix.back();
            for(auto& k : temp){
                k = -1*k;
            }
            matrix.push_back(temp);
        }
        else{
            matrix.push_back( matrix.back() );
        }

    }
    else if( bc["N"] == "PER" ){
        matrix.push_back( matrix[0] );
    }
    else if( bc["N"] == "DIR" ){
        
    }
    else if( bc["N"] == "SOL" ){
        if( var=="MomentumY" ){
            auto temp = matrix.back();
            for(auto& k : temp){
                k = -1*k;
            }
            matrix.push_back(temp);
        }
        else{
            matrix.push_back( matrix.back() );
        }
    }
    else{
        cout << "---ERROR--- Please give correct B.C. for North boundary---" << endl;
    }

    //! South boundary
    if( bc["S"] == "FREE" ){
        matrix.insert(matrix.begin(), matrix[0]);
    }
    else if( bc["S"] == "REF" ){
        if( var=="MomentumY" ){
            auto temp = matrix[0];
            for(auto& k : temp){
                k = -1*k;
            }
            matrix.insert(matrix.begin(), temp);
        }
        else{
            matrix.insert(matrix.begin(), matrix[0]);
        }
    }
    else if( bc["S"] == "PER" ){
        matrix.insert( matrix.begin(), matrix[ matrix.size()-2 ] );
    }
    else if( bc["S"] == "DIR" ){
        
    }
    else if( bc["N"] == "SOL" ){
        if( var=="MomentumY" ){
            auto temp = matrix[0];
            for(auto& k : temp){
                k = -1*k;
            }
            matrix.insert(matrix.begin(), temp);
        }
        else{
            matrix.insert(matrix.begin(), matrix[0]);
        }
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
        if( var=="MomentumX" ){
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
        for(auto& i : matrix){
            i.push_back( i[0] );
        }
    }
    else if( bc["E"] == "DIR" ){
        
    }
    else if( bc["N"] == "SOL" ){
        if( var=="MomentumX" ){
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
        if( var=="MomentumX" ){
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
        for(auto& i : matrix){
            i.insert( i.begin() , i[i.size()-2] );
        }
    }
    else if( bc["W"] == "DIR" ){
        
    }
    else if( bc["N"] == "SOL" ){
        if( var=="MomentumX" ){
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
    else{
        cout << "---ERROR--- Please give correct B.C. for North boundary---" << endl;
    }

}
