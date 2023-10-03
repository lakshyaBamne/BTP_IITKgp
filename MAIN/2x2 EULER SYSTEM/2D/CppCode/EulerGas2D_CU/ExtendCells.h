/*
    * @author - lakshya Bamne (20MA20029)
    * @supervisor - Prof. Naveen Kumar Garg

    ! Implementation of CU scheme for modified Flux functions in Euler's Equations
*/

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>

using namespace std;

#define vd vector<double>
#define vvd vector< vector<double> >
#define umap_ss unordered_map<string,string>

namespace ExtendCells{

    // function to extend the cells of a matrix on different boundaries
    void extend_matrix(vvd& matrix, umap_ss bc, string var);
}

// function to extend the cells of a matrix on different boundaries
void ExtendCells::extend_matrix(vvd& matrix, umap_ss bc, string var){
    // North boundary
    if( bc["N"] == "FREE" ){
        matrix.push_back( matrix.back() );
    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions for North boundary ---" << endl;
    }

    // South boundary
    if( bc["S"] == "FREE" ){
        matrix.insert( matrix.begin() , matrix[0] );
    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions for South boundary ---" << endl;
    }

    // East boundary
    if( bc["E"] == "FREE" ){
        for(auto& i : matrix){
            i.push_back( i.back() );
        }
    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions for East boundary ---" << endl;
    }

    // West boundary
    if( bc["W"] == "FREE" ){
        for(auto& i : matrix){
            i.insert( i.begin() , i[0] );
        }
    }
    else{
        cout << "---ERROR--- Please enter correct boundary conditions for West boundary ---" << endl;
    }
}
