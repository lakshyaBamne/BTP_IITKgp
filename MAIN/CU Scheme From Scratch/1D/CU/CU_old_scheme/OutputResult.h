/*
    ! @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-8 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme - with RH optimization for Euler Equations of Gas Dynamics (CURH Scheme)
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#pragma once

#include<iostream>
#include<vector>
#include<sstream>
#include<fstream>
#include<string>
#include<cmath>

#include "Constants.h"

using namespace std;

namespace cts = Constants;

namespace OutputResult{

    // function to write the given grid vectors to file ComputationalGrid.txt
    void write_grids(int mode, const vector<double>& grid);
    string vector_to_string(const vector<double>& grid);

    // function to write a given 2D matrix to it's corresponding data output file
    void write_conserved_variables(int mode, const vector<double>& U);
}

/*
    Function implementations
*/

string OutputResult::vector_to_string(const vector<double>& grid){
    stringstream ss;

    for(int i=1 ; i<=cts::n ; i++){
        ss << grid[i] << " ";
    }

    return ss.str();
}

void OutputResult::write_grids(int mode, const vector<double>& grid){
    string file_name;

    if( mode == 1 ){
        // CU plot
        file_name = "result1/ComputationalGrid.txt";
    }
    else{
        // Reference plot
        file_name = "result2/ComputationalGrid.txt";
    }

    // convert the vectors to strings
    string strx = vector_to_string(grid);

    // open required file in append mode
    ofstream fout;
    fout.open(file_name, ios::app);

    // write to file
    fout << strx << endl;

    cout << "---LOG--- write successfull in [" << file_name << "] ---" << endl; 

    fout.close();
}



void OutputResult::write_conserved_variables(int mode, const vector<double>& U){

    string file_density;

    if( mode == 1 ){
        file_density = "result1/density.txt";
    }
    else{
        file_density = "result2/density.txt";
    }

    // open the file
    ofstream fout;
    fout.open(file_density, ios::app);

    // write
    stringstream ss;
    for(int i=1 ; i<=cts::n ; i++){
        ss << U[i] << " ";
    }

    fout << ss.str() << endl;

    cout << "---LOG--- write successfull in [" << file_density << "]" << endl;

    fout.close();

}
