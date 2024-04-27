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

using namespace std;

namespace OutputResult{

    // function to write the given grid vectors to file ComputationalGrid.txt
    void write_grids(int mode, const vector<double>& gridX, const vector<double>& gridY);
    string vector_to_string(const vector<double>& grid);

    // function to write a given 2D matrix to it's corresponding data output file
    void write_matrix(int mode, string var, const vector< vector<double> >& matrix);
}

/*
    Function implementations
*/

string OutputResult::vector_to_string(const vector<double>& grid){
    stringstream ss;

    for(auto i : grid){
        ss << i << " ";
    }

    return ss.str();
}

void OutputResult::write_grids(int mode, const vector<double>& gridX, const vector<double>& gridY){
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
    string strx = vector_to_string(gridX);
    string stry = vector_to_string(gridY);

    // open required file in append mode
    ofstream fout;
    fout.open(file_name, ios::app);

    // write to file
    fout << strx << endl;
    fout << stry << endl;

    cout << "---LOG--- write successfull in [" << file_name << "] ---" << endl; 

    fout.close();
}

void OutputResult::write_matrix(int mode, string var, const vector< vector<double> >& matrix){

    string file_name;

    if( mode == 1 ){
        file_name = "result1/" + var + ".txt";
    }
    else{
        file_name = "result2/" + var + ".txt";
    }

    // open the file
    ofstream fout;
    fout.open(file_name, ios::app);

    // write
    for( auto j : matrix ){
        string str = vector_to_string(j);

        fout << str << endl;
    }

    cout << "---LOG--- write successfull in [" << file_name << "]" << endl;

    fout.close();

}
