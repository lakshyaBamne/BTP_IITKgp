/*
    * @author - lakshya Bamne (20MA20029)
    * @supervisor - Prof. Naveen Kumar Garg

    ! Implementation of CU scheme for modified Flux functions in Euler's Equations
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
    void write_grids(int shock_type, const vector<double>& gridX, const vector<double>& gridY);
    string vector_to_string(const vector<double>& grid);

    // function to write a given 2D matrix to it's corresponding data output file
    void write_matrix(int shock_type, string var, const vector< vector<double> >& matrix);
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

void OutputResult::write_grids(int shock_type, const vector<double>& gridX, const vector<double>& gridY){
    stringstream ss;
    ss << "result/" << shock_type << "_ComputationalGrid.txt";

    string file_name = ss.str();

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

void OutputResult::write_matrix(int shock_type, string var, const vector< vector<double> >& matrix){
    stringstream ss;
    ss << "result/" << shock_type << "_" << var << ".txt";

    string file_name = ss.str();

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
