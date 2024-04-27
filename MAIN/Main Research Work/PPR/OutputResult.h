/*
    ! Namespace contains the functions required to output the required files
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>

#include "Constants.h"

using namespace std;

namespace cts = Constants;

namespace OutputResult{
    void export_env(string path);
    void export_env(string path1, string path2);

    // function to convert a vector to a string
    string vector_to_string(const vector<double>& vec);

    // write a vector to an external file
    void output_vector(const vector<double>& vec, string path);

    // write a matrix to an external file
    void output_matrix(const vector< vector<double> >& vec, string path);
}

void OutputResult::export_env(string path){
    ofstream fout;
    fout.open(path, ios::app);

    fout << cts::PROBLEM << endl; 
    fout << cts::BC << endl;
    fout << cts::APPROX << endl;
    fout << cts::FLUX << endl;
    fout << cts::n << endl;
    fout << cts::domx.first << " " << cts::domx.second << endl;

    fout.close();
}

void OutputResult::export_env(string path1, string path2){
    // env1
    ofstream fout;
    fout.open(path1, ios::app);

    fout << cts::PROBLEM << " " << cts::BC << endl;
    fout << cts::APPROX << " " << cts::FLUX << endl;
    fout << cts::n << endl;
    fout << cts::domx.first << " " << cts::domx.second << endl;

    fout.close();

    // env2
    fout;
    fout.open(path2, ios::app);

    fout << cts::PROBLEM << " " << cts::BC << endl;
    fout << cts::APPROX << " " << cts::FLUX << endl;
    fout << cts::nref << endl;
    fout << cts::domx.first << " " << cts::domx.second << endl;

    fout.close();
}


string OutputResult::vector_to_string(const vector<double>& vec){
    stringstream ss;
    for(auto v : vec){
        ss << v << " ";
    }

    return ss.str();
}

void OutputResult::output_vector(const vector<double>& vec, string path){
    // first we need to find the string representation of the vector
    string str = vector_to_string(vec);

    // open the file in append mode
    ofstream fout;
    fout.open(path, ios::app);

    fout << str << endl;

    fout.close();

    cout << "---LOG--- output written to " << path << " ---" << endl;
}

void OutputResult::output_matrix(const vector< vector<double> >& vec, string path){
    int row = vec.size();
    int col = vec[0].size();

    for(int j=col-1 ; j>=0 ; j--){
        // first we need to extract one row out of the grid
        vector<double> one_row;
        for(int i=0 ; i<row ; i++){
            one_row.push_back( vec[i][j] );
        }

        // now output this row
        output_vector(one_row, path);
    }

    cout << "---LOG--- output written to " << path << " ---" << endl;
}

