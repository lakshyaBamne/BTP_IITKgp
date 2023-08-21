/*
    Namespace to output the values of the conserved variables in a text file
    to be printed later
*/
#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>

using namespace std;

namespace OutputResult{
    // Function which takes a vector, a file name and appends the contents
    // of the vector at the end of the file
    void write_vector(vector<double>& nums, string file_name);

    string vector_to_string(vector<double>& nums);
}

void OutputResult::write_vector(vector<double>& nums, string file_name){

    string output = vector_to_string( nums );

    // open the file in append mode
    ofstream fout;
    fout.open(file_name, ios::app);

    if( fout ){
        fout << output << endl;
    }
    else{
        cout << "---ERROR--- Could not write to file ---" << endl;
    }

    fout.close();
}

string OutputResult::vector_to_string(vector<double>& nums){
    stringstream ss;

    for(auto i : nums){
        ss << i << " ";
    }

    return ss.str();
}

