#pragma once

#include<iostream>
#include<string>
#include<fstream>
#include<sstream>

#include "Constants.h"

using namespace std;

namespace cts = Constants;

namespace Utility{
    /*
        ! Nonlinear minmod limiter used in calculating slopes for the linear approximations
    */
    double minmod(double a, double b, double c);

    /*
        ! Function to export the run environment for the running simulation for plotting later
    */
    void export_env();
}

double Utility::minmod(double a, double b, double c){
    if( min( a , min(b,c) ) > 0 ){
        return min(a,min(b,c));
    }
    else if( max(a,max(b,c)) < 0 ){
        return max(a,max(b,c));
    }
    else{
        return 0;
    }
}

void Utility::export_env(){
    string file_name;

    if( cts::MODE == 1 ){
        file_name = "env1/Environment.txt";
    }
    else{
        file_name = "env2/Environment.txt";
    }
    
    stringstream ss;
    
    ofstream fout;
    fout.open(file_name , ios::app);

    ss << cts::MODE << endl;
    ss << cts::PROBLEM << endl;
    ss << cts::dom.first << " " << cts::dom.second << endl;
    ss << cts::time.first << " " << cts::time.second << endl;
    ss << cts::n << endl;

    fout << ss.str() << endl;

    fout.close();
}
