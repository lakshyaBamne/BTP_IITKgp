/*
    Function Implementations for the Utility functions used in the simulation
*/

#pragma once

#include<iostream>
#include<vector>
#include<string>

using namespace std;

namespace Utility{
    double minmod2(double a, double b);
    double minmod3(double a, double b, double c);
}

double Utility::minmod2(double a, double b){
    if( min(a,b) > 0 ){
        return min(a,b);
    }
    else if( max(a,b) < 0 ){
        return max(a,b);
    }
    else{
        return 0;
    }
}

double Utility::minmod3(double a, double b, double c){
    if( min(a,min(b,c)) > 0 ){
        return min(a,min(b,c));
    }
    else if( max(a,max(b,c)) < 0 ){
        return max(a,max(b,c));
    }
    else{
        return 0;
    }
}

