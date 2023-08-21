/*
    Namespace to encapsulate all the Utility functions
*/
#pragma once

#include<iostream>
#include<vector>

using namespace std;

// Namespace definition
namespace Utility{
    // minmod functions with different numbers of arguments
    double minmod(double x, double y);
    double minmod(double x, double y, double z);
}

/*
    Function implementations for the associated functions
*/
double Utility::minmod(double x, double y){
    if( x>0 && y>0 ){
        return min(x,y);
    }
    else if( x<0 && y<0 ){
        return max(x,y);
    }
    else{
        return 0;
    }
}

double Utility::minmod(double x, double y, double z){
    if( x>0 && y>0 && z>0 ){
        return min( x, min(y,z) );
    }
    else if( x<0 && y<0 && z<0 ){
        return max( x, max(y,z) );
    }
    else{
        return 0;
    }
}

