/*
    * @author : Lakshya Bamne(20MA20029) student @IIT Kharagpur, semester-7 (Mathematics and Computing)
    * 2-Dimensional Central Upwind Scheme for Euler Equations of Gas Dynamics
    * Supervisor - Prof. Naveen Kumar Garg (IIT Kharagpur, Dept. of Mathematics)
*/

#pragma once

#include<iostream>
#include<vector>
#include<utility>
#include<unordered_map>
#include<string>

#define vd vector<double>
#define vvd vector< vector<double> >
#define pdd pair<double,double>

using namespace std;

namespace Utility{
    // minmod functions with 2,3,4 arguments respectively
    double minmod(double a, double b);
    double minmod(double a, double b, double c);
    double minmod(double a, double b, double c, double d);
}

/*
    Function implementations
*/

double Utility::minmod(double a, double b){
    if( min(a,b) > 0 ){
        return min(a,b);
    }
    else if( max(a,b) < 0 ){
        return 0;
    }
    else{
        return 0;
    }
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

double Utility::minmod(double a, double b, double c, double d){
    if( min(a,min(b,min(c,d))) > 0 ){
        return min(a,min(b,min(c,d)));
    }
    else if( max(a,max(b,max(c,d))) < 0 ){
        return max(a,max(b,max(c,d)));
    }
    else{
        return 0;
    }
}

