/*
    Function definitions for the Class Domain
*/
#include<iostream>
#include<vector>

#define ll long long int

#include "Domain.h"

using namespace std;

// constructor for the domain
// ! definition for this function is taken directly from the Fortran Code
Domain::Domain(double start, double end, ll num_points){
    // after class variable initialization we need to initialize the domain vector

    x_start = start;
    x_end = end;
    num_grid_points = num_points;

    double dx = (x_end-x_start)/num_grid_points;

    for(int i=1 ; i<=num_grid_points ; i++){
        domain.push_back( x_start + (i-0.5)*(dx) );
    }

}

// utility function to show the domain points after discretizing
void Domain::show_domain(){
    cout << "------------------Discretized domain---------------" << endl;
    for(int i=0 ; i<domain.size() ; i++){
        cout << "x_" << i+1 << " : " << domain[i] << endl;
    }
}

// function to pass the domain vector where required
vector<double> Domain::get_domain(){
    return domain;
}


