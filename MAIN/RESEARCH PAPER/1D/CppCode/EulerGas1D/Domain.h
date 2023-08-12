/*
    Class Domain

    -> using the given value constructs a grid of points on the x-axis
*/

#pragma once

#include<vector>

#define ll long long int

using namespace std;

class Domain{
private:
    // variables used to define the computational domain
    double x_start;
    double x_end;
    ll num_grid_points;
    vector<double> domain;

public:
    // class constructor definition
    Domain(double start, double end, ll num_points);

    // utility function to show the discretized domain vector to the user
    void show_domain();

    // function to pass the domain vector where required
    vector<double> get_domain();
};

