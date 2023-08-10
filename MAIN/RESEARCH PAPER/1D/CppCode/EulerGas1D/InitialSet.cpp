#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>

#define ll long long int

#define GAMMA 1.4

#include "Domain.cpp"

using namespace std;

/*
    Class to initialize the computational domain and the values for 
    the conserved variables according to the given Riemann problem in 1D
*/
class ConservedVariables{
private:
    // vectors representing conserved variables initially according to the riemann problem
    vector<double> density, momentum, energy;

public:

    // function to show the underlying vectors containing the conserved variables
    void show_conserved_variables(){
        cout << "--------------Conserved Variables-----------------" << endl;

        for(int i=0 ; i<density.size() ; i++){
            cout << "Density : " << density[i] << " | Momentum : " << momentum[i] << " | energy : " << energy[i] << " |" << endl; 
        }
    }

    // constructor to the class takes in as argument the discretized computational domain
    // then initializes the Conserved variables using the given 
    ConservedVariables(Domain* domain){
        /*
            To change the problem change the conditions here in the constructor
            ! Later add a more user friendly way to define initial conditions.
        */

        // get the discretized domain points 
        auto x = domain->get_domain();

        for(auto i : x){
            if( i<0.3 ){
                density.push_back(1.4);
            }
            else if( i>0.3 ){
                density.push_back(1.0);
            }

            momentum.push_back( 0.1 * density.back() );
            energy.push_back( (1/(GAMMA-1)) + (density.back()*0.1*0.1)/2 );
        }

    }


};

