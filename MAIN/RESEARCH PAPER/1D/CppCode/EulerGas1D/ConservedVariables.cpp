#include<iostream>
#include<vector>
#include<string>
#include<unordered_map>

#define ll long long int

#include "ConservedVariables.h"
#include "Constants.h"

using namespace std;

/*
    Class to initialize the computational domain and the values for 
    the conserved variables according to the given Riemann problem in 1D
*/
// function to show the underlying vectors containing the conserved variables
void ConservedVariables::show_conserved_variables(){
    cout << "--------------Conserved Variables-----------------" << endl;

    for(int i=0 ; i<conserved_variables[0].size() ; i++){
        cout << "Density : " << conserved_variables[0][i] << " | Momentum : " << conserved_variables[1][i] << " | Energy : " << conserved_variables[2][i] << " |" << endl; 
    }
}

// constructor to the class takes in as argument the discretized computational domain
// then initializes the Conserved variables using the given 
ConservedVariables::ConservedVariables(Domain* domain){
    /*
        To change the problem change the conditions here in the constructor
        ! Later add a more user friendly way to define initial conditions.
    */

    // get the discretized domain points 
    auto x = domain->get_domain();

    // CASE-1 : Moving Contact Wave with Free boundary conditions
    for(auto i : x){
        if( i<0.3 ){
            density.push_back(1.4);
        }
        else if( i>0.3 ){
            density.push_back(1.0);
        }

        momentum.push_back( 0.1 * density.back() );
        energy.push_back( (1/(Constants::GAMMA-1)) + (density.back()*0.1*0.1)/2 );

    }
    // now encapsulate the three vectors into one for easier usage
    conserved_variables.push_back(density);
    conserved_variables.push_back(momentum);
    conserved_variables.push_back(energy);

}

