#include<iostream>
#include<vector>

// Class imports
#include "Domain.h"
#include "ConservedVariables.h"

using namespace std;

int main(){
    /*
        STEP-1

        -> Initialize the domain according to the given Riemann problem
        -> 1) Discretize the computational domain into finite volumes (equal cells)
        -> 2) Initialize the values for the conserved variables on this discretized domain
    */

    // initialize the domain points after discretizing
    // this initialization is done on the heap memory
    Domain* domain = new Domain(0, 1, 200);
    domain->show_domain();

    // after initializing the discretized domain, we can initialize the conserved variables
    ConservedVariables* conserved_variables = new ConservedVariables(domain);
    conserved_variables->show_conserved_variables();

    // now we need to extend the cells to accomodate ghost values
    // which are required in the Numerical Scheme to calculate next iteration values
    

    return 0;
}

