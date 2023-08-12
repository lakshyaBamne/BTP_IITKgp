#pragma once

#include<vector>

#include "ConservedVariables.h"
#include "Constants.h"

using namespace std;

/*
    Class to append ghost values to the computation vectors
    according to a given boundary condition
*/
class ExtendCell{
public:

    // constructor to Extend the cells of a given set of conserved variable vectors
    ExtendCell(ConservedVariables* conserved_variables);


};