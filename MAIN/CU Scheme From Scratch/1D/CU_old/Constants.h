#pragma once

#include<iostream>

using namespace std;

namespace Constants{
    const double GAMMA = 1.4; // ratio of specific heats
    // const double GAMMA = 5/3;
    const double THETA = 1.3; // tuning parameter to control oscillations
    const double CFL = 0.45; // CFL number
    const double EPS = 1.0E-12; //  Epsilon value used in calculating Flux
}

