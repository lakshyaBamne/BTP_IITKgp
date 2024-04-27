#pragma once

#include "Constants.h"

using namespace std;

namespace CTS = Constants;

namespace SecondOrderReconstructions{

    // minmod limiter
    double minmod(double a, double b, double c);

    // function to calculate the reconstruction slopes
    vvd get_slopes(const vvd& U);

    // function to calculate the 2nd order reconstructions
    vvvd get_reconstructions(const vvd& U, const vvd& slx);

    // function to extend the 2nd order reconstructions
    void extend_reconstructions(vvd& Up, vvd& Um);

}

double SecondOrderReconstructions::minmod(double a, double b, double c){
    double max_val = max(a, max(b,c));
    double min_val = min(a, min(b,c));

    if(max_val < 0){
        return max_val;
    }
    else if(min_val > 0){
        return min_val;
    }
    else{
        return 0;
    }
}

vvd SecondOrderReconstructions::get_slopes(const vvd& U){
    vvd slx(3, vd(CTS::N+2));

    for(int i=0 ; i<3 ; i++){
        for(int j=1 ; j<=CTS::N ; j++){
            slx[i][j] = minmod(
                CTS::THETA*(U[i][j] - U[i][j-1]),
                0.5*(U[i][j+1] - U[i][j-1]),
                CTS::THETA*(U[i][j+1] - U[i][j])
            );
        }
    }

    return slx;
}

vvvd SecondOrderReconstructions::get_reconstructions(const vvd& U, const vvd& slx){
    vvd Up(3, vd(CTS::N+1));
    vvd Um(3, vd(CTS::N+1));

    for(int i=0 ; i<3 ; i++){
        // U+
        for(int j=1 ; j<=CTS::N ; j++){
            Up[i][j] = U[i][j] + 0.5*slx[i][j];
        }

        // U-
        for(int j=0 ; j<=CTS::N-1 ; j++){
            Um[i][j] = U[i][j+1] - 0.5*slx[i][j+1];
        }
    }

    // Extend the reconstructions
    extend_reconstructions(Up, Um);

    return {Up, Um};
}

void SecondOrderReconstructions::extend_reconstructions(vvd& Up, vvd& Um){
    if(CTS::BC == "FREE"){
        for(int i=0 ; i<3 ; i++){
            Up[i][0] = Um[i][0];
            Um[i][CTS::N] = Up[i][CTS::N];
        }
    }
    else{
        cout << "---ERROR--- Please select correct boundary conditions---" << endl;
    }
}