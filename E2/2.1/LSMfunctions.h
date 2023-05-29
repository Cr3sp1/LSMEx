#ifndef __LSMfunctions__
#define __LSMfunctions__

#include <armadillo>
#include <iostream>
#include <cmath>

using namespace std;


double block_error(const arma::vec& BlockAvg, int N) {
    if ( (int)BlockAvg.n_elem <= N) {
        cerr << "Error: N is too large." << endl;
        return 0;
    }

    if( N == 0 ) { return 0; }

    // Resizes the input vector to size N
    arma::vec Ai = BlockAvg;
    Ai.resize(N+1);
    // Calculates the squared block averages
    arma::vec squared = square(Ai);

    // Calculates and return the error
    double error = sqrt((mean(squared) - pow(mean(Ai), 2)) / N);
    return error;
}

















#endif