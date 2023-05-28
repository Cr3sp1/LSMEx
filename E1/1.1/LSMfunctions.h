#ifndef __LSMfunctions__
#define __LSMfunctions__

#include <armadillo>
#include <iostream>
#include <cmath>

using namespace std;


double calculate_error(const arma::vec& values, int N) {
    if ( (int)values.n_elem <= N) {
        cerr << "Error: N is too large." << endl;
        return 0;
    }

    if( N == 0 ) { return 0; }

    // Resize the input vector to size N
    arma::vec Ai = values;
    Ai.resize(N+1);
    // Calculate the squared values and the mean
    arma::vec squared = square(Ai);

    // Calculate and return the error
    double error = sqrt((mean(squared) - pow(mean(Ai), 2)) / N);
    return error;
}



#endif