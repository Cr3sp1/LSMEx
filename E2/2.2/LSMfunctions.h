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
    double error = sqrt((arma::mean(squared) - pow(arma::mean(Ai), 2)) / N);
    if( error != error ) {return 0;}
    return error;
}

// Same but for rows
double block_error(const arma::rowvec& BlockAvg, int N) {
    if ( (int)BlockAvg.n_elem <= N) {
        cerr << "Error: N is too large." << endl;
        return 0;
    }

    if( N == 0 ) { return 0; }

    // Resizes the input vector to size N
    arma::vec Ai = BlockAvg.t();
    Ai.resize(N+1);
    // Calculates the squared block averages
    arma::vec squared = square(Ai);

    // Calculates and return the error
    double error = sqrt((arma::mean(squared) - pow(arma::mean(Ai), 2)) / N);
    if( error != error ) {return 0;}
    return error;
}















#endif