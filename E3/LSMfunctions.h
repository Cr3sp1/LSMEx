#ifndef __LSMfunctions__
#define __LSMfunctions__

#include <armadillo>
#include <iostream>
#include <cmath>

using namespace std;
using namespace arma;


double block_error(const vec& BlockAvg, int N) {
    if ( (int)BlockAvg.n_elem <= N) {
        cerr << "Error: N is too large." << endl;
        return 0;
    }

    if( N == 0 ) return 0;

    // Resizes the input vector to size N+1
    vec Ai = BlockAvg;
    Ai.resize(N+1);
    // Calculates the squared block averages
    vec squared = arma::square(Ai);

    // Calculates and return the error
    double error = sqrt((arma::mean(squared) - pow(arma::mean(Ai), 2)) / N);
    if( error != error ) return 0;
    return error;
}

// Same but for rows
double block_error(const rowvec& BlockAvg, int N) {
    if ( (int)BlockAvg.n_elem <= N) {
        cerr << "Error: N is too large." << endl;
        return 0;
    }

    if( N == 0 ) { return 0; }

    // Resizes the input vector to size N
    vec Ai = BlockAvg.t();
    Ai.resize(N+1);
    // Calculates the squared block averages
    vec squared = arma::square(Ai);

    // Calculates and return the error
    double error = sqrt((arma::mean(squared) - pow(arma::mean(Ai), 2)) / N);
    if( error != error ) {return 0;}
    return error;
}


// Builds the error vector
vec block_error(const vec& BlockAvg) {
    int size = BlockAvg.n_elem;
    vec err(size);
    for( int i = 1; i < size; i++ ) { // For cicle starts from 1 so that err(0) remains 0

        // Resizes the input vector to size i+1
        vec Ai = BlockAvg;
        Ai.resize(i+1);
        // Calculates the squared block averages
        vec squared = arma::square(Ai);

        // Calculates and return the error
        double error = sqrt((arma::mean(squared) - pow(arma::mean(Ai), 2)) / i);
        if( error != error ) {  // Only happens if error = nan
            err(i) = 0;
        } else{
            err(i) = error;
        }
    }

    return err;
}

// Builds the average of first i blocks vector
vec tot_avg(const vec& BlockAvg) {

    int size = BlockAvg.n_elem;
    vec tot(size);
    tot(0) = BlockAvg(0);

    for( int i = 1; i < size; i++ ) { 
        tot(i) = ( tot(i-1)*i + BlockAvg(i))/(i+1);
    }
    
    return tot;
}

// Builds and prints in the same file the running average and its error as vectors
void print_avg_err(const vec& BlockAvg, const char* namefile){
    mat mute = arma::join_horiz( tot_avg(BlockAvg), block_error(BlockAvg));
    mute.save(namefile, arma::raw_ascii);
}










#endif