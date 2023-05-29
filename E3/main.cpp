/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <cmath>
#include "random.h"
#include "LSMfunctions.h"
#include "Integral.h"
#include "Functions.h"


using namespace std;
using namespace arma;


int main (int argc, char *argv[]){

    Random rand( "Primes", "seed.in");

    int M = 100000;         // Number of throws
    int N = 100;            // Number of blocks
    int L = M/N;            // Number of throws per block

    double S_0 = 100;         // Asset price at t=0
    double T =1 ;           // Delivery time
    double k = 100;         // Stike price
    double r = 0.1;         // Risk free intrest rate
    double sigma = 0.25;    // Volatility
    int n_step = 100;       // Number of steps in the discretized path

    // Call valued directly
    vec CallDir(N);

    // Put valued directly
    vec PutDir(N); 

    // Call valued throuh discrete path
    vec CallDisc(N);  

    // Put valued throuh discrete path 
    vec PutDisc(N); 


    for( int i = 0; i < N; i++ ){
        for( int j = 0; j < L; j++ ){
            // Evaluate S(T) directly
            double S_T = S_0*exp( (r - pow(sigma, 2)/2)*T + sigma*rand.Gauss(0, T) );
            
            // Evaluate S(T) through discretized path
            double S_i = S_0;
            for( int l = 0; l < n_step; l++ ){
                S_i *= exp( (r - pow(sigma, 2)/2)*T/n_step + sigma*rand.Gauss(0, 1)*sqrt(T/n_step) );
            }

            // Evaluate and add Put and Call prices to the averages
            CallDir(i) += exp(-r*T)*max(0., (S_T-k));
            PutDir(i) += exp(-r*T)*max(0., (k-S_T));
            CallDisc(i) += exp(-r*T)*max(0., (S_i-k));
            PutDisc(i) += exp(-r*T)*max(0., (k-S_i));
        }
        // Normalizes the averages
        CallDir(i) /= L;
        PutDir(i) /= L;
        CallDisc(i) /= L;
        PutDisc(i) /= L;
    }

    print_avg_err( CallDir, "CallDir.dat");
    print_avg_err( PutDir, "PutDir.dat");
    print_avg_err( CallDisc, "CallDisc.dat");
    print_avg_err( PutDisc, "PutDisc.dat");

    return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
