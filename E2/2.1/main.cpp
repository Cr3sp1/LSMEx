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


int main (int argc, char *argv[]){

    IntegralMC Integral( "Primes", "seed.in" );
    Es2I Integrand;
    Es2P Prob;
    Es2Q Quantile;

    int M = 100000;
    int N = 100;
    int L = M/N;


    //Declares the needed vectors
    arma::vec I1Avg(N);  //The average of the first i blocks
    arma::vec I1Ai(N);   //The result of block i 
    arma::vec I1Err(N);  //Error valued on the first i blocks

    //Same but for importance sampling
    arma::vec I2Avg(N);
    arma::vec I2Ai(N);
    arma::vec I2Err(N);


    for( int i = 0; i < N; i++ ) {
        //FIlls the blocks
        I1Ai[i] = Integral.IntegralAve( 0, 1, Integrand, L );
        I2Ai[i] = Integral.IntegralIS( 0, 1, Integrand, Prob, Quantile, L );

        //Fills the Averages
        if( i == 0 ){
            I1Avg[0] = I1Ai[0];
            I2Avg[0] = I2Ai[0];
        } else{
            I1Avg[i] = I1Avg[i-1]*i/(i+1) + I1Ai[i]/(i+1);
            I2Avg[i] = I2Avg[i-1]*i/(i+1) + I2Ai[i]/(i+1);
        }

        //Builds the errors
        I1Err[i] = block_error( I1Ai, i );
        I2Err[i] = block_error( I2Ai, i );

    }

    
    //print everything to file 
    arma::mat Int1 = arma::join_horiz( I1Avg, I1Err );
    Int1.save("Int1.dat", arma::raw_ascii);

    arma::mat Int2 = arma::join_horiz( I2Avg, I2Err );
    Int2.save("Int2.dat", arma::raw_ascii);


    Integral.SaveSeed();
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
