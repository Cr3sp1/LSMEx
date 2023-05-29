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
#include "lez8.h"


using namespace std;
using namespace arma;


int main (int argc, char *argv[]){

    Random rand( "Primes", "seed.in");

    double mu, stepsize, sigma, initialT, temp, alfa;
    int nblocks, blocksize, mcsteps, sasteps, midsteps;
    vec3 app;
    vec psi2(100);  //plots psi2(x), with x in (-3, 3)

    // Reads the input file
    ifstream ReadInput;
    ReadInput.open("input.in"); 

    ReadInput >> nblocks;
    vec Hblocks(nblocks);
    ReadInput >> blocksize;
    mcsteps = nblocks*blocksize;
    ReadInput >> mu;
    ReadInput >> sigma;
    stepsize = 4*sigma;
    ReadInput >> sasteps;
    vec Hsa(sasteps);
    vec Hsa_err(sasteps);
    ReadInput >> midsteps;
    ReadInput >> initialT;
    temp = initialT;
    ReadInput >> alfa;

    ReadInput.close();

    cout << "SA steps: " << sasteps << "; steps at each temperature: " << midsteps << "; MC steps: " << mcsteps << endl;
    cout << "Initial Values: T = " << temp << ", <H> = " << Hval( mu, sigma, mcsteps, stepsize, rand ) << ", mu = " << mu << ", sigma = " << sigma << endl;


    // SA algorithm
    for( int i = 0; i < sasteps; i++ ) {
        for( int j = 0; j < midsteps; j++) {
            SA_Move( mu, sigma, mcsteps, temp, rand);
        }
        app = H_Err_Acc( mu, sigma, nblocks, blocksize, stepsize, rand );
        Hsa[i] = app[0];
        Hsa_err[i] = app[1];
        temp *= alfa;
    }
    
    // Calculates <H> and psi2 for the best mu and sigma
    double x = 0;
    int ix;
    for( int i = 0; i < nblocks; i++ ) {
        for( int j = 0; j < blocksize; j++ ){
            Move( mu, sigma, 4*sigma, x, rand );
            if( -3 < x && x < 3 ) {
                ix = 100*(x+3)/6;
                psi2[ix]++;
            }
            Hblocks[i] += Eloc(x, mu, sigma)/blocksize;
        }
    }
    // Normalizes psi2
    psi2 /= (6*arma::mean(psi2));
    // cout << "1 = " << 6*arma::mean(psi2) << endl;

    //Prints everything
    arma::mat mute = arma::join_horiz( Hsa, Hsa_err);
    mute.save("H_T.dat", arma::raw_ascii);

    print_avg_err( Hblocks, "Hbest.dat");

    vec xpsi(100);
    for( int i = 0; i < 100; i++ ){
        if( i == 0) {
            xpsi[0] = -3;
        } else {
            xpsi[i] = xpsi[i-1] + 6./100;
        }
    }

    arma::mat mute2 = arma::join_horiz( xpsi, psi2 );
    mute2.save("Psi2.dat", arma::raw_ascii);

    vec ms = { mu, sigma, temp };
    ms.save("MuSigmaT.dat", arma::raw_ascii);

    cout << "Final Values: T = " << temp << ", <H> = " << Hsa[sasteps-1] << ", mu = " << mu << ", sigma = " << sigma << endl;

    rand.SaveSeed();
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
