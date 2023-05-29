#ifndef __lez8__
#define __lez8__

#include <armadillo>
#include <iostream>
#include <cmath>
#include "random.h"

using namespace std;
using namespace arma;


// Calculates psi^2
double psi2(double x, double mu, double sigma) {
    double psi = exp(-pow(x-mu, 2)/(2*pow(sigma, 2))) + exp(-pow(x+mu, 2)/(2*pow(sigma, 2)));
    return psi*psi;
}

// Calculates Eloc
double Eloc(double x, double mu, double sigma) {
    double term1 = exp(-pow(x-mu, 2)/( 2*pow(sigma, 2) ));
    double term2 = exp(-pow(x+mu, 2)/( 2*pow(sigma, 2) ));
    double numerator1 = pow(x-mu, 2) - pow(sigma, 2);
    double numerator2 = pow(x+mu, 2) - pow(sigma, 2);
    double psi = term1 + term2;
    double kin = -0.5*(term1 * numerator1+ term2 * numerator2)/(pow(sigma,4)*psi);
    double pot = pow(x,4) - pow(x,2)*5/2;
    return kin + pot;
}

// Moves x and returns 1 if the step was accepted and 0 otherwise
int Move( double mu, double sigma, double stepsize, double &x, Random& rand ) {
    double xnew = x + rand.Rannyu(-stepsize, stepsize);
    if( rand.Rannyu() < psi2(xnew, mu, sigma)/psi2(x, mu, sigma)) {
        x = xnew;
        return 1;
    }
    return 0;
}


// Returns the value of H, its error and the acceptace rate
vec3 H_Err_Acc( double mu, double sigma, int nblocks, int blocksize, double stepsize, Random& rand ) {
    double x = 0;
    int tried = 0;
    int accepted = 0;
    vec3 result;
    vec BlockAv(nblocks);

    for( int i = 0; i < nblocks; i++ ) {       
        for(int j = 0; j < blocksize; j++ ) {
            tried++;
            accepted += Move( mu, sigma, stepsize, x, rand);
            BlockAv[i] += Eloc(x, mu, sigma)/blocksize;
        }
    }

    result[0] = arma::mean(BlockAv);
    result[1] = arma::stddev(BlockAv)/sqrt(nblocks);
    result[2] = (double)accepted/(double)tried;

    return result;
}


// Returns the averages of the blocks of the value of H
vec Hblocks( double mu, double sigma, int nblocks, int blocksize, double stepsize, Random& rand ) {
    double x = 0;
    int tried = 0;
    int accepted = 0;
    vec BlockAv(nblocks);

    for( int i = 0; i < nblocks; i++ ) {       
        for(int j = 0; j < blocksize; j++ ) {
            tried++;
            accepted += Move( mu, sigma, stepsize, x, rand);
            // double xnew = x + rand.Rannyu(-stepsize, stepsize);
            // if( rand.Rannyu() < psi2(xnew, mu, sigma)/psi2(x, mu, sigma)) {
            //     x = xnew;
            //     accepted++;
            // }
            BlockAv[i] += Eloc(x, mu, sigma)/blocksize;
        }
        cout << "Block = " << i << ", x = " << x << endl;
    }

    // cout << "Mu = " << mu << ", sigma = " << sigma << ", <H> = " << arma::mean(BlockAv) << " +- " << arma::stddev(BlockAv)/sqrt(nblocks) << ", Acceptance rate = " << (double)accepted/(double)tried << endl;
    return BlockAv;
}


double Hval( double mu, double sigma, int nstep, double stepsize, Random& rand ) {
    double x = 0;
    double result = 0;

    for( int i = 0; i < nstep; i++ ) {
        // double xnew = x + rand.Rannyu(-stepsize, stepsize);
        // if( rand.Rannyu() < psi2(xnew, mu, sigma)/psi2(x, mu, sigma)) {
        //     x = xnew;
        // }
        Move(mu, sigma, stepsize, x, rand);
        result += Eloc(x, mu, sigma);
    }

    return result/nstep;
}

// Move function for the simulated annealing, returns 1 if it succeeds, 0 othewise
int SA_Move( double& mu, double& sigma, int mcsteps, double temp, Random& rand ){
    double mu_step = 0.1;
    double sigma_step = 0.1;
    double mu_new = mu;
    double sigma_new = sigma;
    if( rand.Rannyu() > 0.5) {
        mu_new += rand.Rannyu(-mu_step, mu_step);
    } else {
        sigma_new += rand.Rannyu(-sigma_step, sigma_step);
    }

    double prob = exp(-( Hval(mu_new, sigma_new, mcsteps, 4*sigma_new, rand) - Hval(mu, sigma, mcsteps, 4*sigma_new, rand) )/temp);
    if( rand.Rannyu() < prob ){
        mu = mu_new;
        sigma = sigma_new;
        return 1;
    } else {
        return 0;
    }
}

#endif