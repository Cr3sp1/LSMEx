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

    int nstep = 100; // Number of steps
    int M = 10000; // Number of throws
    int N = 100;  // Number of blocks

    vec disc( nstep );   // Vector containing the average of all the blocks for each step
    mat disc_i( nstep, N );  // Matrix containing the values of each step for each block
    vec disc_err( nstep ); //Matrix containing the error of each step

    // Same but for continuus case
    vec cont( nstep );   
    mat cont_i( nstep, N );  
    vec cont_err( nstep ); 

    for( int i = 0; i < N; i++ ){       // Iterates through the blocks
        for( int j = 0; j < M/N; j++ ){     // Iterates inside the block
            for( int k = 0; k < nstep; k++ ){       //Iterates through the steps   
                vec disc_k(3);
                vec cont_k(3);
                for( int l = 0; l <= k; l++) {
                    // Generates and adds the discrete step
                    int disc_step = rand.Discrete(0,5); 
                    int sign = -1 + 2*(disc_step%2); // Even = negative step, Odd = positive step
                    disc_k( disc_step/2 ) += sign;

                    // Generates and adds the continuus step
                    vec v = rand.Sphere(3,1);
                    cont_k(0) += v(0);
                    cont_k(1) += v(1);
                    cont_k(2) += v(2);
                }
                // Builds the averages inside the blocks
                disc_i( k, i ) += ( pow( disc_k(0), 2) + pow( disc_k(1), 2) + pow( disc_k(2), 2) )/ ( M/N);
                cont_i( k, i ) += ( pow( cont_k(0), 2) + pow( cont_k(1), 2) + pow( cont_k(2), 2) )/ ( M/N);
            }
        }     
    }

    
    for( int i = 0; i < nstep; i++ ){
        drowvec dis_i_sqrt = arma::sqrt(disc_i.row(i));
        drowvec con_i_sqrt = arma::sqrt(cont_i.row(i));

        // Fills the averages
        disc(i) = arma::mean( dis_i_sqrt);
        cont(i) = arma::mean( con_i_sqrt);

        // Fills the errors
        disc_err(i) = block_error( dis_i_sqrt, N-1 );
        cont_err(i) = block_error( con_i_sqrt, N-1 );
    }
    // disc.print();
    // cont.print();


    // Prints everything to file 
    mat Disc = arma::join_horiz( disc, disc_err);
    Disc.save("disc.dat", arma::raw_ascii);

    mat Cont = arma::join_horiz( cont, cont_err );
    Cont.save("cont.dat", arma::raw_ascii);

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
