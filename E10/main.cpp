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
#include "mpi.h"
#include "random.h"
#include "LSMfunctions.h"
#include "lez10.h"



using namespace std;
using namespace arma;


int main (int argc, char *argv[]){

    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Status stat;

    // Random different for every node
    Random rand( "Primes", "seed.in", rank );
    // "Shared" random that gives the same result in every node
    Random shared( "Primes", "seed.in" );

    bool connected;
    int n_paths, n_gen, n_migrate;
    int n_cities = 50;

    // Reads the input file
    ifstream ReadInput;
    ReadInput.open("input.in"); 
    ReadInput >> connected;
    ReadInput >> n_paths;
    ReadInput >> n_gen;
    ReadInput >> n_migrate;

    // Declares vectors that we will print
    vec Best(n_gen/n_migrate + 1);
    vec BestHalf(n_gen/n_migrate + 1);
    vec Index(n_gen/n_migrate + 1);

    
    // Builds the matrix containing the positions of the cities
    mat Caps = Capitals(50);

    // Builds the distance matrix
    mat Dist = Distances(Caps);
    

    // Builds and sorts the population, columns are the paths, rows are the cities (last city is not included)
    mat Pop = Population( n_cities, n_paths, rand );
    LSort(Pop, Dist);
    mat TotPop = TotalPop( Pop, 0, rank, size);
    LSort(TotPop, Dist);
    // Fills vectors
    if( rank == 0 ){
        Best(0) = lenght(TotPop.col(0), Dist);
        BestHalf(0) = Avlen( TotPop, Dist);
    }



    // Executes the steps
    for( int i = 1; i <= n_gen; i++ ) {
        // Evolves population
        NGE(Pop, Dist, rand);

        if( i%n_migrate == 0 ){
            int ind = i/n_migrate;
            MPI_Barrier(MPI_COMM_WORLD);
            TotPop = TotalPop( Pop, 0, rank, size);

            // Fills vectors
            if( rank == 0 ) {
                // cout << "Migration number " << ind << endl;
                Best(ind) = lenght(TotPop.col(0), Dist);
                BestHalf(ind) = Avlen( TotPop, Dist);
                Index(ind) = i;
            }

            // Mixes them around
            if( connected == 1 ) {
                Switcherout( Pop, Dist, rank, size, shared );
            }
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    TotPop = TotalPop( Pop, 0, rank, size);
    // Prints stuff to file
    if( rank == 0 ) {
        // cout << "Printing" << endl;
        string type = connected ? "con." : "unc.";
        mat mute = arma::join_horiz(Index, Best);
        mute.save( type + "best.dat", arma::raw_ascii);
        mute = arma::join_horiz(Index, BestHalf);
        mute.save( type + "bestHalf.dat", arma::raw_ascii);
        BestPath( TotPop, Dist, Caps ).save( type + "path.dat", arma::raw_ascii);
        // Index.save("Pippo", arma::raw_ascii);
    }

    // if(rank == 0) {
    //     string type = connected ? "con." : "unc.";
    //     mat mute = BestPath( TotPop, Dist, Caps );
    //     mute.save( type + "path.dat", arma::raw_ascii);
    // }
    

    // rand.SaveSeed(rank);
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
