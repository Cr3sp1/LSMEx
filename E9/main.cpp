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
#include "lez9.h"


using namespace std;
using namespace arma;


int main (int argc, char *argv[]){

    Random rand( "Primes", "seed.in");

    int n_cities, n_paths, n_gen;

    // Reads the input file
    ifstream ReadInput;
    ReadInput.open("input.in"); 
    ReadInput >> n_cities;
    ReadInput >> n_paths;
    ReadInput >> n_gen;

    vec L1best_C(n_gen), L1best_S(n_gen), L1avg_C(n_gen), L1avg_S(n_gen);
    
    // Builds the matrices containing the positions of the cities
    mat CitiesC = CityCircle(rand);
    mat CitiesS = CitySquare(rand);
    for( int i = 1; i < n_cities; i++ ){
        CitiesC = arma::join_vert(CitiesC, CityCircle(rand));
        CitiesS = arma::join_vert(CitiesS, CitySquare(rand));
    }

    // Builds the distance matrices
    mat DistC = Distances(CitiesC);
    mat DistS = Distances(CitiesS);

    // Builds and sorts the populations, columns are the paths, rows are the cities (last city is not included)
    mat PopC = Population( n_cities, n_paths, rand );
    LSort(PopC, DistC);
    mat PopS = Population( n_cities, n_paths, rand );
    LSort(PopS, DistS);
    // Fills vectors
    L1best_C(0) = lenght(PopC.col(0), DistC);
    L1best_S(0) = lenght(PopS.col(0), DistS);
    L1avg_C(0) = Avlen( PopC, DistC);
    L1avg_S(0) = Avlen( PopS, DistS);



    // Executes the steps
    for( int i = 0; i < n_gen; i++ ) {
        // Evolves population
        NGE(PopC, DistC, rand);
        NGE(PopS, DistS, rand);

        // Fills vectors
        L1best_C(i) = lenght(PopC.col(0), DistC);
        L1best_S(i) = lenght(PopS.col(0), DistS);
        L1avg_C(i) = Avlen( PopC, DistC);
        L1avg_S(i) = Avlen( PopS, DistS);
    }

    // Prints stuff to file
    L1best_C.save("L1bestC.dat", arma::raw_ascii);
    L1best_S.save("L1bestS.dat", arma::raw_ascii);
    L1avg_C.save("L1avgC.dat", arma::raw_ascii);
    L1avg_S.save("L1avgS.dat", arma::raw_ascii);
    BestPath( PopC, DistC, CitiesC ).save("PathC.dat", arma::raw_ascii);
    BestPath( PopS, DistS, CitiesS ).save("PathS.dat", arma::raw_ascii);
    PopC.save("PopC.dat",  arma::raw_ascii);
    PopS.save("PopS.dat",  arma::raw_ascii);
    CitiesC.save("CitiesC.dat",  arma::raw_ascii);
    CitiesS.save("CitiesS.dat",  arma::raw_ascii);

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
