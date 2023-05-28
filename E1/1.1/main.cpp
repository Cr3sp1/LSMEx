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
#include<cmath>
#include "random.h"
#include "LSMfunctions.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   rnd.SetRandFromFile( "Primes", "seed.in" );

   int M = 1000000;
   int N = 100;
   int L = M/N;

   int M2 = 100;
   int N2 = 10000;

   arma::vec AvgN(N); // media dei primi i+1 blocchi
   arma::vec Ai(N); //media del blocco i
   arma::vec AvgEr(N);

   arma::vec Sig2(N);
   arma::vec Sig2i(N);
   arma::vec Sig2Er(N);

   arma::vec Chi2(M2);
   arma::vec Ni(M2);


   for( int i=0; i<N; i++ ) {

      //runs the block
      for( int j=0; j<L; j++ ) {
         double r = rnd.Rannyu();
         Ai[i] += r/L;
         Sig2i[i] += pow( r - 0.5, 2 )/L;
      }

      //assigns the values
      if( i == 0 ) {
         AvgN[i] = Ai[i];
         Sig2[i] = Sig2i[i];
      } else {
         AvgN[i] = AvgN[i-1]*i/(i+1) + Ai[i]/(i+1);
         Sig2[i] = Sig2[i-1]*i/(i+1) + Sig2i[i]/(i+1);
      }
   }

   //builds the errors
   for( int i = 0; i < N; i++ ) {
      AvgEr[i] = calculate_error( Ai, i );
      Sig2Er[i] = calculate_error( Sig2i, i );
   }

   
   //part 3
   for( int i = 0; i < M2; i++ ) {
      // resets Ni
      Ni.zeros();

      // fills the sub intervals
      for( int j = 0; j < N2; j++ ) {
         //adds 1 to the correct Ni[k];
         int k = rnd.Rannyu()*100;
         Ni[k]++;
      }

      //evaluates Chi2[i]
      for( int j=0; j < M2; j++ ) {
         Chi2[i] += pow( Ni[j] - N2/M2, 2 )/(N2/M2);
      }

   }

   //print everything to file 
   arma::mat Avg = arma::join_horiz( AvgN, AvgEr );
   Avg.save("Avg.dat", arma::raw_ascii);

   arma::mat Sig = arma::join_horiz( Sig2, Sig2Er );
   Sig.save("Sig2.dat", arma::raw_ascii);

   Chi2.save("Chi2.dat", arma::raw_ascii);


   rnd.SaveSeed();
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
