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

using namespace std;


//costruiamo un generatore di cosTheta che non dipenda in alcun modo da PI
double cosTheta( Random rand ){
// genrate a random point in the square (-0.5, 0.5]X(-0.5, 0.5]
double x = rand.Rannyu() - 0.5;
double y = rand.Rannyu() - 0.5;
//troviamo raggio rispetto a origine.
double r = sqrt(pow(x, 2) + pow(y, 2));
//se (x,y) non appartiene alla circonferenza di raggio 1/2 centrata nell'origine riptova.
if( r  >= 0.5 || r ==0 ) { return cosTheta( rand ); };
//ora abbiamo un punto generato uniformemente all'interno della circonferenza di raggio 0.5, anche il suo angolo theta rispetto all'asse delle x è generato uniformemente in (0,2Pi].
// per trovare cosTheta basta invertire x = r cosTheta; 
return x/r;
} 


 
int main (int argc, char *argv[]){


   Random rnd;
   rnd.SetRandFromFile( "Primes", "seed.in" );

   double d = 1; //distanza tra le linee
   double l = 0.7; //lunghezza ago
   double center; //posizione centro dell'ago
   int M = 10000000; //numero di lanci
   int Nh = 0; // numero di successi
   int N = 100; // numero di blocchi
   int L = M/N; //lanci per blocco

   arma::vec Pi(N);
   arma::vec PiI(N);
   arma::vec SigPi(N);


   for( int i = 0; i < N; i++) {

      Nh = 0; // resets Nh

      for( int j = 0; j < L; j++ ) {
         //consideriamo che il centro dell'ago cade in un punto tra due linee con distribuzione piatta.
         center = rnd.Rannyu();
         //lo spazio occupato verticalmente è uguale ad |L*cosTheta|.
         double h = abs( l * cosTheta(rnd) );
         if ( center + h/2 >= d or center - h/2 <= 0 ) {
            Nh++;
         }
      }
      //assegna PiI e Pi
      PiI[i] = 2*l*L/(Nh*d);

      if( i == 0 ){
         Pi[0] = PiI[0];
      } else{
         Pi[i] = ( Pi[i-1]*i + PiI[i] )/(i+1);
      }
   }


   //calculates error
   for( int i = 0; i < N; i++ ) {
      SigPi[i] = calculate_error( PiI, i );
   }

   arma::mat PI = arma::join_horiz( Pi, SigPi);
   PI.save("Pi.dat", arma::raw_ascii);


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
