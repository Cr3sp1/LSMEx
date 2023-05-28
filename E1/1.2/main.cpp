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
   //declare vectors

   arma::vec St1(10000);
   arma::vec St2(10000);
   arma::vec St10(10000);
   arma::vec St100(10000);

   arma::vec Ex1(10000);
   arma::vec Ex2(10000);
   arma::vec Ex10(10000);
   arma::vec Ex100(10000);

   arma::vec Lo1(10000);
   arma::vec Lo2(10000);
   arma::vec Lo10(10000);
   arma::vec Lo100(10000);

   for( int i = 0; i < 10000; i++ ){

      //fill N = 1
      double s = rnd.Rannyu();
      double e = rnd.Exponential(1);
      double l = rnd.CauchyLorentz(0, 1);
      St1[i] = s;
      Ex1[i] = e;
      Lo1[i] = l;

      //fill N = 2
      s /= 2;
      e /= 2;
      l /= 2;
      s += rnd.Rannyu()/2;
      e += rnd.Exponential(1)/2;
      l += rnd.CauchyLorentz(0,1)/2;
      St2[i] = s;
      Ex2[i] = e;
      Lo2[i] = l;

      //fill N = 10
      s /= 5;
      e /= 5;
      l /= 5;
      for( int j = 0; j < 8; j++ ){
         s += rnd.Rannyu()/10;
         e += rnd.Exponential(1)/10;
         l += rnd.CauchyLorentz(0,1)/10;
      }
      St10[i] = s;
      Ex10[i] = e;
      Lo10[i] = l;

      //fill N = 100
      s /= 10;
      e /= 10;
      l /= 10;
      for( int j = 0; j < 90; j++ ){
         s += rnd.Rannyu()/100;
         e += rnd.Exponential(1)/100;
         l += rnd.CauchyLorentz(0,1)/100;
      }
      St100[i] = s;
      Ex100[i] = e;
      Lo100[i] = l;
   }

   //prints all the data
   St1.save( "St1.dat", arma::raw_ascii );
   St2.save( "St2.dat", arma::raw_ascii );
   St10.save( "St10.dat", arma::raw_ascii );
   St100.save( "St100.dat", arma::raw_ascii );

   Ex1.save( "Ex1.dat", arma::raw_ascii );
   Ex2.save( "Ex2.dat", arma::raw_ascii );
   Ex10.save( "Ex10.dat", arma::raw_ascii );
   Ex100.save( "Ex100.dat", arma::raw_ascii );

   Lo1.save( "Lo1.dat", arma::raw_ascii );
   Lo2.save( "Lo2.dat", arma::raw_ascii );
   Lo10.save( "Lo10.dat", arma::raw_ascii );
   Lo100.save( "Lo100.dat", arma::raw_ascii );



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
