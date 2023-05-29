#ifndef __INTEGRAL_H__
#define __INTEGRAL_H__

#include <iostream>
#include <algorithm>
#include "Functions.h"
#include "random.h"

using namespace std;


// 1-dimensional integral
class IntegralMC {

  public:

    IntegralMC(): m_RNG() {}
    IntegralMC( const char* PrimesFile, const char* SeedFile ): m_RNG( PrimesFile, SeedFile) {}
    ~IntegralMC() {}

    void SaveSeed(){
      m_RNG.SaveSeed();
    }

    //Uses the average method
    double IntegralAve( double xmin, double xmax, Function &f, unsigned int n) {
      double sum = 0;
      for( int i = 0; i < (int)n; i++ ) {
        sum += f.Eval( m_RNG.Rannyu( xmin, xmax ) );
      }
      return ( xmax - xmin )*sum/n;
    }

    //Uses the hit or miss method
    double IntegralHoM( double xmin, double xmax, double fmax, Function &f, unsigned int n ) {
      unsigned int nhit = 0;
      for( int i = 0; i < (int)n; i++ ) {
        if ( m_RNG.Rannyu( 0, fmax ) < f.Eval( m_RNG.Rannyu( xmin, xmax ) ) ){
          nhit += 1;
        }
      }
      return ( xmax - xmin )*fmax*nhit/n;
    }

    //Uses the average method with importance sampling
    double IntegralIS( double xmin, double xmax, Function &Integrand, Function &P_x, Function &Quantile, unsigned int n ) {
      double sum = 0;
      for( int i = 0; i < (int)n; i++ ) {
        double x = Quantile.Eval( m_RNG.Rannyu(xmin, xmax));
        sum += Integrand.Eval(x)/P_x.Eval(x);
      }
      return ( xmax - xmin )*sum/n;
    }



  protected:

    Random m_RNG;

};







// n-dimensional integral (work in progress)
class IntegralMCnD {

  public:

    IntegralMCnD(): m_RNG() {}
    IntegralMCnD(const char* PrimesFile, const char* SeedFile ): m_RNG( PrimesFile, SeedFile) {}
    ~IntegralMCnD() {}


    double IntegralAve( vector<double> &xmin, vector<double> &xmax, ScalarFunction &f, unsigned int n ) {
      double sum = 0;
      double vol = 1;
      for( int i = 0; i < (int)n; i++ ) {
        vector<double> x( xmin.size() );
        for( int j = 0; j < (int)xmin.size(); j++){
          x[j]= m_RNG.Rannyu( xmin[j], xmax[j] );
          if( i == 0 ) vol = vol*(xmax[j] - xmin[j]);
        }
        sum += f.Eval( x );
      }
      return vol*sum/n;
    }

    double IntegralHoM( vector<double> &xmin, vector<double> &xmax, double fmax, ScalarFunction &f, unsigned int n ) {
      unsigned int nhit = 0;
      double vol = 1;
      for( int i = 0; i < (int)n; i++ ) {
        vector<double> x( xmin.size() );
        for( int j = 0; j < (int)xmin.size(); j++){
          x[j]= m_RNG.Rannyu( xmin[j], xmax[j] );
          if( i == 0 ) vol = vol*(xmax[j] - xmin[j]);
        }
        if(m_RNG.Rannyu( 0, fmax ) < f.Eval( x ) ) nhit++;
      }
      return vol*fmax*nhit/n;
    }

  protected:

    Random m_RNG;
};




#endif 
