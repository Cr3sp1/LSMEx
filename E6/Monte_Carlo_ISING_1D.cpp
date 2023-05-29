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
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadConf;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs
  type = metro ? "metro" : "gibbs";  // Added

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> reuse;

  ReadInput >> prime;
  cout << "Priming = " << prime << endl;


  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

  


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration 
  // for (int i=0; i<nspin; ++i)
  // {
  //   if(rnd.Rannyu() >= 0.5) s[i] = 1;
  //   else s[i] = -1;
  // }

  // Added
  if(reuse){
    // Reads last generated configuration
    ReadConf.open("config.final");
    for( int i = 0; i < nspin; i++ ){
      ReadConf >> s[i];
    }
    ReadConf.close();
  } else {
    // Initializes all spins as up
    for (int i=0; i<nspin; ++i) {
      s[i] = 1;
    }
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
//Print other observables   // Added
  cout << "Initial heat capacity = " << beta*beta*(walker[ic]-walker[iu]*walker[iu])/nspin << endl;   // Always returns zero
  cout << "Initial susceptibility = " << beta*walker[ix] << endl;
  cout << "Initial magnetization = " << walker[im] << endl;

// Added
Prime(prime);
}


void Move(int metro)
{
  int o;
  // double p, energy_old, energy_new, sm;
  // double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    attempted++;    // Added
    
    if(metro==1) //Metropolis
    {
      // Added
      if( rnd.Rannyu() < exp(-beta*EnDiff(o)) ) {
        s[o]= -s[o];
        accepted++;
      }
    }
    else //Gibbs sampling
    { 
      // Calculates the probability for it to be different from its current value, then either changes it or leaves it the same // Added
      if( rnd.Rannyu() > 1/( 1 + exp(-beta*EnDiff(o)) ) ) {   
        s[o] = -s[o];
      }
      accepted++;
    }
  }
}

double Boltzmann(int sm, int ip)  // sm is spin value, ip spin position
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

// Evalues E_f - E_i from just the position // Added
double EnDiff( int ip ) {
  return 2*(J*s[ip]*( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) + h*s[ip]);
}


void Measure()
{
  // int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];   // Added
  }
  walker[iu] = u;
  // Added
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=14;   // Changed 12 -> 14
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    Ene.open("output.ene."+type,ios::app);      // Changed "output.ene.0" -> "output.ene."+type
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    // Added
    Heat.open("output.heat."+type,ios::app);      
    stima_c = beta*beta*(blk_av[ic]-blk_av[iu]*blk_av[iu]/blk_norm) /(blk_norm*(double)nspin); //Heat Capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_u=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_u << endl;
    Heat.close();

    Mag.open("output.mag."+type,ios::app);    
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_u=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_u << endl;
    Mag.close();

    Chi.open("output.chi."+type,ios::app);
    stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_u=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_u << endl;
    Chi.close();

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf, Ene, Heat, Mag, Chi;  // Added

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  // Prints the final values to file    // Added
  const int wd=14;
  if( h == 0 ) {

    Ene.open("final.ene."+type,ios::app);
    double err = Error(glob_av[iu],glob_av2[iu],nblk);
    Ene << setw(wd) << temp <<  setw(wd)  << glob_av[iu]/(double)nblk << setw(wd) << err << endl;
    Ene.close();

    Heat.open("final.heat."+type,ios::app);
    err = Error(glob_av[ic],glob_av2[ic],nblk);
    Heat << setw(wd) << temp <<  setw(wd)  << glob_av[ic]/(double)nblk << setw(wd) << err << endl;
    Heat.close();

    Chi.open("final.chi."+type,ios::app);
    err = Error(glob_av[ix],glob_av2[ix],nblk);
    Chi << setw(wd) << temp <<  setw(wd)  << glob_av[ix]/(double)nblk << setw(wd) << err << endl;
    Chi.close();

  } else {

    Mag.open("final.mag."+type,ios::app);
    double err = Error(glob_av[im],glob_av2[im],nblk);
    Mag << setw(wd) << temp <<  setw(wd)  << glob_av[im]/(double)nblk << setw(wd) << err << endl;
    Mag.close();

  }

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void Prime(int n){
  for( int i = 0; i < n; i++){
    Move(metro);
  }
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
