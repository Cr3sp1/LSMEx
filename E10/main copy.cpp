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

// Sends a vector using MPI
void SendVec( vec v, int toRank, int tag, MPI_Comm comm ){
    int vsize = v.size();
    MPI_Send(v.memptr(), vsize, MPI_DOUBLE, toRank, tag, comm);
}

// Receives a vector using MPI
vec ReceiveVec( int vsize, int fromRank, int tag, MPI_Comm comm ){
    vec v(vsize);
    MPI_Recv(v.memptr(), vsize, MPI_DOUBLE, fromRank, tag, comm, MPI_STATUS_IGNORE);
    
    return v;
}


int main (int argc, char *argv[]){

    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;

    vec sent = { rank/1., rank/1., rank/1., rank/1., rank/1., rank/1. };
    vec recv;
    if(rank == 0) SendVec( sent, 0, 0, MPI_COMM_WORLD );
    for(int i = 1; i < size-1; i++ ){
        if( i == rank) {
            recv = ReceiveVec(sent.size(), i-1, i-1, MPI_COMM_WORLD );
            SendVec( sent, i+1, i, MPI_COMM_WORLD);
        }
    }
    if (rank == size-1) {
        recv = ReceiveVec(sent.size(), size-2, size-2, MPI_COMM_WORLD );
        SendVec( sent, 0, size-1, MPI_COMM_WORLD);
    }
    if(rank == 0) recv = ReceiveVec(sent.size(), size-1, size-1, MPI_COMM_WORLD );

    cout << "Nodo =" << rank << "/" << size << ", v = " << recv(5) << endl;

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
