#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

void printdata(std::vector<double>data,int N, int Nlocal);


int main(int argc, char **argv)
{
  const int N = std::atoi(argv[1]);
  int np, pid;
  int Nlocal=N/np;
  int tag =0;

    std::vector<double>data(N*Nlocal);
    std::fill(data.begin(),data.end(),0.0);

    for(int ilocal=0; ilocal<Nlocal;++ilocal){

      for(int jlocal = Nlocal*pid < Nlocal*(pid+1);++jlocal;) {
        data[ilocal*(N+1)+jlocal]=1.0;
        break;

        }

    }



    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Status status;



    
    double totaltime = 0;

    if(pid == 0){
      printdata(data,N,Nlocal);

      for(int ipid=1; ipid < np;++ipid){
        double start= MPI_Wtime();

        MPI_Recv(&data[0], N*Nlocal, MPI_DOUBLE, ipid,tag, MPI_COMM_WORLD, &status);

        double end = MPI_Wtime();
        totaltime =end - start;

        printdata(data,N,Nlocal);

        std::cout << 1*sizeof(double) << "\t" << totaltime<< "\t" << 1*sizeof(double)/totaltime/1.0e6 <<std::endl;


      }
    } else {
        MPI_Send(&data[0], N*Nlocal, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }
     MPI_Finalize();
   
    return 0;
}

void printdata(std::vector<double>data,int N, int Nlocal){

    for (int i=0; i<Nlocal;++i){
        for(int j=0; j<N; j++){
          std::cout << data[(i*N)+j]<<"\t";
        }
        std::cout<<"\n";
    }
}
