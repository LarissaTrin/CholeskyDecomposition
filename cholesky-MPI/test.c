#include<mpi.h>
#include<stdio.h>

int main(int argc, char *argv[]){
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //devolve número do processo 
    MPI_Comm_size(MPI_COMM_WORLD, &size); //numero de processos que foram criados
    printf("Eu sou %d de %d\n", rank, size);
    MPI_Finalize();
    return 0;
}