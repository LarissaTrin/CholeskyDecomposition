#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <mpi.h>

void mat_zero(double **x, int n)
{
  int i, j;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      x[i][j] = 0.0;
    }
  }
}

double **mat_new(int n)
{
  int i;
  double **x = malloc(sizeof(double *) * n);
  assert(x != NULL);

  for (i = 0; i < n; i++)
  {
    x[i] = malloc(sizeof(double) * n);
    assert(x[i] != NULL);
  }

  mat_zero(x, n);

  return x;
}

void print_mat(double **x, int n, int p)
{
  printf("p = %d\n",p);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      printf("%f   ",x[i][j]);
    }
    printf("\n");
  }
  
}

int calc_fim(int n_bar, int n, int i, int p, int my_rank)
{
  // caso a razão da quantidade de linhas abaixo do pivo e a quantidade de processo
  // tiver resto diferente de zero, o último processso ficará com as linhas restantes
  int fim;
  if (((((n - 1) - i) % p) != 0) && (my_rank == (p - 1)))
  {
    fim = n - 1;
  }
  else
  {
    fim = (n_bar * (my_rank + 1)) + i;
  }
  return fim;
}

int calc_inicio(int n_bar, int i, int my_rank)
{
  int inicio;
  inicio = ((n_bar * my_rank) + 1) + i;
  return inicio;
}

int main(int argc, char *argv[]){

    int n_bar, inicio, fim;
    int p, my_rank;
    MPI_Status status;


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); //devolve número do processo 
    MPI_Comm_size(MPI_COMM_WORLD, &p); //numero de processos que foram criados
    MPI_Barrier(MPI_COMM_WORLD); /* Timing */
    int n=4;

    double **A=mat_new(n);

    if (my_rank==0)
    {
        int cont =1;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                A[i][j]=cont++;
            }
        }
        // print_mat(A,n,my_rank);
    }
    // 
    print_mat(A,n,my_rank);

    // mandar um vetor de dados em vez de 1 por 1
    // envia valor A para todos os processos OKAY
    // if (my_rank == 0)
    // {
    //     for (int i = 1; i < p; i++)
    //     {
    //         MPI_Send(&(A[0][0]), (n + n)*n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    //     }
    //     // printMat(local_A, n);
    // }
    // else
    // {
    //     MPI_Recv(&(A[0][0]), (n + n)*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    //     // printMat(local_A, n);
    // }

    // mandar um vetor de dados em vez de 1 por 1
    //enviar a coluna inteira atraves do for
    // for (int i = 0; i < n; i++)
    // {
    //     if (my_rank==0)
    //     {
    //         MPI_Send(&A[i][0], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    //     } else {
    //         MPI_Recv(&A[i][0], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    //     }
    // }


    // mandar um vetor de dados em vez de 1 por 1
    //enviar a linha inteira especifica
    // if (my_rank==0)
    // {
    //   MPI_Send(&(A[1][0]), 1*n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    // } else {
    //   MPI_Recv(&(A[1][0]),1*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    // }

    n_bar = round(((n-1)-0) / p);
    inicio = calc_inicio(n_bar, 0, my_rank);
	  // printf("inicio: %d\n",inicio);
	  fim = calc_fim(n_bar, n, 0, p, my_rank);
    // printf("fim: %d\n\n",fim);
    
    // envia a colula em quantidade de linhas especificas
    // if (my_rank!=0){
    //     for (int i = inicio; i <= fim; i++)
    //     {
    //         MPI_Send(&A[i][0], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    //         printf("\nEnviado para 0\n");
    //     }
    // } else {
    //     printf("\nTo em 0\n");
    //     for (int k = 1; k < p; k++)
    //     {
    //         printf("\nTo em 0.1\n");
    //         inicio = calc_inicio(n_bar, 0, k);
	// 		fim = calc_fim(n_bar, n, 0, p, k);
    //         for (int i = inicio; i <= fim; i++){
    //             printf("\nTo em 0.2\n");
    //             MPI_Recv(&A[i][0], 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
    //         }
    //     }
    // }

    // envia a linha em quantidade de linhas especificas
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank!=0){
        // printf("\nENVIANDO para 0\n");
        MPI_Send(&A[0][inicio], (fim-inicio+1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        // printf("\nEnviado para 0\n");
    } else {
        // printf("\nTo em 0\n");
        for (int k = 1; k < p; k++)
        {
            // printf("\nTo em 0.1\n");
            inicio = calc_inicio(n_bar, 0, k);
			      fim = calc_fim(n_bar, n, 0, p, k);
            // printf("\nTo em 0.2\n");
            MPI_Recv(&A[0][inicio], (fim-inicio+1), MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
            // print_mat(A,n,my_rank);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (my_rank==0)
    {
      print_mat(A,n,my_rank);
    }
    
    // print_mat(A,n,my_rank);
    
    MPI_Barrier(MPI_COMM_WORLD); /* Timing */
    MPI_Finalize();
    return 0;
}