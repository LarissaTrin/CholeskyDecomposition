/*****************************************************
 * Site: https://rosettacode.org/wiki/Cholesky_decomposition
 *****************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <mpi.h>

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

void show_matrix(double *A, int n)
{
    FILE *file = fopen("cholesky.out", "w");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            fprintf(file, "%2.5f ", A[i * n + j]);
        fprintf(file, "\n");
    }
    fclose(file);
}

void mat_zero(double *x, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            x[i * n + j] = 0;
        }
    }
}

double *mat_new(int n)
{
    int i;
    i = n * n;
    double *x = malloc(sizeof(double) * i);

    mat_zero(x, n);

    return x;
}

void mat_gen(FILE *file, double *s, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            fscanf(file, "%lf", &s[i * n + j]);
        }
    }
}

void mat_del(double *x)
{
    // free(&x[0]);
    free(x);
}

// void printMat(double **x, int n, int p)
// {
//   printf("p = %d\n", p);
//   for (int i = 0; i < n; i++)
//   {
//     for (int j = 0; j < n; j++)
//     {
//       printf("%f   ", x[i][j]);
//     }
//     printf("\n");
//   }
// }


int main(int argc, char **argv)
{
  printf("Iniciando execucao\n");

  double t_solve;
  struct timeval tstart, tend;

  int n_bar;
  int p, my_rank;
  MPI_Status status;

  int inicio, fim;
  FILE *file;

  int n;

  double s = 0;
  int i, j, k;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  //---------------Get Matrix----------------
  if (my_rank == 0)
  {
    // pega a dimensao da matriz (n) e a matriz A
    file = fopen("cholesky_10000.in", "r");
    fscanf(file, "%d", &n);
    ;
  }
  // envia valor de n para todos os processos
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // OKAY

  double *local_A = mat_new(n);
  if (my_rank == 0)
  {
    mat_gen(file, local_A, n);
    fclose(file);
  }

  //----------------------------------------------
  //------------------Solve Cholesky--------------
  if (my_rank == 0)
  {
    gettimeofday(&tstart, NULL);
  }

  for (i = 0; i < n; i++)
  { // coluna
    if (my_rank == 0) // diagonal
    {
      s = 0;
      for (k = 0; k < i; k++)
      {
        s += local_A[i * n + k] * local_A[i * n + k];
      }
      local_A[i * n + i] = sqrt(local_A[i * n + i] - s);
    }

    MPI_Bcast(local_A, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    n_bar = round(((n - 1) - i) / p);
    inicio = calc_inicio(n_bar, i, my_rank);
    fim = calc_fim(n_bar, n, i, p, my_rank);
    for (j = inicio; j <= fim; j++)
    { // linha
      s = 0;
      for (k = 0; k < i; k++)
      {
        s += local_A[j * n + k] * local_A[i * n + k];
      }
      local_A[j * n + i] = (1.0 / local_A[i * n + i] * (local_A[j * n + i] - s));
      local_A[i * n + j] = local_A[j * n + i];
    }

    // if (inicio <= fim)
    // {
      
    // }

    if (my_rank != 0)
      {
        // printMat(local_A, n, my_rank);
        MPI_Send(&local_A[i * n + inicio], (fim - inicio + 1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      }
      else
      {
        for (int k = 1; k < p; k++)
        {
          inicio = calc_inicio(n_bar, i, k);
          fim = calc_fim(n_bar, n, i, p, k);
          MPI_Recv(&local_A[i * n + inicio], (fim - inicio + 1), MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
          for (int x = inicio; x <= fim; x++)
          {
            local_A[x * n + i] = local_A[i * n + x];
          }
        }
      }
  }


  if (my_rank == 0)
  {
    gettimeofday(&tend, NULL);
    t_solve = (tend.tv_sec - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec) / 1000000.0;
    printf("%f\n", t_solve);

    //---------------Register Matrix---------------
    show_matrix(local_A, n);
    //---------------------------------------------
  }

  // if (my_rank==1)
  // {
  //     show_matrix(local_A, n);
  // }

  printf("\nExecucao finalizada %d\n", my_rank);
  mat_del(local_A);
  MPI_Finalize();
  return 0;
}
