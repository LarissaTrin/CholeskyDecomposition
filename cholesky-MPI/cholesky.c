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

void show_matrix(double **A, int n)
{
  FILE *file = fopen("cholesky.out", "w");
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      fprintf(file, "%2.5f ", A[i][j]);
    fprintf(file, "\n");
  }
  fclose(file);
}

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

double **mat_aloc(int n)
{
  int i;
  double **x = malloc(sizeof(double *) * n);
  assert(x != NULL);

  for (i = 0; i < n; i++)
  {
    x[i] = malloc(sizeof(double) * n);
    assert(x[i] != NULL);
  }

  return x;
}

void mat_gen(FILE *file, double **s, int n)
{
  int i, j;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      fscanf(file, "%lf", &s[i][j]);
    }
  }
}

void mat_del(double **x)
{
  free(x[0]);
  free(x);
}

void printMat(double **x, int n, int p)
{
  printf("p = %d\n", p);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      printf("%f   ", x[i][j]);
    }
    printf("\n");
  }
}

double **matCopy(double **x, double **y, int n)
{
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      x[i][j] = y[i][j];
    }
  }
  return x;
}

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
    file = fopen("cholesky_5000.in", "r");
    fscanf(file, "%d", &n);
    ;
  }
  // envia valor de n para todos os processos
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); // OKAY

  double **local_A = mat_new(n);
  if (my_rank == 0)
  {
    mat_gen(file, local_A, n);
    fclose(file);
  }
  // printf("p = %d\n",my_rank); //OKAY

  // envia valor A para todos os processos OKAY
  if (my_rank == 0)
  {
    for (int i = 1; i < p; i++)
    {
      MPI_Send(&(local_A[0][0]), (n + n) * n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
    // printMat(local_A, n,my_rank);
  }
  else
  {
    MPI_Recv(&(local_A[0][0]), (n + n) * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    // printMat(local_A, n,my_rank);
  }

  //----------------------------------------------
  //------------------Solve Cholesky--------------
  // MPI_Barrier(MPI_COMM_WORLD); /* Timing */
  if (my_rank == 0)
  {
    gettimeofday(&tstart, NULL);
  }

  for (i = 0; i < n; i++)
  { // coluna
    // printf("\ni persent p = %d percent %d = %d\n",i,p,i%p);
    // printf("my_rank = %d\n",my_rank);
    if (my_rank == 0) // diagonal
    {
      s = 0;
      for (k = 0; k < i; k++)
      {
        s += local_A[i][k] * local_A[i][k];
      }
      local_A[i][i] = sqrt(local_A[i][i] - s);
    }

    MPI_Bcast(local_A[i], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // printMat(local_A, n,my_rank);

    // printf("COLUNA: %d   my_rank = %d\n",i,my_rank);
    n_bar = round(((n - 1) - i) / p);
    // printf("n_bar: %d\n",n_bar);
    inicio = calc_inicio(n_bar, i, my_rank);
    // printf("inicio: %d\n",inicio);
    fim = calc_fim(n_bar, n, i, p, my_rank);
    // printf("fim: %d\n\n",fim);
    // MPI_Barrier(MPI_COMM_WORLD); /* Timing */
    if (inicio <= fim)
    {
      // printf("COLUNA: %d   my_rank = %d\n", i, my_rank);
      // printf("     inicio: %d   fim: %d\n", inicio, fim);
      for (j = inicio; j <= fim; j++)
      { // linha
        // if (j < inicio && j > fim) {
        s = 0;
        for (k = 0; k < i; k++)
        {
          s += local_A[j][k] * local_A[i][k];
        }
        local_A[j][i] = (1.0 / local_A[i][i] * (local_A[j][i] - s));
        // local_A[i][j] = local_A[j][i];
        // }

        // MPI_Bcast(local_A[j], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        local_A[i][j] = local_A[j][i];
      }
    }

    if (my_rank != 0)
    {
      // printMat(local_A, n, my_rank);
      for (int a = inicio; a <= fim; a++)
      {
        MPI_Send(&local_A[a][i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&local_A[i][a], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      }
    }
    else
    {
      for (int k = 1; k < p; k++)
      {
        inicio = calc_inicio(n_bar, i, k);
        fim = calc_fim(n_bar, n, i, p, k);
        for (int a = inicio; a <= fim; a++)
        {
          MPI_Recv(&local_A[a][i], 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
          MPI_Recv(&local_A[i][a], 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
        }
      }
      
    }
    // MPI_Barrier(MPI_COMM_WORLD); /* Timing */
    // envia valor A para todos os processos OKAY
    if (my_rank == 0)
    {
      for (int i = 1; i < p; i++)
      {
        MPI_Send(&(local_A[0][0]), (n + n) * n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      }
      // printMat(local_A, n);
    }
    else
    {
      MPI_Recv(&(local_A[0][0]), (n + n) * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
      // printMat(local_A, n);
    }
  }

  // MPI_Barrier(MPI_COMM_WORLD); /* Timing */

  if (my_rank == 0)
  {
    // printMat(local_A, n, my_rank);
    gettimeofday(&tend, NULL);
    t_solve = (tend.tv_sec - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec) / 1000000.0;
    // printf("Solve: %fs    ", t_solve);
    printf("%f\n", t_solve);

    //---------------Register Matrix---------------
    show_matrix(local_A, n);
    //---------------------------------------------
    printf("\nExecucao finalizada\n");
  }
  mat_del(local_A);
  MPI_Finalize();
  return 0;
}
