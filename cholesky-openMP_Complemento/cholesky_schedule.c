/*****************************************************
 * Site: https://rosettacode.org/wiki/Cholesky_decomposition
 *****************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <omp.h>

void cholesky(double** A, int n) {
  double t_diagonal=0;
  double t_parallel=0;
  struct timeval tstart, tend;

  double s = 0;
  int diagonal=0;
  int i,j,k;

  for (i = 0; i < n; i++){ //coluna

    if (diagonal==0)
    {
      s = 0;
      diagonal=1;
      gettimeofday(&tstart, NULL);
      for (k = 0; k < i; k++){
        s += A[i][k] * A[i][k];
      }
      gettimeofday(&tend, NULL);
      A[i][i] = sqrt(A[i][i] - s);
      t_diagonal += (tend.tv_sec - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec) / 1000000.0;
      // printf("%f\n", t_diagonal);
    }
    
    gettimeofday(&tstart, NULL);
    #pragma omp parallel for shared(A,i) private(j,k,s) schedule(dynamic)
    for (j = i+1; j < n; j++) { //linha
      s = 0;
      for (k = 0; k < i; k++){
        s += A[j][k] * A[i][k];
      }
      A[j][i] = (1.0 / A[i][i] * (A[j][i] - s));
      A[i][j] = A[j][i];
    }
    gettimeofday(&tend, NULL);
    t_parallel += (tend.tv_sec - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec) / 1000000.0;
    // printf("%f\n", t_parallel);
    diagonal=0;
  }
  printf("Diagonal: %f\n", t_diagonal);
  printf("Abaixo da diafonal: %f\n", t_parallel);
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

int main()
{
  printf("Iniciando execucao\n");

    double t_solve;

    struct timeval tstart, tend;

    int n;

    //---------------Get Matrix----------------
    gettimeofday(&tstart, NULL);
    FILE *file = fopen("cholesky_5000.in", "r");
    fscanf(file, "%d", &n);

    double **A = mat_new(n);
    mat_gen(file, A, n);
    fclose(file);
    //--------------------------------------------

    //---------------Solve Cholesky---------------
    gettimeofday(&tstart, NULL);
    cholesky(A, n);
    gettimeofday(&tend, NULL);

    t_solve = (tend.tv_sec - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec) / 1000000.0;
    // printf("Solve: %fs    ", t_solve);
    printf("%f\n", t_solve);

    //---------------------------------------------

    //---------------Register Matrix---------------
    gettimeofday(&tstart, NULL);
    show_matrix(A, n);

    mat_del(A);
    //---------------------------------------------

    printf("\nExecucao finalizada\n");
  return 0;
}
