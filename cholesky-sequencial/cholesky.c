/*****************************************************
 * Site: https://rosettacode.org/wiki/Cholesky_decomposition
 *****************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

void cholesky(double **A, int n)
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < (i + 1); j++)
    {
      double s = 0;
      for (int k = 0; k < j; k++)
        s += A[i][k] * A[j][k];
      A[i][j] = (i == j) ? sqrt(A[i][i] - s) : (1.0 / A[j][j] * (A[i][j] - s));
    }
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

  double t_generate;
  double t_solve;
  double t_register;

  struct timeval tstart, tend;

  int n;

  //---------------Generate Matrix----------------
  gettimeofday(&tstart, NULL);
  FILE *file = fopen("cholesky_10000_t.in", "r");
  fscanf(file, "%d", &n);

  double **A = mat_new(n);
  mat_gen(file, A, n);
  fclose(file);
  gettimeofday(&tend, NULL);

  t_generate = (tend.tv_sec - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec) / 1000000.0;
  printf("Generate: %fs    ", t_generate);
  //--------------------------------------------

  //---------------Solve Cholesky---------------
  gettimeofday(&tstart, NULL);
  cholesky(A, n);
  gettimeofday(&tend, NULL);

  t_solve = (tend.tv_sec - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec) / 1000000.0;
  printf("Solve: %fs    ", t_solve);
  //---------------------------------------------

  //---------------Register Matrix---------------
  gettimeofday(&tstart, NULL);
  show_matrix(A, n);

  mat_del(A);
  gettimeofday(&tend, NULL);

  t_register = (tend.tv_sec - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec) / 1000000.0;
  printf("Register: %fs", t_register);
  //---------------------------------------------

  printf("\nExecucao finalizada\n");
  return 0;
}
