#include <cstdio>
#include <ctime>

// i for rows, j for columns
#define C(i, j) C[i + N * j]
#define A(i, j) A[i + N * j]
#define B(i, j) B[i + N * j]

int main() {

  clock_t clocks;

  int N = 1000;
  auto A = new double[N * N];
  auto B = new double[N * N];
  auto C = new double[N * N];

  int i, j, k;

  // Loop nest: IJK
  clocks = clock();
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      for (k = 0; k < N; k++)
        C(i, j) += A(i, k) * B(k, j);

  clocks = clock() - clocks; // timing

  printf("IJK, N = %d, Execution time = %ld ms\n", N, clocks * 1000 / CLOCKS_PER_SEC);

  // Loop nest: IKJ
  clocks = clock();
  for (i = 0; i < N; i++)
    for (k = 0; k < N; k++)
      for (j = 0; j < N; j++)
        C(i, j) += A(i, k) * B(k, j);

  clocks = clock() - clocks; // timing

  printf("IKJ, N = %d, Execution time = %ld ms\n", N, clocks * 1000 / CLOCKS_PER_SEC);

  // Loop nest: JKI
  clocks = clock();
  for (j = 0; j < N; j++)
    for (k = 0; k < N; k++)
      for (i = 0; i < N; i++)
        C(i, j) += A(i, k) * B(k, j);

  clocks = clock() - clocks; // timing

  printf("JKI, N = %d, Execution time = %ld ms\n", N, clocks * 1000 / CLOCKS_PER_SEC);

  // cleanup
  delete[] A;
  delete[] B;
  delete[] C;

  return 0;
}
