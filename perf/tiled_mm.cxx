#include <omp.h>
#include <chrono> // timer
#include <cstdio>
#include <ctime>

// i for rows, j for columns
#define C(i, j) C[i + N * j]
#define A(i, j) A[i + N * j]
#define B(i, j) B[i + N * j]

int main() {

  using namespace std::chrono;
  using Timer = high_resolution_clock;

  int N = 2048;
  auto A = new double[N * N];
  auto B = new double[N * N];
  auto C = new double[N * N];

  // Parameters for tiling
  int nblock_x = 16;
  int nblock_y = 16;

  // Set the number of threads
  int nthreads = 4;
  omp_set_num_threads(nthreads);

  auto start_time = Timer::now(); // record the start time

#pragma omp parallel for
  for (int jj = 0; jj < N; jj += nblock_y)
    for (int ii = 0; ii < N; ii += nblock_x)
      for (int j = jj; j < jj + nblock_y && j < N; j++)
        for (int i = ii; i < ii + nblock_x && i < N; i++)
          for (int k = 0; k < N; k++)
            C(i, j) += A(i, k) * B(k, j);

  auto stop_time = Timer::now(); // record the stop time
  auto diff_time = duration_cast<milliseconds>(stop_time - start_time);

  printf("[%d threads] JKI, N = %d, Execution time = %ld ms\n",
         nthreads, N, diff_time.count());

  // cleanup
  delete[] A;
  delete[] B;
  delete[] C;

  return 0;
}
