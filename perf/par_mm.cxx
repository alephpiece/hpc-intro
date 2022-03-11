#include <omp.h>
#include <chrono> // timer
#include <cstdio>
#include <ctime>

// i for rows, j for columns
#define C(i, j) C[i + N * j]
#define A(i, j) A[i + N * j]
#define B(i, j) B[i + N * j]

void par_mm(double *__restrict__ C,
            double *__restrict__ A,
            double *__restrict__ B,
            int N, int nthreads)
{
  using namespace std::chrono;
  using Timer = high_resolution_clock;

  // Set the number of threads
  omp_set_num_threads(nthreads);

  auto start_time = Timer::now(); // record the start time

#pragma omp parallel for
  for (int j = 0; j < N; j++)
    for (int k = 0; k < N; k++)
      for (int i = 0; i < N; i++)
        C(i, j) += A(i, k) * B(k, j);

  auto stop_time = Timer::now(); // record the stop time
  auto diff_time = duration_cast<milliseconds>(stop_time - start_time);

  printf("[%d threads] JKI, N = %d, Execution time = %ld ms\n",
         nthreads, N, diff_time.count());
}

int main() {

  int nthreads;
  int N = 4000;
  auto A = new double[N * N];
  auto B = new double[N * N];
  auto C = new double[N * N];

  nthreads = 16; // Number of threads = 16
  par_mm(C, A, B, N, nthreads);

  nthreads = 8; // Number of threads = 8
  par_mm(C, A, B, N, nthreads);

  nthreads = 4; // Number of threads = 4
  par_mm(C, A, B, N, nthreads);

  // cleanup
  delete[] A;
  delete[] B;
  delete[] C;

  return 0;
}