#include <hip/hip_runtime.h>
#include <chrono> // timer
#include <cassert>
#include <cstdio>

// i for rows, j for columns
#define C(i, j) C[i + N * j]
#define A(i, j) A[i + N * j]
#define B(i, j) B[i + N * j]

// Innermost loop k
__global__ void mmKernel(double *C, double *A, double *B, int N) {
  int i = threadIdx.x + blockDim.x * blockIdx.x;
  int j = threadIdx.y + blockDim.y * blockIdx.y;

  for (int k = 0; k < N; k++)
    C(i, j) += A(i, k) * B(k, j);
}

int main() {

  using namespace std::chrono;
  using Timer = high_resolution_clock;

  int N = 4096;      // matrix size
  double *A, *B, *C; // matrices
  hipMalloc(&A, N * N * sizeof(double));
  hipMalloc(&B, N * N * sizeof(double));
  hipMalloc(&C, N * N * sizeof(double));

  // Parameters for tiling
  // Each of the thread compute one output
  // grid size * block size == N*N
  dim3 block(16, 16);                    // inner loops i,j
  dim3 grid((N + block.x - 1) / block.x, // outer loops ii,jj
            (N + block.y - 1) / block.y);
  assert(grid.x * block.x == N);
  assert(grid.y * block.y == N);

  printf("Grid size = (%d, %d), Block size = (%d, %d)\n", grid.x, grid.y, block.x, block.y);

  // SIMT
  auto start_time = Timer::now(); // record the start time

  // Kernel launch
  mmKernel<<<grid, block>>>(C, A, B, N);
  hipDeviceSynchronize();

  auto stop_time = Timer::now(); // record the stop time
  auto diff_time = duration_cast<milliseconds>(stop_time - start_time);

  printf("SIMT, N = %d, Execution time = %ld ms\n", N, diff_time.count());

  // cleanup
  hipFree(A);
  hipFree(B);
  hipFree(C);

  return 0;
}
