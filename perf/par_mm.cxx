#include <omp.h>
#include <cstdio>
#include <ctime>

// i for rows, j for columns
#define C(i,j) C[i+N*j]
#define A(i,j) A[i+N*j]
#define B(i,j) B[i+N*j]

int main() {

    clock_t clocks;

    int N = 2000;
    auto A = new double[N*N];
    auto B = new double[N*N];
    auto C = new double[N*N];
    int strip = 8;

    // Set the number of threads
    int nthreads = 4;
    omp_set_num_threads(nthreads);

    // Loop nest: JKI
    clocks = clock();
#pragma omp parallel for
    for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++)
            for (int i = 0; i < N; i++)
                C(i,j) += A(i,k)*B(k,j);

    clocks = clock() - clocks;  // timing

    printf("[%d threads] JKI, N = %d, Execution time = %ld ms\n",
        nthreads, N, clocks * 1000/CLOCKS_PER_SEC);

    // cleanup
    delete [] A;
    delete [] B;
    delete [] C;

    return 0;
}
