
# 稳定性

考虑以下初值问题，

```math
\left\{
    \begin{aligned}
        &\frac{du}{dt}+Au=0,\quad 0\le t \le T, \\
        &u(0) = b.
    \end{aligned}
\right.
```

其中，$`A`$、$`T`$、$`b`$都是常数。

使用有限差分来求解它，把时间导数用一阶差商替换，得到差分格式（一阶前向差分），

```math
\begin{aligned}
&\frac{u_{n+1}-u_{n}}{\Delta t} + Au_{n} = 0, \\
&u_0 = b.
\end{aligned}
```

这个算法只在$`\Delta t \le 2/A`$时才是稳定的，如果不满足这个条件，数值解会随时间增加发散到无穷大。

令$`A=4`$，那么必须要有$`\Delta t \le 1/2`$。从代码运行结果可以看到，第一批结果是正常的，第二批结果完全发散了。下面一段代码演示了这一点。

> 编译代码
>
> 首先把代码保存到文件中，名为`initial-value.cxx`。
>
> 然后使用`g++`编译，生成一个名为`initial-value`的可执行文件。
>
> `g++ initial-value.cxx -o initial-value`

```C
// Initial value problem
//
// ddt(u) + Au = 0,  0 <= t <= T
// u(0) = b
//
// Scheme:
// u_{n+1} - u_{n} + dt*A*u_{n} = 0

#include <vector>
#include <cstdio>

int main() {

    double A = 4;      // coefficient
    double b = 0.5;    // boundary condition
    double T = 10;     // upper limit
    double dt;         // time step
    double u;          // value

    dt = 0.4;
    u = b;
    printf("A = %f, dt = %f, stable? %d\n", A, dt, dt <= 2/A);
    for (double t = 0; t < T; t+=dt) {
        printf("t = %f, u = %f\n", t, u);
        u = u - dt * A * u;
    }

    printf("\n");

    dt = 0.6;
    u = b;
    printf("A = %f, dt = %f, stable? %d\n", A, dt, dt <= 2/A);
    for (double t = 0; t < T; t+=dt) {
        printf("t = %f, u = %f\n", t, u);
        u = u - dt * A * u;
    }

    return 0;
}
```

或者用下面的 Python 代码。

> 新建一个`initial-value.py`文件保存代码。
>
> 运行代码：
>
> `python3 initial-value.py`

```python
# Initial value problem
#
# ddt(u) + Au = 0,  0 <= t <= T
# u(0) = b
#
# Scheme:
# u_{n+1} - u_{n} + dt*A*u_{n} = 0
import numpy as np

A = 4       # coefficient
b = 0.5     # boundary condition
T = 10      # upper limit

dt = 0.4    # time step
u = b
print(f"A = {A}, dt = {dt}, stable? {dt<=2/A}")
for t in np.arange(0, T, dt):
    print(f"t = {t:.3}, u = {u:.5f}")
    u = u - dt*A*u

print();

dt = 0.6
u = b
print(f"A = {A}, dt = {dt}, stable? {dt<=2/A}")
for t in np.arange(0, T, dt):
    print(f"t = {t:.3}, u = {u:.5f}")
    u = u - dt*A*u
```

# 程序性能

考虑一个矩阵乘法，

```math
C=A*B
```

其中的三个矩阵都是$`N\times N`$。

简单的矩阵乘法就是按定义实现，
```C
for (i = 0; i<N; i++)
    for (j = 0; j<N; j++)
        for (k = 0; k<N; k++)
            C[i][j] += A[i][k] * B[k][j];
```

做一次矩阵乘法一共有$N^3$次迭代。下面的代码演示的是，仅仅改变循环嵌套的顺序，就可以让计算时间发生巨大变化。

> 编译代码
>
> 新建文件`mm.cxx`保存下面的代码并编译，生成一个名为`mm`的可执行文件。
>
> `g++ mm.cxx -o mm`

```C
#include <cstdio>
#include <ctime>

// i for rows, j for columns
#define C(i,j) C[i+N*j]
#define A(i,j) A[i+N*j]
#define B(i,j) B[i+N*j]

int main() {

    clock_t clocks;

    long N = 1000;
    auto A = new double[N*N];
    auto B = new double[N*N];
    auto C = new double[N*N];

    long i, j, k;

    // Loop nest: IJK
    clocks = clock();
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            for (k = 0; k < N; k++)
                C(i,j) += A(i,k)*B(k,j);

    clocks = clock() - clocks;  // timing

    printf("IJK, Execution time = %ld ms\n", clocks * 1000/CLOCKS_PER_SEC);

    // Loop nest: IKJ
    clocks = clock();
    for (i = 0; i < N; i++)
        for (k = 0; k < N; k++)
            for (j = 0; j < N; j++)
                C(i,j) += A(i,k)*B(k,j);

    clocks = clock() - clocks;  // timing

    printf("IKJ, Execution time = %ld ms\n", clocks * 1000/CLOCKS_PER_SEC);

    // Loop nest: JKI
    clocks = clock();
    for (j = 0; j < N; j++)
        for (k = 0; k < N; k++)
            for (i = 0; i < N; i++)
                C(i,j) += A(i,k)*B(k,j);

    clocks = clock() - clocks;  // timing

    printf("JKI, Execution time = %ld ms\n", clocks * 1000/CLOCKS_PER_SEC);

    // cleanup
    delete [] A;
    delete [] B;
    delete [] C;

    return 0;
}
```

# 并行计算

## 向量化

并行计算有很多种，其中一种简单的形式是向量化（vectorization），也就是让一个处理器同时做多个标量运算。

例如，把两个长度为$N$的向量相加，下面的代码看起来需要$N$次运算。

```c
for (int i = 0; i < N; i++)
    c[i] = a[i] + b[i];
```

实际上，经过向量化之后（vectorized），有可能只需要$N/4$次运算或者更少，每次运算都是把多个数字当成一个很短的向量一起算。

仍然用前面的矩阵乘法代码`mm.cxx`，不过在编译时我们启用向量化功能。执行代码，可以看到执行时间变短。

```bach
g++ -O3 -fopt-info mm.cxx -o vec_mm
```

> 这里的`-fopt-info-vec`是让编译器打印一条消息，告诉我们有没有成功向量化。

## 多线程

另一种并行的技术是多线程（multithreading），利用多核处理器做并行计算。前面都是利用一个处理器，现在我们要利用多个处理器，因此可以计算更大的矩阵。

例如，下面的代码让 4 个线程同时计算$N=2000$的情况。如果用前面的`mm.cxx`来算，时间会很长，但使用多线程之后，会明显变快。

> 编译代码
>
> 先把代码拷贝到文件`par_mm.cxx`中。
>
> 这里使用了 OpenMP 为我们生成多线程的代码，编译时使用下面的命令。
>
> `g++ -fopenmp par_mm.cxx -o par_mm`

```C
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
    long strip = 8;

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

    printf("JKI, Execution time = %ld ms\n", clocks * 1000/CLOCKS_PER_SEC/nthreads);

    // cleanup
    delete [] A;
    delete [] B;
    delete [] C;

    return 0;
}
```

## 多线程+向量化

多线程通常可以把外层循环并行，向量化通常可以把内层循环并行，这两个可以合起来使用。只要加上两个选项就行。

```bash
g++ -fopenmp -O3 par_mm.cxx -o par_vec_mm
```
