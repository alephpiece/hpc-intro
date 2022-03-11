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

  double A = 4;   // coefficient
  double b = 0.5; // boundary condition
  double T = 10;  // upper limit
  double dt;      // time step
  double u;       // value

  dt = 0.4;
  u = b;
  printf("A = %f, dt = %f, stable? %d\n", A, dt, dt <= 2 / A);
  for (double t = 0; t < T; t += dt)
  {
    printf("t = %f, u = %f\n", t, u);
    u = u - dt * A * u;
  }

  printf("\n");

  dt = 0.6;
  u = b;
  printf("A = %f, dt = %f, stable? %d\n", A, dt, dt <= 2 / A);
  for (double t = 0; t < T; t += dt)
  {
    printf("t = %f, u = %f\n", t, u);
    u = u - dt * A * u;
  }

  return 0;
}
