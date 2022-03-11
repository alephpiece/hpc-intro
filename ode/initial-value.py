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

print()

dt = 0.6
u = b
print(f"A = {A}, dt = {dt}, stable? {dt<=2/A}")
for t in np.arange(0, T, dt):
    print(f"t = {t:.3}, u = {u:.5f}")
    u = u - dt*A*u
