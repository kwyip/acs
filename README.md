[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=kwyip/acs&file=README.md)

**Analog Classical Simulation**
---
This repository contains the following functions:

1. [`toK.m`](https://github.com/kwyip/acs_test/blob/main/toK.m) inputs a Hamiltonian H and outputs the K matrix, i.e, turn H real.

2. [`toX.m`](https://github.com/kwyip/acs_test/blob/main/toX.m) inputs a Hamiltonian H and outputs the X matrix for the ODE solver.

3. [`tovp.m`](https://github.com/kwyip/acs_test/blob/main/tovp.m) inputs a state vector psi and outputs the varphi vector, i.e, turn psi real.

4. [`todvp.m`](https://github.com/kwyip/acs_test/blob/main/todvp.m) inputs a Hamiltonian H, a state vector psi and outputs the velocity vector.
   
5. [`toPhi.m`](https://github.com/kwyip/acs_test/blob/main/toPhi.m) inputs a Hamiltonian H, a state vector psi and outputs the Phi vector for the ODE solver.

6. [`topsi.m`](https://github.com/kwyip/acs_test/blob/main/topsi.m) inputs a Phi vector from the ODE solver and outputs the a state vector psi.
   
7. [`dPhi.m`](https://github.com/kwyip/acs_test/blob/main/dPhi.m) is an ode function for solving a 2nd order ODE.

8. [`qaoa_p2_n8.m`](https://github.com/kwyip/acs_test/blob/main/qaoa_p2_n8.m) is an executable function for an analog classical simulation of QAOA.
