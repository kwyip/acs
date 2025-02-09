function Phi0 = toPhi(H,psi0)
Phi0 = [tovp(psi0); todvp(H,psi0)];