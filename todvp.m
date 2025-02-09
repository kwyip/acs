function varphidot0 = todvp(H,psi0)
varphidot0 = [imag(H)*real(psi0) + real(H)*imag(psi0); -real(H)*real(psi0) + imag(H)*imag(psi0)];