function psi = topsi(Phi)
neval = length(Phi)/4;
psi = Phi(1:neval)+1i*Phi(neval+1:neval*2);