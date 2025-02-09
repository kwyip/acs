function K = toK(H)
K = H^2;
K = [real(K) -imag(K);imag(K) real(K)];