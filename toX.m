function X = toX(H)
K = toK(H);
sizeK = size(K,1);
X = [zeros(sizeK) -eye(sizeK);K zeros(sizeK)];