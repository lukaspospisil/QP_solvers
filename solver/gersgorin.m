function [ normA ] = gersgorin( A )

S = diag(A);
R = sum(abs(A-diag(S)),2);

normA = max(S+R);

end

