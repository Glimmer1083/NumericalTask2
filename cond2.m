function cond2A = cond2(A)
B = A'*A;
[V1, D1] = eig(B);
[V2, D2] = eig(B^(-1));
cond2A = sqrt(max(max(D1))) * sqrt(max(max(D2)));
end