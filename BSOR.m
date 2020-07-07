function [X, k, relerr] = BSOR(A, b, eps, iter, omega, batch_size)
n = length(b);
m = n / batch_size;

X0 = zeros(n, 1);
X = zeros(n, 1);
relerr = zeros(iter+1, 1);
for k = 1:iter
    bx0 = X0(1:batch_size); Aii = A(1:batch_size, 1:batch_size); bii = b(1:batch_size);
    X(1:batch_size) = Aii \ (Aii * bx0 + omega * (bii - sum(A(1:batch_size, 1:n) * X0(1:n), 2)));
    for ii = 2:m
        st = (ii-1) * batch_size+1;
        ed = ii * batch_size;
        bx0 = X0(st:ed); Aii = A(st:ed, st:ed); bii = b(st:ed);
        bx = Aii \ (Aii * bx0 + omega * (bii - sum(A(st:ed, 1:st-1) * X(1:st-1), 2) - ...
            sum(A(st:ed, st:n) * X0(st:n), 2)));
        X(st:ed) = bx;
    end
    relerr(k) = norm(X - X0) / norm(X0);
    if max(abs(X - X0)) < eps
       break
    else
       X0 = X;
    end
end
end