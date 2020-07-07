function [X, k, relerr] = Jacobi(A, b, eps, Iteration)
n = length(b);
X0 = zeros(n, 1);
X = zeros(n, 1);
relerr = zeros(Iteration+1, 1);
for k = 1:Iteration
   for i = 1:n
       X(i) = (b(i) - A(i, 1:i-1) * X0(1:i-1) - A(i, i+1:n) * X0(i+1:n)) / A(i, i);
   end
   relerr(k) = norm(X - X0) / norm(X0);
   if max(abs(X - X0)) < eps
       break
   else
       X0 = X;
   end
end

end