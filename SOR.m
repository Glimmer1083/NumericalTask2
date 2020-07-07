function [X, k, relerr] = SOR(A, b, eps, Iteration, omega)
n = length(b);
X0 = zeros(n, 1);
X = zeros(n, 1);
relerr = zeros(Iteration+1, 1);
for k = 1:Iteration
   for i = 1:n
       X(i) = X0(i) + omega * (b(i) - A(i, 1:i-1) * X(1:i-1) - A(i, i:n) * X0(i:n)) / A(i, i);
   end
   relerr(k) = norm(X - X0) / norm(X0);
   if max(abs(X - X0)) < eps
       break
   else
       X0 = X;
   end
end

end