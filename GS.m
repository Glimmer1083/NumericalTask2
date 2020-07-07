function [X, k, relerr] = GS(A, b, eps, Iteration, mode)
n = length(b);
X0 = zeros(n, 1);
X = zeros(n, 1);
relerr = zeros(Iteration+1, 1);
if mode == 0
    for k = 1:Iteration
        for i = 1:n
           X(i) = (b(i) - A(i, 1:i-1) * X(1:i-1) - A(i, i+1:n) * X0(i+1:n)) / A(i, i);
        end
        relerr(k) = norm(X - X0) / norm(X0);
        if max(abs(X - X0)) < eps
           break
        else
           X0 = X;
        end
    end
elseif mode == 1
    for k = 1:Iteration        
        for i = 1:1:n
            if mod(i, 2) == 0
                X(i) = (b(i) - A(i, 1:i-1) * X0(1:i-1) - A(i, i+1:n) * X0(i+1:n)) / A(i, i);
            end
        end
        for i = 1:1:n
            if mod(i, 2) == 1
                X(i) = (b(i) - A(i, 1:i-1) * X(1:i-1) - A(i, i+1:n) * X(i+1:n)) / A(i, i);
            end
        end
        
        if max(abs(X - X0)) < eps
            break
        else
            X0 = X;
        end
    end
end
end