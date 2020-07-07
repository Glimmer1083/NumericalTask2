function [X] = gauss(A, b, mode)
n = length(b);
X = zeros(n, 1);
% 顺序消元
if mode == 0
    for k = 1:n - 1
        for i = k+1:n
            A(i,k) = A(i,k) / A(k,k);
            A(i,k+1:n) = A(i,k+1:n) - A(i,k) * A(k,k+1:n);
            b(i) = b(i) - A(i,k) * b(k);
        end
    end
% 选最大主元
elseif mode == 1
    for k = 1:n - 1
        [value, index] = max(abs(A(k:n,k))); index = index + k - 1;
        tmpA = A(k,:); A(k,:) = A(index,:); A(index,:) = tmpA;
        tmpb = b(k); b(k) = b(index); b(index) = tmpb;
        for i = k+1:n
            A(i,k) = A(i,k) / A(k,k);
            A(i,k+1:n) = A(i,k+1:n) - A(i,k) * A(k,k+1:n);
            b(i) = b(i) - A(i,k) * b(k);
        end
    end
% 选最小不为0主元
elseif mode == 2
    for k = 1:n - 1
        arr = abs(A(k:n,k));
        [value, index] = min(arr(arr~=0)); index = index + k - 1;
        tmpA = A(k,:); A(k,:) = A(index,:); A(index,:) = tmpA;
        tmpb = b(k); b(k) = b(index); b(index) = tmpb;
        for i = k+1:n
            A(i,k) = A(i,k) / A(k,k);
            A(i,k+1:n) = A(i,k+1:n) - A(i,k) * A(k,k+1:n);
            b(i) = b(i) - A(i,k) * b(k);
        end
    end
end

% 回代
X(n) = b(n) / A(n,n);
for i = n-1:-1:1
   X(i) = (b(i) - A(i, i+1:n) * X(i+1:n)) / A(i, i);
end

end