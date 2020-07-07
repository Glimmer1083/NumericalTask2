
%% 5.1 主元的选取与算法的稳定性      
% 1. 准备变量
clear;
close all;
n = 20;
x = ones(n, 1);
% mode = input('请选择模式: 0 顺序消元, 1 列选最大主元, 2 列选最小主元 ');

% 2. 程序开始
alpha = 6; beta = 8; gamma = 1;
% 定义三对角矩阵A
A = diag(alpha*ones(n, 1)) + diag(beta*ones(n-1, 1), -1) + diag(gamma*ones(n-1, 1), 1);
b = A * x;
C0 = cond(A, 1);
% 顺序消元
X0 = gauss(A, b, 0);
% 列主元
X1 = gauss(A, b, 1);
% 列最小
X2 = gauss(A, b, 2);
fig1 = figure();
plot(1:n, X2 - x, 1:n, X1 - x, 'LineWidth', 2);
legend('模最小', '列主元');
grid on;
title('模最小和列主元的误差比较');
set(gca, 'linewidth', 1.5, 'fontsize', 14, 'fontname', 'Times New Roman');
saveas(fig1, 'fig1.png')

%% 5.2 线性代数方程组的性态与条件数的估计
clear;
close all;
% 定义变量
n = 10;
factor = 1e-6;
alpha = 6; beta = 8; gamma = 1;

% 求扰动后的误差    
A = diag(alpha*ones(n, 1)) + diag(beta*ones(n-1, 1), -1) + diag(gamma*ones(n-1, 1), 1);
x = ones(n, 1);
b = A * x;
delta_A = rand(n) * factor; Abar = A + delta_A;
delta_b = rand(n, 1) * factor; bbar = b + delta_b;
C = cond(A, 1);
xbar = Abar \ bbar;
% 计算误差
actEps = norm(xbar - x, 1) / norm(x, 1);

%% 5.2.2)比较函数condest用的机器时间, 比较cond2
max = 500;
time = zeros(max, 1);
mycond2List = zeros(max, 1);
cond2List = zeros(max, 1);
for N = 1:1:max
    CA = fix(100*randn(N));
    tic
    condestCA = condest(CA);
    toc
    time(N) = toc;
%     mycond2List(N) = cond2(CA);
%     cond2List(N) = cond(CA, 2);
end
for N = 1:1:max
    CA = fix(100*randn(N));
    mycond2List(N) = cond2(CA);
    cond2List(N) = cond(CA, 2);
end
fig2 = figure();
plot(time, 'LineWidth', 2);
grid on;
title('condest所需机器时间');
xlabel('n')
ylabel('t/s')
set(gca, 'linewidth', 1.5, 'fontsize', 14, 'fontname', 'Times New Roman');
saveas(fig2, 'fig2.png')

fig3 = figure();
plot(1:max, abs(mycond2List-cond2List), 'LineWidth', 2);
legend('eig计算', 'cond计算');
grid on;
title('cond2比较');
set(gca, 'linewidth', 1.5, 'fontsize', 14, 'fontname', 'Times New Roman');
saveas(fig3, 'fig3.png')

%% 5.2.3)理论计算误差限
condestA = condest(A);
InvA = condestA / norm(A, 1);
theoryEps = condestA / (1 - InvA * norm(delta_A, 1)) * (norm(delta_A, 1) / norm(A, 1) + norm(delta_b', 1) / norm(b', 1));
diff = theoryEps - actEps;

%% 5.2.4)Hilbert矩阵条件数
cond1H = ones(1,15); cond2H = ones(1,15); condInfH = ones(1,15);
N = 15;
for n = 1:N
    hilbert = hilb(n);           
    cond1H(n) = cond(hilbert,1);
    cond2H(n) = cond(hilbert,2);
    condInfH(n) = cond(hilbert,inf);
end

fig4 = figure();
plot(1:N, cond1H, 1:N, cond2H, 1:N, condInfH, 'LineWidth', 2);
legend('1-范数', '2-范数', '无穷范数');
grid on;
title('Hilbert矩阵的条件数');
set(gca, 'linewidth', 1.5, 'fontsize', 14, 'fontname', 'Times New Roman');
saveas(fig4, 'fig4.png')

%% 6.1 病态的线性方程组的求解
clear;
close all;
N = 20;
H = hilb(N);
x_exact = ones(N, 1);
b = H * x_exact;
x_gauss_1 = gauss(H, b, 0); gauss_err_1 = x_gauss_1 - x_exact;
x_gauss_2 = gauss(H, b, 1); gauss_err_2 = x_gauss_2 - x_exact;
i = 1:N;

[x_jacobi, kj]= Jacobi(H, b, 1e-5, 100);
[x_GS, kgs]= GS(H, b, 1e-7, 10000, 0);
[x_SOR, ksor]= SOR(H, b, 1e-7, 10000, 1.5);

fig5 = figure();
plot(i, gauss_err_1, i, gauss_err_2, i, x_GS - x_exact, i, x_SOR - x_exact, 'LineWidth', 2);
legend('顺序消元', '列主元', 'GS', 'SOR', 'Location', 'northwest')
grid on;
title('Hilbert矩阵解的比较');
set(gca, 'linewidth', 1.5, 'fontsize', 14, 'fontname', 'Times New Roman');
saveas(fig5, 'fig5.png')

%% 6.2 Jacobi、GS、SOR、BSOR收敛速度
N = 10;
h = 1 / (N + 1);
A = zeros(N^2, N^2); b = zeros(N^2, 1); u = zeros(N^2, 1);
% 构造矩阵
for j = 1:N
    for i = 1:N
        origin = -h^2 * ((i*h)^2 + (j*h)^2) * exp(i*h*j*h);
        if i == 1 && j == 1
            b((j-1)*N+i) = 2 + origin;
        elseif i == 1 && j == N
            b((j-1)*N+i) = 1 + exp((i+1)*h) + origin;
        elseif i == N && j == N
            b((j-1)*N+i) = exp((i+1)*h) + exp((j+1)*h) + origin;
        elseif j == 1 && i == N
            b((j-1)*N+i) = 1 + exp((j+1)*h) + origin;
        elseif (i == 1 && (1 < j && j < N))
            b((j-1)*N+i) = 1 + origin;
        elseif (j == 1 && (1 < i && i < N))
            b((j-1)*N+i) = 1 + origin;
        elseif (i == N && (1 < j && j < N))
            b((j-1)*N+i) = exp(j*h) + origin;
        elseif (j == N && (1 < i && i < N))
            b((j-1)*N+i) = exp(i*h) + origin;
        else
            b((j-1)*N+i) = origin;
        end
        A((j-1)*N+i, (j-1)*N+i) = 4;
        if i ~= N
            A((j-1)*N+i, (j-1)*N+i+1) = -1;
        end
        if j ~= N
            A((j-1)*N+i, (j)*N+i) = -1;
        end
        if i ~= 1
            A((j-1)*N+i, (j-1)*N+i-1) = -1;
        end
        if j ~= 1
            A((j-1)*N+i, (j-2)*N+i) = -1;
        end
    end
end
% tu = A \ b;
% rtu = reshape(tu, N, []);
[X,Y] = meshgrid(1/(N+1):1/(N+1):1-1/(N+1));
u_exact = exp(X.*Y);

[tu_Jacobi, kj, err_jacobi]= Jacobi(A, b, 1e-6, 200);
[tu_GS, kgs, err_gs]= GS(A, b, 1e-6, 200, 0);
[tu_SOR, ksor, err_sor]= SOR(A, b, 1e-6, 200, 1.6);
[tu_BSOR, kbsor, err_bsor]= BSOR(A, b, 1e-6, 200, 1.6, 10);

fig6 = figure();
rtu = reshape(tu_BSOR, N, []);
surf(X, Y, rtu);
saveas(fig6, 'fig6.png')
uInf = norm(rtu - u_exact, inf);

fig7 = figure();
xgs = 1:kgs; xj = 1:kj; xsor = 1:ksor; xbsor = 1:kbsor;
semilogy(xgs, err_gs(xgs), '-ro', xj, err_jacobi(xj), '-b+', xsor, err_sor(xsor), '-kd', xbsor, err_bsor(xbsor), '-gd');
legend('GS', 'Jacobi', 'SOR(omega=1.6)', 'BSOR(omega=1.6)')
title('不同解法的收敛速度');
grid on;
set(gca, 'linewidth', 1.5, 'fontsize', 14, 'fontname', 'Times New Roman');
saveas(fig7, 'fig7.png')

%% 求解最佳迭代因子
omega = 1.5;
D = diag(diag(A));
L = D - tril(A);
U = D - triu(A);
L_omega = (D - omega * L) \ ((1 - omega) * D + omega * U);
rho_max_L_omega = zeros(1, 181);
for i = 1:200
    omega = i * 0.01;
    L_omega = (D - omega * L) \ ((1 - omega) * D + omega * U);
    rho_max = max(abs(eig(L_omega)));
    rho_max_L_omega(i) = rho_max;
end
fig8 = figure();
plot(1:200, rho_max_L_omega, 'LineWidth', 2);
title('谱半径随omega的变化');
grid on;
set(gca, 'linewidth', 1.5, 'fontsize', 14, 'fontname', 'Times New Roman');
saveas(fig8, 'fig8.png')