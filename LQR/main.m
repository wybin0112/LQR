% 清空
clc;
clear;
close all;

% pkg load control;

%% 系统定义
A = [0 1; -1 -0.5];
[nA, ~] = size(A);
B = [0; 1];
mB = size(B, 2);
C = [0 0];
D = 0;
%% 系统离散化
Ts = 0.1; % 步长
sys_d = c2d(ss(A,B,C,D), Ts); % 连续系统转离散化
A = sys_d.A;
B = sys_d.B;

% 状态初始化
x0 = [1; 0];
x = x0;

% 输入初始化
u0 = -2;
u = u0;

% 设置权重矩阵
R = [0.01];
Q = [1 0; 0 1];
S = [1 0; 0 1];

% 定义系统运行步数
k_steps = 100;

% 定义x_history零矩阵，存储系统状态
x_history = zeros(nA, k_steps);
x_history(:, 1) = x0;

% 定义u_history零矩阵，存储输入；
u_history = zeros(mB, k_steps);
u_history(:, 1) = u0;

% 计算P矩阵和F矩阵
P0 = S;
%[P, F] = solvePF(A, B, Q, R, P0, k_steps);

% 计算F常值矩阵
F = LQR_solveF(A, B, Q, R, P0);

% 仿真开始，for循环
for i = 1:1:k_steps
    % 计算系统输入
    % u = -F((i-1)*mB+1:i*mB, :)*x;
    u = -F*x; % F常值；
    % 计算系统响应
    x = A*x +B*u;
    % 存储
    x_history(:, i+1) = x;
    u_history(:, i+1) = u;
end

figure(1)
subplot(2, 1, 1)
for i = 1:nA
    plot(x_history(i, :));
    hold on;
end
subplot(2, 1, 2)
for i = 1:mB
    plot(u_history(i, :));
    hold on;
end



















