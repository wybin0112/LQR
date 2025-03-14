% 清空
clc;
clear;
close all;

% pkg load control;

% 定义系统参数
m_sys = 1;
b_sys = 0.5;
k_sys = 1;

%% 系统定义
AA = [0 1; -k_sys/m_sys -b_sys/m_sys];
[nA, ~] = size(AA);
BB = [0; 1/m_sys];
[nB, mB] = size(BB);
C = [0 0];
D = 0;

%% 系统离散化
Ts = 0.1; % 步长
sys_d = c2d(ss(AA,BB,C,D), Ts); % 连续系统转离散化
AA = sys_d.A;
BB = sys_d.B;

Ad = eye(nA); % xd 为常值

%% 增广矩阵形式
Aa = [AA zeros(nA, nA) BB; zeros(nA, nA) Ad zeros(nB,mB); zeros(mB, nA), zeros(mB, nA), eye(mB)];
Ba = [BB; zeros(nA*2-nB, mB); eye(mB)];
Ca = [eye(nA,nA) -eye(nA,nA) zeros(nA, mB)]; 

A = Aa;
[nA, ~] = size(A);
B = Ba;
mB = size(B, 2);
C = Ca;
D = 0;

% 预设期望
xd = [1; 0];
% % 求解ud,采用线性方程中方法，因为B可能不是方阵，无法求逆；
% ud = solve_ud(AA,BB,xd);

% 输入初始化
u0 = 0;
u = u0;

% 状态初始化
x0 = [0; 0];
x = x0;
xa = [x0; xd; u0];

% 设置权重矩阵
R = [0.1]*10;
Q = [1 0; 0 1];
S = [1 0; 0 1];

% 计算增广形式下的权重矩阵
Ra = R;
Qa =Ca'*Q*Ca;
Sa = Ca'*S*Ca;

R = Ra;
Q = Qa;
S = Sa;

% 定义系统运行步数
k_steps = 100;

% 定义x_history零矩阵，存储系统状态
x_history = zeros(size(x0, 1), k_steps);
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
    delta_u = -F*xa; %; % F常值；
    u = u + delta_u;
    % 计算系统响应
    x = AA*x +BB*u;
    xd = Ad*xd;
    xa = [x; xd; u];
    % 存储
    x_history(:, i+1) = x;
    u_history(:, i) = u;

    
    
end

figure(1)
subplot(2, 1, 1)
for i = 1:nA/2
    plot(x_history(i, :));
    hold on;
end
subplot(2, 1, 2)
for i = 1:mB
    plot(u_history(i, :));
    hold on;
end



















