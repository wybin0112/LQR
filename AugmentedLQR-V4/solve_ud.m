function ud = solve_ud(AA,BB,xd)
% 计算稳态控制输入ud,采用线性方程中方法，因为B可能不是方阵，无法求逆；

    ud = mldivide(BB, (eye(size(AA, 1))-AA)*xd);

end