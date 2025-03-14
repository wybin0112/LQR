function [Fn_k] = LQR_solveF(A, B, Q, R, P)
    % 返回常值F矩阵
    Pk_1 = P;
    i = 1;
    % 定义系统稳态误差阈值
    tol = 1e-3;

    max_iter = 200;
    diff = inf;
    Fn_k_pre = inf;

    while(diff>tol)
        Fn_k = inv(B'*Pk_1*B+R)*B'*Pk_1*A;
        Pk_1 = (A-B*Fn_k)'*Pk_1*(A-B*Fn_k) + Fn_k'*R*Fn_k + Q;

        diff = abs(max(Fn_k-Fn_k_pre));
        i = i + 1;
        %if(i>max_iter)
         %   error('max_iter');
        %end
        Fn_k_pre = Fn_k;
    end
end
