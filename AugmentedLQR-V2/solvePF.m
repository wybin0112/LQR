function [Pk, F] = solvePF(A, B, Q, R, P, k_steps)
    Pk_1 = P;
    for i = 1:1:k_steps
        Fn_k = inv(B'*Pk_1*B+R)*B'*Pk_1*A;
        Pk_1 = (A-B*Fn_k)'*Pk_1*(A-B*Fn_k) + Fn_k'*R*Fn_k + Q;
        if(i == 1)
            F = Fn_k
            Pk = Pk_1;
        end
        F = [Fn_k; F];
        Pk = [Pk; Pk_1];
    end
end

