function [H_OMP,Pos_t] = fun_2stage_IRS_CE_OMP(Y,S, A_r,A_t,U_dr,U_dt,NN)
% perform the two stage channel estimation for IRS systems using OMP
% algorithm 
N_r = NN(1);
N_t = NN(2);
% Rx channel
y_r = sum((Y * A_t),2);
[~,Pos_r] = fun_OMP_estimation(y_r,S,A_r, U_dr, N_r, 1); 
% Tx channel
y_t = (sum( (A_r(:,Pos_r)' * Y),1 ))';
[~,Pos_t] = fun_OMP_estimation(y_t, S ,A_t, U_dt, N_t,1); 

U_r = A_r(:,Pos_r);
U_t = A_t(:,Pos_t);
H_ls = (U_r' * U_r)^(-1) * U_r' * Y * ((U_t' * U_t)^(-1) * U_t')';
H_virtual = zeros(N_r, N_t);
for idx_r = 1:length(Pos_r)
    pos_r = Pos_r(idx_r);
    for idx_t = 1:length(Pos_t)
        pos_t = Pos_t(idx_t);
        H_virtual(pos_r,pos_t) = H_ls(idx_r,idx_t);
    end
end
H_OMP = U_dr * H_virtual * U_dt';
end
