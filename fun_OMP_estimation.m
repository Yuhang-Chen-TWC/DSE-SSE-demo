function [H_OMP, Pos_u1] = fun_OMP_estimation(y_vec, Sparse_level,A_observe,U_total,N_r, N_t)
% this function is used to calculate the OMP estimation result
% using different codebooks in A_observe
h_virtual = zeros(size(U_total,2),1);
r_n = y_vec;
Pos_u1 = zeros(1,Sparse_level); 
U1_save = zeros(size(A_observe,1), Sparse_level);
for idx_sparse = 1:Sparse_level
    product = A_observe' * r_n;
    [~,pos] = max(abs(product)); % find  the max column
    U1_save(:,idx_sparse) = A_observe(:,pos); % save the column
    Pos_u1(idx_sparse) = pos; % save the position of the column
    A_observe(:, pos) = zeros(size(A_observe,1),1);
    h_ls = (U1_save(:,1:idx_sparse)' * U1_save(:,1:idx_sparse))^(-1) *  U1_save(:,1:idx_sparse)' * y_vec;
    % estimation of the positions of the channel
    r_n = y_vec - U1_save(:,1:idx_sparse) * h_ls;
end
h_virtual(Pos_u1)= h_ls;
h_est = U_total * h_virtual;
H_OMP = reshape(h_est,N_r,N_t);
end