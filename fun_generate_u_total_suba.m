function [U_total] = fun_generate_u_total_suba(K_r,K_t,N_r, N_t)
% This function is used to generate the subarray based codebook
% size of U_total: N_rN_t * N_r N_t
N_ar = N_r/K_r;
N_at = N_t/K_t;
U_ar = fun_2D_spatial_DFT(N_ar);
U_at = fun_2D_spatial_DFT(N_at);
U_r = [];
U_t = [];
for idx_kr = 1:K_r
    U_r = blkdiag(U_r, U_ar);
end
for idx_kt = 1:K_t
    U_t = blkdiag(U_t, U_at);
end
U_total = kron(U_t, U_r);
end




