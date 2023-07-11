%% a demo of DSE and SSE channel estimation for UM-MIMO system
clear all
clc
N_t = 256;
N_r = 1024;
K_t = 1; % number of subarrays and RF chains Tx
K_r = 4; % number of subarrays and RF chains Rx
T_r = N_r/K_r * 3/4;
T_t = N_t * 3/4;
Aver = 10; 
L_path= 2;
N_at = N_t/K_t;
N_ar = N_r/K_r;
SNR_save = linspace(-20,20,6);
%% generate the HSPM channel
NMSE_SSE = zeros(Aver,length(SNR_save));
NMSE_DSE = zeros(Aver,length(SNR_save));
U_suba_r = fun_generate_u_total_suba(K_r,1,N_r,1); 
U_suba_t = fun_generate_u_total_suba(1,K_t,1,N_t);
NN = [N_r, N_t];
for aver = 1:Aver
    [H_HSPM] = fun_generate_HSPM(N_r, N_t, 1, 1, L_path);
    W = ((rand(N_r, N_r * 3/4) > 0.5) * 2 - 1)/sqrt(N_r);
    F = ((rand(N_t, T_t) > 0.5) * 2 - 1)/sqrt(N_t); 
    A_observe_r = W' * U_suba_r;
    A_observe_t = F' * U_suba_t;
    Neighbor = 5;
    %% Try to use the subarray by sybarray method for channel estimation
    for idx_snr = 1:length(SNR_save)
        SNR = SNR_save(idx_snr);
        Y = W' * awgn(H_HSPM * F,SNR,'measured');
        y_vec = vec(Y);
        Sparse_level_total = L_path * K_t * K_r * 2;
        Sparse_level_LC = L_path * 2;
        [H_SSE] = fun_2stage_IRS_CE_OMP(Y,Sparse_level_total, A_observe_r,A_observe_t,U_suba_r,U_suba_t,NN);

        [H_DSE] = fun_DSE_estimation(Y,Sparse_level_LC,U_suba_r,U_suba_t,W,F, N_r, N_t, K_r, K_t,Neighbor);
        NMSE_SSE(aver, idx_snr) = 10 * log10((norm(H_HSPM - H_SSE, 'f')/norm(H_HSPM, 'f'))^2);
        NMSE_DSE(aver, idx_snr) = 10 * log10((norm(H_HSPM - H_DSE, 'f')/norm(H_HSPM, 'f'))^2);
    end
end 

figure()
plot(SNR_save,mean(NMSE_SSE),'r-^', 'linewidth',1)
hold on
plot(SNR_save,mean(NMSE_DSE),'b-*', 'linewidth',1)
hold on
set(gca,'FontSize',16);
xlabel('SNR (dB)','FontSize',16,'FontName','Times New Roman')
ylabel('Estimation NMSE (dB)','FontSize',16,'FontName','Times New Roman')
hhandle=legend('SSE','DSE');
set(hhandle,'FontSize',16,'FontName','Times New Roman');
grid on;
