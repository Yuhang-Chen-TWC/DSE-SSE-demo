
function [H_HSPM] = fun_generate_HSPM(N_r, N_t, K_r, K_t, L_path)
%% Generate the Tx Rx channel matrix.
% We consider a random HSPM channel
N_atx = sqrt(N_t/K_t);
N_arx = sqrt(N_r/K_r);
N_at = N_t/K_t;
N_ar = N_r/K_r;
d_atx = 0:N_atx-1;
d_arx = 0:N_arx-1;
H_HSPM = zeros(N_r,N_t);
for idx_ltm = 1:L_path
    if idx_ltm == 1
        alpha(idx_ltm) = 1/sqrt(1) * (rand(1) + 1j * rand(1));
    else
        alpha(idx_ltm) = 10^(-0.7) * 1/sqrt(1) * (rand(1) + 1j * rand(1));
        % assume the NLoS path is 7dB lower than the LoS path
    end 
    for idx_kr = 1:K_r
        for idx_kt = 1:K_t
            if idx_kr == 1 && idx_kt == 1
                theta_r(idx_kr,idx_kt,idx_ltm) = rand(1) * pi;
                phi_r(idx_kr,idx_kt,idx_ltm) =  rand(1) * pi;
                theta_t(idx_kr,idx_kt,idx_ltm) = rand(1) * pi;
                phi_t(idx_kr,idx_kt,idx_ltm) = rand(1) * pi;
            else
                theta_r(idx_kr,idx_kt,idx_ltm) = theta_r(1,1,idx_ltm) + 1/N_arx * rand(1) * pi;
                phi_r(idx_kr,idx_kt,idx_ltm) =  phi_r(1,1,idx_ltm) + 1/N_arx * rand(1) * pi;
                theta_t(idx_kr,idx_kt,idx_ltm) = theta_t(1,1,idx_ltm) + 1/N_atx * rand(1) * pi;
                phi_t(idx_kr,idx_kt,idx_ltm) = phi_t(1,1,idx_ltm) + 1/N_atx * rand(1) * pi;
            end
            beta_rx(idx_kr,idx_kt,idx_ltm) = sin(theta_r(idx_kr,idx_kt,idx_ltm)) * cos(phi_r(idx_kr,idx_kt,idx_ltm));
            beta_rz(idx_kr,idx_kt,idx_ltm) = sin(phi_r(idx_kr,idx_kt,idx_ltm));
            beta_tx(idx_kr,idx_kt,idx_ltm) = sin(theta_t(idx_kr,idx_kt,idx_ltm)) * cos(phi_t(idx_kr,idx_kt,idx_ltm));
            beta_tz(idx_kr,idx_kt,idx_ltm) = sin(phi_t(idx_kr,idx_kt,idx_ltm)); 
            
            a_rx = exp(1j * pi * beta_rx(idx_kr,idx_kt,idx_ltm) * d_arx); 
            a_rz = exp(1j * pi * beta_rz(idx_kr,idx_kt,idx_ltm) * d_arx); 
            a_tx = exp(1j * pi * beta_tx(idx_kr,idx_kt,idx_ltm) * d_atx); 
            a_tz = exp(1j * pi * beta_tz(idx_kr,idx_kt,idx_ltm) * d_atx); 
            a_r = kron(a_rx,a_rz);
            a_t = kron(a_tx,a_tz);
            H_HSPM(N_ar*(idx_kr-1)+1:N_ar*idx_kr, N_at*(idx_kt-1)+1:N_at*idx_kt) = H_HSPM(N_ar*(idx_kr-1)+1:N_ar*idx_kr, N_at*(idx_kt-1)+1:N_at*idx_kt) ...
                + alpha(idx_ltm) * a_r' * a_t;
        end
    end
end

end