function [H_DSE] = fun_DSE_estimation(Y,Sparse_level_LC,U_suba_r,U_suba_t,W,F, N_r, N_t, K_r, K_t,Neighbor)
% the SSE CE algorithm
N_at = N_t/K_t;
N_ar = N_r/K_r;
A_r_suba = W' * U_suba_r; 
A_t_suba = F' * U_suba_t;

y_r = sum((Y * F'),2);
U_subasuba_r = U_suba_r(1:N_ar,1:N_ar);

for idx_kr = 1:K_r
    A_suba0 = A_r_suba(:, (idx_kr - 1) * N_ar + 1:idx_kr * N_ar);
    if idx_kr == 1
        A_suba = A_suba0;
%         figure()
%         plot(abs( A_suba' *  y_r))
        [~,Pos_initial_r] = fun_OMP_estimation(y_r, Sparse_level_LC ,A_suba, U_subasuba_r, N_ar, 1);
        Pos_DSE_r = Pos_initial_r;
        % calculate positions of others
        for idx_pos = 1:length(Pos_initial_r)
            Max_pos = Pos_initial_r(idx_pos); 
            for idx_neighbor = 1:Neighbor
                value = Max_pos - (Neighbor-1)/2 + idx_neighbor;
                if value > N_ar
                    value = value - N_ar + 1;
                else
                    if value <= 0
                        value = value + N_ar;
                    end
                end
                Pos_others_r((idx_pos - 1) * Neighbor + idx_neighbor) = value;
            end 
        end
        Pos_others_r = unique(Pos_others_r); % unique positions for others
    else
        A_suba = A_suba0(:,Pos_others_r);
        [~,Pos_res] = fun_OMP_estimation(y_r, Sparse_level_LC ,A_suba, U_subasuba_r, N_ar, 1);
        Pos_res_map = [];
        for idx = 1:length(Pos_res)
            val = (idx_kr - 1) * N_ar + Pos_others_r(Pos_res(idx));
            Pos_res_map(idx) = val;
        end
        Pos_DSE_r = [Pos_DSE_r, Pos_res_map];
    end
end



U_subasuba_t = U_suba_t(1:N_at,1:N_at);
y_t = (sum((W * Y ),1))';
for idx_kt = 1:K_t
    A_suba0 = A_t_suba(:, (idx_kt - 1) * N_at + 1:idx_kt * N_at);
    if idx_kt == 1
        A_suba = A_suba0;
%         figure()
%         plot(abs( A_suba' *  y_t))
        [~,Pos_initial_t] = fun_OMP_estimation(y_t, Sparse_level_LC ,A_suba, U_subasuba_t, N_at, 1);
        Pos_DSE_t = Pos_initial_t;
        % calculate positions of others
        for idx_pos = 1:length(Pos_initial_t)
            Max_pos = Pos_initial_t(idx_pos); 
            for idx_neighbor = 1:Neighbor
                value = Max_pos - (Neighbor-1)/2 + idx_neighbor;
                if value > N_at
                    value = value - N_at + 1;
                else
                    if value <= 0
                        value = value + N_at;
                    end
                end
                Pos_others_t((idx_pos - 1) * Neighbor + idx_neighbor) = value;
            end 
        end
        Pos_others_t = unique(Pos_others_t); % unique positions for others
    else
        A_suba = A_suba0(:,Pos_others_t);
        [~,Pos_res] = fun_OMP_estimation(y_t, Sparse_level_LC ,A_suba, U_subasuba_t, N_at, 1);
        Pos_res_map = [];
        for idx = 1:length(Pos_res)
            val = (idx_kt - 1) * N_at + Pos_others_t(Pos_res(idx));
            Pos_res_map(idx) = val;
        end
        Pos_DSE_t = [Pos_DSE_t, Pos_res_map];
    end
end


U_r = A_r_suba(:,Pos_DSE_r);
U_t = A_t_suba(:,Pos_DSE_t);
H_ls = (U_r' * U_r)^(-1) * U_r' * Y * ((U_t' * U_t)^(-1) * U_t')';
H_virtual = zeros(N_r, N_t);
for idx_r = 1:length(Pos_DSE_r)
    pos_r = Pos_DSE_r(idx_r);
    for idx_t = 1:length(Pos_DSE_t)
        pos_t = Pos_DSE_t(idx_t);
        H_virtual(pos_r,pos_t) = H_ls(idx_r,idx_t);
    end
end
H_DSE = U_suba_r * H_virtual * U_suba_t';
end 