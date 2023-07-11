function U = fun_2D_spatial_DFT(N)
% calculate the 2D spatial DFT to obtain the on-grid channel
% N is the total number of antennas
    N_axis = sqrt(N);
    U_axis = fun_spatial_DFT(N_axis);
    U = kron(U_axis,U_axis);
end
                     