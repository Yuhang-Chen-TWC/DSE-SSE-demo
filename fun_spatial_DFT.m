function U = fun_spatial_DFT(N)
 deta = 1/N;
 for ii= -(N-1)/2:1:(N-1)/2
      U(:,ii+(N+1)/2) = sqrt(1/N)*exp(1i*[0:N-1]*2*pi*deta*ii).'; % spatial DFT matrix n*n
 end
end
