function [x0,h0] = initialization(N,L)
% INITIALIZATION    Initialize the sparse signal and the kernel.
% 
%   N: length of sparse signal
%   L: length of kernel to be estimated.
%

sigma = 1 ; % 1 to 0.5
c = floor(L/2)+1 ; 
gg = 1:L ;
gauss_c = exp(-(gg - c).^2/(2*sigma^2));
h0 = (gauss_c)';
h0 = h0 / sum(h0);


x0 = ones(N,1);% 

end
