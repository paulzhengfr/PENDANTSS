function [x_out] = remove_bias(x_in,K,y,seuil)
%remove_bias remove bias of the result, part of the post processing.
%   
%   x_in: signal to be bias removed.
%   K: kernel vec
%   y: observation signal (without trend)
%   seuil: threshold.
%

index_unb = find(abs(x_in)>=seuil);
K_unb = K(:,index_unb);
x_out = x_in;
x_out(index_unb) = pinv(K_unb'*K_unb)*K_unb'*y;
end

