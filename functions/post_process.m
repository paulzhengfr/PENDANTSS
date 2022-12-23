function [s_out] = post_process(s, P, y, thr)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
index_unb = find(abs(s)>= thr);

index_double_peak_center = [];
bool = 0;
for i = 1: length(index_unb)-2
    if index_unb(i+1) - index_unb(i) == 2 %|| index_unb(i+1) - index_unb(i) == 4
        index_double_peak_center = [index_double_peak_center, floor((index_unb(i+1) + index_unb(i))/2)];
%         if index_unb(i+1) - index_unb(i) == 4
%             bool = 1;
%         end
    elseif index_unb(i+1) - index_unb(i) == 1 && index_unb(i+2) - index_unb(i+1) == 1
        index_double_peak_center = [index_double_peak_center, index_unb(i+1)];
%     elseif index_unb(i+1) - index_unb(i) == 2 && index_unb(i+2) - index_unb(i+1) == 2
%         index_double_peak_center = [index_double_peak_center, index_unb(i+1)];
%         bool = 1;
    end
end
disp(index_double_peak_center)
s_out = s;
if ~isempty(index_double_peak_center)
for ind = index_double_peak_center
    s_out(ind) = s_out(ind) + s_out(ind - 1) + s_out(ind + 1);
    
    s_out(ind - 1) = 0;
    s_out(ind + 1) = 0;
    if bool
    s_out(ind - 2) = 0;
    s_out(ind + 2) = 0;
    end
end
index_unb = find(abs(s_out)>= thr);
K_unb = P(:,index_unb);
s_out(index_unb) = pinv(K_unb'*K_unb)*K_unb'*y;
end
end

