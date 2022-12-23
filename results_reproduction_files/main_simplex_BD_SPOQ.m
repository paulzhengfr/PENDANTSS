%% One realisation of signal
% clear all
% close all
% dataId = 'B1';
% p  = 0.75; q =2;
dataId = [dataset, num2str(noise_ratio*100)];
switch dataId
    case {'A1', 'A0.5'}
        load('data/xA_1_realization_sig01.mat')
        fc = 0.04091;     % fc : cut-off frequency (cycles/sample)
    case 'A2'
        load('data/xA_1_realization_sig02.mat')
        fc = 0.05455;     % fc : cut-off frequency (cycles/sample)
    case {'B1', 'B0.5'}
        load('data/xB_1_realization_sig01.mat')
        fc = 0.04091;     % fc : cut-off frequency (cycles/sample)
    case 'B2'
        load('data/xB_1_realization_sig02.mat')
        fc = 0.0347;     % fc : cut-off frequency (cycles/sample)
%         fc = 0.02273;
end
w = w/ sigma;
sigma = noise_ratio * xmax;
y = xbar+ tbar + w * sigma;
%% Fix Parameters
% Filter parameters
% fc = 0.0347;     % fc : cut-off frequency (cycles/sample) 0.038 0.14
d = 1;    % d : filter order parameter (d = 1 or 2)
theta = 0.5;
B = 50; % 50
iter_max = 2000;
gamma = 1.9 * ones(iter_max,1);%1.9
gamma = sparse(gamma);
% xsi = 1.0*sigma * sqrt(M);  
% itermax_PPXA = 2000;
% prec_PPXA = 1e-7;
iter_o = 1;
[s0, p0] = initialization(N,L);
%% Parameters SPOQ

alpha = 7e-7; %irrelevant.
% load(['param/param_p_', num2str(p), '_q_', num2str(q),'_',dataId,'.mat'])
% beta_i = paramVn(1); eta_i = paramVn(2);
% load(['param/param_Pendants_bestgrid_p_', num2str(p), '_q_', num2str(q),'_',dataId,'.mat'])
load(['param/param_Pendants_simplex_p_', num2str(p), '_q_', num2str(q),'_',dataId,'.mat'])
beta_i = paramVn(1); eta_i = paramVn(2); lambda_i = paramVn(3);
%% Simplex search
% initialization 
% load(['param/param_BD_p_', num2str(p), '_q_', num2str(q),'_',dataId,'_lambda.mat'])
% lambda = 1;
% lambda = 1;
% paramVn = [beta_i;eta_i; lambda];
% paramVn = [lambda];
paramVn = paramVn;

[paramVn,t]=Simplex('init',paramVn, paramVn * 0.5);
iters=0;
for i=1:200
    beta_i =paramVn(1);
    eta_i = paramVn(2);
    lambda_i = paramVn(3); %paramVn
%     lambda_i = paramVn; %paramVn
    if ~test_parameters(p, q, alpha, beta_i, eta_i, B, theta)
        val2mini = 1;
        [paramVn,t]=Simplex(val2mini);
        disp(['iteration ', num2str(i), ' val2mini ', num2str(val2mini)])
        ok=Simplex('converged', 1e-5); 
    else
    [s_res, t_res, p_res,~, converged] = SPOQ_BD_quadraFidel(y, L, p,q,theta,B, iter_o, gamma, iter_max,...
    alpha, beta_i, eta_i,  fc, d, lambda_i,s0, p0);
%     beta_i = paramVn(1); eta_i = paramVn(2); lambda_i = paramVn(3);

%     [s_res, t_res, p_res, cost] = SPOQ_BD_quadraFidel(y, L, p,q,theta,B, iter_o, gamma, iter_max,...
%     alpha, beta_i, eta_i, xsi,itermax_PPXA, prec_PPXA, fc, d,lambda_i);
    P_res = convmtx(p_res,N); 

    snr_res = 20 * log10(norm(sbar) / norm(sbar - s_res));
    thr= sigma / 4;
    s_res_2 = remove_bias(s_res, P_res, y - t_res, thr);
    snr_res2 = 20 * log10(norm(sbar) / norm(sbar - s_res_2));
    s_sup = (sbar ~=0);
    tsnr_res = 20* log10(norm(sbar(s_sup)) / norm(sbar(s_sup) - s_res(s_sup)));
    tsnr_res2 = 20* log10(norm(sbar(s_sup)) / norm(sbar(s_sup) - s_res_2(s_sup)));
    snr_reconstructed = 20 * log10(norm(y)  / norm(y - P_res* s_res_2 - t_res));
    snr_pi = 20 * log10(norm(pi) / norm(pi - p_res));
    diff_t = norm(t_res - tbar, 2)/norm(tbar,2);
    snr_baseline= -20 * log10(diff_t);
    val2mini = - 2* snr_res2 - snr_baseline- snr_pi; %- snr_pi - snr_baseline;
%     end
    disp(['iteration ', num2str(i), ' val2mini ', num2str(val2mini), ' snr ', num2str(snr_res2), ' tsnr ', num2str(tsnr_res2), ' snr_pi ', num2str(snr_pi), ' snr_t ', num2str(snr_baseline)])
    [paramVn,t]=Simplex(val2mini);
%     if ~mod(i, 5)  % every 5 iterations print out the fitted value
%         disp(['beta=', num2str(paramVn(1)),' eta=',num2str(paramVn(2))])
%     end
    ok=Simplex('converged', 1e-5); 
    if ok
        break
    end
    end
end

disp([dataId, ' done'])
disp(paramVn)
% paramVn=Simplex('centroid'); 
% paramVn
if ok
    paramVn=Simplex('centroid'); % obtain the final value if ok
    disp(['valeur ', num2str(val2mini)])
%     save(['param/param_BD_testp_p_',num2str(p), '_q_', num2str(q),'_',dataId,'_lambda.mat'],'paramVn')
    save(['param/param_PENDANTS_simplex2_p_',num2str(p), '_q_', num2str(q),'_',dataId,'.mat'],'paramVn')
else
    disp([dataId,' not ok'])
    paramVn=Simplex('centroid'); % obtain the final value if ok
%     save(['param/param_BD_testp_p_',num2str(p), '_q_', num2str(q),'_',dataId,'_lambda.mat'],'paramVn')
    save(['param/param_PENDANTS_simplex2_p_',num2str(p), '_q_', num2str(q),'_',dataId,'.mat'],'paramVn')
    disp(['valeur ', num2str(val2mini)])
end