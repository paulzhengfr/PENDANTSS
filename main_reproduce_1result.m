%% Load data
addpath(genpath(pwd));
% clear all
% close all
dataset = 'A';
p  = 0.75; q =2;
noise_ratio= 0.01;
dataId = [dataset, num2str(noise_ratio*100)];
% load noise
fileID = fopen('data/noise200','r');
noiseMn_raw = fscanf(fileID,'%f');
fclose(fileID);
% load data
switch dataId
    case {'A1', 'A0.5'}
        load('data/xA_1_realization_sig01.mat')
        fc = 0.04091;     % fc : low-pass filter cut-off frequency (cycles/sample)
    case 'A2'
        load('data/xA_1_realization_sig02.mat')
        fc = 0.05455;     % fc :low-pass filter  cut-off frequency (cycles/sample)
    case {'B1', 'B0.5'}
        load('data/xB_1_realization_sig01.mat')
        fc = 0.04091;     % fc :low-pass filter  cut-off frequency (cycles/sample)
    case 'B2'
        load('data/xB_1_realization_sig02.mat')
        fc = 0.0347;     % fc : low-pass filter cut-off frequency (cycles/sample)
%         fc = 0.02273;
end
w = w/ sigma;
xmax = max(xbar);
sigma = noise_ratio * xmax;
y = xbar+ tbar + w * sigma;

%% Load Algo parameters

d = 1;    % d : low-pass filter order parameter (d = 1 or 2)
theta = 0.5; % trust-region radius reduction factor.
B = 50; % Maximum trust-region iterations.
iter_max = 2000; % maximum iteration of the outer loop.
gamma = 1.9 * ones(iter_max,1);%1.9
gamma = sparse(gamma);
iter_o = 1; % kernel updates at each iteration.


%% Load Hyperparameters
hyperfound=1;
alpha=7e-7;
switch hyperfound
    case 1
        load(['param/param_PENDANTS_simplex2_p_',num2str(p), '_q_', num2str(q),'_',dataId,'.mat'])
beta = paramVn(1); eta = paramVn(2); lambda = paramVn(3);
    case 0
        main_run_simplex_BD_SPOQ; % parameter search with a known signal.
        load(['param/param_PENDANTS_simplex2_p_',num2str(p), '_q_', num2str(q),'_',dataId,'.mat'])
beta = paramVn(1); eta = paramVn(2); lambda = paramVn(3);
end

%% Run the algorithms and get the results
[s0, p0] = initialization(N,L);
[s_res, t_res, p_res,~, converged] = SPOQ_BD_quadraFidel(y, L, p,q,theta,B, iter_o, gamma, iter_max,... 
    alpha, beta, eta,  fc, d, lambda,s0, p0);
P_res = convmtx(p_res,N); 
thr = sigma / 2^2;
s_res_2 = remove_bias(s_res, P_res, y - t_res, thr);

%% Evaluate the results.
res = evaluate_results(s_res_2, t_res, p_res, converged, sbar,...
    tbar, pi, y, thr);
res