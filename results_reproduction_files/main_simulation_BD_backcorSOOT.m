%% data preparation
% num_sample_evaluate = 50;
% label = '200';
% 
% dataId = 'A1';


% penalty = 'backcorSOOT'; % or 'backcorSPOQ'
% p = 1; q = 2;
% p = 0.75;q =2;
% noise
nb_samples = 200;

fileID = fopen('data/noise200','r');
noiseMn_raw = fscanf(fileID,'%f');
fclose(fileID);
noiseMn = reshape(noiseMn_raw, [220,nb_samples]);
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
end
xmax = max(xbar);
sigma = noise_ratio * xmax;
load(['param/param_backcor_',dataId,'.mat'])
ord = 10; threshold = paramVn; fct = 'ah'; 
% Filter parameters
% t_var = linspace(-1,1,L)';
% s = 0.05;
% p0 = 1/sqrt(2*3.1415926) / s .*exp( -(t_var.^2)/(2*s^2) );
% p0 = p0/sum(p0); % other
% s0  = ones(M,1);
% [s0,p0,~,~,~,~] = initialization(sbar,pi);
% epsx = [min(sbar),max(sbar)];epsh = norm(pi);
% hmin = min(pi);
% hmax = max(pi);

fc = [];
d = 1;    % d : filter order parameter (d = 1 or 2)
theta = 0.5;
B = 50; % 50
iter_max = 2000;
gamma = 1.9 * ones(iter_max,1);%1.9
gamma = sparse(gamma);
iter_o = 1;
[s0, p0] = initialization(N,L);
% t_var = linspace(-1,1,L)';
% s = 0.05;
% p0 = 1/sqrt(2*3.1415926) / s .*exp( -(t_var.^2)/(2*s^2) );
% p0 = p0/sum(p0); % other
% s0  = ones(N,1);

alpha = 7e-7; %irrelevant.

% % load(['param/param_BD_SOOT_',dataId,'.mat'])
% % alpha = 7e-7; beta = paramVn(1); eta = paramVn(2);
% % lambda = paramVn(3);
% load(['param/param_p_',num2str(0.75), '_q_',num2str(q), '_',dataId,'.mat'])
% % load(['param/param_p_',num2str(p), '_q_',num2str(q), '_',dataId,'02273.mat'])
% alpha = 7e-7; beta = paramVn(1); eta = paramVn(2);
% load(['param/param_BD_p_',num2str(p), '_q_',num2str(q), '_',dataId,'_lambda.mat'])
% if length(paramVn) == 1
% lambda = paramVn(1);
% else
% beta = paramVn(1); eta = paramVn(2);lambda = paramVn(3);
% end
load(['param/param_backcorSPOQ_simplex2_p_',num2str(p), '_q_', num2str(q),'_',dataId,'.mat'])
beta = paramVn(1); eta = paramVn(2); lambda = paramVn(3);
% load(['param/param_backcorSPOQ_Its_p_',num2str(p), '_q_', num2str(q),'_',dataId,'.mat'])
best_Its = 1;
for id_noise = 1:num_sample_evaluate
    w = sigma * noiseMn(:,id_noise);
    y = xbar+ tbar + w;
    
    % run algorithm
    % backcor
    [t_res,coef,it] = backcor(1:M,y,ord,threshold,fct);
    y_woBaseline = y - t_res;
    fc = [];
    % SOOT
    [s_res, ~, p_res, cost, converged] = SPOQ_BD_quadraFidel(y_woBaseline, L, p,q,theta,B, iter_o, gamma, iter_max,...
    alpha, beta, eta, fc, d, lambda,s0,p0, best_Its);
    P_res = convmtx(p_res,N); 
    fileID = fopen(['result/BD_',label,'_',dataId,'_',penalty,'_n_',num2str(id_noise),'_s'],'w');
    fprintf(fileID,'%5d\n',s_res);
    fclose(fileID);
    fileID = fopen(['result/BD_',label,'_',dataId,'_',penalty,'_n_',num2str(id_noise),'_t'],'w');
    fprintf(fileID,'%5d\n',t_res);
    fclose(fileID);
    fileID = fopen(['result/BD_',label,'_',dataId,'_',penalty,'_n_',num2str(id_noise),'_p'],'w');
    fprintf(fileID,'%5d\n',p_res);
    fclose(fileID);

    
    % Calculate metrics
    
    thr = sigma / 2^2;
s_res_2 = remove_bias(s_res, P_res, y_woBaseline, thr);
    res = evaluate_results(s_res_2, t_res, p_res, converged, sbar,...
    tbar, pi, y, thr);
    writetable(struct2table(res), ['result/resBD_',label,'_',dataId,'_',penalty,'_n_',num2str(id_noise),'.txt'])
    fprintf('noise id %d ends  ',id_noise);
%     disp(['iteration ', num2str(id_noise)])
end
figure