%% Choose case
clear all
close all
dataset = 'B';
noise_ratio = 0.01;

% dataId = 'B1';
% 
penalty = 'SPOQ';
label = 'final_v4';


% penalty = 'SCAD';
% penalty ="backcorSOOT";
% label = '200testq_initsigma0_6'A
% label = 'initsigma1_ampli1_hNormalized';
% label ='test_postprocess';
p = 1; q = 2;
nb_noise = 1;
%% load signal


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
% load noise
fileID = fopen('data/noise200','r');
noiseMn_raw = fscanf(fileID,'%f');
fclose(fileID);
noiseMn = reshape(noiseMn_raw, [220,200]);

for id_noise = 1:nb_noise
    close all
% id_noise = 1;
w = sigma * noiseMn(:,id_noise);
y = xbar+ tbar + w;

%% Plot original figures
% figure(1)
% subplot(321)
% plot(1:N,sbar);
% title('s ');
% xlim([1 N]) 
% ylim([1.5*min(sbar) 1.5*max(sbar)]);
% 
% subplot(322)
% plot(1:M,xbar);
% title('xbar = Pi s');
% xlim([1 N]) 
% ylim([min(y) 1.2*max(y)]);
% 
% subplot(323)
% plot(t_var,pi);
% title('pi');
% %xlim([1 M]) 
% ylim([min(pi) max(pi)])
% 
% subplot(324)
% plot(1:M,w);
% title(['w: Bsnr ', num2str(BSNRinit)]);
% xlim([1 N]) 
% ylim([min(y) 1.2*max(y)]);
% subplot(325)
% plot(1:M,tbar);
% title('t tbar');
% xlim([1 N]) 
% ylim([min(y) 1.2*max(y)]);
% 
% subplot(326)
% plot(1:M,y);
% title('y = Pi s + w + t');
% xlim([1 N]) 
% ylim([min(y) 1.2*max(y)]);
% orient tall% tall
% print(['figure/resBD_orig_', dataId,'_n_', num2str(id_noise),'.pdf'], '-dpdf')

%% Load result


switch penalty
case 'SPOQ'
filename_ID = [label, '_',dataId,'_p_', num2str(p),'_q_', num2str(q), '_n_',num2str(id_noise)];
fileID = fopen(['result/BD_',label, '_',dataId,'_p_', num2str(p),'_q_', num2str(q), '_n_',num2str(id_noise),'_s3'],'r');
s_res_3 = fscanf(fileID,'%f');
fclose(fileID);
fileID = fopen(['result/BD_',label, '_',dataId,'_p_', num2str(p),'_q_', num2str(q), '_n_',num2str(id_noise),'_s'],'r');
s_res = fscanf(fileID,'%f');
fclose(fileID);
fileID = fopen(['result/BD_',label, '_',dataId,'_p_', num2str(p),'_q_', num2str(q), '_n_',num2str(id_noise),'_t'],'r');
t_res =fscanf(fileID,'%f');
fclose(fileID);
fileID = fopen(['result/BD_',label, '_',dataId,'_p_', num2str(p),'_q_', num2str(q), '_n_',num2str(id_noise),'_p'],'r');
p_res = fscanf(fileID,'%f');
fclose(fileID);

otherwise
filename_ID = [label,'_',dataId,'_', penalty, '_n_',num2str(id_noise)];
fileID = fopen(['result/BD_',label,'_',dataId,'_', penalty, '_n_',num2str(id_noise),'_s'],'r');
s_res = fscanf(fileID,'%f');
fclose(fileID);
fileID = fopen(['result/BD_',label,'_',dataId,'_', penalty, '_n_',num2str(id_noise),'_t'],'r');
t_res =fscanf(fileID,'%f');
fclose(fileID);
fileID = fopen(['result/BD_',label,'_',dataId,'_', penalty, '_n_',num2str(id_noise),'_p'],'r');
p_res = fscanf(fileID,'%f');
fclose(fileID);

end
P_res = convmtx(p_res,N); 
thr = sigma / 2^2; % / 2^6;
s_res_2 = remove_bias(s_res, P_res, y - t_res, thr);
%% Plot
% f = figure('visible','off');
% s = get(gca, 'Position');
% set(gca, 'Position', [s(1), s(2), s(3), s(4)*0.8])
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% figure
% title('sparse signals')
create_figure_sparse([sbar, s_res_2], dataset)
% stem(sbar,'--','Marker', 'o', 'MarkerFaceColor','none', 'MarkerEdgeColor','blue','MarkerSize', 6,'LineWidth',1.3)
% hold on
% stem(s_res_2,'--','Marker', 'x', 'MarkerFaceColor','none', 'MarkerEdgeColor','red','MarkerSize', 6,'LineWidth',1.3)
% xlim([0,200])
% ylim([0,30])
% stem(s0,'--', 'MarkerFaceColor','yellow', 'MarkerEdgeColor','none','MarkerSize', 3.5)
% stem(s_res_2,'--', 'MarkerFaceColor','yellow', 'MarkerEdgeColor','none','MarkerSize', 2.5)
% stem(s_res_3,'--', 'MarkerFaceColor','green', 'MarkerEdgeColor','none','MarkerSize', 2.5)
% legend('original signal', 'estimated signal')
% legend boxoff
% grid()
savefig(['figures2/',filename_ID,'_s.fig'])
saveas(gcf,['figures2/',filename_ID,'_s.pdf'])
markerIndex = [1, 6,10, 15:2:20, 24, 30, 37,39:2:60, 65, 70:2:80, 83, 92, 101, 105:3:145, 150, 158, 159:2:180, 190, 200:2:210, 215,220 ];
% [1:4:220]
% figure
Ymat = [tbar, t_res, tbar+xbar, t_res+ P_res*s_res];
create_figure_reconstruct(Ymat, dataset)
% plot(tbar, "-", 'LineWidth',2, 'Color','black')
% hold on
% plot(t_res,"--", 'LineWidth',2, 'Color','black')
% stem([2:8:220],t_res(2:8:220),"LineStyle", "none",'LineWidth',1,'Marker', 'x','MarkerFaceColor','none', 'MarkerEdgeColor','black','MarkerSize', 10)
% plot(tbar+ xbar,"-", 'LineWidth',1, 'Color','black')
% plot(t_res+ P_res*s_res,"--", 'LineWidth',1.5, 'Color','black')
% y_reconstructed = t_res+ P_res*s_res;
% stem(markerIndex,y_reconstructed(markerIndex),"LineStyle", "none", 'LineWidth',1, 'Marker', 'o','MarkerFaceColor','none', 'MarkerEdgeColor','black','MarkerSize', 6)
% plot(tbar, 'LineWidth',1,'Marker', 'o','MarkerFaceColor','none', 'MarkerEdgeColor','blue','MarkerSize', 6, 'MarkerIndices', [2:2:220])
% hold on
% plot(t_res, 'LineWidth',1,'Marker', 'x','MarkerFaceColor','none', 'MarkerEdgeColor','red','MarkerSize', 6, 'MarkerIndices', [1:2:220])
% plot(tbar+ xbar, 'LineWidth',1, 'color', 'blue','Marker', 'o','MarkerFaceColor','none', 'MarkerEdgeColor','blue','MarkerSize', 6, 'MarkerIndices', [2:2:220])
% plot(t_res+ P_res*s_res, 'LineWidth',1, 'color', 'red','Marker', 'x','MarkerFaceColor','none', 'MarkerEdgeColor','red','MarkerSize', 6, 'MarkerIndices', [1:2:220])
% if strcmp(dataset, 'A')
%     ylim([0,7])
% else
%     ylim([0,14])
% end
% 
% % title('baseline')
% % legend( 'orignal', 'estimated')
% % legend boxoff
% grid()
savefig(['figure_tmp/',filename_ID,'_t.fig'])
savefig(['figures2/',filename_ID,'_t.fig'])
saveas(gcf,['figures2/',filename_ID,'_t.pdf'])
% subplot(324)
% plot(P_res*s_res+t_res)
% hold on
% plot(P_res*s_res_2+t_res)
% plot(y)
% legend("reconstructed signal", "observed signal")
% legend boxoff
% grid()
% 
% subplot(325)
% plot(y - P_res* s_res - t_res)
% hold on
% plot(y - P_res* s_res_2 - t_res)
% grid()
% title('Residual')
% ylim([0, max(y)])
% 
% % set(gcf, 'PaperSize', [5 7])
% % s = get(gca, 'Position');
% % set(gca, 'Position', [s(1), s(2), s(3), s(4)*0.8])

create_figure_pi(t_var, [pi, p_res], dataset)
% figure
% plot(t_var,pi, 'LineWidth',1.5);
% hold on 
% plot(t_var, p_res,"--",'LineWidth',1.5);
% % title('convolution kernel');
% % legend('original kernel', 'estimated kernel')
% % legend boxoff
% xlim([-1,1])
% ylim([0,0.3])
% grid()
savefig(['figures2/',filename_ID,'_pi.fig'])
saveas(gcf,['figures2/',filename_ID,'_pi.pdf'])
% orient tall
switch penalty
    case 'SPOQ'
% sgtitle(['Blind deconvolution result with ', penalty,' $p =$ ', num2str(p), ' $q =$ ', num2str(q), ' noise ', num2str(id_noise)],'FontSize', 12)
% sgtitle(['Blind deconvolution result with ', penalty, ' noise ', num2str(id_noise)],'FontSize', 12)

% print(['diag_figure2/',label,'_', dataId,'_p_', num2str(p), '_q_', num2str(q),'_n_', num2str(id_noise),'.pdf'], '-dpdf')
    otherwise
% sgtitle(['Blind deconvolution result with ', penalty, ' noise ', num2str(id_noise)],'FontSize', 12)
% print(['diag_figure2/',label,'_', dataId,'_',penalty, '_n_', num2str(id_noise),'.pdf'], '-dpdf')

end
    
% close(f)
end
