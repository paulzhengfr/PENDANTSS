% dataId= 'A2';
% penalty = "SPOQ";
% 
% % label = 'initsigma1_ampli1_hNormalized';
% label ='test_postprocess';
% num_sample_evaluate = 50;

% p = 0.75;q = 2;
tab = [];



%%

for id_noise = 1:num_sample_evaluate % 200 %10
    switch penalty
        case "SPOQ"
            tab = [tab;readtable(['result/resBD_',label,'_',dataId,'_p_', num2str(p),'_q_', num2str(q), '_n_',num2str(id_noise),'.txt'])];
%     delete(['result/resBD_',dataId,'_p_', num2str(p),'_q_', num2str(q), '_n_',num2str(id_noise),'.txt'])
        otherwise
            tab = [tab;readtable(['result/resBD_',label,'_',dataId,'_', penalty,'_n_',num2str(id_noise),'.txt'])];
%     delete(['result/resBD_',dataId,'_backcorSOOT_n_',num2str(id_noise),'.txt'])
    end
end
switch penalty
    case "SPOQ"
        writetable(tab, ['result/resBD_',label,'_',dataId,'_p_', num2str(p),'_q_', num2str(q),'.txt'])
    otherwise
        writetable(tab, ['result/resBD_',label,'_',dataId,'_', penalty,'.txt'])
end

%% 
% dataId= 'A1';
% p = 0.75;q = 2;
% tab = readtable(['result/resBD_',dataId,'_p_', num2str(p),'_q_', num2str(q),'.txt']);
% tab = readtable(['result/resBD_',dataId,'_SCAD.txt']);
tabA = table2array(tab);
stat.mean = mean(tabA);
stat.std = std(tabA);
stat.min = min(tabA);
stat.max = max(tabA);
stat.median = median(tabA);
stat.mad = mad(tabA,1);

statArray = [stat.mean;stat.std ;stat.min;stat.max;stat.median;stat.mad];
statTab = array2table(statArray);
statTab.Properties.VariableNames = tab.Properties.VariableNames;
statMetrics = fieldnames(stat);
statTab.Properties.RowNames = statMetrics(1:6);

switch penalty
    case "SPOQ"
        writetable(statTab, ['result/stat/stat_resBD_',label,'_',dataId,'_p_', num2str(p),'_q_', num2str(q),'.xlsx'])
    otherwise
        writetable(statTab, ['result/stat/stat_resBD_',label,'_',dataId,'_', penalty, '.xlsx'])
end
