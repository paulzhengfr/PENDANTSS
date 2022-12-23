
%% PENDANTS SPOQ
p = 0.75 ; q = 2;

for p = [0,75, 1]
for dataset = ['A', 'B']
for noise_ratio = [0.005, 0.01]
% dataset = 'A';  
dataId = [dataset, num2str(noise_ratio*100)]; penalty = 'SPOQ';
disp(dataId)
main_simplex_BD_SPOQ ;
% main_plot_simulations;

disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])

clearvars -except label num_sample_evaluate p q;
end
end
end

%% backcorSOOT
% Simulation of backcorSOOT
clearvars -except label num_sample_evaluate p q;
p = 1; q = 2;

for dataset = ['A', 'B']
for noise_ratio = [0.005, 0.01]
% dataset = 'A';
% noise_ratio = 0.005;
penalty = 'backcorSOOT';
dataId = [dataset, num2str(noise_ratio*100)]; 
main_simplex_BD_backcorSOOT ;
% main_plot_simulations; 
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label num_sample_evaluate p q;
end
end

%% backcorSPOQ
clearvars -except label num_sample_evaluate p q;
 p = 0.75; q = 2;
for dataset = ['A', 'B']
for noise_ratio = [0.005, 0.01]
% dataset = 'A';
% noise_ratio = 0.005;
penalty = 'backcorSPOQ';
dataId = [dataset, num2str(noise_ratio*100)]; 
main_simplex_BD_backcorSOOT ;
% main_plot_simulations; 
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label num_sample_evaluate p q;
end
end
