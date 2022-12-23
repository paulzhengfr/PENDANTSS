%% simulation of SPOQ

label = 'final_v4_test';
nb_noise = 5 ;
p = 0.75 ; q = 2;
clearvars -except label nb_noise p q;

%% PENDANTS SPOQ
p = 0.75 ; q = 2;
penalty = 'SPOQ';
dataset = 'A'; noise_ratio = 0.005; dataId = [dataset, num2str(noise_ratio*100)]; 

main_plot_simulations;

disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])

clearvars -except label nb_noise p q;
penalty = 'SPOQ';
dataset = 'A'; noise_ratio = 0.01; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;

disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])

clearvars -except label nb_noise p q;
 
penalty = 'SPOQ';
dataset = 'A';
noise_ratio = 0.02; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])


clearvars -except label nb_noise p q;
penalty = 'SPOQ';
dataset = 'B';
noise_ratio = 0.005; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;
penalty = 'SPOQ';
dataset = 'B';
noise_ratio = 0.01; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;
 
penalty = 'SPOQ';
dataset = 'B';
noise_ratio = 0.02; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;
%% PENDANTS SOOT

p = 1; q = 2;
penalty = 'SPOQ';
dataset = 'A';
noise_ratio = 0.005; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])

clearvars -except label nb_noise p q;

penalty = 'SPOQ';
dataset = 'A';
noise_ratio = 0.01; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;

disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])

clearvars -except label nb_noise p q;
 
penalty = 'SPOQ';
dataset = 'A';
noise_ratio = 0.02; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])


clearvars -except label nb_noise p q;
 penalty = 'SPOQ';
dataset = 'B';
noise_ratio = 0.005; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;

penalty = 'SPOQ';
dataset = 'B';
noise_ratio = 0.01; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;
 
penalty = 'SPOQ';
dataset = 'B';
noise_ratio = 0.02; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;
%% backcorSOOT
% Simulation of backcorSOOT
clearvars -except label nb_noise p q;
p = 1; q = 2;
penalty = 'backcorSOOT';
dataset = 'A';
noise_ratio = 0.005; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations; 
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])


penalty = 'backcorSOOT';
dataset = 'A';
noise_ratio = 0.01; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])

clearvars -except label nb_noise p q;
 
penalty = 'backcorSOOT';
 
dataset = 'A';
noise_ratio = 0.02; dataId = [dataset, num2str(noise_ratio*100)]; 
 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;
 penalty = 'backcorSOOT';
 
dataset = 'B';
noise_ratio = 0.005; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations; 
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;

penalty = 'backcorSOOT';
 
dataset = 'B';
noise_ratio = 0.01; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations; 
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;
 
penalty = 'backcorSOOT';
 
dataset = 'B';
noise_ratio = 0.02; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;
%% backcorSPOQ
clearvars -except label nb_noise p q;
 p = 0.75; q = 2;
 penalty = 'backcorSPOQ';
dataset = 'A';
noise_ratio = 0.005; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])

penalty = 'backcorSPOQ';
dataset = 'A';
noise_ratio = 0.01; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])

clearvars -except label nb_noise p q;
 
penalty = 'backcorSPOQ';
 
dataset = 'A';
noise_ratio = 0.02; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations; 
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;
penalty = 'backcorSPOQ';
 
dataset = 'B';
noise_ratio = 0.005; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;

penalty = 'backcorSPOQ';
 
dataset = 'B';
noise_ratio = 0.01; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;
 
penalty = 'backcorSPOQ';
 
dataset = 'B';
noise_ratio = 0.02; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
%% Simulation of SCAD
clearvars -except label nb_noise p q;
penalty = 'SCAD';
dataset = 'A';
noise_ratio = 0.005; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])

clearvars -except label nb_noise p q;

penalty = 'SCAD';
dataset = 'A';
noise_ratio = 0.01; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])

clearvars -except label nb_noise p q;
 
penalty = 'SCAD';
dataset = 'A';
 noise_ratio = 0.02; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;

penalty = 'SCAD';
dataset = 'B';
noise_ratio = 0.005; dataId = [dataset, num2str(noise_ratio*100)]; 

main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;


penalty = 'SCAD';
dataset = 'B';
 noise_ratio = 0.01; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
clearvars -except label nb_noise p q;
 
penalty = 'SCAD';
dataset = 'B';
 noise_ratio = 0.02; dataId = [dataset, num2str(noise_ratio*100)]; 
main_plot_simulations;
disp(['Simulation of ', penalty,' with ',dataId , ' is accomplished'])
