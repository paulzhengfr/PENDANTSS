M = 220;
rng(2)
noiseMn = randn(M,200);
i = 1;
noise_i = noiseMn(:,i);
fileID = fopen('data/noise200_2','w');
fprintf(fileID,'%9d\n',noiseMn);
fclose(fileID);

%% Read the noise file
fileID = fopen('data/noise','r');
noiseMn_read = fscanf(fileID,'%f');
fclose(fileID);
noiseMn_r = reshape(noiseMn_read, [220,10]);

plot(noiseMn_r(:,1))
