function [x,p,finalsnr,SNR,MAE] = searchHP(p0,x0,I,algorithm,crit,iscvx,iter)


disp('---')
disp('SIMPLEX ALGORITHM');


[Nx,Ny] = size(I);
 

%p = p0;
psimp = [sqrt(p0(1));sqrt(p0(2))];

disp(['Estimation of parameters lambda and delta']);
[psimp,t]=Simplex('init',psimp);
 
 
P = zeros(iter,2);
 
SNR = zeros(iter,1);
MSSIM = SNR;
MAE = SNR;

for k=1:iter
    disp(['ITERATION ',num2str(k)])
    x = algorithm(x0,psimp);
    Irec = reshape(x,Nx,Ny);
    snr = 10*log10(sum(I(:).^2)/sum((I(:)-Irec(:)).^2));
    mae = sum(abs((I(:)-Irec(:))))/(Nx*Ny);
%    mssim = ssim(reshape(x,Nx,Ny), I);

    switch crit
        case 1
            f = -snr;
        case 2
            f = mae;
    end

    [psimp,t]=Simplex(f);
    
%     if(sum(psimp<=0)>0)
%         disp('ENTREE NEGATIVE!!')
%         [psimp,t]=Simplex('init',abs(psimp));
%     end
    
    if(iscvx)   %convex case => unique minimizer
    x0 = x;
    size(x0)
    end
    

    
    SNR(k) = snr;
    MAE(k) = mae;
    
%     figure(11)
%     imagesc(Irec);
%     colormap gray;
%     set(gca,'xtick',[ ])
%     set(gca,'ytick',[ ])
%     title(['k = ',num2str(k),' : SNR = ',num2str(snr),'; MAE = ',num2str(mae)])
%     xlabel(num2str((psimp(:).^2)'));
%     pause(1);
    
    disp(['RESULTS: ',num2str(k),'th iterate :']);
    for j=1:2
        disp(['p(',num2str(j),')^2 = ',num2str(psimp(j)^2)]);
    end
    disp(['SNR = ',num2str(snr),'; MAE = ',num2str(mae)]);
    disp('**********');
    
    %update
    P(k,:) = psimp';
%     p.lambda = psimp(1)^2;
%     p.delta  = psimp(2)^2;
end

psimp=Simplex('centroid'); % obtain the final value.
% p.lambda = psimp(1)^2;
% p.delta  = psimp(2)^2;
 
disp('Final result')
p = psimp;
x = algorithm(x0,psimp);
Irec = reshape(x,Nx,Ny);
finalsnr = 10*log10(sum(I(:).^2)/sum((I(:)-Irec(:)).^2));


% figure(2)
% plot(1:iter,SNR)

end

 