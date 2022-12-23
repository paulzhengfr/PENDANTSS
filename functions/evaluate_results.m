function [res] = evaluate_results(s_res, t_res, p_res, converged, sbar,...
    tbar, pbar, y, thr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    P_res = convmtx(p_res,length(s_res)); 
    snr_res = 20 * log10(norm(sbar) / norm(sbar - s_res));
    res.snr = snr_res;
    
    s_sup = (sbar ~=0);
    tsnr_res = 20* log10(norm(sbar(s_sup)) / norm(sbar(s_sup) - s_res(s_sup)));
    res.tsnr = tsnr_res;
    
    diff_t = norm(t_res - tbar, 2)/norm(tbar,2);
    snr_baseline= -20 * log10(diff_t);
    res.snrT = snr_baseline;
    
    snr_pi = 20 * log10(norm(pbar) / norm(pbar - p_res));
    res.snrPi = snr_pi;
    
    sparsity_res = sum(abs(s_res)>thr);
    res.sparsity= sparsity_res;
    
    sparsity_orig = sum(s_sup);
    res.sparsityOrig = sparsity_orig;
    
    % test post process, but only works well for well-separated peaks
    % dataset
    s_res_3 = post_process(s_res, P_res, y - t_res, thr);
    snr_res3 = 20 * log10(norm(sbar) / norm(sbar - s_res_3));
    tsnr_res3 = 20* log10(norm(sbar(s_sup)) / norm(sbar(s_sup) - s_res_3(s_sup)));
    res.snr3 = snr_res3;
    res.tsnr3 = tsnr_res3;
    
    snr_reconstructed= 20 * log10(norm(y) / norm(y - P_res * s_res - t_res ));
    res.snr_reconstr = snr_reconstructed;
    
    indSpikeOrig = find(sbar);
    indSpikeEstim = find(abs(s_res)>thr);
    indDiffMn = indSpikeOrig - indSpikeEstim';
    [distNear, indNearest] = min(abs(indDiffMn),[],2);
    NbSpikeMissed = sum(distNear ~= 0);
    distMeanSpikeMissed = mean(distNear);

    sidePeak = zeros(10,1);
    for k = 1:10
    %     k= 5;
        indN = indNearest(k);
        indN_l = indN - 3 * (indN > 3);
        indN_u = indN + 3 * (indN < length(indSpikeEstim)-2); 
        sidePeak(k) = sum(abs(indDiffMn(k, indN_l:indN_u)) == 1) + sum(abs(indDiffMn(k, indN_l:indN_u)) == 2);
    end
    NbsidePeak = sum(sidePeak);

    NbPeakElseWhere = length(indSpikeEstim)-NbsidePeak - (sparsity_orig - NbSpikeMissed);
    
    res.NbsidePeak = NbsidePeak;
    res.NbSpikeMissed = NbSpikeMissed;
    res.NbPeakElseWhere = NbPeakElseWhere;
    res.converged = converged;
end

