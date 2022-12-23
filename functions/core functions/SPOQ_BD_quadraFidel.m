function [s_res, t,p_res, cost, converged] = SPOQ_BD_quadraFidel(y, L, p,q,theta,I_iter, J,gamma, iter_max,...
    alpha,beta, eta, fc, d_order,lambda, s0, p0, nb_s_innerloop)
%SPOQ_BD_quadraFidel PENDANTSS main algorithm for joint blind deconvolution and
%detrending.
%
%   y: vector of length N of observed data 
%   L: kernel length to be estimated
%   p, q: parameters of SPOQ, lp / lq
%   theta: factor of trust region size reduction
%   I_iter: maximum number of trust region trial
%   J: number of kernel updates per iteration.
%   gamma: vector of length euqal to iter_max of sequence of positive
%   stepsize between ]0,2[
%   iter_max: maximum of iterations
%   alpha, eta, beta: parameters of SPOQ penalty.
%   fc: cutoff frequency for trend removal.
%   d_order: order of filter for trend removal(usually taken = 1).
%   lambda: weight hyperparameter for SPOQ penalty
%   s0, p0: initialized s and p.
%   nb_s_innerloop: number of iterations on s. usually at = 1.
%return
%   s_res: sparse signal result of PENDANTSS.
%   t: baseline estimated.
%   p_res: kernel estimated.
%   cost: not used here
%   converged: bool for convergence.

if nargin == 17
    nb_s_innerloop = 1;
end

%% initialization
converged = false;
M = length(y);
N = M - L+1;
s_res = s0;
p_res = p0;

P = convmtx(p_res, N); % init convolution matrix of pi
cost = zeros(iter_max,4);
% test if parameters satisfy the convergent conditions
if ~test_parameters(p, q, alpha, beta, eta, I_iter, theta)
    s_res = zeros(N,1);
    t = zeros(M,1);
    warning("hyperparameters error")
    return
end
    

% tic;
% calculate the high pass filter H
[A_filt, B_filt] = BAfilt(d_order, fc, M);
H_mat =  B_filt* inv(A_filt);
H_func = @(x) B_filt*(A_filt\x);
Hy = H_func(y);
%HtH = transpose(H_mat) * H_mat;            

HP = H_mat*P;
HPt = HP';
HPtHP = HP' * HP; % (HP)^T HP
L1 = norm(HPtHP); % 


% SPOQ regularization function
func_lpp = @(x) sum((x.^2 + alpha^2).^(p/2) - alpha^p);
% func_lq = @(x) (eta^q + sum(abs(x).^q)).^(1/q);
% Psi = @(x) 1/p * log(func_lpp(x) + beta^p) -  log(func_lq(x));
delta_lpp= @(x) p * x .* (x.^2 + alpha^2).^(p/2-1);
delta_lqq = @(x) q * sign(x) .* abs(x).^(q-1);
delta_Psi1 = @(x) 1/p * delta_lpp(x) ./ (func_lpp(x) + beta^p);
delta_Psi2 = @(x) 1 / q * delta_lqq(x) ./ (eta^q + sum(abs(x).^q));
delta_Psi = @(x) delta_Psi1(x) - delta_Psi2(x);


% %initialization of loop
rho = 1; % take a value, no importance, init later in the function

% cost(1,2) = (norm(HP*s - Hy,2)- xi) * (norm(HP*s - Hy,2)>xi) +...
%     norm(s.*(s<=0),2);
% cost(1,3)  =  Psi(s);
% cost(1, 1) = cost(1,2) + cost(1,3);
% cost(1,4) = norm(1 - sum(p_res));
% disp(['initialization takes this amount of time'])
% toc;
%% SPOQ BD loop

% tStart = tic;  
for k = 1:iter_max
%     tic;
    gamma_k = gamma(k);
    L1mat = L1 * eye(N);
    % trust region loop
    s_beforeloop = s_res;
    for ii = 1:nb_s_innerloop
    for i = 1:I_iter
        rho =  calc_TRradius(s_res, i,q,I_iter, rho, theta); % trust region radius
        if length(rho) > 1
            warning("length of rho superior than 1")
        end
        % Ak: MM metrics
        A_q_rho = calc_A_rho_p(N,s_res,rho,p, q, alpha,beta,eta, func_lpp);
        Ak = lambda * A_q_rho + L1mat; % forgot lambda
        
        % Find s_k,i
        delta_rho1 = HPtHP * s_res - HPt * Hy; 
        delta_f1 = lambda * delta_Psi(s_res) + delta_rho1;
        x_f = s_res - gamma_k  * (Ak \ delta_f1);  % forward step
        % calculate the proximity operator value: projection to positive
        % values
        s_ki  = prox_g1(x_f);
        if sum(abs(s_ki).^q) >= rho^q
%             disp(['trust region loop breaks at i = ', num2str(i)])
            break
        end
        if i == I_iter
            disp('trust region loop arrives at the maximum iterations')
        end
    end
    
    s_res = s_ki;
    end
    Norm_Sdiff = norm(s_res-s_beforeloop); % norm difference of s update
%     if rem(k,500) == 0
%         disp(['iteration k = ', num2str(k)])
%     end
%     if norm(zk - x,2) < 1
%         x = zk;
%         disp(['general loops break'])
%         break
%     end

    % Testing inner loop for trust region approach.
%     Norm_Sdiff = norm(s_res-s_ki); % norm difference of s update
%     s_res = s_ki;
    
    SMn = convmtx(s_res,L); 
    HS = H_mat*SMn;
    HSt = HS';
    HStHS = HS' * HS;
    L2 = norm(HStHS,2);

    % Kernel updates
    pVn = p_res;
    
    for j = 1:J
        deltaf2 = HStHS * pVn - HSt * Hy;
        ytild = pVn - gamma_k / L2 * deltaf2;
        pVn = prox_g2(ytild, L2);
    end
%     Norm_Pidiff = norm(p_res - pVn);
    p_res = pVn;
    P = convmtx(p_res,N); 
    HP = H_mat*P;
    HPt = HP';
    HPtHP = HP' * HP;
    L1 = norm(HPtHP); % + 0.1    
    
%     cost(k+1,2) = (norm(HP*s - Hy,2)- xi) * (norm(HP*s - Hy,2)>xi) +...
%         norm(s.*(s<=0),2);
%     cost(k+1,4) = norm(0 - sum(p_res));
%     cost(k+1,3)  =  Psi(s);
%     cost(k+1, 1) = cost(k+1,2) + cost(k+1,3)+cost(k+1,4);
    
    if(k>2)
%         CondS_norm = Norm_Xold/sqrt(M) + Norm_Piold/ sqrt(L);
        CondS_norm = Norm_Sdiff/sqrt(M);
        if(CondS_norm < 1e-6)
            converged = true;
            break
        elseif k == iter_max
            warning("Algorithm not converged")
            disp("Algorithm not converged")
        end
    end
%     T(k)= toc;
end

% tMul = sum(T);
% tEnd = toc(tStart);
% disp(['The loop takes ', num2str(tEnd), ' seconds'])
% disp(['Each iteration takes ', num2str(T)])
s_res = s_res;
vec = y - P*s_res;
t = vec - H_func(vec);  % Compute baseline estimated
[s_res, p_res] = correct_shift(s_res,p_res);
if nargout==4
    converged = {};
end
end

%% used functions
function [res] = calc_TRradius(x, i,q,B, rho_prev, theta)
%calc_TRradius calculate the trust region radius
    if i == 1
        res = sum(abs(x).^q);
    elseif i >=2 && i <= B - 1
        res = theta * rho_prev;
    elseif i == B
        res = 0;
    else 
        warning("index when calculating TR radius is erroneous")
    end
end

function [A] = calc_A_rho_p(N,x,rho,p, q, alpha,beta,eta, func_lpp)
%calc_A_rho_p calculate the SPD metric of this step
    qi = (q - 1) / (eta^q + rho^q)^(2/q);
    iden = qi * eye(N);
    diag_term = 1 / (func_lpp(x)+beta^p) * diag((x.^2+alpha^2).^(p/2 -1));
    A = iden + diag_term;
%     A= sparse(A);
end

% --- local function ----
% 
function [A, B] = BAfilt(d, fc, N)
% [A, B] = BAfilt(d, fc, N)
%
% Banded matrices for zero-phase high-pass filter.
% The matrices are 'sparse' data type in MATLAB.
%
% INPUT
%   d  : degree of filter is 2d (use d = 1 or 2)
%   fc : cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   N  : length of signal
if d == 0 || isempty(fc)
    A = eye(N);
    B = A;
else
b1 = [1 -1];
for i = 1:d-1
    b1 = conv(b1, [-1 2 -1]);
end
b = conv(b1, [-1 1]);

omc = 2*pi*fc;
t = ((1-cos(omc))/(1+cos(omc)))^d;

a = 1;
for i = 1:d
    a = conv(a,[1 2 1]);
end
a = b + t*a;

A = spdiags( a(ones(N, 1), :), -d:d, N, N);   % A: Symmetric banded matrix
B = spdiags(b(ones(N, 1), :), -d:d, N, N);    % B: banded matrix

end
end


%% 
function [bool] = test_parameters(p, q, alpha, beta, eta, B, theta)
bool = true;
% check the domain of definition
if p<=0 || p >= 2|| q < 2
    bool = false;
    disp(['false domain of definition of (p,q)'])
elseif alpha <= 0 || beta <= 0 || eta <= 0
    bool = false;
    disp(['false domain of definition of (alpha, beta, eta)'])
elseif q <= 2 && (q ~= 2 || eta^2 * alpha^(p-2) <= beta^p)
    bool = false;
    disp(['0 may not be a local minimizer of \Psi'])
elseif eta^2 < beta^2 * max(8 * alpha^(2-p) / p / (2+p) / beta^(2-p),1 / (2^(p/2) - 1)^(2/p))
    bool = false;
    disp(['0 may not be a global minimizer of \Psi'])
end
if B < 1 || theta <= 0 || theta >=1 
    bool = false;
    disp(['The value of B or theta is erroneous'])
end
end


%% calculate the prox wrt o in Blind deconvolution 
function[zk] = prox_g1(x)
%proxPPXAplus calculate the proximity operator's value by PPXA+ 
% D = HX
% x = o_res (yVn)
zk = x;
zk(zk < 0) = 0;
end

function[zk] = prox_g2(x,L2)
vecx = x;
if size(x,1) > 1
    vecx = transpose(x);
end
zk = projection_metric_simplex(vecx, ones(1,size(vecx,2))/L2, 0);
zk = zk';
end
% function[zk] = prox_g2(x)
% mu_u = 0;
% mu_l = max(x);
% mu_m = (mu_u + mu_l)/2;
% 
% vec_diff = x - mu_m;
% val_mu = sum(vec_diff(vec_diff>0));
% iter = 0;
% while abs(val_mu - 1)>1e-13 && iter < 30000
% if val_mu > 1
%     mu_u = mu_m;
% else
%     mu_l = mu_m;
% end
% mu_m = (mu_u + mu_l)/2;
% 
% vec_diff = x - mu_m;
% val_mu = sum(vec_diff(vec_diff>0));
% iter = iter+1;
% end
% zk = vec_diff;
% zk(vec_diff<0) = 0;
% disp(['yoyo ', num2str(sum(zk)), ' iter ', num2str(iter)])
% end

% function[zk] = prox_g2(x)
% %proxPPXAplus calculate the proximity operator's value by PPXA+ 
% % D = HX
% % x = o_res (yVn)
% L = length(x);
% oneVn = ones(L,1);
% % zk = x + (0- sum(oneVn .* x))/norm(oneVn,2)^2 *oneVn;
%  zk = x + (1- sum(oneVn .* x))/norm(oneVn,2)^2 *oneVn;
%  zk(zk<0)=0;
% end

function[s_res, p_res] = correct_shift(s,p)
[~,i_max] = max(p);
L = length(p);
mid = ceil(L/2);
p_res = p;
s_res = s;
while i_max > mid
    p_res= circshift(p_res,-1);
    s_res = circshift(s_res, 1);
    [~,i_max] = max(p_res);
end
while i_max < mid
    p_res = circshift(p_res, 1);
    s_res = circshift(s_res,-1);
    [~,i_max] = max(p_res);
end

end