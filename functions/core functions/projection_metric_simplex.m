function [ P ] = projection_metric_simplex( A, D, epsilon )
%
% function [ P ] = projection_metric_simplex( A, D )
%--------------------------------------------------------------------------
% Computes the point in the unit simplex conv(e_1,...,e_n) with minimal 
% distance to an arbitrary given point.
% Projection sur le simplex pour la normalisation de l'op?rateur de flou
%
% ENTREE : 
% A : \in R^{p \times n}
% D : \in R^{p \times n}
%
% SORTIE : 
% P : projection (relative a metrique induite par D) de A sur le simplex 
%--------------------------------------------------------------------------


[p,n]=size(A); 

D12 = sqrt(D) ; 
D12inv = 1./D12 ;
Dinv = 1./D;

DA = D12.*A ;

epsilonD = epsilon.*D12;

%projection of p on hyperplane w_1 + ... + w_n = 1
lambda = (sum(A,2)-1)./sum(Dinv,2) ;
P = DA-lambda(:,ones(n,1)).*D12inv;
% % 
subzerodetect = zeros(p,n);
subzerodetect(P<epsilonD) = 1;

while(sum(subzerodetect(:)) > 0)
    % while projection is not in the simplex:
    % project orthogonally on to the corresponding coordinate hyperplane and
    % and project along the corresponding direction on to the simplex
    direction = zeros(p,n);
    direction(P>epsilonD) = 1;
    P(P<epsilonD)=epsilonD(P<epsilonD);
    
    denominator = sum(Dinv.*direction,2);
    lambda = (sum(D12inv.*P,2)-1) ./ denominator;
    
    P = P - lambda(:,ones(n,1)).* (D12inv .* direction);
    subzerodetect = zeros(p,n);
    subzerodetect(P<epsilonD) = 1;
end

P = P.*D12inv; 
