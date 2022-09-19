function [output] = UOPut(S0,K,b,T,r,a,sigma,m)
% Discretely monitored up-and-out put option in Black–Scholes model.
%   
%   [EstPut] = UOPut(S0,K,b,T,r,a,sigma,m)
%
%   Inputs:
%       S0    - Initial stock price
%       K     - Strike Price
%       b     - barrier
%       T     - Years until maturity
%       r     - Interest rate
%       a     - Alpha
%       sigma - Variance
%       m     - Number of simulations
%
%   Output:
%       EstPut - Estimated price for up-and-out put option

n     = round(252*T); % Working-days until maturity
payoffP = zeros(1,m);
pEst    = zeros(1,m);
for j = 1:m
    % Black Scholes model
    [S] = BlackScholes(S0,T,r,sigma);
    
    % Up-and-out barrier 
    if max(S) < b
        payoffP(j) = max(K-S(end), 0); % Put option
    else
        payoffP(j) = 0;
    end
    pEst(j)   = exp(-r*T)*mean(payoffP(1:j)); % Estimated stock price
end
% Confidence interval
[u,est,l]=confFunc(pEst,a);
A1 = [l est u];
A2= {'Lower bound','Estimated price','Upper bound'};
output = array2table(A1,'VariableNames',A2);
end

