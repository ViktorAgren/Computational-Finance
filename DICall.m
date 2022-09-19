function [output] = DICall(S0,V0,K,b,T,r,a,sigma,kappa,theta,rho,m)
% Discretely monitored down-and-in call option under the Heston model 
%   
%   [output] = DICall(S0,V0,K,b,T,r,a,sigma,kappa,theta,rho,m)
%
%   Inputs:
%       S0    - Initial stock price
%       V0    - the initial variance
%       K     - Strike Price
%       b     - barrier
%       T     - Years until maturity
%       r     - Interest rate
%       a     - Alpha
%       sigma - Variance
%       kappa - Rate at which variance reverts to theta
%       theta - Long variance
%       rho   - Volatility of the volatility
%       m     - Number of simulations
%
%   Output:
%       output - Estimated price for down-and-in call option with lower and
%       upper confidence interval
n1       = round(252*T); % Working-days until maturity for the larger sample
n2       = 1/(2*T/n1); % Working-days until maturity for the smaller sample
payoffC1 = zeros(1,m);
payoffC2 = zeros(1,m);
cEst    = zeros(2,m);
cEstp    = zeros(1,m);
for j = 1:m 
    % Heston models
    [S1] = hestonmodel(S0,V0,r,kappa,theta,sigma,rho,T,n1);
    [S2] = hestonmodel(S0,V0,r,kappa,theta,sigma,rho,T,n2);
    
    % Barrier payoff
	if min(S1) < b
        payoffC1(j) = max(S1(end)-K, 0);
    else
        payoffC1(j)=0;
    end
	if min(S2) < b
        payoffC2(j) = max(S2(end)-K, 0);
    else
        payoffC2(j)=0;
    end
    
    % Estimated prices
	cEst(1,j)	= exp(-r*T)*mean(payoffC1(1:j));
    cEst(2,j)	= exp(-r*T)*mean(payoffC2(1:j));
    
    % Romberg extrapolation
    cEstp(j)    = 2*cEst(1,j)-cEst(2,j);
end
% Confidence interval
[u,est,l]=confFunc(cEstp,a);
A1 = [l est u];
A2= {'Lower bound','Estimated price','Upper bound'};
output = array2table(A1,'VariableNames',A2);
end