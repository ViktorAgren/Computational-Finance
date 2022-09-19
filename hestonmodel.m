function [S] = hestonmodel(S0,V0,r,kappa,theta,sigma,rho,T,n)
% Heston model 
%   
%   [S] = hestonmodel(S0,V0,r,kappa,theta,sigma,rho,T,n)
%
%   Inputs:
%       S0    - Initial stock price
%       V0    - the initial volatility
%       r     - Interest rate
%       kappa - Rate at which variance reverts to theta
%       theta - Long variance
%       sigma - Variance
%       rho   - Volatility of the volatility
%       T     - Years until maturity
%       n     - Number of timesteps
%
%   Output:
%       S - Price of the asset at every timestep
dt=T/n;
S     = zeros(1,n+1);
V     = zeros(1,n+1);
S(1)  = S0;
V(1)  = V0;
Z     = randn(2,n+1); % Brownian motion 
Z1    = Z(1,:)*rho+Z(2,:)*sqrt(1-rho*rho); % Correlated Brownian motion
    for i = 1:n
        %S(i+1) = S(i)+(r - V(i)/2)*S(i)*dt+sqrt(abs(V(i)))*S(i)*sqrt(dt)*Z(1,i); 
        S(i+1) = S(i)*exp((r - V(i)/2)*dt+sqrt(abs(V(i)))*sqrt(dt)*Z(1,i)); 
            % Stochastic process for the price of the asset
        V(i+1) = V(i)+kappa*(theta-abs(V(i)))*dt+sigma*sqrt(abs(V(i)))*sqrt(dt)*Z1(1,i);
            % The instantaneous variance
    end
end

