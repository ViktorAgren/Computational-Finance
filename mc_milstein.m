function [Value, STD] = mc_milstein(S0,K,r,sigma,T,N,M,gamma)
%   Compute European call option prices using Milstein scheme and Monte Carlo 
%   Simulation.
%
% [Value, STD] = mc_milstein(S0,K,r,sigma,T,N,M,gamma)
%
% Inputs:
% S0          - Current price of the underlying asset.
% K           - Strike (i.e., exercise) price of the option.
% r           - Annualized continuously compounded risk-free rate of return
%               over the life of the option, expressed as a positive decimal
%               number.
% T           - Time to expiration of the option, expressed in years.
% sigma       - Annualized asset price volatility (i.e., annualized standard
%               deviation of the continuously compounded asset return),
%               expressed as a positive decimal number.
% N           - Number of timesteps
% M           - Number of simulations
% Gamma       - controls the relationship between volatility and price
%               [0,1]
% Output:
% Value       - Price (i.e., value) of a European call option.
% STD         - Standard Deviation

dt = T/N;  
S = S0*ones(M,1);
V = zeros(1,M);
for k=1:M % Loop for number of simulations
    Z  = randn(1,N);
    for i=1:N % Loop for each time step
        S(k) = S(k) + r*S(k)*dt+sigma*S(k)^gamma*sqrt(dt)*Z(i) +...
            0.5*gamma*sigma^2*S(k)^(2*gamma-1)*sqrt(dt)*(Z(i)^2-1);
    end
    V(k) = max(S(k,end)-K,0); %Phi
end
Value = exp(-r*T)*mean(V); %Expected value
STD = exp(-r*T) *std(V);
end