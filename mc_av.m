function [Value, Error] = mc_av(S0,K,r,sigma,T,N,M,gamma,Z)
%   Compute European call option prices using Euler's method and Monte Carlo
%   Simulation with Anthetic Variate.
%
% [Value, Error] = mc_av(S0,K,r,sigma,T,N,M,gamma,Z)
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
% Error       - ...

dt = T/(N/2);  
S1 = S0*ones(M,1);
S2 = S0*ones(M,1);
for i=1:N/2 % loop for each time step
    dW = Z(:,i)*sqrt(dt);
	S1 = S1 + r*S1*dt + sigma*S1.^gamma.*dW;
    S2 = S2 + r*S2*dt + sigma*S2.^gamma.*-dW;
end
Value1 = exp(-r*T)*mean(max(S1-K,0)); %Expected value
Value2 = exp(-r*T)*mean(max(S2-K,0)); %Expected value
Value=(Value1+Value2)/2;
BS = blsprice(S0,K,r,T,sigma);
Error = abs(BS-Value);
end