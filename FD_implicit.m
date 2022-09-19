function [Value,Error] = FD_implicit(S0,Smin,Smax,K,r,sigma,T,gamma,M,N)
% Value of a European Call option and error between the analytical and the implicit
% finite difference method.
%
% [Value,Error] = FD_implicit(S0,Smin,Smax,K,r,sigma,T,gamma,M,N)
%
% Inputs:
% S0       - Initial stock price
% Smin     - Minimum value of S
% Smax     - Maximum value of S
% K        - Strike (i.e., exercise) price of the option.
% r        - Annualized continuously compounded risk-free rate of return
%            over the life of the option, expressed as a positive decimal
%            number.
% sigma    - Annualized asset price volatility (i.e., annualized standard
%            deviation of the continuously compounded asset return),
%            expressed as a positive decimal number.
% T        - Time to expiration of the option, expressed in years.
% M        - Number of time steps
% N        - Number of stock steps
%
% Output: 
% Value    - Price of a European call option
% Error    - Absolute error between the Black-Scholes analytical price and
%            the implicit Forward Difference method
dS = (Smax - Smin)/N;
dt = T/M;
t = 0:dt:T;
S = Smin:dS:Smax;

% Grid in space and time
v=zeros(M+1,N+1);

v(end,:) = max(0,S-K); %V(T,s)=phi(s)
v(:,end) = Smax - K*exp(-r*(T-t)); %V(t,Smax)=S-Ke^-r(T-t)

% tridiagonal coefficients
alpha = 0.5*sigma^2*S.^(2*gamma)/dS^2*dt;
beta  = r*S*dt/(2*dS);

l = -alpha + beta;
d = 1 + r*dt + 2*alpha;
u = -beta - alpha;

% A implicit
A = diag(l(3:N-1),-1)+diag(d(2:N-1))+diag(u(2:N-2),1);

for k = M:-1:1
    v(k,2:N-1) = A\(v(k+1,2:N-1)'- u(N-1)*v(k,N));
end
BS = blsprice(S0, K, r, T, sigma); % Analytical price
Value = interp1(S,v(1,:),S0); % Estimated
Error = abs(Value-BS); % Absolute error
end