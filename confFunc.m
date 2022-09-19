function [upper, est, lower] = confFunc(Data,alpha)
% Confidence interval
%   
%   [upper, est, lower] = confFunc(Data,alpha)
%
%   Inputs:
%       Data     - Simulated prices
%       alpha    - alpha
%
%   Output:
%       upper - Upper bound 
%       est   - Estimated price
%       lower - Lower bound
zAlphas = norminv([alpha/2, 1-alpha/2]);
m       = length(Data);
sEst    = std(Data);
est     = mean(Data);
upper   = est + zAlphas(2) * sEst/sqrt(m);
lower   = est + zAlphas(1) * sEst/sqrt(m);
end