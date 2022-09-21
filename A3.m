clear;clc;
randn('state',0)
S0=14;K=15;r=0.1;sigma=0.25;T=0.5;N=2024;M=1000000;gamma=1;Smin=0;Smax = 4*K;

%% Simulate sample error MC vs MC with Antithetic variates
Mvec=round(logspace(2,5,100));

Error_mc = zeros(size(Mvec,2),1);
Error_mc_antithetic = zeros(size(Mvec,2),1);
for i=1:size(Mvec,2)
    Z=randn(Mvec(i),N);
    [~, Error_mc(i)] = mc(S0,K,r,sigma,T,N,Mvec(i),gamma,Z);
    [~, Error_mc_antithetic(i)] = mc_av(S0,K,r,sigma,T,N/2,Mvec(i),gamma,Z);
end
figure(1)
loglog(Mvec,Error_mc,'bd')
p_mc = polyfit(log(Mvec), log(Error_mc), 1);
hold on
loglog(Mvec,exp(polyval(p_mc,log(Mvec))),'b','LineWidth',2)

loglog(Mvec,Error_mc_antithetic,'rd')
p_mc_antithetic = polyfit(log(Mvec), log(Error_mc_antithetic), 1);
loglog(Mvec,exp(polyval(p_mc_antithetic,log(Mvec))),'r','LineWidth',2)
hold off
xlabel('log(No. of simulations)')
ylabel('log(Error)')
legend('Monte Carlo sample points','Monte Carlo regression line',...
    'Antithetic Variates sample points','Antithetic Variates regression line')
%% Implicit FD Error (stock)
Nt = 10000; Ns = 120*(1:10);
Error_fd = zeros(size(Ns,2),1);
for i=1:size(Ns,2)
    %[~, Error_fd(i)] = FD_implicit(S0,Smin,Smax,K,r,sigma,T,gamma,M_FD,Nvec_FD(i));
    [~, Error_fd(i)] = FDimplicit(S0,Smin,Smax,K,r,sigma,T,gamma,Ns(i),Nt);
end
figure(2)
loglog(Ns,Error_fd,'bd')
hold on
p_fd = polyfit(log(Ns), log(Error_fd), 1);
loglog(Ns,exp(polyval(p_fd,log(Ns))),'b')
hold off
xlabel('log(No. of stock steps)')
ylabel('log(Error)')
xlim([0,Ns(end)])
legend('Implicit FD sample points','Implicit FD regression line')

%% Implicit FD Error (time)
Ns = 1000; Nt_vec = 100*(1:20);
Error_fd = zeros(size(Nt_vec,2),1);
for i=1:size(Nt_vec,2)
    [~, Error_fd(i)] = FD_implicit(S0,Smin,Smax,K,r,sigma,T,gamma,Mvec_FD(i),N_FD);
end
figure(3)
loglog(Nt_vec,Error_fd,'rd')
hold on
p_fd = polyfit(log(Nt_vec), log(Error_fd), 1);
loglog(Nt_vec,exp(polyval(p_fd,log(Nt_vec))),'r')
xlabel('log(No. of time steps)')
ylabel('log(Error)')
legend('Implicit FD sample points','Implicit FD regression line)')
xlim([0,Nt_vec(end)])
hold off
%% Runtime MC
sim           = 100;
N             = 1000;
Mvec          = round(logspace(2,5,100));
Time_mc = zeros(1,sim); 
Time_mc_av = zeros(1,sim);
for i = 1:sim
    tic
    Z          = randn(Mvec(i),N);
    Value_mc   = mc(S0,K,r,sigma,T,N,Mvec(i),gamma,Z);
    Time_mc(i) = toc;
    tic
    Z             = randn(Mvec(i),N/2);
    Value_mc_av   = mc_av(S0,K,r,sigma,T,N/2,Mvec(i),gamma,Z);
    Time_mc_av(i) = toc;
end
figure(4)
loglog(Mvec,Time_mc,'r*')
hold on
p_fd = polyfit(log(Mvec), log(Time_mc), 1);
loglog(Mvec,exp(polyval(p_fd,log(Mvec))),'r')

loglog(Mvec,Time_mc_av,'b*')
p_fd = polyfit(log(Mvec), log(Time_mc_av), 1);
loglog(Mvec,exp(polyval(p_fd,log(Mvec))),'b')
xlabel('log(No. of simulations)')
ylabel('log(Simulation time)')
legend('Monte Carlo','Monte Carlo regression','Antithetic variates',...
    'Antithetic variates regression')
hold off
%% Runtime FD
sim           = 100;
Ns            = 1000;
Nt_time       = (1:sim)*100;
Time_implicit = zeros(1,sim); 
Time_explicit = zeros(1,sim);
for i = 1:sim
    tic
    Value_implicit = FD_implicit(14,0,Smax,K,r,sigma,T,gamma,Nt_time(i),Ns);
    Time_implicit(i) = toc;
    tic
    Value_explicit = FD_explicit(14,0,Smax,K,r,sigma,T,gamma,Nt_time(i),Ns);
    Time_explicit(i) = toc;
end
figure(5)
loglog(Nt_time,Time_implicit,'r*')
hold on
p_fd = polyfit(log(Nt_time), log(Time_implicit), 1);
loglog(Nt_time,exp(polyval(p_fd,log(Nt_time))),'r')

loglog(Nt_time,Time_explicit,'b*')
p_fd = polyfit(log(Nt_time), log(Time_explicit), 1);
loglog(Nt_time,exp(polyval(p_fd,log(Nt_time))),'b')
xlabel('log(No. of stock steps)')
ylabel('log(Simulation time)')
legend('Implicit','Implicit regression','Explicit','Explicit regression')
hold off

%% Gamma increase 
gamma = 0:0.1:1;Ns=1000;Nt=1000;N=10000;M=10000;
Value_FD = zeros(size(gamma,2),1);
Value_mc = zeros(size(gamma,2),1);
Value_mc_av = zeros(size(gamma,2),1);
Value_mc_milstein = zeros(size(gamma,2),1);
Z = randn(N,N);
for i = 1:size(gamma,2) % Error convergence for stock steps
    [Value_FD(i),~] = FD_implicit(S0,Smin,Smax,K,r,sigma,T,gamma(i),Nt,Ns);
    [Value_mc(i),~] = mc(S0,K,r,sigma,T,N,N,gamma(i),Z);
    [Value_mc_av(i), ~] = mc_av(S0,K,r,sigma,T,N,N,gamma(i),Z);
    [Value_mc_milstein(i), ~] = mc_milstein(S0,K,r,sigma,T,N,N,gamma(i));
end

figure(6)
subplot(2,2,1)
plot(gamma,Value_FD)
title('Implicit')
ylabel('Price')
xlabel('Gamma')
subplot(2,2,2)
plot(gamma,Value_mc)
title('MC')
ylabel('Price')
xlabel('Gamma')
subplot(2,2,3)
plot(gamma,Value_mc_av)
title('AV')
ylabel('Price')
xlabel('Gamma')
subplot(2,2,4)
plot(gamma,Value_mc_milstein)
title('Milstein')
ylabel('Price')
xlabel('Gamma')