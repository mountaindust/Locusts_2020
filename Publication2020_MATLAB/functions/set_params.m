%%% A script converting between versions of our parameters

%%% Biological Parameters
N = 7000;%7294.540405;       %total mass
Rplus = 200;%194.487762;     %initial resources
v = 0.04;%0.046614;           %speed of a moving locust (m/s)
lambda = 10^(-5);%10^(-5.002960);    %feeding rate of a stationary locust
% Sobol version of Switching Rate parameters
    % M -> S low resource
    beta = 0.02;%0.015695;
    % S -> M rate: ratio of large R/small R
    etaalpha = 0.8;%0.759527; % < 1
    % M -> S rate: ratio of large R/small R
    thetabeta = 7;%7.215061; % > 1
    % alpha/beta - eta/theta = small R ratio SM/MS - large R ratio SM/MS
    DELTA = 10^(-0.7);%10^(-0.653073);
    % exponential rates
    gamma = 0.03;%0.034318;     %S->M
    delta = 0.005;%0.004681;     %M->S
    
%%% call a function to convert these to actual greek parameter values
[alpha, beta, eta, theta] = fxn_paramRats(beta,thetabeta,etaalpha,DELTA);
    
% %%% for PDE
% %%% Numerical Parameters
% L = 100;     %domain length <---- Increase this to run longer!
% nx = 1000;  %spatial points
% nt = 1000;  %time steps
% %%% Computational grid
% dx = L/nx;  %spatial discretization
% dt = dx/v;  %time step -- chosen to for convenience
% X = (1:nx)*dx;
% 
% %%% Report Info
% tfinal = nt*dt
% 
% [R, S, M, rhomax] = fxn_locustRK_shift(L,nx,nt,N,Rplus,v,lambda,alpha,beta,eta,theta,gamma,delta,video_flag);
% % Returns
% % R,S,M as space,time arrays