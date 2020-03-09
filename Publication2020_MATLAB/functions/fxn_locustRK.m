%%% Split step code for continuum (PDE) Locust Model
%%% Uses Runge-Kutta (4th order) for time derivative
%%% Uses simple shift for the spatial derivative
%%%%%%%%%%% Outputs: %%%%%%%%%%%
%   R = matrix (nx by nt) of resource densities
%   S = matrix (nx by nt) of stationary locust densities
%   M = matrix (nx by nt) of moving locust densities
%   rhomax = vector (nt) of maximum locust densities = max(S(:,t)+M(:,t))
%   Optional outputs
%   c = varargout{1} = speed of traveling wave over last 10 time steps
%   Rm = varargout{2} = min(R(:,end)) = minimum resources left behind
%%%%%%%%%%% Inputs: %%%%%%%%%%%
%%% Numerical parameters
%   L = length of spatial domain
%   nx = number of spatial points
%   nt = number of timesteps
%   video_flag = 1 ==> watch real-time video as matlab solves
%%% Biological Parameters
%   N = 7000;%7294.540405;       %total mass
%   Rplus = 200;%194.487762;     %initial resources
%   v = 0.04;%0.046614;           %speed of a moving locust (m/s)
%   lambda = 10^(-5);%10^(-5.002960);    %feeding rate of a stationary locust
%%%   Switching rates                                %nondimensional name
%     %     ksm = eta-(eta-alpha)*exp(-gamma.*R);
%     %     S -> M parameters in ksm, which is decreasing
%     alpha = .01*1/2;     %1/10*1.5;% 1; %small R limit, must be > eta
%     eta = 0.00001*1/2;       %1/10*.25;%0; %large R limit, must be < alpha   
%     gamma = .01;        %0; % expl rate                           
%     %     kms = theta-(theta-beta)*exp(-delta.*R);
%     %     M -> S parameters in kms, which is increasing
%     beta = .005*1/2;     %1/10*1.45;%0.5;   %small R limit, must be < theta   
%     theta = .04*1/2;     %1/10*1.5;%1.3 ;%large R limit, must be > beta     
%     delta = .001;        %2;%expl rate                         

function [R, S, M, rhomax, varargout] = fxn_locustRK(L,nx,nt,N,Rplus,v,lambda,alpha,beta,eta,theta,gamma,delta,video_flag)

close all
% clear all

% video_flag = 1;

% %Input parameters (the ones we like to vary)
% N = 10000;      %total mass
% Rplus = 200;    %1.5; %initial resources available
% 
% %%% Initialize parameters
% % %   Switching rates                                %nondimensional name
% %     ksm = eta-(eta-alpha)*exp(-gamma.*R);
% %     S -> M parameters in ksm, which is decreasing
%     alpha = .01*1/2;     %1/10*1.5;% 1; %small R limit, should be > beta     %nondim = 1
%     eta = 0.00001*1/2;       %1/10*.25;%0; %large R limit, should be < theta     %nondim = phi = eta/alpha
%     gamma = .01;        %0; % expl rate                           %nondim = 1
% %     kms = theta-(theta-beta)*exp(-delta.*R);
% %   M -> S parameters in kms, which is increasing
%     beta = .005*1/2;     %1/10*1.45;%0.5;   %small R limit, should be < alpha   %nondim = mu = beta/alpha
%     theta = .04*1/2;     %1/10*1.5;%1.3 ;%large R limit, should be > eta     %nondim = psi = theta/alpha
%     delta = .001;        %2;%expl rate                          %nondim = omega = delta/gamma
% % %Biological parameters
%     v = 0.005;       %speed of an individual, needed for dt~dx relation above
%     lambda = 5.0e-6;    %0.01; %feeding rate... so far poorly defined I think...

    
% 'Info for comparison to Sobol ranges'
% 
% beta
% 
% Delta = alpha/beta-eta/theta    % should be > 0
% 
% etaalpha = eta/alpha    % should be < 1
% thetabeta = theta/beta  % should be > 1

 
    
% % Setup computational grid    
% L = 60;             % Domain size    
% nx = 2000; %spatial points
% nt = 1000;  %time steps


dx = L/nx; %spatial gridsize
dt = dx/v;  %time step -- chosen to be dx/v for convenience

X = (1:nx)*dx; 


%for plotting
green =  [  0    0.5000      0];

%%% Initialize variables
% Initialize S, M, R
S = zeros(nx,nt);
M = zeros(nx,nt);
R = zeros(nx,nt);

CM = zeros(nt,1); % initialize vector to hold center-of-mass at each time
Rhotot = zeros(nt,1);

%%% Initial conditions
% R is Rplus
R(1:end,1) = Rplus;
% S and M are Gaussians centered at 10
cen = 10; % X(end)/10; <-- for changing domains
sig = 2.75/sqrt(2);
S(:,1) = N/2*normpdf(X,cen,sig);
M(:,1) = N/2*normpdf(X,cen,sig);

% S(:,1) = 1/sqrt(2*pi*sig^2)*exp(-(X-cen).^2./(2*sig^2));
% M(:,1) = 1/sqrt(2*pi*sig^2)*exp(-(X-cen).^2./(2*sig^2));
% S(:,1) = N/2*S(:,1);
% M(:,1) = N/2*M(:,1);

% fac = 1/3+0.03;
% S(:,1) = exp(-((X-cen)*fac).^2);
% M(:,1) = exp(-((X-cen)*fac).^2);
% % Normalize S and M 
% Stot = sum(S(:,1))*dx;
% Mtot = sum(M(:,1))*dx;
% S(:,1) = N/2*S(:,1)/Stot;
% M(:,1) = N/2*M(:,1)/Mtot;

% Calculate measurables
Stot = sum(S(:,1))*dx;
Mtot = sum(M(:,1))*dx;
initRhotot = Stot+Mtot; %total mass, should = N -- Rhotot should always = N
CM(1) = sum(X'.*(S(:,1)+M(:,1)))*dx/N;  %center of mass
rhomax(1) = max(S(:,1)+M(:,1));  %maximum
Rhotot(1) = initRhotot;

%%% Time Stepping
for t = 2:nt
    [S1, M1, R2] = rk4(S(:,t-1), M(:,t-1), R(:,t-1), dt/2, alpha, eta, gamma, beta, theta, delta, lambda);
    [S2, M2] = locustMove(S1,M1);
    [S(:,t), M(:,t), R(:,t)] = rk4(S2, M2, R2, dt/2, alpha, eta, gamma, beta, theta, delta, lambda);
%     max(abs(S(:,t)-S2))/mean(S2)
%     max(abs(M(:,t)-M2))/mean(M2)
%     max(abs(R(:,t)-R2))/mean(R2)
    
    % Calculate population totals
    Stot = sum(S(:,t))*dx;
    Mtot = sum(M(:,t))*dx;

    Rhotot(t) = Stot+Mtot; %total mass, should = N
    
    CM(t) = sum(X'.*(S(:,t)+M(:,t)))*dx/Rhotot(t); % center of mass
    rhomax(t) = max(S(:,t)+M(:,t));  %maximum
    
    if video_flag == 1 || t == nt
%         if t==nt
    % plot current stationary/moving locusts and resources
    figure(31)
    h=plot( ...X, M(:,t),...
          ...X, S(:,t),...
          X, S(:,t)+M(:,t),...
          X, R(:,t)*max(rhomax)/Rplus);
    axis([X(1) X(end) 0 1.1*max(rhomax)])
    h(end).Color=green;
%     %pause(0.1)
%         end
    end
end

% plot center of mass over time
figure(32)
plot(1:nt, rhomax);
title('peak')
movegui('northeast')
%
% figure(33)
% plot(1:nt, Rhotot)
% title('total mass')
%
figure(34)
plot(diff(CM)/dt)
title('speed')
movegui('northwest')

Rs = 0:.1:Rplus;
[ksm, kms] = fxn_switching(Rs, alpha, eta, gamma, beta, theta, delta);
figure(35)
plot(   Rs, ksm,...
        Rs, kms)
title('switching rates')
axis([0 Rplus min(beta,eta) max(alpha,theta)])
legend('ksm','kms')
movegui('southeast')


%Rhotot(nt) % print off final number of locusts (can compare to initRhotot)

c=(CM(nt)-CM(nt-10))/(10*dt) %print avg speed of last 10 steps
Rm = min(R(:,end));

varargout{1} = c;
varargout{2} = Rm;

tfinal = nt*dt

% Plot K(R) = ksm/c - kms/(v-c)
figure(36)
Kr = ksm/c - kms/(v-c);
plot(Rs,Kr,'r')
hold on
plot(Rs,zeros(length(Rs)),'--k')
Rminus=min(R(:,end));
yl = ylim;
plot([Rminus,Rminus],[yl(1),yl(2)],'-b')
legend('K(R)','K=0','Rminus')
xlabel('R'); ylabel('K(R)')
title('K(R) = ksm/c - kms/(v-c)')
movegui('southwest')

end

%%%%%%%%%%%%%%%%%%%
%%% local functions
%%%%%%%%%%%%%%%%%%%
% NOTE: While S, M, R are matrices above, the same symbols denote a single
% row of each matrix when appear as inputs in the functions below. A single
% row represents the S, M, R variable at the current time step only.

function [Snew, Mnew] = locustMove(S, M)
    % pass in S and M = row vectors (i.e. for current time step)
    % return S and M = row vectors, result of movement
    Snew = S;
    Mnew = [0; M(1:end-1)]; % move movers foward one, zero boundary cond.
end


% function [ksm, kms] = locustRate(R, alpha, eta, gamma, beta, theta, delta)
%     % pass in:
%         % R = row vector, resources for current time step only.
%         % greek letters = scalars, switching parameters
%     % return: ksm,kms = row vectors, switching rates to/from moving/stationary
%     
%     ksm = eta-(eta-alpha)*exp(-gamma.*R);
%     kms = theta-(theta-beta)*exp(-delta.*R);
% end

function [St,Mt,Rt] = f(S, M, R, alpha, eta, gamma, beta, theta, delta, lambda)
        % pass in:
            % S, M, R = row vectors
            % greek letters = scalars, switching parameters
            % lambda = scalar, feeding rate parameter
        % return:
            % St, Mt, Rt = row vectors, time-derivs at specific time step
            
        [ksm, kms] = fxn_switching(R, alpha, eta, gamma, beta, theta, delta);

        St = (-ksm.*S + kms.*M);
        Mt = (ksm.*S - kms.*M);
        Rt = (-lambda*R.*S);
end

function [up,vp,wp] = rk4(u,v,w,h,alpha, eta, gamma, beta, theta, delta, lambda)
    % pass in
        % u = current time S, row vector
        % v = current time M, row vector
        % w = current time R, row vector
        % h = time step
    % RK4 calls f which returns the derivative for S, M, and R.
    %
    % return
        %[up, vp, wp] = [S(:,t), M(:,t), R(:,t)], row vector values at time t
        
        [k1,l1,m1]=f(u,v,w,alpha, eta, gamma, beta, theta, delta, lambda);
		[k2,l2,m2]=f(u+h*k1/2,v+h*l1/2,w+h*m1/2,alpha, eta, gamma, beta, theta, delta, lambda);
		[k3,l3,m3]=f(u+h*k2/2,v+h*l2/2,w+h*m2/2,alpha, eta, gamma, beta, theta, delta, lambda);
		[k4,l4,m4]=f(u+h*k3,v+h*l3,w+h*m3,alpha, eta, gamma, beta, theta, delta, lambda);
		up=u+h*(k1+2*k2+2*k3+k4)/6;
		vp=v+h*(l1+2*l2+2*l3+l4)/6;
        wp=w+h*(m1+2*m2+2*m3+m4)/6;
end