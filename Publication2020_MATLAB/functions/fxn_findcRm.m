% % % % % % % % % % % % % % % % % %
%
% Produce plots of the traveling wave profile
%
% Given the 10 inputs (v,lambda,Rp,N) and the greeks
%
% Produce c, Rm 
%
% Then compute the profile and
%   (width, peak, skew)
% 
% % % % % % % % % % % % % % % % % % 
%     close all
%     clear all

function [Rm, c] = fxn_findcRm(N,Rp,v,lambda,alpha,beta,eta,theta,gamma,delta)

%  ratmax is the switchover point for Rp/Rm from root finding to
%  asymptotic formulation
    ratmax=10^(20);   
% % % %%% Biological Parameters
% %     N = 9502.868652;      % total mass
% %     Rp = 209.351196;   % initial resources
% %     v = 0.041886;         % speed of a moving locust (m/s)
% %     loglam=-5.084778;
% %     lambda =10^(loglam);%feeding rate of a stationary locust
% % % 
% % %   The greeks
% % %
% %     gamma = 0.037436 ;  %S->M
% %     delta = 0.029516;   %M->S
% % %    
% %     Delta = 0.111419; % alpha/beta-eta/theta  (should be >0
% %     eoa = 0.309109;  %  eta/alpha  (should be <1)
% %     tob=9.762344;    %  theta/beta (should be >1)
% %     beta = 0.020937; 
% % %
% %     alpha =beta*Delta*(tob/(tob-eoa));
% % %
% %     eta = eoa*alpha;%	%S->M
% %     theta = tob*beta;   %M->S

%%%%%%%%%%%%%%%%% Test case 1 %%%%%%%%%%%%%%%
%     alpha =2.0; eta = 1.0; gamma = .01;
%     beta = 1; theta = 2; delta = .01
%     v = 0.45;   lambda =1e-5;
%     N = 10000;  Rp = 200;
% %
% %  From PDE code
% %     Rmin is: 142.2886
% %     The speed (c) is: 0.17764
%%%%%%%%%%%%%% test case 2 %%%%%%%%%%%%%%%%%%%%
% alpha =2.0; eta = 1.0;  gamma = .04;
% beta = 4.0; theta = 5;  delta = .01;
% v = 0.45;  lambda =1e-5;
% N = 10000; Rp = 200; 
% 
% % From PDE code
% %     Rmin is: 71.5707
% %     The speed (c) is: 0.079869
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% First solve Rm, then solve for c
%   
%  Redefine function for root finding as a single-variable function
%
    fun = @(rat) rootfind(rat,Rp,N,v,lambda,alpha,eta,gamma,beta,theta,delta);  
%
if fun(ratmax) <0  % Will root finding fail ??
    ratout=fzero(fun,[1+10^(-8),ratmax]);
    Rm =Rp/ratout ;
    c = speed(Rp,Rm,alpha,eta,gamma,beta,theta,delta,v);
    disp('Using zero finding')
else
    logRpRm=Lasymp(Rp,N,v,lambda,alpha,eta,gamma,beta,theta,delta);
    c=casymp(logRpRm,Rp,v,alpha,eta,gamma,beta,theta,delta);
    Rm=Rp/exp(logRpRm);
    disp('Using asymptotics')
end

format long e
disp(['Rm is ',num2str(Rm)]);
disp(['c is ',num2str(c)]);

    
function out = Ival(Rp,Rm,c1,c2,c3)
    out = c1*(log(Rp)-log(Rm)) - (c1-c2)*(expint(c3.*Rm)-expint(c3.*Rp));
end

function out = rootfind(rat,Rp,N,v,lambda,alpha,eta,gamma,beta,theta,delta) 
         out= (N*lambda/v)*Ival(Rp,Rp/rat,theta,beta,delta) ...
                    -log(rat)*Ival(Rp,Rp/rat,eta,alpha,gamma);
end

function c = speed(Rp,Rm,alpha,eta,gamma,beta,theta,delta,v)
     rat = Ival(Rp,Rm,eta,alpha,gamma)./Ival(Rp,Rm,theta,beta,delta);
     c = v*rat./(1+rat);
end

 function lograt=Lasymp(Rp,N,v,lambda,alpha,eta,gamma,beta,theta,delta)
 a=alpha;
 b=(eta-alpha)*Fasymp(Rp,gamma) - (N*lambda/v)*beta;
 chere=-(N*lambda/v)*(theta-beta)*Fasymp(Rp,delta);
 lograt=(-b +sqrt(b^2-4*a*chere))/(2*a);
 end
 
 function F=Fasymp(Rp,gam)
 % euler = 5.772156649015329e-01  % Euler-Mascheroni constant
 euler=double(eulergamma);
 F=euler+log(Rp*gam)+expint(Rp*gam);
 end
 
 function c = casymp(logRpRm,Rp,v,alpha,eta,gamma,beta,theta,delta)
     I1overI2 = (alpha*logRpRm+(eta-alpha)*Fasymp(Rp,gamma))./(beta*logRpRm+(theta-beta)*Fasymp(Rp,delta));
     c = v*I1overI2./(1+I1overI2);
end

end