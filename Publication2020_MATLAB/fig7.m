% Creates figures comparing the geometry of profiles of solutions to the
% ABM and PDE, for varying log(lambda).

close all
clear all

addpath('functions','data')

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
close all
clear all

% fonts
titlesize = 12;
axislabelsize = 12;%14;
ticklabelsize = 10;%12;
subfiglabelsize = 12;

%colors
ksm_color = 1/255*[237, 143, 33];   %gold
kms_color = 1/255*[111,73,179];     %purple
green =  [  0    0.5000      0];
gold = [0.9290, 0.6940, 0.1250];
orange = [0.85 0.325 0.098];
blue = [      0    0.4470    0.7410];
defwhite = [1 1 1];

%linewidth
lwidth = 2;

%%% options
newdata = 0; %1 saves new data, 0 loads old data file and makes new picutre
savePlot = 1; % Save new .eps files

pcheck = 0; % pause after each run to visually check the profile
pcheck_plots = 0;

PDE_flag = 1; % 1 to run the PDE
ABM_flag = 1; % 1 to run the ABM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 7 %%%
number = 7;
filename = sprintf('fig%d',number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% biological parameters
set_params

%%% for PDE
%%% Numerical Parameters
nt = 200000;  % time steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set a range of log(lambda) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Manually run this on three intervals
loglamsNew = -7:-0.1:-8; %-7:0.1:-5.7; %-4:-0.1:-6.3;
lams = 10.^(loglamsNew); 

n = length(lams);


% get an intial condition from old data
load('data/cfShape_9:27.mat')

ind = find(loglams == loglamsNew(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruct Numerical Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% for PDE
    % have X = spatial grid = linspace(-L/2, L/2, nx);
    X = PDE{4,ind};
    nx = length(X);
    L = 2*X(end);
    dx = L/(nx-1);
    if abs(dx - mean(diff(X))) >= 10^-4 % just a small threshold number
        disp('dx value does not match the loaded spatial grid X!')
        pause
    end
    dt = dx/v;
    
%%% for ABM
%     T_end = nt*dt;
%     delta_t = 1;
%     T_steps_end = ceil(T_end/delta_t);
    T_steps_end = nt;
    delta_t = dt;
    delta_x = v*delta_t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% set PDE initial conditions
    Sin = PDE{1,ind};
    Min = PDE{2,ind};
    Rin = PDE{3,ind};
    X = PDE{4,ind};
    
%%% set ABM initial conditions
    agents_out = ABM{1,ind};
    R_out = ABM{2,ind};

    % first shift initial condition
        agents_in = agents_out;
        % location of first/last locust
        first_loc = min(agents_in(:,1));
        last_loc = max(agents_in(:,1));
        % shift so that first locust is at first spatial index
        agents_in(:,1) = agents_in(:,1)-first_loc+1;

        domain_length = T_steps_end;
        R_in = ones(domain_length,1)*Rplus;
        % shift so that first (last_loc-first_loc) values of R come from R_out
        R_in(1:1+last_loc-first_loc) = R_out(first_loc:last_loc);

% Peek at the initial conditions
figure(101)

subplot(2,1,1)
green =  [  0    0.5000      0];
plot(   X, Sin+Min)
ylabel('Locust Density, PDE')
xlabel('Initial Condition')

yyaxis right
plot(X, Rin, 'Color', green)
ylim([0,Rplus*1.1])
set(gca,'ycolor',green)

subplot(2,1,2)
[cnts, edges] = histcounts(agents_in(:,1),'BinMethod','Integers');
    x_ind = floor(edges(2:end));
    x_grid = x_ind*delta_x;
    plot(x_grid,cnts/delta_x)
    %xlim([x_grid(1), x_grid(end)])
    ylabel('Locust Density','Color','k')
    xlabel('Space (not spatial grid, not shifted), INIT')
    
    %add resource curve corresponding to histogram
    yyaxis right
    %resources(min(agents(:,1))-max_idx+mode_pt:max(agents(:,1))-max_idx+mode_pt),...
    plot(x_grid,R_in(x_ind),...
        ':','Color',green,'LineWidth',3);
    ylim([0,Rplus*1.1])
    set(gca,'ycolor',green)
if pcheck == 1
    disp('Peek at the initial conditions, press any key to continue')
    pause
end

clear('widths','peaks','skews','PDE','ABM')
% initialize storage variables
% top row for PDE, bottom row for ABM
widths = zeros(2,n);
peaks = zeros(2,n);
skews = zeros(2,n);

PDE = cell(4,n); % ending values of { S; M; R; X }
ABM = cell(5,n); % { agents_out; R_out; mean_pulse; mean_resources; xmesh }

for k = 1:n
    close all
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Run Script %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if newdata == 1
        
    % set current lambda
    lambda = lams(k);
    disp(['lambda is now ' num2str(lambda)])
    disp(['log10(lambda) is now ' num2str(log10(lambda))])
        
        if PDE_flag == 1

        %%% call PDE
        tic
                [R, S, M, c] = fxn_locustRK_shift(Sin,Min,Rin,L,nx,nt,N,Rplus,v,lambda,alpha,beta,eta,theta,gamma,delta,pcheck_plots);
        disp('time for PDE run was')
        toc
        %%% view end profile to check achieved end of transients and did not run off domain
                if pcheck == 1
                 disp('Check PDE profile and peak, then press a key to continue')
                 pause %until keypress %view end profile
                end
                
        %%% RECORD DATA
                % end locust density profile (stationary + moving)
                rho = S(:,end) + M(:,end);
                mass = sum(rho)*dx;
                disp(['Mass after run is ', num2str(mass)])
                
                % record model outputs
                PDE{1,k} = S(:,end);
                PDE{2,k} = M(:,end);
                PDE{3,k} = R(:,end);
                PDE{4,k} = X;

                %%% compute, record peak
                [m, max_idx] = max(rho);  % m = maximum, max_idx = index of maximimum (location of peak)
                % Interpolate to get a more accurate maximum
                    ym = rho(max_idx-1);
                    y0 = rho(max_idx);
                    yp = rho(max_idx+1);
                    interpmax = (16*y0^2-8*y0*(yp+ym)+(ym-yp)^2)/(16*y0-8*(ym+yp));
                peaks(1,k) = interpmax;

                %%% compute, record standard deviation
                cm = sum(X.*(rho)/N)*dx; % center of mass
                width = sqrt(sum((X-cm).^2.*(rho)/N)*dx); %std of the profile
                widths(1,k) = width;

                %%% compute, record skewness
                skews(1,k) = -sum((X-cm).^3.*(rho)/N)*dx/(widths(1,k))^3;
                % skews = skewness of end profile

        %%% Resize Computational Grid
                L = max([20*width, 50]);         % 10 std should get most of the profile, right?
                dx = L/(nx-1);  %spatial discretization
                dt = dx/v;  %time step -- chosen to for convenience
                Xnew = linspace(-L/2,L/2,nx);
                Xnew = Xnew';

                %%% interpolate
                % interp1(grid, samples on grid,new grid, method, extrapolation)
                Sint = interp1(X,S(:,end),Xnew,'pchip',0); % pchip = shape-preserving piecewise cubic interp
                Mint = interp1(X,M(:,end),Xnew,'pchip',0); % 0 ==> outside of domain grid, assign value 0 
                Rint = interp1(X,R(:,end),Xnew,'pchip',Rplus);
                % Xold = X; %if needed get it from spaces
                X = Xnew;

                % check mass
                % NOTE: mass before interpolation is guaranteed to be N, due to rescaling
                % written into the PDE code.
                mass = sum(Sint+Mint)*dx;
                disp(['Mass after interpolation is ', num2str(mass)])

        %%% reset initial condition for next lambda step
                Sin = Sint*N/mass;
                Min = Mint*N/mass;
                Rin = Rint;
        end
        
        if ABM_flag == 1
            
        %%% call ABM
        tic
            % optional outputs: mean_pulse, mean_resources, xmesh
            %(must change transient detection threshold to obtain)
            [agents_out, R_out, shape_stats, xmesh, mean_pulse, mean_resources] =...
                 fxn_ABM_ave(agents_in,R_in,T_steps_end,delta_t,N,Rplus,v,lambda,alpha,beta,eta,theta,gamma,delta,pcheck_plots);
        disp('time for ABM run was')
        toc
        %%% view mass, speed, resources to ensure that it has
        %%% equilibriated and not run off edge of domain
            if pcheck == 1 
             disp('Check ABM peak, then press a key to continue')
             pause %until keypress %view end profile
            end 
            
        %%% Record Data
            % record model outputs
            ABM{1,k} = agents_out;
            ABM{2,k} = R_out;
            
            ABM{3,k} = mean_pulse;
            ABM{4,k} = mean_resources;
            ABM{5,k} = xmesh;

            % take max by interpolating across three points in the mean pulse
            [m, max_idx] = max(mean_pulse);  %maximum
                ym = mean_pulse(max_idx-1);
                y0 = mean_pulse(max_idx);
                yp = mean_pulse(max_idx+1);
                max_pulse = (16*y0^2-8*y0*(yp+ym)+(ym-yp)^2)/(16*y0-8*(ym+yp));
            %%% record max, scaled by bin size
            peaks(2,k) = max_pulse;
            
            %%% compute width
            CM = sum(xmesh.*mean_pulse/N)*delta_x;
            widths(2,k) = sqrt(sum((xmesh-CM).^2.*mean_pulse/N)*delta_x); % check this!!!
            % don't multiply by delta_x
            
            %%% compute skewness
            skews(2,k) = -1*sum( (xmesh-CM).^3 .* mean_pulse/N )*delta_x/(widths(2,k))^3; % check this!!!
            
            % assign next initial condition
            agents_in = agents_out;
            
            % location of first/last locust
            first_loc = min(agents_in(:,1));
            last_loc = max(agents_in(:,1));
            % shift so that first locust is at first spatial index
            agents_in(:,1) = agents_in(:,1)-first_loc+1;
            
            R_in = ones(domain_length,1)*Rplus;
            % shift so that first (last_loc-first_loc) values of R come from R_out
            R_in(1:1+last_loc-first_loc) = R_out(first_loc:last_loc);
            
            %%% Rescale delta_t according to PDE rescaling
            delta_t = dt;
            delta_x = v*delta_t;
            
        end
        
        % SAVE DATA
    save('data/cfShape.mat','peaks','widths','skews','loglamsNew','PDE','ABM')
    
    else
        %load three data files
        
        %file 1
        load('data/cfShape_from7to8.mat','peaks','widths','skews')
        peaks1 = fliplr(peaks);
        widths1 = fliplr(widths);
        skews1 = fliplr(skews);
        
        %file 2
        load('data/cfShape_from7to57.mat','peaks','widths','skews')
        peaks2 = peaks(:,2:end);
        widths2 = widths(:,2:end);
        skews2 = skews(:,2:end);
        
        %file 3
        load('data/cfShape_from4to63.mat','peaks','widths','skews')
        % only go up to 5.7 ==> first 17 entries
        peaks3 = peaks(:,1:17); peaks3 = fliplr(peaks3);
        widths3 = widths(:,1:17); widths3 = fliplr(widths3);
        skews3 = skews(:,1:17); skews3 = fliplr(skews3);
        
        % concatenate them all together!
        peaks = [peaks1 peaks2 peaks3];
        widths = [widths1 widths2 widths3];
        skews = [skews1 skews2 skews3];
        
        loglamsNew = -8:0.1:-4;
        
        % also load data from Snapshots (created in fig8.m)
        load('data/cfShape_snap_9by27.mat')
        % loads PDE_snap, ABM_snap, peaks_snap, widths_snap, skews_snap
%         S = PDE_snap{1,k};
%         M = PDE_snap{2,k};
%         rho = S + M;
%         R = PDE_snap{3,k};
%         X = PDE_snap{4,k};
%         
%         agents_out = ABM_snap{1,k};
%         R_out = ABM_snap{2,k};
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(number);

fullwidth = 6.5;

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Manual') 
% Setting 'PaperPositionMode' to 'manual' unlinks 'Position' (what gets displayed) and 'PaperPosition' (what gets printed)
pos = [ 0 0  pos(3)*fullwidth/pos(3) pos(4)*fullwidth/pos(3)*1];
set(h,'PaperPosition',pos);
set(h,'Position',pos);
get(h,'Position');

titlesize = 20;
axesfontsize = 16;
mrk = 8;
mrkW = 2;
dotSize=30;
xSize=15;
xWid = 5;
snapClr = gold;

subplot(3,1,1)
    plot(   loglamsNew, peaks(1,:) ,'LineWidth', 3)
    hold on
    plot(   loglamsNew, peaks(2,:) ,'o', 'MarkerSize', mrk,'LineWidth',mrkW);
    plot(   loglams_snap, peaks_snap(1,:), '.','MarkerSize', dotSize,'LineWidth',mrkW,'Color',snapClr);
    plot(   loglams_snap, peaks_snap(2,:), 'x','MarkerSize', xSize,'LineWidth',xWid,'Color',snapClr);
    hold off
set(gca,'FontSize',ticklabelsize)
ylabel('Max Density P','FontSize',axislabelsize)
set(get(gca,'ylabel'),'rotation',90)

%     %remove whitespace
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset;
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];

subplot(3,1,2)
    plot(   loglamsNew, widths(1,:) ,'LineWidth', 3)
    hold on
    plot(   loglamsNew, widths(2,:) ,'o', 'MarkerSize', mrk,'LineWidth',mrkW)
    plot(   loglams_snap, widths_snap(1,:), '.','MarkerSize', dotSize,'LineWidth',mrkW,'Color',snapClr);
    plot(   loglams_snap, widths_snap(2,:), 'x','MarkerSize', xSize,'LineWidth',xWid,'Color',snapClr);
    hold off
set(gca,'FontSize',ticklabelsize)
ylabel('Std Width $W_\sigma$','FontSize',axislabelsize)
set(get(gca,'ylabel'),'rotation',90)

%     %remove whitespace
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset;
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];

subplot(3,1,3)
    plot(   loglamsNew, skews(1,:) ,'LineWidth', 3)
    hold on
    plot(   loglamsNew, skews(2,:) ,'o', 'MarkerSize', mrk,'LineWidth',mrkW)
    plot(   loglams_snap, skews_snap(1,:), '.','MarkerSize', dotSize,'LineWidth',mrkW,'Color',snapClr);
    plot(   loglams_snap, skews_snap(2,:), 'x','MarkerSize', xSize,'LineWidth',xWid,'Color',snapClr); 
    hold off
set(gca,'FontSize',ticklabelsize)
xlabel('$\log_{10}\lambda$','FontSize',axislabelsize)
ylabel('Skewness $\Sigma$','FontSize',axislabelsize)
set(get(gca,'ylabel'),'rotation',90)

%     %remove whitespace
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset;
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(h,[filename, '.eps'],'-depsc')
set(h,'PaperPositionMode','Auto')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Local Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sin, Min, Rin, varargout] = smallLam(L,nx,N,Rplus,v,lambda,alpha,beta,eta,theta,gamma,delta)
    [ksm, kms] = fxn_switching(Rplus, alpha, eta, gamma, beta, theta, delta);
    
    Sbar = kms/(kms+ksm);
    Mbar = ksm/(kms+ksm);
    varargout{1} = Sbar;
    varargout{2} = Mbar;
    
    % speed
    c = v*Mbar;
    % K'(Rplus)
    Kp = gamma*(eta-alpha)*exp(-gamma*Rplus)/c-delta*(theta-beta)*exp(-delta*Rplus)/(v-c);
    
    d = -2*v*ksm/(N*lambda*Rplus*Kp*kms);
    
    dx = L/(nx-1);
    X = linspace(-L/2,L/2,nx);
    % X = -L/2:dx:L/2;
    X = X';
    rho = N/(4*d)*(sech(X/(2*d))).^2;
    
    disp(['Mass ratio of first initial condition is ', num2str(sum(rho)*L/(nx-1)/N)])
    
    % initial locust populations
    Sin = Sbar*rho;
    Min = Mbar*rho;
    
    % initial resources
    Rin = ones(nx,1)*Rplus-(lambda*Rplus*kms/ksm*1/v)*N/2*(1-tanh(X/(2*d)));
    
end

function positions = sample_dist(mean_pulse, Xold, Xnew, N_agents)
% inputs:   mean_pulse = nx by 1 array = density of locusts (unscaled probability density function)
%           Xold = nx by 1 array = old spatial grid
%           Xnew = nx by 1 array = new spatial grid
% outputs:  positions = N by 1 array = list of locations of N locusts

    nx = length(Xold);

    % build a cumulative density distribution
    cumu_dist = @(k) sum(mean_pulse(1:floor(1*k)));
    vals = arrayfun(cumu_dist,1:nx);
    
    cumu_dist_new = interp1(Xold,vals,Xnew,'pchip',0);
    cumu_dist_new(Xnew>max(Xold)) = cumu_dist(nx);
    
    mass = cumu_dist(nx);
    
    cumu_vals = mass*linspace(1/(2*N_agents),(2*N_agents-1)/(2*N_agents),N_agents);
    
    positions = zeros(N_agents,1);
    for k = 1:N_agents
       %positions(k) = first index such that cumu_dist_new > cumu_vals(k)
       positions(k) = find(cumu_dist_new>cumu_vals(k),1);
    end
    
end
