% Creates three figures comparing the profiles of the ABM and PDE, for varying log(lambda).

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
axislabelsize = 12;%16;
ticklabelsize = 10;%12;
subfiglabelsize = 14;

%colors
ksm_color = 1/255*[237, 143, 33];   %gold
kms_color = 1/255*[111,73,179];     %purple
green =  [  0    0.5000      0];
gold = [0.9290, 0.6940, 0.1250];
orange = [0.85 0.325 0.098];
blue = [      0    0.4470    0.7410];
defwhite = [1 1 1];

%linewidth
lwidth = 1;

pcheck_plots = 0;
new_data = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 8 %%%
number = 8;
filename = sprintf('fig%d',number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% biological parameters
set_params


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHOOSE number of time steps %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The initial conditions come from runs out to nt = 200,000
% (This assumes we are loading 'data/cfShape_9:27.mat' on line 47)
nt = 800000; % probably want to go out to 800,000... that way total is nt = 10^6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set a range of log(lambda) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use a subset of that contained in 
%%% 'data/cfShape_9:27.mat'
loglams_snap = [-7.4, -6.3, -4.2];

n = length(loglams_snap);
% initialize storage variables
% top row for PDE, bottom row for ABM
widths_snap = zeros(2,n);
peaks_snap = zeros(2,n);
skews_snap = zeros(2,n);
PDE_snap = cell(4,n);
ABM_snap = cell(5,n);

% get data to be used for
% a) initial conditions from cfShape_DATE.mat
% b) making a big plot after loop
% DATE has format m:dd (month:dayday)
load('data/cfShape_9:27.mat')
% 'loglams' = 1x41 array of values of log10(lambda)
% 'peaks' = 2x41 array of values for maximum of the profile with PDE in top row, ABM in bottom
% 'widths' = 2x41 array of values for standard deviation of the profile with PDE in top row, ABM in bottom
% 'skews' = 2x41 array of values for skewness of the profile with PDE in top row, ABM in bottom
% PDE = 4 x 41 cell, rows described from top to bottom
    % S = nx x 1 array for density of stationary locusts at the last time step of the PDE
    % M = nx x 1 array for density of moving locusts at the last time step
    % R = nx x 1 array for density of resources at the last time step
    % X = nx x 1 array for spatial grid corresponding to the three variables above
% ABM = 5 x 41 cell
    % agents_out = N x 2 array of locust positions and states at the last time step of the ABM
    % R_out = T_steps_end x 1 array of resources remaining on the full spatial grid at the last time step
    % [optional]
    % mean_pulse = ?? x 1 array describing the locust density around the center of mass, time-averaged over the second half of the run
    % mean_resources = ?? x 1 array of the resource density around the center of mass of the locust pulse, time-averaged over the second half
    % xmesh = ?? x 1 array of spatial grid points corresponding to the two variables above
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%% --> Also includes computation because plotting 6 plots with a for loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initialize figure window %%%
h = figure(number);

fullwidth = 6.5;

titlesize = 24;
axesfontsize = 16;
lWidth = 2;
lWidABM = 0.5;
snapClr = gold;

set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Manual') 
% Setting 'PaperPositionMode' to 'manual' unlinks 'Position' (what gets displayed) and 'PaperPosition' (what gets printed)
pos = [ 0 0  pos(3)*fullwidth/pos(3) pos(4)*fullwidth/pos(3)*1.5];
set(h,'PaperPosition',pos);
set(h,'Position',pos);
get(h,'Position');


%sgtitle(['Model Output after $10^6$ time steps'],'FontSize',titlesize)

%%% loop through lambda values %%%
for k = 1:n
    
    lambda = 10^loglams_snap(k);
    
    % get index corresponding to loglam
    kk = find(loglams == loglams_snap(k));
    if length(kk) ~= 1
        disp('More than one index with matching loglam value!')
        pause
    end
    
%%% set PDE initial conditions
    Sin = PDE{1,kk};
    Min = PDE{2,kk};
    Rin = PDE{3,kk};
    X = PDE{4,kk};
    
%%% set ABM initial conditions
    agents_out = ABM{1,kk};
    R_out = ABM{2,kk};
    
    % we will SHIFT THEM TO THE BEGINNING OF A SPATIAL GRID!!!
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reconstruct Numerical Parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% for PDE
    % have X = spatial grid = linspace(-L/2, L/2, nx);
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
%%% call PDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if new_data == 1
    tic
    [R, S, M, c] = fxn_locustRK_shift(Sin,Min,Rin,L,nx,nt,N,Rplus,v,lambda,alpha,beta,eta,theta,gamma,delta,pcheck_plots);
    disp('time for PDE run was')
    toc

    S = S(:,end);
    M = M(:,end);
    R = R(:,end);
%%% record PDE data
    % record model outputs
    PDE_snap{1,k} = S;
    PDE_snap{2,k} = M;
    PDE_snap{3,k} = R;
    PDE_snap{4,k} = X;
    
    % compute peak
    rho = S+M;
    [m, max_idx] = max(rho);  % m = maximum, max_idx = index of maximimum (location of peak)
    % Interpolate to get a more accurate maximum
        ym = rho(max_idx-1);
        y0 = rho(max_idx);
        yp = rho(max_idx+1);
        interpmax = (16*y0^2-8*y0*(yp+ym)+(ym-yp)^2)/(16*y0-8*(ym+yp));
    peaks_snap(1,k) = interpmax;
    
    % compute standard deviation
    cm = sum(X.*(rho)/N)*dx; % center of mass
    width = sqrt(sum((X-cm).^2.*(rho)/N)*dx); %std of the profile
    widths_snap(1,k) = width;
    
    % compute skewness
    skews_snap(1,k) = -sum((X-cm).^3.*(rho)/N)*dx/(widths_snap(1,k))^3;
else
    load('data/cfShape_snap_9by27.mat')
    % loads PDE_snap, ABM_snap, peaks_snap, widths_snap, skews_snap
    S = PDE_snap{1,k};
    M = PDE_snap{2,k};
    rho = S + M;
    R = PDE_snap{3,k};
    X = PDE_snap{4,k};
end % if new_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% call ABM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if new_data == 1
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

    tic
    % optional outputs: mean_pulse, mean_resources, xmesh
    %(must change transient detection threshold to obtain)
    [agents_out, R_out, shape_stats] =...
         fxn_ABM_shift(agents_in,R_in,T_steps_end,delta_t,N,Rplus,v,lambda,alpha,beta,eta,theta,gamma,delta,pcheck_plots);
    disp('time for ABM run was')
    toc   

%%% record ABM data
    % record model outputs
    ABM_snap{1,k} = agents_out;
    ABM_snap{2,k} = R_out;
    
    % record shape statistics
    peaks_snap(2,k) = mean(shape_stats(ceil(nt/2):end,1));
    widths_snap(2,k) = mean(shape_stats(ceil(nt/2):end,2));
    skews_snap(2,k) = mean(shape_stats(ceil(nt/2):end,3));
    
    %%% Save Data
    save('data/cfShape_snap.mat','peaks_snap','widths_snap','skews_snap','loglams_snap','PDE_snap','ABM_snap')

else
    % load('data/cfSHape_snap.mat') in earlier if statement under "call PDE"
    agents_out = ABM_snap{1,k};
    R_out = ABM_snap{2,k};
    % already loaded peaks_snap, widths_snap, skews_snap
end %if new_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For odd reasons, plot ABM first!
    
    ytix = 1900;
    
    %%% Plot ABM %%%
    figure(number)
    subplot(3,2,2*k)
    set(gcf,'Units','Inches');
    pos = get(gcf, 'Position'); % [0.1300    0.5838    0.3347    0.3412] from full figure [0 0 6.5 5];
    set(gcf,'Position', [pos(1),pos(2), pos(3),pos(4)])
    
    CM_idx = round(mean(agents_out(:,1)));
    
    %%% plot the location of the locusts
    [cnts, edges] = histcounts(agents_out(:,1),'BinMethod','Integers');
    x_idx = floor(edges(2:end)); % spatial indices
    buff = 100; %100 %index buffer
    % add buffer onto left and right of spatial indices 
    x_ind = [fliplr(linspace(x_idx(1)-1,x_idx(1)-buff,buff)), x_idx, linspace(x_idx(end)+1,x_idx(end)+buff,buff)];
    x_grid = (x_ind-CM_idx)*delta_x; % make spatial indices into spatial grid, with CM at x = 0
    cnts = [zeros(1,buff), cnts, zeros(1,buff)]; % extend cnts by buffer 
    
    left_dist = CM_idx-min(x_ind);
    right_dist = max(x_ind)-CM_idx;
    
    yheight = 2000;% 1.1*max( max(rho), max(cnts)/delta_x);
    plot(x_grid,cnts/delta_x,... %cnts/delta_x makes this to density
                     ... %cnts makes this locusts per gridpoint
        'Color',orange,'LineWidth',lWidABM) 
    axis([x_grid(1) x_grid(end) 0 yheight])
    yticktop = floor(max(cnts)/delta_x/100)*100;
%     if k == 2
%         yticks([0 yticktop])
%     else
%         yticks([0 yticktop ytix])
%     end
    yticks([])
    %ylabel('Locusts per gridpoint','Color','k')
    %xlabel('Space, centered at CM')
    
    %%% add resource curve corresponding to histogram
    yyaxis right
    %resources(min(agents(:,1))-max_idx+mode_pt:max(agents(:,1))-max_idx+mode_pt),...
    plot(x_grid,R_out(x_ind),...
        ':','Color',green,'LineWidth',lWidth);
    ylim([0,Rplus*1.1])
    yticks([0 200])
    set(gca,'ycolor',green)
    %ylabel('Resource Density','Color','k') %fontsize doesn't stick for some reason...
    set(gca,'FontSize',ticklabelsize)
    if k == 3
        xlabel('Distance from Mean','Fontsize',axislabelsize)
    end
    ylabel('Resource Density','Fontsize',axislabelsize)
    %title(['$\log_{10} \lambda = $ ', num2str(loglams_snap(k))],'Fontsize',titlesize)
    
%     %remove whitespace
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset;
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];
    
%%% Plot PDE %%%
    figure(number)
    subplot(3,2,2*k-1)
    set(gcf,'Units','Inches');
    pos = get(gcf, 'Position'); % [0.1300    0.5838    0.3347    0.3412] from full figure [0 0 6.5 5];
    set(gcf,'Position', [pos(1),pos(2), pos(3), pos(4)])
    
    CM = sum(X.*rho)*dx/N;
    CMidx = find( X > CM,1);
    first_idx = CMidx-left_dist;
    last_idx = CMidx+right_dist;
    
    % in case we have made it outside of where rho, R are defined
    if first_idx < 1
        % disp('first_idx < 1') 
        % pause %used to check
        rho = [zeros(1-first_idx,1); rho];
        R = [R(1)*ones(1-first_idx,1); R];
        last_idx = last_idx+1-first_idx;
        first_idx = 1;
    end
    
    if last_idx > nx
        pause
        disp('last_idx > nx')
        rho = [rho; zeros(last_idx-nx+1,1)];
        R = [R; Rplus*ones(last_idx-nx+1,1)];
        last_idx = length(rho);
    end
    
    plot(   ...X, M,...
            ...X, S,...
            x_grid, rho(first_idx:last_idx),...
            'Color',blue,'LineWidth',lWidth);
    axis([x_grid(1) x_grid(end) 0 yheight])
    if k == 2
        yticks([0 yticktop])
    else
        yticks([0 yticktop ytix])
    end
    
    set(gca,'FontSize',ticklabelsize)
    
    if k == 3
        xlabel('Distance from Mean','Fontsize',axislabelsize)
    end
    ylabel('Locust Density','Fontsize',axislabelsize)
    %ylabel('Locust Density')
    %xlabel('Displacement from CM')
    
    % add resource curve
    yyaxis right
    plot(   x_grid, R(first_idx:last_idx),...
        'Color',green,'LineWidth',lWidth);
    ylim([0,Rplus*1.1])
    yticks([])%yticks([0 200])
    set(gca,'ycolor',green)
    %ylabel('Resource Density','Color','k')
    
%     title(['$\log_{10} \lambda = $ ', num2str(loglams_snap(k))],'Fontsize',axesfontsize)
    
    
%     %remove whitespace
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset;
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];
    
end % for k=1:n

% subfigure labels
str = '$\log_{10} \lambda = -7.4$ ';
annotation('textbox', [.4, .975, 0, 0], 'string', str, 'FontSize', subfiglabelsize,'Interpreter','latex')

str = '$\log_{10} \lambda = -6.3$ ';
annotation('textbox', [.4, .675, 0, 0], 'string', str, 'FontSize', subfiglabelsize,'Interpreter','latex')

str = '$\log_{10} \lambda = -4.2$ ';
annotation('textbox', [.4, .375, 0, 0], 'string', str, 'FontSize', subfiglabelsize,'Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(h,[filename, '.eps'],'-depsc')
set(h,'PaperPositionMode','Auto')

