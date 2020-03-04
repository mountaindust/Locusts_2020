%%%%%%%%%%%%%%
% Run this script for a single ABM realization
% This is supposed to be a cleaned upversion
%%%%%%%%%%%%%%
% varargout = {xmesh_out, mean_pulse_out, mean_resources_out} add mean_resources_out later?
function [agents_out, R_out, shape_stats, varargout] =...
         fxn_ABM_ave(agents_in,R_in,T_steps_end,delta_t,N_agents,R_plus,v,lambda,alpha,beta,eta,theta,gamma,delta,pcheck_plots)

% set(0, 'defaulttextinterpreter', 'latex')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  

%for plotting
green =  [  0    0.5000      0];

%
%  Set up flags
%  movieflag 1=movie, 0=nomovie to save movie
%  savePlots will save .eps of figures
%  PLOT_FLAG Real-time plotting, yes=1, no=0 (will plot last time-step)
%  PLOT_SUMMARY %plot a summary pulse shape with fitted functions after simulation

    movieflag=0;
    savePlots = 0;
    PLOT_FLAG = 0; %plot each integer time as it is calculated
    PLOT_SUMMARY = 0; %

% Number of locusts that we will track as individuals
N_sav=5;

% Initialize a spatial grid with delta_x based on delta_t
% spatial domain length set by input (initial condition) R_in
delta_x = v*delta_t;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Probability function parameters               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p_SM = b - (b-a)*e^{-r_sm*R},  
% p_SM = d - (d-c)*e^{-r_ms*R}

% the greeks
a = alpha;
b = eta;
r_sm = gamma;
c = beta;
d = theta;
r_ms = delta;

% Probability function of going from S -> M per delta_t
%p_SM = @(R) (a-b.*R).*delta_t;
p_SM = @(R) (b - (b-a).*exp(-r_sm.*R)).*delta_t;

% Probability function of going from M -> S per delta_t
%p_MS = @(R) (c+d.*R).*delta_t;
p_MS = @(R) (d - (d-c).*exp(-r_ms.*R)).*delta_t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Resource function                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exp function of S to reduce R by. R(t+1) = R*R_exp(S)
% lambda = above %rate of resourse proportion decrease per stationary agent 
               %per unit time in a given location
R_exp = @(S) exp(-lambda*S*delta_t/delta_x); % needs to depend on delta_t and delta_x!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Plotting and data gathering parameters           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% length of time window after which we begin gathering data
% choose this LARGE enough so that the profile is "near" equilibrium at this value
% choose this SMALL enough so that we average over a "large" number of time steps after this value
t_delay_window_len = round(.5*T_steps_end); % half the full number of time steps
T_steps_ave = T_steps_end-t_delay_window_len;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Simulation initialization %%%%%%%%%%
% agents = N_agents x 2 array
% Columns are x_ind, S/M: Stationary = 0, Moving = 1
agents = agents_in;

% resources = length(R_in) x 1 array
% typically length(R_in) = T_steps_end 
    % an alternative might be < T_steps_end*c/v, or similar (ensures no locust runs off)
% units: resources per 1 meter
resources = R_in;

% final time
% units: seconds
T_end = T_steps_end*delta_t;

%%%%%%%%%% Initialize Stat Arrays %%%%%%%%%%

% initialize end-of-transient values
start_t_idx = 0;
start_x_idx = 0;

% Get an array in which to store spatial mean, median, and max
mean_ary = zeros(T_steps_end,1); % center of mass (up to factor delta_x)
% median_ary = zeros(T_steps_end,1); 
max_ary = zeros(T_steps_end,1); % peak value (up to factor delta_x)

% Get an array in which to store spatial variance
std_ary = zeros(T_steps_end,1); % 'width' = standard devation (up to factor delta_x)

% Get an array in which to store the spatial skewness
skew_ary = zeros(T_steps_end,1); % skewness needs no factor, is dimensionless

% Get an array in which to store total number of moving locusts
M_cnt = zeros(T_steps_end,1);

% Get an cell array in which to store the positions of all locusts
posn = cell(T_steps_ave,2);
% % for alternative
% posn_sum = zeros(7000,1);

% Get a cell array in which to store the resources on some intervals of varying sizes
res = cell(T_steps_ave,1);

% Get array in which to store the positions of several agents for
% examination
sav_idx = [1:N_sav]; %number of agents
sav_loc = zeros(T_steps_end,length(sav_idx));

pltcnt = 0;

%%%%%%%%%% Main Loop %%%%%%%%%%

for tstep = 1:T_steps_end
    
    if mod(tstep,10000) == 0
        disp(['Timestep ', num2str(tstep), ' of ', num2str(T_steps_end)])
    end
    
    %
    %%% change state with some probability %%%
    %
    r_nums = rand(N_agents,1); %get random values for each agent
    
    % Select moving agents, apply transitions with some probability
    % MS_agents = moving agent AND randomly within prob(R(x))
    moving_agents = (agents(:,2) == 1); %record moving agents - will overwrite!
    MS_agents = (moving_agents).*(r_nums <= p_MS(resources(agents(:,1))));
    % If so, subract 1 to make stationary (0)
    agents(:,2) = agents(:,2) - MS_agents;
    
    % Select stationary agents, apply transitions with some probability
    % SM_agents = stationary agent AND randomly within prob(R(x))
    stationary_agents = ~moving_agents;
    SM_agents = (stationary_agents).*(r_nums <= p_SM(resources(agents(:,1))));
    % If so, add 1 to make moving
    agents(:,2) = agents(:,2) + SM_agents;
    
    %
    %%% move the moving locusts by an integer (one), which represents v*delta_t=delta_x %%%
    %
    agents(:,1) = agents(:,1) + agents(:,2);

    %
    %%% Reduce R by an amount based on number of S locally %%%
    %
    if min(agents(:,2)) < 1 %only do this if some are stationary
        % get x-values of stationary agents
        stat_agent_loc = agents(~logical(agents(:,2)),1);
        if max(stat_agent_loc)-min(stat_agent_loc) >= 2^16
            disp('Stationary agents are too spread out to be tabulated using histcounts, "Integers"')
            pause
        end
        [cnts, edges] = histcounts(stat_agent_loc, 'BinMethod', 'Integers');
        values = floor(edges(2:end));
        tbl = [values; cnts]; %values (first row) and counts (second row)
        % change amount of resources at these x values by an amount calculated
        % using the counts, multiplied by time and rate of eating
        resources(tbl(1,:)) = resources(tbl(1,:)).*R_exp(tbl(2,:))';
    end
    % else they are all moving!
    
    %%%%%%%%%%%%%%%%%%%
    %%% Record Data %%%  
    %%%%%%%%%%%%%%%%%%%
    
    % Sanity check
    assert(max(agents(:,2))<=1 && min(agents(:,2))>=0)
    % Record number moving
    M_cnt(tstep) = sum(agents(:,2)); %number moving
    
    % Record position of chosen few
    sav_loc(tstep,:) = agents(sav_idx,1)*delta_x;
    
    % Record Positions and mean position (center of mass)
    % and resources wherever locusts are present
    if tstep > t_delay_window_len
        posn_ind = agents(:,1);
        CM_ind = floor(mean(agents(:,1)));
        posn{tstep-t_delay_window_len,1} = posn_ind;
        posn{tstep-t_delay_window_len,2} = CM_ind;
        
%         % for alternative method
%         posn_sum = posn_sum+posn_ind;

        % Record resources
        res{tstep-t_delay_window_len} = resources(min(agents(:,1)):max(agents(:,1)));
    end
    
    % Basic statistics
    mean_ary(tstep) = mean(agents(:,1)*delta_x); % center of mass
    %%% Width
    std_ary(tstep) = std(agents(:,1)*delta_x);
    
    %%% Peak
    [M, F] = mode(agents(:,1));
    max_ary(tstep) = F/delta_x ;
    
    %%% Skewness
    skew_ary(tstep) = -skewness(agents(:,1));  %dimensionless, so no need for delta_x
    % -1/N_agents*sum((agents(:,1)*delta_x-mean_ary(tstep)).^3)/(std_ary(tstep))^3;
    
end

% Construct the mean profile from posn and the mean resources from res
% tic
Xlist_ind = zeros(2,T_steps_ave);
profiles = cell(1,T_steps_ave);

    for k = 1:T_steps_ave
       %%% make histograms
       posn_shifted = posn{k,1}-posn{k,2};
       [counts, edges] = histcounts(posn_shifted,'BinMethod','Integers');
       Xlist_ind(1:2,k) = [floor(edges(2)); floor(edges(end))];
       profiles{k} = counts;
    end

xmin_ind = min(Xlist_ind(1,:));
xmax_ind = max(Xlist_ind(2,:));
xlen = xmax_ind - xmin_ind + 1;
profiles_sum = zeros(1,xlen);
profiles_num = zeros(1,xlen);
% resources
res_sum = zeros(xlen,1);
res_num = zeros(xlen,1);

    for k = 1:T_steps_ave
       ini = Xlist_ind(1, k) - xmin_ind + 1;
       iniR =  Xlist_ind(1, k) + 1;
       fin = Xlist_ind(2, k) - xmin_ind + 1;
       finR = Xlist_ind(2, k) + 1;
       profiles_sum(ini:fin) = profiles_sum(ini:fin) + profiles{k};
       profiles_num(ini:fin) = profiles_num(ini:fin) + 1;
       
       % resources
       res_sum(ini:fin) = res_sum(ini:fin)+res{k};
       res_num(ini:fin) = res_num(ini:fin) + 1;
    end

%mean_pulse = (profiles_sum ./ profiles_num)/delta_x; % don't average over missing data
mean_pulse = (profiles_sum ./ T_steps_ave)/delta_x; % assume missing data = 0

% resources
mean_resources = (res_sum ./ res_num);

% disp('Constructing the mean profile took ')
% toc

% construct an xgrid!
xmesh = xmin_ind:1:xmax_ind;
xmesh = xmesh*delta_x;

% % for alternative
% posn_ave = posn_sum/T_steps_ave;
% [counts, edges] = histcounts(posn_ave,'BinLimits',[xmesh(1)-delta_x/2,xmesh(end)+delta_x/2],'BinWidth',delta_x);
% xmesh2 = floor(edges(2:end))*delta_x;
% mean_pulse2 = counts/delta_x;


%%% Set output variables
agents_out = agents;
R_out = resources;

shape_stats = [max_ary, std_ary, skew_ary];
varargout{1} = xmesh;
varargout{2} = mean_pulse;
varargout{3} = mean_resources;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time varying stats %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get a time and space grid, leave off first time point (velocity)
tgrid = (2:T_steps_end)*delta_t;
xgrid = (1:length(resources))*delta_x;

% peak saved in max_ary
% width saved in std_ary
% skewness saved in skew_ary

% velocity of the mean
vel_mean = (mean_ary(2:end) - mean_ary(1:end-1))/delta_t;
% this is equal to diff(mean_ary)/delta_t

if pcheck_plots
    axesfontsize = 15;
    labelsize = 20;

    %%%%%%%%%%%%%%%%%%%%%%%
    % Plot Stats - Figure 3
    %%%%%%%%%%%%%%%%%%%%%%%
    h = figure(3);

    pos = get(h, 'Position');
    set(h, 'Position', [round(pos(1)/2),round(pos(1)/2),pos(3)*1,pos(4)*1.5]);
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    subplot(5,1,1)
    plot(tgrid,vel_mean,'LineWidth',2);
    title('Velocity of mean position over time','FontSize',14)
    xlabel('time (s)')
    ylim([0,2*max(vel_mean)])
    'average velocity is'
    mean(vel_mean(floor(0.75*length(vel_mean)):end))

    subplot(5,1,2)
    plot([0 tgrid],max_ary,'LineWidth',2);
    title('Maximum density over time','FontSize',14)
    xlabel('time (s)')
    ylim([0,1.25*max(max_ary)])
    'the final peak density is'
    max_ary(end)

    subplot(5,1,3)
    plot([0 tgrid],std_ary,'LineWidth',2);
    title('Std of position over time','FontSize',14)
    xlabel('time (s)')
    ylim([0,1.25*max(std_ary)])
    
    subplot(5,1,4)
    plot([0 tgrid],skew_ary,'LineWidth',2);
    title('Skewness of profile over time','FontSize',14)
    xlabel('time (s)')
    ylim([min(1.25*min(skew_ary),-.01),max(1.25*max(skew_ary),.01)])
    
    subplot(5,1,5)
    plot(xgrid,R_in,'LineWidth',2,'Color',green) %initial resources
    hold on
    plot(xgrid,R_out,':','LineWidth',2) %final resources
    hold off
    xlim([0,max(agents_out(:,1))*delta_x])
    ylim([min(R_out),R_plus*1.01]);
    title('Resources: initial (green) and final (dotted)','FontSize',14)
    
    if savePlots == 1
        saveas(gcf,'ABM_stats.eps','epsc')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Locust Profiles - Figure 4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = figure(4);

    pos = get(h, 'Position');
    set(h, 'Position', [round(pos(1)/2),round(pos(1)/2),pos(3)*1,pos(4)*1.5]);
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    subplot(2,1,1)
    % plot the mean locust density
    plot(xmesh,mean_pulse)
%     % for alternative
%     hold on
%     plot(xmesh2,mean_pulse2, 'Color', 'k')
%     hold off
    ylabel('Locust Density','Color','k')
    xlabel('Space (not spatial grid, YES shifted)')
    title(['averaged over ' num2str(T_steps_ave) 'time steps'])
    
    % add a mean resource curve
    yyaxis right
    plot(xmesh,mean_resources,'Color',green,'LineWidth',3);
    ylim([min(R_out),R_plus*1.01])
    set(gca,'ycolor',green)
    ylabel('Resource Density','Color','k')
    
    subplot(2,1,2)
    % plot histogram of the final location of the locusts
    [cnts, edges] = histcounts(agents_out(:,1),'BinMethod','Integers');
    x_ind = floor(edges(2:end));
    x_grid = x_ind*delta_x;
    plot(x_grid,cnts/delta_x)
    xlim([x_grid(1), x_grid(end)])
    ylabel('Locust Density','Color','k')
    xlabel('Space (not spatial grid, NOT shifted)')
    
    %add resource curve corresponding to histogram
    yyaxis right
    %resources(min(agents(:,1))-max_idx+mode_pt:max(agents(:,1))-max_idx+mode_pt),...
    plot(x_grid,R_out(x_ind),...
        ':','Color',green,'LineWidth',3);
    ylim([min(R_out),R_plus*1.01])
    set(gca,'ycolor',green)
    ylabel('Resource Density','Color','k') %fontsize doesn't stick for some reason...
    
end