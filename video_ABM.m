%%%%%%%%%%%%%%
% Run this script for a single ABM realization and create a video
%%%%%%%%%%%%%%

close all
clear all

addpath('functions','data')

% set(0, 'defaulttextinterpreter', 'latex')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose how often to plot and record a frame for video

% plots visualizations of ABM every 'plot_frames' timesteps
% choose plot_frames = 5 for a longer video with higher time-resolution
% choose plot_frames = 25 for a faster video
plot_frames = 5;


%%% for plotting
% fonts
titlesize = 16;
axislabelsize = 14;
ticklabelsize = 12;
subfiglabelsize = 12;

% colors
ksm_color = 1/255*[237, 143, 33];   %gold
kms_color = 1/255*[111,73,179];     %purple
green =  [  0    0.5000      0];
gold = [0.9290, 0.6940, 0.1250];
orange = [0.85 0.325 0.098];
blue = [      0    0.4470    0.7410];
maroon = [0.6350    0.0780    0.1840];

% linewidth
lwidth = 1;


%%% biological parameters
set_params
N_agents = N;
R_plus = Rplus;

% % Number of locusts that we will track as individuals
N_sav=5;

% How long to run
T_steps_end = 7500;
% Time step size. Do not increase above 1 or things will break!!
% Also, don't make too small or things will get really noisy...
delta_t = 1;

% Initialize a spatial domain width and delta_x based on delta_t
delta_x = v*delta_t;
domain_length = ceil(50/delta_x); %ceil(T_steps_end/5);
X = delta_x*(1:1:domain_length);

% pause and check plots
pcheck_plots = 0;

% set a seed for the random number generator 
% so that we don't accidentally run off the right hand edge.
rng(2)

% Initial Conditions
agents_in = zeros(N_agents,2); % Columns are x_ind, S/M
% Stationary = 0, Moving = 1

%%% start all locusts in the first 5 spatial points, equally spaced, stationary
% Ndiv5 = floor(N_agents/5);
% for n=1:5
%     agents_in((Ndiv5*(n-1)+1):Ndiv5*n,1) = n; %set spatial locations
% end
% if Ndiv5*5 ~= N_agents %check for remaining locusts
%     agents_in((Ndiv5*5+1):end) = 1; %put them in the first grid point
% end
% agents_in(:,2) = 0; %set all agents to stationary

%%% start locusts arranged in a Gaussian distribution
cen = 5; % X(end)/10; <-- for changing domains
sig = 1;
cen_idx = cen/delta_x;
sig_idx = sig/delta_x;
agents_in(:,1) = ceil(normrnd(cen_idx,sig_idx,N_agents,1));
agents_in(:,2) = 0; %set all agents to stationary

R_in = ones(domain_length,1).*R_plus; %units: resources per delta_x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 101 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number = 101;
filename = sprintf('fig%d',number);

    xaxis = 'Space $x_n$ (m)';
    yaxis = 'Locusts per Gridpoint'; %$S_{n,m}+M_{n,m}$
    %yaxis18 = 'Locusts per meter (1/m)';
    ryaxis = 'Resource Density (g/m)';
    
    axesfontsize = ticklabelsize;
    labelsize = axislabelsize;
    
    lyt = [0 400 800 1200 1600]; % left y-axis ticks
    lyt_stoch = [0 16 32 48 64];

    resColor = green;
    resWidth = 2;
    
    abm_width = 1;
    smooth_width = 2;
    exp_width = 2;
    
% setting up figure
h = figure(number);
dims = [ 0 0 13 6.5];
set(h,'Units','Inches');
% pos = get(h,'Position')
set(h,'PaperPositionMode','Manual') % Setting this to 'manual' unlinks 'Position' (what gets displayed) and 'PaperPosition' (what gets printed)
set(h,'PaperPosition', dims);
set(h,'Position', dims);
get(h,'Position');

hgt = 1;
% ABM plot
h2 = subplot(2,1,2);
pos = [0.1300    0.1100    0.7750    0.3412]; %get(h2, 'Position');
set(h2,'Position', [pos(1),pos(2)*1, pos(3), pos(4)*hgt])
    
        yl = [-4 75];
        ylim(yl);
        xl = xlim;
        hold on
    %     hline = plot(xl,[0,0],'k');
    %     vline = plot([0,0],yl,'k--');
    %     uistack(hline,'bottom')
    %     uistack(vline,'bottom')

        yticks(lyt_stoch)
        set(gca,'FontSize',ticklabelsize)

        xlabel(xaxis,'FontSize',axislabelsize)
        set(gca,'ycolor',orange)
        ylabel(yaxis,'FontSize',axislabelsize,'Color', orange)
        %title(['Histogram of Locusts at Time $t_m = $ ' num2str(0) ' sec'],...
        %    'FontSize',titlesize)

        yyaxis right
        ylim([(yl(1)/yl(2))*R_plus*1.1,R_plus*1.1])
        set(gca,'ycolor',resColor)
        set(gca,'FontSize',ticklabelsize)
        ylabel(ryaxis,'FontSize',axislabelsize,'Color',green)


        hold off
        %remove whitespace
        ax = gca;
        outerpos = ax.OuterPosition;
        ti = ax.TightInset; 
        left = outerpos(1) + ti(1);
        bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3);
        ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];

% track schematic
h1 = subplot(2,1,1);
pos = get(h1, 'Position');
set(h1,'Position', [pos(1),pos(2), pos(3), pos(4)*1]) %adjust this later to make top plot thinner

axis([X(1) X(end) 0 N_agents])
        set(gca,'FontSize',ticklabelsize)
        
        ylabel('Spatial Cross Section (1 m)','Fontsize',axislabelsize)
        xlabel(xaxis,'FontSize',axislabelsize)
        yticks([])
        title(['Visualizations of the Agent-Based Model, Time $t_m = $ ' num2str(0) ' sec'],...
            'FontSize',titlesize)
        
            %remove whitespace
            ax = gca;
            outerpos = ax.OuterPosition;
            %ti = ax.TightInset; 
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Begin Main Script %%% (from fxn_ABM_shift.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Set up flags
%  movieflag 1=movie, 0=nomovie to save movie
%  savePlots will save .eps of figures
%  PLOT_FLAG Real-time plotting, yes=1, no=0 (will plot last time-step)
%  PLOT_SUMMARY %plot a summary pulse shape with fitted functions after simulation

    movieflag = 1;
    savePlots = 0;
    PLOT_FLAG = 0; %plot each integer time as it is calculated
    PLOT_SUMMARY = 0; %

%
% Set up video
% 
  if movieflag
      vid = VideoWriter('video','MPEG-4');
      vid.FrameRate = 20;
      vid.Quality = 75;
      open(vid);
  end    
    
% Number of locusts that we will track as individuals
N_sav=5;

% Initialize a spatial grid with delta_x based on delta_t
% spatial domain length set by input (initial condition) R_in
delta_x = v*delta_t;
X = delta_x*(1:1:domain_length);


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
    
    if mod(tstep,plot_frames) == 0
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Construct Histogram %%%  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        [counts, edges] = histcounts(agents(:,1),'BinMethod','Integers');
        min_ind = min(agents(:,1));
        max_ind = max(agents(:,1));
        histo = zeros(1,domain_length);
        histo(min_ind:max_ind) = counts;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Plot! %%%  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(number)
        
        % ABM plot
        subplot(2,1,2);
%         set(gcf,'Units','Inches');
%         pos = [0.1300    0.1100    0.7750    0.3412]; %get(h2, 'Position');
%         set(gcf,'Position', [pos(1),pos(2)*1, pos(3), pos(4)*hgt])

        histplot = plot(X,histo,'LineWidth',abm_width,'Color',orange);
        % uistack(histplot, 'bottom')
        % yl defined above
        ylim(yl);
        xl = xlim;
        hold on
    %     hline = plot(xl,[0,0],'k');
    %     vline = plot([0,0],yl,'k--');
    %     uistack(hline,'bottom')
    %     uistack(vline,'bottom')

        yticks(lyt_stoch)
        set(gca,'FontSize',ticklabelsize)

        xlabel(xaxis,'FontSize',axislabelsize)
        set(gca,'ycolor',orange)
        ylabel(yaxis,'FontSize',axislabelsize,'Color', orange)
        %title(['Time $t_m = $ ' num2str(tstep) ' sec'],...
        %    'FontSize',titlesize)

        yyaxis right
        resplot = plot(X, resources, ':', 'Color', resColor,'LineWidth',resWidth);
        ylim([(yl(1)/yl(2))*R_plus*1.1,R_plus*1.1])
        set(gca,'ycolor',resColor)
        %set(gca,'FontSize',ticklabelsize)
        ylabel(ryaxis,'FontSize',axislabelsize,'Color',green)


        hold off
            %remove whitespace
            ax = gca;
            outerpos = ax.OuterPosition;
            ti = ax.TightInset; 
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti(2) - ti(4);
            ax.Position = [left bottom ax_width ax_height];
        
        % track schematic
        subplot(2,1,1);
        
        locSize = 4;
        
        agent_lane = [agents(:,1)*delta_x agents(:,1) (1:N_agents)'];
        agent_slane = [agent_lane(~logical(agents(:,2)),1) agent_lane(~logical(agents(:,2)),3)];
        agent_mlane = [agent_lane(logical(agents(:,2)),1) agent_lane(logical(agents(:,2)),3)];
        
        plot(agent_slane(:,1), agent_slane(:,2) ,'.','MarkerSize', 6,'Color',maroon,'DisplayName','stationary, feeding locust');
        hold on
        plot(agent_mlane(:,1), agent_mlane(:,2) ,'.','MarkerSize', 6,'Color',blue,'DisplayName','moving locust')
        legend({},'AutoUpdate','off','location','northeast')
        axis([X(1) X(end) 0 N_agents])
        set(gca,'FontSize',ticklabelsize)
        hold off
        
        ylabel('Spatial Cross Section (1 m)','Fontsize',axislabelsize)
        xlabel(xaxis,'FontSize',axislabelsize)
        yticks([])
        title(['Visualizations of the Agent-Based Model, Time $t_m = $ ' num2str(tstep) ' sec'],...
            'FontSize',titlesize)
        
            %remove whitespace
            ax = gca;
            outerpos = ax.OuterPosition;
            ti2 = ax.TightInset; 
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            ax_width = outerpos(3) - ti(1) - ti(3);
            ax_height = outerpos(4) - ti2(2) - ti2(4);
            ax.Position = [left bottom ax_width ax_height];
            
            
            
        if movieflag
            frame = getframe(h);
            writeVideo(vid,frame);
        end
    end
      
end

if movieflag
        close(vid);
end 
