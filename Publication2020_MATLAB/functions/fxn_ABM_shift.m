%%%%%%%%%%%%%%
% ABM time/space refinement: Stationary/Moving, resources, discrete grid
%
% Run this script for a single ABM realization w/ statistical analysis
%%%%%%%%%%%%%%
%cspeed, cstd, R_minus, Rstd -- old outputs
% agents_out, R_out, shape_stats,  -- old outputs
% varargout = {mean_pulse_out, mean_resources_out, xmesh_out}
function [fig3outs, fig4outs, fig6outs] =...
         fxn_ABM_shift(agents_in,R_in,T_steps_end,delta_t,N_agents,R_plus,v,lambda,alpha,beta,eta,theta,gamma,delta,pcheck_plots)

% set(0, 'defaulttextinterpreter', 'latex')
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
% close all
% clear all

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
%
% Set up video
% 
  if movieflag
      vid = VideoWriter('abmjwb','MPEG-4');
      vid.FrameRate = 20;
      vid.Quality = 75;
      open(vid);
  end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                     Input Parameters                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N_agents = 7000;%7294.540405;       %total mass
% R_plus = 200;%194.487762;     %initial resources
% v = 0.04;%0.046614;           %speed of a moving locust (m/s)
% lambda = 10^(-5);%10^(-5.002960);    %feeding rate of a stationary locust
% % Sobol version of Switching Rate parameters
%     % M -> S low resource
%     beta = 0.02;%0.015695;
%     % S -> M rate: ratio of large R/small R
%     etaalpha = 0.8;%0.759527; % < 1
%     % M -> S rate: ratio of large R/small R
%     thetabeta = 7;%7.215061; % > 1
%     % alpha/beta - eta/theta = small R ratio SM/MS - large R ratio SM/MS
%     DELTA = 10^(-0.7);%10^(-0.653073);
%     % exponential rates
%     gamma = 0.03;%0.034318;     %S->M
%     delta = 0.005;%0.004681;     %M->S
% 
% %%% call a function to convert these to actual greek parameter values
% [alpha, beta, eta, theta] = fxn_paramRats(beta,thetabeta,etaalpha,DELTA);


% Number of locusts that we will track as individuals
N_sav=5;

% Initialize a spatial domain with width and delta_x based on delta_t
% domain_length = ceil(T_steps_end*(.007/v));
% spatial domain length set by initial condition R_init
delta_x = v*delta_t;

% Initial resources. Assign a scalar or an array of length domain_length
% Domain with uniform noise around 1
%R_init = ones(domain_length,1) + rand(domain_length,1).*1.6 - 0.8;
% Domain with oscillations
%R_init = 1 + 0.8*sin((1:domain_length)'.*delta_x./20.*2.*pi);

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

% length of time window for auto detection of end of transients using std
% of wave-speed velocity
t_delay_window_len = round(.5*T_steps_end);
% threshold std of velocity under which we are satisfied the transients
% have passed
vel_std_thrs = 1; %0.5*10^(-1); make this window large so that we use just the second half of the run as arranged by t_delay_window_len
%%% (From some older version) We WANT to never detect the end of transients, so that we never compute mean_pulse, mean_resources

%%% t_delay_window_len and vel_std_thrs must be chosen by examining mean
%%% velocity over the course of a full run. Choose window length long
%%% enough to detect slow increases. Choose threshold large enough to allow
%%% small stochasticity. Balance these!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Simulation initialization %%%%%%%%%%
% agents = zeros(N_agents,2); % Columns are x_ind, S/M
% Stationary = 0, Moving = 1
agents = agents_in;

% % Create array of resources %%% THIS BREAKS IF NON-CONSTANT resources
% resources = ones(domain_length,1).*R_plus; %units: resources per delta_x
resources = R_in;

% number of time steps that need to be run
% T_steps_end = ceil(T_end/delta_t);
T_end = T_steps_end*delta_t;

%%%%%%%%%% Initialize Stat Arrays %%%%%%%%%%

% array for transient detection
t_vel_window = zeros(t_delay_window_len,1);
% initialize end-of-transient values
start_t_idx = 0;
start_x_idx = 0;

% Get an array in which to store spatial mean, median, and max
mean_ary = zeros(T_steps_end,1);
median_ary = zeros(T_steps_end,1);
max_ary = zeros(T_steps_end,1);
% Get an array in which to store spatial variance
std_ary = zeros(T_steps_end,1);
% Get an array in which to store the spatial skewness
skew_ary = zeros(T_steps_end,1);
% Get an array in which to store total number of stationary and moving
M_cnt = zeros(T_steps_end,1);
% Get an array in which to store exp fit params
B_ary = zeros(2,T_steps_end);
% Get array in which to store root mean square error
rmserr = zeros(T_steps_end,1);
% Get array in which to store the positions of several agents for
% examination
sav_idx = [1:N_sav]; %number of agents
sav_loc = zeros(T_steps_end,length(sav_idx));

pltcnt = 0;

%%%%%%%%%%%
%%% Code for Figure 1 %%
%%%%%%%%%%%%

%figure(1)
pos = get(gcf,'OuterPosition');
set(gcf, 'OuterPosition', 1.6*[0  0   560   442]);
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

%%%%%%%%%% Main Loop %%%%%%%%%%

for tstep = 1:T_steps_end
    
    if mod(tstep,1000) == 0
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
        % get x-values of stationary agents, and tabulate into unique
        % values (first column) and counts (second column)
        % tbl = tabulate(agents(~logical(agents(:,2)),1));
        stat_agent_loc = agents(~logical(agents(:,2)),1);
        if max(stat_agent_loc)-min(stat_agent_loc) >= 2^16
            disp('Stationary agents are too spread out to be tabulated using histcounts "Integers"')
            pause
        end
        [counts, edges] = histcounts(stat_agent_loc, 'BinMethod', 'Integers');
        values = floor(edges(2:end));
        % NEW TABLE = values (first row) and counts (second row)
        tbl = [values; counts];
        % change amount of resources at these x values by an amount calculated
        % using the counts, multiplied by time and rate eating
        resources(tbl(1,:)) = resources(tbl(1,:)).*R_exp(tbl(2,:))';
    end
    %else they are all moving!
    
    %%%%%%%%%%%%
    %%% Plot Figure 1 %%% Realtime movie of simulation
    %%%%%%%%%%%%
    
%     % only plot at integer times
    
    % only plot every snapfreq time steps
     snapfreq = ceil(T_end/50);
     if PLOT_FLAG && (mod(tstep,snapfreq)==0 || tstep == 1)%new_int
       pltcnt = pltcnt + 1;
        %figure(1) % repetition is okay
        
        %%% sort out x domain and indices %%%
        
        % get xmin/xmax bounds, incl a little buffer
        
        xmin = max([min(agents(:,1))-5,1]); %5 is a O(1) number to show a bit more than all agents
%         xmin=1;
        xmax = max(agents(:,1))+5;
%         xmax=v*T_end;
        % create local xgrid and scale to put in actual units (instead of index)
        xgrid = (xmin:xmax)*delta_x;
        %%% plot %%%
        %agents histogram
        h1 = subplot(3,1,2);
        %plot histogram of agents
        %histogram(agents(:,1),[xgrid-delta_x/2, xgrid(end)+delta_x/2])
        histogram(delta_x*agents(:,1),'BinMethod','integers',...
            'BinLimits',[xmin*delta_x,xmax*delta_x])
        hold on
        k = round(1/delta_x)+1; %odd to window around center point
        [cnts,edges] = histcounts(delta_x*agents(:,1),...
            [xgrid-delta_x/2, xgrid(end)+delta_x/2]);
        plot(xgrid,movsum(cnts,k),'LineWidth',1.5)
%         ylim([0 500]);
%         xlim([0 25]);
        axis([ xmin*delta_x xmax*delta_x 0 1500])
        title('Histogram of Locusts','FontSize',14)
        hold off
        
        %%%%%%%%%%%%%%%%
        % plot agents position - subplot 1
        %%%%%%%%%%%%%%%%
        
        h2 = subplot(3,1,1);
        plot(agents(:,1)*delta_x,1:N_agents,'.','MarkerSize', 10);
        hold on
%         %also plot a moving sum
%         k = round(1/delta_x)+1; %odd to window around center point
%         [cnts,edges] = histcounts(v*delta_x*agents(:,1),...
%             [xgrid-delta_x/2, xgrid(end)+delta_x/2]);
%         %plot(xgrid,movsum(cnts,k),'LineWidth',1.5)
%         %also plot the position of the first agent on the x-axis
%         %scatter(agents(1,1)*delta_x,0,'filled')

%         xlim([0,25]);
%         ylim([0,length(agents(:,1))])
        axis([xmin*delta_x xmax*delta_x 0 N_agents])
        
        title('Locusts','FontSize',14)
        hold off
        
        %%%%%%%%%%%%%%%%%%%%
        % plot resources - subplot 3
        %%%%%%%%%%%%%%%%%%%%
        
        h3 = subplot(3,1,3);
        %plot a moving mean of resources, taking the average over k points
        %where k scales with delta_x
        
        h=plot(xgrid,resources(xmin:xmax),'LineWidth',2);
%         might want full window of x values later
%         xlim([0,domain_length]); %xlim([0,25])
%         ylim([0,R_plus+10]); %ylim([0,1.6])
        axis([xmin*delta_x xmax*delta_x 0 R_plus+.1*R_plus])
        title('Resources','FontSize',14)
        xlabel('position (m)')
        h(end).Color=green;
        
        if movieflag
            frame = getframe(1);
            writeVideo(vid,frame);
        end

        %make things prettier
%         set(h1,'position',[0.13,0.9,0.7750,0.4]);
%         set(h1,'position',[0.13,0.3,0.7750,0.2]);
%         set(h3,'position',[0.13,0.1,0.7750,0.1]);
        %set(h1,'position',[0.13,0.3838,0.7750,0.5412]);
        %set(h3,'position',[0.13,0.11,0.7750,0.1412]);
        %%%pause(0.1)       
     end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
     %  Record Agent Data    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %%% Record data %%%
    %
    % Sanity check, and record number moving
    assert(max(agents(:,2))<=1 && min(agents(:,2))>=0)
    M_cnt(tstep) = sum(agents(:,2)); %number moving
    
    % Record position of chosen few
    sav_loc(tstep,:) = agents(sav_idx,1)*delta_x;
    
    % Basic statistics
    median_ary(tstep) = median(agents(:,1)*delta_x);
    mean_ary(tstep) = mean(agents(:,1)*delta_x);
    %%% Width
    std_ary(tstep) = std(agents(:,1)*delta_x);
    
    % Calculate max, % mode of locations of pulse = most common location = location of peak
    [max_idx_full, freq] = mode(agents(:,1));
    this_pulse_max = freq/delta_x;
    %%% Peak
    max_ary(tstep) = this_pulse_max;
    
    %%% Skewness
    skew_ary(tstep) = -1/N_agents*sum((agents(:,1)*delta_x-mean_ary(tstep)).^3)/(std_ary(tstep))^3;
    
    %
    %%% Auto detection of end of transients
    %
    if tstep > t_delay_window_len+1 && start_t_idx == 0
        % update window
        t_vel_window = (mean_ary((tstep-t_delay_window_len):tstep) -...
                mean_ary((tstep-1-t_delay_window_len):(tstep-1)))/delta_t;
        % calculate std of vel in window
        vel_std = std(t_vel_window);
        if vel_std < vel_std_thrs % end of transients found!
            start_t_idx = tstep;
            start_x_idx = mode(agents(:,1));
            % initialize pulse shape array
            wdth = (max(agents(:,1))-min(agents(:,1)));
            wdth = max(2*wdth,2000); %10*wdth
            % wdth is width of pulse at end of transients. Since transients
            % are detected using mean velocity, which stabilizes faster
            % than width, we expect width to continue growing. Therefore we multiply width by a factor of 10. 
            if mod(wdth,2) == 0
                wdth = wdth+1; % make it odd
            end
            last_trans_t = tstep-1;
            shape_ary = zeros(T_steps_end-last_trans_t,wdth);
            % Make a shape array for resources
            Rshape_ary = zeros(T_steps_end-last_trans_t,wdth);            
            mode_pt = round((wdth+1)/2); %ignore variable name. This is the median not mode
        end % else do nothing and keep looking next tstep
    end
    
    %
    %%% Post transients, at each step, acquire data on pulse shape %%%
    %
    if start_t_idx ~= 0
        %max_idx = mode(agents(:,1)); % mode of locations of pulse = most common location = location of peak
        mean_idx = round(mean(agents(:,1))); % mean of locations of pulse
        % tbl_x = tabulate(agents(:,1)); % columns: unique vals, counts, percent
        % tbl_x = tbl_x(:,1:2);          % drop third column
        if max(stat_agent_loc)-min(stat_agent_loc) >= 2^16
            disp('Stationary agents are too spread out to be tabulated using histcounts "Integers"')
            pause
        end
        [counts, edges] = histcounts(agents(:,1),'BinMethod','Integers');
        values = floor(edges(2:end));
        % NEW TABLE X = values (first row) and counts (second row)
        tbl_x = [values; counts];
        tbl_x = tbl_x(:,tbl_x(2,:)>0); % ditch the zeros
        min_loc = min(tbl_x(1,:));
        max_loc = max(tbl_x(1,:));
        tbl_x(1,:) = tbl_x(1,:) - mean_idx; %normalize x to the mean location of the band %before: -max_idx (not mean_idx) to align to locust mode
        % Do something here with resources?
        if min(tbl_x(1,:))>-mode_pt && max(tbl_x(1,:)) < mode_pt
            % record the pulse shape in the array, centered around mode_pt
            shape_ary(tstep-last_trans_t,tbl_x(1,:)+mode_pt) = tbl_x(2,:);
            % record resource shape in the array
            Rshape_ary(tstep-last_trans_t,...
                min_loc-mean_idx+mode_pt:max_loc-mean_idx+mode_pt)...
                = resources(min_loc:max_loc);
            % set all values ahead of max_loc = R_plus
            Rshape_ary(tstep-last_trans_t,max_loc-mean_idx+mode_pt+1:end) = R_plus;
        else
            % pulse is too big for shape_ary. throw a warning
            fprintf('Pulse too big for shape_ary! Initialize a bigger array in loop') 
            %to initialize bigger array, increase mode_pt by increasing
            %the variable wdth.
        end
        
    end
    
end

%
%%% Fit exponential curve to average tail shape and plot average pulse %%%
%
if start_t_idx > 0 && start_t_idx < T_steps_end
    shape_ary_mean = mean(shape_ary);
    Rshape_ary_mean = mean(Rshape_ary);
    
    % Find the mode of the mean shape
    [M, max_idx] = max(shape_ary_mean); 
    
    tail_start = find(shape_ary_mean,1); % find location of first nonzero entry
    pulse_end = find(shape_ary_mean,1,'last'); % find location of last nonzero entry
    xmesh = ((-(mode_pt-1):(mode_pt-1))-(max_idx-mode_pt))*delta_x; %shift xmesh back to being centered around the mode
    % get rid of unused parts of the array and calculate std
    mean_pulse = shape_ary_mean(tail_start:pulse_end); % ditch zero entries
    mean_resources = Rshape_ary_mean(tail_start:pulse_end); % look at same segment
    xmesh = xmesh(tail_start:pulse_end);
    
    std_pulse = std(shape_ary(:,tail_start:pulse_end));
    std_resources = std(Rshape_ary(:,tail_start:pulse_end));
    
    pulse_center = max_idx-(tail_start-1); %put pulse center at mode
%     pulse_center = mean_idx; <---- I would like to fix this... :(
    assert(xmesh(pulse_center)==0) % sanity check
    
    % get tail of distribution
    mean_tail = mean_pulse(1:pulse_center);
    xmesh_tail = xmesh(1:pulse_center);
    
    %%% Fit tail to exponential curve %%%
    % y=a*exp(b*x), solve lsq to find a and b
    % mean_pulse may (stochastically) have entries with a value of zero in
    %   the tail. Find these and get rid of them for least squares so that
    %   there will be a feasible solution.
    mean_tail_nonzero = mean_tail(mean_tail~=0);
    xmesh_tail_nonzero = xmesh_tail(mean_tail~=0);
    size(mean_tail_nonzero)
    % On second thought, let's throw out all the really small points
    mean_tail_nonzero = mean_tail(mean_tail>=0.01*M);
    xmesh_tail_nonzero = xmesh_tail(mean_tail>=0.01*M);
    size(mean_tail_nonzero)
    
    Y = log(mean_tail_nonzero');
    X = [ones(length(xmesh_tail_nonzero),1), xmesh_tail_nonzero'];
    B = X\Y; % least squares. B(1)=log(a),  B(2)=b
    % compute root mean square error
    rms_err = sqrt(immse(mean_tail./delta_x,...
        exp(B(1)+B(2).*xmesh_tail)./delta_x));
    fprintf('Least-squares exponential parameters y=a*exp(b*x): a=%g, b=%g\n',...
        B(1)./delta_x,B(2))
    fprintf('Root-mean squared error: %g\n',rms_err)
    
    if PLOT_SUMMARY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot labels for (17) and (18)  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xaxis = 'Displacement from Mean Position (m)';
        %yaxis17 = 'Locusts per gridpoint $S_{n,m}+M_{n,m}$';
        yaxis18 = 'Locusts per meter (1/m)';
        yaxis17 = yaxis18;
        ryaxis = 'Resource Density (g/m)';
        titl17 = 'Time-Average of Outputs from Agent-Based Model, with Fitted Curve';
        titl18 = ['Agent-Based Model Output at $t = $ ' num2str(T_end)];

        axesfontsize = 24;
        labelsize = 24;

        lyt = [0 400 800 1200 1600]; % left y-axis ticks

        resColor = green;
        resWidth = 3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot Agent Data - Figure 17: time-averaged profile %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        h = figure(17);

        pos = get(h, 'Position');
        set(h, 'Position', [round(pos(1)/2),round(pos(1)/2),pos(3)*2,pos(4)*1.25]);
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        hold on
        %overlay a density plot
        y = mean_pulse./delta_x;
        dy = std_pulse./delta_x;
        %errorbar(xmesh,y,dy,'LineWidth',1.5,'DisplayName','locust density')
        d_plot = plot(xmesh,y,...
            'LineWidth',3,'DisplayName','avg locust density');
        %plot logistic fit
    %     plot(xmesh_tail,logistic_func(xmesh_tail,K,r,t0)./delta_x,...
    %         'LineWidth',3,'DisplayName','logistic')
        %plot exponential fit
    %     plot(xmesh_tail,(exp(B(1)+B(2).*xmesh_tail))./delta_x,...
    %         'LineWidth',3,'Color',[0.9290, 0.6940, 0.1250],'DisplayName','exponential')
        legend({},'AutoUpdate','off','location','northwest')
        %add error shading
        plot_color = get(d_plot,'Color');
        f = fill([xmesh';flipud(xmesh')],[y'-dy';flipud(y'+dy')],[.78 .78 .78],'EdgeColor',[.65,.65,.65]);
        %alpha(f,0.3)
        uistack(f,'down',3)
        %add some convience lines
        xl = xlim;
        yl = ylim;
        ylsave = ylim;
        ylim([yl(1)/4, mean_pulse(pulse_center)*1.25/delta_x])
        yl = ylim;
        hline = plot(xl,[0,0],'k');
        vline = plot([0,0],yl,'k--');
        uistack(hline,'bottom')
        uistack(vline,'bottom')

        %add mean resource curve
        yyaxis right
        Ry = mean_resources;
        dRy = std_resources;
        %add error shading
        %plot_color = get(resplot,'Color');
        fR = fill([xmesh';flipud(xmesh')],[Ry'-dRy';flipud(Ry'+dRy')],[.78 .78 .78],'EdgeColor',[.65,.65,.65]);
        %alpha(f,0.3)
    %     uistack(fR,'down',2)
        resplot = plot(xmesh,Ry,...resources(mean_idx-(pulse_center-1):mean_idx+(length(xmesh)-pulse_center)),...
            '-','Color',resColor,'LineWidth',resWidth);
        ylim([(yl(1)/yl(2))*R_plus*1.1,R_plus*1.1])
        set(gca,'ycolor',resColor)
        ylabel(ryaxis,'Color','k','FontSize',labelsize) %fontsize doesn't stick for some reason...

        %make the original plot prettier, add labels
        yyaxis left
        set(gca,'FontSize',axesfontsize)
        xlabel(xaxis,'FontSize',labelsize)
        yticks(lyt)
        ylabel(yaxis17,'FontSize',labelsize)
        title(titl17,...
            'FontSize',labelsize)
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
        if savePlots == 1
            saveas(gcf,'ABM_smooth.eps','epsc')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Plot Agent Data - Figure 18: agent-based profile %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        h = figure(18);

        pos = get(h, 'Position');
        set(h, 'Position', [round(pos(1)/2),round(pos(1)/2),pos(3)*2,pos(4)*1.25]);
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        hold on

        %add some convience lines
    %     xl = xlim;
        yl = ylsave;
        ylim([yl(1)/4, mean_pulse(pulse_center)*1.25./delta_x]);
        yl = ylim;
        hline = plot(xl,[0,0],'k');
        vline = plot([0,0],yl,'k--');
        uistack(hline,'bottom')
        uistack(vline,'bottom')

        %add last histogram as an example
        %add plot of number of locusts in each dx bin
        %k = round(1/delta_x)+1; %number of bins over which to average, small k => rougher curve
        max_idx = mode(agents(:,1)); %center around mode of locations of pulse
        [cnts,edges] = histcounts(delta_x*(agents(:,1)-max_idx),...
                        [xmesh-delta_x/2, xmesh(end)+delta_x/2]);
        cntplot = plot(xmesh,movsum(cnts,1)/delta_x,'LineWidth',1.5, 'Color', [0.85 0.325 0.098]);
        uistack(cntplot, 'bottom')

        %add resource curve corresponding to histogram
        yyaxis right
        clrs = parula(100);
        %resources(min(agents(:,1))-max_idx+mode_pt:max(agents(:,1))-max_idx+mode_pt),...
        resplot = plot(xmesh,resources(max_idx-(pulse_center-1):max_idx+(length(xmesh)-pulse_center)),...
            ':','Color',resColor,'LineWidth',resWidth);
        ylim([(yl(1)/yl(2))*R_plus*1.1,R_plus*1.1])
        set(gca,'ycolor',resColor)
        ylabel(ryaxis,'Color','k') %fontsize doesn't stick for some reason...

        %make the original plot prettier, add labels
        yyaxis left
        set(gca,'FontSize',axesfontsize)
        xlabel(xaxis,'FontSize',labelsize)
        yticks(lyt)
        ylabel(yaxis18,'FontSize',labelsize)
        title(titl18,...
            'FontSize',labelsize)
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

            if savePlots == 1
                saveas(gcf,'ABM_stoch.eps','epsc')
            end
            pause
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Time varying stats %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get a time and space grid, leave off first time point (velocity)
tgrid = (2:T_steps_end)*delta_t;
xgrid = (1:length(resources))*delta_x;

% velocity of the mean
vel_mean = (mean_ary(2:end) - mean_ary(1:end-1))/delta_t;

%
%%% Report summary statistics %%%
%
if start_t_idx > 0 && start_t_idx < T_steps_end
    cspeed = mean(vel_mean(start_t_idx:end));
    fprintf('Mean velocity of mean position after t=%g: %g\n',...
        start_t_idx*delta_t,cspeed);
    cstd = std(vel_mean(start_t_idx:end));
    fprintf('Std of mean velocity after t=%g: %g\n',...
        start_t_idx*delta_t,cstd);
    if min(agents(:,1)) > start_x_idx+5
        last_x_ind = min(agents(:,1))-1;
        R_minus = mean(resources(start_x_idx:last_x_ind));
        fprintf('Mean resources left in interval %g<x<%g behind the wave: %g\n',...
            start_x_idx*delta_x,last_x_ind*delta_x,...
            R_minus);
        Rstd = std(resources(start_x_idx:last_x_ind));
        fprintf('Std deviation of these resources: %g\n',...
            Rstd);
    else
        fprintf('For resource computation, run longer.\n')
    end
else
    fprintf(strcat('End of transients went undetected. Limiting behavior stats unavailable.\n',...
        'Examine the velocity of mean position plot and possibly adjust the threshold value.\n'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set output variables
% additional, averaged shapes above in transient dedection loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shape_stats = [max_ary, std_ary, skew_ary];

fig3outs = cell(1,3);
fig3outs{1} = vel_mean;
fig3outs{2} = sav_loc;
fig3outs{3} = mean_ary;

fig4outs = cell(1,10);
fig4outs{1} = mean_pulse; % mean_pulse_out
fig4outs{2} = mean_resources; % mean_resources_out 
fig4outs{3} = xmesh; % xmesh_out
fig4outs{4} = xmesh_tail;
fig4outs{5} = std_pulse; 
fig4outs{6} = std_resources;
fig4outs{7} = B; % constants for fitted exponential
fig4outs{8} = pulse_center;
fig4outs{9} = agents;
fig4outs{10} = resources;

fig6outs = cell(1,2);
fig6outs{1} = cspeed;
fig6outs{2} = R_minus;

axesfontsize = 15;
labelsize = 20;


%%%%%%%%%%%%%%%%%%%%%%%
% Plot Stats - Figure 3
%%%%%%%%%%%%%%%%%%%%%%%
if pcheck_plots
    % Median is noisy! - use velocity of mean
    h = figure(3);

    pos = get(h, 'Position');
    set(h, 'Position', [round(pos(1)/2),round(pos(1)/2),pos(3)*1,pos(4)*1.5]);
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    subplot(4,1,1)
    plot(tgrid,vel_mean,'LineWidth',2);
    title('Velocity of mean position over time','FontSize',14)
    xlabel('time (s)')
    ylim([0,2*max(vel_mean)])
    'average velocity is'
    mean(vel_mean(floor(0.75*length(vel_mean)):end))

    subplot(4,1,2)
    plot([0 tgrid],max_ary,'LineWidth',2);
    title('Maximum density over time','FontSize',14)
    xlabel('time (s)')
    ylim([0,1.25*max(max_ary)])
    'the final peak density is'
    max_ary(end)

    subplot(4,1,3)
    plot([0 tgrid],std_ary,'LineWidth',2);
    title('Std of position over time','FontSize',14)
    xlabel('time (s)')
    ylim([0,1.25*max(std_ary)])
    
%     hold on
%     for ii = 1:size(sav_loc,2)
%         plot((1:T_steps_end)*delta_t,sav_loc(:,ii)-mean_ary);
%     end
%     hold off
%     xlabel('time (s)')
%     title('Position from pulse mean of selected locusts','FontSize',14)

    % subplot(3,2,3)
    % trunc_xgrid = xgrid(xgrid<=mean_ary(end));
    % plot(trunc_xgrid,R_init(xgrid<=mean_ary(end)))
    % xlabel('position (m)')
    % xlim([xgrid(1),trunc_xgrid(end)])
    % title('Starting resources encountered','FontSize',14)
    subplot(4,1,4)
    plot([0 tgrid],skew_ary,'LineWidth',2);
    title('Skewness of profile over time','FontSize',14)
    xlabel('time (s)')
    ylim([min(1.25*min(skew_ary),-.01),max(1.25*max(skew_ary),.01)])
%     
%     h=plot(xgrid,resources,'LineWidth',2);
%     title('Resources left as a function of space','FontSize',14)
%     xlabel('position (m)')
% %     xlim([0,100]);
%     ylim([0,1.1*R_plus]);
%     h(end).Color=green;
%     'average resources leftover'
%     mean(resources(ceil(0.25*length(resources)):ceil(0.5*length(resources))))
    
    figure(19)
    % plot the location of the locusts
    [cnts, edges] = histcounts(agents(:,1),'BinMethod','Integers');
    x_ind = floor(edges(2:end));
    x_grid = x_ind*delta_x;
    plot(x_grid,cnts/delta_x)
    ylabel('Locust Density','Color','k')
    xlabel('Space (not spatial grid, not shifted)')
    
    %add resource curve corresponding to histogram
    yyaxis right
    %resources(min(agents(:,1))-max_idx+mode_pt:max(agents(:,1))-max_idx+mode_pt),...
    plot(x_grid,resources(x_ind),...
        ':','Color',green,'LineWidth',3);
    ylim([0,R_plus*1.1])
    set(gca,'ycolor',green)
    ylabel('Resource Density','Color','k') %fontsize doesn't stick for some reason...
    
    if savePlots == 1
        saveas(gcf,'ABM_stats.eps','epsc')
    end
    
end

%
%%% Create movie if requested %%%
%
    if movieflag
        close(vid);
    end
end
    
function y = logistic_func(t,K,r,t0)
y=K./(1+exp(-r*(t-t0)));
end

function z = penalty_logistic(params,t,y)
r=params(1);t0=params(2);

h=1./(1+exp(-r*(t-t0)));

numerator  =  -dot(h,y)^2;
denominator =  dot(h,h);
z=numerator/denominator;
end