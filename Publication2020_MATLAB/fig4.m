close all
clear all

addpath('functions','data')

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

% fonts
titlesize = 12;
axislabelsize = 10;
ticklabelsize = 8;
subfiglabelsize = 12;

%colors
ksm_color = 1/255*[237, 143, 33];   %gold
kms_color = 1/255*[111,73,179];     %purple
green =  [  0    0.5000      0];
gold = [0.9290, 0.6940, 0.1250];
orange = [0.85 0.325 0.098];

%linewidth
lwidth = 1;

%new_data = 1 to run the script and overwrite.
%         = 0 to load existing data
new_data = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 4 %%%
number = 4;
filename = sprintf('fig%d',number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% biological parameters
set_params
N_agents = N;
R_plus = Rplus;

% % Number of locusts that we will track as individuals
% N_sav=5;

% How long to run
T_steps_end = 15000;
% Time step size. Do not increase above 1 or things will break!!
% Also, don't make too small or things will get really noisy...
delta_t = 1;

% Initialize a spatial domain width and delta_x based on delta_t
domain_length = T_steps_end;%+20;
delta_x = v*delta_t;
T_end = T_steps_end*delta_t;

% pause and check plots
pcheck_plots = 0;

% Initial Conditions
agents_in = zeros(N_agents,2); % Columns are x_ind, S/M
% Stationary = 0, Moving = 1
% start all locusts in the first 5 spatial points, equally spaced, stationary
Ndiv5 = floor(N_agents/5);
for n=1:5
    agents_in((Ndiv5*(n-1)+1):Ndiv5*n,1) = n; %set spatial locations
end
if Ndiv5*5 ~= N_agents %check for remaining locusts
    agents_in((Ndiv5*5+1):end) = 1; %put them in the first grid point
end
agents_in(:,2) = 0; %set all agents to stationary
R_in = ones(domain_length,1).*R_plus; %units: resources per delta_x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run Script %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if new_data
    [fig3outs, fig4outs,fig6outs] =...
             fxn_ABM_shift(agents_in,R_in,T_steps_end,delta_t,N_agents,R_plus,v,lambda,alpha,beta,eta,theta,gamma,delta,pcheck_plots);
    
    mean_pulse = fig4outs{1}; % mean_pulse_out
    mean_resources = fig4outs{2}; % mean_resources_out 
    xmesh = fig4outs{3}; % xmesh_out
    xmesh_tail = fig4outs{4};
    std_pulse = fig4outs{5}; 
    std_resources = fig4outs{6};
    B = fig4outs{7}; % constants for fitted exponential
    pulse_center = fig4outs{8};
    agents = fig4outs{9};
    resources = fig4outs{10};

    save(['data/' filename], 'fig4outs', 'T_end','delta_x','R_plus')
else
    load(['data/' filename])
    
    mean_pulse = fig4outs{1}; % mean_pulse_out
    mean_resources = fig4outs{2}; % mean_resources_out 
    xmesh = fig4outs{3}; % xmesh_out
    xmesh_tail = fig4outs{4};
    std_pulse = fig4outs{5}; 
    std_resources = fig4outs{6};
    B = fig4outs{7}; % constants for fitted exponential
    pulse_center = fig4outs{8};
    agents = fig4outs{9};
    resources = fig4outs{10};
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(number);


    xaxis = 'Distance from Mode Position (m)';
    yaxis18 = 'Locusts per Gridpoint'; %$S_{n,m}+M_{n,m}$
    yaxis17 = 'Locusts per meter (1/m)';
    %yaxis17 = yaxis18;
    ryaxis = 'Resource Density (g/m)';
    titl17 = 'Average Agent-Based Output, with Fitted Curve';
    titl18 = ['Agent-Based Model Output at $t = $ ' num2str(T_end)];
    
    axesfontsize = ticklabelsize;
    labelsize = axislabelsize;
    
    lyt = [0 400 800 1200 1600]; % left y-axis ticks
    lyt_stoch = [0 16 32 48 64];

    resColor = green;
    resWidth = 2;
    
    abm_width = 1;
    smooth_width = 2;
    exp_width = 2;

pos = get(h, 'Position');
% %set(h, 'Position', [-round(pos(1)/2),round(pos(1)/2),pos(3)*2,pos(4)*1]); %from original
% set(h, 'Position', [-round(pos(1)/2),round(pos(1)/2),pos(3)*2,pos(4)*2.5]); %from original
set(h,'Units','Inches');
% pos = get(h,'Position')
set(h,'PaperPositionMode','Manual') % Setting this to 'manual' unlinks 'Position' (what gets displayed) and 'PaperPosition' (what gets printed)
set(h,'PaperPosition',[ 0 0 6.5 6.5]);
set(h,'Position',[ 0 0 6.5 6.5]);
get(h,'Position');

% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[2*pos(3), pos(4)])
hold on

%%% ABM SMOOTH %%%
h2 = subplot(2,1,2);
pos = get(h2, 'Position');
set(h2,'Position', [pos(1),pos(2), pos(3), pos(4)*1])
hold on

y = mean_pulse./delta_x;
    dy = std_pulse./delta_x;
    %errorbar(xmesh,y,dy,'LineWidth',1.5,'DisplayName','locust density')
    d_plot = plot(xmesh,y,...
        'LineWidth',smooth_width,'DisplayName','avg locust density');
    %plot logistic fit
%     plot(xmesh_tail,logistic_func(xmesh_tail,K,r,t0)./delta_x,...
%         'LineWidth',3,'DisplayName','logistic')
    %plot exponential fit
    plot(xmesh_tail,(exp(B(1)+B(2).*xmesh_tail))./delta_x,'--',...
        'LineWidth',exp_width,'Color',gold,'DisplayName','exponential')
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
    ylim([yl(1)/2/delta_x, mean_pulse(pulse_center)*1.25/delta_x])
    yl = ylim;
    hline = plot(xl,[0,0],'k');
    vline = plot([0,0],yl,'k--');
    uistack(hline,'bottom')
    uistack(vline,'bottom')
    set(gca,'FontSize',axesfontsize)

    %add resource curve corresponding to histogram
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
    ylabel(ryaxis,'Color','k','FontSize',axislabelsize) %fontsize doesn't stick for some reason...
    
    %make the original plot prettier, add labels
    yyaxis left
    set(gca,'FontSize',axesfontsize)
    xlabel(xaxis,'FontSize',axislabelsize)
    yticks(lyt)
    ylabel(yaxis17,'FontSize',axislabelsize)
    title(titl17,...
        'FontSize',titlesize)
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

% subfigure label
annotation('textbox', [0, 0.525, 0, 0], 'string', '\textbf{(B)}', 'FontSize', subfiglabelsize,'Interpreter','latex')

%%% ABM STOCHASTIC %%%
h1 = subplot(2,1,1);
pos = get(h1, 'Position');
set(h1,'Position', [pos(1),pos(2)*1, pos(3), pos(4)*1])

    %add some convience lines
%     xl = xlim;
hold on
    yl = ylsave;
    ylim([yl(1)/2, mean_pulse(pulse_center)*1.25]);
    yl = ylim;
    hline = plot(xl,[0,0],'k');
    vline = plot([0,0],yl,'k--');
    uistack(hline,'bottom')
    uistack(vline,'bottom')

    %add plot of number of locusts in each dx bin
    k = round(1/delta_x)+1; %number of bins over which to average, small k => rougher curve
    max_idx = mode(agents(:,1)); %center around mode of locations of pulse
    [cnts,edges] = histcounts(delta_x*(agents(:,1)-max_idx),...
                    [xmesh-delta_x/2, xmesh(end)+delta_x/2]);
    cntplot = plot(xmesh,movsum(cnts,1),'LineWidth',abm_width, 'Color', orange);
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
    set(gca,'FontSize',ticklabelsize)
    xlabel(xaxis,'FontSize',axislabelsize)
    yticks(lyt_stoch)
    ylabel(yaxis18,'FontSize',axislabelsize)
    title(titl18,...
        'FontSize',titlesize)
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
    


% subfigure label
annotation('textbox', [0, 1, 0, 0], 'string', '\textbf{(A)}', 'FontSize', subfiglabelsize,'Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(h,[filename, '.eps'],'-depsc')
set(h,'PaperPositionMode','Auto')

