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

%linewidth
lwidth = 1;

%new_data = 1 to run the script and overwrite.
%         = 0 to load existing data
new_data = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 3 %%%
number = 3;
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
    [fig3outs, fig4outs, fig6outs] =...
             fxn_ABM_shift(agents_in,R_in,T_steps_end,delta_t,N_agents,R_plus,v,lambda,alpha,beta,eta,theta,gamma,delta,pcheck_plots);
    
    vel_mean = fig3outs{1};
    sav_loc = fig3outs{2};
    mean_ary = fig3outs{3};
    
    tgrid = (2:T_steps_end)*delta_t;
    save(['data/' filename], 'fig3outs','tgrid','T_steps_end','delta_t')
else
    load(['data/' filename])
    
    vel_mean = fig3outs{1};
    sav_loc = fig3outs{2};
    mean_ary = fig3outs{3};
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(number);

pos = get(h, 'Position');
% %set(h, 'Position', [-round(pos(1)/2),round(pos(1)/2),pos(3)*2,pos(4)*1]); %from original
% set(h, 'Position', [-round(pos(1)/2),round(pos(1)/2),pos(3)*2,pos(4)*2.5]); %from original
set(h,'Units','Inches');
% pos = get(h,'Position')
set(h,'PaperPositionMode','Manual') % Setting this to 'manual' unlinks 'Position' (what gets displayed) and 'PaperPosition' (what gets printed)
set(h,'PaperPosition',[ 0 0 6.5 4.5]);
set(h,'Position',[ 0 0 6.5 4.5]);
get(h,'Position')

% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[2*pos(3), pos(4)])
hold on

%%% MEAN SPEED %%%
h1 = subplot(2,1,1);
pos = get(h1, 'Position');
set(h1,'Position', [pos(1),pos(2)*1.025, pos(3), pos(4)*.9])
plot(tgrid,vel_mean,'LineWidth',lwidth);
% axes limits, labels
ylim([0,0.007])
ax = axis;
axis([0 10000 ax(3) ax(4)])
ax = axis;

% add convencience line
hold on
mean_speed = mean(vel_mean(floor(0.75*length(vel_mean)):end));
hline = plot(xlim,[mean_speed,mean_speed],'k--');
hold off

% titles and ticks
set(gca,'FontSize',ticklabelsize)
xticks(round(linspace(ax(1),ax(2),5)))
xlabel('Time (s)','FontSize',axislabelsize)
yticks([0 0.003 round(mean_speed*10000)/10000 0.006])
ylabel('Speed (m/s)','FontSize',axislabelsize)
title('Speed of Center of Mass','FontSize',titlesize)
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

%%% INDIVIDUAL LOCUSTS %%%
h2 = subplot(2,1,2);
pos = get(h2, 'Position');
set(h2,'Position', [pos(1),pos(2), pos(3), pos(4)*.9])
hold on
for ii = 1:size(sav_loc,2)
    plot((1:T_steps_end)*delta_t,sav_loc(:,ii)-mean_ary,'LineWidth',lwidth);
end
hold off
ax = axis;
axis([0 10000 ax(3) ax(4)])
ax = axis;
set(gca,'Box','on')
% add convencience line
hold on
hline = plot(xlim,[0,0],'k--');
hold off

% titles and ticks
set(gca,'FontSize',ticklabelsize)
xticks(round(linspace(ax(1),ax(2),5)))
xlabel('Time (s)','FontSize',axislabelsize)
yticks([-10 -5 0 5])%round(linspace(ax(3),ax(4),5)))
ylabel('Distance from Mean (m)','FontSize',axislabelsize)
title('Paths of Sample Locusts','FontSize',titlesize)
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
annotation('textbox', [0, 0.5, 0, 0], 'string', '\textbf{(B)}', 'FontSize', subfiglabelsize,'Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(h,[filename, '.eps'],'-depsc')
set(h,'PaperPositionMode','Auto')

