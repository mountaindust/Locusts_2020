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
blue = [      0    0.4470    0.7410];
defwhite = [1 1 1];

%linewidth
lwidth = 2;

%new_data = 1 to run the script and overwrite.
%         = 0 to load existing data
new_data = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 6 %%%
number = 6;
filename = sprintf('fig%d',number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% biological parameters
set_params

%%% Numerical Parameters
L = 100;     %domain length
nx = 2000;  %spatial points
nt = 10000;  %time steps % MUST BE a multiple of 3

%%% Computational grid
dx = L/nx;  %spatial discritization
dt = dx/v;  %time step -- chosen to for convenience
X = (1:nx)*dx;

%%% Set a range of N and a few Rps %%%

Nmin = 3500;
Nmax = 16000;

Npts = 50;
Nspace=linspace(Nmin,Nmax,Npts);

Rps = [100 150 200 250 300];
Rpts = length(Rps);

% initialize variables for storage
csolve = zeros(Npts,Rpts);
Rmsolve = zeros(Npts,Rpts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run Script %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if new_data
    
    % solve for corresponding c's, Rm's
    
    for iR = 1:Rpts

        Rp = Rps(iR);

        for iN = 1:Npts
            N = Nspace(iN);

            [Rm, c] = fxn_findcRm(N,Rp,v,lambda,alpha,beta,eta,theta,gamma,delta);

            csolve(iN, iR) = c;
            Rmsolve(iN,iR) = Rm;
        end

    end
    
    % now for the ABM data
    T_end = 15000; %2*dt*nt; % The actual end time, NOT the total number of time steps
    delta_t = 1; %dt;
    PLOT_FLAG = 0;
    
    % Initialize a spatial domain width and delta_x based on delta_t
    T_steps_end = T_end; assert(delta_t==1); % because delta_t = 1;
    domain_length = T_steps_end;%+20;
    delta_x = v*delta_t;

    n = 8; % create n data points via the ABM
    Ns_ABM = linspace(5000,15000,n-3);
    Ns_ABM = [5250  5500  6000  Ns_ABM];
    Rps_ABM = Rps(ceil(Rpts/2)); % data should land on the middle Rp value

    % initialize data arrays
    cs_ABM = zeros(1,n);
    cstd_ABM = zeros(1,n);
    Rs_ABM = zeros(1,n);
    Rstd_ABM = zeros(1,n);
    
    for k = 1:n
        N_agents = Ns_ABM(k);
        R_plus = Rps_ABM;%(k)
        
        %%% initial conditions:
        % start all locusts in the first 5 spatial points, equally spaced, stationary
        agents = zeros(N_agents,2);
        Ndiv5 = floor(N_agents/5);
        for n=1:5
            agents((Ndiv5*(n-1)+1):Ndiv5*n,1) = n; %set spatial locations
        end
        if Ndiv5*5 ~= N_agents %check for remaining locusts
            agents((Ndiv5*5+1):end) = 1; %put them in the first grid point
        end
        agents(:,2) = 0; %set all agents to stationary
        agents_in = agents;
        
        R_in = ones(domain_length,1) * R_plus;
        
        %%% run the script %%%
        pcheck_plots = 0;
        [fig3outs, fig4outs, fig6outs] = fxn_ABM_shift(agents_in,R_in,T_steps_end,delta_t,N_agents,R_plus,v,lambda,alpha,beta,eta,theta,gamma,delta,pcheck_plots);
                
        % extract data from run
        cspeed = fig6outs{1};
        R_minus = fig6outs{2};
        
        % assign to recording variables
        cs_ABM(k) = cspeed;
        %cstd_ABM(k) = cstd;
        Rs_ABM(k) = R_minus;
        %Rstd_ABM(k) = Rstd;
    end
    
    fig6save = cell(2,4);
    % columns are:   [N vals      Rp vals     c vals      Rm vals]
    % first row data via equations from PDE Theory
    fig6save{1,1} = Nspace;
    fig6save{1,2} = Rps;
    fig6save{1,3} = csolve;
    fig6save{1,4} = Rmsolve;
    % second row data via simulation of ABM
    fig6save{2,1} = Ns_ABM;
    fig6save{2,2} = Rps_ABM;
    fig6save{2,3} = cs_ABM;
    fig6save{2,4} = Rs_ABM;
    
    save(['data/' filename], 'fig6save')
    
else
    load(['data/' filename])
%     fig6save = fig6outs;
    % assign variables from fig6save
    % columns are:   [N vals      Rp vals     c vals      Rm vals]
    % first row data via equations from PDE Theory
    Nspace = fig6save{1,1};
    Rps = fig6save{1,2}; 
    csolve = fig6save{1,3};
    Rmsolve = fig6save{1,4};
    % second row data via simulation of ABM
    Ns_ABM = fig6save{2,1};
    Rps_ABM = fig6save{2,2};
    cs_ABM = fig6save{2,3};
    Rs_ABM = fig6save{2,4};
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
pos = [ 0 0  pos(3)*fullwidth/pos(3) pos(4)*fullwidth/pos(3)*1/2];
set(h,'PaperPosition',pos);
set(h,'Position',pos);
get(h,'Position');


levels = linspace(130,250,5);
curvelabel = ticklabelsize;
% linecolor = [0    0.4470    0.7410];
dotsize = 75;
defaultOrange = orange;

%%% C VS MASS %%%
h2 = subplot(1,2,1);
set(h2,'Units','Inches');
pos = get(h2, 'Position'); % [0.1300    0.5838    0.3347    0.3412] from full figure [0 0 6.5 5];
set(h2,'Position', [pos(1),1.25*pos(2), 0.9*pos(3), 0.85*pos(4)])
%set(gcf,'Position',[100 100 scrsz(3) .25*scrsz(3)])

for k = 1:Rpts
        plot(   Nspace, csolve(:,k), 'k',...
                'LineWidth', lwidth...
            );
        idx = find(Nspace >= 15000,1);
        str = sprintf('%.0f',Rps(k));
        text(Nspace(idx),csolve(idx-1,k),str);
        hold on
end

    % beautify
    set(gca,'FontSize',ticklabelsize)
    xlabel('Mass $N$','FontSize',axislabelsize)
    ylabel('Speed $c$','FontSize',axislabelsize)
    %set(findall(gcf,'-property','FontSize'),'Fontsize',fS);
    axis([5000-100 15000+100 4.5*10^(-3) 6.5*10^(-3)])
    
    
% add AMB data points and error bars
hold on
% errorbar(Ns_ABM(4:end), cs_ABM(4:end), cstd_ABM(4:end), 'LineStyle', 'none')
sN = scatter(Ns_ABM(4:end), cs_ABM(4:end),dotsize);
sN.MarkerEdgeColor = defaultOrange;
sN.LineWidth = lwidth;
%sN.MarkerFaceColor = defaultOrange; % for filled circles

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


%%% C VS RMINUS %%%
h3 = subplot(1,2,2);
set(h3,'Units','Inches');
pos = get(h3, 'Position'); %[0.1300    0.1100    0.3347    0.3412];% from full figure [0 0 6.5 5];
set(h3,'Position', [pos(1),1.25*pos(2), 0.9*pos(3), 0.85*pos(4)])
%set(gcf,'Position', [100 100 .45*scrsz(3) .25*scrsz(3)])%[100 100 .45*scrsz(3) .12*scrsz(3)])

for k = 1:Rpts
        plot(   Rmsolve(:,k), csolve(:,k), 'k',...log(Rmsolve(:,k)), csolve(:,k), 'k',...
                'LineWidth', lwidth...
            );
        idx = find(Rmsolve(:,k) <= 0.04,1);
        str = sprintf('%.0f',Rps(k));
        if k == Rpts(end)
            idx2 = find(Rmsolve(:,k) <= 0.02,1);
            text(Rmsolve(idx2,k),csolve(idx2,k),str);
        else
            text(Rmsolve(idx,k),csolve(idx,k),str);
        end
        hold on
end
 
    % beautify
    set(gca,'FontSize',ticklabelsize)
    xlabel('Remaining Resources $R^-$','FontSize',axislabelsize)
    ylabel('Speed $c$','FontSize',axislabelsize)
    %axis([-16 log(0.1) 4.5*10^(-3) 6.5*10^(-3)])
    axis([0 0.04 4.5*10^(-3) 6.5*10^(-3)])
 
% add ABM data points and error bars
% hold on
% errorbar(log(Rs_ABM), cs_ABM, cstd_ABM, 'LineStyle', 'none')
% errorbar(log(Rs_ABM), cs_ABM, exp(Rstd_ABM), 'horizontal','LineStyle', 'none')
% scatter(log(Rs_ABM), cs_ABM, dotsize, 'r','filled')
hold on
% errorbar(Rs_ABM(1:5), cs_ABM(1:5), cstd_ABM(1:5), 'LineStyle', 'none')
% errorbar(Rs_ABM(1:5), cs_ABM(1:5), Rstd_ABM(1:5), 'horizontal','LineStyle', 'none')
sRm = scatter(Rs_ABM(1:5), cs_ABM(1:5),dotsize);
sRm.MarkerEdgeColor = defaultOrange;
sRm.LineWidth = lwidth;

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
annotation('textbox', [.485, 1, 0, 0], 'string', '\textbf{(B)}', 'FontSize', subfiglabelsize,'Interpreter','latex')


%%% Quantifying difference between simulation and theory

% interpolate theory values onto ABM values for N
theory_mass = interp1(Nspace,csolve(:,3),Ns_ABM); % use k=3 for middle curve R^+ = 200
diff = theory_mass-cs_ABM;
percent_mass = diff./theory_mass*100;

disp(['The speed as a function of mass matches to within ' num2str(max(percent_mass)) ' %'])

%interpolate theory values onto ABM values for R^-
theory_Rm = interp1(Rmsolve(:,k), csolve(:,k),Rs_ABM);
diff = theory_Rm-cs_ABM;
percent_Rm = diff./theory_Rm;

disp(['The speed as a function of Rm matches to within ' num2str(max(percent_Rm)) ' %'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(h,[filename, '.eps'],'-depsc')
set(h,'PaperPositionMode','Auto')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Local Functions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = Ival(Rp,Rm,c1,c2,c3)
    out = c1*(log(Rp)-log(Rm)) - (c1-c2)*(expint(c3.*Rm)-expint(c3.*Rp));
end