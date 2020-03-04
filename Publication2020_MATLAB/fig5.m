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
lwidth = 1;

%new_data = 1 to run the script and overwrite.
%         = 0 to load existing data
new_data = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 5 %%%
number = 5;
filename = sprintf('fig%d',number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% biological parameters
set_params

%%% Numerical Parameters
L = 100;     %domain length
nx = 1*2000;  %spatial points
nt = 4*3000;  %time steps % MUST BE a multiple of 3

%%% Computational grid
dx = L/nx;  %spatial discritization
dt = dx/v  %time step -- chosen to for convenience

X = (1:nx)*dx;

%%% Report Info
tfinal = nt*dt

sprSpeed = v*alpha/(alpha+beta)
% 0.004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run Script %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if new_data
    
    %%% run PDE simulation %%%
    video_flag = 0;
    [R, S, M, rhomax, c] = fxn_locustRK(L,nx,nt,N,Rplus,v,lambda,alpha,beta,eta,theta,gamma,delta,video_flag);

    %%% spread out parameters
    %%% chosen so that speed is same for spreading pulse
    % beta stays the same
    theta = beta; % so that kms is constant

    alpha = beta/(v/c-1); % by solving c = v*alpha/(alpha+beta), from Telegrapher's equation
    eta = alpha; % so that ksm is constant

    [Rspr, Sspr, Mspr, rhomaxSpr] = fxn_locustRK(L,nx,nt,N,Rplus,v,lambda,alpha,beta,eta,theta,gamma,delta,video_flag);
    
    fig5outs = cell(2,5); % 2 x 5 cell with first row from pulse, second row from spread.
    fig5outs{1,1} = R;
    fig5outs{1,2} = S;
    fig5outs{1,3} = M;
    fig5outs{1,4} = rhomax;
    fig5outs{1,5} = c;
    
    fig5outs{2,1} = Rspr;
    fig5outs{2,2} = Sspr;
    fig5outs{2,3} = Mspr;
    fig5outs{2,4} = rhomaxSpr;
%     fig5outs{2,5} = c;
    
    save(['data/' filename], 'fig5outs')
else
    load(['data/' filename])
    
    R = fig5outs{1,1};
    S = fig5outs{1,2};
    M = fig5outs{1,3};
    rhomax = fig5outs{1,4};
    c = fig5outs{1,5};
    
    Rspr = fig5outs{2,1};
    Sspr = fig5outs{2,2};
    Mspr = fig5outs{2,3};
    rhomaxSpr = fig5outs{2,4};
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(number);

scfactor = .7;
fullwidth = 6.5;
width = fullwidth*scfactor;

pos = get(h, 'Position');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Manual') % Setting this to 'manual' unlinks 'Position' (what gets displayed) and 'PaperPosition' (what gets printed)
pos = [ 0 0  fullwidth 1.25*fullwidth/2];
set(h,'PaperPosition',pos);
set(h,'Position',pos);
get(h,'Position');

% tick positions
tFifteenth = nt/15;
tSixth = nt/6;
tThird = nt/3;
tmid = nt/2;
t2Third = 2*nt/3;
t5Sixth = 5*nt/6;

tick1 = tSixth;
tick2 = tmid;
tick3 = nt;

% colors
RdepColor = blue;
al = .7;
RindepColor = [0 0 0] + al;

%map color
map = [defwhite; RindepColor; RdepColor];

%Graph features
LWsmall = 1; %line width size
LWbig = 2;
LSrd = '-'; %line style resource dependent
LSri = '--';%line style resource independent
tlabH = .61;


%%% SPREADING PULSES %%%
h2 = subplot(2,2,1);
set(h2,'Units','Inches');
pos = get(h2, 'Position'); % [0.1300    0.5838    0.3347    0.3412] from full figure [0 0 6.5 5];
set(h2,'Position', [pos(1),1.05*pos(2), 1.2*width, width/4])
%set(gcf,'Position',[100 100 scrsz(3) .25*scrsz(3)])

p = plot(   X, S(1:nx,tick1)+M(1:nx,tick1),...
            X, S(1:nx,tick2)+M(1:nx,tick2),...
            X, S(1:nx,tick3)+M(1:nx,tick3),...
            X, Sspr(1:nx,1)+Mspr(1:nx,1),...
            X, Sspr(1:nx,tick1)+Mspr(1:nx,tick1),...
            X, Sspr(1:nx,tick2)+Mspr(1:nx,tick2),...
            X, Sspr(1:nx,tick3)+Mspr(1:nx,tick3),...
            X, S(1:nx,1)+M(1:nx,1)...
        );

    axis([X(1) X(end) 0 1500])
    p(1).Color = RdepColor;
    p(2).Color = RdepColor;
    p(3).Color = RdepColor;
    p(4).Color = RdepColor;
    p(5).Color = RindepColor;
    p(6).Color = RindepColor;
    p(7).Color = RindepColor;
    p(8).Color = RindepColor;
    for i = 1:4
        p(i).LineWidth = LWsmall;
        p(i).LineStyle = LSrd;
    end
    for k = 5:8
        p(k).LineWidth = LWbig;
        p(k).LineStyle = LSri;
    end
    uistack(p(6),'top')

    set(gca,'FontSize',ticklabelsize)
    str = '$t = 0$'; %create textbox
    dim = [.08 tlabH .3 .3];
    annotation('textbox',dim,'String',str,'FitBoxToText','on','edgecolor','none','Interpreter','Latex','Fontsize',axislabelsize)
    %set(findall(gcf,'-property','FontSize'),'Fontsize',fS);
    str = ['$t =\, $' num2str(tick1*dt)]; %create textbox
    dim = [.29 tlabH .3 .3];
    annotation('textbox',dim,'String',str,'FitBoxToText','on','edgecolor','none','Interpreter','Latex','Fontsize',axislabelsize)
    %set(findall(gcf,'-property','FontSize'),'Fontsize',fS);
    str = ['$t =\, $' num2str(tick2*dt)]; %create textbox
    dim = [.5 tlabH .3 .3];
    annotation('textbox',dim,'String',str,'FitBoxToText','on','edgecolor','none','Interpreter','Latex','Fontsize',axislabelsize)
    %set(findall(gcf,'-property','FontSize'),'Fontsize',fS);
    str = ['$t =\, $' num2str(tick3*dt)]; %create textbox
    dim = [.81 tlabH .3 .3];
    annotation('textbox',dim,'String',str,'FitBoxToText','on','edgecolor','none','Interpreter','Latex','Fontsize',axislabelsize)
    %set(findall(gcf,'-property','FontSize'),'Fontsize',fS);
    xlabel('Space ($x$)')
    ylabel('Locust Density ($\rho$)')
    %'location','west'
    leg=legend([p(1) p(5)],'$R$-dependent','$R$-independent','Interpreter','Latex');
    pos = get(leg, 'position');
    pos(1) = pos(1)*.82; % reset x position
    pos(2) = pos(2)*.9; % reset y position
    set(leg,'position',pos); % reset legend position
    %set(leg.BoxFace,'ColorType','truecoloralpha','ColorData',uint8(255*[1;1;1;0.8]))
    %set(lgnd,'color','none')
    title('Locust Density Profiles','FontSize',titlesize)
    
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


%%% SPACETIME PLOT %%%
h3 = subplot(2,2,3);
set(h3,'Units','Inches');
pos = get(h3, 'Position'); %[0.1300    0.1100    0.3347    0.3412];% from full figure [0 0 6.5 5];
set(h3,'Position', [pos(1),pos(2), .45*width, width/4])
%set(gcf,'Position', [100 100 .45*scrsz(3) .25*scrsz(3)])%[100 100 .45*scrsz(3) .12*scrsz(3)])

    lthresh=20;
    SM = S+M;
    SM(find(SM<=lthresh))=0;
    SM(find(SM>lthresh))=1;
    SMspr = Sspr+Mspr;
    SMspr(find(SMspr<=lthresh))=0;
    SMspr(find(SMspr>lthresh))=1;

    imagesc(fliplr(SM+SMspr)')
    colormap(map)
    xl=xlim;
    xlim=([0,xl(2)]);
    set(gca,'FontSize',ticklabelsize)
    xlabel('Space ($x$)','Fontsize',axislabelsize)
    ylabel('Time ($t$)','Fontsize',axislabelsize)
    yticks([1, tick2, 5*tick1, nt]) % just hack it together
%     yticklabels({num2str(tick3*dt),num2str(tick2*dt),num2str(tick1*dt),'0'}) % in reverse order because of how Matlab indexes imagesc()
    yticklabels({'15k', '7.5k', '2.5k', '0'})
    xticks(0:nx/5:nx)
    xticklabels({0:L/5:L})
    % set(findall(gcf,'-property','FontSize'),'Fontsize',fS);
    title(['Locust Density $\rho(t)> $ ' num2str(lthresh)], 'FontSize',titlesize)
    
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
annotation('textbox', [0, 0.48, 0, 0], 'string', '\textbf{(B)}', 'FontSize', subfiglabelsize,'Interpreter','latex')


%%% MAX DENSITY PLOT %%%
h4 = subplot(2,2,4);
set(h4,'Units','Inches');
pos = get(h4, 'Position'); % [0.5703    0.1100    0.3347    0.3412];% from full figure [0 0 6.5 5];
set(h4,'Position', [pos(1),pos(2), 0.45*width, width/4])
%set(gcf,'Position', [100 100 .45*scrsz(3) .25*scrsz(3)])%[100 100 .45*scrsz(3) .12*scrsz(3)])

g = plot(   1:nt, rhomax,...
            1:nt, rhomaxSpr...
            );
    g(1).Color = RdepColor;
    g(2).Color = RindepColor;
    g(1).LineWidth = LWsmall;
    g(2).LineWidth = LWbig;
    g(1).LineStyle = LSrd;
    g(2).LineStyle = LSri;
    uistack(h(1),'top');
    ax = gca;
    ax.XLim = [0 nt];
    ax.YLim = [0 1500];
    xticks([1 tick1 tick2 tick3])
    xticklabels({0, tick1*dt, tick2*dt, tick3*dt})%({0:(nt*dt)/3:nt*dt})
    yl = ylim;
    yticks(yl(1):(yl(2)-yl(1))/2:yl(2))
    
    set(gca,'FontSize',ticklabelsize)
    ylabel('Max Density ($\mathrm{P}$ )','FontSize',axislabelsize)
    xlabel('Time ($t$)','FontSize',axislabelsize)
    title('Pulse Height','FontSize',titlesize)

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
annotation('textbox', [.45, .48, 0, 0], 'string', '\textbf{(C)}', 'FontSize', subfiglabelsize,'Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print(h,[filename, '.eps'],'-depsc')
set(h,'PaperPositionMode','Auto')

