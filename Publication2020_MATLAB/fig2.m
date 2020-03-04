close all
clear all

addpath('functions','data')

set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

titlesize = 12;
axislabelsize = 10;
fontSize = axislabelsize;
ticklabelsize = 8;
subfiglabelsize = 12;

%colors
ksm_color = 1/255*[237, 143, 33];   %gold
kms_color = 1/255*[111,73,179];     %purple

%linewidth
lwidth = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure 2 %%%
number = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_plus = 200;
alpha2 = 0.12; %alpha - not sample choice = 0.0045
eta2 = 0.005; %eta - not sample choice = 0.0036
gamma2 = 0.03;%gamma
beta2 = 0.02;%beta
theta2 = 0.14;%theta
delta2 = 0.015;%delta - not sample choice = 0.005

R_test = linspace(0,R_plus*1.25,200);
[ksm, kms] = fxn_switching(R_test, alpha2, eta2, gamma2, beta2, theta2, delta2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(number);
set(h,'Units','Inches');
pos = get(h, 'Position')
set(h,'PaperPositionMode','Manual') % Setting this to 'manual' unlinks 'Position' (what gets displayed) and 'PaperPosition' (what gets printed)
set(h,'PaperPosition',[ 0 0 4 3]);
set(h,'Position',[ 0 0 4 3]);
get(h,'Position')

psm = plot(R_test,ksm,'DisplayName','$k_{sm}$', 'Color', ksm_color);
hold on
pms = plot(R_test,kms,'DisplayName','$k_{ms}$', 'Color', kms_color);
%vline = plot([max(R_plus),max(R_plus)],[0,1.4],'k--','DisplayName','max R^+');
hold on
yline(eta2, '--k', 'DisplayName', '$\eta$')
hold on
yline(theta2, '-.k', 'DisplayName', '$\theta$')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Touch Up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca, 'FontSize', ticklabelsize)
xlabel('Resources','FontSize',axislabelsize)
ylabel('Switching per unit time','FontSize',axislabelsize)
xlim([0, R_plus*1.25])
ylim([0, 0.15])
hlegend = legend;
set(hlegend.BoxFace, 'ColorType', 'truecoloralpha', 'ColorData', uint8(255*[1;1;1;.7]));
set(hlegend, 'Location', 'east');
set(psm,'LineWidth',lwidth);
set(pms,'LineWidth',lwidth);
title('Switching functions','FontSize',titlesize)

% labels
text(R_plus*0.02, alpha2, '$\alpha$', 'FontSize',fontSize);
text(R_plus*1.15, eta2+0.0125, '$\eta$', 'FontSize',fontSize);
text(125, 0.017, '$k_{sm}$', 'FontSize',fontSize);
text(R_plus*0.02, beta2, '$\beta$', 'FontSize',fontSize);
text(R_plus*1.15, theta2-0.015, '$\theta$', 'FontSize',fontSize);
text(125, 0.11, '$k_{ms}$', 'FontSize',fontSize);

    %remove whitespace
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save Figure %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = sprintf('fig%d',number);
print(h,[filename, '.eps'],'-depsc')
