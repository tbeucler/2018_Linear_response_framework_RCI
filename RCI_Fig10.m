%% RCI_Fig10
% tbeucler - 2017
% Longwave and shortwave linear response matrices derived from RRTM
% simulations

close all; clearvars; fclose('all');

%% 0. Parameters
% Conversion constants
spd = 24*3600; % Number of seconds per day

% Figure characteristics
AMP = 0.15; % Max of linear response matrix [1/day]
fsz = 12; % Fontsize
lw = 1.5; % Linewidth
M0 = 0.01; % Normalization factor to convert LRM in log-space [1/day]

%% 1. Load RCE profiles and LRM

% 1.1 Load RCE thermodynamic profiles from SAM300K simulation
L = load('MAT_DATA/SAM300K_96x96x64_mean_profile.mat');
Lp = numel(L.p); % Size of pressure vector

% 1.2 Interpolate to equal pressure level for eigenvalue analysis
p = linspace(L.p(1),L.p(end),Lp)'; % Interpolated pressure [hPa]
QV = interp1(L.p,L.QV,p,'pchip'); % Specific humidity [kg/kg]
T = interp1(L.p,L.T,p,'pchip'); % Temperature [K]
z = interp1(L.p,L.z,p,'pchip'); % Geopotential height [m]

% 1.3 Load linear response matrices derived from RRTM
load('MAT_DATA/RCI_SAM_RRTM_mult.mat','LRM_LW_mult','LRM_SW_mult');
LRM_LW = squeeze(LRM_LW_mult(5,3,:,:)); % LW linear response matrix
LRM_SW = squeeze(LRM_SW_mult(5,3,:,:)); % SW linear response matrix

% 1.4 Loads restriction for LRM
load('MAT_DATA/RCI_SAM_convection_mult.mat','ibl_mult','itp_mult');
res = ibl_mult(5,3):itp_mult(5,3);

%% 2. Paper Figure 08
F = figure('position',[100 100 800 400]);

%%% 2.1 Longwave linear response matrix
SUB1=subplot(3,2,[1 3]);

TOPLOT = sign(LRM_LW(res,res)).*log10(1+abs(spd*LRM_LW(res,res)/M0)); % Log scale
P = pcolor(p(res),p(res),TOPLOT); % Plot scaled LRM in index-space
P.LineStyle = 'none'; % No black lines in pcolor plot
colormap(redblue); CLIM = log10(1+AMP/M0); caxis([-CLIM CLIM]); % Colormap

% Axis and appearance
xlim([min(p(res)) max(p(res))]); ylim([min(p(res)) max(p(res))]);
ylabel('$p_{i}\ [hPa]$','Interpreter','Latex','Fontsize',fsz);
set(gca,'Ydir','reverse','Xdir','reverse'); % Adjust appearance of axis
axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
set(SUB1,'tickLabelInterpreter','Latex','Box','on','tickDir','out','XTick','');
SUB1_P = SUB1.Position; % Save position of first sub-plot

%%% 2.2 Longwave column-integrated growth-rate
SUB2=subplot(3,2,5);

[GR_LW,~] = LRM_growthrate(spd*LRM_LW',p',res);
plot(p(res),GR_LW,'color',[0 0 0],'Linewidth',lw+1);

% Axis and appearance
grid on; xlim([min(p(res)) max(p(res))]); ylim([-0.45 0.7]);
xlabel('$p_{j}\ [hPa]$','Interpreter','Latex','Fontsize',fsz);
ylabel('$\widehat{M}_{j}\ [d^{-1}]$','Interpreter','Latex','Fontsize',fsz);
set(SUB2,'tickLabelInterpreter','Latex','Box','on','tickDir','out',...
    'Fontsize',fsz,'Xdir','reverse');
SUB2.Position(2) = SUB1_P(2)-SUB2.Position(4);
SUB2.Position(1) = SUB1_P(1); SUB2.Position(3) = SUB1_P(3);
SUB1.Position = SUB1_P; SUB1.FontSize = fsz;

%%% 2.3 Shortwave linear response matrix
SUB3=subplot(3,2,[2 4]);

TOPLOT = sign(LRM_SW(res,res)).*log10(1+abs(spd*LRM_SW(res,res)/M0)); % Log scale
P = pcolor(p(res),p(res),TOPLOT); % Plot scaled LRM in index-space
P.LineStyle = 'none'; % No black line in pcolor plot
colormap(redblue); CLIM = log10(1+AMP/M0); caxis([-CLIM CLIM]); % Colormap

c = colorbar; % Colorbar 
ylabel(c,'$day^{-1}$','Interpreter','Latex','Fontsize',fsz+1); % Colorbar's label
C = str2double(c.TickLabels); % Colorbar's ticks in log-space
c.TickLabels=num2str(M0*sign(C).*(10.^abs(C)-1),'%02.2f'); % Convert them to lin-space
c.TickLabelInterpreter = 'Latex';

% Axis and appearance
xlim([min(p(res)) max(p(res))]); ylim([min(p(res)) max(p(res))]);
set(gca,'Ydir','reverse','Xdir','reverse'); % Adjust appearance of axis
axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
set(SUB3,'tickLabelInterpreter','Latex','Box','on','tickDir','out',...
    'XTick','','YTick','');
SUB3_P = SUB3.Position; % Save position of first sub-plot

%%% 2.4 Shortwave column-integrated growth-rate
SUB4=subplot(3,2,6);

[GR_SW,~] = LRM_growthrate(spd*LRM_SW',p',res);
plot(p(res),GR_SW,'color',[0 0 0],'Linewidth',lw+1);

% Axis and appearance
grid on; xlim([min(p(res)) max(p(res))]); ylim([-0.06 0.06]);
xlabel('$p_{j}\ [hPa]$','Interpreter','Latex','Fontsize',fsz);
ylabel('$\widehat{M}_{j}\ [d^{-1}]$','Interpreter','Latex','Fontsize',fsz);
set(SUB4,'tickLabelInterpreter','Latex','Box','on','tickDir','out',...
    'Fontsize',fsz,'YAxisLocation','right','Xdir','reverse');
SUB4.Position(2) = SUB3_P(2)-SUB4.Position(4);
SUB4.Position(1) = SUB3_P(1); SUB4.Position(3) = SUB3_P(3);
SUB3.Position = SUB3_P; SUB3.FontSize = fsz;

%%% 2.5 Final adjustments

% Shift right panel to the left and colorbar to the right
hold on; xshift = 0.05;
c.Position(1) = c.Position(1)+xshift;
SUB3.Position(1)=SUB3.Position(1)-xshift;
SUB4.Position(1)=SUB4.Position(1)-xshift;
SUB3.Position(3)=SUB1.Position(3);
SUB4.Position(3)=SUB2.Position(3);

% Leading eigenvalue real part
[~,D]=eig(spd*LRM_LW(res,res));[lam,~]=sort(real(diag(D))); LAM_LW = lam(end);
[~,D]=eig(spd*LRM_SW(res,res));[lam,~]=sort(real(diag(D))); LAM_SW = lam(end);
annotation(F,'textbox',[SUB2.Position(1)+0.1*SUB2.Position(3)...
    SUB2.Position(2)+0.05*SUB2.Position(4) 0.1 0.06],'String',...
    ['$\lambda=$',num2str(LAM_LW,'%03.2f'),'$\mathrm{d^{-1}}$'],...
    'LineStyle','none','Interpreter','latex',...
    'FitBoxToText','off','BackgroundColor',[1 1 1],'Fontsize',fsz+1);
annotation(F,'textbox',[SUB4.Position(1)+0.1*SUB4.Position(3)...
    SUB4.Position(2)+0.05*SUB4.Position(4) 0.1 0.06],'String',...
    ['$\lambda=$',num2str(LAM_SW,'%03.2f'),'$\mathrm{d^{-1}}$'],...
    'LineStyle','none','Interpreter','latex',...
    'FitBoxToText','off','BackgroundColor',[1 1 1],'Fontsize',fsz+1);

% Panel identifiers
annotation(F,'textbox',[SUB1.Position(1)-0.2*SUB1.Position(3)...
    SUB1.Position(2)+0.9*SUB1.Position(4) 0.025 0.075],'String','a)',...
    'LineStyle','none','Interpreter','latex',...
    'FitBoxToText','off','BackgroundColor',[1 1 1],'Fontsize',fsz+1);
annotation(F,'textbox',[SUB3.Position(1)-0.1*SUB3.Position(3)...
    SUB3.Position(2)+0.9*SUB3.Position(4) 0.025 0.075],'String','b)',...
    'LineStyle','none','Interpreter','latex',...
    'FitBoxToText','off','BackgroundColor',[1 1 1],'Fontsize',fsz+1);

%%% 2.6 Save plot
% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'PDF_DATA\RCI_Fig10.pdf']);
