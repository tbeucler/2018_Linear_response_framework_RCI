%% RCI_Fig04.m
% tbeucler - 4/11/2018
% Plots the q,T,alpha,O3 profiles used as base RCE profiles
% The thermodynamic profiles are ensemble-mean of square SAM runs
% from Tristan Abbott with a sea surface temperatures of
% 280-305K with an interval of 5K

close all; fclose('all'); clearvars;

%% 0. Parameters

% Thermodynamic constants
cp = 1005.7; % Specific heat capacity of dry air at constant pressure [J/K/kg]
g = 9.8; % Gravity constant [m2/s]
Lv = 2.501e6; % Specific latent heat of vaporization [J/kg]
sig = 5.67e-8; % Stefan-Boltzmann constant (W/m^2/K^4)

% Figure characteristics
fsz = 12; % Fontsize
lw = 1.5; % Linewidth
Cmap = [0.36 0.36 1;0.27 0.27 0.75; 0.18 0.18 0.5;0.5 0 0;...
    0.75 0 0;1 0 0]; % Colormap from cold to warm

%% 1. Prepare thermodynamic profiles

% 1.1 Load the different base states
load('MAT_DATA/RCI_SAM_convection_mult.mat','alpha_mult','QV_mult');

% 1.2 Compute temperature of different base states
T_mult = zeros(6,64); % Initializes T for multiple temperatures
SST_array = linspace(280,305,6); % Surface temperature space
for iSST = 1:6, SST = SST_array(iSST);
    load(['MAT_DATA/RCI',num2str(SST),'K_Modified_base_state.mat'],'T','z','p');
    T_mult(iSST,:) = T;
end

% 1.3 Load Ozone profile (ppmv)
load('MAT_DATA/SAM_O3_ppmv.mat'); % Load Ozone profile (in ppmv)

%% 2. Paper Figure 04
F = figure('position',[100 100 1000 400]);

%%% 2.1 Plot specific humidity profile
SUB1 = subplot(1,4,1);
for iSST = 1:6, P = plot(1e3*squeeze(QV_mult(iSST,3,:)),p,'Linewidth',lw,'color',...
        Cmap(iSST,:)); grid on; hold on;
    if iSST==5, P.LineWidth = 2*lw; end
end
LEG = legend('280K','285K','290K','295K','300K','305K');
set(LEG,'Fontsize',fsz,'Location','Northeast','Interpreter','Latex');
xlim([0 1.1*max(1e3*QV_mult(end,3,:))]); ylim([0 1000]);
xlabel('Specific humidity [g/kg]','Interpreter','Latex','Fontsize',fsz);
ylabel('Pressure [hPa]','Interpreter','Latex','Fontsize',fsz);
set(SUB1,'Ydir','reverse','TickLabelInterpreter','Latex','Box','on','TickDir','out');

%%% 2.2 Plot temperature profile
SUB2 = subplot(1,4,2);
for iSST = 1:6, P = plot(T_mult(iSST,:),p,'Linewidth',lw,'color',...
        Cmap(iSST,:)); grid on; hold on;
    if iSST==5, P.LineWidth = 2*lw; end
end
xlim([0.9*min(T_mult(1,:)) 1.02*max(T_mult(end,:))]); ylim([0 1000]);
xlabel('Temperature [K]','Interpreter','Latex','Fontsize',fsz);
set(SUB2,'Ydir','reverse','TickLabelInterpreter','Latex','Box','on','TickDir','out');

%%% 2.3 Plot alpha profile
SUB3 = subplot(1,4,3);
TOPLOT = squeeze(alpha_mult(:,3,:)); TOPLOT(TOPLOT<0) = 0; % Only plot al>0
for iSST = 1:6, P = plot(TOPLOT(iSST,2:end),p(2:end),'Linewidth',lw,...
        'color',Cmap(iSST,:)); grid on; hold on;
    if iSST==5, P.LineWidth = 2*lw; end
end
line([1 1],[0 1000],'Linewidth',lw,'color','k','Linestyle','--');
xlim([0 5]); ylim([0 1000]);
xlabel('$\alpha $','Interpreter','Latex','Fontsize',fsz);
set(SUB3,'Ydir','reverse','TickLabelInterpreter','Latex','Box','on','TickDir','out');

%%% 2.4 Plot ozone profile
SUB4 = subplot(1,4,4);
plot(log10(1e6*pO3),p,'Linewidth',lw,'color','k'); grid on; hold on;
line([0 0],[0 1000],'Linewidth',lw,'color','k','Linestyle','--');
ylim([0 1000]);
xlabel(['$\mathrm{log_{10}}\left(p_{\mathrm{O3}}\ \left[',...
    '\mathrm{ppmv}\right]\right)$'],'Interpreter','Latex','Fontsize',fsz);
set(SUB4,'Ydir','reverse','TickLabelInterpreter','Latex','Box','on','TickDir','out');

%%% 2.5 Panel indicators
annotation(F,'textbox',[SUB1.Position(1)-0.3*SUB1.Position(3) ...
    SUB1.Position(2)+0.99*SUB1.Position(4) 0.02 0.05],'String',...
    'a)','LineStyle','none','Interpreter','latex','FitBoxToText','off',...
    'Backgroundcolor',[1 1 1],'Fontsize',fsz+1,'color','k');
annotation(F,'textbox',[SUB2.Position(1)-0.3*SUB2.Position(3) ...
    SUB2.Position(2)+0.99*SUB2.Position(4) 0.02 0.05],'String',...
    'b)','LineStyle','none','Interpreter','latex','FitBoxToText','off',...
    'Backgroundcolor',[1 1 1],'Fontsize',fsz+1,'color','k');
annotation(F,'textbox',[SUB3.Position(1)-0.3*SUB3.Position(3) ...
    SUB3.Position(2)+0.99*SUB3.Position(4) 0.02 0.05],'String',...
    'c)','LineStyle','none','Interpreter','latex','FitBoxToText','off',...
    'Backgroundcolor',[1 1 1],'Fontsize',fsz+1,'color','k');
annotation(F,'textbox',[SUB4.Position(1)-0.3*SUB4.Position(3) ...
    SUB4.Position(2)+0.99*SUB4.Position(4) 0.02 0.05],'String',...
    'd)','LineStyle','none','Interpreter','latex','FitBoxToText','off',...
    'Backgroundcolor',[1 1 1],'Fontsize',fsz+1,'color','k');

%%% 2.5 Save plot
% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'PDF_DATA\RCI_Fig04.pdf']);