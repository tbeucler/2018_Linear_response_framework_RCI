%% RCI_Fig09_FigA2
% tbeucler - 2017
% Gray radiation 2-stream Schwarzschild toy-model
% Uses CRSRCE_LWgray.m and CRSRCE_SWgray.m

clearvars; close('all'); fclose('all');

%% 0. Parameters

% Thermodynamic constants
cp = 1005.7; % Specific heat capacity of dry air at constant pressure [J/K/kg]
g = 9.8; % Gravity constant [m2/s]
Lv = 2.501e6; % Specific latent heat of vaporization [J/kg]
sig = 5.67e-8; % Stefan-Boltzmann constant (W/m^2/K^4)

% Conversion constants
spd = 24*3600; % Number of seconds per day

% Figure characteristics
AMP = 0.15; % Max of linear response matrix [1/day]
fsz = 12; % Fontsize
lw = 1.5; % Linewidth
M0 = 0.01; % Normalization factor to convert LRM in log-space [1/day]

% Simulation characteristics
Ts = 300; % Sea surface temperature of the SAM simulation [K]

% Toy model parameters: Free-tropospheric radiative cooling
load('MAT_DATA/RCI_SAM_RRTM_mult.mat','QLW_mult','QSW_mult'); % RRTM file
Qcool = QLW_mult(5,3)-QSW_mult(5,3); % Net radiative cooling (from RRTM)
N_kappaL = 100; % Number of points in longwave absorption space
kappaL_space = 10.^linspace(-2.5,-0.715,N_kappaL); % Longwave absorption space
eps = 0.077; % Ratio of shortwave to longwave absorption
kappaS_space = eps*kappaL_space; % Shortwave absorption space
pbd = 1; % Pressure broadening coefficient for SW radiation calculation [1]
ps = 1e3; % Surface pressure [hPa]

%% 1. Compute properties of RCE

% 1.1 Load RCE thermodynamic profiles from SAM300K simulation
L = load('MAT_DATA/SAM300K_96x96x64_mean_profile.mat');
Lp = numel(L.p); % Size of pressure vector

% 1.2 Interpolate to equal pressure level for eigenvalue analysis
p = linspace(L.p(1),L.p(end),Lp)'; % Interpolated pressure [hPa]
QV = interp1(L.p,L.QV,p,'pchip'); % Specific humidity [kg/kg]
T = interp1(L.p,L.T,p,'pchip'); % Temperature [K]
z = interp1(L.p,L.z,p,'pchip'); % Geopotential height [m]

% 1.3 Compute relative humidity, saturation specific humidity & deficit

RH = zeros(Lp,1); % Relative humidity taking into account liq/sol
% saturation pressure
qsat = RH; % Saturation specific humidity [kg/kg]
for ip = 1:Lp,
    RH(ip) = RH_ls(QV(ip)/(1-QV(ip)),T(ip),p(ip));
    qsat(ip) = q_sat(p(ip),T(ip));
end
qdef = qsat-QV; % Saturation deficit [kg/kg]

[~,ibl]=max(RH(1:20)); ibl=ibl+1; % Top of boundary layer's index
[~,itp]=min(T); itp = itp-5; % Tropopause index
res=ibl:itp; % Index-space for linear response matrix

% 1.4 Compute energetics and alpha coefficient

DSE = cp*T+g*z; % Dry static energy [J/kg]
dDSE_dp = differentiate(DSE,p); % Vertical gradient of DSE [J/kg/hPa]
dq_dp = differentiate(QV,p); % Vert. grad. of moisture [kg/kg/hPa]
dqsat_dp = differentiate(qsat,p); % Vert. grad. of moisture [kg/kg/hPa]
alpha = -Lv*dq_dp./dDSE_dp; % Heating to advection of moisture coefficient

%% 2. Compute LW gray radiative fluxes for kappa_L between 1e-3 and 1 m^2/kg

%%% 2.1 Span absorption space to find the different LW & SW heatings
QLW_space = zeros(N_kappaL,1); tauLWs = QLW_space;
QSW_space = zeros(N_kappaL,1); tauSWs = QSW_space;
for ikappa = 1:N_kappaL, disp([num2str(ikappa),'/',num2str(N_kappaL)]);
    [ ~,FupLW,FdoLW,tauLW ] = CRSRCE_LWgray( p,QV,T,Ts,kappaL_space(ikappa) );
    [ ~,FdoSW,tauSW ] = CRSRCE_SWgray( p,QV,kappaS_space(ikappa),L.S0,pbd );
    QLW_space(ikappa) = FupLW(end)+FdoLW(1)-FupLW(1)-FdoLW(end);
    QSW_space(ikappa) = FdoSW(end)-FdoSW(1);
    tauLWs(ikappa) = tauLW(1); tauSWs(ikappa) = tauSW(1);
end

%%% 2.2 Fit the total radiative cooling to the gray model
[~,imax] = max(QLW_space-QSW_space);
% Optically thin best fit
[~,ithin] = min(abs(Qcool-QLW_space(1:imax)+QSW_space(1:imax)));
% Optically thick best fit
[~,ithick] = min(abs(Qcool-QLW_space(imax+1:end)+QSW_space(imax+1:end)));
ithick = ithick+imax; % Convert back to original index space

%%% 2.3 Save radiative characteristics of best fits
[ QLWthin,FupLWthin,FdoLWthin,tauLWthin ] = ...
    CRSRCE_LWgray( p,QV,T,Ts,kappaL_space(ithin) );
[ QSWthin,FdoSWthin,tauSWthin ] = ...
    CRSRCE_SWgray( p,QV,kappaS_space(ithin),L.S0,pbd );
[ QLWthick,FupLWthick,FdoLWthick,tauLWthick ] = ...
    CRSRCE_LWgray( p,QV,T,Ts,kappaL_space(ithick) );
[ QSWthick,FdoSWthick,tauSWthick ] = ...
    CRSRCE_SWgray( p,QV,kappaS_space(ithick),L.S0,pbd );

%% 3. Paper Figure A2
F = figure('position',[100 100 800 400]);

%%% 3.1 Plot the radiative heating as a function of surf. LW opt. thickness
SUB1 = subplot(1,3,1:2);
plot(tauLWs,QLW_space-QSW_space,'Linewidth',lw,'color','k'); hold on;
plot(tauLWs,QLW_space,'Linewidth',lw,'color',[1 0.25 0.25]); hold on;
plot(tauLWs,QSW_space,'Linewidth',lw,'color',[0.25 0.25 1]); hold on;

% Indicate best fit using dotted line
line([min(tauLWs) max(tauLWs)],[Qcool Qcool],'color','k','Linestyle','--');
hold on; line([tauLWs(ithin) tauLWs(ithin)],[0 1e3*Qcool],...
    'color',[0 0.5 0],'Linestyle','--'); hold on;
line([tauLWs(ithick) tauLWs(ithick)],[0 1e3*Qcool],'color',[1 0.5 0],...
    'Linestyle','--'); hold on;

% General plot characteristics
xlim([min(tauLWs) max(tauLWs)]); ylim([0 1.1*max(QLW_space)]); grid on;
xlabel('$\mathrm{Surface\ optical\ thickness}$',...
    'Fontsize',fsz,'Interpreter','Latex');
ylabel('$\mathrm{Radiative\ cooling\ [W.m^{-2}]}$','Fontsize',fsz,'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex','Box','on','TickDir','out');

%%% 3.2 Plot the optical thicknes profiles for the two best fits
SUB2 = subplot(1,3,3);
plot(log10(tauLWthin),p,'color',[0 0.5 0],'Linewidth',lw); hold on;
plot(log10(tauLWthick),p,'color',[1 0.5 0],'Linewidth',lw); hold on;

% General plot characteristics
ylim([0 1000]); xlim([1.1*min(log10(tauLWthin)) 1.1*max(log10(tauLWthick))]); 
grid on; SUB2.YDir = 'reverse';SUB2.YAxisLocation = 'right';
set(gca,'TickLabelInterpreter','Latex','Box','on','TickDir','out');
xlabel('$\log_{10}\left(\mathrm{Optical\ thickness}\right)$','Fontsize',fsz,'Interpreter','Latex');
ylabel('$\mathrm{Pressure\ [hPa]}$','Fontsize',fsz,'Interpreter','Latex');

%%% 3.3 Final adjustments
annotation(F,'textbox',[SUB1.Position(1)+0.2*SUB1.Position(3) ...
    SUB1.Position(2)+0.95*SUB1.Position(4) 0.11 0.04],'String',...
    '$\mathrm{Longwave\ cooling}$','LineStyle','none','Interpreter','latex',...
    'FitBoxToText','off','Fontsize',fsz+1,'color',[1 0.25 0.25]);
annotation(F,'textbox',[SUB1.Position(1)+0.2*SUB1.Position(3) ...
    SUB1.Position(2)+0.675*SUB1.Position(4) 0.11 0.04],'String',...
    '$\mathrm{Total\ cooling}$','LineStyle','none','Interpreter','latex',...
    'FitBoxToText','off','Fontsize',fsz+1,'color','k');
annotation(F,'textbox',[SUB1.Position(1)+0.2*SUB1.Position(3) ...
    SUB1.Position(2)+0.1*SUB1.Position(4) 0.11 0.04],'String',...
    '$\mathrm{Shortwave\ heating}$','LineStyle','none','Interpreter','latex',...
    'FitBoxToText','off','Fontsize',fsz+1,'color',[0.25 0.25 1]);
annotation(F,'textbox',[SUB1.Position(1)-0.1*SUB1.Position(3) ...
    SUB1.Position(2)+0.97*SUB1.Position(4) 0.04 0.07],'String',...
    'a)','LineStyle','none','Interpreter','latex','FitBoxToText','off',...
    'Backgroundcolor',[1 1 1],'Fontsize',fsz+1,'color','k');
annotation(F,'textbox',[SUB1.Position(1)+1.05*SUB1.Position(3) ...
    SUB1.Position(2)+0.97*SUB1.Position(4) 0.04 0.07],'String',...
    'b)','LineStyle','none','Interpreter','latex','FitBoxToText','off',...
    'Backgroundcolor',[1 1 1],'Fontsize',fsz+1,'color','k');

%%% 3.4 Save plot
% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'PDF_DATA\RCI_FigA2.pdf']);

%% 4. Compute LW and SW linear response matrices using the
% best fit to the total radiative cooling

iabs = ithick; % Choose which absorption index to use for the LRM calculation
LRM_LW = zeros(Lp,Lp); LRM_SW = LRM_LW; % Initialize LRM

for iper = ibl:itp % Loops over perturbation level
    for ires = ibl:itp % Loops over response level
        
        delp = -50*(p(iper+1)-p(iper-1)); % Layer's pressure thickness
        pre = alpha(ires)*kappaL_space(iabs)/Lv*(p(iper)/ps)^pbd; % Pre-factor for both LRM
        
        % LW LRM
        T1 = (FupLWthick(iper)+FdoLWthick(iper)-2*sig*T(iper)^4)*(iper==ires);
        T2 = -FupLWthick(iper)*(iper<ires)*exp(-abs(tauLWthick(iper)-tauLWthick(ires)));
        T3 = sig*T(iper)^4*exp(-abs(tauLWthick(ires)-tauLWthick(iper)))*(iper~=ires);
        T4 = -FdoLWthick(iper)*(iper>ires)*exp(-abs(tauLWthick(ires)-tauLWthick(iper)));
        LRM_LW(ires,iper) = pre*(T1+QV(ires)*kappaL_space(iabs)*delp/g*(p(ires)/ps)^pbd*(T2+T3+T4));
        
        % SW LRM
        T1 = FdoSWthick(iper)*(iper==ires);
        T2 = -L.S0*delp/g*kappaS_space(iabs)*QV(ires)*(p(ires)/ps)^pbd*exp(-tauSWthick(ires))*(iper>ires);
        LRM_SW(ires,iper) = pre*eps*(T1+T2);
        
    end
end

%% 5. Paper Figure 09
F = figure('position',[100 100 800 400]);

%%% 5.1 Longwave linear response matrix
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

%%% 5.2 Longwave column-integrated growth-rate
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

%%% 5.3 Shortwave linear response matrix
SUB3=subplot(3,2,[2 4]);

TOPLOT = sign(LRM_SW(res,res)).*log10(1+abs(spd*LRM_SW(res,res)/M0)); % Log scale
P = pcolor(p(res),p(res),TOPLOT); % Plot scaled LRM in index-space
P.LineStyle = 'none'; % No black lines for pcolor plot
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

%%% 5.4 Shortwave column-integrated growth-rate
SUB4=subplot(3,2,6);

[GR_SW,~] = LRM_growthrate(spd*LRM_SW',p',res);
plot(p(res),GR_SW,'color',[0 0 0],'Linewidth',lw+1);

% Axis and appearance
grid on; xlim([min(p(res)) max(p(res))]); ylim([-0.4 0.4]);
xlabel('$p_{j}\ [hPa]$','Interpreter','Latex','Fontsize',fsz);
ylabel('$\widehat{M}_{j}\ [d^{-1}]$','Interpreter','Latex','Fontsize',fsz);
set(SUB4,'tickLabelInterpreter','Latex','Box','on','tickDir','out',...
    'Fontsize',fsz,'YAxisLocation','right','Xdir','reverse');
SUB4.Position(2) = SUB3_P(2)-SUB4.Position(4);
SUB4.Position(1) = SUB3_P(1); SUB4.Position(3) = SUB3_P(3);
SUB3.Position = SUB3_P; SUB3.FontSize = fsz;

%%% 5.5 Final adjustments

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

%%% 5.6 Save plot
% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'PDF_DATA\RCI_Fig09.pdf']);
