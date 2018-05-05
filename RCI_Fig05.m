%% RCI_Fig05.m
% tbeucler - 2017
% Betts-Miller toy model

clearvars; close('all'); fclose('all');

%% 0. Parameters

% Thermodynamic constants
cp = 1005.7; % Specific heat capacity of dry air at constant pressure [J/K/kg]
g = 9.8; % Gravity constant [m2/s]
Lv = 2.501e6; % Specific latent heat of vaporization [J/kg]

% Toy model parameters
tau_BM = 3/24; % Parameter: Betts-Miller relaxation time [day]

% Figure characteristics
AMP = 5; % Max of linear response matrix [1/day]
fsz = 12; % Fontsize
lw = 1.5; % Linewidth
M0 = 0.1; % Normalization factor to convert LRM in log-space [1/day]

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
for ip = 1:Lp, RH(ip) = RH_ls(QV(ip)/(1-QV(ip)),T(ip),p(ip));
end

[~,ibl]=max(RH(1:20)); ibl=ibl+1; % Top of boundary layer's index
[~,itp]=min(T); % Tropopause index
res=ibl:itp; % Index-space for linear response matrix

% 1.4 Compute energetics and alpha coefficient

DSE = cp*T+g*z; % Dry static energy [J/kg]
dDSE_dp = differentiate(DSE,p); % Vertical gradient of DSE [J/kg/hPa]
dq_dp = differentiate(QV,p); % Vert. grad. of moisture [kg/kg/hPa]
alpha = -Lv*dq_dp./dDSE_dp; % Heating to advection of moisture coefficient

%% 2. Betts-Miller linear convective response (in 1/day)

LRM=zeros(Lp,Lp); % Linear response matrix
for ires=ibl:itp
    for iper=ibl:itp
        LRM(ires,iper) = alpha(ires)/tau_BM*0.5*(p(ires+1)-p(ires-1))/...
            (p(res(end)+1)-p(res(1))); % Convective heating
        if ires==iper
            LRM(ires,iper)=LRM(ires,iper)-1/tau_BM; % Convective drying
        end
    end
end

%% 3. Paper figure 05

F = figure('position',[100 100 800 690]);

%%% 3.1. Central plot (Linear response matrix)
SUB1 = subplot(3,4,[1 2 3 5 6 7]);

TOPLOT = sign(LRM(res,res)).*log10(1+abs(LRM(res,res)/M0)); % Logarithmic scale for LRM
P = pcolor(p(res),p(res),TOPLOT); % Plot scaled LRM in pressure-space
P.LineStyle = 'none'; % No black contours for pcolor plot

colormap(redblue); CLIM = log10(1+AMP/M0); caxis([-CLIM CLIM]); % Colormap
c=colorbar; c.Location = 'westoutside'; % Colorbar
ylabel(c,'$day^{-1}$','Interpreter','Latex','Fontsize',fsz+1); % Colorbar's label
C = str2double(c.TickLabels); % Colorbar's ticks in log-space
c.TickLabels=num2str(M0*sign(C).*(10.^abs(C)-1),'%02.1f'); % Convert them to lin-space

xlim([min(p(res)) max(p(res))]); ylim([min(p(res)) max(p(res))]);

set(gca,'Ydir','reverse','Xdir','reverse'); % Adjust appearance of axis
axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
set(SUB1,'tickLabelInterpreter','Latex','Box','on','tickDir','out',...
    'Xtick','','Ytick','');
SUB1_P = SUB1.Position; % Save position of first sub-plot

%%% 3.2. Bottom plot (Column-integrated growth rate)
SUB2 = subplot(3,4,[10 11]);
[LRM_j,~] = LRM_growthrate(LRM',p',res); % Computes col-int growth rate
plot(p(res),LRM_j,'color',[0 0 0],'Linewidth',lw+0.5); % Plot in pressure-space

% Plot's appearance
grid on; xlim([min(p(res)) max(p(res))]); ylim([-0.5 0.1]);
xlabel('$p_{j}\ [hPa]$','Interpreter','Latex','Fontsize',fsz);
ylabel('$\widehat{M}_{j}\ [d^{-1}]$','Interpreter','Latex','Fontsize',fsz);
set(SUB2,'tickLabelInterpreter','Latex','Box','on','tickDir','out',...
    'Fontsize',fsz,'Xdir','reverse');
SUB2.Position(2) = SUB1_P(2)-SUB2.Position(4);
SUB2.Position(1) = SUB1_P(1); SUB2.Position(3) = SUB1_P(3);
SUB1.Position = SUB1_P;

%%% 3.3. Right-side plot (Normalized eigenvector)
SUB3 = subplot(3,4,[4 8]);

% Eigenvector modulus and maximal eigenvalue real part
[V,D]=eig(LRM(res,res)); [lam,ilam]=sort(real(diag(D))); LAM = lam(end);
EIGV = ((real(V(:,ilam(end)))).^2+(imag(V(:,ilam(end)))).^2).^0.5;
norm = 1/max(EIGV);

plot(norm*EIGV,p(res),'Linewidth',lw+0.5,'color',[0 0 0]); hold on;

% Plot's appearance
set(SUB3,'TickLabelInterpreter','Latex','Box','on','TickDir','out',...
    'Fontsize',fsz,'Ydir','reverse');
SUB3.YAxisLocation='right'; SUB3.XAxisLocation = 'top';
grid on; xlim([min(norm*EIGV) max(norm*EIGV)]); ylim([min(p(res)) max(p(res))]);
xlabel('Normalized eigenvector','Interpreter','Latex','Fontsize',fsz);
ylabel('$p_{i}\ [hPa]$','Interpreter','Latex','Fontsize',fsz);
SUB3.Position(1) = SUB1_P(1)+SUB1_P(3);
SUB3.Position(2) = SUB1_P(2); SUB3.Position(4) = SUB1_P(4);
SUB1.Position = SUB1_P;

%%% 3.4 Final adjustments

% Move all the subplots down
Yshift = 0.1;
SUB3.Position(2) = SUB3.Position(2)-Yshift;
SUB2.Position(2) = SUB2.Position(2)-Yshift;
SUB1.Position(2) = SUB1.Position(2)-Yshift;

% Add leading eigenvalue real part
annotation(F,'textbox',[SUB3.Position(1)+0.1*SUB3.Position(3) ...
    SUB3.Position(2)+0.9*SUB3.Position(4) 0.11 0.04],'String',...
    ['$\lambda=$',num2str(LAM,'%03.2f'),...
    '$\mathrm{d^{-1}}$'],'LineStyle','none','Interpreter','latex',...
    'FitBoxToText','off','BackgroundColor',[1 1 1],'Fontsize',fsz+1);

%%% 3.5 Save figure in PDF format
% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
gcfsavepdf([basedir 'PDF_DATA\RCI_Fig05.pdf']);
