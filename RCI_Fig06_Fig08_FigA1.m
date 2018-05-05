%% RCI_Fig06_Fig08_FigA1.m
% tbeucler - 2017
% Bulk-plume toy model

clearvars; close('all'); fclose('all');

%% 0. Parameters

% Thermodynamic constants
cp = 1005.7; % Specific heat capacity of dry air at constant pressure [J/K/kg]
g = 9.8; % Gravity constant [m2/s]
Lv = 2.501e6; % Specific latent heat of vaporization [J/kg]

% Toy model parameters: Free-tropospheric radiative cooling
load('MAT_DATA/RCI_SAM_RRTM_mult.mat','QLW_mult','QSW_mult');
Qcool = QLW_mult(5,3)-QSW_mult(5,3);

% Figure characteristics
AMP = 5; % Max of linear response matrix [1/day]
fsz = 12; % Fontsize
lw = 1.5; % Linewidth
M0 = 0.1; % Normalization factor to convert LRM in log-space [1/day]

% Conversion constants
sph = 3600; % Number of seconds per hour
spd = sph*24; % Number of seconds per day

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

%% 2. RCE entrainment, detrainment, condensation rates & mass flux profile

% 2.1 Integral I on the LRM restriction domain
I = -cumtrapz(p,dDSE_dp./(Lv*qdef)); I = I(1:itp)-I(ibl);

% 2.2 Convective properties in the bulk-plume toy model
% m [kg/m2/s] e [1/s] d [1/s] c [1/s]
m = Qcool*exp(I)./(Lv*qdef(1:itp));
e = g*Qcool*exp(I).*(dDSE_dp(1:itp)+Lv*dqsat_dp(1:itp))./(Lv^2*qdef(1:itp).^2);
d = g*Qcool*exp(I).*dq_dp(1:itp)./(Lv*qdef(1:itp).^2);
c = -g*Qcool*exp(I).*dDSE_dp(1:itp)./(Lv^2*qdef(1:itp));

%% 3. Convective linear response matrices

% 3.1 Remedy the e<0 singularity
eREM = e; % Initializes the remediated entrainment eREM with e
[ism0,~]=find((e(ibl:itp)<0)-(e(ibl-1:itp-1)<0)==1);
ism0 = ism0+(ibl-2); % Find indices where e goes from positive to negative
[igr0,~]=find((e(ibl:itp)>0)-(e(ibl-1:itp-1)>0)==1);
igr0 = igr0+(ibl-2); % Find indices where e goes from negative to positive
if e(end)<0, igr0 = [igr0 itp-1]; end
% Replaces value of entrainment by value right below
for ind = 1:length(ism0), eREM(ism0(ind):igr0(ind)+1)=e(ism0(ind));
end

% 3.2 Loop over response and perturbation levels to compute LRM
Cvc_moi=zeros(Lp,Lp); Cvc_hea = Cvc_moi;
for i=ibl:itp
    for j=ibl:itp
        
        Sdc = 0; % Sum of detrainment and condensation
        for ip = (j+1):itp
            if d(ip)>0, Sdc = Sdc + d(ip)*qsat(ip)^1; end
            if c(ip)>0, Sdc = Sdc + c(ip); end
        end
        Cvc_moi(i,j) = eREM(j)*(-(i==j)+(i>j)*d(i)*qsat(i)^1*(d(i)>0)/Sdc)+...
            g*m(j)*((i==j-1)-(i==j))/(50*(p(j-1)-p(j+1)));
        Cvc_hea(i,j) = alpha(i)*eREM(j)*((i>j)*c(i)*(c(i)>0)/Sdc);
        if j==itp
            if i==j, Cvc_moi(i,j)=-eREM(j); Cvc_hea(i,j)=0;
            else Cvc_moi(i,j)=0; Cvc_hea(i,j)=0;
            end
        end
        
    end
end
LRM = spd*(Cvc_moi+Cvc_hea); % Linear response matrix 
% = sum of cvc moistening and heating [1/day]

%% 4. Paper Figure 06

figure('position',[100 100 350 700]);

% Plot convective profiles
plot(sph*e,p(1:itp),'Linewidth',lw,'color',[0.25 0 0.75]); hold on;
plot(sph*d,p(1:itp),'Linewidth',lw,'color',[1 0 0]); hold on;
plot(sph*1e3*c,p(1:itp),'Linewidth',lw,'color',[1 0.5 0],'Linestyle','--');
hold on; plot(1e3*m,p(1:itp),'Linewidth',lw,'color',[0 0 0]); hold on;
L = legend('$e\left(p\right)$','$d\left(p\right)$',...
    '$10^{3}c\left(p\right)$','$m\left(p\right)$');
set(L,'Location','NorthEast','Interpreter','Latex','Fontsize',fsz+2);

% Set axis limits
xlim([0 95]); ylim([p(itp) p(1)]);
set(gca,'Ydir','reverse');grid on;

% Bottom axis
axB=gca;
xlabel(axB,['$e\left(p\right),d\left(p\right)\ \&\ 10^{3}c\left(p\right)',...
    '\ \left[\mathrm{hour^{-1}}\right]\ $'],'Interpreter','Latex');
ylabel(axB,'p [hPa]','Interpreter','Latex');
set(axB,'Fontsize',fsz,'TickLabelInterpreter','Latex','Box','on','TickDir','out');

% Top axis
axT=axes('Position',axB.Position,'XAxisLocation','top','Color','none');
xlabel(axT,'$m\left(p\right)\ \left[\mathrm{gm^{-2}s^{-1}}\right]$',...
    'Interpreter','Latex','Fontsize',fsz);
set(axT,'Fontsize',fsz,'TickLabelInterpreter','Latex'); 
axT.XTickLabel=axB.XTickLabel; 
axT.XTickLabel=''; axT.YAxis.Visible='off'; 

% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);

% Save plot
gcfsavepdf([basedir 'PDF_DATA\RCI_Fig06.pdf']);

%% 5. Paper figure 08

F = figure('position',[100 100 800 690]);

%%% 5.1. Central plot (Linear response matrix)
SUB1 = subplot(3,4,[1 2 3 5 6 7]);

TOPLOT = sign(LRM(res,res)).*log10(1+abs(LRM(res,res)/M0)); % Logarithmic scale for LRM
P = pcolor(p(res),p(res),TOPLOT); % Plot scaled LRM in pressure-space
P.LineStyle = 'none'; % No black lines in the pcolor plot

colormap(redblue); CLIM = log10(1+AMP/M0); caxis([-CLIM CLIM]); % Colormap
c = colorbar; c.Location = 'westoutside'; % Colorbar
ylabel(c,'$day^{-1}$','Interpreter','Latex','Fontsize',fsz+1); % Colorbar's label
C = str2double(c.TickLabels); % Colorbar's ticks in log-space
c.TickLabels=num2str(M0*sign(C).*(10.^abs(C)-1),'%02.1f'); % Convert them to lin-space

xlim([min(p(res)) max(p(res))]); ylim([min(p(res)) max(p(res))]);

set(gca,'Ydir','reverse','Xdir','reverse'); % Adjust appearance of axis
axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
set(SUB1,'tickLabelInterpreter','Latex','Box','on','tickDir','out',...
    'Xtick','','Ytick','');
SUB1_P = SUB1.Position; % Save position of first sub-plot

%%% 5.2. Bottom plot (Column-integrated growth rate)
SUB2 = subplot(3,4,[10 11]);
[LRM_j,~] = LRM_growthrate(LRM',p',res); % Computes col-int growth rate
plot(p(res),LRM_j,'color',[0 0 0],'Linewidth',lw+0.5); % Plot in index-space

% Plot's appearance
grid on; xlim([min(p(res)) max(p(res))]); ylim([-2.5 11.5]);
xlabel('$p_{j}\ [hPa]$','Interpreter','Latex','Fontsize',fsz);
ylabel('$\widehat{M}_{j}\ [d^{-1}]$','Interpreter','Latex','Fontsize',fsz);
set(SUB2,'tickLabelInterpreter','Latex','Box','on','tickDir','out',...
    'Fontsize',fsz,'Xdir','reverse');
SUB2.Position(2) = SUB1_P(2)-SUB2.Position(4);
SUB2.Position(1) = SUB1_P(1); SUB2.Position(3) = SUB1_P(3);
SUB1.Position = SUB1_P;

%%% 5.3. Right-side plot (Normalized eigenvector)
SUB3 = subplot(3,4,[4 8]);

% Eigenvector modulus and maximal eigenvalue real part
[V,D]=eig(LRM(res,res));[lam,ilam]=sort(real(diag(D))); LAM = lam(end);
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

%%% 5.4 Final adjustments

% Move all the subplots down
Yshift = 0.1;
SUB3.Position(2) = SUB3.Position(2)-Yshift;
SUB2.Position(2) = SUB2.Position(2)-Yshift;
SUB1.Position(2) = SUB1.Position(2)-Yshift;

% Add leading eigenvalue real part
annotation(F,'textbox',[SUB3.Position(1)+0.1*SUB3.Position(3) ...
    SUB3.Position(2)+0.05*SUB3.Position(4) 0.11 0.04],'String',...
    ['$\lambda=$',num2str(LAM,'%03.2f'),...
    '$\mathrm{d^{-1}}$'],'LineStyle','none','Interpreter','latex',...
    'FitBoxToText','off','BackgroundColor',[1 1 1],'Fontsize',fsz+1);

%%% 5.5 Save figure in PDF format
% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
gcfsavepdf([basedir 'PDF_DATA\RCI_Fig08.pdf']);

%% 6. Paper Figure A1

F = figure('position',[100 100 800 400]);
LRM_moi = spd*Cvc_moi(res,res); LRM_hea = spd*Cvc_hea(res,res);

%%% 6.1 Convective moistening linear response matrix
SUB1=subplot(1,13,1:6);
TOPLOT = sign(LRM_moi).*log10(1+abs(LRM_moi/M0)); % Log scale for LRM
P = pcolor(p(res),p(res),TOPLOT); % Plot scaled LRM in index-space
P.LineStyle = 'none';

colormap(redblue); CLIM = log10(1+AMP/M0); caxis([-CLIM CLIM]); % Colormap
xlabel('$p_{j}\ [hPa]$','Interpreter','Latex','Fontsize',fsz+1);
ylabel('$p_{i}\ [hPa]$','Interpreter','Latex','Fontsize',fsz+1);
xlim([min(p(res)) max(p(res))]); ylim([min(p(res)) max(p(res))]); grid on;

set(gca,'Ydir','reverse','Xdir','reverse'); % Adjust appearance of axis
axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
set(SUB1,'tickLabelInterpreter','Latex','Box','on','tickDir','out');

%%% 6.2 Convective heating linear response matrix
SUB2=subplot(1,13,7:13);
TOPLOT = sign(LRM_hea).*log10(1+abs(LRM_hea/M0)); % Log scale for LRM
P = pcolor(p(res),p(res),TOPLOT); % Plot scaled LRM in index-space
P.LineStyle = 'none';

colormap(redblue); CLIM = log10(1+AMP/M0); caxis([-CLIM CLIM]); % Colormap
c=colorbar; grid on; % Colorbar
C = str2double(c.TickLabels); % Colorbar's ticks in log-space
c.TickLabels=num2str(M0*sign(C).*(10.^abs(C)-1),'%02.1f'); % Convert them to lin-space
ylabel(c,'$day^{-1}$','Interpreter','Latex','Fontsize',fsz+1); % Colorbar's label
xlabel('$p_{j}\ [hPa]$','Interpreter','Latex','Fontsize',fsz+1);
xlim([min(p(res)) max(p(res))]); ylim([min(p(res)) max(p(res))]);

set(gca,'Ydir','reverse','Xdir','reverse'); % Adjust appearance of axis
axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
set(SUB2,'tickLabelInterpreter','Latex','Box','on','tickDir','out','YTick','');

%%% 6.3 Final adjustments

% Move all the subplots up
Yshift = 0.033;
SUB2.Position(2) = SUB2.Position(2)+Yshift;
SUB1.Position(2) = SUB1.Position(2)+Yshift;

% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);

% Save plot
gcfsavepdf([basedir 'PDF_DATA\RCI_FigA1.pdf']);
