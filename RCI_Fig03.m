%% RCI_Fig03.m
% tbeucler - 3/23/2018
% Loads linear response matrices from Kuang (2012)
% Uses the following scripts:
% differentiate.m
% gcfsavepdf.m
% LRM_growthrate.m
% RH_ls.m

close all; fclose('all'); clearvars;

%% 0. Parameters

% Thermodynamic constants
cp = 1005.7; % Specific heat capacity of dry air at constant pressure [J/K/kg]
g = 9.8; % Gravity constant [m2/s]
Lv = 2.501e6; % Specific latent heat of vaporization [J/kg]

% Figure characteristics
AMP = 5; % Max of linear response matrix [1/day]
fsz = 12; % Fontsize
lw = 1.5; % Linewidth
M0 = 0.1; % Normalization factor to convert LRM in log-space [1/day]
res = 1:18; % Index space

%% 1. Load linear matrices from Kuang, 2012
% The data is not provided here but is available from Z. Kuang
Kuang_data = 'MAT_DATA/RCI_Kuang_2012.mat'; % Z. Kuang's reduced data 
if exist(Kuang_data,'file') == 2, load(Kuang_data);

%% 2. Compute alpha = -Lv*dq/dp / ds/dp
s = cp*t+g*z'; % Dry static energy [J/kg]
ds_dp = differentiate(s',p); % Vertical gradient of DSE [J/kg/hPa]
dq_dp = differentiate(1e-3*q',p); % Vert. grad. of moisture [kg/kg/hPa]
alpha = -Lv*dq_dp./ds_dp; % Heating to advection of moisture coefficient [1]

%% 3. Compute convective heating linear response matrix using alpha
% and convective heating from moisture perturbations [K/day per kg/kg]
Cvc_heat = cp/(Lv*1e-3)*dT_dt_from_q_square.*repmat(alpha(1:18),1,18);

%% 4. Paper figure 03
F = figure('position',[100 100 800 690]);
LRM = Cvc_heat+dq_dt_from_q; % Total linear response matrix [1/day]

%%% 4.1. Central plot (Linear response matrix)
SUB1 = subplot(3,4,[1 2 3 5 6 7]);

TOPLOT = sign(LRM).*log10(1+abs(LRM/M0)); % Logarithmic scale for LRM
P = pcolor(p(res),p(res),TOPLOT); % Plot scaled LRM in pressure-space
P.LineStyle = 'none'; % No black lines on pcolor plot

colormap(redblue); CLIM = log10(1+AMP/M0); caxis([-CLIM CLIM]); % Colormap
c=colorbar; c.Location = 'westoutside'; % Colorbar
ylabel(c,'$day^{-1}$','Interpreter','Latex','Fontsize',fsz+1); % Colorbar's label
C = str2double(c.TickLabels); % Colorbar's ticks in log-space
c.TickLabels=num2str(M0*sign(C).*(10.^abs(C)-1),'%02.1f'); % Convert them to lin-space

% Find and mark planetary boundary layer based on relative humidity profile
RH = zeros(numel(p),1); % Relative humidity taking into account liq/sol
% saturation pressure
for ip = 1:numel(p)
    RH(ip) = RH_ls(1e-3*q(ip)/(1-1e-3*q(ip)),t(ip),p(ip)); % RH profile
end
[~,ibl]=max(RH(res)); ibl=ibl+1; % Top of boundary layer's index
line([p(ibl) p(ibl)],[p(min(res)) p(max(res))],'color','k','Linewidth',lw,...
    'Linestyle','--'); hold on;
line([p(min(res)) p(max(res))],[p(ibl) p(ibl)],'color','k','Linewidth',lw,...
    'Linestyle','--'); hold on;
xlim([min(p(res)) max(p(res))]); ylim([min(p(res)) max(p(res))]);

set(gca,'Ydir','reverse','Xdir','reverse'); % Adjust appearance of axis
axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
set(SUB1,'TickLabelInterpreter','Latex','Box','on','TickDir','out',...
    'XTick','','YTick','');
SUB1_P = SUB1.Position; % Save position of first sub-plot

%%% 4.2. Bottom plot (Column-integrated growth rate)
SUB2 = subplot(3,4,[10 11]);
[LRM_j,~] = LRM_growthrate(LRM',p',res); % Computes col-int growth rate
% Plot in p-space and use a dotted line for the boundary layer
plot(p(res),LRM_j,'color',[0 0 0],'Linewidth',lw+0.5,'Linestyle','--'); 
hold on; plot(p(ibl:res(end)),LRM_j(ibl:res(end)),...
    'color','k','Linewidth',lw+0.5);
hold on; line([p(ibl) p(ibl)],[-100 100],'color','k','Linewidth',lw,...
    'Linestyle','--');

% Plot's appearance
grid on; xlim([min(p(res)) max(p(res))]); ylim([-7.5 7.5]);
xlabel('$p_{j}\ [hPa]$','Interpreter','Latex','Fontsize',fsz);
ylabel('$\widehat{M}_{j}\ [d^{-1}]$','Interpreter','Latex','Fontsize',fsz);
set(SUB2,'TickLabelInterpreter','Latex','Box','on','TickDir','out',...
    'Fontsize',fsz,'Xdir','reverse');
SUB2.Position(2) = SUB1_P(2)-SUB2.Position(4);
SUB2.Position(1) = SUB1_P(1); SUB2.Position(3) = SUB1_P(3);
SUB1.Position = SUB1_P;

%%% 4.3. Right-side plot (Normalized eigenvector)
SUB3 = subplot(3,4,[4 8]);

% Eigenvector modulus and maximal eigenvalue real part
[V,D]=eig(LRM);[lam,ilam]=sort(real(diag(D))); LAM = lam(end);
EIGV = ((real(V(:,ilam(end)))).^2+(imag(V(:,ilam(end)))).^2).^0.5;
norm = 1/max(EIGV);

plot(norm*EIGV,p(res),'Linewidth',lw+0.5,'color',[0 0 0],...
    'Linestyle','--'); hold on;
plot(norm*EIGV(ibl:res(end)),p(ibl:res(end)),'Linewidth',lw+0.5,...
    'color','k');
line([-100 100],[p(ibl) p(ibl)],'color','k','Linewidth',lw,...
    'Linestyle','--');

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

%%% 4.4 Final adjustments

% Move all the subplots down
Yshift = 0.1;
SUB3.Position(2) = SUB3.Position(2)-Yshift;
SUB2.Position(2) = SUB2.Position(2)-Yshift;
SUB1.Position(2) = SUB1.Position(2)-Yshift;

% Add leading eigenvalue real part
annotation(F,'textbox',[SUB3.Position(1)+0.1*SUB3.Position(3) ...
    SUB3.Position(2)+0.885*SUB3.Position(4) 0.11 0.05],'String',...
    ['$\lambda=$',num2str(LAM,'%03.2f'),...
    '$\mathrm{d^{-1}}$'],'LineStyle','none','Interpreter','latex',...
    'FitBoxToText','off','BackgroundColor',[1 1 1],'Fontsize',fsz+1);
annotation(F,'textbox',[SUB1.Position(1)+0.01*SUB1.Position(3) ...
    SUB1.Position(2)+0.9*SUB1.Position(4) 0.05 0.05],...
    'FontSize',fsz+1,'Interpreter','latex','String','(a)',...
    'Linestyle','none');
annotation(F,'textbox',[SUB1.Position(1)+1.5*(1-p(ibl)/max(p(res)))*SUB1.Position(3) ...
    SUB1.Position(2)+0.9*SUB1.Position(4) 0.05 0.05],...
    'FontSize',fsz+1,'Interpreter','latex','String','(b)',...
    'Linestyle','none');
annotation(F,'textbox',[SUB1.Position(1)+0.01*SUB1.Position(3) ...
    SUB1.Position(2)+0.5*(1-p(ibl)/max(p(res)))*SUB1.Position(4) 0.05 0.05],...
    'FontSize',fsz+1,'Interpreter','latex','String','(c)',...
    'Linestyle','none');
annotation(F,'textbox',[SUB1.Position(1)+6.8*(1-p(ibl)/max(p(res)))*SUB1.Position(3) ...
    SUB1.Position(2)+0.5*(1-p(ibl)/max(p(res)))*SUB1.Position(4) 0.05 0.05],...
    'FontSize',fsz+1,'Interpreter','latex','String','(d)',...
    'Linestyle','none');

% Add (a),(b),(c),(d) panel indicators

%%% 4.5 Save figure in PDF format
% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
gcfsavepdf([basedir 'PDF_DATA\RCI_Fig03.pdf']);

else disp('The data for this figure was provided by Zhiming Kuang and is not be provided here.');
    disp('The reference is "Kuang (2012) Weakly Forced Mock Walker Cells, JAS".');
    disp('The corresponding authors email is kuang@fas.harvard.edu');
end
