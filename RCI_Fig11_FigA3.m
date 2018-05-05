%% RCI_Fig11_FigA3.m
% tbeucler - 4/24/2018
% Computes the leading eigenvalue real parts from different base states
% The surface temperature goes from 280K to 305K and the saturation deficit
% Optional figure: Activate with option_RH==1 ************************
% brings the RH max to 1 (blue=moist) or the RH min to 0 (orange=dry)
% Appendix figure ****************************************************
% For the standard case, explores the sensitivity to convective parameters

close all; clearvars; fclose('all');

%% 0. Parameters

% 0.1 Figure characteristics
AMP = 20; % Max of leading eigenvalue real part [1/day]
fsz = 12; % Fontsize
mz = 50; % Marker size
lw = 1.5; % Linewidth
M0 = 0.1; % Normalization factor to convert LRM in log-space [1/day]
L0 = 0.033; % Normalization factor to convert lambda in log-space [1/day]
N_cvc = 100; % Number of points in the convective parameter space

% 0.2 Color map from dry to moist
cmap = [1 0.5 0; 0.75 0.5 0.25; 0 0 0; 0.25 0.25 1; 0 0 0.75];

% 0.3 Toy model parameters
tau_BM = 3/24; % Parameter: Betts-Miller relaxation time [day]

% 0.4 Load data
load('MAT_DATA/RCI_SAM_RRTM_mult.mat'); % Radiative response
load('MAT_DATA/RCI_SAM_convection_mult.mat'); % Convective response

% 0.5 Options for figure 11
option_RH = 0; % [1] to display the dependence of lambda on sat. deficit

%% 1. Loop over different base states to get the leading real eigenvalues

% 1.1 Initialization
laRRTM = zeros(6,5); laBM0 = laRRTM; laBP0 = laRRTM;
laBM0pRRTM = laRRTM; laBP0pRRTM = laRRTM;

% 1.2 Loop over base states
SST_array = [280 285 290 295 300 305]; % Surface temperature domain
for iSST = 1:6, SST = SST_array(iSST);
    for imult = 1:5, disp(['SST=',num2str(SST),'K and RH multiplier = ',...
            num2str(imult),'/5']);
        
        % Restriction domain for the LRM
        res = ibl_mult(iSST,imult):itp_mult(iSST,imult);
        
        % Radiative LRM
        RADLRM = 24*3600*(LRM_LW_mult(iSST,imult,res,res)+...
            LRM_SW_mult(iSST,imult,res,res));
        % If imult==1, there is a singularity where q=0, so set to 0
        RADLRM(RADLRM==Inf)=0; RADLRM(RADLRM==-Inf)=0;
        RADLRM(isnan(RADLRM))=0;
        
        % Eigenvalue decomposition
        [l,~] = RCI_lead_eig_real_part(squeeze(RADLRM));
        laRRTM(iSST,imult) = l; % RRTM
        [l,~] = RCI_lead_eig_real_part(squeeze(LRM_BM_mult(...
            iSST,imult,res,res))); laBM0(iSST,imult) = l; % Betts-Miller
        [l,~] = RCI_lead_eig_real_part(squeeze(LRM_BP_mult(iSST,imult,...
            res,res))); laBP0(iSST,imult) = l; % Bulk-plume
        [l,~] = RCI_lead_eig_real_part(squeeze(LRM_BM_mult(iSST,imult,...
            res,res)+RADLRM)); laBM0pRRTM(iSST,imult) = l; % BM+RRTM
        [l,~] = RCI_lead_eig_real_part(squeeze(LRM_BP_mult(iSST,imult,...
            res,res)+RADLRM)); laBP0pRRTM(iSST,imult) = l; % BP+RRTM
        
        % For the reference state, save the corresponding eigenvector
        if SST==300 && imult==3, % Reference case
            [lamRRTM,eigRRTM] = RCI_lead_eig_real_part(squeeze(RADLRM));
            [lamBM0,eigBM0] = RCI_lead_eig_real_part(squeeze(LRM_BM_mult(...
                iSST,imult,res,res))); % Betts-Miller
            [lamBP0,eigBP0] = RCI_lead_eig_real_part(squeeze(LRM_BP_mult(...
                iSST,imult,res,res))); % Bulk-plume
            [lamBM0pRRTM,eigBM0pRRTM] = RCI_lead_eig_real_part(squeeze(...
                LRM_BM_mult(iSST,imult,res,res)+RADLRM)); % BM+RRTM
            [lamBP0pRRTM,eigBP0pRRTM] = RCI_lead_eig_real_part(squeeze(...
                LRM_BP_mult(iSST,imult,res,res)+RADLRM)); % BP+RRTM
            res300 = res; % Save LRM index space
            RADLRM300 = RADLRM; % Save radiative response
        end
        
    end
end

%% 2. Paper figure 11
F = figure('Position',[50 50 1000 500]);

%%% 2.1 Left panel: Leading eigenvalues for different base states
SUB1 = subplot(1,3,[1 2]);

% Scatter plot of leading real eigenvalues
if option_RH == 1,
    for imult = [1 3 5]; % Dry/Original/Moist
        % RRTM
        S = scatter(SST_array,sign(laRRTM(:,imult)).*log10(1+abs(laRRTM(:,imult)/L0)));
        S.MarkerEdgeColor = cmap(imult,:);
        S.SizeData = mz; S.LineWidth = lw; S.Marker = 'square'; hold on;
        % Betts-Miller
        S = scatter(SST_array,sign(laBM0(:,imult)).*log10(1+abs(laBM0(:,imult)/L0)));
        S.MarkerEdgeColor = cmap(imult,:);
        S.SizeData = mz; S.LineWidth = lw; hold on; grid on;
        % Bulk-plume
        S = scatter(SST_array,sign(laBP0(:,imult)).*log10(1+abs(laBP0(:,imult)/L0)));
        S.MarkerEdgeColor = cmap(imult,:);
        S.SizeData = mz; S.LineWidth = lw; S.Marker = '^'; hold on;
        % Betts-Miller + RRTM
        S = scatter(SST_array,sign(laBM0pRRTM(:,imult)).*log10(1+abs(laBM0pRRTM(:,imult)/L0)));
        S.MarkerEdgeColor = cmap(imult,:); S.MarkerFaceColor = cmap(imult,:);
        S.SizeData = mz; S.LineWidth = lw; hold on; grid on;
        % Bulk-plume + RRTM
        S = scatter(SST_array,sign(laBP0pRRTM(:,imult)).*log10(1+abs(laBP0pRRTM(:,imult)/L0)));
        S.MarkerEdgeColor = cmap(imult,:); S.MarkerFaceColor = cmap(imult,:);
        S.SizeData = mz; S.LineWidth = lw; S.Marker = '^'; hold on;
        % Legend
        LEG = legend('$\mathrm{RRTM}$','$\mathrm{BM}$','$\mathrm{BP}$',...
            '$\mathrm{RRTM+BM}$','$\mathrm{RRTM+BP}$');
    end
elseif option_RH == 0, imult = 3; % Original SAM RCE profile
    % RRTM
    S = scatter(SST_array,sign(laRRTM(:,imult)).*log10(1+abs(laRRTM(:,imult)/L0)));
    S.MarkerEdgeColor = [0 0.25 0.75]; S.SizeData = mz; S.LineWidth = lw; hold on;
    % Betts-Miller
    S = scatter(SST_array,sign(laBM0(:,imult)).*log10(1+abs(laBM0(:,imult)/L0)));
    S.MarkerEdgeColor = [0 0 0]; S.SizeData = mz; S.LineWidth = lw; hold on; grid on;
    % Bulk-plume
    S = scatter(SST_array,sign(laBP0(:,imult)).*log10(1+abs(laBP0(:,imult)/L0)));
    S.MarkerEdgeColor = [1 0 0]; S.SizeData = mz; S.LineWidth = lw; hold on;
    % Betts-Miller + RRTM
    S = scatter(SST_array,sign(laBM0pRRTM(:,imult)).*log10(1+abs(laBM0pRRTM(:,imult)/L0)));
    S.MarkerEdgeColor = [0 0 0]; S.MarkerFaceColor = [0 0 0];
    S.SizeData = mz; S.LineWidth = lw; hold on; grid on;
    % Bulk-plume + RRTM
    S = scatter(SST_array,sign(laBP0pRRTM(:,imult)).*log10(1+abs(laBP0pRRTM(:,imult)/L0)));
    S.MarkerEdgeColor = [1 0 0]; S.MarkerFaceColor = [1 0 0];
    S.SizeData = mz; S.LineWidth = lw; hold on;
    % Legend
    LEG = legend('$\mathrm{RRTM}$','$\mathrm{BM}$','$\mathrm{BP}$',...
        '$\mathrm{RRTM+BM}$','$\mathrm{RRTM+BP}$');
end
line([min(SST_array)-2 max(SST_array)+2],[0 0],'color',[0.6 0.6 0.6],...
    'Linewidth',1.5*lw);

% First panel's axis and appearance
xlim([min(SST_array)-2 max(SST_array)+2]);
if option_RH == 0, ylim([-log10(1+abs(3.5/L0)) log10(1+abs(1.5/L0))]);
else ylim([-log10(1+abs(20/L0)) log10(1+abs(3/L0))]);
end
G = gca; Y = str2double(G.YTickLabels);
set(gca,'YTickLabels',num2str(L0*sign(Y).*(10.^abs(Y)-1),'%03.2f'),...
    'Fontsize',fsz+2,'TickLabelInterpreter','Latex');
xlabel('$\mathrm{Surface\ temperature\ \left[K\right]}$','Fontsize',fsz+2,...
    'Interpreter','Latex');
ylabel('$\mathrm{Leading\ eigenvalue\ real\ part\ }\left[\mathrm{day^{-1}}\right]$',...
    'Fontsize',fsz+2,'Interpreter','Latex');
set(LEG,'Interpreter','Latex','Fontsize',fsz,'Location','northoutside',...
    'Orientation','horizontal');

%%% 2.2 Normalized eigenvector for standard convective parameters
SUB2 = subplot(1,3,3);

% Normalize the eigenvectors so that they all have the same column-integral
% and the maximum of all eigenvectors is 1
load('MAT_DATA/RCI300K_Modified_base_state.mat','p');
p = linspace(p(1),p(end),length(p))'; % Interpolated pressure [hPa]
eigBP0_int = -trapz(p(res300),eigBP0)/max(eigBP0);
NBP0 = 1/max(eigBP0);
NRRTM = -eigBP0_int/trapz(p(res300),eigRRTM);
NBM0 = -eigBP0_int/trapz(p(res300),eigBM0);
NBM0pRRTM = -eigBP0_int/trapz(p(res300),eigBM0pRRTM);
NBP0pRRTM = -eigBP0_int/trapz(p(res300),eigBP0pRRTM);

% Plot eigenvectors
plot(NRRTM*eigRRTM,p(res300),'Linewidth',lw,'color',[0 0.25 0.75],'Linestyle','--');
hold on; plot(NBM0*eigBM0,p(res300),'Linewidth',lw,'color',[0 0 0],'Linestyle','--');
hold on; plot(NBP0*eigBP0,p(res300),'Linewidth',lw,'color',[1 0 0],'Linestyle','--');
hold on; plot(NBM0pRRTM*eigBM0pRRTM,p(res300),'Linewidth',lw,'color',[0 0 0]);
hold on; plot(NBP0pRRTM*eigBP0pRRTM,p(res300),'Linewidth',lw,'color',[1 0 0]);

% Legend
L = legend('$\mathrm{RRTM}$','$\mathrm{BM}$','$\mathrm{BP}$',...
    '$\mathrm{RRTM+BM}$','$\mathrm{RRTM+BP}$');
set(L,'Interpreter','Latex','Fontsize',fsz);
L.Position(2) = L.Position(2)-0.3675;

% Axis limit
ylim([min(p(res300)) max(p(res300))]);xlim([-0.05 1]); SUB2.YAxisLocation = 'right';
set(gca,'Ydir','reverse'); grid on;
set(gca,'TickLabelInterpreter','Latex','Box','on','TickDir','out','Fontsize',fsz);
xlabel('Normalized eigenvectors','Interpreter','Latex','Fontsize',fsz+2);
ylabel('$p\ [hPa]$','Interpreter','Latex','Fontsize',fsz+2);
axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');

% Panel indicators
annotation(F,'textbox',[SUB1.Position(1)-0.1*SUB1.Position(3) ...
    SUB2.Position(2)+0.97*SUB2.Position(4) 0.04 0.07],'String',...
    'a)','LineStyle','none','Interpreter','latex','FitBoxToText','off',...
    'Backgroundcolor',[1 1 1],'Fontsize',fsz+1,'color','k');
annotation(F,'textbox',[SUB2.Position(1)-0.2*SUB2.Position(3) ...
    SUB2.Position(2)+0.97*SUB2.Position(4) 0.04 0.07],'String',...
    'b)','LineStyle','none','Interpreter','latex','FitBoxToText','off',...
    'Backgroundcolor',[1 1 1],'Fontsize',fsz+1,'color','k');

%%% 2.3 Save plot
% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
if option_RH==0, gcfsavepdf([basedir 'PDF_DATA\RCI_Fig11.pdf']); end

%% 3. Eigenvalues associated to a wider range of convective parameters
iSST = 5; imult = 3; % Reference values for standard 300K base state
%%% 3.1 Convective parameter spaces
logtau = linspace(-3,1,N_cvc); tau_space = 10.^logtau; % tau_BM space [day]
load('MAT_DATA/RCI_SAM_RRTM_mult.mat');
Qcool = QLW_mult(iSST,imult)-QSW_mult(iSST,imult); % Reference value of radiative cooling
Qcool_space = Qcool*tau_BM./tau_space; % Qcool space [W/m2]
%%% 3.2 Loop over convective parameters and get eigenvalues
lamBM = zeros(N_cvc,1); lamBP = lamBM; lamBMpRRTM = lamBM; lamBPpRRTM = lamBM;
for icvc = 1:N_cvc
    % Betts-Miller alone
    [l,~] = RCI_lead_eig_real_part(squeeze(LRM_BM_mult(iSST,imult,...
        res300,res300)*tau_BM/tau_space(icvc))); lamBM(icvc) = l;
    % Bulk-plume alone
    [l,~] = RCI_lead_eig_real_part(squeeze(LRM_BP_mult(iSST,imult,...
        res300,res300)*Qcool_space(icvc)/Qcool)); lamBP(icvc) = l;
    % Betts-Miller + RRTM
    [l,~] = RCI_lead_eig_real_part(squeeze(LRM_BM_mult(iSST,imult,...
        res300,res300)*Qcool_space(icvc)/Qcool)+squeeze(RADLRM300));
    lamBMpRRTM(icvc) = l;
    % Bulk-plume + RRTM
    [l,~] = RCI_lead_eig_real_part(squeeze(LRM_BP_mult(iSST,imult,...
        res300,res300)*Qcool_space(icvc)/Qcool)+squeeze(RADLRM300));
    lamBPpRRTM(icvc) = l;
end

%% 4. Paper Figure A3
F = figure('position',[100 100 500 400]);

% 4.1 Plot
plot(logtau,0*logtau,'Linewidth',lw,'color',[0.6 0.6 0.6]); hold on;
line([log10(tau_BM) log10(tau_BM)],[-99 99],'Linewidth',lw,...
    'color',[0.6 0.6 0.6]); hold on;
line([1+log10(tau_BM) 1+log10(tau_BM)],[-99 99],'Linewidth',lw,...
    'color',[0.6 0.6 0.6]); hold on;
line([log10(tau_BM)-1 log10(tau_BM)-1],[-99 99],'Linewidth',lw,...
    'color',[0.6 0.6 0.6]); hold on;
BPR = plot(logtau,sign(lamBPpRRTM).*log10(1+abs(lamBPpRRTM)/M0),...
    'Linewidth',lw,'color',[1 0 0]); hold on;
BMR = plot(logtau,sign(lamBMpRRTM).*log10(1+abs(lamBMpRRTM)/M0),...
    'Linewidth',lw,'color',[0 0 0]); hold on;
BP = plot(logtau,sign(lamBP).*log10(1+abs(lamBP)/M0),...
    'Linewidth',lw,'color',[1 0 0],'Linestyle','--'); hold on;
BM = plot(logtau,sign(lamBM).*log10(1+abs(lamBM)/M0),...
    'Linewidth',lw,'color',[0 0 0],'Linestyle','--'); hold on;
RRTM = plot(logtau,sign(lamRRTM)*log10(1+abs(lamRRTM)/M0)*logtau.^0,...
    'Linewidth',lw,'Linestyle','--','color',[0 0.25 0.75]);

% 4.2 Axis limits and labels
YMIN = -2.56; YMAX = 1; ylim([YMIN YMAX]); xlim([min(logtau) max(logtau)]);
G = gca; L = str2double(G.YTickLabels); % Colorbar's ticks in log-space
G.YTickLabels=num2str(M0*sign(L).*(10.^abs(L)-1),'%03.2f'); % Convert them to lin-space
xlabel('$\mathrm{log_{10}}\left(\tau_{\mathrm{BM}}\ \left[\mathrm{day}\right]\right)$',...
    'Fontsize',fsz+2,'Interpreter','Latex');
title('$\widehat{Q}_{\mathrm{BP}} \left[\mathrm{W.m^{-2}}\right]$',...
    'Fontsize',fsz+2,'Interpreter','Latex','color',[1 0 0]);
ylabel('$\mathrm{Leading\ eigenvalue\ real\ part\ }\left[\mathrm{day^{-1}}\right]$',...
    'Fontsize',fsz+2,'Interpreter','Latex');
set(gca,'TickLabelInterpreter','Latex','Box','on','TickDir','out',...
    'Fontsize',fsz); grid on; hold on;

% 4.3 Top-axis labels
Q_labels = {'10000';'1500';'150';'15';'3'};
logtau_Qlab = interp1(Qcool_space,logtau,[10000 10*Qcool Qcool 0.1*Qcool 3]);
for i=1:length(logtau_Qlab)
    T=text(logtau_Qlab(i),YMAX-0.05*(YMAX-YMIN),char(Q_labels(i)),'Color',...
        [1 0 0],'HorizontalAlignment','center','VerticalAlignment',...
        'top','Interpreter','Latex','Fontsize',fsz,'Backgroundcolor',[1 1 1]);
    line([logtau_Qlab(i) logtau_Qlab(i)],[YMAX-0.05*(YMAX-YMIN) YMAX],...
        'Color',[1 0 0],'LineWidth',1);
end

% 4.4 Legend
L = legend([RRTM BM BP BMR BPR],{'$\mathrm{RRTM}$','$\mathrm{BM}$',...
    '$\mathrm{BP}$','$\mathrm{RRTM+BM}$','$\mathrm{RRTM+BP}$'});
set(L,'Interpreter','Latex','Fontsize',fsz,'Location','SouthEast');

% 4.5 Save plot
% Trick to get current directory on different machines
thisfile  = which(mfilename);
basedir = thisfile(1:strfind(thisfile,mfilename)-1);
% Save plot
gcfsavepdf([basedir 'PDF_DATA\RCI_FigA3.pdf']);
