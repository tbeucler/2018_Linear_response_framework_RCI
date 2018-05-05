function [ Qrad,FupLW,FdoLW,tau ] = CRSRCE_LWgray( p,q,T,Ts,kappa )
% LW gray modes - 11/10/2017
%%% INPUT
% p: Pressure profile [hPa]
% q: Specific humidity profile [kg/kg]
% T: Temperature profile [K]
% Ts: Surface temperature [K]
% kappa: Absorption coefficient [m^2/kg]
%%% OUTPUT
% Qrad: Radiative heating [K/d]
% Fup: Upwelling flux [W/m2]
% Fdo: Downwelling flux [W/m2]
% tau: Optical thickness profile [1]

close all;fclose('all');

cp = 1005.7; % Specific heat capacity of dry air at constant pressure [J/K/kg]
g = 9.81; % Gravity constant [m/s2]
Nbl = 5e3; % Boundary layer resolution [def=5e3]
NI = 1e3; % Minimal resolution (max linear resolution is 10 times that) [def=1e3]
sig = 5.67e-8; % Stefan-Boltzmann constant (W/m^2/K^4)
tauc = 1e3; % Surface optical thickness above which the Laplace approx. is made

Np = numel(p); ps = p(1); % Profile resolution and surface pressure [hPa]
Kpd = 24*3600*g*sig*Ts^4/(cp*100*ps); % Conversion factor from [1] to [K/d]

tau = cumtrapz(100*p,kappa*(p/ps).*q/g); tau = tau-tau(end); % Optical thickness
taus = tau(1); % Surface optical thickness
i0 = find(diff(tau)==0,1); % First index for which opt. thickness>0
fdo = exp(tau).*(T/Ts).^4; fup = exp(-tau).*(T/Ts).^4; % Integrands in the Schwarzschild model

if taus<tauc, % Trapezoidal integration
    
    Nres = min(10*NI,max(floor((kappa^1.5)*10^5.8468),NI)); % Necessary resolution (11/7/2017)
    tauI = linspace(tau(1),tau(i0),Nres); % High res optical thickness
    TI = interp1(tau(1:i0),T(1:i0),tauI,'linear','extrap'); % Interpolated temperature
    qI = interp1(tau(1:i0),q(1:i0),tauI,'linear','extrap'); % Interpolated specific humidity
    pI = interp1(tau(1:i0),p(1:i0),tauI,'linear','extrap'); % Interpolated pressure
    fupI = interp1(tau(1:i0),fup(1:i0),tauI,'linear','extrap'); % Interpolated upwelling integrand
    fdoI = interp1(tau(1:i0),fdo(1:i0),tauI,'linear','extrap'); % Interpolated downwelling integrand
    
    Fup = zeros(1,Nres); Fdo = Fup; Fup(1,1) = 1; Fup(1,2) = 1; % Upwelling & downwelling fluxes
    for ip = 1:Nres, u = 1:(ip-1); d = (ip+1):Nres;
        if ip>02, Fup(1,ip) = exp(tauI(ip))*(exp(-tauI(1))-trapz(tauI(u),fupI(u))); end
        if ip<Nres-1, Fdo(1,ip) = -exp(-tauI(ip))*trapz(tauI(d),fdoI(d)); end
    end
    
    Q = Kpd*100*ps*kappa*(pI/ps).*qI/g.*(Fup+Fdo-2*(TI/Ts).^4); % Radiative heating [K/d]
    Qrad = interp1(tauI,Q,tau(1:i0),'linear','extrap'); % Interpolated back to original levels
    Qrad(i0+1:Np) = 0; % Zeroed where tau=0
    FupLW = interp1(tauI,Fup,tau(1:i0),'linear','extrap'); % Upwelling flux [1]
    FupLW(i0+1:Np) = FupLW(i0); FupLW = sig*Ts^4*FupLW; % Zero gradient where tau=0 and converted to [W/m2]
    FdoLW = interp1(tauI,Fdo,tau(1:i0),'linear','extrap'); % Downwelling flux [1]
    FdoLW(i0+1:Np) = 0; FdoLW = sig*Ts^4*FdoLW; % Zeroed where tau=0 and converted to [W/m2]
    
else % Asymptotic expression for large tau with boundary layer correction for tau in [0,tauc/2]
    
    Qfac = Kpd*100*ps*kappa*(p/ps).*q/g; % Pre-factor for radiative heating
    
    Nres = NI; % Sets resolution to minimal resolution
    tauI = linspace(tau(1),tau(i0),Nres); % High res opt. thickness
    TI = interp1(tau(1:i0),T(1:i0).^4/Ts^4,tauI,'linear','extrap'); % Interpolated temperature
    QfacI = interp1(tau(1:i0),Qfac(1:i0),tauI,'linear','extrap'); % Interpolated rad heating pre-factor
    
    Fup = exp(tauI-taus)+TI.*(1-exp(tauI-taus)); % Asymptotic form of the upwelling flux
    Fdo = TI.*(1-exp(-tauI)); % Asymptotic form of the downwelling flux
    QI = QfacI.*(Fup+Fdo-2*TI); % Asymptotic form of the radiative heating
    
    %%% Boundary layer at the top (Small tau) to resolve Fdown %%%
    
    [~,itop] = min(abs(tau-(tauc/2))); if itop==i0, itop=itop-1; end % Index for which tau=tauc/2
    tauT = linspace((tauc/2),0,Nbl); % Optical thickness coordinate in the BL
    fdoT = interp1(tau(itop:i0),fdo(itop:i0),tauT,'linear','extrap'); % Extrapolated downwelling integrand
    FdoT = zeros(1,Nbl); % Initialization of the downwelling BL flux computed below ?
    for ibl = 1:(Nbl-1), d = ibl:Nbl; FdoT(ibl) = -exp(-tauT(ibl))*trapz(tauT(d),fdoT(d)); end
    
    TT = interp1(tau(itop:i0),T(itop:i0).^4/Ts^4,tauT,'linear','extrap'); % Extrapolated temperature in BL
    QfacT = interp1(tau(itop:i0),Qfac(itop:i0),tauT,'linear','extrap'); % Extrapolated rad heating pre-factor
    
    [~,itopI] = min(abs(tauI-(tauc/2))); if itopI==Nres, itopI=itopI-1; end % Intex for which tauI=tauc/2
    FupT = interp1(tauI(itopI:Nres),Fup(itopI:Nres),tauT,'linear','extrap'); % Interpolated upwelling flux in BL
    
    QT = QfacT.*(FupT+FdoT-2*TT); % Radiative heating in BL
    
    %%% Interpolate back to original coordinates %%%
    
    Qrad = interp1(tauI,QI,tau(1:i0),'linear','extrap'); % Interpolates asymptotic solution
    Qrad(itop:i0) = interp1(tauT,QT,tau(itop:i0),'linear','extrap'); % Interpolates BL solution
    ineg = find(Qrad<0,1); Qrad(1:ineg-1)=0; % Sets radiative heating to 0 for high tau
    Qrad(i0+1:Np) = 0; % Sets radiative heating to 0 where tau=0
    
    FupLW = interp1(tauI,Fup,tau(1:i0),'linear','extrap'); % Interpolates asymptotic solution
    FupLW(itop:i0) = interp1(tauT,FupT,tau(itop:i0),'linear','extrap'); % Interpolates BL solution
    FupLW(i0+1:Np) = FupLW(i0); FupLW(1:ineg-1)=1; % Sets upwelling flux to FupLW(i0)
    % where tau=0 and 1 for high tau
    
    FdoLW = interp1(tauI,Fdo,tau(1:i0),'linear','extrap'); % Interpolates asymptotic solution
    FdoLW(itop:i0) = interp1(tauT,FdoT,tau(itop:i0),'linear','extrap'); % Interpolates BL solution
    FdoLW(i0+1:Np) = 0; FdoLW(1:ineg-1)=1; % Sets downwelling flux to 0 where
    % tau=0 and 1 for high tau
    
    FupLW = sig*Ts^4*FupLW; FdoLW = sig*Ts^4*FdoLW; % Conversion to [W/m2]
    
end

end

