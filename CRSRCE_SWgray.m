function [ Qrad,FdoSW,tau ] = CRSRCE_SWgray( p,q,kappa,S0,pbd )
% SW gray modes - 11/10/2017
% One-stream model
%%% INPUT
% p: Pressure profile [hPa]
% q: Specific humidity profile [kg/kg]
% kappa: Absorption coefficient [m^2/kg]
% S0: Average solar constant [W/m^2]
% pbd: Exponent used for pressure broadening
%%% OUTPUT
% Qrad: Radiative heating [K/d]
% Fdo: Downwelling flux [W/m2]
% tau: Optical thickness profile [1]

close all;fclose('all');

cp = 1005.7; % Specific heat capacity of dry air at constant pressure [J/K/kg]
g = 9.81; % Gravity constant [m/s2]
Np = numel(p); ps = p(1); % Profile resolution and surface pressure [hPa]
Kpd = 24*3600*g/(cp*100*ps); % Conversion factor from [W/m2] to [K/d]

tau = cumtrapz(100*p,kappa*(p/ps).^pbd.*q/g); tau = tau-tau(end); % Optical thickness
i0 = find(diff(tau)==0,1); % First index for which opt. thickness>0

Nres = 1e3; % Resolution for the high res optical thickness profile
tauI = linspace(tau(1),tau(i0),Nres); % High res optical thickness
qI = interp1(tau(1:i0),q(1:i0),tauI,'linear','extrap'); % Interpolated specific humidity
pI = interp1(tau(1:i0),p(1:i0),tauI,'linear','extrap'); % Interpolated pressure

Fdo = S0*exp(-tauI); % Downwelling SW flux for the one-stream model
Q = Kpd*100*ps*kappa*(pI/ps).^pbd.*qI/g.*Fdo; % Radiative heating [K/d]

FdoSW = interp1(tauI,Fdo,tau(1:i0),'linear','extrap'); % Downwelling flux [1]
FdoSW(i0+1:Np) = S0; % Set to S0 where tau=0 and converted to [W/m2]

Qrad = interp1(tauI,Q,tau(1:i0),'linear','extrap'); % Interpolated back to original levels
Qrad(i0+1:Np) = 0; % Zeroed where tau=0

end