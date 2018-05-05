function q=q_sat(p,T)

% tbeucler - 2017
% Computes the saturation mixing ratio [1] as a function of pressure [hPa]
% and temperature [K]
% Uses Bolton formula and Dalton's law

Rd = 287.04; % Specific gas constant for dry air [J/kg/K]
Rv = 461.50; % Specific gas constant for water vapor [J/kg/K]
eps = Rd/Rv; % Ratio of the molecular weight of water vapor to the molecular weight of dry air
Tcrit = -30; % Critical temperature to separate between liquid and ice [C]

x = T-273.15; % Temperature in Celcius
e_sat_l = 6.112*exp(17.67*x/(x+243.5));
e_sat_s = exp(23.33086-6111.72784/T+0.15215*log(T));
e_sat = e_sat_l*(x>0)+e_sat_s*(x<Tcrit)+...
    ((Tcrit-x)/Tcrit*e_sat_l+x/Tcrit*e_sat_s)*(x<=0&&x>=Tcrit);
q=eps.*e_sat./(p-(1-eps)*e_sat);