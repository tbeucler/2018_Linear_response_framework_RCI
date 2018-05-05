function QV = RH_inv(RH,T,p)

% Inverses RH to compute specific hum. taking into account liquid and ice
% Computes QV [1] from RH [1], T[K] and p[hPa]
% Uses Bolton formula for liq
% Uses Smithsonian tables (Kerry, 4.4.15)

Rd = 287.04;
Rv = 461.50;
eps = Rd/Rv; % Ratio of the molecular weight of water vapor to the molecular weight of dry air
Tcrit = -30; % Critical temperature to separate between liquid and ice

x = T-273.15; % Temperature in Celcius
e_sat_l = 6.112*exp(17.67*x/(x+243.5));
e_sat_s = exp(23.33086-6111.72784/T+0.15215*log(T));
e_sat = e_sat_l*(x>0)+e_sat_s*(x<Tcrit)+...
    ((Tcrit-x)/Tcrit*e_sat_l+x/Tcrit*e_sat_s)*(x<=0&&x>=Tcrit);
e = RH*e_sat;
QV = eps*e./(p+(eps-1)*e); QV = max(QV,0);


