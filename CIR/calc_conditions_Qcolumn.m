function [rho_pump rho_oven rho_BPR M Q_pump Q_BPR] = calc_conditions_Qcolumn(P,T,Q) 
%calc_conditions determines carbon dioxide conditions for a certain pressure (bar) in different components of the equipment:
%- pump, according to the temperature in the bath (T_bath (ºC));
%- inside the column, according to the temperature in the oven (T_oven (ºC));
%- in the BPR, according to its temperature (T_BPR (ºC)).
%These temperatures should be given as T = [ T_bath T_oven T_BPR].
%Q is the desired volumetric flow rate inside the column (oven) in mL/min.
%The function calculates CO2 densities rho(g/cm3) in the pump, oven and BPR,
%the massic flow rate M(g/mL) and the correspondent volumetric flow rates in
%the pump and BPR:
%[rho_pump rho_oven rho_BPR M Q_pump Q_BPR] = calc_conditions(P,T,Q)


T_bath = T(1);
T_oven = T(2);
T_BPR = T(3);

rho_pump = densityCO2(P/1.01325, T_bath+273.15); % (g/cm3)
rho_oven = densityCO2(P/1.01325, T_oven+273.15); % (g/cm3)
rho_BPR = densityCO2(P/1.01325, T_BPR+273.15); % (g/cm3)

M = Q*rho_oven; %g/mL

Q_pump = M/rho_pump; %mL/min
Q_BPR = M/rho_BPR; %mL/min
if Q_BPR <= 0.1
    disp('DANGER! BPR flow rate should be higher than 0.1 mL/min.')
end

end