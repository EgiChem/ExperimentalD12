function [Q_oven Q_BPR Q_buble] = calc_conditions_Qpump(P,T,Q) 

%calc_conditions determines carbon dioxide conditions for a certain pressure (bar) in different components of the equipment:
%- pump, according to the temperature in the bath (T_bath (ºC));
%- inside the column, according to the temperature in the oven (T_oven (ºC));
%- in the BPR, according to its temperature (T_BPR (ºC)).
%These temperatures should be given as T = [ T_bath T_oven T_BPR].
%Q is the volumetric flow rate in the pump (oven) in mL/min.
%The function calculates CO2 densities rho(g/cm3) in the pump, oven and BPR,
%the massic flow rate M(g/mL) and the correspondent volumetric flow rates in
%the pump and BPR:
%[rho_pump rho_oven rho_BPR M Q_pump Q_BPR] = calc_conditions(P,T,Q)


T_bath = T(1);
T_oven = T(2);
T_BPR = T(3);
T_amb = T(4);

rho_pump = densityCO2(P(1)/1.01325, T_bath+273.15); % (g/cm3)
rho_oven = densityCO2(P(1)/1.01325, T_oven+273.15); % (g/cm3)
rho_BPR = densityCO2(P(1)/1.01325, T_BPR+273.15); % (g/cm3)
rho_buble = densityCO2(P(2)/1.01325, T_amb+273.15); % (g/cm3)

M = Q*rho_pump; %g/mL

Q_oven = M/rho_oven %mL/min
Q_BPR = M/rho_BPR; %mL/min
Q_buble = M/rho_buble;

if Q_BPR <= 0.1
    disp('DANGER! BPR flow rate should be higher than 0.1 mL/min.')
end

%restrições
L = 1118.2; %cm atualizado com ultimo corte a 23/04/2014
id = 0.0523; %cm 
% e = 0.0001; %cm
e = 0; %cm
Ro = id/2; %cm
Dc = 30; %cm

D12 = 1.2*10^-4; %cm2/s

A = pi*Ro^2; %cm2
lamda = (Dc/2)/Ro;

visco = viscosityCO2(T_oven+273.15, rho_oven); % (g/cm.s)
Sc = visco/rho_oven/D12; 
 
        
    Uo = Q_oven/60/A %cm/s
%     Uo = 1.06
%     t = L*A/Q_oven %min
    t = L/Uo %min

    Pex = Uo*L/D12 
    disp('large?')%must be large (>10^4)
    Re = rho_oven*Uo*(2*Ro)/visco
    disp('laminar?')%must be laminar (<2100)
    
    De = Re/sqrt(lamda);
    D = D12+Ro^2*Uo^2/(48*D12);
    
    D/(Uo*L)
    disp('<0.01?')
    De*sqrt(Sc)
    disp('<10?')
    Uo*L/D
    disp('>1000?')
    
    nada='';

end