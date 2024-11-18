function [Q_oven Q_BPR Q_Pump] = calc_conditions_Qbolha(P2,T,Q) 

%P2 = [Pexperimental Patm].
%T = [ T_bolha T_oven T_BPR T_Pump].
%Q is the volumetric flow rate measured in the soap bubble flow meter in mL/min.

%[rho_pump rho_oven rho_BPR M Q_pump Q_BPR] = calc_conditions(P,T,Q)

P=P2(1); %pressão da experiencia
Patm = P2(2); %pressão atmosferica

T_bolha = T(1);
T_oven = T(2);
T_BPR = T(3);
T_Pump = T(4);

rho_bolha = densityCO2(Patm/1.01325, T_bolha+273.15); % (g/cm3)
rho_oven = densityCO2(P/1.01325, T_oven+273.15); % (g/cm3)
rho_BPR = densityCO2(P/1.01325, T_BPR+273.15); % (g/cm3)
rho_pump = densityCO2(P/1.01325, T_Pump+273.15); % (g/cm3)
M = Q*rho_bolha; %g/mL

Q_oven = M/rho_oven; %mL/min
Q_BPR = M/rho_BPR; %mL/min
Q_Pump = M/rho_pump; %mL/min

if Q_BPR <= 0.1
    disp('DANGER! BPR flow rate should be higher than 0.1 mL/min.')
end

%restrições
L = 1118.2; %cm
id = 0.053; %cm
e = 0.0001; %cm
Ro = id/2-e; %cm
Dc = 30; %cm

% D12 = 10^-4; %cm2/s
% 
A = pi*Ro^2; %cm2
% lamda = (Dc/2)/Ro;
% 
% visco = viscosityCO2(T_oven+273.15, rho_oven); % (g/cm.s)
% Sc = visco/rho_oven/D12; 
%  
%         
    Uo = Q_oven/60/A %cm/s
%     t = L*A/Q_oven %min
%     
%     Pex = Uo*L/D12 
%     disp('large?')%must be large (>10^4)
%     Re = rho_oven*Uo*(2*Ro)/visco
%     disp('laminar?')%must be laminar (<2100)
%     
%     De = Re/sqrt(lamda);
%     D = D12+Ro^2*Uo^2/(48*D12);
%     
%     D/(Uo*L)
%     disp('<0.01?')
%     De*sqrt(Sc)
%     disp('<10?')
%     Uo*L/D
%     disp('>1000?')
    

end