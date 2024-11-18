function retrictions(P,T,Q)

% Q is the flowrate inside the column

rho = densityCO2(P/1.01325, T+273.15); % (g/cm3)
visco = viscosityCO2(T+273.15, rho); % (g/cm.s)

L = 1688.7; %cm
% d = 0.053;; %cm
% e = 0.0001; %cm
% Ro = id/2-e; %cm
id = 0.0516;; %cm
Ro = id/2; %cm
Dc = 30; %cm

D12 = 10^-4; %cm2/s

A = pi*Ro^2; %cm2
lamda = (Dc/2)/Ro;

Sc = visco/rho/D12; 
 
        
    Uo = Q/60/A %cm/s
    t = L*A/Q %min
    
    Pex = Uo*L/D12 
    disp('large?')%must be large (>10^4)
    Re = rho*Uo*(2*Ro)/visco
    disp('laminar?')%must be laminar (<2100)
    
    De = Re/sqrt(lamda);
    D = D12+Ro^2*Uo^2/(48*D12);
    
    D/(Uo*L)
    disp('<0.01?')
    De*sqrt(Sc)
    disp('<10?')
    Uo*L/D
    disp('>1000?')
    
%     F = De*sqrt(Sc);
end
    