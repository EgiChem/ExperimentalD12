function Qmax = calc_Qmax(P,T) %mL/min


P = 175; % bar
T = 40; % ºC

% P(bar), T(ºC)
rho = densityCO2(P/1.01325, T+273.15); % (g/cm3)
visco = viscosityCO2(T+273.15, rho); % (g/cm.s) 

L = 1688.7; %cm
id = 0.0516; %cm
e = 0.0001; %cm
Ro = id/2; %cm
Dc = 30; %cm

D12 = 10^-4; %cm2/s

A = pi*Ro^2; %cm2
lamda = (Dc/2)/Ro;

Sc = visco/rho/D12; 
 
[Q, fval] = fminsearch(@(Q)fobj(Q,A,L,D12,rho,visco,lamda,Ro,Sc),0.05)


    Uo = Q/60/A; %cm/s
    t = L*A/Q; %min
    
    Pex = Uo*L/D12; %must be large (>10^4)
    Re = rho*Uo*(2*Ro)/visco; %must be laminar (<2100)
    
    De = Re/sqrt(lamda);
    D = D12+Ro^2*Uo^2/(48*D12);

    D/(Uo*L)
    De*sqrt(Sc) 
    Uo*L/D
    
    De*sqrt(Sc)
    disp('<10?')


end


function F = fobj(Q,A,L,D12,rho,visco,lamda,Ro,Sc)
        
    Uo = Q/60/A; %cm/s
    t = L*A/Q; %min
    
    Pex = Uo*L/D12; %must be large (>10^4)
    Re = rho*Uo*(2*Ro)/visco; %must be laminar (<2100)
    
    De = Re/sqrt(lamda);
    D = D12+Ro^2*Uo^2/(48*D12);
    
    if D/(Uo*L)>0.01 || De*sqrt(Sc)>10 || Uo*L/D<1000
       Q = 10^-3; 
    end
    Q;
    F = 1/Q;
    
end

