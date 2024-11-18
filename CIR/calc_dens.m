function rho = calc_dens(P,T) 
%calc_dens determines carbon dioxide density for a certain pressure and temperature:
%rho = calc_dens(P,T),
%with P(bar) and T(ºC)
rho = densityCO2(P/1.01325, T+273.15); % (g/cm3)

end