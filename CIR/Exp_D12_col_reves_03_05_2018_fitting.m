function [Resultados_CIR_fitt D12Calc_variance time_max_absorbance] = Exp_D12_col_reves_03_05_2018_fitting(FileName,FileNumber, altura)

global R_col L_col Dados_Experimentais h_pico

%% Programa de c�lculo das difusividades a partir dos pontos experimentais (por fitting e por an�lise da largura do pico ou "m�todo dos momentos")

% NOTAS IMPORTANTES : 
%
%                     Escolha do m�todo de c�lculo das �reas (Metodo_Integral): 
%                       - 1 para o m�todo dos trap�zios,
%                       - 2 para o m�todo de Simpson 1/3,
%                       - 3 para o m�todo de Simpson 3/8.
%                ------------------------------------------------
% 
%                     Escolha da forma a partir do qual se obt�m a velocidade linear dentro da coluna de difus�o (Veloc_Calc):
%                       - 11 para a bomba;
%                       - 22 para o medidor de bolha de sab�o.

Metodo_Integral = 1; 
Veloc_Calc = 11;


%% Dados da coluna

d_enrolamento = 30.0; % Di�metro do enrolamento (cm)
diam_col = (0.526)/10; % Di�metro interno da coluna (cm) 
L_col = (16.887-0.028)*100; % Comprimento da coluna (cm)
espess_revestimento = 1E-4; % (cm) 
% R_col = (diam_col)/2 - espess_revestimento; %(cm) Raio da coluna
R_col = (diam_col)/2;%(cm) Raio da coluna


Dl = 0.1e-6; %difus�o no s�lido

% h_abertuta_pico = 0.607; %60.7 % da altura do pico
% fator_h = 1; % fator para 60.7 % da altura do pico

h_abertuta_pico = 0.5; %50.0 % da altura do pico
fator_h = 5.545/4 ; % fator para 50.0 % da altura do pico

%% Informa��o acerca do soluto (massa, ou volume, ou concentra��o do soluto injectado
vol_inject = 0.1*1E-3; % mL

%% Leitura dos dados de um ficheiro de dados em txt (tempo vs Absorv�ncia) 

[data_matrix, text_data, alldata] = xlsread(...
    FileName,FileNumber); 

Dados_Experimentais.wavelength = data_matrix(1,3); % (nm)
Dados_Experimentais.solvente = text_data{2,4};
Dados_Experimentais.solvente_formula = text_data{3,4};
Dados_Experimentais.soluto = text_data{2,5};
Dados_Experimentais.Pexp = data_matrix(1,6); %(bar)
Dados_Experimentais.Texp = data_matrix(1,7)+ 273.15; % (K)
Dados_Experimentais.Texp_bomba = data_matrix(1,8)+ 273.15; % (K)
Dados_Experimentais.Pexp_bomba = data_matrix(1,6); % (bar)

Dados_Experimentais.caudal_bomba = data_matrix(1,9); % (mL/min)
Dados_Experimentais.Pamb = data_matrix(1,10); %(bar)
Dados_Experimentais.Tbolha = data_matrix(1,11)+ 273.15; % (K)
Dados_Experimentais.caudal_bolha = data_matrix(1,12); % (mL/min)


Dados_Experimentais.tempo = data_matrix(1:end,1); % tempo (s)
Dados_Experimentais.absorbance = data_matrix(1:end,2); % Unidades de Absorv�ncia (AU)

h_pico=altura;

% % smoothing signal
% Dados_Experimentais.absorbance= sgolayfilt(Dados_Experimentais.absorbance,2,9);


%% Figura 1 - Absorv�ncia vs tempo
figure(1)
subplot(1,3,1)
plot(Dados_Experimentais.tempo, Dados_Experimentais.absorbance, '*k')
xlabel('{\itt} (s)')
ylabel('Absorbance (AU)')

%% ADICIONAR
% % % % % % Title([Dados_Experimentais.soluto, ' in ', Dados_Experimentais.solvente_formula, ' - ', num2str(Dados_Experimentais.Texp), ...
% % % % % %     ' K / ', num2str(Dados_Experimentais.Pexp/10), ' MPa (\lambda = ',...
% % % % % %     num2str(Dados_Experimentais.wavelength),' nm)'])
%% --------------------------------------------------------------

%% Tratamento dos dados experimentais

% Dete��o do tempo correpondente � absorv�ncia m�xima

max_absorbance = max(Dados_Experimentais.absorbance); 

for i = 1:length(Dados_Experimentais.tempo)
    if Dados_Experimentais.absorbance(i) == max_absorbance
        time_max_absorbance = Dados_Experimentais.tempo(i);
        posicao_max_absorbance = i;
    end 
end


%% Verifica�ao da simetria do pico obtido (atrav�s do par�metro S10)
 
% Determina�ao da exist�ncia de um ponto experimenal a 10 % da altura do pico (1� Parte)
% ---(indice correspondente � abs max)
n1 = 0;
for j1 = 1:posicao_max_absorbance
    if Dados_Experimentais.absorbance(j1) == h_pico*max_absorbance 
        S10_1_absorbance = Dados_Experimentais.absorbance(j1);
        S10_1_time = Dados_Experimentais.tempo(j1);
        n1 = n1 + 1;
    end 
end
% ---
npontos_int = 3; %pontos usados para cado lado na interpola��o

if n1 ==0 %se n�o existir nenhum ponto correspondente a 0.1*maxabs
    dif_absorb1 = abs(Dados_Experimentais.absorbance(1:posicao_max_absorbance)- h_pico*max_absorbance); %diferen�a entre cada ponto e o maxabs
    min_dif_absorb1 = min(dif_absorb1);
    
% ---(indice correspondente � diferen�a minima, ie, ponto + proximo de 0.1*maxabs)    
    for j1 = 1:posicao_max_absorbance
        if min_dif_absorb1 == dif_absorb1(j1)
            absorbance_1 = Dados_Experimentais.absorbance(j1);
            tempo_1 = Dados_Experimentais.tempo(j1);
            posicao_min_1 = j1;
        end
    end
% ---
    
    if absorbance_1 < h_pico*max_absorbance % se o ponto + proximo esta abaixo
    
            absorbance_interpolar_1 = Dados_Experimentais.absorbance(posicao_min_1-(npontos_int-1):posicao_min_1+npontos_int); %pontos a usar na fun��o interpola��o: 2 para tr�s e 3 para a frente
            tempo_intepolar_1 = Dados_Experimentais.tempo(posicao_min_1-(npontos_int-1):posicao_min_1+npontos_int); %igual para o y (tempo)
    
    elseif absorbance_1 > h_pico*max_absorbance % se o ponto + proximo esta acima
        
            absorbance_interpolar_1 = Dados_Experimentais.absorbance(posicao_min_1-npontos_int:posicao_min_1+(npontos_int-1)); %pontos a usar na fun�ao de interpola�ao: 3 para tr�s e 2 para a frente
            tempo_intepolar_1 = Dados_Experimentais.tempo(posicao_min_1-npontos_int:posicao_min_1+(npontos_int-1));%igual para o y (tempo)
       
    end
    
    p_1 = polyfit(absorbance_interpolar_1,tempo_intepolar_1,2); %fun��o de interpola��o: polin�mio de 2� grau
    S10_1_time = polyval(p_1,h_pico*max_absorbance); %tempo correspondente a 0.1*maxabs
    
    tempo_I_ver = polyval(p_1,linspace(min(absorbance_interpolar_1),max(absorbance_interpolar_1), 500)); %para ver no grafico
    
    hold on
    figure(1)
    plot(tempo_intepolar_1, absorbance_interpolar_1, '*r', tempo_I_ver, linspace(min(absorbance_interpolar_1),max(absorbance_interpolar_1), 500), '-r', S10_1_time, h_pico*max_absorbance, 'ob') %representa��o do 0.1*maxabs, do ponto exp + proximo e da f. interpola��o
    
    
end

% --- tudo igual mas para a outra metade da curva

% Determina�ao da exist�ncia de um ponto experimenal a 10 % da altura do pico (2� Parte)
n2 = 0;
for j2 = posicao_max_absorbance+1:length(Dados_Experimentais.tempo)
    if Dados_Experimentais.absorbance(j2) == h_pico*max_absorbance 
        S10_2_absorbance = Dados_Experimentais.absorbance(j2);
        S10_2_time = Dados_Experimentais.tempo(j2);
        n2 = n2 + 1;
    end 
end


if n2 ==0
    dif_absorb2 = abs(Dados_Experimentais.absorbance(posicao_max_absorbance+1:length(Dados_Experimentais.tempo))- h_pico*max_absorbance);
    min_dif_absorb2 = min(dif_absorb2);

    
    m2 = 0;
    
    for j2 = posicao_max_absorbance+1:length(Dados_Experimentais.tempo)
        m2 = m2+1;
        if min_dif_absorb2 == dif_absorb2(m2)
            absorbance_2 = Dados_Experimentais.absorbance(j2);
            tempo_2 = Dados_Experimentais.tempo(j2);
            posicao_min_2 = j2;
        end
    end
        
    if absorbance_2 < h_pico*max_absorbance
    
            absorbance_interpolar_2 = Dados_Experimentais.absorbance(posicao_min_2-npontos_int:posicao_min_2+(npontos_int-1));
            tempo_intepolar_2 = Dados_Experimentais.tempo(posicao_min_2-npontos_int:posicao_min_2+(npontos_int-1));
    
    elseif absorbance_2 > h_pico*max_absorbance
        
            absorbance_interpolar_2 = Dados_Experimentais.absorbance(posicao_min_2-(npontos_int-1):posicao_min_2+npontos_int);
            tempo_intepolar_2 = Dados_Experimentais.tempo(posicao_min_2-(npontos_int-1):posicao_min_2+npontos_int);
       
    end    
    
    p_2 = polyfit(absorbance_interpolar_2,tempo_intepolar_2,2);
    S10_2_time = polyval(p_2,h_pico*max_absorbance);
    
    tempo_II_ver = polyval(p_2,linspace(min(absorbance_interpolar_2),max(absorbance_interpolar_2), 500));
    
    hold on
    figure(1)
    subplot(1,3,1)
    plot(tempo_intepolar_2, absorbance_interpolar_2, '*r', tempo_II_ver, linspace(min(absorbance_interpolar_2),max(absorbance_interpolar_2), 500),'-r', S10_2_time, h_pico*max_absorbance, 'ob')
    xlabel('{\itt} (s)')
    ylabel('Absorbance (AU)')
end    
% ---    

% Express�o para o c�lculo da simetria do pico
S10 = (S10_2_time-time_max_absorbance)/(time_max_absorbance-S10_1_time);



if Veloc_Calc == 11
    
 % Determina��o da velocidade linear do solvente no interior da coluna (neste caso e espec�fico para CO2)
    CO2_dens_bomba = densityCO2(Dados_Experimentais.Pexp_bomba/1.01325, Dados_Experimentais.Texp_bomba); % (g/cm3)
    CO2_dens_estufa = densityCO2(Dados_Experimentais.Pexp/1.01325, Dados_Experimentais.Texp); % (g/cm3)
    
    caudal_mass_CO2 = Dados_Experimentais.caudal_bomba*CO2_dens_bomba; % (g/mL)
     
    Q_estufa = caudal_mass_CO2/CO2_dens_estufa; % (mL/min)
    
    Area_SR_coluna = pi*R_col^2; %(cm2)
    
    u0_exp = Q_estufa/60/Area_SR_coluna; % (cm/s)
    
elseif Veloc_Calc == 22
    
  % Determina��o da velocidade linear do solvente no interior da coluna (neste caso e espec�fico para CO2)
    CO2_dens_bolha = densityCO2(Dados_Experimentais.Pamb/1.01325, Dados_Experimentais.Tbolha); % (g/cm3)
    CO2_dens_estufa = densityCO2(Dados_Experimentais.Pexp/1.01325, Dados_Experimentais.Texp); % (g/cm3)
    
    caudal_mass_CO2 = Dados_Experimentais.caudal_bolha*CO2_dens_bolha; % (g/mL)
    
    Q_estufa = caudal_mass_CO2/CO2_dens_estufa; % (mL/min)
    
    Area_SR_coluna = pi*R_col^2; % (cm2)
   
    u0_exp = Q_estufa/60/Area_SR_coluna; % (cm/s)

elseif Veloc_Calc ~= 11 || Veloc_Calc ~= 22

return

end


%% Optimiza�ao
options = optimset;
% Modify options setting
options = optimset(options,'Display' ,'notify');
options = optimset(options,'MaxIter' ,50000);
options = optimset(options,'MaxFunEvals' ,100000);
options = optimset(options,'TolX' ,1e-10);
options = optimset(options,'TolFun',1e-15);
options = optimset(options,'Diagnostics' ,'on');

% Par�metros iniciais (U e a)
% params_ini = [0.2  0.5];
t0=L_col/u0_exp; %s
tmax = find(Dados_Experimentais.absorbance==max_absorbance); %indice correspondente � abs maxima
k0 = (Dados_Experimentais.tempo(tmax)-t0)/t0; %ko estimado
U0 = u0_exp/(k0+1); %estimativa inicial
params_ini = [U0  0.5];


[params_optim, fval,exitflag,output] = fminsearch(@(params) D12_optim_CIR_fun(params, S10_1_time, S10_2_time, Metodo_Integral), params_ini, options);

epsilon = fval;

U_optim = params_optim(1);
a_optim = params_optim(2);


% Dois tempos que definem o intervalo de integra�ao
t1 = S10_1_time; % tempo para o qual a altura do pico e 10% do valor m�ximo
t2 = S10_2_time; % tempo para o qual a altura do pico e 10% do valor m�ximo (NOTA: t1<t2)


n = 0;

for i = 1:length(Dados_Experimentais.tempo)
    if Dados_Experimentais.tempo(i)>= t1 & Dados_Experimentais.tempo(i)<= t2 %se estiverem entre t1 e t2 (acima 10% altura pico)
    n = n+1;
    time_integral(n) = Dados_Experimentais.tempo(i); %tempos usados para fzer o integral
    absorbance_integral(n) = Dados_Experimentais.absorbance(i); %absorbancias usadas para fzer o integral
    end
end


aviso1 = 0;
aviso2 = 0;

for p = 1:length(Dados_Experimentais.tempo)
    if t1 == Dados_Experimentais.tempo(p)
       aviso1 = 1; % Significa que existe ponto experimental;
    elseif t2 ==  Dados_Experimentais.tempo(p)
       aviso2 = 1; % Significa que existe ponto experimental;
    end
end


%% - Resultados a partir do c�lculo da velocidade linear a partir das condi��es de opera�� da BOMBA 
  
k_optim_CIR_1 = u0_exp/U_optim - 1;

variavel_1 = (1 + 6*k_optim_CIR_1 + 11*k_optim_CIR_1^2)/(1 + k_optim_CIR_1);
variavel_2 = (1 + 6*k_optim_CIR_1 + 11*k_optim_CIR_1^2)/(1 + k_optim_CIR_1)^2;

termo_D12_1 = variavel_1*(R_col^2*U_optim^2)/(24*a_optim);
termo_D12_2 = 1 + (1 - variavel_2*R_col^2*U_optim^2/(12*a_optim^2))^(1/2);

D12_optim_CIR = termo_D12_1/termo_D12_2; % Valor de Difusividade (cm2/s)


Condicao_11 = (1 + 6*k_optim_CIR_1 + 11*k_optim_CIR_1^2).*(R_col^2*U_optim^2)/(48*D12_optim_CIR^2); % NOTA: Tem de ser >> 1

Resultados_CIR_fitt.S10 = S10;
Resultados_CIR_fitt.u0_exp = u0_exp; % (cm/s) velocidade linear dentro da coluna - calculadas a partir das condi�oes operatorias da bomba
Resultados_CIR_fitt.max_absorbance=max_absorbance;
Resultados_CIR_fitt.D12 = D12_optim_CIR; % (cm^2/s) Valor da difusividade experimental
Resultados_CIR_fitt.k_optim = k_optim_CIR_1; % coeficiente de reten��o 
Resultados_CIR_fitt.epsilon = epsilon*100;
Resultados_CIR_fitt.a = a_optim;
Resultados_CIR_fitt.U = U_optim;


ana = '';


%%
if aviso1 == 1 && aviso2 == 1
    
    % Pontos experimenais qe entram no c�lculo do integral entre os valores t1 e t2
    absorbance_integral_final = absorbance_integral;
    time_integral_final = time_integral;
    
       
elseif aviso1 == 1 &&  aviso2 ~= 1
    
    % Pontos necess�rios para calcular o integral entre os valores t1 e t2
    absorbance_integral_final = [absorbance_integral h_pico*max(Dados_Experimentais.absorbance)];
    time_integral_final = [time_integral t2];
        
elseif aviso1 ~= 1 &&  aviso2 == 1
    
    % Pontos necess�rios para calcular o integral entre os valores t1 e t2   
    absorbance_integral_final = [h_pico*max(Dados_Experimentais.absorbance) absorbance_integral];
    time_integral_final = [t1 time_integral];
        
elseif aviso1 ~= 1 &&  aviso2 ~= 1
    
    % Pontos necess�rios para calcular o integral entre os valores t1 e t2
   absorbance_integral_final = [h_pico*max(Dados_Experimentais.absorbance) absorbance_integral h_pico*max(Dados_Experimentais.absorbance)];
   time_integral_final = [t1 time_integral t2];
   
end 
    

if Metodo_Integral == 1
    % Curva de absorv�ncia experimental (normalizada) COMPLETA 
    Area_Absorbance = trapz(time_integral_final,absorbance_integral_final); % (AU*s) -> M�todo dos trap�zios
    Absor_norm = Dados_Experimentais.absorbance./Area_Absorbance; % [s^(-1)]   
    
    % Normaliza��o da curva calculada COMPLETA [(normaliza��o segundo a express�o (Eq. 24) que � fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
    Ca_norm_num = 1./(4*pi*a_optim.*time_integral_final).^(1/2).*exp(-(L_col-U_optim.*time_integral_final).^2./(4.*a_optim.*time_integral_final)); % cm^(-1)
    % �rea da express�o anterior - integral da express�o anterior em ordem ao tempo.
    Area_Ca_norm_denom = trapz(time_integral_final,Ca_norm_num); % (s.cm^(-1)) -> M�todo dos trap�zios
    

elseif Metodo_Integral == 2 % Simpson's 1/3 rule
    % Curva de absorv�ncia experimental (normalizada) COMPLETA  
    intervalo_t = linspace(min(time_integral_final), max(time_integral_final), 3000);
    pp = spline(time_integral_final,absorbance_integral_final);
    Interp_abs_exp_det = ppval(pp, intervalo_t);
    
    Area_Absorbance = Simpson_Method(intervalo_t,Interp_abs_exp_det,[],'1/3'); % Simpson's 3/8 rule
    Absor_norm = Dados_Experimentais.absorbance./Area_Absorbance; % [s^(-1)]
    
    % Normaliza��o da curva calculada COMPLETA [(normaliza��o segundo a express�o (Eq. 24) que � fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
    Ca_norm_num = 1./(4*pi*a_optim.*intervalo_t).^(1/2).*exp(-(L_col-U_optim.*intervalo_t).^2./(4.*a_optim.*intervalo_t)); % cm^(-1)
    % �rea da express�o anterior - integral da express�o anterior em ordem ao tempo.
    Area_Ca_norm_denom = Simpson_Method(intervalo_t,Ca_norm_num,[],'1/3'); % Simpson's 1/3 rule
        

elseif Metodo_Integral == 3  % Simpson's 3/8 rule
    intervalo_t = linspace(min(time_integral_final), max(time_integral_final), 3000);
    pp = spline(time_integral_final,absorbance_integral_final);
    Interp_abs_exp_det = ppval(pp, intervalo_t);
    
    Area_Absorbance = Simpson_Method(intervalo_t,Interp_abs_exp_det,[],'3/8'); % Simpson's 3/8 rule
    Absor_norm = Dados_Experimentais.absorbance./Area_Absorbance; % [s^(-1)] 
    
    % Normaliza��o da curva calculada COMPLETA [(normaliza��o segundo a express�o (Eq. 24) que � fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
    Ca_norm_num = 1./(4*pi*a_optim.*intervalo_t).^(1/2).*exp(-(L_col-U_optim.*intervalo_t).^2./(4.*a_optim.*intervalo_t)); % cm^(-1)
    % �rea da express�o anterior - integral da express�o anterior em ordem ao tempo.
    Area_Ca_norm_denom = Simpson_Method(intervalo_t,Ca_norm_num,[],'3/8'); % Simpson's 3/8 rule

    
end


% Representa��o da concentra��o calculada, mas em termos de absorv�ncia, com mais detalhe (uso de mais pontos em t) 
tempo_detalhado = min(Dados_Experimentais.tempo):1:max(Dados_Experimentais.tempo);
Ca_norm_num_alldata = 1./(4*pi*a_optim.*tempo_detalhado).^(1/2).*exp(-(L_col-U_optim.*tempo_detalhado).^2./(4.*a_optim.*tempo_detalhado)); % g/cm3 -> pela express�o (Eq. 25) dada no artigo Kong et al (J. Chromatogr. A, 1035, 177)  

Absorvance_calc = Ca_norm_num_alldata.*Area_Absorbance./Area_Ca_norm_denom; % (AU)



 figure(1)
 subplot(1,3,3)
plot(Dados_Experimentais.tempo./60, Dados_Experimentais.absorbance, 'ok', tempo_detalhado/60, Absorvance_calc, '-r')
xlabel('{\itt} (s)')
ylabel('Absorbance (AU)')
%% ADICIONAR
title([Dados_Experimentais.soluto, ' in ', Dados_Experimentais.solvente, ' - ', ...
    num2str(Dados_Experimentais.Texp), ' K / ', num2str(Dados_Experimentais.Pexp/10), ' MPa (\lambda = ',...
    num2str(Dados_Experimentais.wavelength),' nm)'])
%% ------------------------------------------------
legend1 = legend('Experimental', ['Fitting: {\itD}_{12} = ', num2str(D12_optim_CIR*10^4, 5), '\times10^{-4} (cm^2/s)', ', \epsilon = ', num2str(fval*100,4), '%']);

set(legend1,'EdgeColor',[1 1 1]);
%hold on


Resultados_CIR_fitt.Area_Absorbance = Area_Absorbance;


Resultados_CIR_fitt;

%% An�lise da Largura do pico

% Determina�ao de 60.7% da altura do pico em ambas as partes (1� parte) 
p1 = 0;
for q1 = 1:posicao_max_absorbance
    if Dados_Experimentais.absorbance(q1) == h_abertuta_pico*max_absorbance
        W_time_607_1 = Dados_Experimentais.absorbance(q1);  
        p1 = p1 + 1;
   end
end

if p1 == 0
    
    difernca_607_1 = abs(Dados_Experimentais.absorbance(1:posicao_max_absorbance)- h_abertuta_pico*max_absorbance);
    min_difernca_607_1 = min(difernca_607_1);

     for q1 = 1:posicao_max_absorbance
        if min_difernca_607_1 == difernca_607_1(q1)
            absorbance_607_1 = Dados_Experimentais.absorbance(q1);
            tempo_607_1 = Dados_Experimentais.tempo(q1);
            posicao_min_607_1 = q1;
        end
     end   
    
    if absorbance_607_1 < h_abertuta_pico*max_absorbance
       absorbance_interpolar_607_1 = Dados_Experimentais.absorbance(posicao_min_607_1-2:posicao_min_607_1+3);
       tempo_intepolar_607_1 = Dados_Experimentais.tempo(posicao_min_607_1-2:posicao_min_607_1+3);
    
    elseif absorbance_607_1 > h_abertuta_pico*max_absorbance
        absorbance_interpolar_607_1 = Dados_Experimentais.absorbance(posicao_min_607_1-3:posicao_min_607_1+2);
        tempo_intepolar_607_1 = Dados_Experimentais.tempo(posicao_min_607_1-3:posicao_min_607_1+2);
    end    
     
    p_607_1 = polyfit(absorbance_interpolar_607_1,tempo_intepolar_607_1,3);
    time_607_1 = polyval(p_607_1,h_abertuta_pico*max_absorbance);

    tempo_607_ver_1 = polyval(p_607_1, linspace(min(absorbance_interpolar_607_1),max(absorbance_interpolar_607_1), 500));

%     figure(1)
    subplot(2,2,1); plot(time_607_1, h_abertuta_pico*max_absorbance, 'sg',tempo_607_ver_1,linspace(min(absorbance_interpolar_607_1),max(absorbance_interpolar_607_1), 500), '-g')
    
end 

% Determina�ao de 60.7% da altura do pico em ambas as partes (2� parte) 

p2 = 0;
for q2 = posicao_max_absorbance+1:length(Dados_Experimentais.tempo)
    if Dados_Experimentais.absorbance(q2) == h_abertuta_pico*max_absorbance
        W_time_607_2 = Dados_Experimentais.absorbance(q2);  
        p2 = p2 + 1;
   end
end


if p2 == 0
    
    difernca_607_2 = abs(Dados_Experimentais.absorbance(posicao_max_absorbance+1:length(Dados_Experimentais.tempo))- h_abertuta_pico*max_absorbance);
    min_difernca_607_2 = min(difernca_607_2);

    m3 = 0;
     for q2 = posicao_max_absorbance+1:length(Dados_Experimentais.tempo)
         
        m3 = m3 + 1;
        if min_difernca_607_2 == difernca_607_2(m3)
            absorbance_607_2 = Dados_Experimentais.absorbance(q2);
            tempo_607_2 = Dados_Experimentais.tempo(q2);
            posicao_min_607_2 = q2;
        end
     end   
    

    if absorbance_607_2 < h_abertuta_pico*max_absorbance
       absorbance_interpolar_607_2 = Dados_Experimentais.absorbance(posicao_min_607_2-3:posicao_min_607_2+2);
       tempo_intepolar_607_2 = Dados_Experimentais.tempo(posicao_min_607_2-3:posicao_min_607_2+2);
    
    elseif absorbance_607_2 > h_abertuta_pico*max_absorbance
        absorbance_interpolar_607_2 = Dados_Experimentais.absorbance(posicao_min_607_2-2:posicao_min_607_2+3);
        tempo_intepolar_607_2 = Dados_Experimentais.tempo(posicao_min_607_2-2:posicao_min_607_2+3);
    end    
     
    p_607_2 = polyfit(absorbance_interpolar_607_2,tempo_intepolar_607_2,3);
    time_607_2 = polyval(p_607_2,h_abertuta_pico*max_absorbance);
    

    tempo_607_ver_2 = polyval(p_607_2, linspace(min(absorbance_interpolar_607_2),max(absorbance_interpolar_607_2), 500));

%     figure(1)
    hold on
    
    subplot(2,2,1); plot(time_607_2, h_abertuta_pico*max_absorbance, 'sg',tempo_607_ver_2,linspace(min(absorbance_interpolar_607_2),max(absorbance_interpolar_607_2), 500), '-g')
    
end 


W_time_607 = (time_607_2-time_607_1)/2;
S10_607 = (time_607_2-time_max_absorbance)/(time_max_absorbance-time_607_1);


% NOTA neste caso uso a velocidade optimizada, mas devo sempre usar a % velocidade experimental obtida
Analise_pico.u0_exp_soluto = L_col/time_max_absorbance;

Analise_pico.H = Analise_pico.u0_exp_soluto ^2*W_time_607^2/(fator_h*L_col);
% Analise_pico.H = u0_exp^2*W_time_607^2/L_col/5.545;

Analise_pico.k = -1 + u0_exp/Analise_pico.u0_exp_soluto;

%sem difus�o no solido
Analise_pico.D12_calc = u0_exp/12/(Analise_pico.k+1) * ( 3*Analise_pico.H*(1+Analise_pico.k) - 3^(1/2)*( Analise_pico.H^2*(3*Analise_pico.k^2 +6*Analise_pico.k +3)...
    -11*R_col^2*Analise_pico.k^2 - 6*R_col^2*Analise_pico.k-R_col^2  )^(1/2) ); 

%com difus�o no solido

Analise_pico.H_s = Analise_pico.H-2/3*espess_revestimento^2*u0_exp*Analise_pico.k/Dl/(1+Analise_pico.k)^2;
Analise_pico.D12_calc_c_s = u0_exp/12/(Analise_pico.k+1) * ( 3*Analise_pico.H_s*(1+Analise_pico.k) - 3^(1/2)*( Analise_pico.H_s^2*(3*Analise_pico.k^2 +6*Analise_pico.k +3)...
    -11*R_col^2*Analise_pico.k^2 - 6*R_col^2*Analise_pico.k-R_col^2  )^(1/2) ); 

Analise_pico.u0_opt = Analise_pico.D12_calc/R_col*sqrt(48)*(1+Analise_pico.k)/(11*Analise_pico.k^2+ 6*Analise_pico.k+1)^(1/2);
Analise_pico.erro = (Analise_pico.D12_calc-Resultados_CIR_fitt.D12)/Resultados_CIR_fitt.D12*100;

% Compila��o da Informa��o (M�todo an�lise da Largura do Pico) 
D12Calc_variance.S10_607 = S10_607;
D12Calc_variance.H = Analise_pico.H;
D12Calc_variance.u0_exp_soluto= Analise_pico.u0_exp_soluto;
D12Calc_variance.u0_opt= Analise_pico.u0_opt;
D12Calc_variance.D12_calc = Analise_pico.D12_calc;
D12Calc_variance.k_calc=Analise_pico.k;
% D12Calc_variance.H_s=Analise_pico.H_s;
% D12Calc_variance.D12_calc_c_s=Analise_pico.D12_calc_c_s;
D12Calc_variance.erro = Analise_pico.erro;


Fim = '';



end 



%% Fun��o de Optimiza��o dos par�metros
function Fobj = D12_optim_CIR_fun(params, S10_1_time, S10_2_time, Metodo_Integral)

global R_col L_col Dados_Experimentais m_real_inject h_pico

U = params(1);
a = params(2);


% Dois tempos que definem o intervalo de integra�ao
t1 = S10_1_time; % tempo para o qual a altura do pico e 10% do valor m�ximo
t2 = S10_2_time; % tempo para o qual a altura do pico e 10% do valor m�ximo (NOTA: t1<t2)


n = 0;

for i = 1:length(Dados_Experimentais.tempo)
    if Dados_Experimentais.tempo(i)>= t1 & Dados_Experimentais.tempo(i)<= t2 %se estiver entre t1 e t2 8acima 10% da altura do pico)
    n = n+1;
    time_integral(n) = Dados_Experimentais.tempo(i); %tempos a usar para as areas/integrais
    absorbance_integral(n) = Dados_Experimentais.absorbance(i); %abs a usar para as areas/integrais
    end
end


aviso1 = 0;
aviso2 = 0;

for p = 1:length(Dados_Experimentais.tempo) %correr os dados que se vao usar  para o integral
    if t1 == Dados_Experimentais.tempo(p)
       aviso1 = 1; % Significa que existe ponto experimental;
    elseif t2 ==  Dados_Experimentais.tempo(p)
       aviso2 = 1; % Significa que existe ponto experimental;
    end
end



% Curva calculada (normalizada) 
% K = D12_optim + u0^2*R_col^2/(48*D12_optim);


if aviso1 == 1 && aviso2 == 1
    
    % Pontos experimenais qe entram no c�lculo do integral entre os valores t1 e t2
    absorbance_integral_final = absorbance_integral;
    time_integral_final = time_integral;
    
       
elseif aviso1 == 1 &&  aviso2 ~= 1 %se n existe o t2 experimntal, acrescenta-se esse ponto
    
    % Pontos necess�rios para calcular o integral entre os valores t1 e t2
    absorbance_integral_final = [absorbance_integral h_pico*max(Dados_Experimentais.absorbance)];
    time_integral_final = [time_integral t2];
        
elseif aviso1 ~= 1 &&  aviso2 == 1 %se n existe o t1 experimntal, acrescenta-se esse ponto
    
    % Pontos necess�rios para calcular o integral entre os valores t1 e t2   
    absorbance_integral_final = [h_pico*max(Dados_Experimentais.absorbance) absorbance_integral];
    time_integral_final = [t1 time_integral];
        
elseif aviso1 ~= 1 &&  aviso2 ~= 1 %se n existem o t1 e o t2 experimntais, acrescentam-se esses pontos
    
    % Pontos necess�rios para calcular o integral entre os valores t1 e t2
   absorbance_integral_final = [h_pico*max(Dados_Experimentais.absorbance) absorbance_integral h_pico*max(Dados_Experimentais.absorbance)];
   time_integral_final = [t1 time_integral t2];
   
end

%so pontos experimentais (sem os t1 e t2 "inventados" qd n existem):
absorbance_exp = absorbance_integral;
time_exp = time_integral;

% Curva de absorv�ncia experimental (normalizada)
if Metodo_Integral == 1
    Area_Absorbance = trapz(time_integral_final,absorbance_integral_final); % (AU*s)
elseif Metodo_Integral == 2
    Area_Absorbance = Simpson_Method(time_integral_final,absorbance_integral_final,[],'1/3');
elseif Metodo_Integral == 3
    Area_Absorbance = Simpson_Method(time_integral_final,absorbance_integral_final,[],'3/8');
end


 Absor_norm = absorbance_exp./Area_Absorbance; % [s^(-1)]   absorvancia (pts so acima de 10%) normalizada pela area da curva acima de 10% da altura do pico


% Curva de absorv�ncia "experimental" (normalizada) COMPLETA  
intervalo_t = linspace(min(time_integral_final), max(time_integral_final), 3000);
% pp = spline(time_integral_final,absorbance_integral_final); %spline cubica
% Interp_abs_exp_det = ppval(pp, intervalo_t); %Absorvancia interpolada para os 3000 pontos igualmente espa�ados entre t1 e t2


% Normaliza��o da curva calculada [(normaliza��o segundo a express�o (Eq. 24) que � fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
Ca_norm_num = 1./(4*pi*a.*time_exp).^(1/2).*exp(-(L_col-U.*time_exp).^2./(4*a.*time_exp)); % cm^(-1)
% �rea da express�o anterior - integral da express�o anterior em ordem ao tempo.


if Metodo_Integral == 1
       
    % Normaliza��o da curva calculada COMPLETA [(normaliza��o segundo a express�o (Eq. 24) que � fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
    Ca_norm_num_trap = 1./(4*pi*a.*intervalo_t).^(1/2).*exp(-(L_col-U.*intervalo_t).^2./(4.*a.*intervalo_t)); % cm^(-1)
    % �rea da express�o anterior - integral da express�o anterior em ordem ao tempo.
    Area_Ca_norm_denom = trapz(intervalo_t,Ca_norm_num_trap); % (s.cm^(-1)) -> M�todo dos trap�zios
        
elseif Metodo_Integral == 2
    
    Ca_norm_num_S13 = 1./(4*pi*a.*intervalo_t).^(1/2).*exp(-(L_col-U.*intervalo_t).^2./(4*a.*intervalo_t)); % cm^(-1)
    
    Area_Ca_norm_denom = Simpson_Method(intervalo_t,Ca_norm_num_S13,[],'1/3'); % Simpson's 3/8 rule

elseif Metodo_Integral == 3
        
    Ca_norm_num_S38 = 1./(4*pi*a.*intervalo_t).^(1/2).*exp(-(L_col-U.*intervalo_t).^2./(4*a.*intervalo_t)); % cm^(-1)
    
    Area_Ca_norm_denom = Simpson_Method(intervalo_t,Ca_norm_num_S38,[],'3/8'); % Simpson's 3/8 rule
    
end
                                                    
% Normaliza��o da express�o da concentra��o (�rea normalizada)
Ca_norm = Ca_norm_num/Area_Ca_norm_denom; % [s^(-1)] (Ca_numerador so para os mm t_exp entre 1 e t2; area � q foi calculda com 3000 pontos entre t1 e t2)

figure(1)
subplot(1,3,2)
plot(time_exp, Absor_norm, '*k', time_exp, Ca_norm, '-k')
xlabel('{\itt} (s)')
ylabel('{\itC}_{2,exp,norm}, {\itC}_{2,calc,norm} (s^{-1})')


% Defini��o da fun��o objectivo
Termo1 = (Ca_norm-Absor_norm).^2;
Termo2 = Absor_norm.^2;


Fobj1 = trapz(time_exp,Termo1); %integral ou area da fun�ao "termo 1" entre t1 e t2 (ou o ponto + proximo)
Fobj2 = trapz(time_exp, Termo2); %integral ou area da fun�ao "termo 2" entre t1 e t2 (ou o ponto + proximo)

Fobj = (Fobj1/Fobj2)^(1/2);


end



