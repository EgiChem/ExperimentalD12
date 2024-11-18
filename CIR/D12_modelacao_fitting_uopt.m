function D12_modelacao_fitting()


clear all 
clc
format short

global R_col L_col Dados_Experimentais perc time_max_absorbance

% m_real_inject
%% Programa de cálculo das difusividades a partir dos pontos experimentais (por fitting e por análise da largura do pico ou "método dos momentos")

% NOTA IMPORTANTE : Escolha do método de cálculo das áreas: 
%                       - 1 para o método dos trapézios,
%                       - 2 para o método de Simpson 1/3,
%                       - 3 para o método de Simpson 3/8.


Metodo_Integral = 1; 

percentagem = 10; %
perc = percentagem/100;
%% 


coluna = 3;
%% Dados da coluna

% L_col = 17.681*100; % Comprimento da coluna (cm)
d_enrolamento = 30.0; % Diâmetro do enrolamento (cm)

% diam_col = (0.516)/10; % Diâmetro interno da coluna (cm) 
% R_col = (diam_col)/2; %(cm) Raio da coluna

if coluna == 1
    % L_col = 16.887*100; % Comprimento da coluna (cm)
    L_col = 2125; % Comprimento da coluna (cm)
elseif coluna == 3
%         L_col = 1035; 
    L_col = 1024.3; %cm atualizado com ultimo corte a 23/04/2014
end   
% diam_col = (0.53)/10; % Diâmetro interno da coluna (cm) 
diam_col = 0.0522; % Diâmetro interno da coluna (cm) 
% diam_col = 0.0556; % Diâmetro interno da coluna (cm) 
% espess_revestimento = 1E-4; % (cm)
% espess_revestimento = 0; % (cm) 
R_col = (diam_col)/2; %(cm) Raio da coluna

%% Informação acerca do soluto (massa, ou volume, ou concentração do soluto injectado

vol_inject = 0.1*1E-3; % mL
%%
percentagem = 10; %
perc = percentagem/100;

%% Leitura dos dados de um ficheiro de dados em txt (tempo vs Absorvância) 

[data_matrix, text_data, alldata] = xlsread(...
    '2015_02_24ethanol_em_CO2_250bar_50ºC_coluna3.xlsx',8); 

% Q_dado = 0.13;
Dados_Experimentais.soluto = text_data{2,5};
Dados_Experimentais.solvente = text_data{2,4};
Dados_Experimentais.solvente_formula = text_data{3,4};
%%
Dados_Experimentais.Pexp = data_matrix(1,6); %(bar)
Dados_Experimentais.Pamb = data_matrix(1,10); %(bar)
Dados_Experimentais.Texp = data_matrix(1,7)+ 273.15; % (K)
Dados_Experimentais.Tbolha = data_matrix(1,11)+ 273.15; % (K)
Dados_Experimentais.wavelength = data_matrix(1,3); % (nm)
Dados_Experimentais.Texp_bomba = data_matrix(1,8)+ 273.15; % (K)
Dados_Experimentais.caudal_bomba = data_matrix(1,9); % (mL/min)
% Dados_Experimentais.caudal_bomba = Q_dado;
Dados_Experimentais.caudal_bolha = data_matrix(1,12); % (mL/min)
% Dados_Experimentais.tempo =  double(uint32(data_matrix(1:end,1))); % tempo (s)
Dados_Experimentais.tempo = data_matrix(1:end,1); % tempo (s)
Dados_Experimentais.absorbance = data_matrix(1:end,2); % Unidades de Absorvância (AU)


%% Determinação da velocidade linear do solvente no interior da coluna (neste caso e específico para CO2)
CO2_dens_bolha = densityCO2(Dados_Experimentais.Pamb/1.01325, Dados_Experimentais.Tbolha); % (g/cm3)
CO2_dens_bomba = densityCO2(Dados_Experimentais.Pexp/1.01325, Dados_Experimentais.Texp_bomba); % (g/cm3)
CO2_dens_estufa = densityCO2(Dados_Experimentais.Pexp/1.01325, Dados_Experimentais.Texp); % (g/cm3)


caudal_mass_CO2 = Dados_Experimentais.caudal_bomba*CO2_dens_bomba; %g/mL
caudal_mass_CO2_2 = Dados_Experimentais.caudal_bolha*CO2_dens_bolha; %g/mL


Q_estufa = caudal_mass_CO2/CO2_dens_estufa; %mL/min
Q_estufa_2 = caudal_mass_CO2_2/CO2_dens_estufa; %mL/min

Area_SR_coluna = pi*R_col^2; %cm2


u0_exp_bomba = Q_estufa/60/Area_SR_coluna %cm/s (determinada a partir do caudal na bomba com Tbanho)
u0_exp_2 = Q_estufa_2/60/Area_SR_coluna %cm/s (determinada pelo medidor de bolha de sabão)






% Figura 1 - Absorvância vs tempo
figure(1)
plot(Dados_Experimentais.tempo, Dados_Experimentais.absorbance, '*k')
xlabel('{\itt} (s)')
ylabel('Absorbance (AU)')
Title([Dados_Experimentais.soluto, ' in ', Dados_Experimentais.solvente_formula, ' - ', num2str(Dados_Experimentais.Texp), ...
    ' K / ', num2str(Dados_Experimentais.Pexp/10), ' MPa (\lambda = ',...
    num2str(Dados_Experimentais.wavelength),' nm)'])



%% Tratamento dos dados experimentais

% Deteção do tempo correpondente à absorvância máxima

max_absorbance = max(Dados_Experimentais.absorbance); 

for i = 1:length(Dados_Experimentais.tempo)
    if Dados_Experimentais.absorbance(i) == max_absorbance
        time_max_absorbance = Dados_Experimentais.tempo(i);
        posicao_max_absorbance = i;
    end 
end


%% Verificaçao da simetria do pico obtido (através do parâmetro S10)


 
% Determinaçao da existência de um ponto experimenal a 10 % da altura do pico (1ª Parte)
% ---(indice correspondente á abs max)
n1 = 0;
for j1 = 1:posicao_max_absorbance
    if Dados_Experimentais.absorbance(j1) == perc*max_absorbance 
        S10_1_absorbance = Dados_Experimentais.absorbance(j1);
        S10_1_time = Dados_Experimentais.tempo(j1);
        n1 = n1 + 1;
    end 
end
% ---
npontos_int = 3; %pontos usados para cado lado na interpolação

if n1 ==0 %se não existir nenhum ponto correspondente a 0.1*maxabs
    dif_absorb1 = abs(Dados_Experimentais.absorbance(1:posicao_max_absorbance)- perc*max_absorbance); %diferença entre cada ponto e o maxabs
    min_dif_absorb1 = min(dif_absorb1);
    
% ---(indice correspondente á diferença minima, ie, ponto + proximo de 0.1*maxabs)    
    for j1 = 1:posicao_max_absorbance
        if min_dif_absorb1 == dif_absorb1(j1)
            absorbance_1 = Dados_Experimentais.absorbance(j1);
            tempo_1 = Dados_Experimentais.tempo(j1);
            posicao_min_1 = j1;
        end
    end
% ---
    
    if absorbance_1 < perc*max_absorbance % se o ponto + proximo esta abaixo
    
            absorbance_interpolar_1 = Dados_Experimentais.absorbance(posicao_min_1-(npontos_int-1):posicao_min_1+npontos_int); %pontos a usar na função interpolação: 2 para trás e 3 para a frente
            tempo_intepolar_1 = Dados_Experimentais.tempo(posicao_min_1-(npontos_int-1):posicao_min_1+npontos_int); %igual para o y (tempo)
    
    elseif absorbance_1 > perc*max_absorbance % se o ponto + proximo esta acima
        
            absorbance_interpolar_1 = Dados_Experimentais.absorbance(posicao_min_1-npontos_int:posicao_min_1+(npontos_int-1)); %pontos a usar na funçao de interpolaçao: 3 para trás e 2 para a frente
            tempo_intepolar_1 = Dados_Experimentais.tempo(posicao_min_1-npontos_int:posicao_min_1+(npontos_int-1));%igual para o y (tempo)
       
    end
    
    p_1 = polyfit(absorbance_interpolar_1,tempo_intepolar_1,2); %função de interpolação: polinómio de 2º grau
    S10_1_time = polyval(p_1,perc*max_absorbance); %tempo correspondente a 0.1*maxabs
    
    tempo_I_ver = polyval(p_1,linspace(min(absorbance_interpolar_1),max(absorbance_interpolar_1), 500)); %para ver no grafico
    
    hold on
    figure(1)
    plot(tempo_intepolar_1, absorbance_interpolar_1, '*r', tempo_I_ver, linspace(min(absorbance_interpolar_1),max(absorbance_interpolar_1), 500), '-r', S10_1_time, perc*max_absorbance, 'ob') %representação do 0.1*maxabs, do ponto exp + proximo e da f. interpolação
    
    
end

% --- tudo igual mas para a outra metade da curva

% Determinaçao da existência de um ponto experimenal a 10 % da altura do pico (2ª Parte)
n2 = 0;
for j2 = posicao_max_absorbance+1:length(Dados_Experimentais.tempo)
    if Dados_Experimentais.absorbance(j2) == perc*max_absorbance 
        S10_2_absorbance = Dados_Experimentais.absorbance(j2);
        S10_2_time = Dados_Experimentais.tempo(j2);
        n2 = n2 + 1;
    end 
end


if n2 ==0
    dif_absorb2 = abs(Dados_Experimentais.absorbance(posicao_max_absorbance+1:length(Dados_Experimentais.tempo))- perc*max_absorbance);
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
        
    if absorbance_2 < perc*max_absorbance
    
            absorbance_interpolar_2 = Dados_Experimentais.absorbance(posicao_min_2-npontos_int:posicao_min_2+(npontos_int-1));
            tempo_intepolar_2 = Dados_Experimentais.tempo(posicao_min_2-npontos_int:posicao_min_2+(npontos_int-1));
    
    elseif absorbance_2 > perc*max_absorbance
        
            absorbance_interpolar_2 = Dados_Experimentais.absorbance(posicao_min_2-(npontos_int-1):posicao_min_2+npontos_int);
            tempo_intepolar_2 = Dados_Experimentais.tempo(posicao_min_2-(npontos_int-1):posicao_min_2+npontos_int);
       
    end    
    
    p_2 = polyfit(absorbance_interpolar_2,tempo_intepolar_2,2);
    S10_2_time = polyval(p_2,perc*max_absorbance);
    
    tempo_II_ver = polyval(p_2,linspace(min(absorbance_interpolar_2),max(absorbance_interpolar_2), 500));
    
    hold on
    figure(1)
    plot(tempo_intepolar_2, absorbance_interpolar_2, '*r', tempo_II_ver, linspace(min(absorbance_interpolar_2),max(absorbance_interpolar_2), 500),'-r', S10_2_time, perc*max_absorbance, 'ob')
    xlabel('{\itt} (s)')
    ylabel('Absorbance (AU)')
end    
% ---    

% Expressão para o cálculo da simetria do pico
S10 = (S10_2_time-time_max_absorbance)/(time_max_absorbance-S10_1_time);



%% Optimizaçao
options = optimset;
% Modify options setting
options = optimset(options,'Display' ,'notify');
options = optimset(options,'MaxIter' ,50000);
options = optimset(options,'MaxFunEvals' ,100000);
options = optimset(options,'TolX' ,1e-10);
options = optimset(options,'TolFun',1e-15);
options = optimset(options,'Diagnostics' ,'on');

% Parâmetros iniciais (U e a)55
% params_ini = [0.2  0.5];
% t0=L_col/u0_exp; %s
% tmax = find(Dados_Experimentais.absorbance==max_absorbance); %indice correspondente à abs maxima
% k0 = (Dados_Experimentais.tempo(tmax)-t0)/t0; %ko estimado
% U0 = u0_exp/(k0+1); %estimativa inicial
U0 = L_col/time_max_absorbance;
params_ini = [U0  0.2];


[params_optim, fval,exitflag,output] = fminsearch(@(params) D12_optim_CIR(params, S10_1_time, S10_2_time, Metodo_Integral), params_ini, options);



ana = 'toto';


u0_optim = U0;
u0_optim = params_optim(1);
K_optim = params_optim(2);


% Dois tempos que definem o intervalo de integraçao
t1 = S10_1_time; % tempo para o qual a altura do pico e 10% do valor máximo
t2 = S10_2_time; % tempo para o qual a altura do pico e 10% do valor máximo (NOTA: t1<t2)


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



% Curva calculada (normalizada) 
S10
epsilon = fval*100;
D12_optim = (K_optim - sqrt(K_optim^2-R_col^2*u0_optim^2/12))/2;
t=time_max_absorbance/60
u0_exp_bomba;
L1=u0_exp_bomba*time_max_absorbance;
u0_optim
L2=u0_optim*time_max_absorbance;

a='';
% k_optim_CIR = 0;

% k_optim_CIR = u0_exp/U_optim - 1
% 
% variavel_1 = (1 + 6*k_optim_CIR + 11*k_optim_CIR^2)/(1 + k_optim_CIR);
% variavel_2 = (1 + 6*k_optim_CIR + 11*k_optim_CIR^2)/(1 + k_optim_CIR)^2;
% 
% 
% termo_D12_1 = variavel_1*(R_col^2*U_optim^2)/(24*a_optim);
% termo_D12_2 = 1 + (1 - variavel_2*R_col^2*U_optim^2/(12*a_optim^2))^(1/2);
% 
% D12_optim_CIR = termo_D12_1/termo_D12_2
% 
% (1+6*k_optim_CIR+11*k_optim_CIR^2)*R_col^2*U_optim^2/(48*D12_optim_CIR^2)
% disp('>>1?')
% 
% S10
% 
% epsilon = fval.*100
% 
% u0_exp

%com a outra velocidade:

% k_optim_CIR_2 = u0_exp_2/U_optim - 1;
% 
% 
% variavel_1_2 = (1 + 6*k_optim_CIR_2 + 11*k_optim_CIR_2^2)/(1 + k_optim_CIR_2);
% variavel_2_2 = (1 + 6*k_optim_CIR_2 + 11*k_optim_CIR_2^2)/(1 + k_optim_CIR_2)^2;
% 
% 
% termo_D12_1_2 = variavel_1_2*(R_col^2*U_optim^2)/(24*a_optim);
% termo_D12_2_2 = 1 + (1 - variavel_2_2*R_col^2*U_optim^2/(12*a_optim^2))^(1/2);
% 
% D12_optim_CIR_2 = termo_D12_1_2/termo_D12_2_2;
% 
% 
% u0_exp_2;

% Condição 
% Condicao_1 = (1 + 6*k_optim_CIR + 11*k_optim_CIR^2).*(R_col^2*U_optim^2)/(48*D12_optim_CIR^2); % NOTA: Tem de ser >> 1

%% restrições


%%
if aviso1 == 1 && aviso2 == 1
    
    % Pontos experimenais qe entram no cálculo do integral entre os valores t1 e t2
    absorbance_integral_final = absorbance_integral;
    time_integral_final = time_integral;
    
       
elseif aviso1 == 1 &&  aviso2 ~= 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2
    absorbance_integral_final = [absorbance_integral perc*max(Dados_Experimentais.absorbance)];
    time_integral_final = [time_integral t2];
        
elseif aviso1 ~= 1 &&  aviso2 == 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2   
    absorbance_integral_final = [perc*max(Dados_Experimentais.absorbance) absorbance_integral];
    time_integral_final = [t1 time_integral];
        
elseif aviso1 ~= 1 &&  aviso2 ~= 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2
   absorbance_integral_final = [perc*max(Dados_Experimentais.absorbance) absorbance_integral perc*max(Dados_Experimentais.absorbance)];
   time_integral_final = [t1 time_integral t2];
   
end 
    

if Metodo_Integral == 1
    % Curva de absorvância experimental (normalizada) COMPLETA 
    Area_Absorbance = trapz(time_integral_final,absorbance_integral_final); % (AU*s) -> Método dos trapézios
    Absor_norm = Dados_Experimentais.absorbance./Area_Absorbance; % [s^(-1)]   
    
    % Normalização da curva calculada COMPLETA [(normalização segundo a expressão (Eq. 24) que é fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
    Ca_norm_num = 1./(4*pi*K_optim.*time_integral_final).^(1/2).*exp(-(L_col-u0_optim.*time_integral_final).^2./(4.*K_optim.*time_integral_final)); % cm^(-1)
    % Área da expressão anterior - integral da expressão anterior em ordem ao tempo.
    Area_Ca_norm_denom = trapz(time_integral_final,Ca_norm_num); % (s.cm^(-1)) -> Método dos trapézios
    

elseif Metodo_Integral == 2 % Simpson's 1/3 rule
    % Curva de absorvância experimental (normalizada) COMPLETA  
    intervalo_t = linspace(min(time_integral_final), max(time_integral_final), 3000);
    pp = spline(time_integral_final,absorbance_integral_final);
    Interp_abs_exp_det = ppval(pp, intervalo_t);
    
    Area_Absorbance = Simpson_Method(intervalo_t,Interp_abs_exp_det,[],'1/3'); % Simpson's 3/8 rule
    Absor_norm = Dados_Experimentais.absorbance./Area_Absorbance; % [s^(-1)]
    
    % Normalização da curva calculada COMPLETA [(normalização segundo a expressão (Eq. 24) que é fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
    Ca_norm_num = 1./(4*pi*K_optim.*time_integral_final).^(1/2).*exp(-(L_col-u0_optim.*time_integral_final).^2./(4.*K_optim.*time_integral_final)); % cm^(-1)
    % Área da expressão anterior - integral da expressão anterior em ordem ao tempo.
    Area_Ca_norm_denom = Simpson_Method(intervalo_t,Ca_norm_num,[],'1/3'); % Simpson's 1/3 rule
        

elseif Metodo_Integral == 3  % Simpson's 3/8 rule
    intervalo_t = linspace(min(time_integral_final), max(time_integral_final), 3000);
    pp = spline(time_integral_final,absorbance_integral_final);
    Interp_abs_exp_det = ppval(pp, intervalo_t);
    
    Area_Absorbance = Simpson_Method(intervalo_t,Interp_abs_exp_det,[],'3/8'); % Simpson's 3/8 rule
    Absor_norm = Dados_Experimentais.absorbance./Area_Absorbance; % [s^(-1)] 
    
    % Normalização da curva calculada COMPLETA [(normalização segundo a expressão (Eq. 24) que é fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
    Ca_norm_num = 1./(4*pi*K_optim.*time_integral_final).^(1/2).*exp(-(L_col-u0_optim.*time_integral_final).^2./(4.*K_optim.*time_integral_final)); % cm^(-1)
    % Área da expressão anterior - integral da expressão anterior em ordem ao tempo.
    Area_Ca_norm_denom = Simpson_Method(intervalo_t,Ca_norm_num,[],'3/8'); % Simpson's 3/8 rule

    
end


% Representação da concentração calculada, mas em termos de absorvância, com mais detalhe (uso de mais pontos em t) 
tempo_detalhado = min(Dados_Experimentais.tempo):1:max(Dados_Experimentais.tempo);
Ca_norm_num_alldata = 1./(4*pi*K_optim.*tempo_detalhado).^(1/2).*exp(-(L_col-u0_optim.*tempo_detalhado).^2./(4.*K_optim.*tempo_detalhado)); % g/cm3 -> pela expressão (Eq. 25) dada no artigo Kong et al (J. Chromatogr. A, 1035, 177)  

Absorvance_calc = Ca_norm_num_alldata.*Area_Absorbance./Area_Ca_norm_denom; % (AU)



figure(3)
plot(Dados_Experimentais.tempo, Dados_Experimentais.absorbance, '*k', tempo_detalhado, Absorvance_calc, '-r')
xlabel('{\itt} (s)')
ylabel('Absorbance (AU)')
Title([Dados_Experimentais.soluto, ' in ', Dados_Experimentais.solvente, ' - ', ...
    num2str(Dados_Experimentais.Texp), ' K / ', num2str(Dados_Experimentais.Pexp/10), ' MPa (\lambda = ',...
    num2str(Dados_Experimentais.wavelength),' nm)'])
legend1 = legend('Experimental', ['Fitting: {\itD}_{12} = ', num2str(D12_optim*10^4, 5), '\times10^{-4} (cm^2/s)', ', \epsilon = ', num2str(fval*100,4), '%']);

set(legend1,'EdgeColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1]);
hold on


%%


Fim = '';



end 



%% Função de Optimização dos parâmetros
function Fobj = D12_optim_CIR(params, S10_1_time, S10_2_time, Metodo_Integral)

global R_col L_col Dados_Experimentais m_real_inject perc time_max_absorbance

U = params(1);
% U = L_col/time_max_absorbance;
a = params(2);


% Dois tempos que definem o intervalo de integraçao
t1 = S10_1_time; % tempo para o qual a altura do pico e 10% do valor máximo
t2 = S10_2_time; % tempo para o qual a altura do pico e 10% do valor máximo (NOTA: t1<t2)


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
    
    % Pontos experimenais qe entram no cálculo do integral entre os valores t1 e t2
    absorbance_integral_final = absorbance_integral;
    time_integral_final = time_integral;
    
       
elseif aviso1 == 1 &&  aviso2 ~= 1 %se n existe o t2 experimntal, acrescenta-se esse ponto
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2
    absorbance_integral_final = [absorbance_integral perc*max(Dados_Experimentais.absorbance)];
    time_integral_final = [time_integral t2];
        
elseif aviso1 ~= 1 &&  aviso2 == 1 %se n existe o t1 experimntal, acrescenta-se esse ponto
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2   
    absorbance_integral_final = [perc*max(Dados_Experimentais.absorbance) absorbance_integral];
    time_integral_final = [t1 time_integral];
        
elseif aviso1 ~= 1 &&  aviso2 ~= 1 %se n existem o t1 e o t2 experimntais, acrescentam-se esses pontos
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2
   absorbance_integral_final = [perc*max(Dados_Experimentais.absorbance) absorbance_integral perc*max(Dados_Experimentais.absorbance)];
   time_integral_final = [t1 time_integral t2];
   
end

%so pontos experimentais (sem os t1 e t2 "inventados" qd n existem):
absorbance_exp = absorbance_integral;
time_exp = time_integral;

% Curva de absorvância experimental (normalizada)
if Metodo_Integral == 1
    Area_Absorbance = trapz(time_integral_final,absorbance_integral_final); % (AU*s)
elseif Metodo_Integral == 2
    Area_Absorbance = Simpson_Method(time_integral_final,absorbance_integral_final,[],'1/3');
elseif Metodo_Integral == 3
    Area_Absorbance = Simpson_Method(time_integral_final,absorbance_integral_final,[],'3/8');
end


 Absor_norm = absorbance_exp./Area_Absorbance; % [s^(-1)]   absorvancia (pts so acima de 10%) normalizada pela area da curva acima de 10% da altura do pico


% Curva de absorvância "experimental" (normalizada) COMPLETA  
intervalo_t = linspace(min(time_integral_final), max(time_integral_final), 3000);
% pp = spline(time_integral_final,absorbance_integral_final); %spline cubica
% Interp_abs_exp_det = ppval(pp, intervalo_t); %Absorvancia interpolada para os 3000 pontos igualmente espaçados entre t1 e t2



% Normalização da curva calculada [(normalização segundo a expressão (Eq. 24) que é fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
Ca_norm_num = 1./(4*pi*a.*time_exp).^(1/2).*exp(-(L_col-U.*time_exp).^2./(4*a.*time_exp)); % cm^(-1)
% Área da expressão anterior - integral da expressão anterior em ordem ao tempo.


if Metodo_Integral == 1
    
%     Normalização da curva calculada [(normalização segundo a expressão (Eq. 24) que é fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
    Ca_norm_num_trap = 1./sqrt(4*pi*a.*intervalo_t).*exp(-(L_col-U.*intervalo_t).^2./(4*a.*intervalo_t)); % cm^(-1)
%    e a area pelo metodo dos trapézios (Area_Ca_norm_denom)??
    Area_Ca_norm_denom = trapz(intervalo_t,Ca_norm_num_trap);
    
elseif Metodo_Integral == 2
    
    Ca_norm_num_S13 = 1./(4*pi*a.*intervalo_t).^(1/2).*exp(-(L_col-U.*intervalo_t).^2./(4*a.*intervalo_t)); % cm^(-1)
    
    Area_Ca_norm_denom = Simpson_Method(intervalo_t,Ca_norm_num_S13,[],'1/3'); % Simpson's 3/8 rule

elseif Metodo_Integral == 3
        
    Ca_norm_num_S38 = 1./(4*pi*a.*intervalo_t).^(1/2).*exp(-(L_col-U.*intervalo_t).^2./(4*a.*intervalo_t)); % cm^(-1)
    
    Area_Ca_norm_denom = Simpson_Method(intervalo_t,Ca_norm_num_S38,[],'3/8'); % Simpson's 3/8 rule
    
end
                                                    
% Normalização da expressão da concentração (Área normalizada)
    Ca_norm = Ca_norm_num/Area_Ca_norm_denom; % [s^(-1)] (Ca_numerador so para os mm t_exp entre 1 e t2; area é q foi calculda com 3000 pontos entre t1 e t2)

figure(2)
plot(time_exp, Absor_norm, '*k', time_exp, Ca_norm, '-k')
xlabel('{\itt} (s)')
ylabel('{\itC}_{2,exp,norm}, {\itC}_{2,calc,norm} (s^{-1})')


% Definição da função objectivo
Termo1 = (Ca_norm-Absor_norm).^2;
Termo2 = Absor_norm.^2;


Fobj1 = trapz(time_exp,Termo1); %integral ou area da funçao "termo 1" entre t1 e t2 (ou o ponto + proximo)

Fobj2 = trapz(time_exp, Termo2); %integral ou area da funçao "termo 2" entre t1 e t2 (ou o ponto + proximo)

Fobj = (Fobj1/Fobj2)^(1/2);


ana = '';

end
