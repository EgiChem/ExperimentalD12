function [Resultados_CIR_moment] = Exp_D12_col_reves_03_05_2018_moment_method(FileName,FileNumber, altura)


global R_col L_col Dados_Experimentais h_pico

%% Programa de cálculo das difusividades a partir dos pontos experimentais (por fitting e por análise da largura do pico ou "método dos momentos")

% NOTAS IMPORTANTES : 
%
%                     Escolha do método de cálculo das áreas (Metodo_Integral): 
%                       - 1 para o método dos trapézios,
%                       - 2 para o método de Simpson 1/3,
%                       - 3 para o método de Simpson 3/8.
%                ------------------------------------------------
% 
%                     Escolha da forma a partir do qual se obtém a velocidade linear dentro da coluna de difusão (Veloc_Calc):
%                       - 11 para a bomba;
%                       - 22 para o medidor de bolha de sabão.

Metodo_Integral = 1; 
Veloc_Calc = 11;


%% Dados da coluna

d_enrolamento = 30.0; % Diâmetro do enrolamento (cm)
diam_col = (0.526)/10; % Diâmetro interno da coluna (cm) 
L_col = (16.887-0.028)*100; % Comprimento da coluna (cm)
espess_revestimento = 1E-4; % (cm) 
% R_col = (diam_col)/2 - espess_revestimento; %(cm) Raio da coluna
R_col = (diam_col)/2;%(cm) Raio da coluna


%% Informação acerca do soluto (massa, ou volume, ou concentração do soluto injectado
vol_inject = 0.1*1E-3; % mL

%% Leitura dos dados de um ficheiro de dados em txt (tempo vs Absorvância) 

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
Dados_Experimentais.absorbance = data_matrix(1:end,2); % Unidades de Absorvância (AU)

h_pico=altura;

% % smoothing signal
% Dados_Experimentais.absorbance= sgolayfilt(Dados_Experimentais.absorbance,2,9);


%% Figura 1 - Absorvância vs tempo
figure(2)
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
    if Dados_Experimentais.absorbance(j1) == h_pico*max_absorbance 
        S10_1_absorbance = Dados_Experimentais.absorbance(j1);
        S10_1_time = Dados_Experimentais.tempo(j1);
        n1 = n1 + 1;
    end 
end
% ---
npontos_int = 3; %pontos usados para cado lado na interpolação

if n1 ==0 %se não existir nenhum ponto correspondente a 0.1*maxabs
    dif_absorb1 = abs(Dados_Experimentais.absorbance(1:posicao_max_absorbance)- h_pico*max_absorbance); %diferença entre cada ponto e o maxabs
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
    
    if absorbance_1 < h_pico*max_absorbance % se o ponto + proximo esta abaixo
    
            absorbance_interpolar_1 = Dados_Experimentais.absorbance(posicao_min_1-(npontos_int-1):posicao_min_1+npontos_int); %pontos a usar na função interpolação: 2 para trás e 3 para a frente
            tempo_intepolar_1 = Dados_Experimentais.tempo(posicao_min_1-(npontos_int-1):posicao_min_1+npontos_int); %igual para o y (tempo)
    
    elseif absorbance_1 > h_pico*max_absorbance % se o ponto + proximo esta acima
        
            absorbance_interpolar_1 = Dados_Experimentais.absorbance(posicao_min_1-npontos_int:posicao_min_1+(npontos_int-1)); %pontos a usar na funçao de interpolaçao: 3 para trás e 2 para a frente
            tempo_intepolar_1 = Dados_Experimentais.tempo(posicao_min_1-npontos_int:posicao_min_1+(npontos_int-1));%igual para o y (tempo)
       
    end
    
    p_1 = polyfit(absorbance_interpolar_1,tempo_intepolar_1,2); %função de interpolação: polinómio de 2º grau
    S10_1_time = polyval(p_1,h_pico*max_absorbance); %tempo correspondente a 0.1*maxabs
    
    tempo_I_ver = polyval(p_1,linspace(min(absorbance_interpolar_1),max(absorbance_interpolar_1), 500)); %para ver no grafico
    
    hold on
    figure(2)
    subplot(1,3,1)
    plot(tempo_intepolar_1, absorbance_interpolar_1, '*r', tempo_I_ver, linspace(min(absorbance_interpolar_1),max(absorbance_interpolar_1), 500), '-r', S10_1_time, h_pico*max_absorbance, 'ob') %representação do 0.1*maxabs, do ponto exp + proximo e da f. interpolação
    
    
end

% --- tudo igual mas para a outra metade da curva

% Determinaçao da existência de um ponto experimenal a 10 % da altura do pico (2ª Parte)
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
    figure(2)
    subplot(1,3,1)
    plot(tempo_intepolar_2, absorbance_interpolar_2, '*r', tempo_II_ver, linspace(min(absorbance_interpolar_2),max(absorbance_interpolar_2), 500),'-r', S10_2_time, h_pico*max_absorbance, 'ob')
    xlabel('{\itt} (s)')
    ylabel('Absorbance (AU)')
end    
% ---    

% Expressão para o cálculo da simetria do pico
S10 = (S10_2_time-time_max_absorbance)/(time_max_absorbance-S10_1_time);


if Veloc_Calc == 11
    
 % Determinação da velocidade linear do solvente no interior da coluna (neste caso e específico para CO2)
    CO2_dens_bomba = densityCO2(Dados_Experimentais.Pexp_bomba/1.01325, Dados_Experimentais.Texp_bomba); % (g/cm3)
    CO2_dens_estufa = densityCO2(Dados_Experimentais.Pexp/1.01325, Dados_Experimentais.Texp); % (g/cm3)
    
    caudal_mass_CO2 = Dados_Experimentais.caudal_bomba*CO2_dens_bomba; % (g/mL)
     
    Q_estufa = caudal_mass_CO2/CO2_dens_estufa; % (mL/min)
    
    Area_SR_coluna = pi*R_col^2; %(cm2)
    
    u0_exp = Q_estufa/60/Area_SR_coluna; % (cm/s)
    
elseif Veloc_Calc == 22
    
  % Determinação da velocidade linear do solvente no interior da coluna (neste caso e específico para CO2)
    CO2_dens_bolha = densityCO2(Dados_Experimentais.Pamb/1.01325, Dados_Experimentais.Tbolha); % (g/cm3)
    CO2_dens_estufa = densityCO2(Dados_Experimentais.Pexp/1.01325, Dados_Experimentais.Texp); % (g/cm3)
    
    caudal_mass_CO2 = Dados_Experimentais.caudal_bolha*CO2_dens_bolha; % (g/mL)
    
    Q_estufa = caudal_mass_CO2/CO2_dens_estufa; % (mL/min)
    
    Area_SR_coluna = pi*R_col^2; % (cm2)
   
    u0_exp = Q_estufa/60/Area_SR_coluna; % (cm/s)

elseif Veloc_Calc ~= 11 || Veloc_Calc ~= 22

return

end

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




%%
if aviso1 == 1 && aviso2 == 1
    
    % Pontos experimenais qe entram no cálculo do integral entre os valores t1 e t2
    absorbance_integral_final = absorbance_integral;
    time_integral_final = time_integral;
    
       
elseif aviso1 == 1 &&  aviso2 ~= 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2
    absorbance_integral_final = [absorbance_integral h_pico*max(Dados_Experimentais.absorbance)];
    time_integral_final = [time_integral t2];
        
elseif aviso1 ~= 1 &&  aviso2 == 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2   
    absorbance_integral_final = [h_pico*max(Dados_Experimentais.absorbance) absorbance_integral];
    time_integral_final = [t1 time_integral];
        
elseif aviso1 ~= 1 &&  aviso2 ~= 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2
   absorbance_integral_final = [h_pico*max(Dados_Experimentais.absorbance) absorbance_integral h_pico*max(Dados_Experimentais.absorbance)];
   time_integral_final = [t1 time_integral t2];
   
end 
    

if Metodo_Integral == 1
    % Curva de absorvância experimental (normalizada) COMPLETA 
    Area_Absorbance = trapz(time_integral_final,absorbance_integral_final); % (AU*s) -> Método dos trapézios
    Absor_norm = Dados_Experimentais.absorbance./Area_Absorbance; % [s^(-1)]   
    
   
elseif Metodo_Integral == 2 % Simpson's 1/3 rule
    % Curva de absorvância experimental (normalizada) COMPLETA  
    intervalo_t = linspace(min(time_integral_final), max(time_integral_final), 3000);
    pp = spline(time_integral_final,absorbance_integral_final);
    Interp_abs_exp_det = ppval(pp, intervalo_t);
    
    Area_Absorbance = Simpson_Method(intervalo_t,Interp_abs_exp_det,[],'1/3'); % Simpson's 3/8 rule
    Absor_norm = Dados_Experimentais.absorbance./Area_Absorbance; % [s^(-1)]
            

elseif Metodo_Integral == 3  % Simpson's 3/8 rule
    intervalo_t = linspace(min(time_integral_final), max(time_integral_final), 3000);
    pp = spline(time_integral_final,absorbance_integral_final);
    Interp_abs_exp_det = ppval(pp, intervalo_t);
    
    Area_Absorbance = Simpson_Method(intervalo_t,Interp_abs_exp_det,[],'3/8'); % Simpson's 3/8 rule
    Absor_norm = Dados_Experimentais.absorbance./Area_Absorbance; % [s^(-1)] 
  
    
end


% Representação da concentração calculada, mas em termos de absorvância, com mais detalhe (uso de mais pontos em t) 
tempo_detalhado = min(Dados_Experimentais.tempo):1:max(Dados_Experimentais.tempo);


%% Cálculo (Metodo dos momentos)

Absor_norm_exp = absorbance_integral_final./Area_Absorbance; % [s^(-1)]   absorvancia (pts so acima de 10%) normalizada pela area da curva acima de 10% da altura do pico

Absor_norm_exp_i=Absor_norm_exp-Absor_norm_exp(1);

figure(2)
subplot(1,3,2)
plot(time_integral_final,Absor_norm_exp)

%momento 0
momento_0 = trapz(time_integral_final,Absor_norm_exp);

%momento 1
t_abs = time_integral_final.* Absor_norm_exp;
momento_1_numerador = trapz(time_integral_final, t_abs);
momento_1 = momento_1_numerador./momento_0; %tmed

%momento 2
t2_abs = (time_integral_final-momento_1).^2.* Absor_norm_exp;
momento_2_numerador = trapz(time_integral_final,t2_abs );
momento2 = momento_2_numerador./momento_0; %sig^2

% calculo k

% expressão (Eq. 38) que é fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)
alfa=momento2./momento_1.^2;
k_moment_method=2*(2-alfa)./(3+(1+4*alfa).^0.5)  .*  u0_exp*momento_1./L_col -1;


%calculo D12

%metodo 1
% expressão (Eq. 39) que é fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)

beta=(2*alfa - 1 + (1+4*alfa)^0.5) / (4*(2-alfa));
gama=(1 + 6*k_moment_method + 11*k_moment_method^2) * R_col^2 / (48*L_col^2 *(1+k_moment_method)^2);

D12_moment_method_1 = 2*gama/ (beta+(beta.^2-4*gama)^0.5) *L_col * u0_exp;


%metodo 2
% através do 1º momento (Eq. 7d, 28, 36)

t_med=momento_1;
U_mm=u0_exp/(1+k_moment_method);
a_mm=L_col*U_mm/2 * (U_mm/L_col*t_med-1); 
% a_mm=D12_moment_method_1/(1+k_moment_method) + (1+6*k_moment_method+11*k_moment_method^2 )/(1+k_moment_method) * (R_col^2*U_mm^2)/(48*D12_moment_method_1)

% variavel_1_mm = (1 + 6*k_moment_method + 11*k_moment_method^2)/(1 + k_moment_method);
% variavel_2_mm = (1 + 6*k_moment_method + 11*k_moment_method^2)/(1 + k_moment_method)^2;
% 
% termo_D12_mm_1 = variavel_1_mm*(R_col^2*U_mm^2)/(24*a_mm);
% termo_D12_mm_2 = 1 + (1 - variavel_2_mm*R_col^2*U_mm^2/(12*a_mm^2))^(1/2);
% 
% D12_moment_method_2 = termo_D12_mm_1/termo_D12_mm_2; % Valor de Difusividade (cm2/s


% t=L_col/U_mm * (1+2*a_mm/(L_col*U_mm))
% sig=2*(L_col/U_mm)^2 * a_mm/(L_col*U_mm) * (1+4*a_mm/(L_col*U_mm))
% Calculo da obs normalizada

%trapézio
    % Normalização da curva calculada COMPLETA [(normalização segundo a expressão (Eq. 24) que é fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
    Ca_norm_num_mm = 1./(4*pi*a_mm.*time_integral_final).^(1/2).*exp(-(L_col-U_mm.*time_integral_final).^2./(4.*a_mm.*time_integral_final)); % cm^(-1)
    % Área da expressão anterior - integral da expressão anterior em ordem ao tempo.
    Area_Ca_norm_denom_mm = trapz(time_integral_final,Ca_norm_num_mm); % (s.cm^(-1)) -> Método dos trapézios
    
Ca_norm_num_alldata_mm = 1./(4*pi*a_mm.*tempo_detalhado).^(1/2).*exp(-(L_col-U_mm.*tempo_detalhado).^2./(4.*a_mm.*tempo_detalhado)); % g/cm3 -> pela expressão (Eq. 25) dada no artigo Kong et al (J. Chromatogr. A, 1035, 177)
Absorvance_calc_mm = Ca_norm_num_alldata_mm.*Area_Absorbance./Area_Ca_norm_denom_mm; % (AU)


% Compilação da Informação metodo dos momentos

Resultados_CIR_moment.Area=Area_Absorbance ;
Resultados_CIR_moment.a=a_mm ;
Resultados_CIR_moment.U=U_mm ;
Resultados_CIR_moment.D12_1=D12_moment_method_1 ;
% Resultados_CIR_moment.D12_2=D12_moment_method_2;
Resultados_CIR_moment.k=k_moment_method ;


% gráfico

figure(2)
subplot(1,3,3)
plot(Dados_Experimentais.tempo./60, Dados_Experimentais.absorbance, 'ok', tempo_detalhado/60, Absorvance_calc_mm, '-g')
xlabel('{\itt} (s)')
ylabel('Absorbance (AU)')
%% ADICIONAR
title([Dados_Experimentais.soluto, ' in ', Dados_Experimentais.solvente, ' - ', ...
    num2str(Dados_Experimentais.Texp), ' K / ', num2str(Dados_Experimentais.Pexp/10), ' MPa (\lambda = ',...
    num2str(Dados_Experimentais.wavelength),' nm)'])
%% ------------------------------------------------
legend1 = legend('Experimental', ['Moment Method: {\itD}_{12} = ', num2str(D12_moment_method_1*10^4, 5), '\times10^{-4} (cm^2/s)']);

set(legend1,'EdgeColor',[1 1 1]);


Resultados_CIR_moment;
Fim = '';



end 
