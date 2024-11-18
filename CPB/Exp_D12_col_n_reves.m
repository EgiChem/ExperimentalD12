function [D12Calc_fitting, D12Calc_variance, time_max_absorbance, S10_1_time, S10_2_time, tres] = Exp_D12_col_n_reves(FileName,soluto, solvente, solvente_formula, altura, poly_order_10_1, poly_order_10_2, poly_order_607_1, poly_order_607_2)

global R_col L_col Dados_Experimentais altura tres

%% Dados da coluna
d_enrolamento = 30; % Diâmetro do enrolamento (cm)
L_col = 1118.2-13.140-3.48; % Comprimento da coluna (cm)
diam_col = 0.522/10; % Diâmetro interno da coluna (cm) 

R_col = diam_col/2;%*1.027; %(cm) Raio da coluna

%% Informação acerca do soluto (massa, ou volume, ou concentração do soluto injectado

vol_inject = 0.1*1E-3; % mL

%% Leitura dos dados de um ficheiro de dados em txt (tempo vs Absorvância) 
% [data_matrix, text_data, alldata] = xlsread(FileName,FileNumber);

[system, date, Inject_time, Wavelenght, T, P, x2mol, x2vol, Q, Tbath, data] = readpeak(FileName);

%% Dados Experimentais 
Dados_Experimentais.soluto = soluto;
Dados_Experimentais.solvente = solvente;
Dados_Experimentais.solvente_formula = solvente_formula;
%%
Dados_Experimentais.wavelength = Wavelenght; % (nm)
Dados_Experimentais.Pexp = P; %(bar)
Dados_Experimentais.Texp = T; % (K)
% Dados_Experimentais.Texp_bomba = data_matrix(1,9)+ 273.15; % (K)
Dados_Experimentais.caudal_bomba = Q; % (mL/min)
Dados_Experimentais.tempo = data(1:end,1); % tempo (s)
Dados_Experimentais.absorbance = data(1:end,3); % Unidades de Absorvância (uAU)



% Figura 1 - Absorvância vs tempo
figure(1)
subplot(2,2,1); plot(Dados_Experimentais.tempo, Dados_Experimentais.absorbance, '*k')
xlabel('{\itt} (s)')
ylabel('Absorbance (AU)')
title([Dados_Experimentais.soluto, ' in ', Dados_Experimentais.solvente_formula, ' - ', ...
    num2str(Dados_Experimentais.Texp), ' K / ', num2str(Dados_Experimentais.Pexp), ' bar (\lambda = ',...
    num2str(Dados_Experimentais.wavelength),' nm)'])
axis([min(Dados_Experimentais.tempo)-10 max(Dados_Experimentais.tempo)+10  -0.001 max(Dados_Experimentais.absorbance)+0.01])

hold on

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
n1 = 0;
for j1 = 1:posicao_max_absorbance
    if Dados_Experimentais.absorbance(j1) == altura*max_absorbance 
        S10_1_absorbance = Dados_Experimentais.absorbance(j1);
        S10_1_time = Dados_Experimentais.tempo(j1);
        n1 = n1 + 1;
    end 
end


if n1 ==0
    dif_absorb1 = abs(Dados_Experimentais.absorbance(1:posicao_max_absorbance)- altura*max_absorbance);
    min_dif_absorb1 = min(dif_absorb1);
    
    
    for j1 = 1:posicao_max_absorbance
        if min_dif_absorb1 == dif_absorb1(j1)
            absorbance_1 = Dados_Experimentais.absorbance(j1);
            tempo_1 = Dados_Experimentais.tempo(j1);
            posicao_min_1 = j1;
        end
    end
    
    if absorbance_1 < altura*max_absorbance
    
            absorbance_interpolar_1 = Dados_Experimentais.absorbance(posicao_min_1-2:posicao_min_1+3);
            tempo_intepolar_1 = Dados_Experimentais.tempo(posicao_min_1-2:posicao_min_1+3);
    
    elseif absorbance_1 > altura*max_absorbance
        
            absorbance_interpolar_1 = Dados_Experimentais.absorbance(posicao_min_1-3:posicao_min_1+2);
            tempo_intepolar_1 = Dados_Experimentais.tempo(posicao_min_1-3:posicao_min_1+2);
       
    end
    
    p_1 = polyfit(absorbance_interpolar_1,tempo_intepolar_1,poly_order_10_1);
    S10_1_time = polyval(p_1,altura*max_absorbance);
    
    tempo_I_ver = polyval(p_1,linspace(min(absorbance_interpolar_1),max(absorbance_interpolar_1), 500));
    
    hold on
%     figure(1)
    subplot(2,2,1); plot(tempo_intepolar_1, absorbance_interpolar_1, '*r', tempo_I_ver, linspace(min(absorbance_interpolar_1),max(absorbance_interpolar_1), 500), '-r', S10_1_time, altura*max_absorbance, 'ob')
    
    
end


% Determinaçao da existência de um ponto experimenal a 10 % da altura do pico (2ª Parte)
n2 = 0;
for j2 = posicao_max_absorbance+1:length(Dados_Experimentais.tempo)
    if Dados_Experimentais.absorbance(j2) == altura*max_absorbance 
        S10_2_absorbance = Dados_Experimentais.absorbance(j2);
        S10_2_time = Dados_Experimentais.tempo(j2);
        n2 = n2 + 1;
    end 
end


if n2 ==0
    dif_absorb2 = abs(Dados_Experimentais.absorbance(posicao_max_absorbance+1:length(Dados_Experimentais.tempo))- altura*max_absorbance);
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
        
    if absorbance_2 < altura*max_absorbance
    
            absorbance_interpolar_2 = Dados_Experimentais.absorbance(posicao_min_2-3:posicao_min_2+2);
            tempo_intepolar_2 = Dados_Experimentais.tempo(posicao_min_2-3:posicao_min_2+2);
    
    elseif absorbance_2 > altura*max_absorbance
        
            absorbance_interpolar_2 = Dados_Experimentais.absorbance(posicao_min_2-2:posicao_min_2+3);
            tempo_intepolar_2 = Dados_Experimentais.tempo(posicao_min_2-2:posicao_min_2+3);
       
    end    
    
    p_2 = polyfit(absorbance_interpolar_2,tempo_intepolar_2,poly_order_10_2);
    S10_2_time = polyval(p_2,altura*max_absorbance);
    
    tempo_II_ver = polyval(p_2,linspace(min(absorbance_interpolar_2),max(absorbance_interpolar_2), 500));
    
    hold on
%     figure(1)
    subplot(2,2,1); plot(tempo_intepolar_2, absorbance_interpolar_2, '*r', tempo_II_ver, linspace(min(absorbance_interpolar_2),max(absorbance_interpolar_2), 500),'-r', S10_2_time, altura*max_absorbance, 'ob')
    xlabel('{\itt} (s)')
    ylabel('Absorbance (AU)')
end    
    

% Expressão para o cálculo da simetria do pico
S10 = (S10_2_time-time_max_absorbance)/(time_max_absorbance-S10_1_time);



%% Optimizaçao
options = optimset;
% Modify options setting
options = optimset(options,'Display' ,'notify');
options = optimset(options,'MaxIter' ,50000);
options = optimset(options,'MaxFunEvals' ,100000);
options = optimset(options,'TolX' ,1e-8);
options = optimset(options,'TolFun',1e-10);
options = optimset(options,'Diagnostics' ,'on');

% Parâmetros iniciais (u0 e K)
params_ini = [1  0.1]; %BZ

[params_optim, fval,exitflag,output] = fminsearch(@(params) D12_optim_Taylor(params, S10_1_time, S10_2_time), params_ini, options);

u0_optim = params_optim(1);
K_optim = params_optim(2);


% Dois tempos que definem o intervalo de integraçao
t1 = S10_1_time; % tempo para o qual a altura do pico e 10% do valor máximo
t2 = S10_2_time; % tempo para o qual a altura do pico e 10% do valor máximo (NOTA: t1<t2)


n = 0;

for i = 1:length(Dados_Experimentais.tempo)
    if Dados_Experimentais.tempo(i)>= t1 & Dados_Experimentais.tempo(i)<= t2 
    n = n+1;
    time_integral(n) = Dados_Experimentais.tempo(i);
    absorbance_integral(n) = Dados_Experimentais.absorbance(i);
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



% D12 Calculado pelo método de ajuste da curva 
D12_optim1 = (K_optim -sqrt(K_optim^2 - R_col^2*u0_optim^2/12))/2;

% D12_optim2 = K_optim/2*(1-(1-u0_optim^2*R_col^2/(12*K_optim^2))^(1/2));
% D12_optim3 = (R_col^2*u0_optim^2)/(24*K_optim)/(1 + (1-R_col^2*u0_optim^2/(12*K_optim^2))^(1/2));



if aviso1 == 1 && aviso2 == 1
    
    % Pontos experimenais qe entram no cálculo do integral entre os valores t1 e t2
    absorbance_integral_final = absorbance_integral;
    time_integral_final = time_integral;
    
       
elseif aviso1 == 1 &&  aviso2 ~= 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2
    absorbance_integral_final = [absorbance_integral altura*max(Dados_Experimentais.absorbance)];
    time_integral_final = [time_integral t2];
        
elseif aviso1 ~= 1 &&  aviso2 == 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2   
    absorbance_integral_final = [altura*max(Dados_Experimentais.absorbance) absorbance_integral];
    time_integral_final = [t1 time_integral];
        
elseif aviso1 ~= 1 &&  aviso2 ~= 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2
   absorbance_integral_final = [altura*max(Dados_Experimentais.absorbance) absorbance_integral altura*max(Dados_Experimentais.absorbance)];
   time_integral_final = [t1 time_integral t2];
   
end 
    
    
% Curva de absorvância experimental (normalizada) COMPLETA 
Area_Absorbance = trapz(time_integral_final,absorbance_integral_final); % (AU*s)
Absor_norm = Dados_Experimentais.absorbance./Area_Absorbance; % [s^(-1)]   


% Normalização da curva calculada COMPLETA [(normalização segundo a expressão (Eq. 24) que é fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
Ca_norm_num = 1./(4*pi*K_optim.*time_integral_final).^(1/2).*exp(-(L_col-u0_optim.*time_integral_final).^2./(4.*K_optim.*time_integral_final)); % cm^(-1)
% Área da expressão anterior - integral da expressão anterior em ordem ao tempo.
Area_Ca_norm_denom = trapz(time_integral_final,Ca_norm_num); % (s.cm^(-1))


% Representação da concentração calculada, mas em termos de absorvância, com mais detalhe (uso de mais pontos em t) 
tempo_detalhado = min(Dados_Experimentais.tempo):1:max(Dados_Experimentais.tempo);
Ca_norm_num_alldata = 1./(4*pi*K_optim.*tempo_detalhado).^(1/2).*exp(-(L_col-u0_optim.*tempo_detalhado).^2./(4.*K_optim.*tempo_detalhado)); % g/cm3 -> pela expressão (Eq. 25) dada no artigo Kong et al (J. Chromatogr. A, 1035, 177)  

Absorvance_calc = Ca_norm_num_alldata.*Area_Absorbance./Area_Ca_norm_denom;

% figure(3)
subplot(2,2,3); plot(Dados_Experimentais.tempo, Dados_Experimentais.absorbance, '*k', tempo_detalhado, Absorvance_calc, '-r')
xlabel('{\itt} (s)')
ylabel('Absorbance (AU)')
title([Dados_Experimentais.soluto, ' in ', Dados_Experimentais.solvente_formula, ' - ', ...
    num2str(Dados_Experimentais.Texp), ' K / ', num2str(Dados_Experimentais.Pexp), ' bar (\lambda = ',...
    num2str(Dados_Experimentais.wavelength),' nm)'])
legend1 = legend('Experimental', ['Fitting: {\itD}_{12} = ', num2str(D12_optim1*10^4, 5), ' x 10^{-4} (cm^2/s)', ', \epsilon = ', num2str(fval*100,4), '%']);
set(legend1,'EdgeColor',[1 1 1]);
axis([min(Dados_Experimentais.tempo)-10 max(Dados_Experimentais.tempo)+10  -0.001 max(Dados_Experimentais.absorbance)+0.01])
hold on


% Compilação da Informação (Método por ajuste da curva - Fitting)
D12Calc_fitting.S10_10 = S10;
D12Calc_fitting.u0_optim = u0_optim; % (cm/s)
D12Calc_fitting.AreaPeak = Area_Absorbance;
D12Calc_fitting.maxAbs = max_absorbance;
D12Calc_fitting.D12_Calc = D12_optim1; % (cm2/s)

%% 
%% Análise da Largura do pico

% Determinaçao de 60.7% da altura do pico em ambas as partes (1ª parte) 
p1 = 0;
for q1 = 1:posicao_max_absorbance
    if Dados_Experimentais.absorbance(q1) == 0.607*max_absorbance
        W_time_607_1 = Dados_Experimentais.absorbance(q1);  
        p1 = p1 + 1;
   end
end

if p1 == 0
    
    difernca_607_1 = abs(Dados_Experimentais.absorbance(1:posicao_max_absorbance)- 0.607*max_absorbance);
    min_difernca_607_1 = min(difernca_607_1);

     for q1 = 1:posicao_max_absorbance
        if min_difernca_607_1 == difernca_607_1(q1)
            absorbance_607_1 = Dados_Experimentais.absorbance(q1);
            tempo_607_1 = Dados_Experimentais.tempo(q1);
            posicao_min_607_1 = q1;
        end
     end   
    
    if absorbance_607_1 < 0.607*max_absorbance
       absorbance_interpolar_607_1 = Dados_Experimentais.absorbance(posicao_min_607_1-2:posicao_min_607_1+3);
       tempo_intepolar_607_1 = Dados_Experimentais.tempo(posicao_min_607_1-2:posicao_min_607_1+3);
    
    elseif absorbance_607_1 > 0.607*max_absorbance
        absorbance_interpolar_607_1 = Dados_Experimentais.absorbance(posicao_min_607_1-3:posicao_min_607_1+2);
        tempo_intepolar_607_1 = Dados_Experimentais.tempo(posicao_min_607_1-3:posicao_min_607_1+2);
    end    
     
    p_607_1 = polyfit(absorbance_interpolar_607_1,tempo_intepolar_607_1,poly_order_607_1);
    time_607_1 = polyval(p_607_1,0.607*max_absorbance);

    tempo_607_ver_1 = polyval(p_607_1, linspace(min(absorbance_interpolar_607_1),max(absorbance_interpolar_607_1), 500));

%     figure(1)
    subplot(2,2,1); plot(time_607_1, 0.607*max_absorbance, 'sg',tempo_607_ver_1,linspace(min(absorbance_interpolar_607_1),max(absorbance_interpolar_607_1), 500), '-g')
    
end 

% Determinaçao de 60.7% da altura do pico em ambas as partes (2ª parte) 

p2 = 0;
for q2 = posicao_max_absorbance+1:length(Dados_Experimentais.tempo)
    if Dados_Experimentais.absorbance(q2) == 0.607*max_absorbance
        W_time_607_2 = Dados_Experimentais.absorbance(q2);  
        p2 = p2 + 1;
   end
end


if p2 == 0
    
    difernca_607_2 = abs(Dados_Experimentais.absorbance(posicao_max_absorbance+1:length(Dados_Experimentais.tempo))- 0.607*max_absorbance);
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
    

    if absorbance_607_2 < 0.607*max_absorbance
       absorbance_interpolar_607_2 = Dados_Experimentais.absorbance(posicao_min_607_2-3:posicao_min_607_2+2);
       tempo_intepolar_607_2 = Dados_Experimentais.tempo(posicao_min_607_2-3:posicao_min_607_2+2);
    
    elseif absorbance_607_2 > 0.607*max_absorbance
        absorbance_interpolar_607_2 = Dados_Experimentais.absorbance(posicao_min_607_2-2:posicao_min_607_2+3);
        tempo_intepolar_607_2 = Dados_Experimentais.tempo(posicao_min_607_2-2:posicao_min_607_2+3);
    end    
     
    p_607_2 = polyfit(absorbance_interpolar_607_2,tempo_intepolar_607_2,poly_order_607_2);
    time_607_2 = polyval(p_607_2,0.607*max_absorbance);
    

    tempo_607_ver_2 = polyval(p_607_2, linspace(min(absorbance_interpolar_607_2),max(absorbance_interpolar_607_2), 500));

%     figure(1)
    hold on
    
    subplot(2,2,1); plot(time_607_2, 0.607*max_absorbance, 'sg',tempo_607_ver_2,linspace(min(absorbance_interpolar_607_2),max(absorbance_interpolar_607_2), 500), '-g')
    
end 


W_time_607 = (time_607_2-time_607_1)/2;
S10_607 = (time_607_2-time_max_absorbance)/(time_max_absorbance-time_607_1);


% NOTA neste caso uso a velocidade optimizada, mas devo sempre usar a % velocidade experimental obtida
Analise_pico.u0_exp = L_col/time_max_absorbance;

Analise_pico.H = Analise_pico.u0_exp^2*W_time_607^2/L_col;
% Analise_pico.H = Analise_pico.u0_exp^2*W_time_607^2/L_col/5.545;

Analise_pico.D12_calc = Analise_pico.u0_exp/4*(Analise_pico.H-sqrt(Analise_pico.H^2-R_col^2/3));

Analise_pico.u_opt = sqrt(48)*Analise_pico.D12_calc/R_col;

% Compilação da Informação (Método análise da Largura do Pico)
D12Calc_variance.S10_607 = S10_607;
D12Calc_variance.H = Analise_pico.H;
D12Calc_variance.u0_exp = Analise_pico.u0_exp;
D12Calc_variance.D12_calc = Analise_pico.D12_calc;
D12Calc_variance.u_opt=Analise_pico.u_opt;



% Curva calculada (normalizada) 
Analise_pico.K = Analise_pico.D12_calc + Analise_pico.u0_exp^2*R_col^2/(48*Analise_pico.D12_calc);

% Normalização da curva calculada COMPLETA [(normalização segundo a expressão (Eq. 24) que é fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
Analise_pico.Ca_norm_num_all = 1./(4*pi*Analise_pico.K.*time_integral_final).^(1/2).*exp(-(L_col-Analise_pico.u0_exp.*time_integral_final).^2./(4.*Analise_pico.K.*time_integral_final)); % cm^(-1)
% Área da expressão anterior - integral da expressão anterior em ordem ao tempo.
Analise_pico.Area_Ca_norm_all = trapz(time_integral_final,Analise_pico.Ca_norm_num_all); % (s.cm^(-1))


Analise_pico.Ca_n_num_all_detal = 1./(4*pi*Analise_pico.K.*tempo_detalhado).^(1/2).*exp(-(L_col-Analise_pico.u0_exp.*tempo_detalhado).^2./(4.*Analise_pico.K.*tempo_detalhado)); % g/cm3 -> pela expressão (Eq. 25) dada no artigo Kong et al (J. Chromatogr. A, 1035, 177)  

Analise_pico.Absorvance_calc = Analise_pico.Ca_n_num_all_detal.*Area_Absorbance./Analise_pico.Area_Ca_norm_all;


% figure(4)
subplot(2,2,4); plot(Dados_Experimentais.tempo, Dados_Experimentais.absorbance, '*k')
title([Dados_Experimentais.soluto, ' in ', Dados_Experimentais.solvente_formula, ' - ', ...
    num2str(Dados_Experimentais.Texp), ' K / ', num2str(Dados_Experimentais.Pexp), ' bar (\lambda = ',...
    num2str(Dados_Experimentais.wavelength),' nm)'])

hold on
plot(tempo_detalhado, Analise_pico.Absorvance_calc, '-g')
xlabel('{\itt} (s)')
ylabel('Absorbance (AU)')
legend1 = legend('Experimental', ['Calc: {\itD}_{12} = ', num2str(Analise_pico.D12_calc*10^4, 5), '\times10^{-4} (cm^2/s)']);
% legend1 = legend('Experimental', 'Calculated');
set(legend1,'EdgeColor',[1 1 1]);
axis([min(Dados_Experimentais.tempo)-10 max(Dados_Experimentais.tempo)+10  -0.001 max(Dados_Experimentais.absorbance)+0.01])

Fim = '';



end 



%% Função de Optimização dos parâmetros
function Fobj = D12_optim_Taylor(params, S10_1_time, S10_2_time)

global R_col L_col Dados_Experimentais m_real_inject altura tres

u0 = params(1);
K = params(2);


% Dois tempos que definem o intervalo de integraçao
t1 = S10_1_time; % tempo para o qual a altura do pico e 10% do valor máximo
t2 = S10_2_time; % tempo para o qual a altura do pico e 10% do valor máximo (NOTA: t1<t2)


n = 0;

for i = 1:length(Dados_Experimentais.tempo)
    if Dados_Experimentais.tempo(i)>= t1 & Dados_Experimentais.tempo(i)<= t2 
    n = n+1;
    time_integral(n) = Dados_Experimentais.tempo(i);
    absorbance_integral(n) = Dados_Experimentais.absorbance(i);
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
% K = D12_optim + u0^2*R_col^2/(48*D12_optim);


if aviso1 == 1 && aviso2 == 1
    
    % Pontos experimenais qe entram no cálculo do integral entre os valores t1 e t2
    absorbance_integral_final = absorbance_integral;
    time_integral_final = time_integral;
    
       
elseif aviso1 == 1 &&  aviso2 ~= 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2
    absorbance_integral_final = [absorbance_integral altura*max(Dados_Experimentais.absorbance)];
    time_integral_final = [time_integral t2];
        
elseif aviso1 ~= 1 &&  aviso2 == 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2   
    absorbance_integral_final = [altura*max(Dados_Experimentais.absorbance) absorbance_integral];
    time_integral_final = [t1 time_integral];
        
elseif aviso1 ~= 1 &&  aviso2 ~= 1
    
    % Pontos necessários para calcular o integral entre os valores t1 e t2
   absorbance_integral_final = [altura*max(Dados_Experimentais.absorbance) absorbance_integral altura*max(Dados_Experimentais.absorbance)];
   time_integral_final = [t1 time_integral t2];
   
end 
    
    
% Curva de absorvância experimental (normalizada) 
Area_Absorbance = trapz(time_integral_final,absorbance_integral_final); % (AU*s)
Absor_norm = absorbance_integral_final./Area_Absorbance; % [s^(-1)]   


% Normalização da curva calculada [(normalização segundo a expressão (Eq. 24) que é fornecido no artigo Kong et al. (J. Chromatogr. A, 1035, 177)]  
Ca_norm_num = 1./(4*pi*K.*time_integral_final).^(1/2).*exp(-(L_col-u0.*time_integral_final).^2./(4.*K.*time_integral_final)); % cm^(-1)
% Área da expressão anterior - integral da expressão anterior em ordem ao tempo.
Area_Ca_norm_denom = trapz(time_integral_final,Ca_norm_num); % (s.cm^(-1))
% Normalização da expressão da concentração (Área normalizada)
Ca_norm = Ca_norm_num./Area_Ca_norm_denom; % [s^(-1)] 

time_exp = time_integral_final;
tet=time_exp.*Ca_norm;
tres=trapz(time_exp,tet);
% Normalização da expressão da concentração (Área normalizada)
Ca_norm = Ca_norm_num/Area_Ca_norm_denom; % [s^(-1)] (Ca_numerador so para os mm t_exp entre 1 e t2; area é q foi calculda com 3000 pontos entre t1 e t2)

% figure(2)
subplot(2,2,2); plot(time_exp, Absor_norm, '*k', time_exp, Ca_norm, '-k')
xlabel('{\itt} (s)')
ylabel('{\itC}_{2,exp,norm}, {\itC}_{2,calc,norm} (s^{-1})')
                                                       

% Definição da função objectivo
Termo1 = (Ca_norm-Absor_norm).^2;
Termo2 = Absor_norm.^2;


Fobj1 = trapz(time_integral_final,Termo1);

Fobj2 = trapz(time_integral_final, Termo2);

Fobj = (Fobj1/Fobj2)^(1/2);

end



