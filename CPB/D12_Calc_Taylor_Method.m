function D12_Calc_Taylor_Method()

clear all 
clc
format short
close all
%D:\OneDrive\OneDrive - Universidade de Aveiro\DAD\Testes_old_style\Testes_eucalptol_210503\210511_1630\Resultados Matlab
% % Diretorio onde serao lidos os resultados experimentais (Diretorio e nome do ficheiro)
% FileName = 'C:\Users\bruno\Desktop\DAD\Testes_old_style\Testes_eucalptol_210423\teste2_eucalyptol_210423_1300_analise (250nm).xlsx';
% FileNumber = 3  % Número do ficheiro a ler dentro do mesmo livro de Excel
% sistema 
soluto = 'Gallic acid';
solvente = 'SCCO2_EtOH';
solvente_formula = 'CC';

Diretorio_gelal = 'Z:\OneDrive - Universidade de Aveiro\DAD\Testes_old_style\CPB\210715_Eucalyptol_CO2_etoh\210705_1200\Resultados\CH 220 nm\';
Diretorio_resultados = '';
lista = 'list_of_files.txt';
peak_no = 3;
altura_pico = 0.10;

poly_order_10_1 = 3; % should be 3, decrease only if fitting is poor; aplies to fitting method
poly_order_10_2 = 3; % should be 3, decrease only if fitting is poor; aplies to fitting method
poly_order_607_1 = 3; % should be 3, decrease only if fitting is poor; aplies to graphical method
poly_order_607_2 = 3; % should be 3, decrease only if fitting is poor; aplies to graphical method

%% ler nome dos ficheiros

Filelist = strcat(Diretorio_gelal, Diretorio_resultados, lista);

fid = fopen(Filelist);
tline = 1;
n = 1; 
while tline~=-1
    tline = fgets(fid);
    if ~ischar(tline)
        break
    end
    if n == peak_no
       File = tline;
       n = n + 1;
    else
        n = n + 1;
    end
end

FileName = strcat(Diretorio_gelal, Diretorio_resultados, File);

%% calcular D12      

[D12Calc_fitting, D12Calc_variance, time_max_absorbance, S10_1_time, S10_2_time, tres] = Exp_D12_col_n_reves(FileName,soluto, solvente, solvente_formula, altura_pico, poly_order_10_1, poly_order_10_2, poly_order_607_1, poly_order_607_2)

BrunoZezere_AnaMagalhaes ='';

end 