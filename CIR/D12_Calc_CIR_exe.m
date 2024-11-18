function D12_Calc_CIR_exe()

clear all 
clc
close all
format short

% Diretorio onde serao lidos os resultados experimentais (Diretorio e nome do ficheiro)
FileName = 'D:\OneDrive - Universidade de Aveiro\DAD\Testes_old_style\CIR\Testes_EtOH_Acetone_210618\210618_1500\Testes_210618_1500_200_270nm.xlsx';
FileNumber = 3;  % Número do ficheiro a ler dentro do mesmo livro de Excel
altura_pico_fitting = 0.1;
altura_pico_moment_method = 0.001; %standart 0.001

[D12_calc_fitting D12Calc_variance time_retention] = Exp_D12_col_reves_03_05_2018_fitting(FileName, FileNumber, altura_pico_fitting)

% [Resultados_CIR_moment] = Exp_D12_col_reves_03_05_2018_moment_method(FileName,FileNumber, altura_pico_moment_method)

BZ ='';

