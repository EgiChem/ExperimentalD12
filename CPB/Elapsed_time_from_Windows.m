function [Elapsed_time_Windows Elapsed_time_Windows_real] = Elapsed_time_from_Windows(Elapsed_time, Windows_time_n_date) 

NDP = numel(Windows_time_n_date);
Elapsed_time_Windows_real = zeros(NDP-1,1);
Elapsed_time_Windows = zeros(NDP-1,1);
Elapsed_time_start = Elapsed_time(1);
Key1 = ' ';
Key2 = ':';
for n = 2:NDP
    scan = char(Windows_time_n_date(n));
    index_Key1 = strfind(scan, Key1);
    scan(1:index_Key1) = [];
    index_Key2 = strfind(scan, Key2);
    
    horas = str2num(scan(1:(index_Key2(1)-1)));
    minutos = str2num(scan((index_Key2(1)+1):index_Key2(2)-1));
    segundos = str2num(scan((index_Key2(2)+1):end));
    
    Elapsed_time_Windows_real(n-1) = horas*60*60 + minutos*60 + segundos;
    Elapsed_time_Windows(n-1) = Elapsed_time_Windows_real(n-1)-Elapsed_time_Windows_real(1)+Elapsed_time_start; 
end

