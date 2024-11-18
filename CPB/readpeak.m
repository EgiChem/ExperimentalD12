function [system, date, Inject_time, Wavelenght, T, P, x2mol, x2vol, Q, Tbath, data] = readpeak(Filepath)

%load peaks from txt test


% %fid = fopen('D:\PhD\DAD\2021_03_22_test3_Peak_numer_1_of_channel_1_Wavelenght_5_nm.txt');
% fid = fopen('Z:\OneDrive - Universidade de Aveiro\DAD\2021_04_22_test3_sofia_Peak_numer_3_of_channel_1_Wavelenght_5_nm.txt');

fid = fopen(Filepath);

system = fgetl(fid);
date = fgetl(fid);
Inject_time = fgetl(fid);
Wavelenght = fgetl(fid);
T = fgetl(fid);
P = fgetl(fid);
x2mol = fgetl(fid);
x2vol = fgetl(fid);
Q = fgetl(fid);
Tbath = fgetl(fid);

%% extract values from file
Inject_time(strfind(Inject_time, '=')) = [];
Key   = 'time';
Index = strfind(Inject_time, Key);
Inject_time = sscanf(Inject_time(Index(1) + length(Key):end), '%g', 1); % in s

Wavelenght(strfind(Wavelenght, '=')) = [];
Key   = 'Wavelenght';
Index = strfind(Wavelenght, Key);
Wavelenght = sscanf(Wavelenght(Index(1) + length(Key):end), '%g', 1); % in s

T(strfind(T, '=')) = [];
Key   = 'T';
Index = strfind(T, Key);
T = sscanf(T(Index(1) + length(Key):end), '%g', 1); % in s

P(strfind(P, '=')) = [];
Key   = 'P';
Index = strfind(P, Key);
P = sscanf(P(Index(1) + length(Key):end), '%g', 1); % in s

x2mol(strfind(x2mol, '=')) = [];
Key   = 'xcosolvent';
Index = strfind(x2mol, Key);
x2mol = sscanf(x2mol(Index(1) + length(Key):end), '%g', 1); % in s

x2vol(strfind(x2vol, '=')) = [];
Key   = 'xcosolvent';
Index = strfind(x2vol, Key);
x2vol = sscanf(x2vol(Index(1) + length(Key):end), '%g', 1); % in s

Q(strfind(Q, '=')) = [];
Key   = 'Q';
Index = strfind(Q, Key);
Q = sscanf(Q(Index(1) + length(Key):end), '%g', 1); % in s

Tbath(strfind(Tbath, '=')) = [];
Key   = 'Tbanho';
Index = strfind(Tbath, Key);
Tbath = sscanf(Tbath(Index(1) + length(Key):end), '%g', 1); % in s

tline = fgets(fid);

%%

tline_num = 1;
n = 1; 

while tline_num~=-1
    tline_num = fgets(fid);
    if ~ischar(tline_num)
        break
    end
    
    data(n,:) = str2num(tline_num); % colun 1 - time; column 2 - original signal; column 3 - modified signal
        
    n = n + 1;
end

