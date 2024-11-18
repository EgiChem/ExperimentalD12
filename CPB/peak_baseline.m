clc
clear all
close all
set(0,'DefaultFigureWindowStyle','docked')

%% User input
    pathtoread = 'Z:\OneDrive - Universidade de Aveiro\DAD\teste_1_eucalyptol_212204 - Cópia.xlsx';
    pathtosave = 'Z:\OneDrive - Universidade de Aveiro\DAD\';
    date = '2021_04_22';
    system = 'eucalyptol_ethanol';
    T = 313.15;  %K
    P = 1;     %bar
    x2mol = 0;  %cosolvent (mol/mol)
    x2vol = 0;  %cosolvent (v/v)
    Q = 0.15;      % caudal mL/min
    Tbath = 21+273.15; % T do banho em K
    
    % Injection time in seconds
    inject_time(1) = 60*60+3*60+49;
    inject_time(2) = 60*60+16*60+29;
    inject_time(3) = 60*60+29*60+16;
    inject_time(4) = 60*60+50*60+50;
  
    
    retention_time = 1200; % time for the peak compleatly come out in s
    
    % Wavelenght range from 190 nm to 700 nm, the values must be pair
    Wave_min = 190;
    Wave_max = 194;

%% Loading data and organizing stuff

[allpeakdata_num, allpeakdata_str, allpeakdata_cell] = xlsread(pathtoread);

Elapsed_time = allpeakdata_num(:,2);

Channel = allpeakdata_num(:,3:6);

Channel_w = allpeakdata_num(:,7:10);

Absolute_absorbance = allpeakdata_num(:,11:end);

%% Calculos

% Channel

for n = 1:4
    
    k = 1;
    while k <= numel(inject_time)
        
             
        injection_b_index = find(Elapsed_time==inject_time(k));
        injection_e_index = find(Elapsed_time==(inject_time(k)+retention_time));
        time = Elapsed_time(injection_b_index:injection_e_index)-inject_time(k) ;
        signal = Channel(injection_b_index:injection_e_index,n);
        Wavelenght = Channel_w(injection_b_index,n);
        
        Figuretitle = sprintf('Peak numer %s Wavelenght %s nm (CH)', num2str(k), num2str(Wavelenght));  
        
        figure(1)
        plot (time, signal, 'b.', [min(time) max(time)], [0 0], '--k')
        hold on
        xlabel('time (s)')
        ylabel('Abs (uAU)')
        xlim([min(time) max(time)])
        title(Figuretitle)
        [t(1), abs(1)] = ginput(1);
        plot(t(1),abs(1),'xr')
        [t(2), abs(2)] = ginput(1);
        plot(t(2),abs(2),'xr')
        [t(3), abs(3)] = ginput(1);
        plot(t(3),abs(3),'xb')
        [t(4), abs(4)] = ginput(1);
        hold off

        t1 = round(t(1))+1;
        t2 = round(t(2))+1;
        t3 = round(t(3))+1;
        t4 = round(t(4))+1;
        
        try
            timefit = [time(t1:t2); time(t3:t4)];
            timepeak = time(t1: t4);
            signalfit = [signal(t1:t2); signal(t3:t4)];
            signalpeak = signal(t1: t4);
        catch
             warning('\n ERROR: If input error try again\n')
             display(lasterror)
             fprintf('\n')
             [t, abs] = ginput(4);
    
             t1 = round(t(1))+1;
             t2 = round(t(2))+1;
             t3 = round(t(3))+1;
             t4 = round(t(4))+1;
             timefit = [time(t1:t2); time(t3:t4)];
             timepeak = time(t1: t4);
             signalfit = [signal(t1:t2); signal(t3:t4)];
             signalpeak = signal(t1: t4);
        end

        [B,BINT,R,RINT] = regress(signalfit,[timefit, ones( numel(signalfit), 1)]);

        m = B(1);
        b = B(2);

        baseline = m.*timepeak + b;

        signalbased = signalpeak - baseline;
        
        figure(2)
        plot( timepeak, signalpeak, '.b', timepeak, baseline, '--r', timepeak, signalbased, '.m', [0 timepeak(end)], [0 0], 'k-')
        xlabel('time (s)')
        ylabel('Abs (uAU)')
        xlim([min(timepeak) max(timepeak)]);
        legend('original peak', 'baseline fitting', 'final peak')
        legend('location', 'northeastoutside')
        title(Figuretitle)
        
        Cont = input('continue ??? (yes (press 1)/ no (press any other number))');
            if Cont==1
                Figurename = sprintf('%s_%s_Peak_numer_%s_of_channel_%s_Wavelenght_%s_nm.png', date, system, num2str(k), num2str(n), num2str(Wavelenght));
                Figurename_save = strcat(pathtosave,Figurename);
                saveas(gcf,Figurename_save)
        
                Filename = sprintf('%s_%s_Peak_numer_%s_of_channel_%s_Wavelenght_%s_nm.txt', date, system, num2str(k), num2str(n), num2str(Wavelenght));
                Filename_save = strcat(pathtosave,Filename);
                fileID = fopen(Filename_save,'w');
                fprintf(fileID,'%3s\n', system);
                fprintf(fileID,'%3s\n', date);
                fprintf(fileID,'Injection time= %3s s\n', num2str(inject_time(k)));
                fprintf(fileID,'Wavelenght= %3s nm\n', num2str(Wavelenght));
                fprintf(fileID,'T= %3.2f K\n', T);
                fprintf(fileID,'P= %3.2f bar\n', P);
                fprintf(fileID,'xcosolvent= %1.4f (mol/mol)\n', x2mol);
                fprintf(fileID,'xcosolvent= %1.4f (v/v)\n', x2vol);
                fprintf(fileID,'Q= %1.5f mL/min\n', Q);
                fprintf(fileID,'Tbanho= %1.2f K\n', Tbath);
                fprintf(fileID,'%6s %20s %20s\n','time (s)','Signal original', 'Signal based');
                    for fuck = 1:numel(timepeak)
                        fprintf(fileID,'%6.2f %20.15f %20.15f\n', timepeak(fuck), signalpeak(fuck), signalbased(fuck));
                    end
                fclose(fileID);
                k = k + 1;
                fprintf('saved as %s \n', Filename)
            else
                % Cont = input('continue ??? (yes (press 1)/ no (press any other number))');
            end

%%           
    % Função para calcular D12 aqui????           
%%
        close all 
        clear  injection_b_index  injection_e_index time signal  Wavelenght t t1 t2 t3 t4 timefit timepeak signalfit signalpeak B m b baseline signalbased Figuretitle Figurename Filename 
    end 
end


% convert wavelenght to index

Wave_min_index =  (Wave_min-190)/2+1
Wave_max_index =  (Wave_max-190)/2+1


for n=Wave_min_index:Wave_max_index
    
    Wavelenght = (n-1)*2+190;
    
    k = 1;
    while k <= numel(inject_time)
        
        injection_b_index = find(Elapsed_time==inject_time(k));
        injection_e_index = find(Elapsed_time==(inject_time(k)+retention_time));
        time = Elapsed_time(injection_b_index:injection_e_index)-inject_time(k) ;
        signal = Absolute_absorbance(injection_b_index:injection_e_index,n);
        
        Figuretitle = sprintf('Peak numer %s Wavelenght %s nm (DAD)', num2str(k), num2str(Wavelenght));  
        
        figure(1)
        plot (time, signal, 'b.', [min(time) max(time)], [0 0], '--k')
        hold on
        xlabel('time (s)')
        ylabel('Abs (uAU)')
        xlim([min(time) max(time)])
        title(Figuretitle)
        [t(1), abs(1)] = ginput(1);
        plot(t(1),abs(1),'xr')
        [t(2), abs(2)] = ginput(1);
        plot(t(2),abs(2),'xr')
        [t(3), abs(3)] = ginput(1);
        plot(t(3),abs(3),'xb')
        [t(4), abs(4)] = ginput(1);
        hold off

        t1 = round(t(1))+1;
        t2 = round(t(2))+1;
        t3 = round(t(3))+1;
        t4 = round(t(4))+1;
        
        try
            timefit = [time(t1:t2); time(t3:t4)];
            timepeak = time(t1: t4);
            signalfit = [signal(t1:t2); signal(t3:t4)];
            signalpeak = signal(t1: t4);
        catch
             warning('\n ERROR: If input error try again\n')
             display(lasterror)
             fprintf('\n')
             [t, abs] = ginput(4);
    
             t1 = round(t(1))+1;
             t2 = round(t(2))+1;
             t3 = round(t(3))+1;
             t4 = round(t(4))+1;
             timefit = [time(t1:t2); time(t3:t4)];
             timepeak = time(t1: t4);
             signalfit = [signal(t1:t2); signal(t3:t4)];
             signalpeak = signal(t1: t4);
        end

        [B,BINT,R,RINT] = regress(signalfit,[timefit, ones( numel(signalfit), 1)]);

        m = B(1);
        b = B(2);

        baseline = m.*timepeak + b;

        signalbased = signalpeak - baseline;

        plot( timepeak, signalpeak, '.b', timepeak, baseline, '--r', timepeak, signalbased, '.m', [0 timepeak(end)], [0 0], 'k-')
        xlabel('time (s)')
        ylabel('Abs (uAU)')
        xlim([min(timepeak) max(timepeak)]);
        legend('original peak', 'baseline fitting', 'final peak')
        legend('location', 'northeastoutside')
        title(Figuretitle)
        
         Cont = input('continue ??? (yes (press 1)/ no (press any other number))')
            if Cont==1

                Figurename = sprintf('%s_%s_Peak_numer_%s_of_Wavelenght_%s_nm.png', date, system, num2str(k),  num2str(Wavelenght));
                Figurename_save = strcat(pathtosave,Figurename);
                saveas(gcf,Figurename_save)
        
                Filename = sprintf('%s_%s_Peak_numer_%s_of_Wavelenght_%s_nm.txt', date, system, num2str(k),  num2str(Wavelenght));
                Filename_save = strcat(pathtosave,Filename);
                fileID = fopen(Filename_save,'w');
                fprintf(fileID,'%3s\n', system);
                fprintf(fileID,'%3s\n', date);
                fprintf(fileID,'Injection time= %3s s\n', num2str(inject_time(k)));
                fprintf(fileID,'Wavelenght= %3s nm\n', num2str(Wavelenght));
                fprintf(fileID,'T= %3.2f K\n', T);
                fprintf(fileID,'P= %3.2f bar\n', P);
                fprintf(fileID,'xcosolvent= %1.4f (mol/mol)\n', x2mol);
                fprintf(fileID,'xcosolvent= %1.4f (v/v)\n', x2vol);
                fprintf(fileID,'Q= %1.5f mL/min\n', Q);
                fprintf(fileID,'Tbanho= %1.2f K\n', Tbath);
                fprintf(fileID,'%6s %12s %12s\n','time (s)','Signal original', 'Signal based');
                    for fuck = 1:numel(timepeak)
                        fprintf(fileID,'%6.2f %12.15f %12.15f\n',timepeak(fuck), signalpeak(fuck), signalbased(fuck));
                    end
                fclose(fileID);
                k = k + 1;
                fprintf('saved as %s', Filename)
            else
                % Cont = input('continue ??? (yes (press 1)/ no (press any other number))')
            end   
            
        close all 
        clear  injection_b_index  injection_e_index time signal  Wavelenght t t1 t2 t3 t4 timefit timepeak signalfit signalpeak B m b baseline signalbased Figuretitle Figurename Filename
    end 
    
    clear Wavelenght
end