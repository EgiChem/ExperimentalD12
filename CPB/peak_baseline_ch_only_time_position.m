clc
clear all
close all
set(0,'DefaultFigureWindowStyle','docked')

% Channel only

%% User input
    pathtoread = 'Z:\OneDrive - Universidade de Aveiro\DAD\Testes_old_style\CPB\210715_Eucalyptol_CO2_etoh\210705_1200\210705_1200.xlsx';
    pathtosave = 'Z:\OneDrive - Universidade de Aveiro\DAD\Testes_old_style\CPB\210715_Eucalyptol_CO2_etoh\210705_1200\Resultados\CH 220 nm\';
    date = '2021_07_14';
    system = 'eucalyptol_etoh_w_duplicates';
    T = 303.15;  %K
    P = 1;     %bar
    x2mol = 0;  %cosolvent (mol/mol)
    x2vol = 0;  %cosolvent (v/v)
    Q = 0.150;      % caudal mL/min
    Tbath = 21+273.15; % T do banho em K
    
    % Injection time in seconds
     inject_time(1) = 4*60+43;
     inject_time(2) = 23*60+10;
     inject_time(3) = 40*60+13;



    start_peak = 500; 
    retention_time = 1500; % time for the peak compleatly come out in s

    Channel_no = 3;
    
%% Loading data and organizing stuff

[allpeakdata_num, allpeakdata_str, allpeakdata_cell] = xlsread(pathtoread);

Elapsed_time = allpeakdata_num(:,2);

Channel = allpeakdata_num(:,3:6);

Channel_w = allpeakdata_num(:,7:10);

Absolute_absorbance = allpeakdata_num(:,11:end);
 
% [Elapsed_time_Windows, Elapsed_time_Windows_real] = Elapsed_time_from_Windows(Elapsed_time, allpeakdata_str(:,1)); 
% Elapsed_time = Elapsed_time_Windows;
% % fix = Elapsed_time - Elapsed_time_Windows;
% % Elapsed_time = Elapsed_time - fix;




%% Calculos

% Channel

for p=1:numel(Channel_no)
    
    n = Channel_no(p);
    k = 1;
    while k <= numel(inject_time)
        
        if k == numel(inject_time)
            injection_e_index = numel(Elapsed_time);
        else
            injection_e_index = find(Elapsed_time==(inject_time(k)+retention_time));
        end
        
        injection_b_index = find(Elapsed_time==inject_time(k));
        time = Elapsed_time(injection_b_index:injection_e_index)-inject_time(k) ;
        signal = Channel(injection_b_index:injection_e_index,n);
        Wavelenght = Channel_w(injection_b_index,n);
        
%         [time_new signal_new] = check_duplicates(time, signal);
%         time = time_new;
%         signal = signal_new;
        
        Figuretitle = sprintf('Peak numer %s Wavelenght %s nm (CH)', num2str(k), num2str(Wavelenght));  
        
        figure(1)
        plot (time, signal, 'g.', [min(time) max(time)], [0 0], '--k')

        xlabel('time (s)')
        ylabel('Abs (uAU)')
        xlim([start_peak retention_time])
        set(gca,'XMinorTick','on','YMinorTick','on', 'FontSize', 15)
        title(Figuretitle)
%         hold on
%         [t(1), abs(1)] = ginput(1);
%         plot(t(1),abs(1),'xr')
%         [t(2), abs(2)] = ginput(1);
%         plot(t(2),abs(2),'xr')
%         [t(3), abs(3)] = ginput(1);
%         plot(t(3),abs(3),'xb')
%         [t(4), abs(4)] = ginput(1);
%         hold off
        
        
        t(1) = input('tempo 1 = ');
        t(2) = input('tempo 2 = ');
        t(3) = input('tempo 3 = ');
        t(4) = input('tempo 4 = ');
        
        
        
        t1 = find(time==round(t(1)));
        t2 = find(time==round(t(2)));
        t3 = find(time==round(t(3)));
        t4 = find(time==round(t(4)));
        
        while isempty(t1)==1
              t(1) = t(1) + 1;
              t1 = find(time==round(t(1)));
        end  
        while isempty(t2)==1
              t(2) = t(2) - 1;
              t2 = find(time==round(t(2)));
        end  
        while isempty(t3)==1
              t(3) = t(3) + 1;
              t3 = find(time==round(t(3)));
        end 
        while isempty(t4)==1
              t(4) = t(4) - 1;
              t4 = find(time==round(t(4)));
        end 
     
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
    
             t1 = find(time==round(t(1)));
             t2 = find(time==round(t(2)));
             t3 = find(time==round(t(3)));
             t4 = find(time==round(t(4)));
             
             while isempty(t1)==1
              t(1) = t(1) + 1;
              t1 = find(time==round(t(1)));
        end  
             while isempty(t2)==1
              t(2) = t(2) - 1;
              t2 = find(time==round(t(2)));
             end  
             while isempty(t3)==1
              t(3) = t(3) + 1;
              t3 = find(time==round(t(3)));
             end 
             while isempty(t4)==1
              t(4) = t(4) - 1;
              t4 = find(time==round(t(4)));
        end 
             
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
        plot( timepeak, signalbased, '.m', [0 timepeak(end)], [0 0], 'k-')
        xlabel('time (s)')
        ylabel('Abs (uAU)')
        xlim([min(timepeak) max(timepeak)]);
        legend('final peak', 'baseline fitting')
        legend('location', 'northeastoutside')
        title('Peak Based')
        set(gca,'XMinorTick','on','YMinorTick','on', 'FontSize', 15)
        
%         figure(3)
%         plot( timepeak, signalpeak, '.g', timepeak, baseline, '--r', timepeak, signalbased, '.m', [0 timepeak(end)], [0 0], 'k-')
%         xlabel('time (s)')
%         ylabel('Abs (uAU)')
%         xlim([min(timepeak) max(timepeak)]);
%         legend('original peak', 'baseline fitting', 'final peak')
%         legend('location', 'northeastoutside')
%         title(Figuretitle)
%         set(gca,'XMinorTick','on','YMinorTick','on', 'FontSize', 15)
        
        Cont = input('continue ??? (yes (press 1)/ no (press any other number))');
            if Cont==1
                
                Figurename = sprintf('%s_%s_Peak_numer_%s_of_channel_%s_Wavelenght_%s_nm.png', date, system, num2str(k), num2str(n), num2str(Wavelenght));
                Figurename_save = strcat(pathtosave,Figurename);
                saveas(gcf,Figurename_save)
        
                Filename = sprintf('%s_%s_Peak_numer_%s_of_channel_%s_Wavelenght_%s_nm.txt', date, system, num2str(k), num2str(n), num2str(Wavelenght));
                List_of_txt(k,:) = Filename;
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
                %
            end
  
%%
        close all 
        clear  injection_b_index  injection_e_index time signal  Wavelenght t t1 t2 t3 t4 timefit timepeak signalfit signalpeak B m b baseline signalbased Figuretitle Figurename Filename 

        fileID = fopen(strcat(pathtosave,'list_of_files.txt'),'w');
    for r=1:(k-1)
        write = List_of_txt(r,:);
        fprintf(fileID,'%s\n', write);
    end
    fclose(fileID);
    
    end 
end

display('Complete')


