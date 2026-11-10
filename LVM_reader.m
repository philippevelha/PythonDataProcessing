%% File reader for measurement done using PGC program
%Calculates the Power Spectral DEnsity
% If multiple files chosen creates a graph with the variation of mean
clear all variables
close all
ScanReady =1; %flag that tell if the file needs to be analyzed with textscan
DialogTitle ='Choose file(s) to start with';
[FileName,PathName,FilterIndex] = uigetfile('*.lvm','E:\MeasurementsDATA\',DialogTitle,'MultiSelect', 'on');
disp(FileName)
if isequal(FilterIndex,0)
    disp('User selected Cancel')
else
    disp(['User selected: ', fullfile(FileName)])
end
cd(PathName)
Fs = 1000;
keyword = 'nm'; %'T' for temperature
if iscell(FileName) %Means there is a multiple selection
    numFiles = length(FileName);
    for ii = 1:numFiles
        DATA = importdata(FileName{ii},'\t',23);%jumps the first 23 lines
        comments = DATA.textdata; %All the comment are put in this variable
        Measurement = DATA.data;
        sizM = size(Measurement);
        for jj = 1:sizM(2)/2
            eval(['X',num2str(jj),' = Measurement(:,',num2str(jj*2-1),');'])
            eval(['Y',num2str(jj),' = Measurement(:,',num2str(jj*2),');'])
        end
        
        
        indnan = isnan(Y1);%Remove non valid values
        Y1(indnan) = [];%Remove non valid values
        d1 = designfilt('lowpassiir','FilterOrder',12, ...
            'HalfPowerFrequency',0.0001,'DesignMethod','butter');
        y = filtfilt(d1,Y1);
        %plot(y)
        M(ii) = mean(Y1);
        minstd=10;
        for jj=1:round(length(Y1)/100)
            if minstd>std(Y1(1+100*(jj-1):100*jj))
                minstd=std(Y1(1+100*(jj-1):100*jj));
            end
        end
        ST(ii)=minstd;
        indT = strfind(FileName{ii},keyword);
        indC = strfind(FileName{ii},'C');
        %         Temperature(ii) = str2double(FileName{ii}(indT+1:indC(2)-1));
        indF = strfind(FileName{ii},'kHz');
        if isempty(indF)
            Fs = 1000;
        else
            if str2double(FileName{ii}(indF-1)) == 0
                Fs = 10000;
            else
                Fs = str2double(FileName{ii}(indF-1))*1000;
            end
        end
        %% Calculation of the Power Spectral Density
        N = length(Y1-y);
        xdft = fft(Y1-y);
        xdft = xdft(1:N/2+1);
        psdx = (1/(Fs*N)) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        freq = 0:Fs/N:Fs/2;
        %freq = 0:(2*pi)/N:pi;
        
        semilogx(freq,10*log10(psdx))
        hold on
        
        grid on
        title('Periodogram Using FFT')
        xlabel('Frequency (Hz)')
        ylabel('Power/Frequency (dB/Hz)')
    end
    legend(FileName)
    
    figure
    errorbar(Temperature, M,ST)
    
    figure
    plot(Temperature,ST)
else
    if ScanReady
        fileID = fopen(FileName);
        flag_exit = 1;
        skipheader(fileID);
        X1=[];
        Y1 = [];
        X2 = [];
        Y2 = [];
        X3 = [];
        Y3 = [];
        X4 = [];
        Y4 = [];
        lastvalue = [0,0,0,0,0,0,0,0];
        while flag_exit
            skipheader(fileID);%Do it twice AT THEN BEGINING
            
            headers = textscan(fileID,'%s %s %s %s %s %s %s %s %s',1); %CONTAINS HEADERS OF THE DATA
            
            DATA = textscan(fileID,'%f %f %f %f %f %f %f %f %*[^\n]');%read data columns remove possible extra at the end of the line
            X1 = [X1;DATA{1}+lastvalue(1)]; Y1 = [Y1;DATA{2}]; %concatenate data
            X2 = [X2;DATA{3}+lastvalue(3)]; Y2 = [Y2;DATA{4}];
            X3 = [X3;DATA{5}+lastvalue(5)]; Y3 = [Y3;DATA{6}];
            X4 = [X4;DATA{7}+lastvalue(7)]; Y4 = [Y4;DATA{8}];
            lastvalue = [X1(end),Y1(end),X2(end),Y2(end),X3(end),Y3(end),X4(end),Y4(end)];
            C_text = textscan(fileID,'%s',1);
            if isempty(C_text{1})
                %reached the end
                flag_exit = 0;
            else
                %continue
            end
        end
        fclose(fileID);
    else
        DATA = importdata(FileName,'\t',23);%jumps the first 23 lines
        
        comments = DATA.textdata; %All the comment are put in this variable
        Measurement = DATA.data;
        
        X1 = Measurement(:,1);
        Y1 = Measurement(:,2);
        X2 = Measurement(:,3);
        Y2 = Measurement(:,4);
    end
    Y1 = unwrap(Y1);
    Y2 = unwrap(Y2);
    figure;plot(X1,Y1)
    hold on
    plot(X2,Y2)
    xlabel('time [s]')
    ylabel('phase [°]')
    title(FileName)
    indnan = isnan(Y1);
    Y1(indnan) = [];
    X1(indnan) = [];
    N = length(Y1);
    xdft = fft(Y1);
    xdft = xdft(1:N/2+1);
    psdx = (1/(2*pi*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:(2*pi)/N:pi;
    figure
    semilogx(freq/pi,10*log10(psdx))
    grid on
    axis([0 1 -80 20])
    title('Periodogram Using FFT')
    xlabel('Normalized Frequency (\times\pi rad/sample)')
    ylabel('Power/Frequency (dB/rad/sample)')
end

function skipheader(fileID)
flag_exit = 1;
while flag_exit %changes to zero on finding exit condition
    C_text = textscan(fileID,'%s',1);
    if strcmp(C_text{1}{1},'***End_of_Header***');flag_exit = 0; end
end
end