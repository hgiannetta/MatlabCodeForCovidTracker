% -------------------------------------
%   script for Spaceapp 2020 Coronavirus
%   by Hernan Giannetta 
%   email: hgiannetta@hotmail.com
%   30/05/2020
% -------------------------------------

loadByFiles=0; %loadByFiles=1 or loadByDabaBase=0

%--------------
clc;
if loadByFiles
    clear all;
    loadByFiles=1;
end
close all;

%VectorFilesTORead=1:length(FolderFilesNames)
VectorFilesTORead=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]

plotRawData=0; % Ploting the raw data
plotTestAData=0; % Ploting the raw data
plotTestBData=0; % Ploting the raw data
plotTestCData=1; % Ploting the raw data
plotTestDData=1; % Ploting the raw data
plotTestEData=1; % Ploting the raw data
plotTestFData=1; % Ploting the raw data
plotTestGData=1; % Ploting the raw data
plotTestHData=1; % Ploting the raw data
plotTestiData=1; % Ploting the raw data
plotTestjData=1; % Ploting the raw data
%-------------

lb5=LibAv05;
CurVar.titulo=  'Script de compilación global de cada imagen';
CurVar.autor=   'Hernan Giannetta';
CurVar.fecha=   'Jun 2017';
% Load Lib
VarInit=lb5.loadInitVariables(CurVar.titulo,CurVar.autor,CurVar.fecha); 
CurVar.message = lb5.DispLibVersion();  
CurVar.fechaInit=clock; 
lb5.DispCommWindwPortada(CurVar.titulo,CurVar.autor,CurVar.fecha);

% Variables
% experimental data
DirFileNameToConvert= './SweepData';
dhindex=0; %figure handler
dh=0; %figure handler

%------------------------------------
% DB tag for everi Signal measured
%------------------------------------
% lleno texto a leer
i=0;

i=i+1;
legendText{i}='01-(Temp 35.6°C) With Congestion(1-20Khz)'; %ImpedanciaOidoTempCorporal35_6CconMoco_Sweep30s_1_20000.ogg
i=i+1;
legendText{i}='02-(Temp 35.6°C) With Congestion(1-20Khz)'; %ImpedanciaOidoTempCorporal35_6CconMoco_Sweep30sde1a20000.ogg
i=i+1;
legendText{i}='03-(Temp 36.2°C) Without Congestion(1-20Khz)'; %ImpedanciaOidoTempCorporal35_6CconMoco_Sweep30s_1_20000.ogg
i=i+1;
legendText{i}='04-(Temp 36.2°C) Without Congestion(1-20Khz)'; %ImpedanciaOidoTempCorporal35_6CconMoco_Sweep30s_1_20000.ogg
i=i+1;
legendText{i}='05-(Temp 36.1°C) With Congestion (1-22Khz)';  % '39-With  Congestion (1-22Khz)';   
i=i+1;
legendText{i}='06-(Temp 36.1°C) With Congestion (1-22Khz)';   %'40-With Congestion (1-21.5Khz)';  
i=i+1;
legendText{i}='07-(Temp 36.1°C) Without Congestion (1-22Khz)';   %'41-Without Congestion (1-22Khz)';   
i=i+1;
legendText{i}='08-(Temp 35.7°C) Calibration Signal (1-22Khz)';   %':008_SweepImp_Temp35_7_Calibración_30s_de1a22000.ogg
i=i+1;
legendText{i}='09-(Temp 35.7°C) Calibration Signal (1-22Khz)';   %':008_SweepImp_Temp35_7_Calibración_30s_de1a22000.ogg
i=i+1;
legendText{i}='10-(Temp 35.7°C) With Congestion (1-22Khz)';   %':008_SweepImp_Temp35_7_Calibración_30s_de1a22000.ogg
i=i+1;
legendText{i}='11-(Temp 35.7°C) With Congestion (22-1Khz)';   %':008_SweepImp_Temp35_7_Calibración_30s_de1a22000.ogg
i=i+1;
legendText{i}='12-(Temp 35.7°C) Without Congestion (22-1Khz)';   %':008_SweepImp_Temp35_7_Calibración_30s_de1a22000.ogg
i=i+1;
legendText{i}='13-(Temp 35.7°C) Calibration Signal (22-1Khz)';   %':008_SweepImp_Temp35_7_Calibración_30s_de1a22000.ogg
i=i+1;
legendText{i}='14-(Temp 35.7°C) Without Congestion (1-22Khz)';   %':008_SweepImp_Temp35_7_Calibración_30s_de1a22000.ogg
i=i+1;
legendText{i}='15-(Temp 35.7°C) Calibration Signal (1-22Khz)';   %':008_SweepImp_Temp35_7_Calibración_30s_de1a22000.ogg


%------------------------------------
% Main
%------------------------------------

if (loadByFiles)
    
    DB.filenameAtributes.AllFilesDir=dir(DirFileNameToConvert);
     j=0;
     for i=1:length(DB.filenameAtributes.AllFilesDir)
         % print mensaje
         LibAv05.DispCommWindwText(strcat...
             ('Fill data struct (',num2str(i),') : ',DB.filenameAtributes.AllFilesDir(i).name));
         if(not(isdir(DB.filenameAtributes.AllFilesDir(i).name))) 
            j=j+1;
            FolderFilesNames{j}=DB.filenameAtributes.AllFilesDir(i).name;
         end                 
     end
     
     for i=VectorFilesTORead
        LibAv05.DispCommWindwText(strcat...
                 ('Opening data file: (',num2str(i),') : ',FolderFilesNames{i}));
        [y,Fs]=audioread(fullfile(DirFileNameToConvert,...
            FolderFilesNames{i}));   
        DB.RawAudioData{i}.info = audioinfo(fullfile(DirFileNameToConvert,...
            FolderFilesNames{i}));
        DB.RawAudioData{i}.FileName=FolderFilesNames{i};
        DB.RawAudioData{i}.data_y= y;  
           t = 0:seconds(1/DB.RawAudioData{i}.info.SampleRate):seconds(length(y)/Fs);
        DB.RawAudioData{i}.data_t = t(1:end-1);
        DB.RawAudioData{i}.data_x= 0:1/length(y):length(y)/...
                                DB.RawAudioData{i}.info.SampleRate;  
        DB.RawAudioData{i}.Fs = Fs;  
        DB.filenameAtributes.filename{i}=FolderFilesNames{i};
        DB.legendText{i}=legendText{i};
     end
else
%     LibAv05.DispCommWindwText(strcat('Loading DB please wait'));
%     load('DBExpCurv.mat','DC' )
%     LibAv05.DispCommWindwText(strcat('Data loading done '));
%     
%     DB.RawAudioData.info = DC.info;
%     DB.RawAudioData.FileName=DC.FileName;
%     DB.RawAudioData.data_y= DC.data_y;   
%     DB.RawAudioData.data_t = DC.data_t;
%     DB.RawAudioData.data_x= DC.data_x;  
%     DB.RawAudioData.Fs = DC.Fs;  
%     DB.filenameAtributes.filename=DC.filename;
    
end


if (plotRawData)
 %graficando raw data
            for i=1:length(DB.RawAudioData)
                dhindex=i;
                dh(dhindex)=figure;
                set(dh(dhindex),'Name', DB.RawAudioData{i}.FileName);
                ax1{i} = subplot(1,2,1); % top subplot
                plot(ax1{i},DB.RawAudioData{i}.data_t,DB.RawAudioData{i}.data_y)
                xlabel(ax1{i},'Time')
                ylabel(ax1{i},'Audio')
%                 title(ax1{i},DB.RawAudioData{i}.FileName)
                % Calculo el FFT                
                N = length(DB.RawAudioData{i}.data_y); 
                dt = 1/DB.RawAudioData{i}.Fs;
                t = dt*(0:N-1)';
                datafft=fft(DB.RawAudioData{i}.data_y);
                datafft_abs=abs(datafft/N);
                datafft_abs=datafft_abs(1:N/2+1);
                freq=Fs*(0:N/2)/N;
                ax2{i} = subplot(1,2,2); % top subplot
                plot(ax2{i},freq,datafft_abs);
                xlabel(ax2{i},'f[Hz]');
                ylabel(ax2{i},'Amp');
%                 title(ax2{i},strcat('FFT:',DB.RawAudioData{i}.FileName));
                title(ax2{i},strcat('FFT:',DB.RawAudioData{i}.FileName));


%                 plot (DB.RawAudioData{i}.data_x,DB.RawAudioData{i}.data_y)
            end
end      
            
            
% TEst A
if (plotTestAData)
VectorToRead_Comp_1=[1 3];    
dhindex=dhindex+1;
dh(dhindex)=figure;         
k=0;
   VectorToRead=VectorToRead_Comp_1;
            %graficando raw data
            for i=VectorToRead
                
                % Señal temporal
                k=k+1;
                ax1{i} = subplot(length(VectorToRead),2,k); % top subplot
                plot(ax1{i},DB.RawAudioData{i}.data_t,DB.RawAudioData{i}.data_y)
                xlabel(ax1{i},'Time')
                ylabel(ax1{i},'Audio')
                title(ax1{i},DB.legendText{i})
                
                % Calculo el FFT   
                k=k+1;
                N = length(DB.RawAudioData{i}.data_y); 
                dt = 1/DB.RawAudioData{i}.Fs;
                t = dt*(0:N-1)';
                datafft=fft(DB.RawAudioData{i}.data_y);
                datafft_abs=abs(datafft/N);
                datafft_abs=datafft_abs(1:N/2+1);
                freq=Fs*(0:N/2)/N;
                ax2{i} = subplot(length(VectorToRead),2,k); % top subplot
                plot(ax2{i},freq,datafft_abs);
                xlabel(ax2{i},'f[Hz]');
                ylabel(ax2{i},'Amp');
                title(ax2{i},strcat(DB.legendText{i}));

            end
end


%TEst B
if(plotTestBData)
   VectorToRead_Comp_1=[5 6 7]; 
    dhindex=dhindex+1;
   dh(dhindex)=figure;          
   k=0;
   VectorToRead=VectorToRead_Comp_1;
            %graficando raw data
            for i=VectorToRead
                
                % Señal temporal
                k=k+1;
                ax1{i} = subplot(length(VectorToRead),2,k); % top subplot
                plot(ax1{i},DB.RawAudioData{i}.data_t,DB.RawAudioData{i}.data_y)
                xlabel(ax1{i},'Time')
                ylabel(ax1{i},'Audio')
                title(ax1{i},DB.legendText{i})
                
                % Calculo el FFT   
                k=k+1;
                N = length(DB.RawAudioData{i}.data_y); 
                dt = 1/DB.RawAudioData{i}.Fs;
                t = dt*(0:N-1)';
                datafft=fft(DB.RawAudioData{i}.data_y);
                datafft_abs=abs(datafft/N);
                datafft_abs=datafft_abs(1:N/2+1);
                freq=Fs*(0:N/2)/N;
                ax2{i} = subplot(length(VectorToRead),2,k); % top subplot
                plot(ax2{i},freq,datafft_abs);
                xlabel(ax2{i},'f[Hz]');
                ylabel(ax2{i},'Amp');
                title(ax2{i},strcat(DB.legendText{i}));

            end
end


%TestC
if (plotTestCData)
            
%----------------------
FigureText='Ploting Calibratios Signals';
%----------------------
VectorToRead_Comp_1=[8 9 13 15]; 
dhindex=dhindex+1;
dh(dhindex)=figure;         
set(dh(dhindex),'Name',FigureText);
k=0;
VectorToRead=VectorToRead_Comp_1;
            %graficando raw data
            for i=VectorToRead
                
                % Señal temporal
                k=k+1;
                ax1{i} = subplot(length(VectorToRead),2,k); % top subplot
                plot(ax1{i},DB.RawAudioData{i}.data_t,DB.RawAudioData{i}.data_y)
                xlabel(ax1{i},'Time')
                ylabel(ax1{i},'Audio')
                title(ax1{i},DB.legendText{i})
                
                % Calculo el FFT   
                k=k+1;
                N = length(DB.RawAudioData{i}.data_y); 
                dt = 1/DB.RawAudioData{i}.Fs;
                t = dt*(0:N-1)';
                datafft=fft(DB.RawAudioData{i}.data_y);
                datafft_abs=abs(datafft/N);
                datafft_abs=datafft_abs(1:N/2+1);
                freq=Fs*(0:N/2)/N;
                ax2{i} = subplot(length(VectorToRead),2,k); % top subplot
                plot(ax2{i},freq,datafft_abs);
                xlabel(ax2{i},'f[Hz]');
                ylabel(ax2{i},'Amp');
                title(ax2{i},strcat(DB.legendText{i}));

            end
            
end

% TEst D
if(plotTestDData)
    %----------------------
    FigureText='Ploting Test1 Signals';
    %----------------------
    VectorToRead_Comp_1=[8 9 10 11]; 
    dhindex=dhindex+1;
    dh(dhindex)=figure;         
    set(dh(dhindex),'Name',FigureText);
    k=0;
    VectorToRead=VectorToRead_Comp_1;
            %graficando raw data
            for i=VectorToRead
                
                % Señal temporal
                k=k+1;
                ax1{i} = subplot(length(VectorToRead),2,k); % top subplot
                plot(ax1{i},DB.RawAudioData{i}.data_t,DB.RawAudioData{i}.data_y)
                xlabel(ax1{i},'Time')
                ylabel(ax1{i},'Audio')
                title(ax1{i},DB.legendText{i})
                
                % Calculo el FFT   
                k=k+1;
                N = length(DB.RawAudioData{i}.data_y); 
                dt = 1/DB.RawAudioData{i}.Fs;
                t = dt*(0:N-1)';
                datafft=fft(DB.RawAudioData{i}.data_y);
                datafft_abs=abs(datafft/N);
                datafft_abs=datafft_abs(1:N/2+1);
                freq=Fs*(0:N/2)/N;
                ax2{i} = subplot(length(VectorToRead),2,k); % top subplot
                plot(ax2{i},freq,datafft_abs);
                xlabel(ax2{i},'f[Hz]');
                ylabel(ax2{i},'Amp');
                title(ax2{i},strcat(DB.legendText{i}));

            end
            
end


% TEst E
if(plotTestEData)
    %----------------------
    FigureText='Ploting Test2 Signals';
    %----------------------
    VectorToRead_Comp_1=[12 13 14 15]; 
    dhindex=dhindex+1;
    dh(dhindex)=figure;         
    set(dh(dhindex),'Name',FigureText);
    k=0;
    VectorToRead=VectorToRead_Comp_1;
                %graficando raw data
            for i=VectorToRead
                
                % Señal temporal
                k=k+1;
                ax1{i} = subplot(length(VectorToRead),2,k); % top subplot
                plot(ax1{i},DB.RawAudioData{i}.data_t,DB.RawAudioData{i}.data_y)
                xlabel(ax1{i},'Time')
                ylabel(ax1{i},'Audio')
                title(ax1{i},DB.legendText{i})
                                
                % Calculo el FFT   
                k=k+1;
                N = length(DB.RawAudioData{i}.data_y); 
                dt = 1/DB.RawAudioData{i}.Fs;
                t = dt*(0:N-1)';
                datafft=fft(DB.RawAudioData{i}.data_y);
                datafft_abs=abs(datafft/N);
                datafft_abs=datafft_abs(1:N/2+1);
                freq=Fs*(0:N/2)/N;
                ax2{i} = subplot(length(VectorToRead),2,k); % top subplot
                plot(ax2{i},freq,datafft_abs);
                xlabel(ax2{i},'f[Hz]');
                ylabel(ax2{i},'Amp');
                title(ax2{i},strcat(DB.legendText{i}));

            end
end



            
% TEst F
if(plotTestFData)         
    %----------------------
    FigureText='Ploting Test3 Diferencia FFT';
    %----------------------
    VectorToRead=[12 13]; 
    dhindex=dhindex+1;
    dh(dhindex)=figure;         
    set(dh(dhindex),'Name',FigureText);
    k=0;
    %VectorToRead=VectorToRead_Comp_1;
    %graficando raw data
    for i=VectorToRead    
        k=k+1;
        frec=DB.RawAudioData{i}.data_t;
        data_y=DB.RawAudioData{i}.data_y;
        L_Comp_FFT =lb5.Compute_FFT(frec,data_y);
%         N = length(data_y); 
%         dt = 1/L_Comp_FFT.Fs;
%         t = dt*(0:N-1)';
%         L_Comp_FFT.FrecCorrected=L_Comp_FFT.Fs*(0:N/2)/N;        
        L_Compute_FFT{k}=L_Comp_FFT;
    end
    
    k=0;
    for i=1:length(L_Compute_FFT)

        % los FFT data 
%         Dt.PlotData.Figure_B.f(:,i)=L_Compute_FFT{i}.f';
% %         Dt.PlotData.Figure_B.f(:,i)=L_Compute_FFT{i}.FrecCorrected;
%         % FFT amplitude magnitude
%         Dt.PlotData.Figure_B.P1(:,i)=L_Compute_FFT{i}.P1;
%         % FFT fase
%         Dt.PlotData.Figure_B.P3(:,i)=L_Compute_FFT{i}.P3;

            %--------------------------------
        % DATOS FIltro Resampling  Sabinsky - Golay
        DataInputY_Mod=L_Compute_FFT{i}.P1;
        DataInputY_Phase=unwrap(L_Compute_FFT{i}.P3);
        DataInputX=L_Compute_FFT{i}.f';

        % FIltro Resampling  Sabinsky - Golay Parametros
        InputX.DataInputX=DataInputX;
        InputX.fsData = length(DataInputX); % Sampling frequency
        InputX.fsResamp=InputX.fsData*3; %resampling frec
        InputX.filtOrder=25;
        InputX.filtPoints=499;        
        InputX.EnablePlot_FuncInt=0;   
%         InputX.xlim=VariablesProc.Figure_B.xlim;
%         InputX.ylim=VariablesProc.Figure_B.ylim;
        %TestMAResamp = lb5.ResampSabinGolay(DataInputY,InputX);
        
        
        % Calculo FIltro Resampling  Sabinsky - Golay          
        FFT_Filt_Mod=lb5.ResampSabinGolay(DataInputY_Mod,InputX); 
        FFT_Filt_Phase=lb5.ResampSabinGolay(DataInputY_Phase,InputX); 
%         % FFT amplitude magnitude
%         Dt.PlotData.Figure_B.P1(:,i)= FFT_Filt_Mod.Resamp_Y;
%         % FFT fase
%         Dt.PlotData.Figure_B.P3(:,i)=  FFT_Filt_Phase.Resamp_Y;  
        
        
            % Señal temporal
        k=k+1;
        ax1{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax1{i},DataInputX,DataInputY_Mod);
%         plot(ax1{i},FFT_Filt_Mod.Resamp_X,FFT_Filt_Mod.Resamp_Y);        
        xlabel(ax1{i},'f[Hz]')
        ylabel(ax1{i},'Amp.')
        title(ax1{i},DB.legendText{i})

        % Calculo el FFT   
        k=k+1;
        ax2{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax2{i},DataInputX,DataInputY_Phase);
%         plot(ax2{i},FFT_Filt_Phase.Resamp_X,FFT_Filt_Phase.Resamp_Y);
        
        xlabel(ax2{i},'f[Hz]');
        ylabel(ax2{i},'Fase [deg]');
        title(ax2{i},strcat(DB.legendText{i}));
    end
          
        DataInputY_Mod=L_Compute_FFT{1}.P1(1:length(L_Compute_FFT{2}.P1),1) ./ L_Compute_FFT{2}.P1;
        DataInputY_Phase=unwrap(L_Compute_FFT{1}.P3(1:length(L_Compute_FFT{2}.P3),1))-unwrap(L_Compute_FFT{2}.P3);
        DataInputX=L_Compute_FFT{2}.f';
        
        DIfIndex=1;
        SignalDiferencia(DIfIndex).DataInputY_Mod=DataInputY_Mod;
        SignalDiferencia(DIfIndex).DataInputY_Phase=DataInputY_Phase;
        SignalDiferencia(DIfIndex).DataInputX=DataInputX;

        k=k+1;
        ax1{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax1{i},DataInputX,DataInputY_Mod);
%         plot(ax1{i},FFT_Filt_Mod.Resamp_X,FFT_Filt_Mod.Resamp_Y);        
        xlabel(ax1{i},'f[Hz]')
        ylabel(ax1{i},'Amp.')
        title(ax1{i},DB.legendText{i})

        % Calculo el FFT   
        k=k+1;
        ax2{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax2{i},DataInputX,DataInputY_Phase);
%         plot(ax2{i},FFT_Filt_Phase.Resamp_X,FFT_Filt_Phase.Resamp_Y);  
        xlabel(ax2{i},'f[Hz]');
        ylabel(ax2{i},'Fase [deg]');
        title(ax2{i},strcat(DB.legendText{i}));
        

end   
            




% TEst G
if(plotTestGData)         
    %----------------------
    FigureText='Ploting Test4 Diferencia FFT';
    %----------------------
    VectorToRead=[14 15]; 
    dhindex=dhindex+1;
    dh(dhindex)=figure;         
    set(dh(dhindex),'Name',FigureText);
    k=0;
    %VectorToRead=VectorToRead_Comp_1;
    %graficando raw data
    for i=VectorToRead    
        k=k+1;
        frec=DB.RawAudioData{i}.data_t;
        data_y=DB.RawAudioData{i}.data_y;
        L_Comp_FFT =lb5.Compute_FFT(frec,data_y);
%         N = length(data_y); 
%         dt = 1/L_Comp_FFT.Fs;
%         t = dt*(0:N-1)';
%         L_Comp_FFT.FrecCorrected=L_Comp_FFT.Fs*(0:N/2)/N;        
        L_Compute_FFT{k}=L_Comp_FFT;
    end
    
    k=0;
    for i=1:length(L_Compute_FFT)

        % los FFT data 
%         Dt.PlotData.Figure_B.f(:,i)=L_Compute_FFT{i}.f';
% %         Dt.PlotData.Figure_B.f(:,i)=L_Compute_FFT{i}.FrecCorrected;
%         % FFT amplitude magnitude
%         Dt.PlotData.Figure_B.P1(:,i)=L_Compute_FFT{i}.P1;
%         % FFT fase
%         Dt.PlotData.Figure_B.P3(:,i)=L_Compute_FFT{i}.P3;

            %--------------------------------
        % DATOS FIltro Resampling  Sabinsky - Golay
        DataInputY_Mod=L_Compute_FFT{i}.P1;
        DataInputY_Phase=unwrap(L_Compute_FFT{i}.P3);
        DataInputX=L_Compute_FFT{i}.f';

        % FIltro Resampling  Sabinsky - Golay Parametros
        InputX.DataInputX=DataInputX;
        InputX.fsData = length(DataInputX); % Sampling frequency
        InputX.fsResamp=InputX.fsData*3; %resampling frec
        InputX.filtOrder=25;
        InputX.filtPoints=499;        
        InputX.EnablePlot_FuncInt=0;   
%         InputX.xlim=VariablesProc.Figure_B.xlim;
%         InputX.ylim=VariablesProc.Figure_B.ylim;
        %TestMAResamp = lb5.ResampSabinGolay(DataInputY,InputX);
        
        
        % Calculo FIltro Resampling  Sabinsky - Golay          
        FFT_Filt_Mod=lb5.ResampSabinGolay(DataInputY_Mod,InputX); 
        FFT_Filt_Phase=lb5.ResampSabinGolay(DataInputY_Phase,InputX); 
%         % FFT amplitude magnitude
%         Dt.PlotData.Figure_B.P1(:,i)= FFT_Filt_Mod.Resamp_Y;
%         % FFT fase
%         Dt.PlotData.Figure_B.P3(:,i)=  FFT_Filt_Phase.Resamp_Y;  
        
        
            % Señal temporal
        k=k+1;
        ax1{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax1{i},DataInputX,DataInputY_Mod);
%         plot(ax1{i},FFT_Filt_Mod.Resamp_X,FFT_Filt_Mod.Resamp_Y);        
        xlabel(ax1{i},'f[Hz]')
        ylabel(ax1{i},'Amp.')
        title(ax1{i},DB.legendText{i})

        % Calculo el FFT   
        k=k+1;
        ax2{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax2{i},DataInputX,DataInputY_Phase);
%         plot(ax2{i},FFT_Filt_Phase.Resamp_X,FFT_Filt_Phase.Resamp_Y);
        
        xlabel(ax2{i},'f[Hz]');
        ylabel(ax2{i},'Fase [deg]');
        title(ax2{i},strcat(DB.legendText{i}));
    end
          
        DataInputY_Mod=L_Compute_FFT{1}.P1  ./ L_Compute_FFT{2}.P1(1:length(L_Compute_FFT{1}.P1),1);
        DataInputY_Phase=unwrap(L_Compute_FFT{1}.P3)-unwrap(L_Compute_FFT{2}.P3(1:length(L_Compute_FFT{1}.P3),1));
        DataInputX=L_Compute_FFT{1}.f';
        
        DIfIndex=2;
        SignalDiferencia(DIfIndex).DataInputY_Mod=DataInputY_Mod;
        SignalDiferencia(DIfIndex).DataInputY_Phase=DataInputY_Phase;
        SignalDiferencia(DIfIndex).DataInputX=DataInputX;

        k=k+1;
        ax1{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax1{i},DataInputX,DataInputY_Mod);
%         plot(ax1{i},FFT_Filt_Mod.Resamp_X,FFT_Filt_Mod.Resamp_Y);        
        xlabel(ax1{i},'f[Hz]')
        ylabel(ax1{i},'Amp.')
        title(ax1{i},DB.legendText{i})

        % Calculo el FFT   
        k=k+1;
        ax2{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax2{i},DataInputX,DataInputY_Phase);
%         plot(ax2{i},FFT_Filt_Phase.Resamp_X,FFT_Filt_Phase.Resamp_Y);  
        xlabel(ax2{i},'f[Hz]');
        ylabel(ax2{i},'Fase [deg]');
        title(ax2{i},strcat(DB.legendText{i}));       
end   



% TEst G
if(plotTestGData)         
    %----------------------
    FigureText='Ploting Test4 Diferencia FFT';
    %----------------------
    VectorToRead=[8 10]; 
    dhindex=dhindex+1;
    dh(dhindex)=figure;         
    set(dh(dhindex),'Name',FigureText);
    k=0;
    %VectorToRead=VectorToRead_Comp_1;
    %graficando raw data
    for i=VectorToRead    
        k=k+1;
        frec=DB.RawAudioData{i}.data_t;
        data_y=DB.RawAudioData{i}.data_y;
        L_Comp_FFT =lb5.Compute_FFT(frec,data_y);
%         N = length(data_y); 
%         dt = 1/L_Comp_FFT.Fs;
%         t = dt*(0:N-1)';
%         L_Comp_FFT.FrecCorrected=L_Comp_FFT.Fs*(0:N/2)/N;        
        L_Compute_FFT{k}=L_Comp_FFT;
    end
    
    k=0;
    for i=1:length(L_Compute_FFT)

        % los FFT data 
%         Dt.PlotData.Figure_B.f(:,i)=L_Compute_FFT{i}.f';
% %         Dt.PlotData.Figure_B.f(:,i)=L_Compute_FFT{i}.FrecCorrected;
%         % FFT amplitude magnitude
%         Dt.PlotData.Figure_B.P1(:,i)=L_Compute_FFT{i}.P1;
%         % FFT fase
%         Dt.PlotData.Figure_B.P3(:,i)=L_Compute_FFT{i}.P3;

            %--------------------------------
        % DATOS FIltro Resampling  Sabinsky - Golay
        DataInputY_Mod=L_Compute_FFT{i}.P1;
        DataInputY_Phase=unwrap(L_Compute_FFT{i}.P3);
        DataInputX=L_Compute_FFT{i}.f';

        % FIltro Resampling  Sabinsky - Golay Parametros
        InputX.DataInputX=DataInputX;
        InputX.fsData = length(DataInputX); % Sampling frequency
        InputX.fsResamp=InputX.fsData*3; %resampling frec
        InputX.filtOrder=25;
        InputX.filtPoints=499;        
        InputX.EnablePlot_FuncInt=0;   
%         InputX.xlim=VariablesProc.Figure_B.xlim;
%         InputX.ylim=VariablesProc.Figure_B.ylim;
        %TestMAResamp = lb5.ResampSabinGolay(DataInputY,InputX);
        
        
        % Calculo FIltro Resampling  Sabinsky - Golay          
        FFT_Filt_Mod=lb5.ResampSabinGolay(DataInputY_Mod,InputX); 
        FFT_Filt_Phase=lb5.ResampSabinGolay(DataInputY_Phase,InputX); 
%         % FFT amplitude magnitude
%         Dt.PlotData.Figure_B.P1(:,i)= FFT_Filt_Mod.Resamp_Y;
%         % FFT fase
%         Dt.PlotData.Figure_B.P3(:,i)=  FFT_Filt_Phase.Resamp_Y;  
        
        
            % Señal temporal
        k=k+1;
        ax1{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax1{i},DataInputX,DataInputY_Mod);
%         plot(ax1{i},FFT_Filt_Mod.Resamp_X,FFT_Filt_Mod.Resamp_Y);        
        xlabel(ax1{i},'f[Hz]')
        ylabel(ax1{i},'Amp.')
        title(ax1{i},DB.legendText{i})

        % Calculo el FFT   
        k=k+1;
        ax2{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax2{i},DataInputX,DataInputY_Phase);
%         plot(ax2{i},FFT_Filt_Phase.Resamp_X,FFT_Filt_Phase.Resamp_Y);
        
        xlabel(ax2{i},'f[Hz]');
        ylabel(ax2{i},'Fase [deg]');
        title(ax2{i},strcat(DB.legendText{i}));
    end
          
        DataInputY_Mod=L_Compute_FFT{1}.P1  ./ L_Compute_FFT{2}.P1(1:length(L_Compute_FFT{1}.P1),1);
        DataInputY_Phase=unwrap(L_Compute_FFT{1}.P3)-unwrap(L_Compute_FFT{2}.P3(1:length(L_Compute_FFT{1}.P3),1));
        DataInputX=L_Compute_FFT{1}.f';
        
        DIfIndex=3;
        SignalDiferencia(DIfIndex).DataInputY_Mod=DataInputY_Mod;
        SignalDiferencia(DIfIndex).DataInputY_Phase=DataInputY_Phase;
        SignalDiferencia(DIfIndex).DataInputX=DataInputX;

        k=k+1;
        ax1{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax1{i},DataInputX,DataInputY_Mod);
%         plot(ax1{i},FFT_Filt_Mod.Resamp_X,FFT_Filt_Mod.Resamp_Y);        
        xlabel(ax1{i},'f[Hz]')
        ylabel(ax1{i},'Amp.')
        title(ax1{i},DB.legendText{i})

        % Calculo el FFT   
        k=k+1;
        ax2{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax2{i},DataInputX,DataInputY_Phase);
%         plot(ax2{i},FFT_Filt_Phase.Resamp_X,FFT_Filt_Phase.Resamp_Y);  
        xlabel(ax2{i},'f[Hz]');
        ylabel(ax2{i},'Fase [deg]');
        title(ax2{i},strcat(DB.legendText{i}));       
end   






% TEst i
if(plotTestiData)         
    %----------------------
    FigureText='Ploting Test4 Diferencia FFT';
    %----------------------
    VectorToRead=[9 11]; 
    dhindex=dhindex+1;
    dh(dhindex)=figure;         
    set(dh(dhindex),'Name',FigureText);
    k=0;
    %VectorToRead=VectorToRead_Comp_1;
    %graficando raw data
    for i=VectorToRead    
        k=k+1;
        frec=DB.RawAudioData{i}.data_t;
        data_y=DB.RawAudioData{i}.data_y;
        L_Comp_FFT =lb5.Compute_FFT(frec,data_y);
%         N = length(data_y); 
%         dt = 1/L_Comp_FFT.Fs;
%         t = dt*(0:N-1)';
%         L_Comp_FFT.FrecCorrected=L_Comp_FFT.Fs*(0:N/2)/N;        
        L_Compute_FFT{k}=L_Comp_FFT;
    end
    
    k=0;
    for i=1:length(L_Compute_FFT)

        % los FFT data 
%         Dt.PlotData.Figure_B.f(:,i)=L_Compute_FFT{i}.f';
% %         Dt.PlotData.Figure_B.f(:,i)=L_Compute_FFT{i}.FrecCorrected;
%         % FFT amplitude magnitude
%         Dt.PlotData.Figure_B.P1(:,i)=L_Compute_FFT{i}.P1;
%         % FFT fase
%         Dt.PlotData.Figure_B.P3(:,i)=L_Compute_FFT{i}.P3;

            %--------------------------------
        % DATOS FIltro Resampling  Sabinsky - Golay
        DataInputY_Mod=L_Compute_FFT{i}.P1;
        DataInputY_Phase=unwrap(L_Compute_FFT{i}.P3);
        DataInputX=L_Compute_FFT{i}.f';

        % FIltro Resampling  Sabinsky - Golay Parametros
        InputX.DataInputX=DataInputX;
        InputX.fsData = length(DataInputX); % Sampling frequency
        InputX.fsResamp=InputX.fsData*3; %resampling frec
        InputX.filtOrder=25;
        InputX.filtPoints=499;        
        InputX.EnablePlot_FuncInt=0;   
%         InputX.xlim=VariablesProc.Figure_B.xlim;
%         InputX.ylim=VariablesProc.Figure_B.ylim;
        %TestMAResamp = lb5.ResampSabinGolay(DataInputY,InputX);
        
        
        % Calculo FIltro Resampling  Sabinsky - Golay          
        FFT_Filt_Mod=lb5.ResampSabinGolay(DataInputY_Mod,InputX); 
        FFT_Filt_Phase=lb5.ResampSabinGolay(DataInputY_Phase,InputX); 
%         % FFT amplitude magnitude
%         Dt.PlotData.Figure_B.P1(:,i)= FFT_Filt_Mod.Resamp_Y;
%         % FFT fase
%         Dt.PlotData.Figure_B.P3(:,i)=  FFT_Filt_Phase.Resamp_Y;  
        
        
            % Señal temporal
        k=k+1;
        ax1{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax1{i},DataInputX,DataInputY_Mod);
%         plot(ax1{i},FFT_Filt_Mod.Resamp_X,FFT_Filt_Mod.Resamp_Y);        
        xlabel(ax1{i},'f[Hz]')
        ylabel(ax1{i},'Amp.')
        title(ax1{i},DB.legendText{i})

        % Calculo el FFT   
        k=k+1;
        ax2{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax2{i},DataInputX,DataInputY_Phase);
%         plot(ax2{i},FFT_Filt_Phase.Resamp_X,FFT_Filt_Phase.Resamp_Y);
        
        xlabel(ax2{i},'f[Hz]');
        ylabel(ax2{i},'Fase [deg]');
        title(ax2{i},strcat(DB.legendText{i}));
    end
          
        DataInputY_Mod=L_Compute_FFT{1}.P1(1:length(L_Compute_FFT{2}.P1),1)  ./ L_Compute_FFT{2}.P1;
        DataInputY_Phase=unwrap(L_Compute_FFT{1}.P3(1:length(L_Compute_FFT{2}.P3),1))-unwrap(L_Compute_FFT{2}.P3);
        DataInputX=L_Compute_FFT{2}.f';
        
        DIfIndex=4;
        SignalDiferencia(DIfIndex).DataInputY_Mod=DataInputY_Mod;
        SignalDiferencia(DIfIndex).DataInputY_Phase=DataInputY_Phase;
        SignalDiferencia(DIfIndex).DataInputX=DataInputX;

        k=k+1;
        ax1{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax1{i},DataInputX,DataInputY_Mod);
%         plot(ax1{i},FFT_Filt_Mod.Resamp_X,FFT_Filt_Mod.Resamp_Y);        
        xlabel(ax1{i},'f[Hz]')
        ylabel(ax1{i},'Amp.')
        title(ax1{i},DB.legendText{i})

        % Calculo el FFT   
        k=k+1;
        ax2{i} = subplot(length(VectorToRead)+1,2,k); % top subplot
        plot(ax2{i},DataInputX,DataInputY_Phase);
%         plot(ax2{i},FFT_Filt_Phase.Resamp_X,FFT_Filt_Phase.Resamp_Y);  
        xlabel(ax2{i},'f[Hz]');
        ylabel(ax2{i},'Fase [deg]');
        title(ax2{i},strcat(DB.legendText{i}));       
end   





% TEst j
if(plotTestjData)  
    
    %----------------------
    FigureText='Ploting Test5 CompDiferencia FFT';
    %----------------------
 
    dhindex=dhindex+1;
    dh(dhindex)=figure;         
    set(dh(dhindex),'Name',FigureText);
    k=0;
     %     12 14 10 11
    n=12; 
    GrafTitle{1}=DB.legendText{n};
      n=14; 
    GrafTitle{2}=DB.legendText{n};
      n=10; 
    GrafTitle{3}=DB.legendText{n};
      n=11; 
    GrafTitle{4}=DB.legendText{n};
    
    
    for i=1:length(SignalDiferencia)
        
        DataInputY_Mod=SignalDiferencia(i).DataInputY_Mod;
        DataInputY_Phase=SignalDiferencia(i).DataInputY_Phase;
        DataInputX=SignalDiferencia(i).DataInputX *2.2/80;
                     k=k+1;
        ax1{i} = subplot(length(SignalDiferencia),2,k); % top subplot
        plot(ax1{i},DataInputX,DataInputY_Mod);
%         plot(ax1{i},FFT_Filt_Mod.Resamp_X,FFT_Filt_Mod.Resamp_Y);        
        xlabel(ax1{i},'f[Hz]')
        ylabel(ax1{i},'Amp.')
        title(ax1{i},DB.legendText{i})
        title(ax1{i},GrafTitle{i});
        grid on
        ylim([0 1000])

        % Calculo el FFT   
        k=k+1;
        ax2{i} = subplot(length(SignalDiferencia),2,k); % top subplot
        plot(ax2{i},DataInputX,DataInputY_Phase);
%         plot(ax2{i},FFT_Filt_Phase.Resamp_X,FFT_Filt_Phase.Resamp_Y);  
        xlabel(ax2{i},'f[Hz]');
        ylabel(ax2{i},'Fase [deg]');
        title(ax2{i},GrafTitle{i});
%         title(ax2{i},strcat(DB.legendText{i}));  % 
        grid on
        
    end
    

    
end

