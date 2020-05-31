% Libreria A: Rutinas para el fiteo de modelos de VOx.
% 
% HG 2016
%------------------------------------------------------------------------
% Historial
%----------
% Creación  de la libreria Agosto 2016
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Main
%------------------------------------------------------------------------

classdef LibAv05
 
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % COnstantes como propiedades (Properties — Storing Class Data)
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Const_Boltzmann_eV
        %-------------------------------------------------------------       
        properties(Constant)
            Const_Boltzmann_eV  = 8.6173324e-5; %8,6173324(78)×10-5 	eV?K?1
        end % end function   
 
        %-------------------------------------------------------------
        % Const_Kelvin
        %------------------------------------------------------------- 
        properties(Constant)
            Const_Kelvin  = 273.15; %Celsius 	[K] = [°C] + 273.15
        end     % end function   
        
        
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Metodos
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
        methods(Static)  

        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Rutinas Init
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function Init_______________________l() 
        end % end function    

        %-------------------------------------------------------------
        % loadInitVariables(titulo,autor,fecha)
        %-------------------------------------------------------------          
        function out=loadInitVariables(titulo,autor,fecha)
            disp('Cargo estructura de variables iniciales '); 
            
            %cargo estructura de variables iniciales
            out.titulo=titulo;
            out.autor=autor;
            out.fecha=fecha;
            
            %estructura de handler de imagenes
            out.ImageHandler=0;
         
            %CurrentFolder
            out.CurrentFolder=pwd;
            
            % lb4 figure exchange variable figure handler
            out.Figurehandler{1}=0;
                       
            
            % estructura de directorios
            Data='.\Output\Data';
            i=1;out.DirFile{i}=Data;
            out.DirFileData=Data;
            Working='.\Output\Working';
            out.DirFileTemp=Working;
            i=i+1;out.DirFile{i}=Working;
            Images='.\Output\Images'; 
            out.DirFileImages=Images;
            i=i+1;out.DirFile{i}=Images;
            ImagesPng='.\Output\ImagesPng';
            out.DirFileImagesPng=ImagesPng;
            i=i+1;out.DirFile{i}=ImagesPng;
            ImagesPdf='.\Output\ImagesPdf';
            out.DirFileImagesPdf=ImagesPdf;
            i=i+1;out.DirFile{i}=ImagesPdf;
            Temp='.\Output\Temp';              
            out.DirFileTemp=Temp;
            i=i+1;out.DirFile{i}=Temp;
            
            
        end % end function  
        
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Rutinas de Manejo de Archivos
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function Archivos_______________________l() 
        end % end function    

         %-------------------------------------------------------------
        % Rutinas Apertura Archivos de Datos
        %-------------------------------------------------------------
         function [RawData,filenameAtributes] =ArchiRetrieveData(DefaultFileName) 
                 
            [filenameAtributes.filename, ...
                filenameAtributes.pathname, ...
                filenameAtributes.filterindex] = uigetfile( ...
                {'*.mat','MAT-files (*.mat)'; ...
                '*.txt','txt-files (*.txt)'; ...    
                '*.mdl','Models (*.mdl)'; ...
                '*.*',  'All Files (*.*)'}, ...
                'Pick a file', ...
                'MultiSelect', 'off',...
                DefaultFileName);
            
            RawData=load(fullfile(filenameAtributes.pathname,...
                filenameAtributes.filename));
         
         end % end function    

        %-------------------------------------------------------------
        % Rutinas ArchiSaveFigure(handler,filename) 
        %-------------------------------------------------------------
         function ArchiSaveFigure(handler,filename) 
             
             saveas(handler,filename,'fig') ; 
                % guardando archivo png
             saveas(handler,filename,'png') ;
                        
         end  % end function   
         
         %-------------------------------------------------------------
        % Rutinas ArchiSaveFigurePdf(handler,filename) 
        %-------------------------------------------------------------
         function ArchiSaveFigurePdf(handler,VarInit,filename) 
             
             
             GeneralFontSize=14;
% %              matlabVersion=version ('-release');
% %              if matlabVersion~='2016a'
% %              end    
%                  %handler=gcf;
%                  %FigAxes=get(gcf,'CurrentAxes');
%                  FigAxes=get(handler,'CurrentAxes');             
%                  FigPaperPosition=get(handler,'PaperPosition');
%                  %FigOuterPosition=get(handler,'OuterPosition');  
%                  FigTightInset=get(FigAxes,'TightInset');  
% 
%                  ax_width = FigPaperPosition(3) - FigTightInset(1) - FigTightInset(3);
%                  ax_height = FigPaperPosition(4) - FigTightInset(2) - FigTightInset(4);
%     %             ax_width = FigPaperPosition(3) ;
%      %            ax_height = FigPaperPosition(4) ;
% 
%                  %set(handler,'PaperPositionMode','auto');
%                  %set(handler,'Position',FigOuterPosition);
%                  %set(handler,'PaperSize',[FigPaperPosition(3) FigPaperPosition(4)]);
%                  set(handler,'PaperSize',[ax_width ax_height]);      
            
            % recorto Imagen
            %output=handler;            
            FigurePosition=get(handler,'position');
            set(handler,'PaperPositionMode','auto');
            
            FigurePaperPosition=get(handler,'PaperPosition');
            FigureOuterPosition=get(handler,'OuterPosition');
            FigureCurrentAxes=get(handler,'CurrentAxes');
            FigureTightInset=get(FigureCurrentAxes,'TightInset');
            

            ax_width = FigurePaperPosition(3) - FigureTightInset(1) - FigureTightInset(3);
            ax_height = FigurePaperPosition(4) - FigureTightInset(2) - FigureTightInset(4);
            
            
            % seteo nuevo paperzize
            set(handler,'PaperSize',[ax_width ax_height]);
            
            % seteo nuevo ancho de linea
            %gca_local=get(handler,'CurrentAxes');
            %gca_childens=get(handler,'Children');
            %gca_allchildens=allchild(gca_childens);
            
            % tomo la cantidad de childrens
            %gca_childens=get(handler,'Children');
            % get el la el current axes de figura actual
%             get_CurrentAxes=get(handler,'CurrentAxes');
%             gca_childens=get(handler,'Children');
%             get_CurrentAxes_childens=get(get_CurrentAxes,'Children');


               %_-----------------------
               % paraincluir men el archivo
                %                GeneralFontSize=14;
                %             set(gca,'FontSize',GeneralFontSize);
                %             set(gca,'Title',text('String','' ));
                %             %get(gca,'Children');
                %             set(get(gca,'Children'),'LineWidth',2);
               
               %-----------------

            AllLineWidth=2;
            get_CurrentFigChilds=get(handler,'Children');
            get_CurrentAxes=get(handler,'CurrentAxes');
            get_CurrentAxes_childens=get(get_CurrentAxes,'Children');
            %exist('get_CurrentAxes_childens(i).LineWidth','class')
            % FOntSize
            
            set(get_CurrentAxes,'FontSize',GeneralFontSize);
            set(get_CurrentAxes,'Title',text('String','' ));
%             ylabel('Population','FontSize',12,...
%        'FontWeight','bold','Color','r')
            % asigno la dimensión a cada linea
            for i=1:length(get_CurrentAxes_childens)                
                get_CurrentAxes_childens(i).LineWidth=AllLineWidth; 
            end
            
               
                %gca_childens(i).LineWidth=AllLineWidth;
                %set(gca_childens(i).LineWidth,AllLineWidth);
                %gca_local=get(gca_childens,'CurrentAxes');
                %set(gca_local,'LineWidth',2);            
%             get_CurrentAxes=get(gcf,'CurrentAxes');
%             get_CurrentAxes_childens=get(get_CurrentAxes,'Children')
            %get_CurrentAxes_childens(3).LineWidth=3
            %gca_childens=get(gcf,'Children')
 
%             gca_childens=get(gcf,'Children')
%             gca_Child=get(gcf,'Child')
%             gca_childensCurrentAxes=(gca_childens(2),'CurrentAxes')
            
            
            %gca_childens=get(gcf,'DefaultLineLineWidth')            
            %gca_childensaxes=set(get(gcf,'DefaultLineLineWidth'),2)
%            gca_childensaxes=set(handler,'DefaultLineLineWidth',2)



            % save figure
             FigFilename=strcat(VarInit.DirFileImages,'\',filename);             
             saveas(handler,FigFilename,'fig') ; 
                % guardando archivo png
             FigFilename=strcat(VarInit.DirFileImagesPng,'\',filename);    
             saveas(handler,FigFilename,'png') ;
                 % guardando archivo pdf
             FigFilename=strcat(VarInit.DirFileImagesPdf,'\',filename);    
             saveas(handler,FigFilename,'pdf') ;                       
             
             
         end  % end function   
         
         
        %-------------------------------------------------------------
        % Rutinas ImportDataFromFileFolder
        % importo archivos .MAT o .CSV en Workspace
        %-------------------------------------------------------------
         function [out]=ImportDataFromFileFolder(DirFileNameToConvert) 
         
            % lectura de archivos csv o mat            
            [out.filenameAtributes.filename, ...
                out.filenameAtributes.pathname, ...
                out.filenameAtributes.filterindex] = uigetfile( ...
                {'*.csv','txt-files (*.csv)'; ... 
                '*.mat','MAT-files (*.mat)'; ...                
                '*.txt','txt-files (*.txt)'; ...    
                '*.mdl','Models (*.mdl)'; ...
                '*.*',  'All Files (*.*)'}, ...
                'Pick a file', ...
                'MultiSelect', 'off',...
                DirFileNameToConvert,...
                'MultiSelect','on');
            
            % recupero data en formato de dos columnas
            for i=1:length(out.filenameAtributes.filename)
            out.RawData{i}=importdata(fullfile(out.filenameAtributes.pathname,...
                out.filenameAtributes.filename{i}));            
            end
                        
         end  % end function 
         
         %-------------------------------------------------------------
        % Rutinas ImportDataAllFilesInFolder
        % importo archivos  en Workspace
        %-------------------------------------------------------------
         function [out]=ImportDataAllFilesInFolder(DirFileNameToConvert) 
                
             % imprimo mensaje
             LibAv05.DispCommWindwText(strcat('Importando datos:',DirFileNameToConvert));
             %a=mfilename;
             %b=depfun;
             %LibAv05.DispCommWindwText(a);
             % extraigo la lista de archivos a abrir
             [out.filenameAtributes.AllFilesDir]=dir(DirFileNameToConvert);
             j=1;
             for i=1:length(out.filenameAtributes.AllFilesDir)
                 if(not(isdir(out.filenameAtributes.AllFilesDir(i).name)))
                 FolderFilesNames{j}=out.filenameAtributes.AllFilesDir(i).name;
                 j=j+1;
                 end                 
             end
             
            % recupero data en formato de dos columnas
            for i=1:length(FolderFilesNames)
            out.RawData{i}=importdata(fullfile(DirFileNameToConvert,...
                FolderFilesNames{i}));  
            out.filenameAtributes.filename{i}=FolderFilesNames{i};
            end
                        
         end  % end function 
         
         
          %-------------------------------------------------------------
        % CreateWorkingFolder
        % COn esta rutina creo las subcarpetas de trabajo y de Salida
        %-------------------------------------------------------------
         function [out]=CreateWorkingAndOutputFolder(DirNameToCreate) 
                
             % imprimo mensaje
             LibAv05.DispCommWindwText(strcat('Crando el Directorio:',DirNameToCreate));
             
             % creo folder OutputFolderNK
            if not(exist(DirNameToCreate,'dir'))
                out.mk=mkdir(DirNameToCreate);
            else
                out.rm=rmdir(DirNameToCreate,'s');
                out.mk=mkdir(DirNameToCreate);
            end
                        
         end  % end function 
         
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Rutinas de FFT
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function FFT_______________________l() 
        end  % end function 
                
        %-------------------------------------------------------------
        % Rutinas Compute_FFT
        % Compute the FFT
        %-------------------------------------------------------------
        % Data corresponde a Raw data
        function [out]=Compute_FFT(RawData_X,RawData_Y) 
 
            out.X_data=RawData_X;
            out.Y_data=RawData_Y;

            out.Fs = length(out.X_data); % Sampling frequency
            out.T = 1/out.Fs; % Sampling period
            out.L = length(out.X_data);% Length of signal
            out.t = (0:out.L-1)*out.T;        % Time vector
    
            % %Compute the Fourier transform of the signal.
            % 
            out.X = (out.Y_data);
            out.Y = fft(out.X);
            % 
            % %Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
            % 
            out.P2 = abs(out.Y/out.L);
            out.P1 = out.P2(1:out.L/2+1);
            out.P1(2:end-1) = 2*out.P1(2:end-1);            
            % %Define the frequency domain f and plot the single-sided
            % amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. 
            % On average, longer signals produce better frequency approximations.
        
            out.f = out.Fs*(0:(out.L/2))/out.L;
            % %
            %Compute phase            
            % obtengo la fase y aplico Unwrap d ela fase
            %out.P4 = unwrap(angle(out.Y/out.L)); 
            out.P4 = angle(out.Y/out.L); 
            out.P3 = out.P4(1:out.L/2+1);
            out.P3(2:end-1) = 2*out.P3(2:end-1);

            %Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, because of the added noise. On average, longer signals produce better frequency approximations.

            %out.f = out.Fs*(0:(out.L/2))/out.L;
       
        end  % end function   
        
          
  
        %-------------------------------------------------------------
        % imchi=kkimbook2(omega,rechi,alpha)
        % Rutinas para obtener la parte Imaginaria por Kramers–Kronig
        % Kramers–Kronig Relations in Optical Materials Research
        % V. Lucarini J.J. Saarinen K.-E. Peiponen E.M. Vartiainen
        % 2005
        %-------------------------------------------------------------
        
        function imchi=kkimbook2(omega,rechi,alpha)
            %The program inputs are the vector of the frequency (or energy)
            %components, the vector of the real part of the susceptibility
            %under examination, and the value of the moment considered.
            %The two vectors must have the same length 
            %and the frequency vector omega must be equispaced. 
            %If not, apply MATLAB functions such as interp.
            %If rechi is the real part of a linear susceptibility, 
            %alpha must be 0. 
            %If rechi is the real part of the nth 
            %harmonic generation susceptibility, alpha=0,1,..2n. 
            %If rechi is the real part of a pump and probe
            %susceptibility, alpha=0 or 1.
            %This files accompanies the book 
            %"Kramers-Kronig Relations in Optical Materials Research"
            %by Lucarini, V., Saarinen, J.J., Peiponen, K.-E., Vartiainen, E.M. 
            %Springer, Heidelberg, 2005
            %where the theory and applications are fully developed.
            %The output is the estimate of the imaginary part as obtained
            %with K-K relations.
            %This software is distributed under the GNU licence agreement
            %by Valerio Lucarini
            %email: lucarini@alum.mit.edu
            %University of Camerino
            %Department of Mathematics and Computer Science
            %Camerino, Italy

            if size(omega,1)>size(omega,2);
            omega=omega';
            end; 
            if size(rechi,1)>size(rechi,2);
            rechi=rechi';
            end;
            %Here the program rearranges the two vectors so that,
            %whichever their initial shape, they become row vectors.

            g=size(omega,2);
            %Size of the vectors.%

            imchi=zeros(size(rechi));
            %The output is initialized.

            a=zeros(size(rechi));
            b=zeros(size(rechi));
            %Two vectors for intermediate calculations are initialized

            deltaomega=omega(2)-omega(1);
            %Here we compute the frequency (or energy) interval

            j=1;
            beta1=0;
            for k=2:g;
            b(1)=beta1+rechi(k)*omega(k)^(2*alpha)/(omega(k)^2-omega(1)^2);
            beta1=b(1);
            end;
            imchi(1)=-2/pi*deltaomega*b(1)*omega(1)^(1-2*alpha);
            %First element of the output: the principal part integration
            %is computed by excluding the first element of the input

            j=g;
            alpha1=0;
            for k=1:g-1;
            a(g)=alpha1+rechi(k)*omega(k)^(2*alpha)/(omega(k)^2-omega(g)^2);
            alpha1=a(g);
            end;
            imchi(g)=-2/pi*deltaomega*a(g)*omega(g)^(1-2*alpha);
            %Last element of the output: the principal part integration
            %is computed by excluding the last element of the input.

            for j=2:g-1; 
            %Loop on the inner components of the output vector.
            alpha1=0;
            beta1=0;
            for k=1:j-1;
            a(j)=alpha1+rechi(k)*omega(k)^(2*alpha)/(omega(k)^2-omega(j)^2);
            alpha1=a(j);
            end;
            for k=j+1:g;
            b(j)=beta1+rechi(k)*omega(k)^(2*alpha)/(omega(k)^2-omega(j)^2);
            beta1=b(j);
            end;
            imchi(j)=-2/pi*deltaomega*(a(j)+b(j))*omega(j)^(1-2*alpha);
            %Last element of the output: the principal part integration
            %is computed by excluding the last element of the input
        end;
        
        end;
        
        
        %-------------------------------------------------------------
        % imchi=kkrebook2(omega,rechi,alpha)
        % Rutinas para obtener la parte Imaginaria por Kramers–Kronig
        % Kramers–Kronig Relations in Optical Materials Research
        % V. Lucarini J.J. Saarinen K.-E. Peiponen E.M. Vartiainen
        % 2005
        %-------------------------------------------------------------
        
        function rechi=kkrebook2(omega,imchi,alpha)
            %The program inputs are the vector of the frequency
            %(or energy) components, the vector of the imaginary
            %part of the susceptibility under examination, and
            %the value of the moment considered. 
            %The two vectors must have the same length 
            %and the frequency vector omega must be equispaced. 
            %If not, apply MATLAB functions such as interp.
            %If imchi is the imaginary part of a linear susceptibility, 
            %alpha must be 0. 
            %If imchi is the imaginary part of the nth 
            %harmonic generation susceptibility, alpha=0,1,..2n. 
            %If imchi is the imaginary part of a pump and probe
            %susceptibility, alpha=0 or 1.
            %This files accompanies the book 
            %"Kramers-Kronig Relations in Optical Materials Research"
            %by Lucarini, V., Saarinen, J.J., Peiponen, K.-E., Vartiainen, E.M. 
            %Springer, Heidelberg, 2005
            %where the theory and applications are fully developed.
            %The output is the estimate of the real part as obtained
            %with K-K relations.
            %This software is distributed under the GNU licence agreement
            %by Valerio Lucarini
            %email: lucarini@alum.mit.edu
            %University of Camerino
            %Department of Mathematics and Computer Science
            %Camerino, Italy

            if size(omega,1)>size(omega,2);
            omega=omega';
            end; 
            if size(imchi,1)>size(imchi,2);
            imchi=imchi';
            end;
            %Here the program rearranges the two vectors so that,
            %whichever their initial shape, they become row vectors.
            g=size(omega,2);
            %Size of the vectors.%
            rechi=zeros(size(imchi));
            %The output is initialized.
            a=zeros(size(imchi));
            b=zeros(size(imchi));
            %Two vectors for intermediate calculations are initialized
            deltaomega=omega(2)-omega(1);
            %Here we compute the frequency (or energy) interval
            j=1;
            beta1=0;
            for k=2:g;
            b(1)=beta1+imchi(k)*omega(k)^(2*alpha+1)/(omega(k)^2-omega(1)^2);
            beta1=b(1);
            end;
            rechi(1)=2/pi*deltaomega*b(1)*omega(1)^(-2*alpha);
            %First element of the output: the principal part integration
            %is computed by excluding the first element of the input
            j=g; 
            alpha1=0; 
            for k=1:g-1;
            a(g)=alpha1+imchi(k)*omega(k)^(2*alpha+1)/(omega(k)^2-omega(g)^2);
            alpha1=a(g);
            end;
            rechi(g)=2/pi*deltaomega*a(g)*omega(g)^(-2*alpha);
            %Last element of the output: the principal part integration
            %is computed by excluding the last element of the input
            for j=2:g-1; 
            %Loop on the inner components of the output vector.
            alpha1=0;
            beta1=0;
            for k=1:j-1;
            a(j)=alpha1+imchi(k)*omega(k)^(2*alpha+1)/(omega(k)^2-omega(j)^2);
            alpha1=a(j);
            end;
            for k=j+1:g;
            b(j)=beta1+imchi(k)*omega(k)^(2*alpha+1)/(omega(k)^2-omega(j)^2);
            beta1=b(j);
            end;
            rechi(j)=2/pi*deltaomega*(a(j)+b(j))*omega(j)^(-2*alpha);
            
            end;
            %Last element of the output: the principal part integration
            %is computed by excluding the last element of
            
        end
        
        
        
        

        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Rutinas de Ploteo
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function Plot_______________________l() 
        end  % end function                
  
        %-------------------------------------------------------------
        %Rutinas de  PlotStantard
        % Ploteo basico de X-Y
        %-------------------------------------------------------------
        function [handler]=PlotStantard(PlotST_X,PlotST_Y,PlotST_Title...
                           ,PlotST_xlabel,PlotST_ylabel,PlotSTColorSimb,PlotSTType)
            %default values
            %PlotSTType: tipo de grafigo :log en eje Y o linear        
            if ~exist('PlotST_Title','var'), PlotST_Title = ['Plot Title:',inputname(2),' vs ',inputname(1)]; end
            if ~exist('PlotST_xlabel','var'), PlotST_xlabel = ['X= ',inputname(1)]; end
            if ~exist('PlotST_ylabel','var'), PlotST_ylabel = ['Y= ',inputname(2)]; end
            if ~exist('PlotSTColorSimb','var'), PlotSTColorSimb = 'b'; end
            if ~exist('PlotSTType.YScale','var'), PlotSTType.YScale = 'linear'; end
            % main
            handler=plot(PlotST_X,PlotST_Y,PlotSTColorSimb); 
            set(gca,'Title',text('String',PlotST_Title,'Color',PlotSTColorSimb));
            set(gca,'xlabel',text('String',PlotST_xlabel)) ;
            set(gca,'ylabel',text('String',PlotST_ylabel)) ;
            %set(gca,'YScale',text('String',PlotSTType.YScale)) ;
            set(gca,'YScale',PlotSTType.YScale) ;
            
            disp(['PlotStantard:','Ploting',inputname(2),' vs ',inputname(1) ]);
            
        end% end function   
         
        %-------------------------------------------------------------
        % Rutinas de Plot_AnotateByGetPosition
        %-------------------------------------------------------------
        function output_txt =PlotAnotateByGetPosition(Axehandler,txtLabelArray,posX,posY,Properties)   
        %function  myfunction(obj,event_obj)
        % Display the position of the data cursor
        % obj          Currently not used (empty)
        % event_obj    Handle to event object
        % output_txt   Data cursor text string (string or cell array of strings).

        %PlotSTType: tipo de grafigo :log en eje Y o linear
        if ~exist('Properties','var'), Properties.FontSize = '7'; end
            
        % handlerds
        
        set(Axehandler,'Units','normalized')
        p1 = get(Axehandler,'Position');
        %[left bottom width height]
        
        posx00=p1(1);
        posy00=p1(2);
        width00=p1(3);
        height00=p1(4);
        
        getCurrentAxexlim=get(Axehandler,'xlim');
        getCurrentAxeylim=get(Axehandler,'ylim');
        
        % anotate
        for i=1 : length(posX) 
                
                CoodY=(posY(i)-getCurrentAxeylim(1))/...
                            (getCurrentAxeylim(2)-getCurrentAxeylim(1));
                        
                CoodX=(posX(i)-getCurrentAxexlim(1))/...
                            (getCurrentAxexlim(2)-getCurrentAxexlim(1));
                                        
                NormalizedCoodY=  CoodY*height00; 
                NormalizedCoodX=  CoodX*width00;
                
                txtToVis =txtLabelArray(i)  ;     
                % Create the textarrow object: 
                output_txt =  annotation('textbox',[posx00+NormalizedCoodX,...
                                                    posy00+NormalizedCoodY,...
                                                    .07 .03],...
                                    'String',txtToVis,...
                                    'FontSize',9,...
                                    'FitBoxToText','on',...
                                    'Interpreter','none',...
                                    'LineStyle','none');
               
        end % ned for

        
        end  % end function     

        %-------------------------------------------------------------
        % PlotComputeMeyerNeldelRule
        %-------------------------------------------------------------
        % Calcula, fitea y plot del ser de datos aportados
        % parametros:   
        %               p1: Datos x
        %               p1: Datos y
        %               p1: string labels
        function output =PlotComputeMeyerNeldelRule(LocalDataX,LocalDataY,txtLabel,PlotFit,MNR_Properties)   
            %default values
            %PlotSTType: tipo de grafigo :log en eje Y o linear
            if ~exist('PlotFit','var'), PlotFit = 'no'; end    
            if ~exist('MNR_Properties','var'), MNR_Properties.FontSize = '7'; end    
            
            %plot
            hold on
            %output.Fig_handler=figure;
            output.PlotST=LibAv05.PlotStantard(LocalDataX,LocalDataY,'s');
            %anotate labels
            LibAv05.PlotAnotateByGetPosition(gca,txtLabel,LocalDataX,LocalDataY,MNR_Properties) ;
     
            if strcmp(PlotFit,'yes')
                % calculo el fiteo
                fitType=fittype('poly1');
                %ExcludedData=excludedata(Dt.MNR.Data_X,Dt.MNR.Data_Y,'indices',Dt.MNR.indiceExclude);
                [output.fitcfun,output.fitgof,output.fitoutput] =fit(LocalDataX',...
                                                        LocalDataY',fitType);
%                                                     ,...
%                                                         'Exclude',Dt.MNR.ExcludedData 
                output.fitformula=formula(output.fitcfun);
                output.fitcoeffvalues=coeffvalues(output.fitcfun);
                output.fitrsquare=output.fitgof.rsquare;
                
                    
                output.PlotFit=plot(output.fitcfun,'r')    ;  
                [output.legend_h,output.legend_strings]=legend...
                        ('Data',strcat('(fit) y=',output.fitformula,' ',...
                        ' - Coef=[',num2str(output.fitcoeffvalues),'], ',...
                        ' (R^2=',num2str(output.fitrsquare),')'),...
                        'Location','Best');
                set(output.legend_h, 'FontSize',8, 'Color', 'none');
                
            end
            
            % fin hold off
            hold off
                        
        end  % end function  
        
        
        %-------------------------------------------------------------
        % PlotingIFG
        %-------------------------------------------------------------
        
        % parametros:   
        %               Input 1:  
        function output =PlotingIFG(Input1,Input2,Input3,Input4)   
        
            
        subplot(3,1,1);        
        output.Sub1 =plot(Input1.x,Input1.y);
        %title('Interferogram')
        title(strcat('Interferogram (',Input2.filenameAtributes.filename{Input4},num2str(Input4),')'));
        xlabel('t ');
        ylabel('X(t)'); 
        legend(Input1.legend,...
                    'FontSize',Input3.FondSize,...
                    'Location','best');  
        %xlim(Input3.xlim);
        %ylim(Input3.ylim); 
        
        
        subplot(3,1,2);
        output.Sub2 =plot(Input1.f,Input1.P1);
        title('Single-Sided Amplitude Spectrum of X(t)');
        xlabel('f (cm^{-1})');
        ylabel('|P1(f)|');
        legend(Input1.legend,...
                    'FontSize',Input3.FondSize,...
                    'Location','best'); 
        xlim(Input3.xlim);
        ylim(Input3.ylim);        
        
        subplot(3,1,3);
        output.Sub3 =plot(Input1.f,Input1.P3);
        title('Single-Sided Phase Spectrum of X(t)');
        xlabel('f (cm^{-1})');
        ylabel('\phi(f)[Rad]');
        legend(Input1.legend,...
                    'FontSize',Input3.FondSize,...
                    'Location','best'); 
        xlim(Input3.xlim);
        ylim(Input3.ylim); 
        
        end
        
        %-------------------------------------------------------------
        % PlotingIFG
        %-------------------------------------------------------------
        
        % parametros:   
        %               Input 1:  
        function output =Ploting_nk(Input1,Input2,Input3,Input4)   
        
            
         
            subplot(3,1,1);
            output.Sub1 =plot(Input1.Nhu,Input1.ReflectanceEnergy);
            title('Single-Sided Reflectance Energy Spectrum ');
            xlabel('f (cm^{-1})');
            ylabel('Reflectance Energy');
            legend(Input1.legend,...
                        'FontSize',Input3.FondSize,...
                        'Location','best');                          
            xlim(Input3.xlim);
            ylim(Input3.ylim); 
        
            subplot(3,1,2);
            output.Sub1 =plot(Input1.Nhu,Input1.RefractiveIndex_n);
            title('Single-Sided Refractive Index n Spectrum ');
            set(gca,'XDir','reverse');
            xlabel('f (cm^{-1})');
            ylabel('Refractive Index (n)');
            legend(Input1.legend,...
                        'FontSize',Input3.FondSize,...
                        'Location','best'); 
            xlim(Input3.xlim);
            ylim(Input3.ylim); 

            subplot(3,1,3);
            output.Sub2 =plot(Input1.Nhu,Input1.AbsortionCoefficien_k);
            title('Single-Sided Absortion Coefficient k Spectrum');
            xlabel('f (cm^{-1})');
            ylabel('Absortion Coefficien (k)');
            legend(Input1.legend,...
                        'FontSize',Input3.FondSize,...
                        'Location','best'); 
            xlim(Input3.xlim);
            ylim(Input3.ylim);                
        end
        
        
        %-------------------------------------------------------------
        % Ploting_nk_KramerKroneg
        %-------------------------------------------------------------
        
        % parametros:   
        %               Input 1:  
        function output =Ploting_nk_KramerKroneg(Input1,Input2,Input3,Input4)   
        
            
         
            subplot(3,1,1);
            output.Sub1 =plot(Input1.Nhu,Input1.ReflectanceEnergy);
            title('Single-Sided Reflectance Energy Spectrum ');
            xlabel('f (cm^{-1})');
            ylabel('Reflectance Energy');
            legend(Input1.legend,...
                        'FontSize',Input3.FondSize,...
                        'Location','best');        

            subplot(3,1,2);
            output.Sub1 =plot(Input1.Nhu,Input1.RefractiveIndex_n);
            title('Single-Sided Refractive Index n Spectrum ');
            xlabel('f (cm^{-1})');
            ylabel('Refractive Index (n)');
            legend(Input1.legend,...
                        'FontSize',Input3.FondSize,...
                        'Location','best'); 


            subplot(3,1,3);
            output.Sub2 =plot(Input1.Nhu,Input1.KramerKronegImchi);
            title('Single-Sided Absortion Coefficient k Spectrum');
            xlabel('f (cm^{-1})');
            ylabel('Absortion Coefficien (k)');
            legend(Input1.legend,...
                        'FontSize',Input3.FondSize,...
                        'Location','best'); 
                
        end
        
        
        %-------------------------------------------------------------
        % Ploting n - K por Kramer-Kroneg
        %-------------------------------------------------------------
        
        % parametros:   
        %               Input 1:  
        function output =Ploting_nk_KK(Input1,Input2,Input3,Input4,Input5)   
        
            if ~exist('Input5','var'), Input5.ReversePlot = 0; end
         
            subplot(3,1,1);
            output.Sub1 =plot(Input1.Nhu,Input1.ReflectanceEnergy);
            title('Single-Sided Reflectance Energy Spectrum ');
            % eje X invertido
            if (Input5.ReversePlot), set(gca,'XDir','reverse'); end
            xlabel('f (cm^{-1})');
            ylabel('Reflectance Energy');
            legend(Input1.legend,...
                        'FontSize',Input3.FondSize,...
                        'Location','best');                          
            xlim(Input3.xlim);
            ylim(Input3.ylim); 
        
            subplot(3,1,2);
            output.Sub1 =plot(Input1.Nhu,Input1.RefractiveIndex_n);
            title('Single-Sided Refractive Index n Spectrum ');
            % eje X invertido
            if (Input5.ReversePlot), set(gca,'XDir','reverse'); end
            xlabel('f (cm^{-1})');
            ylabel('Refractive Index (n)');
            legend(Input1.legend,...
                        'FontSize',Input3.FondSize,...
                        'Location','best'); 
            xlim(Input3.xlim);
            ylim(Input3.ylim); 

            subplot(3,1,3);
            output.Sub2 =plot(Input1.Nhu,Input1.AbsortionCoefficien_k);
            % eje X invertido
            if (Input5.ReversePlot), set(gca,'XDir','reverse'); end
            title('Single-Sided Absortion Coefficient k Spectrum');
            xlabel('f (cm^{-1})');
            ylabel('Absortion Coefficien (k)');
            legend(Input1.legend,...
                        'FontSize',Input3.FondSize,...
                        'Location','best'); 
            xlim(Input3.xlim);
            ylim(Input3.ylim);                
        end
        
        
        %-------------------------------------------------------------
        % Ploting n - K por Kramer-Kroneg
        %-------------------------------------------------------------
        
        % parametros:   
        %               Input 1:  
        function output =Ploting_nkOnly_KK(Input1,Input2,Input3,Input4,Input5,Input6,Input7)   
        
            if ~exist('Input5','var'), Input5.ReversePlot = 0; end;
            
            %figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])       
         
            iV=Input4.initVector;
            eV=Input4.endVector;
            
            output = figure;
            GeneralFontSize=10;
            %------------------
            % ploting refractive index n
            subplot(2,1,1);
            hold on
            % plot dots
            
            plot(Input1.Nhu(:,iV:eV),Input1.RefractiveIndex_n(:,iV:eV),'.');
            % plot fit
            plot(Input2.DataXPeakfit(:,iV:eV),Input2.DataYbaseToAdd(:,iV:eV),'-','LineWidth', 2);
            
            %title('Espectro Indice de refracción (n)');
            % eje X invertido
            if (Input5.ReversePlot), set(gca,'XDir','reverse'); end
            xlabel('f (cm^{-1})');
            ylabel('Indice de refracción (n)');
            legend(Input1.legend(:,iV:eV),...
                        'FontSize',Input3.FondSize,...
                        'Location','northeast'); 
                              %_-----------------------
               % paraincluir men el archivo
                            %GeneralFontSize=12;
                            set(gca,'FontSize',GeneralFontSize);
                %             set(gca,'Title',text('String','' ));
                            %get(gca,'Children');
                            set(get(gca,'Children'),'LineWidth',2);
               
               %-----------------
            xlim(Input7.xlim);
            ylim(Input3.ylim); 
            hold off
            %------------------
            % ploting extinction coefficient  k
            subplot(2,1,2);
            hold on
            plot(Input1.Nhu(:,iV:eV),Input1.AbsortionCoefficien_k(:,iV:eV),'.');
            % plot fit
            plot(Input6.DataXPeakfit(:,iV:eV),Input6.DataYbaseToAdd(:,iV:eV),'-','LineWidth', 2);
            % eje X invertido
            if (Input5.ReversePlot), set(gca,'XDir','reverse'); end
            %title('Espectro Coeficiente de absorcion (k)');
                               %_-----------------------
               % paraincluir men el archivo
                            %GeneralFontSize=12;
                            set(gca,'FontSize',GeneralFontSize);
                %             set(gca,'Title',text('String','' ));
                            %get(gca,'Children');
                            set(get(gca,'Children'),'LineWidth',2);
               
               %-----------------
            xlabel('f (cm^{-1})');
            ylabel('Coef. de absorción (k)');
            legend(Input1.legend(:,iV:eV),...
                        'FontSize',Input3.FondSize,...
                        'Location','northeast'); 
            xlim(Input7.xlim);
            ylim(Input7.ylim);        
            hold off
        end
        
        
        
        %-------------------------------------------------------------
        % Ploting_with_Legend
        %-------------------------------------------------------------
        
        % parametros:   
        %               Input 1:  
        function output =Ploting_with_Legend(Input1,Input2)   
                       
        output.Fig1 =plot(Input1.dataX,Input1.dataY);
        title(Input2.Name);
        %xlabel('f (cm^{-1})');
        %ylabel('Refractive Index (n)');
        legend(Input1.legend,...
                    'FontSize',Input2.FondSize,...
                    'Location','best'); 
        %xlim(Input2.xlim);
        
                
        end
        

        
%         % The function plotting figure inside figure (main and inset) from 2 existing figures.
%         % inset_size is the fraction of inset-figure size, default value is 0.35
%         % The outputs are the axes-handles of both.
%         % 
%         % An examle can found in the file: inset_example.m
%         % 
%         % Moshe Lindner, August 2010 (C).
%         
%         function [h_main, h_inset]=Plot_inset(result_handler,...
%                             main_handle, inset_handle,inset_size)
% 
%         if nargin==3
%             inset_size=0.35;
%         end
% 
%         inset_size=inset_size*.7;
%         %figure
%         new_fig=result_handler;%gcf;
%         main_fig = findobj(main_handle,'Type','axes');
%         h_main = copyobj(main_fig,new_fig);
%         set(h_main,'Position',get(main_fig,'Position'))
%         inset_fig = findobj(inset_handle,'Type','axes');
%         h_inset = copyobj(inset_fig,new_fig);
%         ax=get(main_fig,'Position');
%         set(h_inset,'Position', [.7*ax(1)+ax(3)-inset_size .5*ax(2)+ax(4)-inset_size inset_size inset_size])
% 
%         end

        %-------------------------------------------------------------
        % Ploting_Plot_inset
        %-------------------------------------------------------------
        % The function plotting figure inside figure (main and inset) from 2 existing figures.
        % inset_size is the fraction of inset-figure size, default value is 0.35
        % The outputs are the axes-handles of both.
        % 
        % An examle can found in the file: inset_example.m
        % 
        % Moshe Lindner, August 2010 (C).        
        
        function [new_fig_h, h_main, h_inset]=Plot_inset(result_handle,...
                main_handle, inset_handle,inset_size)

            if nargin==2
                inset_size=0.35;
            end

            inset_size=inset_size*.7;
            %figure
            new_fig_h=result_handle;
            main_fig = findobj(main_handle,'Type','axes');
            h_main = copyobj(main_fig,new_fig_h);
            get_main_fig_position=get(main_fig,'Position');
            set(h_main,'Position',get_main_fig_position)
            inset_fig = findobj(inset_handle,'Type','axes');
            h_inset = copyobj(inset_fig,new_fig_h);
            get_main_fig_position2=get(main_fig,'Position');
            ax=get_main_fig_position2;
            set(h_inset,'Position', [.7*ax(1)+ax(3)-inset_size .5*ax(2)+...
                ax(4)-inset_size inset_size inset_size])
        end

        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Rutinas de Dialogo e ingreso de datos
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function Diag_______________________l() 
        end % end function   
    
        
        %-----------------------------------------------------
        % DiagInputBox
        %-----------------------------------------------------
        function OutString=DiagInputBox(tret) 
            

            %prompt={'Enter the matrix size for x^2:',...
            %        'Enter the colormap name:'};

            InputBox.prompt={'name',...
                    'type',...
                    'xlimitesall',...
                    'xlimitesMin',...
                    'xlimitesMax',...
                    'DBFiles.dir_files',...
                    'DBFiles.regExprSubDataFile',...
                    'PlotCurvas.plotAllIndex',...
                    'PlotCurvas.IndexINIT',...
                    'PlotCurvas.IndexEND'};    


            InputBox.name='InputBox';
            InputBox.numlines=1;
            InputBox.defaultanswer={
            RegExpFit(RegExpFitindex).name,...
            num2str(RegExpFit(RegExpFitindex).type) ,...
            num2str(RegExpFit(RegExpFitindex).xlimitesall),...
            num2str(RegExpFit(RegExpFitindex).xlimites(1)) ,...
            num2str(RegExpFit(RegExpFitindex).xlimites(2)) ,...
            num2str(RegExpFit(RegExpFitindex).DBFiles.dir_files),...
            num2str(RegExpFit(RegExpFitindex).DBFiles.regExprSubDataFile),...
            num2str(RegExpFit(RegExpFitindex).PlotCurvas.plotAllIndex),...
            num2str(RegExpFit(RegExpFitindex).PlotCurvas.IndexINIT),...
            num2str(RegExpFit(RegExpFitindex).PlotCurvas.IndexEND)};

            InputBox.options.Resize='on';
            InputBox.options.WindowStyle='normal';
            InputBox.options.Interpreter='tex';
            %[answer ]=inputdlg(prompt,name,numlines,defaultanswer);
            InputBox.answer=inputdlg(InputBox.prompt,InputBox.name,InputBox.numlines,...
                        InputBox.defaultanswer,InputBox.options);


            if  isempty(InputBox.answer)
                %break;
            else
            RegExpFit(RegExpFitindex).name=char(InputBox.answer(1));
            RegExpFit(RegExpFitindex).type=str2double(InputBox.answer(2));
            RegExpFit(RegExpFitindex).xlimitesall=str2double(InputBox.answer(3));
            RegExpFit(RegExpFitindex).xlimites(1)=str2double(InputBox.answer(4));
            RegExpFit(RegExpFitindex).xlimites(2)=str2double(InputBox.answer(5));
            RegExpFit(RegExpFitindex).DBFiles.dir_files=str2double(InputBox.answer(6));
            RegExpFit(RegExpFitindex).DBFiles.regExprSubDataFile=str2double(InputBox.answer(7));
            RegExpFit(RegExpFitindex).PlotCurvas.plotAllIndex=str2double(InputBox.answer(8));
            RegExpFit(RegExpFitindex).PlotCurvas.IndexINIT=str2double(InputBox.answer(9));
            RegExpFit(RegExpFitindex).PlotCurvas.IndexEND=str2double(InputBox.answer(10));
            end

            
        end % end function   
        
    
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Rutinas de Visualización y Display
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function Disp_______________________l() 
        end % end function   
        
        %-------------------------------------------------------------
        % LibVersion
        %-------------------------------------------------------------
        function message = DispLibVersion()
            message=mfilename('fullpath');
            %disp(['Using Lib Vers.: '  message ]);
            LibAv05.DispCommWindwTitle(['Using Lib: ' message])
        end % end function   

        %-------------------------------------------------------------
        % DispCommWindwText : visualizo mensaje comun
        %-------------------------------------------------------------          
        function DispCommWindwText(text) 
            disp(text );            
        end % end function   
        
        %-------------------------------------------------------------
        % DispCommWindwMess : visualizo mensaje comun
        %-------------------------------------------------------------          
        function DispCommWindwMess(text)
            disp(' ');
            disp(text );            
        end % end function   
         
        %-------------------------------------------------------------
        % DispCommWindwTitle: VIsualizo Titulo
        %-------------------------------------------------------------          
        function DispCommWindwTitle(text)% Disp
            disp(' ');
            disp('@----------------------------------------------------------------------');
            disp(['@ ', text  ]);
            disp('@----------------------------------------------------------------------');
        end % end function   
        
        %-------------------------------------------------------------
        % DispCommWindwPortada
        %-------------------------------------------------------------          
        function DispCommWindwPortada(titulo,autor,fecha)
            disp(' ');
            disp('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
            disp('@----------------------------------------------------------------------');
            disp(['@ ' titulo]);
            disp(['@ ' autor]);
            disp(['@ ' fecha]);
            disp('@----------------------------------------------------------------------');
            disp('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@');
             
        end % end function   

        %-------------------------------------------------------------
        % DispCommWindwScriptEndElapced
        %-------------------------------------------------------------          
        function DispCommWindwScriptEndElapced(Scriptinit,Scriptend) 
            
            strDateTimeInit=LibAv05.ConvFecha(Scriptinit);
            strDateTimeEnd=LibAv05.ConvFecha(Scriptend);
            strDateTimeElapced=LibAv05.ConvFecha(Scriptend-Scriptinit);
            disp(' ');                   
            disp('@----------------------------------------------------------------------');
            disp('@ Ending Cript:' );
            disp(['@ INIT:' strDateTimeInit]);
            disp(['@ END:' strDateTimeEnd]); 
            disp(['@ ELAP:' strDateTimeElapced]); 
            disp('@----------------------------------------------------------------------'); 
        end % end function   
        
       
        
        
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Rutinas de Filtrado
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function Filtering_______________________l() 
        end  % end function          
        
        %-------------------------------------------------------------
        % FIltro Resampling  Sabinsky - Golay
        %------------------------------------------------------------- 
        % DataY = dato a resalplear
        % fsData = Frec de la señal  
        % fsResamp= Frecuencia de resampleo
        % filtOrder= orden del filtro FIR Sabinsky - Golay
        % filtPoints=cantidad de puntos del filtro . There are roughly 1000 / 60 =
        % 16.667 samples in a complete cycle of 60 Hz when sampled at
        %1000 Hz. Let's attempt to "round up" and use a 17-point filter. 
        %This will give us maximal filtering at a fundamental frequency 
        %of 1000 Hz / 17 = 58.82 Hz. There are roughly 1000 / 60 = 16.667 samples 
        %in a complete cycle of 60 Hz when sampled at 1000 Hz. 
        %%Let's attempt to "round up" and use a 17-point filter. 
        %This will give us maximal filtering at a fundamental frequency 
        %of 1000 Hz / 17 = 58.82 Hz. 
        function [out]=ResampSabinGolay(DataY,InputX) 
            
             vResamp = resample(DataY, InputX.fsResamp, InputX.fsData);
             %out.Resamp_X = ((0:numel(vResamp)-1) / InputX.fsResamp)';
             out.Resamp_Y= sgolayfilt(vResamp,InputX.filtOrder,InputX.filtPoints);
             out.Resamp_Y= resample(out.Resamp_Y, InputX.fsData,InputX.fsResamp);
             out.Resamp_X= InputX.DataInputX;
             % plot results
             if InputX.EnablePlot_FuncInt
                out.Plot_h=figure;
                plot(InputX.DataInputX, [DataY out.Resamp_Y]);
                legend('DataInputY','DataInputY Filtered and resampled' ,'location','best');
                xlim(InputX.xlim);
                ylim(InputX.ylim);
             end    
            
        end
        
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Rutinas de Manejo y conversiónDatos
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function Conv_______________________l() 
        end  % end function   
        
        %-------------------------------------------------------------
        % Rutinas de Manejo y conversiónDatos de fecha de Fix Now a str
        %------------------------------------------------------------- 
        function strDateTime=ConvFecha(DateTimeNow) 
            a=fix(DateTimeNow);
            strDateTime=['(' num2str(a(3)) '/' ...
                    num2str(a(2)) '/'...
                    num2str(a(1)) ')' ...
                    ' - ' num2str(a(4),'%10.0f') ':' ...
                    num2str(a(5),'%10.0f') ':' ...
                    num2str(a(6),'%10.0f')];
        end % end function   
       
        %-------------------------------------------------------------
        % Rutinas de conversión de temperatura Centigrado en Kelvin
        %------------------------------------------------------------- 
        function dataOut=ConvCentigrad2Kelvin(dataIn) 
            dataOut=273.15+dataIn;                
        end % fin de rutina función
        
        
        
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Rutinas Para dibidir o quebrar  graficos
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function breakyaxes_______________________l() 
        end  % end function   

        
        
        % breakyaxes splits data in an axes so that data is in a low and high pane.
        %
        %   breakYAxes(splitYLim) splitYLim is a 2 element vector containing a range
        %   of y values from splitYLim(1) to splitYLim(2) to remove from the axes.
        %   They must be within the current yLimis of the axes.
        %
        %   breakYAxes(splitYLim,splitHeight) splitHeight is the distance to 
        %   seperate the low and high side.  Units are the same as 
        %   get(AX,'uints') default is 0.015
        % 
        %   breakYAxes(splitYLim,splitHeight,xOverhang) xOverhang stretches the 
        %   axis split graphic to extend past the top and bottom of the plot by
        %   the distance set by XOverhang.  Units are the same as get(AX,'units')
        %   default value is 0.015
        %
        %   breakYAxes(AX, ...) performs the operation on the axis specified by AX
        %
        function breakInfo = breakyaxis(varargin)

            %Validate Arguements
            if nargin < 1 || nargin > 4
               error('Wrong number of arguements'); 
            end

            if isscalar(varargin{1}) && ishandle(varargin{1})
                mainAxes = varargin{1};
                argOffset = 1;
                argCnt = nargin - 1;
                if ~strcmp(get(mainAxes,'Type'),'axes')
                   error('Handle object must be Type Axes'); 
                end
            else
                mainAxes = gca;
                argOffset = 0;
                argCnt = nargin;
            end

            if (strcmp(get(mainAxes,'XScale'),'log'))
                error('Log X Axes are not supported'); 
            end

            if (argCnt < 3)
                xOverhang = 0.015;
            else
                xOverhang = varargin{3 + argOffset};
                if  numel(xOverhang) ~= 1 || ~isreal(xOverhang) || ~isnumeric(xOverhang)
                    error('XOverhang must be a scalar number');
                elseif (xOverhang < 0)
                    error('XOverhang must not be negative');
                end
                xOverhang = double(xOverhang);
            end

            if (argCnt < 2)
                splitHeight = 0.015;
            else
                splitHeight = varargin{2 + argOffset};
                if  numel(xOverhang) ~= 1 || ~isreal(xOverhang) || ~isnumeric(xOverhang)
                    error('splitHeight must be a scalar number');
                elseif (xOverhang < 0)
                    error('splitHeight must not be negative');
                end
                splitHeight = double(splitHeight);
            end

            splitYLim = varargin{1 + argOffset};
            if numel(splitYLim) ~= 2 || ~isnumeric(splitYLim) || ~isreal(xOverhang)
               error(splitYLim,'Must be a vector length 2');
            end
            splitYLim = double(splitYLim);

            mainYLim = get(mainAxes,'YLim');
            if (any(splitYLim >= mainYLim(2)) || any(splitYLim <= mainYLim(1)))
               error('splitYLim must be in the range given by get(AX,''YLim'')');
            end

            mainPosition = get(mainAxes,'Position');
            if (splitHeight > mainPosition(3) ) 
               error('Split width is too large') 
            end

            %We need to create 4 axes
            % lowAxes - is used for the low y axis and low pane data
            % highAxes - is used to the high y axis and high pane data
            % annotationAxes - is used to display the x axis and title
            % breakAxes - this is an axes with the same size and position as main
            %   is it used to draw a seperator between the low and high side


            %Grab Some Parameters from the main axis (e.g the one we are spliting)
            mainYLim = get(mainAxes,'YLim');
            mainXLim = get(mainAxes,'XLim');
            mainPosition = get(mainAxes,'Position');
            mainParent = get(mainAxes,'Parent');
            mainHeight = mainPosition(4); %Positions have the format [low bottom width height]
            %mainYRange = mainYLim(2) - mainYLim(1);
            mainFigure = get(mainAxes,'Parent');
            mainXColor = get(mainAxes,'XColor');
            mainLineWidth = get(mainAxes,'LineWidth');
            figureColor = get(mainFigure,'Color');
            mainXTickLabelMode = get(mainAxes,'XTickLabelMode');
            mainYLabel = get(mainAxes,'YLabel');
            mainYDir = get(mainAxes,'YDir');
            mainLayer = get(mainAxes,'Layer');

            %Save Main Axis Z Order
            figureChildren = get(mainFigure,'Children');
            zOrder = find(figureChildren == mainAxes);

            %Calculate where axesLow and axesHigh will be layed on screen
            %And their respctive YLimits
            lowYLimTemp = [mainYLim(1) splitYLim(1)];
            highYLimTemp = [splitYLim(2) mainYLim(2)];

            lowYRangeTemp = lowYLimTemp(2) - lowYLimTemp(1);
            highYRangeTemp = highYLimTemp(2) - highYLimTemp(1);

            lowHeightTemp = lowYRangeTemp / (lowYRangeTemp + highYRangeTemp) * (mainHeight - splitHeight);
            highHeightTemp = highYRangeTemp / (lowYRangeTemp + highYRangeTemp) * (mainHeight - splitHeight);

            lowStretch = (lowHeightTemp + splitHeight/2) / lowHeightTemp;
            lowYRange = lowYRangeTemp * lowStretch;
            lowHeight = lowHeightTemp * lowStretch;

            highStretch = (highHeightTemp + splitHeight/2) / highHeightTemp;
            highYRange = highYRangeTemp * highStretch;
            highHeight = highHeightTemp * highStretch;

            lowYLim = [mainYLim(1) mainYLim(1)+lowYRange];
            highYLim = [mainYLim(2)-highYRange mainYLim(2)];

            if (strcmp(mainYDir, 'normal')) 
                lowPosition = mainPosition;
                lowPosition(4) = lowHeight; 

                highPosition = mainPosition;    %(!!!) look here for position indices!
                highPosition(2) = mainPosition(2) + lowHeight;
                highPosition(4) = highHeight;
            else
                %Low Axis will actually go on the high side a vise versa
                highPosition = mainPosition;
                highPosition(4) = highHeight; 

                lowPosition = mainPosition;
                lowPosition(2) = mainPosition(2) + highHeight;
                lowPosition(4) = lowHeight;
            end

            %Create the Annotations layer, if the Layer is top, draw the axes on
            %top (e.g. after) drawing the low and high pane
            if strcmp(mainLayer,'bottom')
                annotationAxes = CreateAnnotaionAxes(mainAxes,mainParent)
            end

            %Create and position the lowAxes. Remove all X Axis Annotations, the 
            %title, and a potentially offensive tick mark 
            lowAxes = copyobj(mainAxes,mainParent);
            set(lowAxes,'Position', lowPosition, ...
                'YLim', lowYLim, ... 
                'XLim', mainXLim, ...
                'XGrid' ,'off', ...
                'XMinorGrid', 'off', ...
                'XMinorTick','off', ...
                'XTick', [], ...
                'XTickLabel', [], ...
                'box','off');
            if strcmp(mainLayer,'bottom')
                set(lowAxes,'Color','none');
            end
            delete(get(lowAxes,'XLabel')); 
            delete(get(lowAxes,'YLabel'));
            delete(get(lowAxes,'Title'));

            if strcmp(mainXTickLabelMode,'auto')
                yTick =  get(lowAxes,'YTick');
                set(lowAxes,'YTick',yTick(1:(end-1)));
            end

            %Create and position the highAxes. Remove all X Axis annotations, the 
            %title, and a potentially offensive tick mark 
            highAxes = copyobj(mainAxes,mainParent);
            set(highAxes,'Position', highPosition, ...
                'YLim', highYLim, ...
                'XLim', mainXLim, ...
                'XGrid' ,'off', ...
                'XMinorGrid', 'off', ...
                'XMinorTick','off', ...
                'XTick', [], ...
                'XTickLabel', [], ...
                'box','off');
            if strcmp(mainLayer,'bottom') %(!!!) is it only about layers?
                set(highAxes,'Color','none');
            end
            delete(get(highAxes,'XLabel')); 
            delete(get(highAxes,'YLabel'));
            delete(get(highAxes,'Title'));

            if strcmp(mainXTickLabelMode,'auto')
                yTick =  get(highAxes,'YTick');
                set(highAxes,'YTick',yTick(2:end));
            end

                %Create the Annotations layer, if the Layer is top, draw the axes on
            %top (e.g. after) drawing the low and high pane
            if strcmp(mainLayer,'top')
                annotationAxes = CreateAnnotaionAxes(mainAxes,mainParent);
                set(annotationAxes, 'Color','none');
            end

            %Create breakAxes, remove all graphics objects and hide all annotations
            breakAxes = copyobj(mainAxes,mainParent);
            children = get(breakAxes,'Children');
            for i = 1:numel(children)
               delete(children(i)); 
            end

            set(breakAxes,'Color','none');
            %Stretch the breakAxes horizontally to cover the vertical axes lines
            orignalUnits = get(breakAxes,'Units');
            set(breakAxes,'Units','Pixel');
            breakPosition = get(breakAxes,'Position');
            nudgeFactor = get(breakAxes,'LineWidth');
            breakPosition(3) = breakPosition(3) +  nudgeFactor;
            set(breakAxes,'Position',breakPosition);
            set(breakAxes,'Units',orignalUnits);

            %Stretch the breakAxes horizontally to create an overhang for sylistic
            %effect
            breakPosition = get(breakAxes,'Position');
            breakPosition(1) = breakPosition(1) - xOverhang;
            breakPosition(3) = breakPosition(3) +  2*xOverhang;
            set(breakAxes,'Position',breakPosition);

            %Create a sine shaped patch to seperate the 2 sides
            breakYLim = [mainPosition(2) mainPosition(2)+mainPosition(4)];
            set(breakAxes,'ylim',breakYLim);
            theta = linspace(0,2*pi,100);
            xPoints = linspace(mainXLim(1),mainXLim(2),100);
            amp = splitHeight/2 * 0.9;
            yPoints1 = amp * sin(theta) + mainPosition(2) + lowHeightTemp;
            yPoints2 = amp * sin(theta) + mainPosition(2) + mainPosition(4) - highHeightTemp;
            patchPointsY = [yPoints1 yPoints2(end:-1:1) yPoints1(1)];
            patchPointsX = [xPoints  xPoints(end:-1:1)  xPoints(1)];
            patch(patchPointsX,patchPointsY ,figureColor,'EdgeColor',figureColor,'Parent',breakAxes); %use of pathc(!!!)?

            %Create A Line To Delineate the low and high edge of the patch
            line('yData',yPoints1,'xdata',xPoints,'Parent',breakAxes,'Color',mainXColor,'LineWidth',mainLineWidth);
            line('yData',yPoints2,'xdata',xPoints,'Parent',breakAxes,'Color',mainXColor,'LineWidth',mainLineWidth);

            set(breakAxes,'Visible','off');

            %Make the old main axes invisiable
            invisibleObjects = RecursiveSetVisibleOff(mainAxes);

            %Preserve the z-order of the figure
            uistack([lowAxes highAxes breakAxes annotationAxes],'down',zOrder-1)

            %Set the rezise mode to position so that we can dynamically change the
            %size of the figure without screwing things up
            set([lowAxes highAxes breakAxes annotationAxes],'ActivePositionProperty','Position');

            %Playing with the titles labels etc can cause matlab to reposition
            %the axes in some cases.  Mannually force the position to be correct. 
            set([breakAxes annotationAxes],'Position',mainPosition);

            %Save the axes so we can unbreak the axis easily
            breakInfo = struct();
            breakInfo.lowAxes = lowAxes;
            breakInfo.highAxes = highAxes;
            breakInfo.breakAxes = breakAxes;
            breakInfo.annotationAxes = annotationAxes;
            breakInfo.invisibleObjects = invisibleObjects;
        end

        function list = RecursiveSetVisibleOff(handle) 
            list = [];
            list = SetVisibleOff(handle,list);

        end 

        function list = SetVisibleOff(handle, list)
            if (strcmp(get(handle,'Visible'),'on'))
                set(handle,'Visible','off');
                list = [list handle];
            end

            children = get(handle,'Children');
            for i = 1:numel(children)
                list = SetVisibleOff(children(i),list);
            end
        end

        function annotationAxes = CreateAnnotaionAxes(mainAxes,mainParent)

            %Create Annotation Axis, Remove graphics objects, YAxis annotations
            %(except YLabel) and make background transparent
            annotationAxes = copyobj(mainAxes,mainParent);

            set(annotationAxes,'XLimMode','Manual');

            children = get(annotationAxes,'Children');
            for i = 1:numel(children)
               delete(children(i)); 
            end

            %Save the yLabelpostion because it will move when we delete yAxis
            %ticks
            yLabel = get(annotationAxes,'YLabel');
            yLabelPosition = get(yLabel,'Position');

            set(annotationAxes,'YGrid' ,'off', ...
                'YMinorGrid', 'off', ...
                'YMinorTick','off', ...
                'YTick', [], ...
                'YTickLabel', []);

            %Restore the pevious label postition
            set(yLabel,'Position',yLabelPosition);
        end


        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Rutinas Para exportar graficos en Latex-Tikz
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function Matlab2tikz_______________________l() 
        end  % end function   


        function M2T_matlab2tikz(varargin)
            %MATLAB2TIKZ    Save figure in native LaTeX (TikZ/Pgfplots).
            %   MATLAB2TIKZ() saves the current figure as LaTeX file.
            %   MATLAB2TIKZ comes with several options that can be combined at will.
            %
            %   MATLAB2TIKZ(FILENAME,...) or MATLAB2TIKZ('filename',FILENAME,...)
            %   stores the LaTeX code in FILENAME.
            %
            %   MATLAB2TIKZ('filehandle',FILEHANDLE,...) stores the LaTeX code in the file
            %   referenced by FILEHANDLE. (default: [])
            %
            %   MATLAB2TIKZ('figurehandle',FIGUREHANDLE,...) explicitly specifies the
            %   handle of the figure that is to be stored. (default: gcf)
            %
            %   MATLAB2TIKZ('colormap',DOUBLE,...) explicitly specifies the colormap to be
            %   used. (default: current color map)
            %
            %   MATLAB2TIKZ('strict',BOOL,...) tells MATLAB2TIKZ to adhere to MATLAB(R)
            %   conventions wherever there is room for relaxation. (default: false)
            %
            %   MATLAB2TIKZ('strictFontSize',BOOL,...) retains the exact font sizes
            %   specified in MATLAB for the TikZ code. This goes against normal LaTeX
            %   practice. (default: false)
            %
            %   MATLAB2TIKZ('showInfo',BOOL,...) turns informational output on or off.
            %   (default: true)
            %
            %   MATLAB2TIKZ('showWarnings',BOOL,...) turns warnings on or off.
            %   (default: true)
            %
            %   MATLAB2TIKZ('imagesAsPng',BOOL,...) stores MATLAB(R) images as (lossless)
            %   PNG files. This is more efficient than storing the image color data as TikZ
            %   matrix. (default: true)
            %
            %   MATLAB2TIKZ('externalData',BOOL,...) stores all data points in external
            %   files as tab separated values (TSV files). (default: false)
            %
            %   MATLAB2TIKZ('dataPath',CHAR, ...) defines where external data files
            %   and/or PNG figures are saved. It can be either an absolute or a relative
            %   path with respect to your MATLAB work directory. By default, data files are
            %   placed in the same directory as the TikZ output file. To place data files
            %   in your MATLAB work directory, you can use '.'. (default: [])
            %
            %   MATLAB2TIKZ('relativeDataPath',CHAR, ...) tells MATLAB2TIKZ to use the
            %   given path to follow the external data files and PNG files. This is the
            %   relative path from your main LaTeX file to the data file directory.
            %   By default the same directory is used as the output (default: [])
            %
            %   MATLAB2TIKZ('height',CHAR,...) sets the height of the image. This can be
            %   any LaTeX-compatible length, e.g., '3in' or '5cm' or '0.5\textwidth'.  If
            %   unspecified, MATLAB2TIKZ tries to make a reasonable guess.
            %
            %   MATLAB2TIKZ('width',CHAR,...) sets the width of the image.
            %   If unspecified, MATLAB2TIKZ tries to make a reasonable guess.
            %
            %   MATLAB2TIKZ('noSize',BOOL,...) determines whether 'width', 'height', and
            %   'scale only axis' are specified in the generated TikZ output. For compatibility with the
            %   tikzscale package set this to true. (default: false)
            %
            %   MATLAB2TIKZ('extraCode',CHAR or CELLCHAR,...) explicitly adds extra code
            %   at the beginning of the output file. (default: [])
            %
            %   MATLAB2TIKZ('extraCodeAtEnd',CHAR or CELLCHAR,...) explicitly adds extra
            %   code at the end of the output file. (default: [])
            %
            %   MATLAB2TIKZ('extraAxisOptions',CHAR or CELLCHAR,...) explicitly adds extra
            %   options to the Pgfplots axis environment. (default: [])
            %
            %   MATLAB2TIKZ('extraColors', {{'name',[R G B]}, ...} , ...) adds
            %   user-defined named RGB-color definitions to the TikZ output.
            %   R, G and B are expected between 0 and 1. (default: {})
            %
            %   MATLAB2TIKZ('extraTikzpictureOptions',CHAR or CELLCHAR,...)
            %   explicitly adds extra options to the tikzpicture environment. (default: [])
            %
            %   MATLAB2TIKZ('encoding',CHAR,...) sets the encoding of the output file.
            %
            %   MATLAB2TIKZ('floatFormat',CHAR,...) sets the format used for float values.
            %   You can use this to decrease the file size. (default: '%.15g')
            %
            %   MATLAB2TIKZ('maxChunkLength',INT,...) sets maximum number of data points
            %   per \addplot for line plots (default: 4000)
            %
            %   MATLAB2TIKZ('parseStrings',BOOL,...) determines whether title, axes labels
            %   and the like are parsed into LaTeX by MATLAB2TIKZ's parser.
            %   If you want greater flexibility, set this to false and use straight LaTeX
            %   for your labels. (default: true)
            %
            %   MATLAB2TIKZ('LibAv05.M2T_parseStringsAsMath',BOOL,...) determines whether to use TeX's
            %   math mode for more characters (e.g. operators and figures). (default: false)
            %
            %   MATLAB2TIKZ('showHiddenStrings',BOOL,...) determines whether to show
            %   strings whose were deliberately hidden. This is usually unnecessary, but
            %   can come in handy for unusual plot types (e.g., polar plots). (default:
            %   false)
            %
            %   MATLAB2TIKZ('interpretTickLabelsAsTex',BOOL,...) determines whether to
            %   interpret tick labels as TeX. MATLAB(R) doesn't allow to do that in R2014a
            %   or before. In R2014b and later, please set the "TickLabelInterpreter"
            %   property of the relevant axis to get the same effect. (default: false)
            %
            %   MATLAB2TIKZ('arrowHeadSize', FLOAT, ...) allows to resize the arrow heads
            %   in quiver plots by rescaling the arrow heads by a positive scalar. (default: 10)
            %
            %   MATLAB2TIKZ('tikzFileComment',CHAR,...) adds a custom comment to the header
            %   of the output file. (default: '')
            %
            %   MATLAB2TIKZ('LibAv05.M2T_addLabels',BOOL,...) add labels to plots: using Tag property
            %   or automatic names (where applicable) which make it possible to refer to
            %   them using \ref{...} (e.g., in the caption of a figure). (default: false)
            %
            %   MATLAB2TIKZ('standalone',BOOL,...) determines whether to produce
            %   a standalone compilable LaTeX file. Setting this to true may be useful for
            %   taking a peek at what the figure will look like. (default: false)
            %
            %   MATLAB2TIKZ('checkForUpdates',BOOL,...) determines whether to automatically
            %   check for updates of LibAv05.M2T_matlab2tikz. (default: true (if not using git))
            %
            %   MATLAB2TIKZ('semanticLineWidths',CELLMATRIX,...) allows you to customize
            %   the mapping of semantic "line width" values.
            %   A valid entry is an Nx2 cell array:
            %     - the first column contains the semantic names,
            %     - the second column contains the corresponding line widths in points.
            %   The entries you provide are used in addition to the pgf defaults:
            %     {'ultra thin', 0.1; 'very thin' , 0.2; 'thin', 0.4; 'semithick', 0.6;
            %      'thick'     , 0.8; 'very thick', 1.2; 'ultra thick', 1.6}
            %   or a single "NaN" can be provided to turn off this feature alltogether.
            %   If you specify the default names, their mapping will be overwritten.
            %   Inside your LaTeX document, you are responsible to make sure these TikZ
            %   styles are properly defined.
            %   (Default: NaN)
            %
            %   Example
            %      x = -pi:pi/10:pi;
            %      y = tan(sin(x)) - sin(tan(x));
            %      plot(x,y,'--rs');
            %      LibAv05.M2T_matlab2tikz('myfile.tex');
            %
            %   See also: cleanfigure

            %% Check if we are in MATLAB or Octave.
            minimalVersion = struct('MATLAB', struct('name','2014a', 'num',[8 3]), ...
                                    'Octave', struct('name','3.8', 'num',[3 8]));
            LibAv05.M2T_checkDeprecatedEnvironment(minimalVersion);

            m2t.args            = []; % For command line arguments
            m2t.current         = []; % For currently active objects
            m2t.transform       = []; % For hgtransform groups
            m2t.pgfplotsVersion = [1,3];
            m2t.about.name      = 'LibAv05.M2T_matlab2tikz';
            m2t.about.version   = '1.1.0';
            m2t.about.years     = '2008--2016';
            m2t.about.website   = 'http://www.mathworks.com/matlabcentral/fileexchange/22022-LibAv05.M2T_matlab2tikz-LibAv05.M2T_matlab2tikz';
            m2t.about.github    = 'https://github.com/LibAv05.M2T_matlab2tikz/LibAv05.M2T_matlab2tikz';
            m2t.about.wiki      = [m2t.about.github '/wiki'];
            m2t.about.issues    = [m2t.about.github '/issues'];
            m2t.about.develop   = [m2t.about.github '/tree/develop'];
            VCID = M2T_VersionControlIdentifier();
            m2t.about.versionFull = strtrim(sprintf('v%s %s', m2t.about.version, VCID));

            m2t.tol = 1.0e-15; % numerical tolerance (e.g. used to LibAv05.M2T_test equality of doubles)

            % the actual contents of the TikZ file go here
            m2t.content = struct('name',     '', ...
                                 'comment',  [], ...
                                 'options',  {LibAv05.M2T_opts_new()}, ...
                                 'content',  {cell(0)}, ...
                                 'children', {cell(0)});
            m2t.preamble = sprintf(['\\usepackage[T1]{fontenc}\n', ...
                                    '\\usepackage[utf8]{inputenc}\n', ...
                                    '\\usepackage{pgfplots}\n', ...
                                    '\\usepackage{grffile}\n', ...
                                    '\\pgfplotsset{compat=newest}\n', ...
                                    '\\usetikzlibrary{plotmarks}\n', ...
                                    '\\usetikzlibrary{arrows.meta}\n', ...
                                    '\\usepgfplotslibrary{patchplots}\n', ...
                                    '\\usepackage{amsmath}\n']);

            %% scan the options
            ipp = m2tInputParser;

            ipp = ipp.addOptional(ipp, 'filename',   '', @(x) LibAv05.M2T_filenameValidation(x,ipp));
            ipp = ipp.addOptional(ipp, 'filehandle', [], @LibAv05.M2T_filehandleValidation);

            ipp = ipp.addParamValue(ipp, 'figurehandle', get(0,'CurrentFigure'), @ishandle);
            ipp = ipp.addParamValue(ipp, 'colormap', [], @isnumeric);
            ipp = ipp.addParamValue(ipp, 'strict', false, @islogical);
            ipp = ipp.addParamValue(ipp, 'strictFontSize', false, @islogical);
            ipp = ipp.addParamValue(ipp, 'showInfo', true, @islogical);
            ipp = ipp.addParamValue(ipp, 'showWarnings', true, @islogical);
            ipp = ipp.addParamValue(ipp, 'checkForUpdates', isempty(VCID), @islogical);

            ipp = ipp.addParamValue(ipp, 'semanticLineWidths', NaN, @LibAv05.M2T_isValidSemanticLineWidthDefinition);

            ipp = ipp.addParamValue(ipp, 'encoding' , '', @ischar);
            ipp = ipp.addParamValue(ipp, 'standalone', false, @islogical);
            ipp = ipp.addParamValue(ipp, 'tikzFileComment', '', @ischar);
            ipp = ipp.addParamValue(ipp, 'extraColors', {}, @LibAv05.M2T_isColorDefinitions);
            ipp = ipp.addParamValue(ipp, 'extraCode', {}, @LibAv05.M2T_isCellOrChar);
            ipp = ipp.addParamValue(ipp, 'extraCodeAtEnd', {}, @LibAv05.M2T_isCellOrChar);
            ipp = ipp.addParamValue(ipp, 'extraAxisOptions', {}, @LibAv05.M2T_isCellOrChar);
            ipp = ipp.addParamValue(ipp, 'extraTikzpictureOptions', {}, @LibAv05.M2T_isCellOrChar);
            ipp = ipp.addParamValue(ipp, 'floatFormat', '%.15g', @ischar);
            ipp = ipp.addParamValue(ipp, 'automaticLabels', false, @islogical);
            ipp = ipp.addParamValue(ipp, 'LibAv05.M2T_addLabels', false, @islogical);
            ipp = ipp.addParamValue(ipp, 'showHiddenStrings', false, @islogical);
            ipp = ipp.addParamValue(ipp, 'height', '', @ischar);
            ipp = ipp.addParamValue(ipp, 'width' , '', @ischar);
            ipp = ipp.addParamValue(ipp, 'imagesAsPng', true, @islogical);
            ipp = ipp.addParamValue(ipp, 'externalData', false, @islogical);
            ipp = ipp.addParamValue(ipp, 'dataPath', '', @ischar);
            ipp = ipp.addParamValue(ipp, 'relativeDataPath', '', @ischar);
            ipp = ipp.addParamValue(ipp, 'noSize', false, @islogical);
            ipp = ipp.addParamValue(ipp, 'arrowHeadSize', 10, @(x) x>0);

            % Maximum chunk length.
            % TeX parses files line by line with a buffer of size buf_size. If the
            % plot has too many data points, pdfTeX's buffer size may be exceeded.
            % As a work-around, the plot is split into several smaller chunks.
            %
            % What is a "large" array?
            % TeX parser buffer is buf_size=200 000 char on Mac TeXLive, let's say
            % 100 000 to be on the safe side.
            % 1 point is represented by 25 characters (estimation): 2 coordinates (10
            % char), 2 brackets, comma and white space, + 1 extra char.
            % That gives a magic arbitrary number of 4000 data points per array.
            ipp = ipp.addParamValue(ipp, 'maxChunkLength', 4000, @isnumeric);

            % By default strings like axis labels are parsed to match the appearance of
            % strings as closely as possible to that generated by MATLAB.
            % If the user wants to have particular strings in the LibAv05.M2T_matlab2tikz output that
            % can't be generated in MATLAB, they can disable string parsing. In that case
            % all strings are piped literally to the LaTeX output.
            ipp = ipp.addParamValue(ipp, 'parseStrings', true, @islogical);

            % In addition to regular string parsing, an additional stage can be enabled
            % which uses TeX's math mode for more characters like figures and operators.
            ipp = ipp.addParamValue(ipp, 'LibAv05.M2T_parseStringsAsMath', false, @islogical);

            % As opposed to titles, axis labels and such, MATLAB(R) does not interpret tick
            % labels as TeX. LibAv05.M2T_matlab2tikz retains this behavior, but if it is desired to
            % interpret the tick labels as TeX, set this option to true.
            ipp = ipp.addParamValue(ipp, 'interpretTickLabelsAsTex', false, @islogical);

            %% deprecated parameters (will auto-generate warnings upon parse)
            ipp = ipp.addParamValue(ipp, 'relativePngPath', '', @ischar);
            ipp = ipp.deprecateParam(ipp, 'relativePngPath', 'relativeDataPath');
            ipp = ipp.deprecateParam(ipp, 'automaticLabels', 'LibAv05.M2T_addLabels');

            %% Finally parse all the arguments
            ipp = ipp.parse(ipp, varargin{:});
            m2t.args = ipp.Results; % store the input arguments back into the m2t data struct

            %% Inform users of potentially dangerous options
            LibAv05.M2T_warnAboutParameter(m2t, 'LibAv05.M2T_parseStringsAsMath', @(opt)(opt==true), ...
                ['This may produce undesirable string output. For full control over output\n', ...
                 'strings please set the parameter "parseStrings" to false.']);
            LibAv05.M2T_warnAboutParameter(m2t, 'noSize', @(opt)(opt==true), ...
                 'This may impede both axes sizing and placement!');
            LibAv05.M2T_warnAboutParameter(m2t, 'imagesAsPng', @(opt)(opt==false), ...
                 ['It is highly recommended to use PNG data to store images.\n', ...
                  'Make sure to set "imagesAsPng" to true.']);

            %% Do some global initialization
            m2t.color = M2T_configureColors(m2t.args.extraColors);
            m2t.semantic.LineWidth = M2T_configureSemanticLineWidths(m2t.args.semanticLineWidths);

            % define global counter variables
            m2t.count.pngFile     = 0; % number of PNG files
            m2t.count.tsvFile     = 0; % number of TSV files
            m2t.count.autolabel   = 0; % number of automatic labels
            m2t.count.plotyylabel = 0; % number of plotyy labels

            %% shortcut
            m2t.ff = m2t.args.floatFormat;

            %% add global elements
            if isempty(m2t.args.figurehandle)
                error('LibAv05.M2T_matlab2tikz:figureNotFound','MATLAB figure not found.');
            end
            m2t.current.gcf = m2t.args.figurehandle;
            if m2t.args.colormap
                m2t.current.colormap = m2t.args.colormap;
            else
                m2t.current.colormap = get(m2t.current.gcf, 'colormap');
            end

            %% handle output file handle/file name
            [m2t, fid, fileWasOpen] = M2T_openFileForOutput(m2t);

            % By default, reference the PNG (if required) from the TikZ file
            % as the file path of the TikZ file itself. This works if the MATLAB script
            % is executed in the same folder where the TeX file sits.
            if isempty(m2t.args.relativeDataPath)
                if ~isempty(m2t.args.relativePngPath)
                    %NOTE: eventually break backwards compatibility of relative PNG path
                    m2t.relativeDataPath = m2t.args.relativePngPath;
                    LibAv05.M2T_userWarning(m2t, ['Using "relativePngPath" for "relativeDataPath".', ...
                        ' This will stop working in a future release.']);
                else
                    m2t.relativeDataPath = m2t.args.relativeDataPath;
                end
            else
                m2t.relativeDataPath = m2t.args.relativeDataPath;
            end
            if isempty(m2t.args.dataPath)
                m2t.dataPath = fileparts(m2t.tikzFileName);
            else
                m2t.dataPath = m2t.args.dataPath;
            end

            %% print some version info to the screen
            LibAv05.M2T_userInfo(m2t, ['(To disable info messages, pass [''showInfo'', false] to LibAv05.M2T_matlab2tikz.)\n', ...
                '(For all other options, type ''help LibAv05.M2T_matlab2tikz''.)\n']);

            LibAv05.M2T_userInfo(m2t, '\nThis is %s %s.\n', m2t.about.name, m2t.about.versionFull)

            % In Octave, put a new line and some spaces in between the URLs for clarity.
            % In MATLAB this is not necessary, since the URLs get (shorter) descriptions.
            sep = M2T_switchMatOct('', sprintf('\n '));
            versionInfo = ['The laLibAv05.M2T_test developments can be retrieved from %s.\n', ...
                           'You can find more documentation on %s and %s.\n', ...
                           'If you encounter bugs or want a new feature, go to %s.\n', ...
                           'Please visit %s to rate %s or download the stable release.\n'];
            LibAv05.M2T_userInfo(m2t, versionInfo, ...
                     LibAv05.M2T_clickableUrl(m2t.about.develop, 'our development branch'), ...
                     [sep LibAv05.M2T_clickableUrl(m2t.about.github, 'our GitHub page') sep], ...
                     [sep LibAv05.M2T_clickableUrl(m2t.about.wiki, 'wiki')], ...
                     [sep LibAv05.M2T_clickableUrl(m2t.about.issues, 'our issue tracker')],...
                     [LibAv05.M2T_clickableUrl(m2t.about.website, 'FileExchange') sep],...
                     m2t.about.name);

            %% Save the figure as TikZ to file
            m2t = M2T_saveToFile(m2t, fid, fileWasOpen);

            %% Check for a new LibAv05.M2T_matlab2tikz version outside version control
            if m2t.args.checkForUpdates
                m2tUpdater(m2t.about, m2t.args.showInfo);
            end

        end
        % ==============================================================================
        function [m2t, counterValue] = M2T_incrementGlobalCounter(m2t, counterName)
            % Increments a global counter value and returns its value
            m2t.count.(counterName) = m2t.count.(counterName) + 1;
            counterValue = m2t.count.(counterName);
        end
        % ==============================================================================
        function colorConfig = M2T_configureColors(extraColors)
            % Sets the global color options for LibAv05.M2T_matlab2tikz
            colorConfig = struct();

            % Set the color resolution.
            colorConfig.depth     = 48; %[bit] RGB color depth (typical values: 24, 30, 48)
            colorConfig.precision = 2^(-colorConfig.depth/3);
            colorConfig.format    = sprintf('%%0.%df',ceil(-log10(colorConfig.precision)));

            % The following color RGB-values which will need to be defined:
            %
            %   - 'extraNames' contains their designated names,
            %   - 'extraSpecs' their RGB specifications.
            [colorConfig.extraNames, colorConfig.extraSpecs] = ...
                LibAv05.M2T_dealColorDefinitions(extraColors);
        end
        % ==============================================================================
        function [m2t, fid, fileWasOpen] = M2T_openFileForOutput(m2t)
            % opens the output file and/or show a dialog to select one
            if ~isempty(m2t.args.filehandle)
                fid         = m2t.args.filehandle;
                fileWasOpen = true;
                if ~isempty(m2t.args.filename)
                    LibAv05.M2T_userWarning(m2t, ...
                        'File handle AND file name for output given. File handle used, file name discarded.')
                end
                m2t.tikzFileName = fopen(fid);
            else
                fid         = [];
                fileWasOpen = false;
                % set filename
                if ~isempty(m2t.args.filename)
                    filename = m2t.args.filename;
                else
                    [filename, pathname] = uiputfile({'*.tex;*.tikz'; '*.*'}, 'Save File');
                    filename = fullfile(pathname, filename);
                end
                m2t.tikzFileName = filename;
            end

        end
        % ==============================================================================
        function l = M2T_filenameValidation(x, p)
            % is the filename argument NOT another keyword?
            l = ischar(x) && ~any(strcmp(x,p.Parameters)); %FIXME: See #471
        end
        % ==============================================================================
        function l = M2T_filehandleValidation(x)
            % is the filehandle the handle to an opened file?
            l = isnumeric(x) && any(x==fopen('all'));
        end
        % ==============================================================================
        function bool = M2T_isCellOrChar(x)
            bool = iscell(x) || ischar(x);
        end
        % ==============================================================================
        function bool = M2T_isRGBTuple(color)
            % Returns true when the color is a valid RGB tuple
            bool = numel(color) == 3 && ...
                   all(isreal(color)) && ...
                   all( 0<=color & color<=1 ); % this also disallows NaN entries
        end
        % ==============================================================================
        function bool = M2T_isColorDefinitions(colors)
            % Returns true when the input is a cell array of color definitions, i.e.
            %  a cell array with in each cell a cell of the form {'name', [R G B]}
            isValidEntry = @(e)( iscell(e) && ischar(e{1}) && LibAv05.M2T_isRGBTuple(e{2}) );

            bool = iscell(colors) && all(cellfun(isValidEntry, colors));
        end
        % ==============================================================================
        function bool = M2T_isValidSemanticLineWidthDefinition(defMat)
            % Returns true when the input is a cell array of shape Nx2 and
            % contents in each column a set of string and numerical value as needed
            % for the semanticLineWidth option.
            bool = iscell(defMat) && size(defMat, 2) == 2; % Nx2 cell array
            bool = bool && all(cellfun(@ischar   , defMat(:,1))); % first column: names
            bool = bool && all(cellfun(@isnumeric, defMat(:,2))); % second column: line width in points

            % alternatively: just 1 NaN to remove the defaults
            bool = bool || (numel(defMat)==1 && isnan(defMat));
        end
        % ==============================================================================
        function fid = M2T_fileOpenForWrite(m2t, filename)
            % Set the encoding of the output file.
            % Currently only MATLAB supports different encodings.
            fid = -1;

            [filepath] = fileparts(filename);
            if ~exist(filepath,'dir') && ~isempty(filepath)
                mkdir(filepath);
            end

            switch getEnvironment()
                case 'MATLAB'
                    fid = fopen(filename, 'w', ...
                                'native', m2t.args.encoding);
                case 'Octave'
                    fid = fopen(filename, 'w');
                otherwise
                    errorUnknownEnvironment();
            end

            if fid == -1
                error('LibAv05.M2T_matlab2tikz:fileOpenError', ...
                    'Unable to open file ''%s'' for writing.', filename);
            end
        end
        % ==============================================================================
        function path = M2T_TeXpath(path)
            path = strrep(path, filesep, '/');
            % TeX uses '/' as a file separator (as UNIX). Windows, however, uses
            % '\' which is not supported by TeX as a file separator
        end
        % ==============================================================================
        function m2t = M2T_saveToFile(m2t, fid, fileWasOpen)
            % Save the figure as TikZ to a file. All other routines are called from here.

            % get all axes handles
            [m2t, axesHandles] = M2T_findPlotAxes(m2t, m2t.current.gcf);

            % Turn around the handles vector to make sure that plots that appeared
            % first also appear first in the vector. This makes sure the z-order of
            % superimposed axes is respected and is fundamental for plotyy.
            axesHandles = axesHandles(end:-1:1);

            % Alternative Positioning of axes.
            % Select relevant Axes and draw them.
            [m2t, axesBoundingBox] = M2T_getRelevantAxes(m2t, axesHandles);

            m2t.axesBoundingBox = axesBoundingBox;
            m2t.axes = {};
            for relevantAxesHandle = m2t.relevantAxesHandles(:)'
                m2t = M2T_drawAxes(m2t, relevantAxesHandle);
            end

            % Handle color bars.
            for cbar = m2t.cbarHandles(:)'
                m2t = M2T_handleColorbar(m2t, cbar);
            end

            % Draw annotations
            m2t = M2T_drawAnnotations(m2t);

            % Add all axes containers to the file contents.
            for axesContainer = m2t.axes
                m2t.content = M2T_addChildren(m2t.content, axesContainer);
            end

            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % actually print the stuff
            minimalPgfplotsVersion = M2T_formatPgfplotsVersion(m2t.pgfplotsVersion);

            m2t.content.comment = sprintf('This file was created by %s.\n', m2t.about.name);

            if m2t.args.showInfo
                % disable this info if showInfo=false
                m2t.content.comment = [m2t.content.comment, ...
                    sprintf(['\n',...
                    'The laLibAv05.M2T_test updates can be retrieved from\n', ...
                    '  %s\n', ...
                    'where you can also make suggestions and rate %s.\n'], ...
                    m2t.about.website, m2t.about.name ) ...
                    ];
            end

            LibAv05.M2T_userInfo(m2t, 'You will need pgfplots version %s or newer to compile the TikZ output.',...
                minimalPgfplotsVersion);

            % Add custom comment.
            if ~isempty(m2t.args.tikzFileComment)
                m2t.content.comment = [m2t.content.comment, ...
                    sprintf('\n%s\n', m2t.args.tikzFileComment)
                    ];
            end

            m2t.content.name = 'tikzpicture';

            % Add custom TikZ options if any given.
            m2t.content.options = LibAv05.M2T_LibAv05.M2T_opts_append_userdefined(m2t.content.options, ...
                                           m2t.args.extraTikzpictureOptions);

            m2t.content.colors = M2T_generateColorDefinitions(m2t.color);

            % Open file if was not open
            if ~fileWasOpen
                fid = M2T_fileOpenForWrite(m2t, m2t.tikzFileName);
                finally_fclose_fid = onCleanup(@() fclose(fid));
            end

            % Finally print it to the file
            LibAv05.M2T_addComments(fid, m2t.content.comment);
            LibAv05.M2T_addStandalone(m2t, fid, 'preamble');
            LibAv05.M2T_addCustomCode(fid, '', m2t.args.extraCode, '');
            LibAv05.M2T_addStandalone(m2t, fid, 'begin');

            LibAv05.M2T_printAll(m2t, m2t.content, fid); % actual plotting happens here

            LibAv05.M2T_addCustomCode(fid, '\n', m2t.args.extraCodeAtEnd, '');

            LibAv05.M2T_addStandalone(m2t, fid, 'end');
        end
        % ==============================================================================
        function M2T_addStandalone(m2t, fid, part)
            % writes a part of a standalone LaTeX file definition
            if m2t.args.standalone
                switch part
                    case 'preamble'
                        fprintf(fid, '\\documentclass[tikz]{standalone}\n%s\n',  m2t.preamble);
                    case 'begin'
                        fprintf(fid, '\\begin{document}\n');
                    case 'end'
                        fprintf(fid, '\n\\end{document}');
                    otherwise
                        error('m2t:unknownStandalonePart', ...
                              'Unknown standalone part "%s"', part);
                end
            end
        end
        % ==============================================================================
        function str = M2T_generateColorDefinitions(colorConfig)
            % Output the color definitions to LaTeX
            str   = '';
            names = colorConfig.extraNames;
            specs = colorConfig.extraSpecs;
            ff    = colorConfig.format;

            if ~isempty(names)
                colorDef = cell(1, length(names));
                for k = 1:length(names)
                    % Append with '%' to avoid spacing woes in LaTeX
                    FORMAT      = ['\\definecolor{%s}{rgb}{' ff ',' ff ',' ff '}%%\n'];
                    colorDef{k} = sprintf(FORMAT, names{k}, specs{k});
                end
                str = m2tstrLibAv05.M2T_join([colorDef, sprintf('%%\n')], '');
            end
        end
        % ==============================================================================
        function [m2t, axesHandles] = M2T_findPlotAxes(m2t, fh)
            % find axes handles that are not legends/colorbars
            % store detected legends and colorbars in 'm2t'
            % fh            figure handle
            axesHandles = findall(fh, 'type', 'axes');

            % Remove all legend handles, as they are treated separately.
            if ~isempty(axesHandles)
                % TODO fix for octave
                tagKeyword = M2T_switchMatOct('Tag', 'tag');
                % Find all legend handles. This is MATLAB-only.
                m2t.legendHandles = findall(fh, tagKeyword, 'legend');
                m2t.legendHandles = m2t.legendHandles(:)';
                idx               = ~ismember(axesHandles, m2t.legendHandles);
                axesHandles       = axesHandles(idx);
            end

            % Remove all colorbar handles, as they are treated separately.
            if ~isempty(axesHandles)
                colorbarKeyword = M2T_switchMatOct('Colorbar', 'colorbar');
                % Find all colorbar handles. This is MATLAB-only.
                cbarHandles = findall(fh, tagKeyword, colorbarKeyword);
                % Octave also finds text handles here; no idea why. Filter.
                m2t.cbarHandles = [];
                for h = cbarHandles(:)'
                    if any(strcmpi(get(h, 'Type'),{'axes','colorbar'}))
                        m2t.cbarHandles = [m2t.cbarHandles, h];
                    end
                end
                m2t.cbarHandles = m2t.cbarHandles(:)';
                idx             = ~ismember(axesHandles, m2t.cbarHandles);
                axesHandles     = axesHandles(idx);
            else
                m2t.cbarHandles = [];
            end

            % Remove scribe layer holding annotations (MATLAB < R2014b)
            m2t.scribeLayer = findall(axesHandles, 'Tag','scribeOverlay');
            idx             = ~ismember(axesHandles, m2t.scribeLayer);
            axesHandles     = axesHandles(idx);
        end
        % ==============================================================================
        function M2T_addComments(fid, comment)
            % prints TeX comments to file stream |fid|
            if ~isempty(comment)
                newline = sprintf('\n');
                newlineTeX = sprintf('\n%%');
                fprintf(fid, '%% %s\n', strrep(comment, newline, newlineTeX));
            end
        end
        % ==============================================================================
        function M2T_addCustomCode(fid, before, code, after)
            if ~isempty(code)
                fprintf(fid, before);
                if ischar(code)
                    code = {code};
                end
                if iscellstr(code)
                    for str = code(:)'
                        fprintf(fid, '%s\n', str{1});
                    end
                else
                    error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_saveToFile', 'Need str or cellstr.');
                end
                fprintf(fid,after);
            end
        end
        % ==============================================================================
        function [m2t, pgfEnvironments] = M2T_handleAllChildren(m2t, h)
            % Draw all children of a graphics object (if they need to be drawn).
            % #COMPLEX: mainly a switch-case
            str = '';
            children = allchild(h);

            % prepare cell array of pgfEnvironments
            pgfEnvironments = cell(1, numel(children));
            envCounter      = 1;

            % It's important that we go from back to front here, as this is
            % how MATLAB does it, too. Significant for patch (contour) plots,
            % and the order of plotting the colored patches.
            for child = children(end:-1:1)'

                % Check if object has legend. Some composite objects need to determine
                % their status at the root level. For detailed explanations check
                % LibAv05.M2T_getLegendEntries().
                % TODO: could move this check into LibAv05.M2T_drawHggroup. Need to verify how
                % hgtransform behaves though. (priority - LOW)
                m2t = M2T_hasLegendEntry(m2t,child);

                switch char(get(child, 'Type'))
                    % 'axes' environments are treated separately.

                    case 'line'
                        [m2t, str] = M2T_drawLine(m2t, child);

                    case 'patch'
                        [m2t, str] = M2T_drawPatch(m2t, child);

                    case 'image'
                        [m2t, str] = M2T_drawImage(m2t, child);

                    case {'hggroup', 'matlab.graphics.primitive.Group', ...
                          'scatter', 'bar', 'stair', 'stem' ,'errorbar', 'area', ...
                          'quiver','contour'}
                        [m2t, str] = M2T_drawHggroup(m2t, child);

                    case 'hgtransform'
                        % From http://www.mathworks.de/de/help/matlab/ref/hgtransformproperties.html:
                        % Matrix: 4-by-4 matrix
                        %   Transformation matrix applied to hgtransform object and its
                        %   children. The hgtransform object applies the transformation
                        %   matrix to all its children.
                        % More information at http://www.mathworks.de/de/help/matlab/creating_plots/group-objects.html.
                        m2t.transform = get(child, 'Matrix');
                        [m2t, str] = M2T_handleAllChildren(m2t, child);
                        m2t.transform = [];

                    case 'surface'
                        [m2t, str] = M2T_drawSurface(m2t, child);

                    case 'text'
                        [m2t, str] = M2T_drawVisibleText(m2t, child);

                    case 'rectangle'
                        [m2t, str] = M2T_drawRectangle(m2t, child);

                    case 'histogram'
                        [m2t, str] = M2T_drawHistogram(m2t, child);

                    case guitypes()
                        % don't do anything for GUI objects and their children
                        str = '';

                    case 'light'
                        % These objects are not supported and should not/cannot be
                        % supported by LibAv05.M2T_matlab2tikz or pgfplots.

                    case ''
                        % No children found for handle. (It has only a title and/or
                        % labels). Carrying on as if nothing happened

                    otherwise
                        error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_handleAllChildren',                 ...
                            'I don''t know how to handle this object: %s\n', ...
                            get(child, 'Type'));

                end

                % A composite object might nest LibAv05.M2T_handleAllChildren calls that can
                % modify the m2t.currentHandleHasLegend value. Re-instate the
                % legend status. For detailed explanations check LibAv05.M2T_getLegendEntries().
                m2t                          = M2T_hasLegendEntry(m2t,child);
                [m2t, legendLabel, labelRef] = M2T_addPlotyyReference(m2t, child);
                legendInfo                   = M2T_addLegendInformation(m2t, child);
                % Add labelRef BEFORE next plot to preserve color order
                str = M2T_join(m2t, {labelRef, str, legendLabel, legendInfo}, '');

                % append the environment
                pgfEnvironments{envCounter} = str;
                envCounter = envCounter +1;
            end
        end
        % ==============================================================================
        function [m2t, label, labelRef] = M2T_addPlotyyReference(m2t, h)
            % Create labelled references to legend entries of the main plotyy axis

            % This ensures we are either on the main or secondary axis
            label    = '';
            labelRef = '';
            if ~LibAv05.M2T_isAxisPlotyy(m2t.current.gca)
                return
            end

            % Get current label counter

            if LibAv05.M2T_hasPlotyyReference(m2t,h)
                % Label the plot to later reference it. Only legend entries on the main
                % plotyy axis will have a label
                [m2t, labelNum] = M2T_incrementGlobalCounter(m2t, 'plotyylabel');
                label = sprintf('\\label{%s}\n\n', LibAv05.M2T_plotyylabelName(labelNum));

            elseif m2t.currentHandleHasLegend && ~isempty(m2t.axes{end}.PlotyyReferences)
                % We are on the secondary axis.

                % We have produced a number of labels we can refer to so far.
                % Also, here we have a number of references that are to be recorded.
                % So, we make the last references (assuming the other ones have been
                % realized already)
                nReferences  = numel(m2t.axes{end}.PlotyyReferences);
                nLabels      = m2t.count.plotyylabel;

                % This is the range of labels, corresponding to the references
                labelRange   = (nLabels-nReferences+1):nLabels;

                labelRef = cell(1, numel(labelRange));
                % Create labelled references to legend entries of the main axis
                for iRef = 1:nReferences
                    ref            = m2t.axes{end}.PlotyyReferences(iRef);
                    lString        = M2T_getLegendString(m2t,ref);
                    labelRef{iRef} = sprintf('\\addlegendimage{/pgfplots/refstyle=%s}\n\\addlegendentry{%s}\n',...
                                          LibAv05.M2T_plotyylabelName(labelRange(iRef)), lString);
                end
                labelRef = M2T_join(m2t, labelRef, '');

                % Clear plotyy references. Ensures that references are created only once
                m2t.axes{end}.PlotyyReferences = [];
            else
                % Do nothing: it's gonna be a legend entry.
                % Not a label nor a referenced entry from the main axis.
            end
        end
        % ==============================================================================
        function label = M2T_plotyylabelName(num)
            % creates a LaTeX label for a plotyy trace
            label = sprintf('plotyyref:leg%d', num);
        end
        % ==============================================================================
        function legendInfo = M2T_addLegendInformation(m2t, h)
            % Add the actual legend string

            legendInfo = '';
            if ~m2t.currentHandleHasLegend
                return
            end
            legendString = M2T_getLegendString(m2t,h);

            % We also need a legend alignment option to make multiline
            % legend entries work. This is added by default in LibAv05.M2T_getLegendOpts().
            legendInfo = sprintf('\\addlegendentry{%s}\n\n', legendString);
        end
        % ==============================================================================
        function data = M2T_applyHgTransform(m2t, data)
            if ~isempty(m2t.transform)
                R = m2t.transform(1:3,1:3);
                t = m2t.transform(1:3,4);
                n = size(data, 1);
                data = data * R' + kron(ones(n,1), t');
            end
        end
        % ==============================================================================
        function m2t = M2T_drawAxes(m2t, handle)
            % Input arguments:
            %    handle.................The axes environment handle.

            LibAv05.M2T_assertRegularAxes(handle);

            % Initialize empty environment.
            % Use a struct instead of a custom subclass of hgsetget (which would
            % facilitate writing clean code) as structs are more portable (old MATLAB(R)
            % versions, GNU Octave).
            m2t.axes{end+1} = struct('handle',   handle, ...
                                   'name',     '', ...
                                   'comment',  [], ...
                                   'options',  {LibAv05.M2T_opts_new()}, ...
                                   'content',  {cell(0)}, ...
                                   'children', {cell(0)});

            % update gca
            m2t.current.gca = handle;

            % Check if axis is 3d
            % In MATLAB, all plots are treated as 3D plots; it's just the view that
            % makes 2D plots appear like 2D.
            m2t.axes{end}.is3D = isAxis3D(handle);

            % Flag if axis contains barplot
            m2t.axes{end}.barAddedAxisOption = false;

            % Get legend entries
            m2t.axes{end}.LegendHandle  = M2T_getAssociatedLegend(m2t, handle);
            m2t.axes{end}.LegendEntries = M2T_getLegendEntries(m2t);
            m2t = M2T_getPlotyyReferences(m2t, handle);

            m2t = M2T_retrievePositionOfAxes(m2t, handle);

            m2t = M2T_addAspectRatioOptionsOfAxes(m2t, handle);

            % Axis direction
            for axis = 'xyz'
                m2t.([axis 'AxisReversed']) = ...
                    strcmpi(get(handle,[upper(axis),'Dir']), 'reverse');
            end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Add color scaling
            CLimMode = get(handle,'CLimMode');
            if strcmpi(CLimMode,'manual') || ~isempty(m2t.cbarHandles)
                clim = caxis(handle);
                m2t = M2T_m2t_addAxisOption(m2t, 'point meta min', sprintf(m2t.ff, clim(1)));
                m2t = M2T_m2t_addAxisOption(m2t, 'point meta max', sprintf(m2t.ff, clim(2)));
            end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Recurse into the children of this environment.
            [m2t, childrenEnvs] = M2T_handleAllChildren(m2t, handle);
            m2t.axes{end} = M2T_addChildren(m2t.axes{end}, childrenEnvs);
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % The rest of this is handling axes options.
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Get other axis options (ticks, axis color, label,...).
            % This is set here such that the axis orientation indicator in m2t is set
            % before -- if ~LibAv05.M2T_isVisible(handle) -- the handle's children are called.
            [m2t, xopts] = M2T_getAxisOptions(m2t, handle, 'x');
            [m2t, yopts] = M2T_getAxisOptions(m2t, handle, 'y');

            m2t.axes{end}.options = M2T_opts_merge(m2t.axes{end}.options, xopts, yopts);

            m2t = M2T_add3DOptionsOfAxes(m2t, handle);

            if ~LibAv05.M2T_isVisible(handle)
                % Setting hide{x,y} axis also hides the axis labels in Pgfplots whereas
                % in MATLAB, they may still be visible. Instead use the following.
                m2t = M2T_m2t_addAxisOption(m2t, 'axis line style', '{draw=none}');
                m2t = M2T_m2t_addAxisOption(m2t, 'ticks', 'none');
                %    % An invisible axes container *can* have visible children, so don't
                %    % immediately bail out here.
                %    children = allchild(handle);
                %    for child = children(:)'
                %        if LibAv05.M2T_isVisible(child)
                %            % If the axes contain something that's visible, add an invisible
                %            % axes pair.
                %            m2t.axes{end}.name = 'axis';
                %            m2t.axes{end}.options = {m2t.axes{end}.options{:}, ...
                %                                               'hide x axis', 'hide y axis'};
                %            NOTE: getTag was removed in 76d260d12e615602653d6f7b357393242b2430b3
                %            m2t.axes{end}.comment = getTag(handle);
                %            break;
                %        end
                %    end
                %    % recurse into the children of this environment
                %    [m2t, childrenEnvs] = M2T_handleAllChildren(m2t, handle);
                %    m2t.axes{end} = M2T_addChildren(m2t.axes{end}, childrenEnvs);
                %    return
            end
            m2t.axes{end}.name = 'axis';

            m2t = M2T_drawBackgroundOfAxes(m2t, handle);
            m2t = M2T_drawTitleOfAxes(m2t, handle);
            m2t = M2T_drawBoxAndLineLocationsOfAxes(m2t, handle);
            m2t = M2T_drawGridOfAxes(m2t, handle);
            m2t = M2T_drawLegendOptionsOfAxes(m2t);

            m2t.axes{end}.options = LibAv05.M2T_LibAv05.M2T_opts_append_userdefined(m2t.axes{end}.options, ...
                                                          m2t.args.extraAxisOptions);
        end
        % ==============================================================================
        function m2t = M2T_drawGridOfAxes(m2t, handle)
            % Draws the grids of an axis
            options = M2T_opts_new();

            % Check for major/minor grids
            hasGrid = [LibAv05.M2T_isOn(get(handle, 'XGrid'));
                       LibAv05.M2T_isOn(get(handle, 'YGrid'));
                       LibAv05.M2T_isOn(get(handle, 'ZGrid')) && isAxis3D(handle)];

            hasMinorGrid = [LibAv05.M2T_isOn(get(handle, 'XMinorGrid'));
                            LibAv05.M2T_isOn(get(handle, 'YMinorGrid'));
                            LibAv05.M2T_isOn(get(handle, 'ZMinorGrid')) && isAxis3D(handle)];

            xyz = {'x', 'y', 'z'};

            % Check for local grid options
            % NOTE: for individual axis color options see the pfgmanual under
            % major x grid style
            for i=1:3
                if hasGrid(i)
                    grid    = [xyz{i}, 'majorgrids'];
                    options = M2T_opts_add(options, grid);
                end
                if hasMinorGrid(i)
                    grid    = [xyz{i}, 'minorgrids'];
                    options = M2T_opts_add(options, grid);
                end
            end

            % Check for global grid options
            if any(hasGrid)
                gridOpts = M2T_opts_new();
                % Get the line style and translate it to pgfplots
                [gridLS, isDefault] = M2T_getAndCheckDefault(...
                    'axes', handle, 'GridLineStyle', ':');
                if ~isDefault || m2t.args.strict
                    gridOpts = M2T_opts_add(gridOpts, LibAv05.M2T_translateLineStyle(gridLS));
                end

                % Get the color of the grid and translate it to pgfplots usable
                % values
                [gridColor, defaultColor] = M2T_getAndCheckDefault(...
                    'axes', handle, 'GridColor', [0.15, 0.15, 0.15]);
                if ~defaultColor
                    [m2t, gridColor] = M2T_getColor(m2t, handle, gridColor, 'patch');
                    gridOpts = M2T_opts_add(gridOpts, gridColor);
                end

                % Get the alpha of the grid and translate it to pgfplots
                [gridAlpha, defaultAlpha] = M2T_getAndCheckDefault(...
                    'axes', handle, 'GridAlpha', 0.1);
                if ~defaultAlpha
                    gridOpts = M2T_opts_add(gridOpts, 'opacity', num2str(gridAlpha));
                end

                if ~isempty(gridOpts)
                    options = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(options, 'grid style', gridOpts);
                end
            end

            if any(hasMinorGrid)
                minorGridOpts = M2T_opts_new();
                % Get the line style and translate it to pgfplots
                [minorGridLS, isDefault] = M2T_getAndCheckDefault(...
                    'axes', handle, 'MinorGridLineStyle', ':');
                if ~isDefault || m2t.args.strict
                    minorGridOpts = M2T_opts_add(minorGridOpts, LibAv05.M2T_translateLineStyle(minorGridLS));
                end

                % Get the color of the grid and translate it to pgfplots usable
                % values
                [minorGridColor, defaultColor] = M2T_getAndCheckDefault(...
                    'axes', handle, 'MinorGridColor', [0.1, 0.1, 0.1]);
                if ~defaultColor
                    [m2t, minorGridColor] = M2T_getColor(m2t, handle, minorGridColor, 'patch');
                    minorGridOpts = M2T_opts_add(minorGridOpts, minorGridColor);
                end

                % Get the alpha of the grid and translate it to pgfplots
                [minorGridAlpha, defaultAlpha] = M2T_getAndCheckDefault(...
                    'axes', handle, 'MinorGridAlpha', 0.1);
                if ~defaultAlpha
                    minorGridOpts = M2T_opts_add(minorGridOpts, 'opacity', num2str(minorGridAlpha));
                end

                if ~isempty(minorGridOpts)
                    options = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(options, 'minor grid style', minorGridOpts);
                end
            end

            if ~any(hasGrid) && ~any(hasMinorGrid)
                % When specifying 'axis on top', the axes stay above all graphs (which is
                % default MATLAB behavior), but so do the grids (which is not default
                % behavior).
                %TODO: use proper grid ordering
                if m2t.args.strict
                    options = M2T_opts_add(options, 'axis on top');
                end
                % FIXME: axis background, axis grid, main, axis ticks, axis lines, axis tick labels, axis descriptions, axis foreground
            end

            m2t.axes{end}.options = M2T_opts_merge(m2t.axes{end}.options, options);
        end
        % ==============================================================================
        function m2t = M2T_add3DOptionsOfAxes(m2t, handle)
            % adds 3D specific options of an axes object
            if isAxis3D(handle)
                [m2t, zopts]        = M2T_getAxisOptions(m2t, handle, 'z');
                m2t.axes{end}.options = M2T_opts_merge(m2t.axes{end}.options, zopts);

                VIEWFORMAT = ['{' m2t.ff '}{' m2t.ff '}'];
                m2t = M2T_m2t_addAxisOption(m2t, 'view', sprintf(VIEWFORMAT, get(handle, 'View')));
            end
        end
        % ==============================================================================
        function legendhandle = M2T_getAssociatedLegend(m2t, axisHandle)
            % Get legend handle associated with current axis

            legendhandle = [];
            env = getEnvironment();
            switch env
                case 'Octave'
                    % Make sure that m2t.legendHandles is a row vector.
                    for lhandle = m2t.legendHandles(:)'
                        ud = get(lhandle, 'UserData');
                        % Empty if no legend and multiple handles if plotyy
                        if ~isempty(ud) && any(axisHandle == ud.handle)
                            legendhandle = lhandle;
                            break
                        end
                    end
                case 'MATLAB'
                    legendhandle = legend(axisHandle);
            end

            % NOTE: there is a BUG in HG1 and Octave. Setting the box off sets the
            % legend visibility off too. We assume the legend is visible if it has
            % a visible child.
            isInvisibleHG2 = M2T_isHG2() && ~LibAv05.M2T_isVisible(legendhandle);
            isInvisibleHG1orOctave = (~LibAv05.M2T_isHG2() || strcmpi(env,'Octave')) &&...
                ~LibAv05.M2T_LibAv05.M2T_isVisibleContainer(legendhandle);

            % Do not return the handle if legend is invisible
            if isInvisibleHG1orOctave || isInvisibleHG2;
                legendhandle = [];
            end
        end
        % ==============================================================================
        function entries = M2T_getLegendEntries(m2t)
            % Retrieve the handles of the objects that have a legend entry

            % Non-composite objects are straightforward, e.g. line, and have the
            % legend entry at their same level, hence we return their handle.
            %
            % Hggroups behave differently depending on the environment and we might
            % return the handle to the hgroot or to one of its children:
            %   1) Matlab places the legend entry at the hgroot.
            %
            %      Usually, the decision to place the legend is either unchanged from
            %      the first call to LibAv05.M2T_handleAllChildrena(axis) or delegated to a
            %      specialized drawing routine, e.g. LibAv05.M2T_drawContour(), if the group has to
            %      be drawn atomically. In this case, the legend entry stays with the
            %      hgroot.
            %
            %      If the hggroup is a pure container like in a bodeplot, i.e. the
            %      `type` is not listed in LibAv05.M2T_drawHggroup(), a nested call to
            %      LibAv05.M2T_handleAllChildren(hgroot) follows. But, this second call cannot detect
            %      legend entries on the children. Hence, we pass down the legend entry
            %      from the hgroot to its first child.
            %
            %   2) Octave places the entry with one of the children of the hgroot.
            %      Hence, most of the hggroups are correctly dealt by a nested
            %      LibAv05.M2T_handleAllChildren() call which detects the entry on the child.
            %      However, when we can guess the type of hggroup with
            %      LibAv05.M2T_guessOctavePlotType(), the legend entry should be placed at the root
            %      level, hence we bubble it up from the child to the hgroot.

            entries = [];
            legendHandle = m2t.axes{end}.LegendHandle;

            if isempty(legendHandle)
                return
            end

            switch getEnvironment()
                case 'Octave'
                    % See set(hlegend, "deletefcn", {@deletelegend2, ca, [], [], t1, hplots}); in legend.m
                    delfun  = get(legendHandle,'deletefcn');
                    entries = delfun{6};

                    % Bubble-up legend entry properties from child to hggroup root
                    % for guessable objects
                    for ii = 1:numel(entries)
                        child = entries(ii);
                        anc   = ancestor(child,'hggroup');
                        if isempty(anc) % not an hggroup
                            continue
                        end
                        cl = M2T_guessOctavePlotType(anc);
                        if ~strcmpi(cl, 'unknown') % guessable hggroup, then bubble-up
                            legendString = get(child,'displayname');
                            set(anc,'displayname',legendString);
                            entries(ii) = anc;
                        end
                    end

                case 'MATLAB'
                    % Undocumented property (exists at least since 2008a)
                    entries = get(legendHandle,'PlotChildren');

                    % Take only the first child from a pure hggroup (e.g. bodeplots)
                    for ii = 1:numel(entries)
                        entry     = entries(ii);
                        % Note that class() is not supported in Octave
                        isHggroupClass = strcmpi(class(handle(entry)),'hggroup');
                        if isHggroupClass
                            children    = get(entry, 'Children');
                            firstChild  = children(1);
                            if isnumeric(firstChild)
                                firstChild = handle(firstChild);
                            end
                            % Inherits DisplayName from hggroup root
                            set(firstChild, 'DisplayName', get(entry, 'DisplayName'));
                            entries(ii) = firstChild;
                        end
                    end
            end
        end
        % ==============================================================================
        function m2t = M2T_getPlotyyReferences(m2t,axisHandle)
            % Retrieve references to legend entries of the main plotyy axis
            %
            % A plotyy plot has a main and a secondary axis. The legend is associated
            % with the main axis and hence m2t will only include the legend entries
            % that belong to the \axis[] that has a legend.
            %
            % One way to include the legend entries from the secondary axis (in the
            % same legend) is to first label the \addplot[] and then reference them.
            % See https://tex.stackexchange.com/questions/42697/42752#42752
            %
            % However, in .tex labels should come before they are referenced. Hence,
            % we actually label the legend entries from the main axis and swap the
            % legendhandle to the secondary axis.
            %
            % The legend will not be plotted with the main \axis[] and the labelled
            % legend entries will be skipped until the secondary axis. Then, they will
            % be listed before any legend entry from the secondary axis.

            % Retrieve legend handle
            if LibAv05.M2T_isAxisMain(axisHandle)
                legendHandle = m2t.axes{end}.LegendHandle;
            else
                legendHandle = M2T_getAssociatedLegend(m2t,LibAv05.M2T_getPlotyyPeer(axisHandle));
                m2t.axes{end}.LegendHandle = legendHandle;
            end

            % Not a plotyy axis or no legend
            if ~LibAv05.M2T_isAxisPlotyy(axisHandle) || isempty(legendHandle)
                m2t.axes{end}.PlotyyReferences = [];

            elseif LibAv05.M2T_isAxisMain(axisHandle)
                % Mark legend entries of the main axis for labelling
                legendEntries = m2t.axes{end}.LegendEntries;
                ancAxes       = ancestor(legendEntries,'axes');
                idx           = ismember([ancAxes{:}], axisHandle);
                m2t.axes{end}.PlotyyReferences = legendEntries(idx);

                % Ensure no legend is created on the main axis
                m2t.axes{end}.LegendHandle = [];
            else
                % Get legend entries associated to secondary plotyy axis. We can do
                % this because we took the legendhandle from the peer (main axis)
                legendEntries = M2T_getLegendEntries(m2t);
                ancAxes       = ancestor(legendEntries,'axes');
                if iscell(ancAxes)
                    ancAxes = [ancAxes{:}];
                end
                idx = ismember(double(ancAxes), axisHandle);
                m2t.axes{end}.LegendEntries = legendEntries(idx);

                % Recover referenced legend entries of the main axis
                m2t.axes{end}.PlotyyReferences = legendEntries(~idx);
            end
        end
        % ==============================================================================
        function bool = M2T_isAxisMain(h)
            % Check if it is the main axis e.g. in a plotyy plot

            if ~LibAv05.M2T_isAxisPlotyy(h)
                bool = true;
                return % an axis not constructed by plotyy is always(?) a main axis
            end

            % If it is a Plotyy axis
            switch getEnvironment()
                case 'Octave'
                    plotyyAxes = get(h, '__plotyy_axes__');
                    bool       = find(plotyyAxes == h) == 1;

                case 'MATLAB'
                    bool = ~isempty(getappdata(h, 'LegendPeerHandle'));
            end
        end
        % ==============================================================================
        function bool = M2T_isAxisPlotyy(h)
            % Check if handle is a plotyy axis

            switch getEnvironment()
                case 'Octave'
                    % Cannot LibAv05.M2T_test hidden property with isfield(), is always false
                    try
                        get(h, '__plotyy_axes__');
                        bool = true;
                    catch
                        bool = false;
                    end

                case 'MATLAB'
                    bool = ~isempty(getappdata(h, 'graphicsPlotyyPeer'));
            end
        end
        % ==============================================================================
        function peer = M2T_getPlotyyPeer(axisHandle)
            % Get the other axis coupled in plotyy plots

            switch getEnvironment()
                case 'Octave'
                    plotyyAxes = get(axisHandle, '__plotyy_axes__');
                    peer       = setdiff(plotyyAxes, axisHandle);

                case 'MATLAB'
                    peer = getappdata(axisHandle, 'graphicsPlotyyPeer');
            end
        end
        % ==============================================================================
        function legendString = M2T_getLegendString(m2t, h)
            % Retrieve the legend string for the given handle
            str         = M2T_getOrDefault(h, 'displayname', '');
            interpreter = get(m2t.axes{end}.LegendHandle,'interpreter');

            % HG1: autogenerated legend strings, i.e. data1,..., dataN, do not populate
            % the 'displayname' property. Go through 'userdata'
            if isempty(str)
                ud  = get(m2t.axes{end}.LegendHandle,'userdata');
                idx = ismember(ud.handles, h);
                str = ud.lstrings{idx};
            end

            % split string to cell, if newline character '\n' (ASCII 10) is present
            delimeter    = sprintf('\n');
            str          = regexp(str, delimeter, 'split');
            str          = M2T_prettyPrint(m2t, str, interpreter);
            legendString = M2T_join(m2t, str, '\\');
        end
        % ==============================================================================
        function [m2t, bool] = M2T_hasLegendEntry(m2t, h)
            % Check if the handle has a legend entry and track its legend status in m2t
            legendEntries = m2t.axes{end}.LegendEntries;
            if isnumeric(h)
                legendEntries = double(legendEntries);
            end

            % Should not have a legend reference
            bool = any(ismember(h, legendEntries)) && ~LibAv05.M2T_hasPlotyyReference(m2t,h);
            m2t.currentHandleHasLegend = bool;
        end
        % ==============================================================================
        function bool = M2T_hasPlotyyReference(m2t,h)
            % Check if the handle has a legend reference
            plotyyReferences = m2t.axes{end}.PlotyyReferences;
            if isnumeric(h)
                plotyyReferences = double(plotyyReferences);
            end

            bool = any(ismember(h, plotyyReferences));
        end
        % ==============================================================================
        function m2t = M2T_retrievePositionOfAxes(m2t, handle)
            % This retrieves the position of an axes and stores it into the m2t data
            % structure

            pos = M2T_getAxesPosition(m2t, handle, m2t.args.width, ...
                                  m2t.args.height, m2t.axesBoundingBox);
            % set the width
            if (~m2t.args.noSize)
                % optionally prevents setting the width and height of the axis
                m2t = M2T_setDimensionOfAxes(m2t, 'width',  pos.w);
                m2t = M2T_setDimensionOfAxes(m2t, 'height', pos.h);

                m2t = M2T_m2t_addAxisOption(m2t, 'at', ...
                        ['{(' LibAv05.M2T_formatDim(pos.x.value, pos.x.unit) ','...
                              LibAv05.M2T_formatDim(pos.y.value, pos.y.unit) ')}']);
                % the following is general MATLAB behavior:
                m2t = M2T_m2t_addAxisOption(m2t, 'scale only axis');
            end
        end
        % ==============================================================================
        function m2t = M2T_setDimensionOfAxes(m2t, widthOrHeight, dimension)
            % sets the dimension "name" of the current axes to the struct "dim"
            m2t = M2T_m2t_addAxisOption(m2t, widthOrHeight, ...
                    LibAv05.M2T_formatDim(dimension.value, dimension.unit));
        end
        % ==============================================================================
        function m2t = M2T_addAspectRatioOptionsOfAxes(m2t, handle)
            % Set manual aspect ratio for current axes
            % TODO: deal with 'axis image', 'axis square', etc. (#540)
            if strcmpi(get(handle, 'DataAspectRatioMode'), 'manual') ||...
               strcmpi(get(handle, 'PlotBoxAspectRatioMode'), 'manual')
                % we need to set the plot box aspect ratio
                if m2t.axes{end}.is3D
                    % Note: set 'plot box ratio' for 3D axes to avoid bug with
                    % 'scale mode = uniformly' (see #560)
                    aspectRatio = M2T_getPlotBoxAspectRatio(handle);
                    m2t = M2T_m2t_addAxisOption(m2t, 'plot box ratio', ...
                                      LibAv05.M2T_formatAspectRatio(m2t, aspectRatio));
                end
            end
        end
        % ==============================================================================
        function m2t = M2T_drawBackgroundOfAxes(m2t, handle)
            % draw the background color of the current axes
            backgroundColor = get(handle, 'Color');
            if ~LibAv05.M2T_isNone(backgroundColor) && LibAv05.M2T_isVisible(handle)
                [m2t, col] = M2T_getColor(m2t, handle, backgroundColor, 'patch');
                m2t = M2T_m2t_addAxisOption(m2t, 'axis background/.style', sprintf('{fill=%s}', col));
            end
        end
        % ==============================================================================
        function m2t = M2T_drawTitleOfAxes(m2t, handle)
            % processes the title of an axes object
            [m2t, m2t.axes{end}.options] = M2T_getTitle(m2t, handle, m2t.axes{end}.options);
        end
        % ==============================================================================
        function [m2t, opts] = M2T_getTitle(m2t, handle, opts)
            % gets the title and its markup from an axes/colorbar/...
            [m2t, opts] = LibAv05.M2T_LibAv05.M2T_getTitleOrLabel_(m2t, handle, opts, 'Title');
        end
        function [m2t, opts] = M2T_getLabel(m2t, handle, opts, tikzKeyword)
            % gets the label and its markup from an axes/colorbar/...
            [m2t, opts] = LibAv05.M2T_LibAv05.M2T_getTitleOrLabel_(m2t, handle, opts, 'Label', tikzKeyword);
        end
        function [m2t, opts] = M2T_getAxisLabel(m2t, handle, axis, opts)
            % convert an {x,y,z} axis label to TikZ
            labelName = [upper(axis) 'Label'];
            [m2t, opts] = LibAv05.M2T_LibAv05.M2T_getTitleOrLabel_(m2t, handle, opts, labelName);
        end
        function [m2t, opts] = M2T_getTitleOrLabel_(m2t, handle, opts, labelKind, tikzKeyword)
            % gets a string element from an object
            if ~exist('tikzKeyword', 'var') || isempty(tikzKeyword)
                tikzKeyword = lower(labelKind);
            end
            object = get(handle, labelKind);

            str = get(object, 'String');
            if ~isempty(str)
                interpreter = get(object, 'Interpreter');
                str = M2T_prettyPrint(m2t, str, interpreter);
                [m2t, style] = M2T_getFontStyle(m2t, object);
                if length(str) > 1 %multiline
                    style = M2T_opts_add(style, 'align', 'center');
                end
                if ~isempty(style)
                    opts = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(opts, [tikzKeyword ' style'], style);
                end
                str = M2T_join(m2t, str, '\\[1ex]');
                opts =  LibAv05.M2T_opts_add(opts, tikzKeyword, sprintf('{%s}', str));
            end
        end
        % ==============================================================================
        function m2t = M2T_drawBoxAndLineLocationsOfAxes(m2t, h)
            % draw the box and axis line location of an axes object
            isBoxOn       = M2T_isOn(get(h, 'box'));
            xLoc          = get(h, 'XAxisLocation');
            yLoc          = get(h, 'YAxisLocation');
            isXaxisBottom = strcmpi(xLoc,'bottom');
            isYaxisLeft   = strcmpi(yLoc,'left');

            % Only flip the labels to the other side if not at the default
            % left/bottom positions
            if isBoxOn
                if ~isXaxisBottom
                    m2t = M2T_m2t_addAxisOption(m2t, 'xticklabel pos','right');
                end
                if ~isYaxisLeft
                    m2t = M2T_m2t_addAxisOption(m2t, 'yticklabel pos','right');
                end

                % Position axes lines (strips the box)
            else
                m2t = M2T_m2t_addAxisOption(m2t, 'axis x line*', xLoc);
                m2t = M2T_m2t_addAxisOption(m2t, 'axis y line*', yLoc);
                if m2t.axes{end}.is3D
                    % There's no such attribute as 'ZAxisLocation'.
                    % Instead, the default seems to be 'left'.
                    m2t = M2T_m2t_addAxisOption(m2t, 'axis z line*', 'left');
                end
            end
        end
        % ==============================================================================
        function m2t = M2T_drawLegendOptionsOfAxes(m2t)
            legendHandle = m2t.axes{end}.LegendHandle;
            if isempty(legendHandle)
                return
            end

            [m2t, key, legendOpts] = M2T_getLegendOpts(m2t, legendHandle);
            m2t = M2T_m2t_addAxisOption(m2t, key, legendOpts);
        end
        % ==============================================================================
        function m2t = M2T_handleColorbar(m2t, handle)
            if isempty(handle)
                return;
            end

            % Find the axes environment that this colorbar belongs to.
            parentAxesHandle = double(get(handle,'axes'));
            parentFound = false;
            for k = 1:length(m2t.axes)
                if m2t.axes{k}.handle == parentAxesHandle
                    k0 = k;
                    parentFound = true;
                    break;
                end
            end
            if parentFound
                m2t.axes{k0}.options = M2T_opts_append(m2t.axes{k0}.options, ...
                    LibAv05.M2T_matlab2pgfplotsColormap(m2t, m2t.current.colormap), []);
                % Append cell string.
                m2t.axes{k0}.options = cat(1, m2t.axes{k0}.options, ...
                                            LibAv05.M2T_LibAv05.M2T_getColorbarOptions(m2t, handle));
            else
                warning('LibAv05.M2T_matlab2tikz:parentAxesOfColorBarNotFound',...
                        'Could not find parent axes for color bar. Skipping.');
            end
        end
        % ==============================================================================
        function [m2t, options] = M2T_getAxisOptions(m2t, handle, axis)
            LibAv05.M2T_assertValidAxisSpecifier(axis);

            options = M2T_opts_new();
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % axis colors
            [color, isDfltColor] = M2T_getAndCheckDefault('Axes', handle, ...
                                                      [upper(axis),'Color'], [ 0 0 0 ]);
            if ~isDfltColor || m2t.args.strict
                [m2t, col] = M2T_getColor(m2t, handle, color, 'patch');
                if LibAv05.M2T_isOn(get(handle, 'box'))
                    % If the axes are arranged as a box, make sure that the individual
                    % axes are drawn as four separate paths. This makes the alignment
                    % at the box corners somewhat less nice, but allows for different
                    % axis styles (e.g., colors).
                    options = M2T_opts_add(options, 'separate axis lines');
                end
                % set color of axis lines
                options = ...
                    LibAv05.M2T_opts_add(options, ...
                    ['every outer ', axis, ' axis line/.append style'], ...
                    ['{', col, '}']);
                % set color of tick labels
                options = ...
                    LibAv05.M2T_opts_add(options, ...
                    ['every ',axis,' tick label/.append style'], ...
                    ['{font=\color{',col,'}}']);
                % set color of ticks
                options = ...
                    LibAv05.M2T_opts_add(options, ...
                    ['every ',axis,' tick/.append style'], ...
                    ['{',col,'}']);
            end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % handle the orientation
            isAxisReversed = strcmpi(get(handle,[upper(axis),'Dir']), 'reverse');
            m2t.([axis 'AxisReversed']) = isAxisReversed;
            if isAxisReversed
                options = M2T_opts_add(options, [axis, ' dir'], 'reverse');
            end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            axisScale = M2T_getOrDefault(handle, [upper(axis) 'Scale'], 'lin');
            if strcmpi(axisScale, 'log');
                options = M2T_opts_add(options, [axis,'mode'], 'log');
            end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % get axis limits
            options = M2T_setAxisLimits(m2t, handle, axis, options);
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % get ticks along with the labels
            [options] = M2T_getAxisTicks(m2t, handle, axis, options);
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % get axis label
            [m2t, options] = M2T_getAxisLabel(m2t, handle, axis, options);
        end
        % ==============================================================================
        function [options] = M2T_getAxisTicks(m2t, handle, axis, options)
            % Return axis tick marks Pgfplots style. Nice: Tick lengths and such
            % details are taken care of by Pgfplots.
            LibAv05.M2T_assertValidAxisSpecifier(axis);

            keywordTickMode = [upper(axis), 'TickMode'];
            tickMode = get(handle, keywordTickMode);
            keywordTick = [upper(axis), 'Tick'];
            ticks = get(handle, keywordTick);

            % hidden properties are not caught by LibAv05.M2T_hasProperties
            isDatetimeTicks = M2T_isAxisTicksDateTime(handle, axis);

            if isempty(ticks)
                % If no ticks are present, we need to enforce this in any case.
                pgfTicks = '\empty';
            elseif strcmpi(tickMode, 'auto') && ~m2t.args.strict && ~isDatetimeTicks
                % Let pgfplots decide if the tickmode is auto or conversion is not
                % strict and we are not dealing with datetime ticks
                pgfTicks = [];
            else % strcmpi(tickMode,'manual') || m2t.args.strict
                pgfTicks = M2T_join(m2t, cellstr(num2str(ticks(:))), ', ');
            end

            keywordTickLabelMode = [upper(axis), 'TickLabelMode'];
            tickLabelMode = get(handle, keywordTickLabelMode);
            if strcmpi(tickLabelMode, 'auto') && ~m2t.args.strict && ~isDatetimeTicks
                pgfTickLabels = [];
            else % strcmpi(tickLabelMode,'manual') || m2t.args.strict
                % HG2 allows to set 'TickLabelInterpreter'.
                % HG1 tacitly uses the interpreter 'none'.
                % See http://www.mathworks.com/matlabcentral/answers/102053#comment_300079
                fallback    = M2T_defaultTickLabelInterpreter(m2t);
                interpreter = M2T_getOrDefault(handle, 'TickLabelInterpreter', fallback);
                keywordTickLabel = [upper(axis), 'TickLabel'];
                tickLabels = cellstr(get(handle, keywordTickLabel));
                tickLabels = M2T_prettyPrint(m2t, tickLabels, interpreter);

                keywordScale = [upper(axis), 'Scale'];
                isAxisLog = strcmpi(LibAv05.M2T_getOrDefault(handle,keywordScale, 'lin'), 'log');
                [pgfTicks, pgfTickLabels] = ...
                    LibAv05.M2T_matlabTicks2pgfplotsTicks(m2t, ticks, tickLabels, isAxisLog, tickLabelMode);
            end

            keywordMinorTick = [upper(axis), 'MinorTick'];
            hasMinorTicks = M2T_isOn(LibAv05.M2T_getOrDefault(handle, keywordMinorTick, 'off'));
            tickDirection = M2T_getOrDefault(handle, 'TickDir', 'in');

            options = M2T_setAxisTicks(m2t, options, axis, pgfTicks, pgfTickLabels, ...
                hasMinorTicks, tickDirection, isDatetimeTicks);

            options = M2T_setAxisTickLabelStyle(options, axis, handle);
        end
        % ==============================================================================
        function options = M2T_setAxisTickLabelStyle(options, axis, handle)
            % determine the style of tick labels
            %TODO: translate the style of tick labels fully (font?, weight, ...)
            kwRotation = [upper(axis), 'TickLabelRotation'];
            rotation = M2T_getOrDefault(handle, kwRotation, 0);
            if rotation ~= 0
                options = M2T_opts_add(options, [axis, 'ticklabel style'], ...
                                            sprintf('{rotate=%d}', rotation));
            end
        end
        % ==============================================================================
        function interpreter = M2T_defaultTickLabelInterpreter(m2t)
            % determines the default tick label interpreter
            % This is only relevant in HG1/Octave. In HG2, we use the interpreter
            % set in the object (not the global default).
            if m2t.args.interpretTickLabelsAsTex
                interpreter = 'tex';
            else
                interpreter = 'none';
            end
        end
        % ==============================================================================
        function isDatetimeTicks = M2T_isAxisTicksDateTime(handle, axis)
            % returns true when the axis has DateTime ticks
            try
                % Get hidden properties of the datetime axes manager
                dtsManager = get(handle, 'DatetimeDurationPlotAxesListenersManager');
                oldState   = warning('off','MATLAB:structOnObject');
                dtsManager = struct(dtsManager);
                warning(oldState);

                isDatetimeTicks = dtsManager.([upper(axis) 'DateTicks']) == 1;
            catch
                isDatetimeTicks = false;
            end
        end
        % ==============================================================================
        function options = M2T_setAxisTicks(m2t, options, axis, ticks, tickLabels,hasMinorTicks, tickDir,isDatetimeTicks)
            % set ticks options

            % According to http://www.mathworks.com/help/techdoc/ref/axes_props.html,
            % the number of minor ticks is automatically determined by MATLAB(R) to
            % fit the size of the axis. Until we know how to extract this number, use
            % a reasonable default.
            matlabDefaultNumMinorTicks = 3;
            if ~isempty(ticks)
                options = M2T_opts_add(options, [axis,'tick'], sprintf('{%s}', ticks));
            end
            if ~isempty(tickLabels)
                options = M2T_opts_add(options, ...
                    [axis,'ticklabels'], sprintf('{%s}', tickLabels));
            end
            if hasMinorTicks
                options = M2T_opts_add(options, [axis,'minorticks'], 'true');
                if m2t.args.strict
                    options = M2T_opts_add(options, ...
                        sprintf('minor %s tick num', axis), ...
                        sprintf('{%d}', matlabDefaultNumMinorTicks));
                end
            end

            if strcmpi(tickDir,'out')
                options = M2T_opts_add(options, 'tick align', 'outside');
            elseif strcmpi(tickDir,'both')
                options = M2T_opts_add(options, 'tick align', 'center');
            end

            if isDatetimeTicks
                options = M2T_opts_add(options, ['scaled ' axis ' ticks'], 'false');
            end
        end
        % ==============================================================================
        function M2T_assertValidAxisSpecifier(axis)
            % assert that axis is a valid axis specifier
            if ~ismember(axis, {'x','y','z'})
                error('LibAv05.M2T_matlab2tikz:illegalAxisSpecifier', ...
                      'Illegal axis specifier "%s".', axis);
            end
        end
        % ==============================================================================
        function M2T_assertRegularAxes(handle)
            % assert that the (axes) object specified by handle is a regular axes and not a
            % colorbar or a legend
            tag = lower(get(handle,'Tag'));
            if ismember(tag,{'colorbar','legend'})
                error('LibAv05.M2T_matlab2tikz:notARegularAxes', ...
                      ['The object "%s" is not a regular axes object. ' ...
                       'It cannot be handled with LibAv05.M2T_drawAxes!'], handle);
            end
        end
        % ==============================================================================
        function options = M2T_setAxisLimits(m2t, handle, axis, options)
            % set the upper/lower limit of an axis
            limits = get(handle, [upper(axis),'Lim']);
            if isfinite(limits(1))
                options = M2T_opts_add(options, [axis,'min'], sprintf(m2t.ff, limits(1)));
            end
            if isfinite(limits(2))
                options = M2T_opts_add(options, [axis,'max'], sprintf(m2t.ff, limits(2)));
            end
        end
        % ==============================================================================
        function bool = M2T_isVisibleContainer(axisHandle)
            if ~LibAv05.M2T_isVisible(axisHandle)
                % An invisible axes container *can* have visible children, so don't
                % immediately bail out here. Also it *can* have a visible title,
                % labels or children

                bool = false;
                for prop = {'Children', 'Title', 'xlabel', 'ylabel', 'ZLabel'}
                    property = prop{1};
                    if strcmpi(property, 'Children')
                        children = allchild(axisHandle);
                    elseif isprop(axisHandle, property)
                        children = get(axisHandle, property);
                    else
                        continue; % don't check non-existent properties
                    end
                    for child = children(:)'
                        if LibAv05.M2T_isVisible(child)
                            bool = true;
                            return;
                        end
                    end
                end
            else
                bool = true;
            end
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawLine(m2t, h)
            % Returns the code for drawing a regular line and error bars.
            % This is an extremely common operation and takes place in most of the
            % not too fancy plots.
            str = '';
            if ~LibAv05.M2T_isLineVisible(h)
                return; % there is nothing to plot
            end

            % Color
            color         = get(h, 'Color');
            [m2t, xcolor] = M2T_getColor(m2t, h, color, 'patch');
            % Line and marker options
            [m2t, lineOptions]   = M2T_getLineOptions(m2t, h);
            [m2t, markerOptions] = M2T_getMarkerOptions(m2t, h);

            drawOptions = M2T_opts_new();
            drawOptions = M2T_opts_add(drawOptions, 'color', xcolor);
            drawOptions = M2T_opts_merge(drawOptions, lineOptions, markerOptions);

            % Check for "special" lines, e.g.:
            if strcmpi(get(h, 'Tag'), 'zplane_unitcircle')
                [m2t, str] = M2T_specialDrawZplaneUnitCircle(m2t, drawOptions);
                return
            end

            % build the data matrix
            data       = M2T_getXYZDataFromLine(m2t, h);
            yDeviation = M2T_getYDeviations(h);
            if ~isempty(yDeviation)
                data = [data, yDeviation];
            end

            % Check if any value is infinite/NaN. In that case, add appropriate option.
            m2t = M2T_jumpAtUnboundCoords(m2t, data);

            [m2t, dataString]  = M2T_writePlotData(m2t, data, drawOptions);
            [m2t, labelString] = M2T_addLabel(m2t, h);

            str = [dataString, labelString];
        end
        % ==============================================================================
        function [m2t, str] = M2T_specialDrawZplaneUnitCircle(m2t, drawOptions)
            % Draw unit circle and axes.

            % TODO Don't hardcode "10", but extract from parent axes of |h|
            opts = M2T_opts_print(drawOptions);
            str  = [sprintf('\\draw[%s] (axis cs:0,0) circle[radius=1];\n',  opts), ...
                    sprintf('\\draw[%s] (axis cs:-10,0)--(axis cs:10,0);\n', opts), ...
                    sprintf('\\draw[%s] (axis cs:0,-10)--(axis cs:0,10);\n', opts)];
        end
        % ==============================================================================
        function bool = M2T_isLineVisible(h)
            % check if a line object is actually visible (has markers and so on)

            lineStyle     = get(h, 'LineStyle');
            lineWidth     = get(h, 'LineWidth');
            marker        = M2T_getOrDefault(h, 'Marker','none');
            hasLines      = ~LibAv05.M2T_isNone(lineStyle) && lineWidth > 0;
            hasMarkers    = ~LibAv05.M2T_isNone(marker);
            hasDeviations = ~isempty(LibAv05.M2T_getYDeviations(h));

            bool = M2T_isVisible(h) && (hasLines || hasMarkers || hasDeviations);
        end
        % ==============================================================================
        function [m2t, str] = M2T_writePlotData(m2t, data, drawOptions)
            % actually writes the plot data to file
            str = '';

            is3D = m2t.axes{end}.is3D;
            if is3D
                % Don't try to be smart in parametric 3d plots: Just plot all the data.
                [m2t, table, tableOptions] = M2T_makeTable(m2t, {'','',''}, data);

                % Print out
                drawOpts = M2T_opts_print(drawOptions);
                tabOpts  = M2T_opts_print(tableOptions);
                str      = sprintf('\\addplot3 [%s]\n table[%s] {%s};\n ', ...
                                   drawOpts, tabOpts, table);
            else
                % split the data into logical chunks
                dataCell = M2T_splitLine(m2t, data);

                % plot them
                strPart = cell(1, length(dataCell));
                for k = 1:length(dataCell)
                    % If the line has a legend string, make sure to only include a legend
                    % entry for the *last* occurrence of the plot series.
                    % Hence the condition k<length(xDataCell).
                    %if ~isempty(m2t.legendHandles) && (~m2t.currentHandleHasLegend || k < length(dataCell))
                    if ~m2t.currentHandleHasLegend || k < length(dataCell)
                        % No legend entry found. Don't include plot in legend.
                        hiddenDrawOptions = M2T_maybeShowInLegend(false, drawOptions);
                        opts = M2T_opts_print(hiddenDrawOptions);
                    else
                        opts = M2T_opts_print(drawOptions);
                    end

                    [m2t, Part] = M2T_plotLine2d(m2t, opts, dataCell{k});
                    strPart{k} = Part;
                end
                strPart = M2T_join(m2t, strPart, '');
                str = [str, strPart];
            end
        end
        % ==============================================================================
        function [data] = M2T_getXYZDataFromLine(m2t, h)
            % Retrieves the X, Y and Z (if appropriate) data from a Line object
            %
            % First put them all together in one multiarray.
            % This also implicitly makes sure that the lengths match.
            try
                xData = get(h, 'XData');
                yData = get(h, 'YData');
            catch
                % Line annotation
                xData = get(h, 'X');
                yData = get(h, 'Y');
            end
            is3D  = m2t.axes{end}.is3D;
            if ~is3D
                data = [xData(:), yData(:)];
            else
                zData = get(h, 'ZData');
                data = M2T_applyHgTransform(m2t, [xData(:), yData(:), zData(:)]);
            end
        end
        % ==============================================================================
        function [m2t, labelCode] = M2T_addLabel(m2t, h)
            % conditionally add a LaTeX label after the current plot
            labelCode = '';

            if m2t.args.automaticLabels||m2t.args.LibAv05.M2T_addLabels
                lineTag = get(h,'Tag');
                if ~isempty(lineTag)
                    labelName = sprintf('%s', lineTag);
                else
                    [pathstr, name] = fileparts(m2t.args.filename); %#ok
                    labelName = sprintf('addplot:%s%d', name, m2t.count.autolabel);
                    [m2t] = M2T_incrementGlobalCounter(m2t, 'autolabel');
                    % TODO: First increment the counter, then use it such that the
                    % pattern is the same everywhere
                end
                labelCode = sprintf('\\label{%s}\n', labelName);
                LibAv05.M2T_userWarning(m2t, 'Automatically added label ''%s'' for line plot.', labelName);
            end
        end
        % ==============================================================================
        function [m2t,str] = M2T_plotLine2d(m2t, opts, data)
            errorbarMode = (size(data,2) == 4); % is (optional) yDeviation given?

            errorBar = '';
            if errorbarMode
                m2t      = M2T_needsPgfplotsVersion(m2t, [1,9]);
                errorBar = sprintf('plot [error bars/.cd, y dir = both, y explicit]\n');
            end

            % Convert to string array then cell to call sprintf once (and no loops).
            [m2t, table, tableOptions] = M2T_makeTable(m2t, repmat({''}, size(data,2)), data);
            if errorbarMode
                tableOptions = M2T_opts_add(tableOptions, 'y error plus index', '2');
                tableOptions = M2T_opts_add(tableOptions, 'y error minus index', '3');
            end

            % Print out
            tabOpts = M2T_opts_print(tableOptions);
            str     = sprintf('\\addplot [%s]\n %s table[%s]{%s};\n',...
                              opts, errorBar, tabOpts, table);
        end
        % ==============================================================================
        function dataCell = M2T_splitLine(m2t, data)
            % TeX parses files line by line with a buffer of size buf_size. If the
            % plot has too many data points, pdfTeX's buffer size may be exceeded.
            % As a work-around, split the xData, yData into several chunks of data
            % for each of which an \addplot will be generated.

            % Get the length of the data array and the corresponding chung size
            %TODO: scale `maxChunkLength` with the number of columns in the data array
            len         = size(data, 1);
            chunkLength = m2t.args.maxChunkLength;
            chunks      = chunkLength * ones(ceil(len/chunkLength), 1);
            if mod(len, chunkLength) ~=0
                chunks(end) = mod(len, chunkLength);
            end

            % Cut the data into chunks
            dataCell = mat2cell(data, chunks);

            % Add an extra (overlap) point to the data stream otherwise the line
            % between two data chunks would be broken. Technically, this is only
            % needed when the plot has a line connecting the points, but the
            % additional cost when there is no line doesn't justify the added
            % complexity.
            for i=1:length(dataCell)-1
                dataCell{i}(end+1,:) = dataCell{i+1}(1,:);
            end
        end
        % ==============================================================================
        function [m2t, lineOpts] = M2T_getLineOptions(m2t, h)
            % Gathers the line options.
            lineOpts = M2T_opts_new();

            % Get the options from the handle
            lineWidth = get(h, 'LineWidth');

            % Get the line style and check whether it is the default one
            [lineStyle, isDefaultLS] = M2T_getAndCheckDefault('Line', h, 'LineStyle', '-');

            if ~isDefaultLS && ~LibAv05.M2T_isNone(lineStyle) && (lineWidth > m2t.tol)
                lineOpts = M2T_opts_add(lineOpts, LibAv05.M2T_translateLineStyle(lineStyle));
            end

            % Take over the line width in any case when in strict mode. If not, don't add
            % anything in case of default line width and effectively take Pgfplots'
            % default.
            % Also apply the line width if no actual line is there; the markers make use
            % of this, too.
            matlabDefaultLineWidth = 0.5;
            if ~isempty(m2t.semantic.LineWidth)
                if ismember(lineWidth, [m2t.semantic.LineWidth{:,2}])
                    semStrID = lineWidth == [m2t.semantic.LineWidth{:,2}];
                    lineOpts = M2T_opts_add(lineOpts, m2t.semantic.LineWidth{semStrID,1});
                else
                    warning('LibAv05.M2T_matlab2tikz:semanticLineWidthNotFound',...
                        ['No semantic correspondance for lineWidth of ''%f'' found.'...
                        'Falling back to explicit export in points.'], lineWidth);
                    lineOpts = M2T_opts_add(lineOpts, 'line width', sprintf('%.1fpt', lineWidth));
                end
            elseif m2t.args.strict || ~abs(lineWidth-matlabDefaultLineWidth) <= m2t.tol
                lineOpts = M2T_opts_add(lineOpts, 'line width', sprintf('%.1fpt', lineWidth));
            end

            % print no lines
            if LibAv05.M2T_isNone(lineStyle) || lineWidth==0
                lineOpts = M2T_opts_add(lineOpts, 'draw', 'none');
            end
        end
        % ==============================================================================
        function list = M2T_configureSemanticLineWidths(semanticLineWidths)
            % Defines the default semantic options of pgfplots and updates it when applicable

            if isnan(semanticLineWidths)
                % Remove the list
                list = {};
                return;
            end

            % Pgf/TikZ defaults (see pgfmanual 3.0.1a section 15.3.1 / page 166)
            list = {'ultra thin',  0.1;
                    'very thin',   0.2;
                    'thin',        0.4;
                    'semithick',   0.6;
                    'thick',       0.8;
                    'very thick',  1.2;
                    'ultra thick', 1.6 };

            % Update defaults or append the user provided setting
            for ii = 1:size(semanticLineWidths, 1)
                % Check for redefinitions of defaults
                [isOverride, idx] = ismember(semanticLineWidths{ii, 1}, list{:, 1})
                if isOverride
                    list{idx, 2} = semanticLineWidths{ii, 2};
                else
                    list{end+1} = semanticLineWidths{ii, :};
                end
            end
        end
        % ==============================================================================
        function [m2t, drawOptions] = M2T_getMarkerOptions(m2t, h)
            % Handles the marker properties of a line (or any other) plot.
            drawOptions = M2T_opts_new();

            marker = M2T_getOrDefault(h, 'Marker', 'none');

            if ~LibAv05.M2T_isNone(marker)
                markerSize = get(h, 'MarkerSize');
                lineStyle  = get(h, 'LineStyle');
                lineWidth  = get(h, 'LineWidth');

                [tikzMarkerSize, isDefault] = ...
                    LibAv05.M2T_LibAv05.M2T_translateMarkerSize(m2t, marker, markerSize);

                % take over the marker size in any case when in strict mode;
                % if not, don't add anything in case of default marker size
                % and effectively take Pgfplots' default.
                if m2t.args.strict || ~isDefault
                    drawOptions = M2T_opts_add(drawOptions, 'mark size', ...
                                           sprintf('%.1fpt', tikzMarkerSize));
                end

                markOptions = M2T_opts_new();
                % make sure that the markers get painted in solid (and not dashed)
                % if the 'lineStyle' is not solid (otherwise there is no problem)
                if ~strcmpi(lineStyle, 'solid')
                    markOptions = M2T_opts_add(markOptions, 'solid');
                end

                % get the marker color right
                markerInfo = M2T_getMarkerInfo(m2t, h, markOptions);

                [m2t, markerInfo.options] = M2T_setColor(m2t, h, markerInfo.options, 'fill', markerInfo.FaceColor);

                if ~strcmpi(markerInfo.EdgeColor,'auto')
                    [m2t, markerInfo.options] = M2T_setColor(m2t, h, markerInfo.options, '', markerInfo.EdgeColor);
                else
                    if isprop(h,'EdgeColor')
                        color = get(h, 'EdgeColor');
                    else
                        color = get(h, 'Color');
                    end
                    [m2t, markerInfo.options] = M2T_setColor(m2t, h, markerInfo.options, '', color);
                end

                % add it all to drawOptions
                drawOptions = M2T_opts_add(drawOptions, 'mark', markerInfo.tikz);

                if ~isempty(markOptions)
                    drawOptions = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(drawOptions, 'mark options', ...
                                               markerInfo.options);
                end
            end
        end
        % ==============================================================================
        function [tikzMarkerSize, isDefault] = ...
            M2T_translateMarkerSize(m2t, matlabMarker, matlabMarkerSize)
            % The markersizes of Matlab and TikZ are related, but not equal. This
            % is because
            %
            %  1.) MATLAB uses the MarkerSize property to describe something like
            %      the diameter of the mark, while TikZ refers to the 'radius',
            %  2.) MATLAB and TikZ take different measures (e.g. the
            %      edge of a square vs. its diagonal).
            if(~ischar(matlabMarker))
                error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_LibAv05.M2T_translateMarkerSize',                      ...
                    'Variable matlabMarker is not a string.');
            end

            if(~isnumeric(matlabMarkerSize))
                error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_LibAv05.M2T_translateMarkerSize',                      ...
                    'Variable matlabMarkerSize is not a numeral.');
            end

            % 6pt is the default MATLAB marker size for all markers
            defaultMatlabMarkerSize = 6;
            isDefault = abs(matlabMarkerSize(1)-defaultMatlabMarkerSize)<m2t.tol;
            % matlabMarkerSize can be vector data, use first index to check the default
            % marker size. When the script also handles different markers together with
            % changing size and color, the LibAv05.M2T_test should be extended to a vector norm, e.g.
            % sqrt(e^T*e) < tol, where e=matlabMarkerSize-defaultMatlabMarkerSize

            switch (matlabMarker)
                case 'none'
                    tikzMarkerSize = [];
                case {'+','o','x','*','p','pentagram','h','hexagram'}
                    % In MATLAB, the marker size refers to the edge length of a
                    % square (for example) (~diameter), whereas in TikZ the
                    % distance of an edge to the center is the measure (~radius).
                    % Hence divide by 2.
                    tikzMarkerSize = matlabMarkerSize(:) / 2;
                case '.'
                    % as documented on the Matlab help pages:
                    %
                    % Note that MATLAB draws the point marker (specified by the '.'
                    % symbol) at one-third the specified size.
                    % The point (.) marker type does not change size when the
                    % specified value is less than 5.
                    %
                    tikzMarkerSize = matlabMarkerSize(:) / 2 / 3;
                case {'s','square'}
                    % Matlab measures the diameter, TikZ half the edge length
                    tikzMarkerSize = matlabMarkerSize(:) / 2 / sqrt(2);
                case {'d','diamond'}
                    % MATLAB measures the width, TikZ the height of the diamond;
                    % the acute angle (at the top and the bottom of the diamond)
                    % is a manually measured 75 degrees (in TikZ, and MATLAB
                    % probably very similar); use this as a base for calculations
                    tikzMarkerSize = matlabMarkerSize(:) / 2 / atan(75/2 *pi/180);
                case {'^','v','<','>'}
                    % for triangles, matlab takes the height
                    % and tikz the circumcircle radius;
                    % the triangles are always equiangular
                    tikzMarkerSize = matlabMarkerSize(:) / 2 * (2/3);
                otherwise
                    error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_LibAv05.M2T_translateMarkerSize',                   ...
                        'Unknown matlabMarker ''%s''.', matlabMarker);
            end
        end
        % ==============================================================================
        function [tikzMarker, markOptions] = ...
            M2T_translateMarker(m2t, matlabMarker, markOptions, faceColorToggle)
            % Translates MATLAB markers to their Tikz equivalents
            % #COMPLEX: inherently large switch-case
            if ~ischar(matlabMarker)
                error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_translateMarker:MarkerNotAString',...
                    'matlabMarker is not a string.');
            end

            switch (matlabMarker)
                case 'none'
                    tikzMarker = '';
                case '+'
                    tikzMarker = '+';
                case 'o'
                    if faceColorToggle
                        tikzMarker = '*';
                    else
                        tikzMarker = 'o';
                    end
                case '.'
                    tikzMarker = '*';
                case 'x'
                    tikzMarker = 'x';
                otherwise  % the following markers are only available with PGF's
                    % plotmarks library
                    LibAv05.M2T_signalDependency(m2t, 'tikzlibrary', 'plotmarks');
                    hasFilledVariant = true;
                    switch (matlabMarker)

                        case '*'
                            tikzMarker = 'asterisk';
                            hasFilledVariant = false;

                        case {'s','square'}
                            tikzMarker = 'square';

                        case {'d','diamond'}
                            tikzMarker = 'diamond';

                        case '^'
                            tikzMarker = 'triangle';

                        case 'v'
                            tikzMarker = 'triangle';
                            markOptions = M2T_opts_add(markOptions, 'rotate', '180');

                        case '<'
                            tikzMarker = 'triangle';
                            markOptions = M2T_opts_add(markOptions, 'rotate', '90');

                        case '>'
                            tikzMarker = 'triangle';
                            markOptions = M2T_opts_add(markOptions, 'rotate', '270');

                        case {'p','pentagram'}
                            tikzMarker = 'star';

                        case {'h','hexagram'}
                            LibAv05.M2T_userWarning(m2t, 'MATLAB''s marker ''hexagram'' not available in TikZ. Replacing by ''star''.');
                            tikzMarker = 'star';

                        otherwise
                            error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_translateMarker:unknownMatlabMarker',...
                                'Unknown matlabMarker ''%s''.',matlabMarker);
                    end
                    if faceColorToggle && hasFilledVariant
                        tikzMarker = [tikzMarker '*'];
                    end
            end
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawPatch(m2t, handle)
            % Draws a 'patch' graphics object (as found in contourf plots, for example).
            %
            str = '';
            if ~LibAv05.M2T_isVisible(handle)
                return
            end

            % This is for a quirky workaround for stacked bar plots.
            m2t.axes{end}.nonbarPlotsPresent = true;

            % Each row of the faces matrix represents a distinct patch
            % NOTE: pgfplot uses zero-based indexing into vertices and interpolates
            % counter-clockwise
            Faces    = get(handle,'Faces')-1;
            Vertices = get(handle,'Vertices');

            % 3D vs 2D
            is3D = m2t.axes{end}.is3D;
            if is3D
                columnNames = {'x', 'y', 'z'};
                plotCmd     = 'addplot3';
                Vertices    = M2T_applyHgTransform(m2t, Vertices);
            else
                columnNames = {'x', 'y'};
                plotCmd     = 'addplot';
                Vertices    = Vertices(:,1:2);
            end

            % Process fill, edge colors and shader
            [m2t,patchOptions, s] = M2T_shaderOpts(m2t,handle,'patch');

            % Return empty axes if no face or edge colors
            if LibAv05.M2T_isNone(s.plotType)
                return
            end

            % -----------------------------------------------------------------------
            % gather the draw options
            % Make sure that legends are shown in area mode.
            drawOptions = M2T_opts_add(LibAv05.M2T_opts_new,'area legend');
            verticesTableOptions = M2T_opts_new();

            % Marker options
            [m2t, markerOptions] = M2T_getMarkerOptions(m2t, handle);
            drawOptions          = M2T_opts_merge(drawOptions, markerOptions);

            % Line options
            [m2t, lineOptions] = M2T_getLineOptions(m2t, handle);
            drawOptions = M2T_opts_merge(drawOptions, lineOptions);

            % If the line is not visible, set edgeColor to none. Otherwise pgfplots
            % draws it by default
            if ~LibAv05.M2T_isLineVisible(handle)
                s.edgeColor = 'none';
            end

            % No patch: if one patch and single face/edge color
            isFaceColorFlat = isempty(strfind(LibAv05.M2T_opts_get(patchOptions, 'shader'),'interp'));
            if size(Faces,1) == 1 && s.hasOneEdgeColor && isFaceColorFlat
                ptType = '';
                cycle  = M2T_conditionallyCyclePath(Vertices);

                [m2t, drawOptions] = M2T_setColor(m2t, handle, drawOptions, 'draw', ...
                                                 s.edgeColor, 'none');
                [m2t, drawOptions] = M2T_setColor(m2t, handle, drawOptions, 'fill', ...
                                                 s.faceColor);

                [drawOptions] = M2T_opts_copy(patchOptions, 'draw opacity', drawOptions);
                [drawOptions] = M2T_opts_copy(patchOptions, 'fill opacity', drawOptions);

            else % Multiple patches

                % Patch table type
                ptType      = 'patch table';
                cycle       = '';
                drawOptions = M2T_opts_add(drawOptions,'table/row sep','crcr');
                % TODO: is the above "crcr" compatible with pgfplots 1.12 ?
                % TODO: is a "patch table" externalizable?

                % Enforce 'patch' or cannot use 'patch table='
                if strcmpi(s.plotType,'mesh')
                    drawOptions = M2T_opts_add(drawOptions,'patch');
                end
                drawOptions = M2T_opts_add(drawOptions,s.plotType); % Eventually add mesh, but after patch!

                drawOptions = M2T_getPatchShape(m2t, handle, drawOptions, patchOptions);

                [m2t, drawOptions, Vertices, Faces, verticesTableOptions, ptType, ...
                 columnNames] = LibAv05.M2T_LibAv05.M2T_setColorsOfPatches(m2t, handle, drawOptions, ...
                   Vertices, Faces, verticesTableOptions, ptType, columnNames, ...
                   isFaceColorFlat, s);
            end

            drawOptions = M2T_maybeShowInLegend(m2t.currentHandleHasLegend, drawOptions);
            m2t = M2T_jumpAtUnboundCoords(m2t, Faces(:));

            % Add Faces table
            if ~isempty(ptType)
                [m2t, facesTable] = M2T_makeTable(m2t, repmat({''},1,size(Faces,2)), Faces);
                drawOptions = M2T_opts_add(drawOptions, ptType, sprintf('{%s}', facesTable));
            end

            % Plot the actual data.
            [m2t, verticesTable, tableOptions] = M2T_makeTable(m2t, columnNames, Vertices);
            tableOptions = M2T_opts_merge(tableOptions, verticesTableOptions);

            % Print out
            drawOpts = M2T_opts_print(drawOptions);
            tabOpts  = M2T_opts_print(tableOptions);
            str = sprintf('\n\\%s[%s]\ntable[%s] {%s}%s;\n',...
                          plotCmd, drawOpts, tabOpts, verticesTable, cycle);
        end
        % ==============================================================================
        function [m2t, drawOptions, Vertices, Faces, verticesTableOptions, ptType, ...
                 columnNames] = M2T_setColorsOfPatches(m2t, handle, drawOptions, ...
                   Vertices, Faces, verticesTableOptions, ptType, columnNames, isFaceColorFlat, s)
            % this behemoth does the color setting for patches

            % TODO: this function can probably be split further, just look at all those
            % parameters being passed.

            fvCData   = get(handle,'FaceVertexCData');
            rowsCData = size(fvCData,1);

            % We have CData for either all faces or vertices
            if rowsCData > 1

                % Add the color map
                m2t = M2T_m2t_addAxisOption(m2t, LibAv05.M2T_matlab2pgfplotsColormap(m2t, m2t.current.colormap));

                % Determine if mapping is direct or scaled
                CDataMapping = get(handle,'CDataMapping');
                if strcmpi(CDataMapping, 'direct')
                    drawOptions = M2T_opts_add(drawOptions, 'colormap access','direct');
                end

                % Switch to face CData if not using interpolated shader
                isVerticesCData = rowsCData == size(Vertices,1);
                if isFaceColorFlat && isVerticesCData
                    % Take first vertex color (see FaceColor in Patch Properties)
                    fvCData         = fvCData(Faces(:,1)+ 1,:);
                    rowsCData       = size(fvCData,1);
                    isVerticesCData = false;
                end

                % Point meta as true color CData, i.e. RGB in [0,1]
                if size(fvCData,2) == 3
                    % Create additional custom colormap
                    m2t.axes{end}.options(end+1,:) = ...
                        {LibAv05.M2T_matlab2pgfplotsColormap(m2t, fvCData, 'patchmap'), []};
                    drawOptions = M2T_opts_append(drawOptions, 'colormap name','patchmap');

                    % Index into custom colormap
                    fvCData = (0:rowsCData-1)';
                end

                % Add pointmeta data to vertices or faces
                if isVerticesCData
                    columnNames{end+1}   = 'c';
                    verticesTableOptions = M2T_opts_add(verticesTableOptions, 'point meta','\thisrow{c}');
                    Vertices             = [Vertices, fvCData];
                else
                    ptType = 'patch table with point meta';
                    Faces  = [Faces fvCData];
                end

            else
                % Scalar FaceVertexCData, i.e. one color mapping for all patches,
                % used e.g. by Octave in drawing barseries

                [m2t,xFaceColor] = M2T_getColor(m2t, handle, s.faceColor, 'patch');
                drawOptions      = M2T_opts_add(drawOptions, 'fill', xFaceColor);
            end
        end
        % ==============================================================================
        function [drawOptions] = M2T_maybeShowInLegend(showInLegend, drawOptions)
            % sets the appropriate options to show/hide the plot in the legend
            if ~showInLegend
                % No legend entry found. Don't include plot in legend.
                drawOptions = M2T_opts_add(drawOptions, 'forget plot');
            end
        end
        % ==============================================================================
        function [m2t, options] = M2T_setColor(m2t, handle, options, property, color, noneValue)
            % assigns the MATLAB color of the object identified by "handle" to the LaTeX
            % property stored in the options array. An optional "noneValue" can be provided
            % that is set when the color == 'none' (if it is omitted, the property will not
            % be set).
            % TODO: probably this should be integrated with LibAv05.M2T_getAndCheckDefault etc.
            if LibAv05.M2T_opts_has(options,property) && LibAv05.M2T_isNone(LibAv05.M2T_opts_get(options,property))
                return
            end
            if ~LibAv05.M2T_isNone(color)
                [m2t, xcolor] = M2T_getColor(m2t, handle, color, 'patch');
                if ~isempty(xcolor)
                    % this may happen when color == 'flat' and CData is Nx3, e.g. in
                    % scatter plot or in patches
                    if isempty(property)
                        options = M2T_opts_add(options, xcolor);
                    else
                        options = M2T_opts_add(options, property, xcolor);
                    end
                end
            else
                if exist('noneValue','var')
                    options = M2T_opts_add(options, property, noneValue);
                end
            end
        end
        % ==============================================================================
        function drawOptions = M2T_getPatchShape(m2t, h, drawOptions, patchOptions)
            % Retrieves the shape options (i.e. number of vertices) of patch objects
            % Depending on the number of vertices, patches can be triangular, rectangular
            % or polygonal
            % See pgfplots 1.12 manual section 5.8.1 "Additional Patch Types" and the
            % patchplots library
            vertexCount = size(get(h, 'Faces'), 2);

            switch vertexCount
                case 3 % triangle (default)
                    % do nothing special

                case 4 % rectangle
                    drawOptions = M2T_opts_add(drawOptions,'patch type', 'rectangle');

                otherwise % generic polygon
                    LibAv05.M2T_userInfo(m2t, '\nMake sure to load \\usepgfplotslibrary{patchplots} in the preamble.\n');

                    % Default interpolated shader,not supported by polygon, to faceted
                    isFaceColorFlat = isempty(strfind(LibAv05.M2T_opts_get(patchOptions, 'shader'),'interp'));
                    if ~isFaceColorFlat
                        % NOTE: check if pgfplots supports this (or specify version)
                        LibAv05.M2T_userInfo(m2t, '\nPgfplots does not support interpolation for polygons.\n Use patches with at most 4 vertices.\n');
                        patchOptions = M2T_opts_remove(patchOptions, 'shader');
                        patchOptions = M2T_opts_add(patchOptions, 'shader', 'faceted');
                    end

                    % Add draw options
                    drawOptions = M2T_opts_add(drawOptions, 'patch type', 'polygon');
                    drawOptions = M2T_opts_add(drawOptions, 'vertex count', ...
                                                        sprintf('%d', vertexCount));
            end

            drawOptions = M2T_opts_merge(drawOptions, patchOptions);
        end
        % ==============================================================================
        function [cycle] = M2T_conditionallyCyclePath(data)
            % returns "--cycle" when the path should be cyclic in pgfplots
            % Mostly, this is the case UNLESS the data record starts or ends with a NaN
            % record (i.e. a break in the path)
            if any(~isfinite(data([1 end],:)))
                cycle = '';
            else
                cycle = '--cycle';
            end
        end
        % ==============================================================================
        function m2t = M2T_jumpAtUnboundCoords(m2t, data)
            % signals the axis to allow discontinuities in the plot at unbounded
            % coordinates (i.e. Inf and NaN).
            % See also pgfplots 1.12 manual section 4.5.13 "Interrupted Plots".
            if any(~isfinite(data(:)))
                m2t = M2T_needsPgfplotsVersion(m2t, [1 4]);
                m2t = M2T_m2t_addAxisOption(m2t, 'unbounded coords', 'jump');
            end
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawImage(m2t, handle)
            str = '';
            if ~LibAv05.M2T_isVisible(handle)
                return
            end

            % read x-, y-, and color-data
            xData = get(handle, 'XData');
            yData = get(handle, 'YData');
            cData = get(handle, 'CData');

            if (m2t.args.imagesAsPng)
                [m2t, str] = M2T_imageAsPNG(m2t, handle, xData, yData, cData);
            else
                [m2t, str] = M2T_imageAsTikZ(m2t, handle, xData, yData, cData);
            end

            % Make sure that the axes are still visible above the image.
            m2t = M2T_m2t_addAxisOption(m2t, 'axis on top');
        end
        % ==============================================================================
        function [m2t, str] = M2T_imageAsPNG(m2t, handle, xData, yData, cData)
            [m2t, fileNum] = M2T_incrementGlobalCounter(m2t, 'pngFile');
            % ------------------------------------------------------------------------
            % draw a png image
            [pngFileName, pngReferencePath] = M2T_externalFilename(m2t, fileNum, '.png');

            % Get color indices for indexed images and truecolor values otherwise
            if ndims(cData) == 2 %#ok don't use ismatrix (cfr. #143)
                [m2t, colorData] = M2T_cdata2colorindex(m2t, cData, handle);
            else
                colorData = cData;
            end

            m = size(cData, 1);
            n = size(cData, 2);

            alphaData = M2T_normalizedAlphaValues(m2t, get(handle,'AlphaData'), handle);
            if numel(alphaData) == 1
                alphaData = alphaData(ones(size(colorData(:,:,1))));
            end
            [colorData, alphaData] = M2T_flipImageIfAxesReversed(m2t, colorData, alphaData);

            % Write an indexed or a truecolor image
            hasAlpha = true;
            if isfloat(alphaData) && all(alphaData(:) == 1)
                alphaOpts = {};
                hasAlpha = false;
            else
                alphaOpts = {'Alpha', alphaData};
            end
            if (ndims(colorData) == 2) %#ok don't use ismatrix (cfr. #143)
                if size(m2t.current.colormap, 1) <= 256 && ~hasAlpha
                    % imwrite supports maximum 256 values in a colormap (i.e. 8 bit)
                    % and no alpha channel for indexed PNG images.
                    imwrite(colorData, m2t.current.colormap, ...
                        pngFileName, 'png');
                else % use true-color instead
                    imwrite(ind2rgb(colorData, m2t.current.colormap), ...
                        pngFileName, 'png', alphaOpts{:});
                end
            else
                imwrite(colorData, pngFileName, 'png', alphaOpts{:});
            end
            % -----------------------------------------------------------------------
            % dimensions of a pixel in axes units
            if n == 1
                xlim = get(m2t.current.gca, 'xlim');
                xw = xlim(2) - xlim(1);
            else
                xw = (xData(end)-xData(1)) / (n-1);
            end
            if m == 1
                ylim = get(m2t.current.gca, 'ylim');
                yw = ylim(2) - ylim(1);
            else
                yw = (yData(end)-yData(1)) / (m-1);
            end

            opts = M2T_opts_new();
            opts = M2T_opts_add(opts, 'xmin', sprintf(m2t.ff, xData(1  ) - xw/2));
            opts = M2T_opts_add(opts, 'xmax', sprintf(m2t.ff, xData(end) + xw/2));
            opts = M2T_opts_add(opts, 'ymin', sprintf(m2t.ff, yData(1  ) - yw/2));
            opts = M2T_opts_add(opts, 'ymax', sprintf(m2t.ff, yData(end) + yw/2));

            % Print out
            drawOpts = M2T_opts_print(opts);
            str      = sprintf('\\addplot [forget plot] graphics [%s] {%s};\n', ...
                               drawOpts, pngReferencePath);

            LibAv05.M2T_userInfo(m2t, ...
                ['\nA PNG file is stored at ''%s'' for which\n', ...
                'the TikZ file contains a reference to ''%s''.\n', ...
                'You may need to adapt this, depending on the relative\n', ...
                'locations of the master TeX file and the included TikZ file.\n'], ...
                pngFileName, pngReferencePath);
        end
        % ==============================================================================
        function [m2t, str] = M2T_imageAsTikZ(m2t, handle, xData, yData, cData)
            % writes an image as raw TikZ commands (STRONGLY DISCOURAGED)

            % set up cData
            if ndims(cData) == 3
                cData = cData(end:-1:1,:,:);
            else
                cData = cData(end:-1:1,:);
            end

            % Generate uniformly distributed X, Y, although xData and yData may be
            % non-uniform.
            % This is MATLAB(R) behavior.
            [X, hX] = M2T_constructUniformXYDataForImage(xData, size(cData, 2));
            [Y, hY] = M2T_constructUniformXYDataForImage(yData, size(cData, 1));
            [m2t, xcolor] = M2T_getColor(m2t, handle, cData, 'image');

            % The following section takes pretty long to execute, although in
            % principle it is discouraged to use TikZ for those; LaTeX will take
            % forever to compile.
            % Still, a bug has been filed on MathWorks to allow for one-line
            % sprintf'ing with (string+num) cells (Request ID: 1-9WHK4W);
            % <http://www.mathworks.de/support/service_requests/Service_Request_Detail.do?ID=183481&filter=&sort=&statusorder=0&dateorder=0>.
            % An alternative approach could be to use 'surf' or 'patch' of pgfplots
            % with inline tables.
            str = '';
            m = length(X);
            n = length(Y);
            imageString = cell(1, m);
            for i = 1:m
                subString = cell(1, n);
                for j = 1:n
                    subString{j} = sprintf(['\t\\fill [%s] ', ...
                                    '(axis cs:', m2t.ff,',', m2t.ff,') rectangle ', ...
                                    '(axis cs:', m2t.ff,',',m2t.ff,');\n'], ...
                                    xcolor{n-j+1,i}, ...
                                    X(i)-hX/2, Y(j)-hY/2, ...
                                    X(i)+hX/2, Y(j)+hY/2);
                end
                imageString{i} = M2T_join(m2t, subString, '');
            end
            str = M2T_join(m2t, [str, imageString], '');
        end
        function [XY, delta] = M2T_constructUniformXYDataForImage(XYData, expectedLength)
            % Generate uniformly distributed X, Y, although xData/yData may be
            % non-uniform. Dimension indicates the corresponding dimension in the cData matrix.
            switch length(XYData)
                case 2 % only the limits given; common for generic image plots
                    delta = 1;
                case expectedLength % specific x/y-data is given
                    delta = (XYData(end)-XYData(1)) / (length(XYData)-1);
                otherwise
                    error('LibAv05.M2T_drawImage:arrayLengthMismatch', ...
                          'CData length (%d) does not match X/YData length (%d).', ...
                          expectedLength, length(XYData));
            end
            XY = XYData(1):delta:XYData(end);
        end
        % ==============================================================================
        function [colorData, alphaData] = M2T_flipImageIfAxesReversed(m2t, colorData, alphaData)
            % flip the image if reversed
            if m2t.xAxisReversed
                colorData = colorData(:, end:-1:1, :);
                alphaData = alphaData(:, end:-1:1);
            end
            if ~m2t.yAxisReversed % y-axis direction is reversed normally for images, flip otherwise
                colorData = colorData(end:-1:1, :, :);
                alphaData = alphaData(end:-1:1, :);
            end
        end
        % ==============================================================================
        function alpha = M2T_normalizedAlphaValues(m2t, alpha, handle)
            alphaDataMapping = M2T_getOrDefault(handle, 'AlphaDataMapping', 'none');
            switch lower(alphaDataMapping)
                case 'none'  % no rescaling needed
                case 'scaled'
                    ALim = get(m2t.current.gca, 'ALim');
                    AMax = ALim(2);
                    AMin = ALim(1);
                    if ~isfinite(AMax)
                        AMax = max(alpha(:)); %NOTE: is this right?
                    end
                    alpha = (alpha - AMin)./(AMax - AMin);
                case 'direct'
                    alpha = ind2rgb(alpha, get(m2t.current.gcf, 'Alphamap'));
                otherwise
                    error('LibAv05.M2T_matlab2tikz:UnknownAlphaMapping', ...
                          'Unknown alpha mapping "%s"', alphaMapping);
            end

            if isfloat(alpha) %important, alpha data can have integer type which should not be scaled
                alpha = min(1,max(alpha,0)); % clip at range [0, 1]
            end
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawContour(m2t, h)
            if LibAv05.M2T_isHG2()
                [m2t, str] = LibAv05.M2T_LibAv05.M2T_drawContourHG2(m2t, h);
            else
                % Save legend state for the contour group
                hasLegend = m2t.currentHandleHasLegend;

                % Plot children patches
                children  = allchild(h);
                N         = numel(children);
                str       = cell(N,1);
                for ii = 1:N
                    % Plot in reverse order
                    child          = children(N-ii+1);
                    isContourLabel = strcmpi(get(child,'type'),'text');
                    if isContourLabel
                        [m2t, str{ii}] = M2T_drawText(m2t,child);
                    else
                        [m2t, str{ii}] = M2T_drawPatch(m2t,child);
                    end

                    % Only first child can be in the legend
                    m2t.currentHandleHasLegend = false;
                end
                str = strcat(str,sprintf('\n'));
                str = [str{:}];

                % Restore group's legend state
                m2t.currentHandleHasLegend = hasLegend;
            end
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawContourHG2(m2t, h)
            str = '';
            if ~LibAv05.M2T_isVisible(h)
                return
            end

            % Retrieve ContourMatrix
            contours = get(h,'ContourMatrix')';
            [istart, nrows] = M2T_findStartOfContourData(contours);

            % Scale negative contours one level down (for proper coloring)
            Levels    = contours(istart,1);
            LevelList = get(h,'LevelList');
            ineg      = Levels < 0;
            if any(ineg) && min(LevelList) < min(Levels)
                [idx,pos] = ismember(Levels, LevelList);
                idx       = idx & ineg;
                contours(istart(idx)) = LevelList(pos(idx)-1);
            end

            % Draw a contour group (MATLAB R2014b and newer only)
            isFilled = M2T_isOn(get(h,'Fill'));
            if isFilled
                [m2t, str] = M2T_drawFilledContours(m2t, h, contours, istart, nrows);
            else
                % Add colormap
                cmap = m2t.current.colormap;
                m2t = M2T_m2t_addAxisOption(m2t, LibAv05.M2T_matlab2pgfplotsColormap(m2t, cmap));

                % Contour table in Matlab format
                plotOptions = M2T_opts_new();
                plotOptions = M2T_opts_add(plotOptions,'contour prepared');
                plotOptions = M2T_opts_add(plotOptions,'contour prepared format','matlab');

                % Labels
                if LibAv05.M2T_isOff(get(h,'ShowText'))
                    plotOptions = M2T_opts_add(plotOptions,'contour/labels','false');
                end

                % Get line properties
                [m2t, lineOptions] = M2T_getLineOptions(m2t, h);

                % Check for special color settings
                [lineColor, isDefaultColor] = M2T_getAndCheckDefault('contour', h, 'LineColor', 'flat');
                if ~isDefaultColor
                    [m2t, lineOptions] = M2T_setColor(m2t, h, lineOptions, 'contour/draw color', lineColor, 'none');
                end

                % Merge the line options with the contour plot options
                plotOptions = M2T_opts_merge(plotOptions, lineOptions);

                % Make contour table
                [m2t, table, tableOptions] = M2T_makeTable(m2t, {'',''}, contours);

                % Print out
                plotOpts = M2T_opts_print(plotOptions);
                tabOpts  = M2T_opts_print(tableOptions);
                str      = sprintf('\\addplot[%s] table[%s] {%%\n%s};\n', ...
                                   plotOpts, tabOpts, table);
            end
        end
        % ==============================================================================
        function [istart, nrows] = M2T_findStartOfContourData(contours)
            % Index beginning of contour data (see contourc.m for details)
            nrows  = size(contours,1);
            istart = false(nrows,1);
            pos    = 1;
            while pos < nrows
                istart(pos) = true;
                pos         = pos + contours(pos, 2) + 1;
            end
            istart = find(istart);
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawFilledContours(m2t, h, contours, istart, nrows)
            % Loop each contour and plot a filled region
            %
            % NOTE:
            % - we cannot plot from inner to outer contour since the last
            % filled area will cover the inner regions. Therefore, we need to
            % invert the plotting order in those cases.
            % - we need to distinguish between contour groups. A group is
            % defined by inclusion, i.e. its members are contained within one
            % outer contour. The outer contours of two groups cannot include
            % each other.
            str = '';
            if ~LibAv05.M2T_isVisible(h)
                return
            end

            % Split contours in cell array
            cellcont = mat2cell(contours, diff([istart; nrows+1]));
            ncont    = numel(cellcont);

            % Determine contour groups and the plotting order.
            % The ContourMatrix lists the contours in ascending order by level.
            % Hence, if the lowest (first) contour contains any others, then the
            % group will be a peak. Otherwise, the group will be a valley, and
            % the contours will have to be plotted in reverse order, i.e. from
            % highest (largest) to lowest (narrowest).

            %FIXME: close the contours over the border of the domain, see #723.
            order = NaN(ncont,1);
            ifree = true(ncont,1);
            from  = 1;
            while any(ifree)
                % Select peer with lowest level among the free contours, i.e.
                % those which do not belong to any group yet
                pospeer = find(ifree,1,'first');
                peer    = cellcont{pospeer};
                igroup  = false(ncont,1);

                % Loop through all contours
                for ii = 1:numel(cellcont)
                    if ~ifree(ii), continue, end

                    curr = cellcont{ii};
                    % Current contour contained in the peer
                    if inpolygon(curr(2,1),curr(2,2), peer(2:end,1),peer(2:end,2))
                        igroup(ii) = true;
                        isinverse  = false;
                        % Peer contained in the current
                    elseif inpolygon(peer(2,1),peer(2,2),curr(2:end,1),curr(2:end,2))
                        igroup(ii) = true;
                        isinverse  = true;
                    end
                end
                % Order members of group according to the inclusion principle
                nmembers = nnz(igroup ~= 0);
                if isinverse
                    order(igroup) = nmembers+from-1:-1:from;
                else
                    order(igroup) = from:nmembers+from-1;
                end

                % Continue numbering
                from  = from + nmembers;
                ifree = ifree & ~igroup;
            end

            % Reorder the contours
            cellcont(order,1) = cellcont;

            % Add zero level fill
            xdata = get(h,'XData');
            ydata = get(h,'YData');
            %FIXME: determine the contour at the zero level not just its bounding box
            % See also: #721
            zerolevel = [0,          4;
                min(xdata(:)), min(ydata(:));
                min(xdata(:)), max(ydata(:));
                max(xdata(:)), max(ydata(:));
                max(xdata(:)), min(ydata(:))];
            cellcont = [zerolevel; cellcont];

            % Plot
            columnNames = {'x','y'};
            for ii = 1:ncont + 1
                drawOptions = M2T_opts_new();

                % Get fill color
                zval          = cellcont{ii}(1,1);
                [m2t, xcolor] = M2T_getColor(m2t,h,zval,'image');
                drawOptions   = M2T_opts_add(drawOptions,'fill',xcolor);

                % Get line properties
                lineColor = get(h, 'LineColor');

                [m2t, drawOptions] = M2T_setColor(m2t, h, drawOptions, 'draw', lineColor, 'none');

                [m2t, lineOptions] = M2T_getLineOptions(m2t, h);
                drawOptions = M2T_opts_merge(drawOptions, lineOptions);

                % Toggle legend entry
                hasLegend   = ii == 1 && m2t.currentHandleHasLegend;
                drawOptions = M2T_maybeShowInLegend(hasLegend, drawOptions);

                % Print table
                [m2t, table, tableOptions] = M2T_makeTable(m2t, columnNames, cellcont{ii}(2:end,:));

                % Print out
                drawOpts = M2T_opts_print(drawOptions);
                tabOpts  = M2T_opts_print(tableOptions);
                str      = sprintf('%s\\addplot[%s] table[%s] {%%\n%s};\n', ...
                                   str, drawOpts, tabOpts, table);
            end
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawHggroup(m2t, h)
            % Continue according to the plot type. Since the function `handle` is
            % not available in Octave, the plot type will be guessed or the fallback type
            % 'unknown' used.
            % #COMPLEX: big switch-case
            switch getEnvironment()
                case 'MATLAB'
                    cl = class(handle(h));

                case 'Octave'
                    % Function `handle` is not yet implemented in Octave
                    % Consequently the plot type needs to be guessed. See #645.
                    cl = M2T_guessOctavePlotType(h);

                otherwise
                    errorUnknownEnvironment();
            end

            switch(cl)
                case {'specgraph.barseries', 'matlab.graphics.chart.primitive.Bar'}
                    % hist plots and friends
                    [m2t, str] = M2T_drawBarseries(m2t, h);

                case {'specgraph.stemseries', 'matlab.graphics.chart.primitive.Stem'}
                    % stem plots
                    [m2t, str] = M2T_drawStemSeries(m2t, h);

                case {'specgraph.stairseries', 'matlab.graphics.chart.primitive.Stair'}
                    % stair plots
                    [m2t, str] = M2T_drawStairSeries(m2t, h);

                case {'specgraph.areaseries', 'matlab.graphics.chart.primitive.Area'}
                    % scatter plots
                    [m2t,str] = M2T_drawAreaSeries(m2t, h);

                case {'specgraph.quivergroup', 'matlab.graphics.chart.primitive.Quiver'}
                    % quiver arrows
                    [m2t, str] = M2T_drawQuiverGroup(m2t, h);

                case {'specgraph.errorbarseries', 'matlab.graphics.chart.primitive.ErrorBar'}
                    % error bars
                    [m2t,str] = M2T_drawErrorBars(m2t, h);

                case {'specgraph.scattergroup','matlab.graphics.chart.primitive.Scatter'}
                    % scatter plots
                    [m2t,str] = M2T_drawScatterPlot(m2t, h);

                case {'specgraph.contourgroup', 'matlab.graphics.chart.primitive.Contour'}
                    [m2t,str] = M2T_drawContour(m2t, h);

                case {'hggroup', 'matlab.graphics.primitive.Group'}
                    % handle all those the usual way
                    [m2t, str] = M2T_handleAllChildren(m2t, h);

                case 'unknown'
                    % Octave only: plot type could not be determined
                    % Fall back to basic plotting
                    [m2t, str] = M2T_handleAllChildren(m2t, h);

                otherwise
                    LibAv05.M2T_userWarning(m2t, 'Don''t know class ''%s''. Default handling.', cl);
                    try
                        m2tBackup = m2t;
                        [m2t, str] = M2T_handleAllChildren(m2t, h);
                    catch ME
                        LibAv05.M2T_userWarning(m2t, 'Default handling for ''%s'' failed. Continuing as if it did not occur. \n Original Message:\n %s', cl, getReport(ME));
                        [m2t, str] = deal(m2tBackup, ''); % roll-back
                    end
            end
        end
        % ==============================================================================
        % Function `handle` is not yet implemented in Octave.
        % Consequently the plot type needs to be guessed. See #645.
        % If the type can not be determined reliably, 'unknown' will be set.
        function cl = M2T_guessOctavePlotType(h)
            % scatter plots
            if LibAv05.M2T_hasProperties(h, {'marker','sizedata','cdata'}, {})
                cl = 'specgraph.scattergroup';

                % error bars
            elseif LibAv05.M2T_hasProperties(h, {'udata','ldata'}, {})
                cl = 'specgraph.errorbarseries';

                % quiver plots
            elseif LibAv05.M2T_hasProperties(h, {'udata','vdata'}, {'ldata'})
                cl = 'specgraph.quivergroup';

                % bar plots
            elseif LibAv05.M2T_hasProperties(h, {'bargroup','barwidth', 'barlayout'}, {})
                cl = 'specgraph.barseries';
                % unknown plot type
            else
                cl = 'unknown';
            end
        end
        % ==============================================================================
        function bool = M2T_hasProperties(h, fieldsExpectedPresent, fieldsExpectedAbsent)
            % Check if object has all of the given properties (case-insensitive).
            % h                     handle to object (e.g. `gcf` or `gca`)
            % fieldsExpectedPresent cell array of strings with property names to be present
            % fieldsExpectedPresent cell array of strings with property names to be absent
            fields = lower(fieldnames(get(h)));
            present = all(ismember(lower(fieldsExpectedPresent), fields));
            absent = ~any(ismember(lower(fieldsExpectedAbsent), fields));
            bool = present && absent;
        end
        % ==============================================================================
        function m2t = M2T_drawAnnotations(m2t)
            % Draws annotation in Matlab (Octave not supported).

            % In HG1 annotations are children of an invisible axis called scribeOverlay.
            % In HG2 annotations are children of annotationPane object which does not
            % have any axis properties. Hence, we cannot simply handle it with a
            % LibAv05.M2T_drawAxes() call.

            % Octave
            if strcmpi(getEnvironment,'Octave')
                return
            end

            % Get annotation handles
            if LibAv05.M2T_isHG2
                annotPanes   = findall(m2t.current.gcf,'Tag','scribeOverlay');
                children = allchild(annotPanes);
                %TODO: is this dead code?
                if iscell(children)
                    children = [children{:}];
                end
                annotHandles = findall(children,'Visible','on');
            else
                annotHandles = findall(m2t.scribeLayer,'-depth',1,'Visible','on');
            end

            % There are no anotations
            if isempty(annotHandles)
                return
            end

            % Create fake simplified axes overlay (no children)
            warning('off', 'LibAv05.M2T_matlab2tikz:NoChildren')
            scribeLayer = axes('Units','Normalized','Position',[0,0,1,1],'Visible','off');
            m2t         = M2T_drawAxes(m2t, scribeLayer);
            warning('on', 'LibAv05.M2T_matlab2tikz:NoChildren')

            % Plot in reverse to preserve z-ordering and assign the converted
            % annotations to the converted fake overlay
            for ii = numel(annotHandles):-1:1
                m2t = LibAv05.M2T_LibAv05.M2T_drawAnnotationsHelper(m2t,annotHandles(ii));
            end

            % Delete fake overlay graphics object
            delete(scribeLayer)
        end
        % ==============================================================================
        function m2t = M2T_drawAnnotationsHelper(m2t,h)
            % Get class name
            try
                cl = class(handle(h));
            catch
                cl = 'unknown';
            end

            switch cl

                % Line
                case {'scribe.line', 'matlab.graphics.shape.Line'}
                    [m2t, str] = M2T_drawLine(m2t, h);

                    % Ellipse
                case {'scribe.scribeellipse','matlab.graphics.shape.Ellipse'}
                    [m2t, str] = M2T_drawEllipse(m2t, h);

                    % Arrows
                case {'scribe.arrow', 'scribe.doublearrow'}%,...
                      %'matlab.graphics.shape.Arrow', 'matlab.graphics.shape.DoubleEndArrow'}
                    % Annotation: single and double Arrow, line
                    % TODO:
                    % - write a drawArrow(). Handle all info info directly
                    %   without using LibAv05.M2T_handleAllChildren() since HG2 does not have
                    %   children (so no shortcut).
                    % - It would be good if drawArrow() was callable on a
                    %   matlab.graphics.shape.TextArrow object to draw the arrow
                    %   part.
                    [m2t, str] = M2T_handleAllChildren(m2t, h);

                    % Text box
                case {'scribe.textbox','matlab.graphics.shape.TextBox'}
                    [m2t, str] = M2T_drawText(m2t, h);

                    % Tetx arrow
                case {'scribe.textarrow'}%,'matlab.graphics.shape.TextArrow'}
                    % TODO: rewrite LibAv05.M2T_LibAv05.M2T_drawTextarrow. Handle all info info directly
                    %       without using LibAv05.M2T_handleAllChildren() since HG2 does not
                    %       have children (so no shortcut) as used for
                    %       scribe.textarrow.
                    [m2t, str] = LibAv05.M2T_LibAv05.M2T_drawTextarrow(m2t, h);

                    % Rectangle
                case {'scribe.scriberect', 'matlab.graphics.shape.Rectangle'}
                    [m2t, str] = M2T_drawRectangle(m2t, h);

                otherwise
                    LibAv05.M2T_userWarning(m2t, 'Don''t know annotation ''%s''.', cl);
                    return
            end

            % Add annotation to scribe overlay
            m2t.axes{end} = M2T_addChildren(m2t.axes{end}, str);
        end
        % ==============================================================================
        function [m2t,str] = M2T_drawSurface(m2t, h)

            [m2t, opts, s] = M2T_shaderOpts(m2t, h,'surf');
            tableOptions = M2T_opts_new();

            % Allow for empty surf
            if LibAv05.M2T_isNone(s.plotType)
                str = '';
                return
            end

            [dx, dy, dz, numrows] = M2T_getXYZDataFromSurface(h);
            m2t = M2T_jumpAtUnboundCoords(m2t, [dx(:); dy(:); dz(:)]);

            [m2t, opts] = M2T_addZBufferOptions(m2t, h, opts);

            % Check if 3D
            is3D = m2t.axes{end}.is3D;
            if is3D
                columnNames = {'x','y','z','c'};
                plotCmd     = 'addplot3';
                data        = M2T_applyHgTransform(m2t, [dx(:), dy(:), dz(:)]);
            else
                columnNames = {'x','y','c'};
                plotCmd     = 'addplot';
                data        = [dx(:), dy(:)];
            end

            % There are several possibilities of how colors are specified for surface
            % plots:
            %    * explicitly by RGB-values,
            %    * implicitly through a color map with a point-meta coordinate,
            %    * implicitly through a color map with a given coordinate (e.g., z).
            %

            % Check if we need extra CData.
            CData = get(h, 'CData');
            if length(size(CData)) == 3 && size(CData, 3) == 3

                % Create additional custom colormap
                nrows = size(data,1);
                CData = reshape(CData, nrows,3);
                m2t.axes{end}.options(end+1,:) = ...
                    {LibAv05.M2T_matlab2pgfplotsColormap(m2t, CData, 'patchmap'), []};

                % Index into custom colormap
                color = (0:nrows-1)';

                tableOptions = M2T_opts_add(tableOptions, 'colormap name','surfmap');
            else
                opts = M2T_opts_add(opts,LibAv05.M2T_matlab2pgfplotsColormap(m2t, m2t.current.colormap),'');
                % If NaNs are present in the color specifications, don't use them for
                % Pgfplots; they may be interpreted as strings there.
                % Note:
                % Pgfplots actually does a better job than MATLAB in determining what
                % colors to use for the patches. The circular LibAv05.M2T_test case on
                % http://www.mathworks.de/de/help/matlab/ref/pcolor.html, for example
                % yields a symmetric setup in Pgfplots (and doesn't in MATLAB).
                needsPointmeta = any(xor(isnan(dz(:)), isnan(CData(:)))) ...
                    || any(abs(CData(:) - dz(:)) > 1.0e-10);
                if needsPointmeta
                    color = CData(:);
                else
                    color = dz(:);      % Fallback on the z-values, especially if 2D view
                end
            end
            tableOptions = M2T_opts_add(tableOptions, 'point meta','\thisrow{c}');

            data = [data, color];

            % Add mesh/rows=<num rows> for specifying the row data instead of empty
            % lines in the data list below. This makes it possible to reduce the
            % data writing to one single sprintf() call.
            opts = M2T_opts_add(opts,'mesh/rows',sprintf('%d', numrows));

            % Print the addplot options
            str = sprintf('\n\\%s[%%\n%s,\n%s]', plotCmd, s.plotType, LibAv05.M2T_opts_print(opts));

            % Print the data
            [m2t, table, tabOptsExtra] = M2T_makeTable(m2t, columnNames, data);
            tableOptions = M2T_opts_merge(tabOptsExtra, tableOptions);
            tabOpts = M2T_opts_print(tableOptions);

            % Here is where everything is put together
            str = sprintf('%s\ntable[%s] {%%\n%s};\n', ...
                          str, tabOpts, table);

            % TODO:
            % - remove grids in spectrogram by either removing grid command
            %   or adding: 'grid=none' from/in axis options
            % - handling of huge data amounts in LaTeX.

            [m2t, labelString] = M2T_addLabel(m2t, h);
            str = [str, labelString];
        end
        % ==============================================================================
        function [m2t, opts] = M2T_addZBufferOptions(m2t, h, opts)
            % Enforce 'z buffer=sort' if shader is flat and is a 3D plot. It is to
            % avoid overlapping e.g. sphere plots and to properly mimic Matlab's
            % coloring of faces.
            % NOTE:
            % - 'z buffer=sort' is computationally more expensive for LaTeX, we
            %   could try to avoid it in some default situations, e.g. when dx and
            %   dy are rank-1-matrices.
            % - hist3D plots should not be z-sorted or the highest bars will cover
            %   the shorLibAv05.M2T_test one even if positioned in the back
            isShaderFlat = isempty(strfind(LibAv05.M2T_opts_get(opts, 'shader'), 'interp'));
            isHist3D     = strcmpi(get(h,'tag'), 'hist3');
            is3D         = m2t.axes{end}.is3D;
            if is3D && isShaderFlat && ~isHist3D
                opts = M2T_opts_add(opts, 'z buffer', 'sort');
                % Pgfplots 1.12 contains a bug fix that fixes legend entries when
                % 'z buffer=sort' has been set. So, it's  easier to always require that
                % version anyway. See #504 for more information.
                m2t = M2T_needsPgfplotsVersion(m2t, [1,12]);
            end
        end
        % ==============================================================================
        function [dx, dy, dz, numrows] = M2T_getXYZDataFromSurface(h)
            % retrieves X, Y and Z data from a Surface plot. The data gets returned in a
            % wastefull format where the dimensions of these data vectors is equal, akin
            % to the format used by meshgrid.
            dx = get(h, 'XData');
            dy = get(h, 'YData');
            dz = get(h, 'ZData');
            [numcols, numrows] = size(dz);

            % If dx or dy are given as vectors, convert them to the (wasteful) matrix
            % representation first. This makes sure we can treat the data with one
            % single sprintf() command below.
            if isvector(dx)
                dx = ones(numcols,1) * dx(:)';
            end
            if isvector(dy)
                dy = dy(:) * ones(1,numrows);
            end
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawVisibleText(m2t, handle)
            % Wrapper for LibAv05.M2T_drawText() that only draws visible text

            % There may be some text objects floating around a MATLAB figure which are
            % handled by other subfunctions (labels etc.) or don't need to be handled at
            % all.
            % The HandleVisibility says something about whether the text handle is
            % visible as a data structure or not. Typically, a handle is hidden if the
            % graphics aren't supposed to be altered, e.g., axis labels.  Most of those
            % entities are captured by LibAv05.M2T_matlab2tikz one way or another, but sometimes they
            % are not. This is the case, for example, with polar plots and the axis
            % descriptions therein.  Also, Matlab treats text objects with a NaN in the
            % position as invisible.
            if any(isnan(get(handle, 'Position')) | isnan(get(handle, 'Rotation'))) ...
                    || LibAv05.M2T_isOff(get(handle, 'Visible')) ...
                    || (LibAv05.M2T_isOff(get(handle, 'HandleVisibility')) && ...
                        ~m2t.args.showHiddenStrings)

                str = '';
                return;
            end

            [m2t, str] = M2T_drawText(m2t, handle);

        end
        % ==============================================================================
        function [m2t, str] = M2T_drawText(m2t, handle)
            % Adding text node anywhere in the axes environment.
            % Not that, in Pgfplots, long texts get cut off at the axes. This is
            % Different from the default MATLAB behavior. To fix this, one could use
            % /pgfplots/after end axis/.code.

            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % get required properties
            content     = get(handle, 'String');
            Interpreter = get(handle, 'Interpreter');
            content     = M2T_prettyPrint(m2t, content, Interpreter);
            % Concatenate multiple lines
            content = M2T_join(m2t, content, '\\');
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % translate them to pgf style
            style = M2T_opts_new();

            bgColor = get(handle,'BackgroundColor');
            [m2t, style] = M2T_setColor(m2t, handle, style, 'fill', bgColor);

            style = M2T_getXYAlignmentOfText(handle, style);

            style = M2T_getRotationOfText(m2t, handle, style);

            [m2t, fontStyle] = M2T_getFontStyle(m2t, handle);
            style = M2T_opts_merge(style, fontStyle);

            EdgeColor = get(handle, 'EdgeColor');
            [m2t, style] = M2T_setColor(m2t, handle, style, 'draw', EdgeColor);

            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % plot the thing
            [m2t, posString] = M2T_getPositionOfText(m2t, handle);

            styleOpts = M2T_opts_print(style);
            str       = sprintf('\\node[%s]\nat %s {%s};\n', ...
                                styleOpts, posString, content);
        end
        % ==============================================================================
        function [style] = M2T_getXYAlignmentOfText(handle, style)
            % sets the horizontal and vertical alignment options of a text object
            VerticalAlignment   = get(handle, 'VerticalAlignment');
            HorizontalAlignment = get(handle, 'HorizontalAlignment');

            horizontal = '';
            vertical   = '';
            switch VerticalAlignment
                case {'top', 'cap'}
                    vertical = 'below';
                case {'baseline', 'bottom'}
                    vertical = 'above';
            end
            switch HorizontalAlignment
                case 'left'
                    horizontal = 'right';
                case 'right'
                    horizontal = 'left';
            end
            alignment = strtrim(sprintf('%s %s', vertical, horizontal));
            if ~isempty(alignment)
                style = M2T_opts_add(style, alignment);
            end

            % Set 'align' option that is needed for multiline text
            style = M2T_opts_add(style, 'align', HorizontalAlignment);
        end
        % ==============================================================================
        function [style] = M2T_getRotationOfText(m2t, handle, style)
            % Add rotation, if existing
            defaultRotation = 0.0;
            rot = M2T_getOrDefault(handle, 'Rotation', defaultRotation);
            if rot ~= defaultRotation
                style = M2T_opts_add(style, 'rotate', sprintf(m2t.ff, rot));
            end
        end
        % ==============================================================================
        function [m2t,posString] = M2T_getPositionOfText(m2t, h)
            % makes the tikz position string of a text object
            pos   = get(h, 'Position');
            units = get(h, 'Units');
            is3D  = m2t.axes{end}.is3D;

            % Deduce if text or textbox
            type = get(h,'type');
            if isempty(type) || strcmpi(type,'hggroup')
                type = get(h,'ShapeType'); % Undocumented property valid from 2008a
            end

            switch type
                case 'text'
                    if is3D
                        pos  = M2T_applyHgTransform(m2t, pos);
                        npos = 3;
                    else
                        pos  = pos(1:2);
                        npos = 2;
                    end
                case {'textbox','textboxshape'}
                    % TODO:
                    %   - size of the box (e.g. using node attributes minimum width / height)
                    %   - Alignment of the resized box
                    pos  = pos(1:2);
                    npos = 2;

                otherwise
                    error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_drawText', 'Unrecognized text type: %s.', type);
            end

            % Format according to units
            switch units
                case 'normalized'
                    type    = 'rel axis cs:';
                    fmtUnit = '';
                case 'data'
                    type    = 'axis cs:';
                    fmtUnit = '';
                    % Let Matlab do the conversion of any unit into cm
                otherwise
                    type    = '';
                    fmtUnit = 'cm';
                    if ~strcmpi(units, 'centimeters')
                        % Save old pos, set units to cm, query pos, reset
                        % NOTE: cannot use copyobj since it is buggy in R2014a, see
                        %       http://www.mathworks.com/support/bugreports/368385
                        oldPos = get(h, 'Position');
                        set(h,'Units','centimeters')
                        pos    = get(h, 'Position');
                        pos    = pos(1:npos);
                        set(h,'Units',units,'Position',oldPos)
                    end
            end
            posString = cell(1,npos);
            for ii = 1:npos
                posString{ii} = M2T_formatDim(pos(ii), fmtUnit);
            end

            posString = sprintf('(%s%s)',type,LibAv05.M2T_join(m2t,posString,','));
            m2t = M2T_disableClippingInCurrentAxes(m2t, pos);

        end
        % ==============================================================================
        function m2t = M2T_disableClippingInCurrentAxes(m2t, pos)
            % Disables clipping in the current axes if the `pos` vector lies outside
            % the limits of the axes.
            xlim  = M2T_getOrDefault(m2t.current.gca, 'xlim',[-Inf +Inf]);
            ylim  = M2T_getOrDefault(m2t.current.gca, 'ylim',[-Inf +Inf]);
            zlim  = M2T_getOrDefault(m2t.current.gca, 'ZLim',[-Inf +Inf]);
            is3D  = m2t.axes{end}.is3D;

            xOutOfRange =          pos(1) < xlim(1) || pos(1) > xlim(2);
            yOutOfRange =          pos(2) < ylim(1) || pos(2) > ylim(2);
            zOutOfRange = is3D && (pos(3) < zlim(1) || pos(3) > zlim(2));
            if xOutOfRange || yOutOfRange || zOutOfRange
                m2t = M2T_m2t_addAxisOption(m2t, 'clip', 'false');
            end
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawRectangle(m2t, h)
            str = '';

            % there may be some text objects floating around a Matlab figure which
            % are handled by other subfunctions (labels etc.) or don't need to be
            % handled at all
            if ~LibAv05.M2T_isVisible(h) || LibAv05.M2T_isOff(get(h, 'HandleVisibility'))
                return;
            end

            % TODO handle Curvature = [0.8 0.4]

            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Get draw options.
            [m2t, lineOptions] = M2T_getLineOptions(m2t, h);
            [m2t, lineOptions] = M2T_getRectangleFaceOptions(m2t, h, lineOptions);
            [m2t, lineOptions] = M2T_getRectangleEdgeOptions(m2t, h, lineOptions);
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            pos = M2T_pos2dims(get(h, 'Position'));
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % plot the thing
            lineOpts = M2T_opts_print(lineOptions);
            str = sprintf(['\\draw[%s] (axis cs:',m2t.ff,',',m2t.ff, ')', ...
                           ' rectangle (axis cs:',m2t.ff,',',m2t.ff,');\n'], ...
                           lineOpts, pos.left, pos.bottom, pos.right, pos.top);
        end
        % ==============================================================================
        function [m2t, drawOptions] = M2T_getRectangleFaceOptions(m2t, h, drawOptions)
            % draws the face (i.e. fill) of a Rectangle
            faceColor    = get(h, 'FaceColor');
            isAnnotation = strcmpi(get(h,'type'),'rectangleshape') || ...
                           strcmpi(LibAv05.M2T_getOrDefault(h,'ShapeType',''),'rectangle');
            isFlatColor  = strcmpi(faceColor, 'flat');
            if ~(LibAv05.M2T_isNone(faceColor) || (isAnnotation && isFlatColor))
                [m2t, xFaceColor] = M2T_getColor(m2t, h, faceColor, 'patch');
                drawOptions = M2T_opts_add(drawOptions, 'fill', xFaceColor);
            end
        end
        % ==============================================================================
        function [m2t, drawOptions] = M2T_getRectangleEdgeOptions(m2t, h, drawOptions)
            % draws the edges of a rectangle
            edgeColor = get(h, 'EdgeColor');
            lineStyle = get(h, 'LineStyle');
            if LibAv05.M2T_isNone(lineStyle) || LibAv05.M2T_isNone(edgeColor)
                drawOptions = M2T_opts_add(drawOptions, 'draw', 'none');
            else
                [m2t, drawOptions] = M2T_setColor(m2t, h, drawOptions, 'draw', edgeColor);
            end
        end
        % ==============================================================================
        function [m2t,opts,s] = M2T_shaderOpts(m2t, handle, selectedType)
            % SHADEROPTS Returns the shader, fill and draw options for patches, surfs and meshes
            %
            %   SHADEROPTS(M2T, HANDLE, SELECTEDTYPE) Where SELECTEDTYPE should either
            %   be 'surf' or 'patch'
            %
            %
            %   [...,OPTS, S] = SHADEROPTS(...)
            %       OPTS is a M by 2 cell array with Key/Value pairs
            %       S is a struct with fields, e.g. 'faceColor', to be re-used by the
            %       caller

            % Initialize
            opts              = M2T_opts_new;
            s.hasOneEdgeColor = false;
            s.hasOneFaceColor = false;

            % Get relevant Face and Edge color properties
            s.faceColor = get(handle, 'FaceColor');
            s.edgeColor = get(handle, 'EdgeColor');

            if LibAv05.M2T_isNone(s.faceColor) && LibAv05.M2T_isNone(s.edgeColor)
                s.plotType        = 'none';
                s.hasOneEdgeColor = true;
            elseif LibAv05.M2T_isNone(s.faceColor)
                s.plotType        = 'mesh';
                s.hasOneFaceColor = true;
                [m2t, opts, s]    = LibAv05.M2T_LibAv05.M2T_shaderOptsMesh(m2t, handle, opts, s);
            else
                s.plotType     = selectedType;
                [m2t, opts, s] = LibAv05.M2T_LibAv05.M2T_shaderOptsSurfPatch(m2t, handle, opts, s);
            end
        end
        % ==============================================================================
        function [m2t, opts, s] = M2T_shaderOptsMesh(m2t, handle, opts, s)

            % Edge 'interp'
            if strcmpi(s.edgeColor, 'interp')
                opts = M2T_opts_add(opts,'shader','flat');

                % Edge RGB
            else
                s.hasOneEdgeColor = true;
                [m2t, xEdgeColor] = M2T_getColor(m2t, handle, s.edgeColor, 'patch');
                opts              = M2T_opts_add(opts,'color',xEdgeColor);
            end
        end
        % ==============================================================================
        function [m2t, opts, s] = M2T_shaderOptsSurfPatch(m2t, handle, opts, s)
            % gets the shader options for surface patches

            % Set opacity if FaceAlpha < 1 in MATLAB
            s.faceAlpha = get(handle, 'FaceAlpha');
            if isnumeric(s.faceAlpha) && s.faceAlpha ~= 1.0
                opts = M2T_opts_add(opts,'fill opacity',sprintf(m2t.ff,s.faceAlpha));
            end

            % Set opacity if EdgeAlpha < 1 in MATLAB
            s.edgeAlpha = get(handle, 'EdgeAlpha');
            if isnumeric(s.edgeAlpha) && s.edgeAlpha ~= 1.0
                opts = M2T_opts_add(opts,'draw opacity',sprintf(m2t.ff,s.edgeAlpha));
            end

            if LibAv05.M2T_isNone(s.edgeColor) % Edge 'none'
                [m2t, opts, s] = LibAv05.M2T_LibAv05.M2T_LibAv05.M2T_shaderOptsSurfPatchEdgeNone(m2t, handle, opts, s);

            elseif strcmpi(s.edgeColor, 'interp') % Edge 'interp'
                [m2t, opts, s] = LibAv05.M2T_LibAv05.M2T_LibAv05.M2T_shaderOptsSurfPatchEdgeInterp(m2t, handle, opts, s);

            elseif strcmpi(s.edgeColor, 'flat') % Edge 'flat'
                [m2t, opts, s] = LibAv05.M2T_LibAv05.M2T_LibAv05.M2T_shaderOptsSurfPatchEdgeFlat(m2t, handle, opts, s);

            else % Edge RGB
                [m2t, opts, s] = LibAv05.M2T_LibAv05.M2T_LibAv05.M2T_shaderOptsSurfPatchEdgeRGB(m2t, handle, opts, s);
            end
        end
        % ==============================================================================
        function [m2t, opts, s] = M2T_shaderOptsSurfPatchEdgeNone(m2t, handle, opts, s)
            % gets the shader options for surface patches without edges
            s.hasOneEdgeColor = true; % consider void as true
            if strcmpi(s.faceColor, 'flat')
                opts = M2T_opts_add(opts,'shader','flat');
            elseif strcmpi(s.faceColor, 'interp');
                opts = M2T_opts_add(opts,'shader','interp');
            else
                s.hasOneFaceColor = true;
                [m2t,xFaceColor]  = M2T_getColor(m2t, handle, s.faceColor, 'patch');
                opts              = M2T_opts_add(opts,'fill',xFaceColor);
            end
        end
        function [m2t, opts, s] = M2T_shaderOptsSurfPatchEdgeInterp(m2t, handle, opts, s)
            % gets the shader options for surface patches with interpolated edge colors
            if strcmpi(s.faceColor, 'interp')
                opts = M2T_opts_add(opts,'shader','interp');
            elseif strcmpi(s.faceColor, 'flat')
                opts = M2T_opts_add(opts,'shader','faceted');
            else
                s.hasOneFaceColor = true;
                [m2t,xFaceColor]  = M2T_getColor(m2t, handle, s.faceColor, 'patch');
                opts              = M2T_opts_add(opts,'fill',xFaceColor);
            end
        end
        function [m2t, opts, s] = M2T_shaderOptsSurfPatchEdgeFlat(m2t, handle, opts, s)
            % gets the shader options for surface patches with flat edge colors, i.e. the
            % vertex color
            if strcmpi(s.faceColor, 'flat')
                opts = M2T_opts_add(opts,'shader','flat corner');
            elseif strcmpi(s.faceColor, 'interp')
                LibAv05.M2T_warnFacetedInterp(m2t);
                opts = M2T_opts_add(opts,'shader','faceted interp');
            else
                s.hasOneFaceColor = true;
                opts              = M2T_opts_add(opts,'shader','flat corner');
                [m2t,xFaceColor]  = M2T_getColor(m2t, handle, s.faceColor, 'patch');
                opts              = M2T_opts_add(opts,'fill',xFaceColor);
            end
        end
        function [m2t, opts, s] = M2T_shaderOptsSurfPatchEdgeRGB(m2t, handle, opts, s)
            % gets the shader options for surface patches with fixed (RGB) edge color
            s.hasOneEdgeColor = true;
            [m2t, xEdgeColor] = M2T_getColor(m2t, handle, s.edgeColor, 'patch');
            if isnumeric(s.faceColor)
                s.hasOneFaceColor = true;
                [m2t, xFaceColor] = M2T_getColor(m2t, handle, s.faceColor, 'patch');
                opts              = M2T_opts_add(opts,'fill',xFaceColor);
                opts              = M2T_opts_add(opts,'faceted color',xEdgeColor);
            elseif strcmpi(s.faceColor,'interp')
                LibAv05.M2T_warnFacetedInterp(m2t);
                opts = M2T_opts_add(opts,'shader','faceted interp');
                opts = M2T_opts_add(opts,'faceted color',xEdgeColor);
            else
                opts = M2T_opts_add(opts,'shader','flat corner');
                opts = M2T_opts_add(opts,'draw',xEdgeColor);
            end
        end
        % ==============================================================================
        function M2T_warnFacetedInterp(m2t)
            % warn the user about the space implications of "shader=faceted interp"
            LibAv05.M2T_userWarning(m2t, ...
                ['A 3D plot with "shader = faceted interp" is being produced.\n', ...
                'This may produce big and sluggish PDF files.\n', ...
                'See %s and Section 4.6.6 of the pgfplots manual for workarounds.'], ...
                LibAv05.M2T_issueUrl(m2t, 693, true));
        end
        % ==============================================================================
        function url = M2T_issueUrl(m2t, number, forOutput)
            % Produces the URL for an issue report in the GitHub repository.
            % When the `forOutput` flag is set, this format the URL for printing to the
            % MATLAB terminal.
            if ~exist('forOutput','var') || isempty(forOutput)
                forOutput = false;
            end
            url = sprintf('%s/%d', m2t.about.issues, number);
            if forOutput
                url = M2T_clickableUrl(url, sprintf('#%d', number));
            end
        end
        % ==============================================================================
        function url = M2T_clickableUrl(url, title)
            % Produce a clickable URL for outputting to the MATLAB terminal
            if ~exist('title','var') || isempty(title)
                title = url;
            end
            switch getEnvironment()
                case 'MATLAB'
                    url = sprintf('<a href="%s">%s</a>', url, title);
                case 'Octave'
                    % just use the URL and discard the title since Octave doesn't
                    % support HTML tags in its output.
                otherwise
                    errorUnknownEnvironment();
            end
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawScatterPlot(m2t, h)
            % DRAWSCATTERPLOT Draws a scatter plot
            %
            % A scatter plot is a plot containing only markers and where the
            % size and/or color of each marker can be changed independently.
            %
            % References for TikZ code:
            %  - http://tex.stackexchange.com/questions/197270/how-to-plot-scatter-points-using-pgfplots-with-color-defined-from-table-rgb-valu
            %  - http://tex.stackexchange.com/questions/98646/multiple-different-meta-for-marker-color-and-marker-size
            %
            % See also: scatter
            str = '';
            if ~LibAv05.M2T_isVisible(h)
                return; % there is nothing to plot
            end

            dataInfo   = M2T_getDataInfo(h, 'X','Y','Z','C','Size');
            markerInfo = M2T_getMarkerInfo(m2t, h);

            if isempty(dataInfo.C) && strcmpi(getEnvironment(), 'Octave')
                dataInfo.C = get(h, 'MarkerEdgeColor');
            end

            %TODO: check against LibAv05.M2T_getMarkerOptions() for duplicated code

            dataInfo.Size = M2T_tryToMakeScalar(dataInfo.Size, m2t.tol);

            % Rescale marker size (not definitive, follow discussion in #316)
            % Prescale marker size for octave
            if strcmpi(getEnvironment(), 'Octave')
                dataInfo.Size = dataInfo.Size.^2/2;
            end
            dataInfo.Size = LibAv05.M2T_LibAv05.M2T_translateMarkerSize(m2t, markerInfo.style, sqrt(dataInfo.Size)/2);

            drawOptions = M2T_opts_new();

            %% Determine if we are drawing an actual scatter plot
            hasDifferentSizes  = numel(dataInfo.Size) ~= 1;
            hasDifferentColors = numel(dataInfo.C)    ~= 3;
            isaScatter         = hasDifferentSizes || hasDifferentColors;
            if isaScatter
                drawOptions = M2T_opts_add(drawOptions, 'scatter');
            end
            %TODO: we need to set the scatter source
            drawOptions = M2T_opts_add(drawOptions, 'only marks');
            drawOptions = M2T_opts_add(drawOptions, 'mark', markerInfo.tikz);

            if length(dataInfo.C) == 3
                % gets options specific to scatter plots with a single color
                % No special treatment for the colors or markers are needed.
                % All markers have the same color.
                [m2t, xcolor, hasFaceColor] = LibAv05.M2T_LibAv05.M2T_getColorOfMarkers(m2t, h, 'MarkerFaceColor', dataInfo.C);
                [m2t, ecolor, hasEdgeColor] = LibAv05.M2T_LibAv05.M2T_getColorOfMarkers(m2t, h, 'MarkerEdgeColor', dataInfo.C);

                if length(dataInfo.Size) == 1;
                    drawOptions = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(drawOptions, 'mark options', ...
                                               markerInfo.options);
                    drawOptions = M2T_opts_add(drawOptions, 'mark size', ...
                        sprintf('%.4fpt', dataInfo.Size)); % FIXME: investigate whether to use `m2t.ff`
                    if hasEdgeColor
                        drawOptions = M2T_opts_add(drawOptions, 'draw', ecolor);
                    else
                        drawOptions = M2T_opts_add(drawOptions, 'color', xcolor); %TODO: why do we even need this one?
                    end
                    if hasFaceColor
                        drawOptions = M2T_opts_add(drawOptions, 'fill', xcolor);
                    end
                else % if changing marker size but same color on all marks
                    markerOptions = M2T_opts_new();
                    markerOptions = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(markerOptions, 'mark options', ...
                                                 markerInfo.options);
                    if hasEdgeColor
                        markerOptions = M2T_opts_add(markerOptions, 'draw', ecolor);
                    else
                        markerOptions = M2T_opts_add(markerOptions, 'draw', xcolor);
                    end
                    if hasFaceColor
                        markerOptions = M2T_opts_add(markerOptions, 'fill', xcolor);
                    end
                    % for changing marker size, the 'scatter' option has to be added
                    drawOptions = M2T_opts_add(drawOptions, 'color', xcolor);
                    drawOptions = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(drawOptions, 'mark options', ...
                                               markerInfo.options);

                    if ~hasFaceColor
                        drawOptions = M2T_opts_add(drawOptions, ...
                            'scatter/use mapped color', xcolor);
                    else
                        drawOptions = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(drawOptions, ...
                            'scatter/use mapped color', markerOptions);
                    end
                end
            elseif size(dataInfo.C,2) == 3
                % scatter plots with each marker a different RGB color (not yet supported!)
                LibAv05.M2T_userWarning(m2t, 'Pgfplots cannot handle RGB scatter plots yet.');
                % TODO Get this in order as soon as Pgfplots can do "scatter rgb".
                % See e.g. http://tex.stackexchange.com/questions/197270 and #433
            else
                % scatter plot where the colors are set using a color map
                markerOptions = M2T_opts_new();
                markerOptions = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(markerOptions, 'mark options', ...
                                             markerInfo.options);
                if markerInfo.hasEdgeColor && markerInfo.hasFaceColor
                    [m2t, ecolor] = M2T_getColor(m2t, h, markerInfo.EdgeColor, 'patch');
                    markerOptions = M2T_opts_add(markerOptions, 'draw', ecolor);
                else
                    markerOptions = M2T_opts_add(markerOptions, 'draw', 'mapped color');
                end
                if markerInfo.hasFaceColor
                    markerOptions = M2T_opts_add(markerOptions, 'fill', 'mapped color');
                end

                if numel(dataInfo.Size) == 1
                    drawOptions = M2T_opts_add(drawOptions, 'mark size', ...
                        sprintf('%.4fpt', dataInfo.Size)); % FIXME: investigate whether to use `m2t.ff` 
                else
                    %TODO: warn the user about this. It is not currently supported.
                end

                drawOptions = M2T_opts_add(drawOptions, 'scatter src', 'explicit');
                drawOptions = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(drawOptions, 'scatter/use mapped color', ...
                                           markerOptions);
                % Add color map.
                m2t = M2T_m2t_addAxisOption(m2t, LibAv05.M2T_matlab2pgfplotsColormap(m2t, m2t.current.colormap), []);
            end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Plot the thing.
            [env, data, metaPart, columns] = M2T_organizeScatterData(m2t, dataInfo);

            if hasDifferentSizes
                drawOptions = M2T_opts_append(drawOptions, 'visualization depends on', ...
                    '{\thisrow{size} \as \perpointmarksize}');
                drawOptions = M2T_opts_add(drawOptions, ...
                    'scatter/@pre marker code/.append style', ...
                    '{/tikz/mark size=\perpointmarksize}');
            end

            % The actual printing.
            [m2t, table, tableOptions] = M2T_makeTable(m2t, columns, data);
            tableOptions = M2T_opts_merge(tableOptions, metaPart);

            % Print
            drawOpts = M2T_opts_print(drawOptions);
            tabOpts  = M2T_opts_print(tableOptions);
            str      = sprintf('\\%s[%s] table[%s]{%s};\n',...
                               env, drawOpts, tabOpts, table);
        end
        % ==============================================================================
        function dataInfo = M2T_getDataInfo(h, varargin)
            % retrieves the "*Data  fields from a HG object
            % When no names are specified, it assumes 'X','Y','Z' is requested
            if nargin == 1
                fields = {'X','Y','Z'};
            else
                fields = varargin;
            end
            dataInfo = struct();
            for iField = 1:numel(fields)
                name            = fields{iField};
                dataInfo.(name) = get(h, [name 'Data']);
            end
        end
        % ==============================================================================
        function value = M2T_tryToMakeScalar(value, tolerance)
            % make a vector into a scalar when all its components are equal
            if ~exist('tolerance','var')
                tolerance = 0; % do everything perfectly
            end
            if all(abs(value - value(1)) <= tolerance)
                value = value(1);
            end
        end
        % ==============================================================================
        function marker = M2T_getMarkerInfo(m2t, h, markOptions)
            % gets marker-related options as a struct
            if ~exist('markOptions','var') || isempty(markOptions)
                markOptions = M2T_opts_new();
            end
            marker                        = struct();
            marker.style                  = get(h, 'Marker');
            marker.FaceColor              = get(h, 'MarkerFaceColor');
            marker.EdgeColor              = get(h, 'MarkerEdgeColor');
            marker.hasFaceColor           = ~LibAv05.M2T_isNone(marker.FaceColor);
            marker.hasEdgeColor           = ~LibAv05.M2T_isNone(marker.EdgeColor);
            [marker.tikz, marker.options] = M2T_translateMarker(m2t, marker.style, ...
                                                    markOptions, marker.hasFaceColor);
        end
        % ==============================================================================
        function [env, data, metaOptions, columns] = M2T_organizeScatterData(m2t, dataInfo)
            % reorganizes the {X,Y,Z,S} data into a single matrix
            metaOptions = M2T_opts_new();


            xData = dataInfo.X;
            yData = dataInfo.Y;
            zData = dataInfo.Z;
            cData = dataInfo.C;
            sData = dataInfo.Size;

            % add the actual data
            if ~m2t.axes{end}.is3D
                env     = 'addplot';
                columns = {'x','y'};
                data    = [xData(:), yData(:)];
            else
                env     = 'addplot3';
                columns = {'x','y','z'};
                data    = M2T_applyHgTransform(m2t, [xData(:), yData(:), zData(:)]);
            end

            % add marker sizes
            if length(sData) ~= 1
                columns = [columns, {'size'}];
                data    = [data, sData(:)];
            end

            % add color data
            if length(cData) == 3
                % If size(cData,1)==1, then all the colors are the same and have
                % already been accounted for above.

            elseif size(cData,2) == 3
                %TODO Hm, can't deal with this?
                %[m2t, col] = M2T_rgb2colorliteral(m2t, cData(k,:));
                %str = strcat(str, sprintf(' [%s]\n', col));
                columns = [columns, {'R','G','B'}];
                data    = [data, cData(:,1), cData(:,2), cData(:,3)];
            else
                columns = [columns, {'color'}];
                metaOptions = M2T_opts_add(metaOptions, 'meta', 'color');
                data = [data, cData(:)];
            end
        end
        % ==============================================================================
        function [m2t, xcolor, hasColor] = M2T_getColorOfMarkers(m2t, h, name, cData)
            color = get(h, name);
            hasColor = ~LibAv05.M2T_isNone(color);
            if hasColor && ~strcmpi(color,'flat');
                [m2t, xcolor] = M2T_getColor(m2t, h, color, 'patch');
            else
                [m2t, xcolor] = M2T_getColor(m2t, h, cData, 'patch');
            end
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawHistogram(m2t, h)
            % Takes care of plots like the ones produced by MATLAB's histogram function.
            % The main pillar is Pgfplots's '{x,y}bar' plot.
            %
            % TODO Get rid of code duplication with 'LibAv05.M2T_drawAxes'.

            % Do nothing if plot is invisible
            str = '';
            if ~LibAv05.M2T_isVisible(h)
                return;
            end

            % Init drawOptions
            drawOptions = M2T_opts_new();

            % Data
            binEdges = get(h, 'BinEdges');
            binValue = get(h, 'Values');
            data     = [binEdges(:), [binValue(:); binValue(end)]];

            % Check for orientation of the bars
            isHorizontal = ~strcmpi(get(h, 'Orientation'), 'vertical');
            if isHorizontal
                drawOptions = M2T_opts_add(drawOptions, 'xbar interval');
                data        = fliplr(data);
            else
                drawOptions = M2T_opts_add(drawOptions, 'ybar interval');
            end

            % Get the draw options for the bars
            [m2t, drawOptions] = M2T_getPatchDrawOptions(m2t, h, drawOptions);

            % Make table
            [m2t, table, tableOptions] = M2T_makeTable(m2t, {'x','y'},data);

            % Print out
            drawOpts = M2T_opts_print(drawOptions);
            tabOpts  = M2T_opts_print(tableOptions);
            str      = sprintf('\\addplot[%s] table[%s] {%s};\n', ...
                               drawOpts, tabOpts, table);
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawBarseries(m2t, h)
            % Takes care of plots like the ones produced by MATLAB's bar function.
            % The main pillar is Pgfplots's '{x,y}bar' plot.
            %
            % TODO Get rid of code duplication with 'LibAv05.M2T_drawAxes'.

            % Do nothing if plot is invisible
            str = '';
            if ~LibAv05.M2T_isVisible(h)
                return;
            end

            % Init drawOptions
            drawOptions = M2T_opts_new();

            % Check for orientation of the bars and their layout
            isHorizontal = M2T_isOn(get(h, 'Horizontal'));
            if isHorizontal
                barType = 'xbar';
            else
                barType = 'ybar';
            end

            % Get the draw options for the layout
            [m2t, drawOptions] = M2T_setBarLayoutOfBarSeries(m2t, h, barType, drawOptions);

            % Get the draw options for the bars
            [m2t, drawOptions] = M2T_getPatchDrawOptions(m2t, h, drawOptions);

            % Add 'log origin = infty' if BaseValue differs from zero (log origin=0 is
            % the default behaviour since Pgfplots v1.5).
            baseValue = get(h, 'BaseValue');
            if baseValue ~= 0.0
                m2t = M2T_m2t_addAxisOption(m2t, 'log origin', 'infty');
                %TODO: wait for pgfplots to implement other base values (see #438)
            end

            % Generate the tikz table
            xData = get(h, 'XData');
            yData = get(h, 'YData');
            if isHorizontal
                [yDataPlot, xDataPlot] = deal(xData, yData); % swap values
            else
                [xDataPlot, yDataPlot] = deal(xData, yData);
            end
            [m2t, table, tableOptions] = M2T_makeTable(m2t, '', xDataPlot, '', yDataPlot);

            % Print out
            drawOpts = M2T_opts_print(drawOptions);
            tabOpts  = M2T_opts_print(tableOptions);
            str      = sprintf('\\addplot[%s] table[%s] {%s};\n', ...
                               drawOpts, tabOpts, table);
            % Add a baseline if appropriate
            [m2t, baseline] = M2T_drawBaseline(m2t,h,isHorizontal);
            str             = [str, baseline];
        end
        % ==============================================================================
        function BarWidth = M2T_getBarWidthInAbsolutUnits(h)
            % determines the width of a bar in a bar plot
            XData = get(h,'XData');
            BarWidth = get(h, 'BarWidth');
            if length(XData) > 1
                BarWidth = min(diff(XData)) * BarWidth;
            end
        end
        % ==============================================================================
        function [m2t, drawOptions] = M2T_setBarLayoutOfBarSeries(m2t, h, barType, drawOptions)
            % sets the options specific to a bar layour (grouped vs stacked)
            barlayout = get(h, 'BarLayout');

            switch barlayout
                case 'grouped'  % grouped bar plots

                    % Get number of bars series and bar series id
                    [numBarSeries, barSeriesId] = M2T_getNumBarAndId(h);

                    % Maximum group width relative to the minimum distance between two
                    % x-values. See <MATLAB>/toolbox/matlab/specgraph/makebars.m
                    maxGroupWidth = 0.8;
                    if numBarSeries == 1
                        groupWidth = 1.0;
                    else
                        groupWidth = min(maxGroupWidth, numBarSeries/(numBarSeries+1.5));
                    end

                    % Calculate the width of each bar and the center point shift as in
                    % makebars.m
                    % Get the shifts of the bar centers.
                    % In case of numBars==1, this returns 0,
                    % In case of numBars==2, this returns [-1/4, 1/4],
                    % In case of numBars==3, this returns [-1/3, 0, 1/3],
                    % and so forth.
                    % assumedBarWidth = groupWidth/numBarSeries; % assumption
                    % barShift        = (barSeriesId - 0.5) * assumedBarWidth - groupWidth/2;
                    % FIXME #785: The previous version of barshift lead to
                    % regressions, as the bars were stacked.
                    % Instead remove the calculation of barShift and add x/ybar to
                    % the axis so that pgf determines it automatically.

                    % From http://www.mathworks.com/help/techdoc/ref/bar.html:
                    % bar(...,width) sets the relative bar width and controls the
                    % separation of bars within a group. The default width is 0.8, so if
                    % you do not specify X, the bars within a group have a slight
                    % separation. If width is 1, the bars within a group touch one
                    % another. The value of width must be a scalar.
                    assumedBarWidth = groupWidth/numBarSeries; % assumption
                    barWidth = M2T_getBarWidthInAbsolutUnits(h) * assumedBarWidth;

                    % Bar type
                    drawOptions = M2T_opts_add(drawOptions, barType);

                    % Bar width
                    drawOptions = M2T_opts_add(drawOptions, 'bar width', LibAv05.M2T_formatDim(barWidth, ''));

                    % The bar shift auto feature was introduced in pgfplots 1.13
                    m2t = M2T_needsPgfplotsVersion(m2t, [1,13]);
                    m2t = M2T_m2t_addAxisOption(m2t, 'bar shift auto');
                case 'stacked' % stacked plots
                    % Pass option to parent axis & disallow anything but stacked plots
                    % Make sure this happens exactly *once*.

                    if ~m2t.axes{end}.barAddedAxisOption;
                        barWidth = M2T_getBarWidthInAbsolutUnits(h);
                        m2t = M2T_m2t_addAxisOption(m2t, 'bar width', LibAv05.M2T_formatDim(barWidth,''));
                        m2t.axes{end}.barAddedAxisOption = true;
                    end

                    % Somewhere between pgfplots 1.5 and 1.8 and starting
                    % again from 1.11, the option {x|y}bar stacked can be applied to
                    % \addplot instead of the figure and thus allows to combine stacked
                    % bar plots and other kinds of plots in the same axis.
                    % Thus, it is advisable to use pgfplots 1.11. In older versions, the
                    % plot will only contain a single bar series, but should compile fine.
                    m2t = M2T_needsPgfplotsVersion(m2t, [1,11]);
                    drawOptions = M2T_opts_add(drawOptions, [barType ' stacked']);
                otherwise
                    error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_drawBarseries', ...
                        'Don''t know how to handle BarLayout ''%s''.', barlayout);
            end
        end
        % ==============================================================================
        function [numBarSeries, barSeriesId] = M2T_getNumBarAndId(h)
            % Get number of bars series and bar series id
            prop         = M2T_switchMatOct('BarPeers', 'bargroup');
            bargroup     = get(h, prop);
            numBarSeries = numel(bargroup);

            if LibAv05.M2T_isHG2
                % In HG2, BarPeers are sorted in reverse order wrt HG1
                bargroup = bargroup(end:-1:1);

            elseif strcmpi(getEnvironment, 'MATLAB')
                % In HG1, h is a double but bargroup a graphic object. Cast h to a
                % graphic object
                h = handle(h);

            else
                % In Octave, the bargroup is a replicated cell array. Pick first
                if iscell(bargroup)
                    bargroup = bargroup{1};
                end
            end

            % Get bar series Id
            [dummy, barSeriesId] = ismember(h, bargroup);
        end
        % ==============================================================================
        function [m2t,str] = M2T_drawBaseline(m2t,hparent,isVertical)
            % DRAWBASELINE Draws baseline for bar and stem plots
            %
            % Notes:
            %   - In HG2, the baseline is a specific object child of a bar or stem
            %     plot. So, LibAv05.M2T_handleAllChildren() won't find a line in the axes to plot as
            %     the baseline.
            %   - The baseline is horizontal for vertical bar and stem plots and is
            %     vertical for horixontal barplots. The ISVERTICAL input refers to the
            %     baseline.
            %   - We do not plot baselines with a BaseValue different from 0 because
            %     pgfplots does not support shifts in the BaseValue, e.g. see #438.
            %     We either implement our own data shifting or wait for pgfplots.

            if ~exist('isVertical','var')
                isVertical = false;
            end

            str = '';
            baseValue = get(hparent, 'BaseValue');
            if LibAv05.M2T_isOff(get(hparent,'ShowBaseLine')) || ~LibAv05.M2T_isHG2() || baseValue ~= 0
                return
            end

            hBaseLine = get(hparent,'BaseLine');

            % Line options of the baseline
            [m2t, lineOptions] = M2T_getLineOptions(m2t, hparent);
            color              = get(hBaseLine, 'Color');
            [m2t, lineColor]   = M2T_getColor(m2t, hBaseLine, color, 'patch');

            drawOptions = M2T_opts_new();
            drawOptions = M2T_opts_add(drawOptions, 'forget plot');
            drawOptions = M2T_opts_add(drawOptions, 'color', lineColor);
            drawOptions = M2T_opts_merge(drawOptions, lineOptions);

            % Get data
            if isVertical
                xData = repmat(baseValue,1,2);
                yData = get(m2t.current.gca,'ylim');
            else
                xData = get(m2t.current.gca,'xlim');
                yData = repmat(baseValue,1,2);
            end

            [m2t, table, tableOptions] = M2T_makeTable(m2t, '', xData, '', yData);

            % Print out
            drawOpts = M2T_opts_print(drawOptions);
            tabOpts  = M2T_opts_print(tableOptions);
            str      = sprintf('\\addplot[%s] table[%s] {%s};\n', ...
                               drawOpts, tabOpts, table);
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawAreaSeries(m2t, h)
            % Takes care of MATLAB's area plots.
            %
            % TODO Get rid of code duplication with 'LibAv05.M2T_drawAxes'.

            % Do nothing if plot is invisible
            str = '';
            if ~LibAv05.M2T_isVisible(h)
                return;
            end

            % Init drawOptions
            drawOptions = M2T_opts_new();

            % Get the draw options for the bars
            [m2t, drawOptions] = M2T_getPatchDrawOptions(m2t, h, drawOptions);

            if ~isfield(m2t, 'addedAreaOption') || isempty(m2t.addedAreaOption) || ~m2t.addedAreaOption
                % Add 'area style' to axes options.
                m2t = M2T_m2t_addAxisOption(m2t, 'area style');
                m2t = M2T_m2t_addAxisOption(m2t, 'stack plots', 'y');
                m2t.addedAreaOption = true;
            end

            % Toggle legend entry
            drawOptions = M2T_maybeShowInLegend(m2t.currentHandleHasLegend, drawOptions);

            % Generate the tikz table
            xData = get(h, 'XData');
            yData = get(h, 'YData');
            [m2t, table, tableOptions] = M2T_makeTable(m2t, '', xData, '', yData);

            % Print out
            drawOpts = M2T_opts_print(drawOptions);
            tabOpts  = M2T_opts_print(tableOptions);
            str      = sprintf('\\addplot[%s] table[%s]{%s}\n\\closedcycle;\n',...
                               drawOpts, tabOpts, table);
            %TODO: shouldn't this be "\addplot[] table[] {}" instead?
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawStemSeries(m2t, h)
            [m2t, str] = M2T_drawStemOrStairSeries_(m2t, h, 'ycomb');

            % TODO: handle baseplane with stem3()
            if m2t.axes{end}.is3D
                return
            end
            [m2t, baseline] = M2T_drawBaseline(m2t,h);
            str             = [str, baseline];
        end
        function [m2t, str] = M2T_drawStairSeries(m2t, h)
            [m2t, str] = M2T_drawStemOrStairSeries_(m2t, h, 'const plot');
        end
        function [m2t, str] = M2T_drawStemOrStairSeries_(m2t, h, plotType)

            % Do nothing if plot is invisible
            str = '';
            if ~LibAv05.M2T_isLineVisible(h)
                return % nothing to plot!
            end

            % deal with draw options
            color = get(h, 'Color');
            [m2t, plotColor] = M2T_getColor(m2t, h, color, 'patch');

            [m2t, lineOptions]   = M2T_getLineOptions(m2t, h);
            [m2t, markerOptions] = M2T_getMarkerOptions(m2t, h);

            drawOptions = M2T_opts_new();
            drawOptions = M2T_opts_add(drawOptions, plotType);
            drawOptions = M2T_opts_add(drawOptions, 'color', plotColor);
            drawOptions = M2T_opts_merge(drawOptions, lineOptions, markerOptions);

            % Toggle legend entry
            drawOptions = M2T_maybeShowInLegend(m2t.currentHandleHasLegend, drawOptions);

            drawOpts = M2T_opts_print(drawOptions);

            % Generate the tikz table
            xData = get(h, 'XData');
            yData = get(h, 'YData');
            if m2t.axes{end}.is3D
                % TODO: account for hgtransform
                zData = get(h, 'ZData');
                [m2t, table, tableOptions] = M2T_makeTable(m2t, '', xData, '', yData, '', zData);
                % Print out
                tabOpts  = M2T_opts_print(tableOptions);
                str = sprintf('\\addplot3 [%s]\n table[%s] {%s};\n ', ...
                                 drawOpts, tabOpts, table);
            else
                [m2t, table, tableOptions] = M2T_makeTable(m2t, '', xData, '', yData);
                % Print out
                tabOpts  = M2T_opts_print(tableOptions);
                str = sprintf('\\addplot[%s] table[%s] {%s};\n', ...
                                 drawOpts, tabOpts, table);
            end

        end
        % ==============================================================================
        function [m2t, str] = M2T_drawQuiverGroup(m2t, h)
            % Takes care of MATLAB's quiver plots.
            str = '';

            [x,y,z,u,v,w] = M2T_getAndRescaleQuivers(m2t,h);
            is3D = m2t.axes{end}.is3D;

            % prepare output
            if is3D
                name = 'addplot3';
            else % 2D plotting
                name = 'addplot';
            end

            variables = {'x', 'y', 'z', 'u', 'v', 'w'};
            data = NaN(numel(x),6);
            data(:,1) = x;
            data(:,2) = y;
            data(:,3) = z;
            data(:,4) = u;
            data(:,5) = v;
            data(:,6) = w;

            if ~is3D
                data(:,[3 6]) = []; % remove Z-direction
                variables([3 6]) = [];
            end

            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % gather the arrow options
            showArrowHead = get(h, 'ShowArrowHead');
            if ~LibAv05.M2T_isLineVisible(h)  && ~showArrowHead
                return
            end

            plotOptions = M2T_opts_new();
            if showArrowHead
                plotOptions = M2T_opts_add(plotOptions, '-Straight Barb');
                LibAv05.M2T_signalDependency(m2t, 'tikzlibrary', 'arrows.meta');
            else
                plotOptions = M2T_opts_add(plotOptions, '-');
            end

            % Append the arrow style to the TikZ options themselves.
            color = get(h, 'Color');
            [m2t, lineOptions] = M2T_getLineOptions(m2t, h);
            [m2t, arrowcolor] = M2T_getColor(m2t, h, color, 'patch');
            plotOptions = M2T_opts_add(plotOptions, 'color', arrowcolor);
            plotOptions = M2T_opts_merge(plotOptions, lineOptions);

            % Define the quiver settings
            quiverOptions = M2T_opts_new();
            quiverOptions = M2T_opts_add(quiverOptions, 'u', '\thisrow{u}');
            quiverOptions = M2T_opts_add(quiverOptions, 'v', '\thisrow{v}');
            if is3D
                quiverOptions = M2T_opts_add(quiverOptions, 'w', '\thisrow{w}');
                arrowLength = '{sqrt((\thisrow{u})^2+(\thisrow{v})^2+(\thisrow{w})^2)}';
            else
                arrowLength = '{sqrt((\thisrow{u})^2+(\thisrow{v})^2)}';
            end
            plotOptions = M2T_opts_add(plotOptions, 'point meta', arrowLength);
            plotOptions = M2T_opts_add(plotOptions, 'point meta min', '0');

            if showArrowHead
                arrowHeadOptions = M2T_opts_new();

                % In MATLAB (HG1), the arrow head is constructed to have an angle of
                % approximately 18.263 degrees in 2D as can be derived from the
                % |quiver| function.
                % In 3D, the angle is no longer constant but it is approximately
                % the same as for 2D quiver plots. So let's make our life easy.
                % |LibAv05.M2T_test/examples/example_quivers.m| covers the calculations.
                arrowHeadOptions = M2T_opts_add(arrowHeadOptions, 'angle''', '18.263');

                %TODO: scale the arrows more rigorously to match MATLAB behavior
                % Currently, this is quite hard to do, since the size of the arrows
                % is defined in pgfplots in absolute units (here we specify that those
                % should be scaled up/down according to the data), while the data itself
                % is in axis coordinates (or some scaled variant). I.e. we need the
                % physical dimensions of the axis to compute the correct scaling!
                %
                % There is a "MaxHeadSize" property that plays a role.
                % MaxHeadSize is said to be relative to the length of the quiver in the
                % MATLAB documentation. However, in practice, there seems to be a SQRT
                % involved somewhere (e.g. if u.^2 + v.^2 == 2, all MHS values >
                % 1/sqrt(2) are capped to 1/sqrt(2)).
                %
                % NOTE: `set(h, 'MaxHeadSize')` is bugged in HG1 (not in HG2 or Octave)
                % according to http://www.mathworks.com/matlabcentral/answers/96754

                LibAv05.M2T_userInfo(m2t, ['Please change the "arrowHeadSize" option', ...
                    ' if the size of the arrows is incorrect.']);
                arrowHeadSize = sprintf(m2t.ff, abs(m2t.args.arrowHeadSize));

                % Write out the actual scaling for TikZ.
                % `\pgfplotspointsmetatransformed` is in the range [0, 1000], so
                % divide by this span (as is done in the pgfplots manual) to normalize
                % the arrow head size. First divide to avoid overflows.
                arrowHeadOptions = M2T_opts_add(arrowHeadOptions, 'scale', ...
                    ['{' arrowHeadSize '/1000*\pgfplotspointmetatransformed}']);

                headStyle = ['-{Straight Barb[' LibAv05.M2T_opts_print(arrowHeadOptions) ']}'];
                quiverOptions = M2T_opts_add(quiverOptions, 'every arrow/.append style', ...
                                         ['{' headStyle '}']);
            end
            plotOptions = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(plotOptions, 'quiver', quiverOptions);

            [m2t, table, tableOptions] = M2T_makeTable(m2t, variables, data);

            % Print out
            plotOpts = M2T_opts_print(plotOptions);
            tabOpts  = M2T_opts_print(tableOptions);
            str      = sprintf('\\%s[%s]\n table[%s] {%s};\n', ...
                               name, plotOpts, tabOpts, table);
        end
        % ==============================================================================
        function [x,y,z,u,v,w] = M2T_getAndRescaleQuivers(m2t, h)
            % get and rescale the arrows from a quivergroup object
            x = get(h, 'XData');
            y = get(h, 'YData');
            z = M2T_getOrDefault(h, 'ZData', []);

            u = get(h, 'UData');
            v = get(h, 'VData');
            w = M2T_getOrDefault(h, 'WData', []);

            is3D = m2t.axes{end}.is3D;
            if ~is3D
                z = 0;
                w = 0;
            end

            % MATLAB uses a scaling algorithm to determine the size of the arrows.
            % Before R2014b, the processed coordinates were available. This is no longer
            % the case, so we have to re-implement it. In MATLAB it is implemented in
            % the |quiver3|  (and |quiver|) function.
            if any(size(x)==1)
                nX = sqrt(numel(x)); nY = nX;
            else
                [nY, nX] = size(x);
            end
            range  = @(xyzData)(max(xyzData(:)) - min(xyzData(:)));
            euclid = @(x,y,z)(sqrt(x.^2 + y.^2 + z.^2));
            dx = range(x)/nX;
            dy = range(y)/nY;
            dz = range(z)/max(nX,nY);
            dd = euclid(dx, dy, dz);
            if dd > 0
                vectorLength = euclid(u/dd,v/dd,w/dd);
                maxLength = max(vectorLength(:));
            else
                maxLength = 1;
            end
            if LibAv05.M2T_isOn(LibAv05.M2T_getOrDefault(h, 'AutoScale', 'on'))
                scaleFactor = M2T_getOrDefault(h,'AutoScaleFactor', 0.9) / maxLength;
            else
                scaleFactor = 1;
            end
            x = x(:).'; u = u(:).'*scaleFactor;
            y = y(:).'; v = v(:).'*scaleFactor;
            z = z(:).'; w = w(:).'*scaleFactor;
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawErrorBars(m2t, h)
            % Takes care of MATLAB's error bar plots.
            % Octave's error bar plots are handled as well.
            [m2t, str] = M2T_drawLine(m2t, h);
            % Even though this only calls |LibAv05.M2T_drawLine|, let's keep this wrapper
            % such that the code is easier to read where it is called.
        end
        % ==============================================================================
        function [yDeviations] = M2T_getYDeviations(h)
            % Retrieves upper/lower uncertainty data

            upDev = M2T_getOrDefault(h, 'UData', []);
            loDev = M2T_getOrDefault(h, 'LData', []);

            yDeviations = [upDev(:), loDev(:)];
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawEllipse(m2t, handle)
            % Takes care of MATLAB's ellipse annotations.

            drawOptions = M2T_opts_new();

            p = get(handle,'position');
            radius = p([3 4]) / 2;
            center = p([1 2]) + radius;

            color = get(handle, 'Color');
            [m2t, xcolor] = M2T_getColor(m2t, handle, color, 'patch');
            [m2t, lineOptions] = M2T_getLineOptions(m2t, handle);

            filling = get(handle, 'FaceColor');

            % Has a filling?
            if LibAv05.M2T_isNone(filling)
                drawOptions = M2T_opts_add(drawOptions, xcolor);
                drawCommand = '\draw';
            else
                [m2t, xcolorF] = M2T_getColor(m2t, handle, filling, 'patch');
                drawOptions = M2T_opts_add(drawOptions, 'draw', xcolor);
                drawOptions = M2T_opts_add(drawOptions, 'fill', xcolorF);

                drawCommand = '\filldraw';
            end
            drawOptions = M2T_opts_merge(drawOptions, lineOptions);

            opt = M2T_opts_print(drawOptions);

            str = sprintf('%s [%s] (axis cs:%g,%g) ellipse [x radius=%g, y radius=%g];\n', ...
                drawCommand, opt, center, radius);
        end
        % ==============================================================================
        function [m2t, str] = M2T_drawTextarrow(m2t, handle)
            % Takes care of MATLAB's textarrow annotations.

            % LibAv05.M2T_handleAllChildren to draw the arrow
            [m2t, str] = M2T_handleAllChildren(m2t, handle);

            % LibAv05.M2T_handleAllChildren ignores the text, unless hidden strings are shown
            if ~m2t.args.showHiddenStrings
                child = findall(handle, 'type', 'text');
                [m2t, str{end+1}] = M2T_drawText(m2t, child);
            end
        end
        % ==============================================================================
        function [m2t, drawOptions] = M2T_getPatchDrawOptions(m2t, h, drawOptions)
            % Determines the reoccurring draw options usually applied when drawing
            % a patch/area/bar. These include EdgeColor, LineType, FaceColor/Alpha

            % Get object for color;
            if ~isempty(allchild(h))
                % quite oddly, before MATLAB R2014b this value is stored in a child
                % patch and not in the object itself
                obj = allchild(h);
            else % R2014b and newer
                obj = h;
            end

            % Get the object type
            type = get(h, 'Type');

            % Face Color (inside of area)
            faceColor          = get(obj, 'FaceColor');
            [m2t, drawOptions] = M2T_setColor(m2t, h, drawOptions, 'fill', faceColor, 'none');

            % FaceAlpha (Not applicable for MATLAB2014a/b)
            faceAlpha = M2T_getOrDefault(h, 'FaceAlpha', 'none');
            if ~LibAv05.M2T_isNone(faceColor) && isnumeric(faceAlpha) && faceAlpha ~= 1.0
                drawOptions = M2T_opts_add(drawOptions, 'fill opacity', sprintf(m2t.ff,faceAlpha));
            end

            % Define linestyle
            [lineStyle, isDefaultLS] = M2T_getAndCheckDefault(type, h, 'LineStyle', '-');
            if LibAv05.M2T_isNone(lineStyle)
                drawOptions = M2T_opts_add(drawOptions, 'draw', 'none');
            elseif ~isDefaultLS
                drawOptions = M2T_opts_add(drawOptions, LibAv05.M2T_translateLineStyle(lineStyle));
            end

            % Check for the edge color. Only plot it if it is different from the
            % face color and if there is a linestyle
            edgeColor = get(h, 'EdgeColor');
            if ~LibAv05.M2T_isNone(lineStyle) && ~LibAv05.M2T_isNone(edgeColor) && ~strcmpi(edgeColor,faceColor)
                [m2t, drawOptions] = M2T_setColor(m2t, h, drawOptions, 'draw', edgeColor, 'none');
            end

            % Add 'area legend' to the options as otherwise the legend indicators
            % will just highlight the edges.
            if strcmpi(type, 'bar') || strcmpi(type, 'histogram')
                drawOptions = M2T_opts_add(drawOptions, 'area legend');
            end
        end
        % ==============================================================================
        function out = M2T_linearFunction(X, Y)
            % Return the linear function that goes through (X[1], Y[1]), (X[2], Y[2]).
            out = @(x) (Y(2,:)*(x-X(1)) + Y(1,:)*(X(2)-x)) / (X(2)-X(1));
        end
        % ==============================================================================
        function matlabColormap = M2T_pgfplots2matlabColormap(points, rgb, numColors)
            % Translates a Pgfplots colormap to a MATLAB color map.
            matlabColormap = zeros(numColors, 3);
            % Point indices between which to interpolate.
            I = [1, 2];
            f = M2T_linearFunction(points(I), rgb(I,:));
            for k = 1:numColors
                x = (k-1)/(numColors-1) * points(end);
                if x > points(I(2))
                    I = I + 1;
                    f = M2T_linearFunction(points(I), rgb(I,:));
                end
                matlabColormap(k,:) = f(x);
            end
        end
        % ==============================================================================
        function pgfplotsColormap = M2T_matlab2pgfplotsColormap(m2t, matlabColormap, name)
            % Translates a MATLAB color map into a Pgfplots colormap.

            if nargin < 3 || isempty(name), name = 'mymap'; end

            % First check if we could use a default Pgfplots color map.
            % Unfortunately, MATLAB and Pgfplots color maps will never exactly coincide
            % except to the most simple cases such as blackwhite. This is because of a
            % slight incompatibility of Pgfplots and MATLAB colormaps:
            % In MATLAB, indexing goes from 1 through 64, whereas in Pgfplots you can
            % specify any range, the default ones having something like
            % (0: red, 1: yellow, 2: blue).
            % To specify this exact color map in MATLAB, one would have to put 'red' at
            % 1, blue at 64, and yellow in the middle of the two, 32.5 that is.
            % Not really sure how MATLAB rounds here: 32, 33? Anyways, it will be
            % slightly off and hence not match the Pgfplots color map.
            % As a workaround, build the MATLAB-formatted colormaps of Pgfplots default
            % color maps, and check if matlabColormap is close to it. If yes, take it.

            % For now, comment out the color maps which haven't landed yet in Pgfplots.
            pgfmaps = { %struct('name', 'colormap/autumn', ...
                %       'points', [0,1], ...
                %       'values', [[1,0,0];[1,1,0]]), ...
                %struct('name', 'colormap/bled', ...
                %       'points', 0:6, ...
                %       'values', [[0,0,0];[43,43,0];[0,85,0];[0,128,128];[0,0,170];[213,0,213];[255,0,0]]/255), ...
                %struct('name', 'colormap/bright', ...
                %       'points', 0:7, ...
                %       'values', [[0,0,0];[78,3,100];[2,74,255];[255,21,181];[255,113,26];[147,213,114];[230,255,0];[255,255,255]]/255), ...
                %struct('name', 'colormap/bone', ...
                %       'points', [0,3,6,8], ...
                %       'values', [[0,0,0];[84,84,116];[167,199,199];[255,255,255]]/255), ...
                %struct('name', 'colormap/cold', ...
                %       'points', 0:3, ...
                %       'values', [[0,0,0];[0,0,1];[0,1,1];[1,1,1]]), ...
                %struct('name', 'colormap/copper', ...
                %       'points', [0,4,5], ...
                %       'values', [[0,0,0];[255,159,101];[255,199,127]]/255), ...
                %struct('name', 'colormap/copper2', ...
                %       'points', 0:4, ...
                %       'values', [[0,0,0];[68,62,63];[170,112,95];[207,194,138];[255,255,255]]/255), ...
                %struct('name', 'colormap/hsv', ...
                %       'points', 0:6, ...
                %       'values', [[1,0,0];[1,1,0];[0,1,0];[0,1,1];[0,0,1];[1,0,1];[1,0,0]]), ...
                struct('name', 'colormap/hot', ...
                'points', 0:3, ...
                'values', [[0,0,1];[1,1,0];[1,0.5,0];[1,0,0]]), ... % TODO check this
                struct('name', 'colormap/hot2', ...
                'points', [0,3,6,8], ...
                'values', [[0,0,0];[1,0,0];[1,1,0];[1,1,1]]), ...
                struct('name', 'colormap/jet', ...
                'points', [0,1,3,5,7,8], ...
                'values', [[0,0,128];[0,0,255];[0,255,255];[255,255,0];[255,0,0];[128,0,0]]/255), ...
                struct('name', 'colormap/blackwhite', ...
                'points', [0,1], ...
                'values', [[0,0,0];[1,1,1]]), ...
                struct('name', 'colormap/bluered', ...
                'points', 0:5, ...
                'values', [[0,0,180];[0,255,255];[100,255,0];[255,255,0];[255,0,0];[128,0,0]]/255), ...
                struct('name', 'colormap/cool', ...
                'points', [0,1,2], ...
                'values', [[255,255,255];[0,128,255];[255,0,255]]/255), ...
                struct('name', 'colormap/greenyellow', ...
                'points', [0,1], ...
                'values', [[0,128,0];[255,255,0]]/255), ...
                struct('name', 'colormap/redyellow', ...
                'points', [0,1], ...
                'values', [[255,0,0];[255,255,0]]/255), ...
                struct('name', 'colormap/violet', ...
                'points', [0,1,2], ...
                'values', [[25,25,122];[255,255,255];[238,140,238]]/255) ...
                };

            % The tolerance is a subjective matter of course.
            % Some figures:
            %    * The norm-distance between MATLAB's gray and bone is 6.8e-2.
            %    * The norm-distance between MATLAB's jet and Pgfplots's jet is 2.8e-2.
            %    * The norm-distance between MATLAB's hot and Pgfplots's hot2 is 2.1e-2.
            tol = 5.0e-2;

            for map = pgfmaps
                numColors = size(matlabColormap, 1);
                mmap = M2T_pgfplots2matlabColormap(map{1}.points, map{1}.values, numColors);
                alpha = norm(matlabColormap - mmap) / sqrt(numColors);
                if alpha < tol
                    LibAv05.M2T_userInfo(m2t, 'Found %s to be a pretty good match for your color map (||diff||=%g).', ...
                        map{1}.name, alpha);
                    pgfplotsColormap = map{1}.name;
                    return
                end
            end

            % Build a custom color map.
            % Loop over the data, stop at each spot where the linear
            % interpolation is interrupted, and set a color mark there.
            m = size(matlabColormap, 1);
            steps = [1, 2];
            % A colormap with a single color is valid in MATLAB but an error in
            % pgfplots. Repeating the color produces the desired effect in this
            % case.
            if m==1
                colors=[matlabColormap(1,:);matlabColormap(1,:)];
            else
                colors = [matlabColormap(1,:); matlabColormap(2,:)];
                f = M2T_linearFunction(steps, colors);
                k = 3;
                while k <= m
                    if norm(matlabColormap(k,:) - f(k)) > 1.0e-10
                        % Add the previous step to the color list
                        steps(end) = k-1;
                        colors(end,:) = matlabColormap(k-1,:);
                        steps = [steps, k];
                        colors = [colors; matlabColormap(k,:)];
                        f = M2T_linearFunction(steps(end-1:end), colors(end-1:end,:));
                    end
                    k = k+1;
                end
                steps(end) = m;
                colors(end,:) = matlabColormap(m,:);
            end

            % Get it in Pgfplots-readable form.
            unit = 'pt';
            colSpecs = cell(length(steps), 1);
            for k = 1:length(steps)
                x = steps(k)-1;
                colSpecs{k} = sprintf('rgb(%d%s)=(%g,%g,%g)', x, unit, colors(k,:));
            end
            pgfplotsColormap = sprintf('colormap={%s}{[1%s] %s}',name, unit, LibAv05.M2T_join(m2t, colSpecs, '; '));
        end
        % ==============================================================================
        function [m2t, fontStyle] = M2T_getFontStyle(m2t, handle)
            fontStyle = '';
            if strcmpi(get(handle, 'FontWeight'),'Bold')
                fontStyle = sprintf('%s\\bfseries',fontStyle);
            end
            if strcmpi(get(handle, 'FontAngle'), 'Italic')
                fontStyle = sprintf('%s\\itshape',fontStyle);
            end
            if ~all(get(handle, 'Color')==0)
                color = get(handle, 'Color');
                [m2t, col] = M2T_getColor(m2t, handle, color, 'patch');
                fontStyle = sprintf('%s\\color{%s}', fontStyle, col);
            end
            if m2t.args.strictFontSize
                fontSize  = get(handle,'FontSize');
                fontUnits = M2T_matlab2texUnits(get(handle,'FontUnits'), 'pt');
                fontStyle = sprintf('\\fontsize{%d%s}{1em}%s\\selectfont',fontSize,fontUnits,fontStyle);
            else
                % don't try to be smart and "translate" MATLAB font sizes to proper LaTeX
                % ones: it cannot be done. LaTeX uses semantic sizes (e.g. \small)
                % whose actual dimensions depend on the document style, context, ...
            end

            if ~isempty(fontStyle)
                fontStyle = M2T_opts_add(LibAv05.M2T_opts_new, 'font', fontStyle);
            else
                fontStyle = M2T_opts_new();
            end
        end
        % ==============================================================================
        function axisOptions = M2T_getColorbarOptions(m2t, handle)

            % begin collecting axes options
            axisOptions = M2T_opts_new();
            cbarStyleOptions = M2T_opts_new();

            [cbarTemplate, cbarStyleOptions] = LibAv05.M2T_LibAv05.M2T_getColorbarPosOptions(handle, ...
                                                        cbarStyleOptions);

            % axis label and direction
            if LibAv05.M2T_isHG2
                % VERSION: Starting from R2014b there is only one field `label`.
                % The colorbar's position determines, if it should be a x- or y-label.

                if strcmpi(cbarTemplate, 'horizontal')
                    labelOption = 'xlabel';
                else
                    labelOption = 'ylabel';
                end
                [m2t, cbarStyleOptions] = M2T_getLabel(m2t, handle, cbarStyleOptions, labelOption);

                % direction
                dirString = get(handle, 'Direction');
                if ~strcmpi(dirString, 'normal') % only if not 'normal'
                    if strcmpi(cbarTemplate, 'horizontal')
                        dirOption = 'x dir';
                    else
                        dirOption = 'y dir';
                    end
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, dirOption, dirString);
                end

                % TODO HG2: colorbar ticks and colorbar tick labels

            else
                % VERSION: Up to MATLAB R2014a and OCTAVE
                [m2t, xo] = M2T_getAxisOptions(m2t, handle, 'x');
                [m2t, yo] = M2T_getAxisOptions(m2t, handle, 'y');
                xyo = M2T_opts_merge(xo, yo);
                xyo = M2T_opts_remove(xyo, 'xmin','xmax','xtick','ymin','ymax','ytick');

                cbarStyleOptions = M2T_opts_merge(cbarStyleOptions, xyo);
            end

            % title
            [m2t, cbarStyleOptions] = M2T_getTitle(m2t, handle, cbarStyleOptions);

            if m2t.args.strict
                % Sampled colors.
                numColors = size(m2t.current.colormap, 1);
                axisOptions = M2T_opts_add(axisOptions, 'colorbar sampled');
                cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'samples', ...
                    sprintf('%d', numColors+1));

                if ~isempty(cbarTemplate)
                    LibAv05.M2T_userWarning(m2t, ...
        -               'Pgfplots cannot deal with more than one colorbar option yet.');
                    %FIXME: can we get sampled horizontal color bars to work?
                    %FIXME: sampled colorbars should be inferred, not by using strict!
                end
            end

            % Merge them together in axisOptions.
            axisOptions = M2T_opts_add(axisOptions, strtrim(['colorbar ', cbarTemplate]));

            if ~isempty(cbarStyleOptions)
                axisOptions = LibAv05.M2T_LibAv05.M2T_opts_addSubOpts(axisOptions, ...
                                           'colorbar style', cbarStyleOptions);
            end

            % do _not_ handle colorbar's children
        end
        % ==============================================================================
        function [cbarTemplate, cbarStyleOptions] = M2T_getColorbarPosOptions(handle, cbarStyleOptions)
            % set position, ticks etc. of a colorbar
            loc = get(handle, 'Location');
            cbarTemplate = '';

            switch lower(loc) % case insensitive (MATLAB: CamelCase, Octave: lower case)
                case 'north'
                    cbarTemplate = 'horizontal';
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'at',...
                        '{(0.5,0.97)}');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'anchor',...
                        'north');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions,...
                        'xticklabel pos', 'lower');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'width',...
                        '0.97*\pgfkeysvalueof{/pgfplots/parent axis width}');
                case 'south'
                    cbarTemplate = 'horizontal';
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'at',...
                        '{(0.5,0.03)}');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'anchor', ...
                        'south');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, ...
                        'xticklabel pos','upper');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'width',...
                        '0.97*\pgfkeysvalueof{/pgfplots/parent axis width}');
                case 'east'
                    cbarTemplate = 'right';
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'at',...
                        '{(0.97,0.5)}');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'anchor', ...
                        'east');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, ...
                        'xticklabel pos','left');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'width',...
                        '0.97*\pgfkeysvalueof{/pgfplots/parent axis width}');
                case 'west'
                    cbarTemplate = 'left';
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'at',...
                        '{(0.03,0.5)}');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'anchor',...
                        'west');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions,...
                        'xticklabel pos', 'right');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'width',...
                        '0.97*\pgfkeysvalueof{/pgfplots/parent axis width}');
                case 'eastoutside'
                    %cbarTemplate = 'right';
                case 'westoutside'
                    cbarTemplate = 'left';
                case 'northoutside'
                    % TODO move to top
                    cbarTemplate = 'horizontal';
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'at',...
                        '{(0.5,1.03)}');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'anchor',...
                        'south');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions,...
                        'xticklabel pos', 'upper');
                case 'southoutside'
                    cbarTemplate = 'horizontal';
                case 'manual'
                    origUnits = get(handle,'Units');
                    assocAxes = get(handle,'Axes');
                    origAxesUnits = get(assocAxes,'Units');
                    set(handle,'Units','centimeters');        % Make sure we have
                    set(assocAxes,'Units','centimeters');     % same units
                    cbarDim = M2T_pos2dims(get(handle,'Position'));
                    cbarAxesDim = M2T_pos2dims(get(assocAxes,'Position'));
                    set(handle,'Units',origUnits);            % Restore original
                    set(assocAxes,'Units',origAxesUnits);     % units

                    center = @(dims) (dims.left + dims.right)/2;
                    centerCbar = center(cbarDim);
                    centerAxes = center(cbarAxesDim);

                    % Cases of colorbar axis locations (in or out) depending on center
                    % of colorbar relative to the center it's associated axes.
                    % According to matlab manual (R2016a) colorbars with Location 'manual'
                    % can only be vertical.
                    axisLoc = M2T_getOrDefault(handle, 'AxisLocation', 'out');
                    if centerCbar < centerAxes
                        if strcmp(axisLoc,'in')
                            cbarTemplate = 'right';
                        else
                            cbarTemplate = 'left';
                        end
                    else
                        if strcmp(axisLoc,'in')
                            cbarTemplate = 'left';
                        else
                            cbarTemplate = 'right';
                        end
                    end

                    % Using positions relative to associated axes
                    calcRelPos = @(pos1,pos2,ext2) (pos1-pos2)/ext2; 
                    cbarRelPosX = calcRelPos(cbarDim.left,cbarAxesDim.left,cbarAxesDim.width);
                    cbarRelPosY = calcRelPos(cbarDim.bottom,cbarAxesDim.bottom,cbarAxesDim.height);
                    cbarRelHeight = cbarDim.height/cbarAxesDim.height;

                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'anchor',...
                        'south west');
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'at',...
                        ['{(' LibAv05.M2T_formatDim(cbarRelPosX) ','...
                              LibAv05.M2T_formatDim(cbarRelPosY) ')}']);
                    cbarStyleOptions = M2T_opts_add(cbarStyleOptions, 'height',...
                        [LibAv05.M2T_formatDim(cbarRelHeight),...
                        '*\pgfkeysvalueof{/pgfplots/parent axis height}']);

                otherwise
                    error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_getColorOptions:unknownLocation',...
                        'LibAv05.M2T_LibAv05.M2T_getColorbarOptions: Unknown ''Location'' %s.', loc)
            end
        end
        % ==============================================================================
        function [m2t, xcolor] = M2T_getColor(m2t, handle, color, mode)
            % Handles MATLAB colors and makes them available to TikZ.
            % This includes translation of the color value as well as explicit
            % definition of the color if it is not available in TikZ by default.
            %
            % The variable 'mode' essentially determines what format 'color' can
            % have. Possible values are (as strings) 'patch' and 'image'.

            % check if the color is straight given in rgb
            % -- notice that we need the extra NaN LibAv05.M2T_test with respect to the QUIRK
            %    below
            if LibAv05.M2T_isRGBTuple(color)
                % everything alright: rgb color here
                [m2t, xcolor] = M2T_rgb2colorliteral(m2t, color);
            else
                switch lower(mode)
                    case 'patch'
                        [m2t, xcolor] = M2T_patchcolor2xcolor(m2t, color, handle);
                    case 'image'

                        m = size(color,1);
                        n = size(color,2);
                        xcolor = cell(m, n);

                        if ndims(color) == 3
                            for i = 1:m
                                for j = 1:n
                                    [m2t, xc] = M2T_rgb2colorliteral(m2t, color(i,j, :));
                                    xcolor{i, j} = xc;
                                end
                            end
                        elseif ndims(color) <= 2
                            [m2t, colorindex] = M2T_cdata2colorindex(m2t, color, handle);
                            for i = 1:m
                                for j = 1:n
                                    [m2t, xc] = M2T_rgb2colorliteral(m2t, m2t.current.colormap(colorindex(i,j), :));
                                    xcolor{i, j} = xc;
                                end
                            end
                        else
                            error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_getColor:image:colorDims',...
                                'Image color data cannot have more than 3 dimensions');
                        end
                    otherwise
                        error(['LibAv05.M2T_matlab2tikz:LibAv05.M2T_getColor', ...
                            'Argument ''mode'' has illegal value ''%s''.'], ...
                            mode);
                end
            end
        end
        % ==============================================================================
        function [m2t, xcolor] = M2T_patchcolor2xcolor(m2t, color, patchhandle)
            % Transforms a color of the edge or the face of a patch to an xcolor literal.
            if isnumeric(color)
                [m2t, xcolor] = M2T_rgb2colorliteral(m2t, color);
            elseif ischar(color)
                switch color
                    case 'flat'
                        cdata  = M2T_getCDataWithFallbacks(patchhandle);
                        color1 = cdata(1,1);
                        % RGB cdata
                        if ndims(cdata) == 3 && all(size(cdata) == [1,1,3])
                            [m2t,xcolor] = M2T_rgb2colorliteral(m2t, cdata);
                            % All same color
                        elseif all(isnan(cdata) | abs(cdata-color1)<1.0e-10)
                            [m2t, colorindex] = M2T_cdata2colorindex(m2t, color1, patchhandle);
                            [m2t, xcolor] = M2T_rgb2colorliteral(m2t, m2t.current.colormap(colorindex, :));
                        else
                            % Don't return anything meaningful and count on the caller
                            % to make something of it.
                            xcolor = [];
                        end

                    case 'auto'
                        try
                            color = get(patchhandle, 'Color');
                        catch
                            % From R2014b use an undocumented property if Color is
                            % not present
                            color = get(patchhandle, 'AutoColor');
                        end
                        [m2t, xcolor] = M2T_rgb2colorliteral(m2t, color);

                    case 'none'
                        % Before, we used to throw an error here. However, probably this
                        % is not necessary and actually harmful (#739).
                        xcolor = 'none';

                    otherwise
                        error('LibAv05.M2T_matlab2tikz:anycolor2rgb:UnknownColorModel',...
                        'Don''t know how to handle the color model ''%s''.',color);
                end
            else
                error('LibAv05.M2T_patchcolor2xcolor:illegalInput', ...
                    'Input argument ''color'' not a string or numeric.');
            end
        end
        % ==============================================================================
        function cdata = M2T_getCDataWithFallbacks(patchhandle)
            % Looks for CData at different places
            cdata = M2T_getOrDefault(patchhandle, 'CData', []);

            if isempty(cdata) || ~isnumeric(cdata)
                child = allchild(patchhandle);
                cdata = get(child, 'CData');
            end
            if isempty(cdata) || ~isnumeric(cdata)
                % R2014b+: CData is implicit by the ordering of the siblings
                siblings = allchild(get(patchhandle, 'Parent'));
                cdata = find(siblings(end:-1:1)==patchhandle);
            end
        end
        % ==============================================================================
        function [m2t, colorindex] = M2T_cdata2colorindex(m2t, cdata, imagehandle)
            % Transforms a color in CData format to an index in the color map.
            % Only does something if CDataMapping is 'scaled', really.

            if ~isnumeric(cdata) && ~islogical(cdata)
                error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_cdata2colorindex:unknownCDataType',...
                    'Don''t know how to handle CData ''%s''.',cdata);
            end

            axeshandle = m2t.current.gca;

            % -----------------------------------------------------------------------
            % For the following, see, for example, the MATLAB help page for 'image',
            % section 'Image CDataMapping'.
            try
                mapping = get(imagehandle, 'CDataMapping');
            catch
                mapping = 'scaled';
            end
            switch mapping
                case 'scaled'
                    % need to scale within clim
                    % see MATLAB's manual page for caxis for details
                    clim = get(axeshandle, 'clim');
                    m = size(m2t.current.colormap, 1);
                    colorindex = zeros(size(cdata));
                    idx1 = cdata <= clim(1);
                    idx2 = cdata >= clim(2);
                    idx3 = ~idx1 & ~idx2;
                    colorindex(idx1) = 1;
                    colorindex(idx2) = m;
                    % cdata may be of type uint8. Convert to double to avoid
                    % getting binary indices
                    colorindex(idx3) = fix(double(cdata(idx3)-clim(1)) / (clim(2)-clim(1)) *m) ...
                        + 1;
                case 'direct'
                    % direct index
                    colorindex = cdata;

                otherwise
                    error('LibAv05.M2T_matlab2tikz:anycolor2rgb:unknownCDataMapping',...
                        'Unknown CDataMapping ''%s''.',cdatamapping);
            end
        end
        % ==============================================================================
        function [m2t, key, legendOpts] = M2T_getLegendOpts(m2t, handle)
            lStyle = M2T_opts_new();

            lStyle = M2T_getLegendPosition(m2t, handle, lStyle);
            lStyle = M2T_getLegendOrientation(m2t, handle, lStyle);
            lStyle = M2T_getLegendEntryAlignment(m2t, handle, lStyle);

            % If the plot has 'legend boxoff', we have the 'not visible'
            % property, so turn off line and background fill.
            if ~LibAv05.M2T_isVisible(handle) || LibAv05.M2T_isOff(get(handle,'box'))
                lStyle = M2T_opts_add(lStyle, 'fill', 'none');
                lStyle = M2T_opts_add(lStyle, 'draw', 'none');
            else
                % handle colors
                [edgeColor, isDfltEdge] = M2T_getAndCheckDefault('Legend', handle, ...
                                                             'EdgeColor', [1 1 1]);
                if LibAv05.M2T_isNone(edgeColor)
                    lStyle = M2T_opts_add(lStyle, 'draw', 'none');

                elseif ~isDfltEdge
                    [m2t, col] = M2T_getColor(m2t, handle, edgeColor, 'patch');
                    lStyle = M2T_opts_add(lStyle, 'draw', col);
                end

                [fillColor, isDfltFill] = M2T_getAndCheckDefault('Legend', handle, ...
                                                             'Color', [1 1 1]);
                if LibAv05.M2T_isNone(fillColor)
                    lStyle = M2T_opts_add(lStyle, 'fill', 'none');

                elseif ~isDfltFill
                    [m2t, col] = M2T_getColor(m2t, handle, fillColor, 'patch');
                    lStyle = M2T_opts_add(lStyle, 'fill', col);
                end
            end
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            key = 'legend style';
            legendOpts = M2T_opts_print(lStyle);
            legendOpts = ['{', legendOpts, '}'];
            %TODO: just pass out the `lStyle` instead of `legendOpts`
        end
        % ==============================================================================
        function [lStyle] = M2T_getLegendOrientation(m2t, handle, lStyle)
            % handle legend orientation
            ori = get(handle, 'Orientation');
            switch lower(ori)
                case 'horizontal'
                    numLegendEntries = sprintf('%d',length(get(handle, 'String')));
                    lStyle = M2T_opts_add(lStyle, 'legend columns', numLegendEntries);

                case 'vertical'
                    % Use default.
                otherwise
                    LibAv05.M2T_userWarning(m2t, [' Unknown legend orientation ''',ori,'''' ...
                        '. Choosing default (vertical).']);
            end
        end
        % ==============================================================================
        function [lStyle] = M2T_getLegendPosition(m2t, handle, lStyle)
            % handle legend location
            % #COMPLEX: just a big switch-case
            loc  = get(handle, 'Location');
            dist = 0.03;  % distance to to axes in normalized coordinates
            % MATLAB(R)'s keywords are camel cased (e.g., 'NorthOutside'), in Octave
            % small cased ('northoutside'). Hence, use lower() for uniformity.
            switch lower(loc)
                case 'northeast'
                    return % don't do anything in this (default) case
                case 'northwest'
                    position = [dist, 1-dist];
                    anchor   = 'north west';
                case 'southwest'
                    position = [dist, dist];
                    anchor   = 'south west';
                case 'southeast'
                    position = [1-dist, dist];
                    anchor   = 'south east';
                case 'north'
                    position = [0.5, 1-dist];
                    anchor   = 'north';
                case 'east'
                    position = [1-dist, 0.5];
                    anchor   = 'east';
                case 'south'
                    position = [0.5, dist];
                    anchor   = 'south';
                case 'west'
                    position = [dist, 0.5];
                    anchor   = 'west';
                case 'northoutside'
                    position = [0.5, 1+dist];
                    anchor = 'south';
                case 'southoutside'
                    position = [0.5, -dist];
                    anchor = 'north';
                case 'eastoutside'
                    position = [1+dist, 0.5];
                    anchor = 'west';
                case 'westoutside'
                    position = [-dist, 0.5];
                    anchor = 'east';
                case 'northeastoutside'
                    position = [1+dist, 1];
                    anchor = 'north west';
                case 'northwestoutside'
                    position = [-dist, 1];
                    anchor = 'north east';
                case 'southeastoutside'
                    position = [1+dist, 0];
                    anchor = 'south west';
                case 'southwestoutside'
                    position = [-dist, 0];
                    anchor = 'south east';
                case 'none'
                    legendPos = get(handle, 'Position');
                    unit = get(handle, 'Units');
                    if isequal(unit, 'normalized')
                        position = legendPos(1:2);
                    else
                        % Calculate where the legend is located w.r.t. the axes.
                        axesPos = get(m2t.current.gca, 'Position');
                        axesUnit = get(m2t.current.gca, 'Units');
                        % Convert to legend unit
                        axesPos = M2T_convertUnits(axesPos, axesUnit, unit);
                        % By default, the axes position is given w.r.t. to the figure,
                        % and so is the legend.
                        position = (legendPos(1:2)-axesPos(1:2)) ./ axesPos(3:4);
                    end
                    anchor = 'south west';
                case {'best','bestoutside'}
                    % TODO: Implement these.
                    % The position could be determined by means of 'Position' and/or
                    % 'OuterPosition' of the legend handle; in fact, this could be made
                    % a general principle for all legend placements.
                    LibAv05.M2T_userWarning(m2t, [sprintf(' Option ''%s'' not yet implemented.',loc),         ...
                        ' Choosing default.']);
                    return % use defaults

                otherwise
                    LibAv05.M2T_userWarning(m2t, [' Unknown legend location ''',loc,''''           ...
                        '. Choosing default.']);
                    return % use defaults
            end

            % set legend position
            %TODO: shouldn't this include units?
            lStyle = M2T_opts_add(lStyle, 'at',  sprintf('{(%s,%s)}', ...
                                LibAv05.M2T_formatDim(position(1)), LibAv05.M2T_formatDim(position(2))));
            lStyle = M2T_opts_add(lStyle, 'anchor', anchor);

        end
        % ==============================================================================
        function [lStyle] = M2T_getLegendEntryAlignment(m2t, handle, lStyle)
            % determines the text and picture alignment inside a legend
            textalign = '';
            pictalign = '';
            switch getEnvironment
                case 'Octave'
                    % Octave allows to change the alignment of legend text and
                    % pictograms using legend('left') and legend('right')
                    textpos = get(handle, 'textposition');
                    switch lower(textpos)
                        case 'left'
                            % pictogram right of flush right text
                            textalign = 'right';
                            pictalign = 'right';
                        case 'right'
                            % pictogram left of flush left text (default)
                            textalign = 'left';
                            pictalign = 'left';
                        otherwise
                            LibAv05.M2T_userWarning(m2t, ...
                                ['Unknown legend text position ''',...
                                textpos, '''. Choosing default.']);
                    end
                case 'MATLAB'
                    % does not specify text/pictogram alignment in legends
                otherwise
                    errorUnknownEnvironment();
            end

            % set alignment of legend text and pictograms, if available
            if ~isempty(textalign) && ~isempty(pictalign)
                lStyle = M2T_opts_add(lStyle, 'legend cell align', textalign);
                lStyle = M2T_opts_add(lStyle, 'align', textalign);
                lStyle = M2T_opts_add(lStyle, 'legend plot pos', pictalign);
            else
                % Make sure the entries are flush left (default MATLAB behavior).
                % This is also import for multiline legend entries: Without alignment
                % specification, the TeX document won't compile.
                % 'legend plot pos' is not set explicitly, since 'left' is default.
                lStyle = M2T_opts_add(lStyle, 'legend cell align', 'left');
                lStyle = M2T_opts_add(lStyle, 'align', 'left');
            end
        end
        % ==============================================================================
        function [pTicks, pTickLabels] = ...
            M2T_matlabTicks2pgfplotsTicks(m2t, ticks, tickLabels, isLogAxis, tickLabelMode)
            % Converts MATLAB style ticks and tick labels to pgfplots style (if needed)
            if isempty(ticks)
                pTicks      = '\empty';
                pTickLabels = [];
                return
            end

            % set ticks + labels
            pTicks = M2T_join(m2t, num2cell(ticks), ',');

            % if there's no specific labels, return empty
            if isempty(tickLabels) || (length(tickLabels)==1 && isempty(tickLabels{1}))
                pTickLabels = '\empty';
                return
            end

            % sometimes tickLabels are cells, sometimes plain arrays
            % -- unify this to cells
            if ischar(tickLabels)
                tickLabels = strtrim(mat2cell(tickLabels,                  ...
                    ones(size(tickLabels,1), 1), ...
                    size(tickLabels, 2)          ...
                    ) ...
                    );
            end

            ticks = M2T_removeSuperfluousTicks(ticks, tickLabels);

            isNeeded = M2T_isTickLabelsNecessary(m2t, ticks, tickLabels, isLogAxis);

            pTickLabels = M2T_formatPgfTickLabels(m2t, isNeeded, tickLabels, ...
                isLogAxis, tickLabelMode);
        end
        % ==============================================================================
        function bool = M2T_isTickLabelsNecessary(m2t, ticks, tickLabels, isLogAxis)
            % Check if tickLabels are really necessary (and not already covered by
            % the tick values themselves).
            bool = false;

            k = find(ticks ~= 0.0, 1); % get an index with non-zero tick value
            if isLogAxis || isempty(k) % only a 0-tick
                scalingFactor = 1;
            else
                % When plotting axis, MATLAB might scale the axes by a factor of ten,
                % say 10^n, and plot a 'x 10^k' next to the respective axis. This is
                % common practice when the tick marks are really large or small
                % numbers.
                % Unfortunately, MATLAB doesn't contain the information about the
                % scaling anywhere in the plot, and at the same time the {x,y}TickLabels
                % are given as t*10^k, thus no longer corresponding to the actual
                % value t.
                % Try to find the scaling factor here. This is then used to check
                % whether or not explicit {x,y}TickLabels are really necessary.
                s = str2double(tickLabels{k});
                scalingFactor = ticks(k)/s;
                % check if the factor is indeed a power of 10
                S = log10(scalingFactor);
                if abs(round(S)-S) > m2t.tol
                    scalingFactor = 1.0;
                end
            end

            for k = 1:min(length(ticks),length(tickLabels))
                % Don't use str2num here as then, literal strings as 'pi' get
                % legally transformed into 3.14... and the need for an explicit
                % label will not be recognized. str2double returns a NaN for 'pi'.
                if isLogAxis
                    s = 10^(str2double(tickLabels{k}));
                else
                    s = str2double(tickLabels{k});
                end
                if isnan(s)  ||  abs(ticks(k)-s*scalingFactor) > m2t.tol
                    bool = true;
                    return;
                end
            end
        end
        % ==============================================================================
        function pTickLabels = M2T_formatPgfTickLabels(m2t, plotLabelsNecessary, ...
                tickLabels, isLogAxis, tickLabelMode)
            % formats the tick labels for pgfplots
            if plotLabelsNecessary
                for k = 1:length(tickLabels)
                    % Turn tickLabels from cells containing a cell into
                    % cells containing strings
                    if isnumeric(tickLabels{k})
                        tickLabels(k) = num2str(tickLabels{k});
                    elseif iscell(tickLabels{k})
                        tickLabels(k) = tickLabels{k};
                    end
                    % If the axis is logscaled, MATLAB does not store the labels,
                    % but the exponents to 10
                    if isLogAxis && strcmpi(tickLabelMode,'auto')
                        tickLabels{k} = sprintf('$10^{%s}$', str);
                    end
                end
                tickLabels = cellfun(@(l)(sprintf('{%s}',l)), tickLabels, ...
                    'UniformOutput', false);
                pTickLabels = M2T_join(m2t, tickLabels, ',');
            else
                pTickLabels = [];
            end
        end
        % ==============================================================================
        function ticks = M2T_removeSuperfluousTicks(ticks, tickLabels)
            % What MATLAB does when the number of ticks and tick labels is not the same,
            % is somewhat unclear. Cut of the first entries to fix bug
            %     https://github.com/LibAv05.M2T_matlab2tikz/LibAv05.M2T_matlab2tikz/issues/161,
            m = length(ticks);
            n = length(tickLabels);
            if n < m
                ticks = ticks(m-n+1:end);
            end
        end
        % ==============================================================================
        function tikzLineStyle = M2T_translateLineStyle(matlabLineStyle)
            if(~ischar(matlabLineStyle))
                error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_translateLineStyle:NotAString',...
                    'Variable matlabLineStyle is not a string.');
            end

            switch (matlabLineStyle)
                case 'none'
                    tikzLineStyle = '';
                case '-'
                    tikzLineStyle = 'solid';
                case '--'
                    tikzLineStyle = 'dashed';
                case ':'
                    tikzLineStyle = 'dotted';
                case '-.'
                    tikzLineStyle = 'dashdotted';
                otherwise
                    error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_translateLineStyle:UnknownLineStyle',...
                        'Unknown matlabLineStyle ''%s''.', matlabLineStyle);
            end
        end
        % ==============================================================================
        function [m2t, table, opts] = M2T_makeTable(m2t, varargin)
            %   [m2t,table,opts] = M2T_makeTable(m2t, 'name1', data1, 'name2', data2, ...)
            %   [m2t,table,opts] = M2T_makeTable(m2t, {'name1','name2',...}, {data1, data2, ...})
            %   [m2t,table,opts] = M2T_makeTable(m2t, {'name1','name2',...}, [data1(:), data2(:), ...])
            %
            %  Returns m2t structure, formatted table and table options.
            %  When all the names are empty, no header is printed
            [variables, data] = M2T_parseInputsForTable_(varargin{:});
            opts = M2T_opts_new();

            COLSEP = sprintf('\t');
            if m2t.args.externalData
                ROWSEP = sprintf('\n');
            else
                ROWSEP = sprintf('\\\\\n');
                opts = M2T_opts_add(opts, 'row sep','crcr');
            end

            nColumns = numel(data);
            nRows    = cellfun(@numel, data);
            if ~all(nRows==nRows(1))
                error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_makeTableDifferentNumberOfRows',...
                    'Different data lengths [%s].', num2str(nRows));
            end
            nRows = nRows(1);

            FORMAT = repmat({m2t.ff}, 1, nColumns);
            FORMAT(cellfun(@LibAv05.M2T_isCellOrChar, data)) = {'%s'};
            FORMAT = M2T_join(m2t, FORMAT, COLSEP);
            if all(cellfun(@isempty, variables))
                header = {};
            else
                header = {LibAv05.M2T_join(m2t, variables, COLSEP)};
            end

            table = cell(nRows,1);
            for iRow = 1:nRows
                thisData = cell(1,nColumns);
                for jCol = 1:nColumns
                    thisData{1,jCol} = data{jCol}(iRow);
                end
                table{iRow} = sprintf(FORMAT, thisData{:});
            end
            table = lower(table); % convert NaN and Inf to lower case for TikZ
            table = [LibAv05.M2T_join(m2t, [header;table], ROWSEP) ROWSEP];

            if m2t.args.externalData
                % output data to external file
                [m2t, fileNum] = M2T_incrementGlobalCounter(m2t, 'tsvFile');
                [filename, latexFilename] = M2T_externalFilename(m2t, fileNum, '.tsv');

                % write the data table to an external file
                fid = M2T_fileOpenForWrite(m2t, filename);
                finally_fclose_fid = onCleanup(@() fclose(fid));

                fprintf(fid, '%s', table);

                % put the filename in the TikZ output
                table = latexFilename;
            else
                % output data with "%newline" prepended for formatting consistency
                % do NOT prepend another newline in the output: LaTeX will crash.
                table = sprintf('%%\n%s', table);
            end
        end
        % ==============================================================================
        function [variables, data] = M2T_parseInputsForTable_(varargin)
            % parse input arguments for |LibAv05.M2T_makeTable|
            if numel(varargin) == 2 % cell syntax
                variables = varargin{1};
                data      = varargin{2};
                if ischar(variables)
                    % one variable, one data vector -> (cell, cell)
                    variables = {variables};
                    data      = {data};
                elseif iscellstr(variables) && ~iscell(data)
                    % multiple variables, one data matrix -> (cell, cell) by column
                    data = num2cell(data, 1);
                end
            else % key-value syntax
                variables = varargin(1:2:end-1);
                data      = varargin(2:2:end);
            end
        end
        % ==============================================================================
        function [path, texpath] = M2T_externalFilename(m2t, counter, extension)
            % generates a file name for an external data file and its relative TeX path

            [dummy, name] = fileparts(m2t.tikzFileName); %#ok
            baseFilename  = [name '-' num2str(counter) extension];
            path    = fullfile(m2t.dataPath, baseFilename);
            texpath = M2T_TeXpath(fullfile(m2t.relativeDataPath, baseFilename));
        end
        % ==============================================================================
        function [names,definitions] = M2T_dealColorDefinitions(mergedColorDefs)
            if isempty(mergedColorDefs)
                mergedColorDefs = {};
            end
            [names,definitions] = cellfun(@(x)(deal(x{:})),  mergedColorDefs, ...
                'UniformOutput', false);
        end
        % ==============================================================================
        function [m2t, colorLiteral] = M2T_rgb2colorliteral(m2t, rgb)
            % Translates an rgb value to an xcolor literal
            %
            % Possible outputs:
            %  - xcolor literal color, e.g. 'blue'
            %  - mixture of 2 previously defined colors, e.g. 'red!70!green'
            %  - a newly defined color, e.g. 'mycolor10'

            % Take a look at xcolor.sty for the color definitions.
            % In xcolor.sty some colors are defined in CMYK space and approximated
            % crudely for RGB color space. So it is better to redefine those colors
            % instead of using xcolor's:
            %    'cyan' , 'magenta', 'yellow', 'olive'
            %    [0,1,1], [1,0,1]  , [1,1,0] , [0.5,0.5,0]

            xcolColorNames = {'white', 'black', 'red', 'green', 'blue', ...
                              'brown', 'lime', 'orange', 'pink', ...
                              'purple', 'teal', 'violet', ...
                              'darkgray', 'gray', 'lightgray'};
            xcolColorSpecs = {[1,1,1], [0,0,0], [1,0,0], [0,1,0], [0,0,1], ...
                              [0.75,0.5,0.25], [0.75,1,0], [1,0.5,0], [1,0.75,0.75], ...
                              [0.75,0,0.25], [0,0.5,0.5], [0.5,0,0.5], ...
                              [0.25,0.25,0.25], [0.5,0.5,0.5], [0.75,0.75,0.75]};

            colorNames = [xcolColorNames, m2t.color.extraNames];
            colorSpecs = [xcolColorSpecs, m2t.color.extraSpecs];

            %% check if rgb is a predefined color
            for kColor = 1:length(colorSpecs)
                Ck = colorSpecs{kColor}(:);
                if max(abs(Ck - rgb(:))) < m2t.color.precision
                    colorLiteral = colorNames{kColor};
                    return % exact color was predefined
                end
            end

            %% check if the color is a linear combination of two already defined colors
            for iColor = 1:length(colorSpecs)
                for jColor = iColor+1:length(colorSpecs)
                    Ci = colorSpecs{iColor}(:);
                    Cj = colorSpecs{jColor}(:);

                    % solve color mixing equation `Ck = p * Ci + (1-p) * Cj` for p
                    p  = (Ci-Cj) \ (rgb(:)-Cj);
                    p  = round(100*p)/100;  % round to a percentage
                    Ck = p * Ci + (1-p)*Cj; % approximated mixed color

                    if p <= 1 && p >= 0 && max(abs(Ck(:) - rgb(:))) < m2t.color.precision
                        colorLiteral = sprintf('%s!%d!%s', colorNames{iColor}, round(p*100), ...
                            colorNames{jColor});
                        return % linear combination found
                    end
                end
            end

            %% Define colors that are not a linear combination of two known colors
            colorLiteral = sprintf('mycolor%d', length(m2t.color.extraNames)+1);
            m2t.color.extraNames{end+1} = colorLiteral;
            m2t.color.extraSpecs{end+1} = rgb;
        end
        % ==============================================================================
        function newstr = M2T_join(m2t, cellstr, delimiter)
            % This function LibAv05.M2T_joins a cell of strings to a single string (with a
            % given delimiter in between two strings, if desired).
            %
            % Example of usage:
            %              LibAv05.M2T_join(m2t, cellstr, ',')
            newstr = m2tstrLibAv05.M2T_join(cellstr, delimiter, m2t.ff);
        end
        % ==============================================================================
        function [width, height, unit] = M2T_getNaturalFigureDimension(m2t)
            % Returns the size of figure (in inch)
            % To stay compatible with getNaturalAxesDimensions, the unit 'in' is
            % also returned.

            % Get current figure size
            figuresize = get(m2t.current.gcf, 'Position');
            figuresize = figuresize([3 4]);
            figureunit = get(m2t.current.gcf, 'Units');

            % Convert Figure Size
            unit = 'in';
            figuresize = M2T_convertUnits(figuresize, figureunit, unit);

            % Split size into width and height
            width  = figuresize(1);
            height = figuresize(2);

        end
        % ==============================================================================
        function dimension = M2T_getFigureDimensions(m2t, widthString, heightString)
            % Returns the physical dimension of the figure.

            [width, height, unit] = M2T_getNaturalFigureDimension(m2t);

            % get the natural width-height ration of the plot
            axesWidthHeightRatio = width / height;
            % check LibAv05.M2T_matlab2tikz arguments
            if ~isempty(widthString)
                width = M2T_extractValueUnit(widthString);
            end
            if ~isempty(heightString)
                height = M2T_extractValueUnit(heightString);
            end

            % prepare the output
            if ~isempty(widthString) && ~isempty(heightString)
                dimension.x.unit  = width.unit;
                dimension.x.value = width.value;
                dimension.y.unit  = height.unit;
                dimension.y.value = height.value;
            elseif ~isempty(widthString)
                dimension.x.unit  = width.unit;
                dimension.x.value = width.value;
                dimension.y.unit  = width.unit;
                dimension.y.value = width.value / axesWidthHeightRatio;
            elseif ~isempty(heightString)
                dimension.y.unit  = height.unit;
                dimension.y.value = height.value;
                dimension.x.unit  = height.unit;
                dimension.x.value = height.value * axesWidthHeightRatio;
            else % neither width nor height given
                dimension.x.unit  = unit;
                dimension.x.value = width;
                dimension.y.unit  = unit;
                dimension.y.value = height;
            end
        end
        % ==============================================================================
        function position = M2T_getAxesPosition(m2t, handle, widthString, heightString, axesBoundingBox)
            % Returns the physical position of the axes. This includes - in difference
            % to the Dimension - also an offset to shift the axes inside the figure
            % An optional bounding box can be used to omit empty borders.

            % Deal with optional parameter
            if nargin < 4
                axesBoundingBox = [0 0 1 1];
            end

            % First get the whole figures size
            figDim = M2T_getFigureDimensions(m2t, widthString, heightString);

            % Get the relative position of the axis
            relPos = M2T_getRelativeAxesPosition(m2t, handle, axesBoundingBox);

            position.x.value = relPos(1) * figDim.x.value;
            position.x.unit  = figDim.x.unit;
            position.y.value = relPos(2) * figDim.y.value;
            position.y.unit  = figDim.y.unit;
            position.w.value = relPos(3) * figDim.x.value;
            position.w.unit  = figDim.x.unit;
            position.h.value = relPos(4) * figDim.y.value;
            position.h.unit  = figDim.y.unit;
        end
        % ==============================================================================
        function [position] = M2T_getRelativeAxesPosition(m2t, axesHandles, axesBoundingBox)
            % Returns the relative position of axes within the figure.
            % Position is an (n,4) matrix with [minX, minY, width, height] for each
            % handle. All these values are relative to the figure size, which means
            % that [0, 0, 1, 1] covers the whole figure.
            % It is possible to add a second parameter with the relative coordinates of
            % a bounding box around all axes of the figure (see LibAv05.M2T_getRelevantAxes()). In
            % this case, relative positions are rescaled so that the bounding box is
            % [0, 0, 1, 1]

            % Get Figure Dimension
            [figWidth, figHeight, figUnits] = M2T_getNaturalFigureDimension(m2t);

            % Initialize position
            position = zeros(numel(axesHandles), 4);
            % Iterate over all handles
            for i = 1:numel(axesHandles)
                axesHandle = axesHandles(i);
                axesPos = get(axesHandle, 'Position');
                axesUnits = get(axesHandle, 'Units');
                if isequal(lower(axesUnits), 'normalized')
                    % Position is already relative
                    position(i,:) = axesPos;
                else
                    % Convert figure size into axes units
                    figureSize = M2T_convertUnits([figWidth, figHeight], figUnits, axesUnits);
                    % Figure size into axes units to get the relative size
                    position(i,:) = axesPos ./ [figureSize, figureSize];

                end

                if strcmpi(get(axesHandle, 'DataAspectRatioMode'), 'manual') ...
                        || strcmpi(get(axesHandle, 'PlotBoxAspectRatioMode'), 'manual')

                    if strcmpi(get(axesHandle,'Projection'),'Perspective')
                        LibAv05.M2T_userWarning(m2t,'Perspective projections are not currently supported')
                    end

                    % project vertices of 3d plot box (this results in 2d coordinates in
                    % an absolute coordinate system that is scaled proportionally by
                    % Matlab to fit the axes position box)
                    switch getEnvironment()
                        case 'MATLAB'
                            projection = view(axesHandle);

                        case 'Octave'
                            % Unfortunately, Octave does not have the full `view`
                            % interface implemented, but the projection matrices are
                            % available: http://octave.1599824.n4.nabble.com/Implementing-view-td3032041.html

                            projection = get(axesHandle, 'x_viewtransform');

                        otherwise
                            errorUnknownEnvironment();
                    end


                    vertices = projection * [0, 1, 0, 0, 1, 1, 0, 1;
                        0, 0, 1, 0, 1, 0, 1, 1;
                        0, 0, 0, 1, 0, 1, 1, 1;
                        1, 1, 1, 1, 1, 1, 1, 1];

                    % each of the columns of vertices represents a vertex of the 3D axes
                    % but we only need their XY coordinates
                    verticesXY = vertices([1 2], :);

                    % the size of the projected plot box is limited by the long diagonals
                    % The matrix A determines the connectivity, e.g. the first diagonal runs from vertices(:,3) -> vertices(:,4)
                    A = [ 0,  0,  0, -1, +1,  0,  0,  0;
                        0,  0, -1,  0,  0, +1,  0,  0;
                        0, -1,  0,  0,  0,  0, +1,  0;
                        -1,  0,  0,  0,  0,  0,  0, +1];
                    diagonals = verticesXY * A';
                    % each of the columns of this matrix contains a the X and Y distance of a diagonal
                    dimensions = max(abs(diagonals), [], 2);

                    % find limiting dimension and adjust position
                    aspectRatio = dimensions(2) * figWidth / (dimensions(1) * figHeight);
                    axesAspectRatio = position(i,4) / position(i,3);
                    if aspectRatio > axesAspectRatio
                        newWidth = position(i,4) / aspectRatio;
                        % Center Axis
                        offset = (position(i,3) - newWidth) / 2;
                        position(i,1) = position(i,1) + offset;
                        % Store new width
                        position(i,3) = newWidth;
                    else
                        newHeight = position(i,3) * aspectRatio;
                        offset = (position(i,4) - newHeight) / 2;
                        position(i,2) = position(i,2) + offset;
                        % Store new height
                        position(i,4) = newHeight;
                    end
                end
            end

            %% Rescale if axesBoundingBox is given
            if exist('axesBoundingBox','var')
                % shift position so that [0, 0] is the lower left corner of the
                % bounding box
                position(:,1) = position(:,1) - axesBoundingBox(1);
                position(:,2) = position(:,2) - axesBoundingBox(2);
                % Recale
                position(:,[1 3]) = position(:,[1 3]) / max(axesBoundingBox([3 4]));
                position(:,[2 4]) = position(:,[2 4]) / max(axesBoundingBox([3 4]));
            end
        end
        % ==============================================================================
        function aspectRatio = M2T_getPlotBoxAspectRatio(axesHandle)
            limits = axis(axesHandle);
            if any(isinf(limits))
                aspectRatio = get(axesHandle,'PlotBoxAspectRatio');
            else
                % DataAspectRatio has priority
                dataAspectRatio = get(axesHandle,'DataAspectRatio');
                nlimits         = length(limits)/2;
                limits          = reshape(limits, 2, nlimits);
                aspectRatio     = abs(limits(2,:) - limits(1,:))./dataAspectRatio(1:nlimits);
                aspectRatio     = aspectRatio/min(aspectRatio);
            end
        end
        % ==============================================================================
        function texUnits = M2T_matlab2texUnits(matlabUnits, fallbackValue)
            switch matlabUnits
                case 'pixels'
                    texUnits = 'px'; % only in pdfTex/LuaTeX
                case 'centimeters'
                    texUnits = 'cm';
                case 'characters'
                    texUnits = 'em';
                case 'points'
                    texUnits = 'pt';
                case 'inches'
                    texUnits = 'in';
                otherwise
                    texUnits = fallbackValue;
            end
        end
        % ==============================================================================
        function dstValue = M2T_convertUnits(srcValue, srcUnit, dstUnit)
            % Converts values between different units.
            %   srcValue stores a length (or vector of lengths) in srcUnit.
            % The resulting dstValue is the converted length into dstUnit.
            %
            % Currently supported units are: in, cm, px, pt

            % Use tex units, if possible (to make things simple)
            srcUnit = M2T_matlab2texUnits(lower(srcUnit),lower(srcUnit));
            dstUnit = M2T_matlab2texUnits(lower(dstUnit),lower(dstUnit));

            if isequal(srcUnit, dstUnit)
                dstValue = srcValue;
                return % conversion to the same unit => factor = 1
            end

            units  = {srcUnit, dstUnit};
            factor = ones(1,2);
            for ii = 1:numel(factor) % Same code for srcUnit and dstUnit
                % Use inches as intermediate unit
                % Compute the factor to convert an inch into another unit
                switch units{ii}
                    case 'cm'
                        factor(ii) = 2.54;
                    case 'px'
                        factor(ii) = get(0, 'ScreenPixelsPerInch');
                    case 'in'
                        factor(ii) = 1;
                    case 'pt'
                        factor(ii) = 72;
                    otherwise
                        warning('MATLAB2TIKZ:UnknownPhysicalUnit',...
                        'Can not convert unit ''%s''. Using conversion factor 1.', units{ii});
                end
            end

            dstValue = srcValue * factor(2) / factor(1);
        end
        % ==============================================================================
        function out = M2T_extractValueUnit(str)
            % Decompose m2t.args.width into value and unit.

            % Regular expression to match '4.12cm', '\figurewidth', ...
            fp_regex = '[-+]?\d*\.?\d*(?:e[-+]?\d+)?';
            pattern = strcat('(', fp_regex, ')?', '(\\?[a-zA-Z]+)');

            [dummy,dummy,dummy,dummy,t,dummy] = regexp(str, pattern, 'match'); %#ok

            if length(t)~=1
                error('getAxesDimensions:illegalLength', ...
                    'The width string ''%s'' could not be decomposed into value-unit pair.', str);
            end

            if length(t{1}) == 1
                out.value = 1.0; % such as in '1.0\figurewidth'
                out.unit  = strtrim(t{1}{1});
            elseif length(t{1}) == 2 && isempty(t{1}{1})
                % MATLAB(R) does this:
                % length(t{1})==2 always, but the first field may be empty.
                out.value = 1.0;
                out.unit  = strtrim(t{1}{2});
            elseif length(t{1}) == 2
                out.value = str2double(t{1}{1});
                out.unit  = strtrim(t{1}{2});
            else
                error('getAxesDimensions:illegalLength', ...
                    'The width string ''%s'' could not be decomposed into value-unit pair.', str);
            end
        end
        % ==============================================================================
        function str = M2T_escapeCharacters(str)
            % Replaces "%" and "\" with respectively "%%" and "\\"
            str = strrep(str, '%' , '%%');
            str = strrep(str, '\' , '\\');
        end
        % ==============================================================================
        function bool = M2T_isNone(value)
            % Checks whether a value is 'none'
            bool = strcmpi(value, 'none');
        end
        % ==============================================================================
        function bool = M2T_isOn(value)
            % Checks whether a value is 'on'
            bool = strcmpi(value, 'on');
        end
        % ==============================================================================
        function bool = M2T_isOff(value)
            % Checks whether a value is 'off'.
            % Note that some options are not be solely an on/off boolean, such that `LibAv05.M2T_isOn`
            % and LibAv05.M2T_isOff don't always return the complement of each other and such that we
            % need both functions to check the value.
            % E.g. `set(0, 'HandleVisibility')` allows the value 'callback'.
            bool = strcmpi(value, 'off');
        end
        % ==============================================================================
        function val = M2T_getOrDefault(handle, key, default)
            % gets the value or returns the default value if no such property exists
            if all(isprop(handle, key))
                val = get(handle, key);
            else
                val = default;
            end
        end
        % ==============================================================================
        function val = M2T_getFactoryOrDefault(type, key, fallback)
            % get factory default value for a certain type of HG object
            % this CANNOT be done using |LibAv05.M2T_getOrDefault| as |isprop| doesn't work for
            % factory/default settings. Hence, we use a more expensive try-catch instead.
            try
                groot = 0;
                val = get(groot, ['Factory' type key]);
            catch
                val = fallback;
            end
        end
        % ==============================================================================
        function [val, isDefault] = M2T_getAndCheckDefault(type, handle, key, default)
            % gets the value from a handle of certain type and check the default values
            default   = M2T_getFactoryOrDefault(type, key, default);
            val       = M2T_getOrDefault(handle, key, default);
            isDefault = isequal(val, default);
        end
        % ==============================================================================
        function bool = M2T_isVisible(handles)
            % Determines whether an object is actually visible or not.
            bool = M2T_isOn(get(handles,'Visible'));
            % There's another handle property, 'HandleVisibility', that is unrelated
            % to the "physical" visibility of an object. Rather, it sets whether an
            % object should be visitable by |findobj|. Hence, it is often switched off
            % for non-data objects such as custom axes/grid objects.
        end
        % ==============================================================================
        function [m2t, axesBoundingBox] = M2T_getRelevantAxes(m2t, axesHandles)
            % Returns relevant axes. These are defines as visible axes that are no
            % colorbars. Function 'LibAv05.M2T_findPlotAxes()' ensures that 'axesHandles' does not
            % contain colorbars. In addition, a bounding box around all relevant Axes is
            % computed. This can be used to avoid undesired borders.
            % This function is the remaining code of alignSubPlots() in the alternative
            % positioning system.

            % List only visible axes
            N   = numel(axesHandles);
            idx = false(N,1);
            for ii = 1:N
                idx(ii) = LibAv05.M2T_LibAv05.M2T_isVisibleContainer(axesHandles(ii));
            end
            % Store the relevant axes in m2t to simplify querying e.g. positions
            % of subplots
            m2t.relevantAxesHandles = axesHandles(idx);

            % Compute the bounding box if width or height of the figure are set by
            % parameter
            if ~isempty(m2t.args.width) || ~isempty(m2t.args.height)
                % TODO: check if relevant Axes or all Axes are better.
                axesBoundingBox = M2T_getRelativeAxesPosition(m2t, m2t.relevantAxesHandles);
                % Compute second corner from width and height for each axes
                axesBoundingBox(:,[3 4]) = axesBoundingBox(:,[1 2]) + axesBoundingBox(:,[3 4]);
                % Combine axes corners to get the bounding box
                axesBoundingBox = [min(axesBoundingBox(:,[1 2]),[],1), max(axesBoundingBox(:,[3 4]), [], 1)];
                % Compute width and height of the bounding box
                axesBoundingBox(:,[3 4]) = axesBoundingBox(:,[3 4]) - axesBoundingBox(:,[1 2]);
            else
                % Otherwise take the whole figure as bounding box => lengths are
                % not changed in tikz
                axesBoundingBox = [0, 0, 1, 1];
            end
        end
        % ==============================================================================
        function M2T_userInfo(m2t, message, varargin)
            % Display usage information.
            if m2t.args.showInfo
                mess = sprintf(message, varargin{:});

                mess = strrep(mess, sprintf('\n'), sprintf('\n *** '));
                fprintf(' *** %s\n', mess);
            end
        end
        % ==============================================================================
        function M2T_userWarning(m2t, message, varargin)
            % Drop-in replacement for warning().
            if m2t.args.showWarnings
                warning('LibAv05.M2T_matlab2tikz:LibAv05.M2T_userWarning', message, varargin{:});
            end
        end
        % ==============================================================================
        function M2T_signalDependency(m2t, dependencyType, name)
            % Signals an (optional) dependency to the user
            switch lower(dependencyType)
                case 'tikzlibrary'
                    message = 'Make sure to add "\\usetikzlibrary{%s}" to the preamble.';
                otherwise
                    message = 'Please make sure to load the "%s" dependency';
            end
            LibAv05.M2T_userInfo(m2t, message, name);
        end
        % ==============================================================================
        function M2T_warnAboutParameter(m2t, parameter, isActive, message)
            % warn the user about the use of a dangerous parameter
            line = ['\n' repmat('=',1,80) '\n'];
            if isActive(m2t.args.(parameter))
                LibAv05.M2T_userWarning(m2t, [line, 'You are using the "%s" parameter.\n', ...
                                  message line], parameter);
            end
        end
        % ==============================================================================
        function parent = M2T_addChildren(parent, children)
            if isempty(children)
                return;
            elseif iscell(children)
                for k = 1:length(children)
                    parent = M2T_addChildren(parent, children{k});
                end
            else
                if isempty(parent.children)
                    parent.children = {children};
                else
                    parent.children = [parent.children children];
                end
            end
        end
        % ==============================================================================
        function M2T_printAll(m2t, env, fid)
            if isfield(env, 'colors') && ~isempty(env.colors)
                fprintf(fid, '%s', env.colors);
            end

            if isempty(env.options)
                fprintf(fid, '\\begin{%s}\n', env.name);
            else
                fprintf(fid, '\\begin{%s}[%%\n%s\n]\n', env.name, ...
                        LibAv05.M2T_opts_print(env.options, sprintf(',\n')));
            end

            for item = env.content
                fprintf(fid, '%s', char(item));
            end

            for k = 1:length(env.children)
                if ischar(env.children{k})
                    fprintf(fid, LibAv05.M2T_escapeCharacters(env.children{k}));
                else
                    fprintf(fid, '\n');
                    LibAv05.M2T_printAll(m2t, env.children{k}, fid);
                end
            end

            % End the tikzpicture environment with an empty comment and no newline
            % so no additional space is generated after the tikzpicture in TeX.
            if strcmp(env.name, 'tikzpicture') % LaTeX is case sensitive
                fprintf(fid, '\\end{%s}%%', env.name);
            else
                fprintf(fid, '\\end{%s}\n', env.name);
            end
        end
        % ==============================================================================
        function c = M2T_prettyPrint(m2t, strings, interpreter)
            % Some resources on how MATLAB handles rich (TeX) markup:
            % http://www.mathworks.com/help/techdoc/ref/text_props.html#String
            % http://www.mathworks.com/help/techdoc/creating_plots/f0-4741.html#f0-28104
            % http://www.mathworks.com/help/techdoc/ref/text_props.html#Interpreter
            % http://www.mathworks.com/help/techdoc/ref/text.html#f68-481120

            % If the user set the LibAv05.M2T_matlab2tikz parameter 'parseStrings' to false, no
            % parsing of strings takes place, thus making the user 100% responsible.
            if ~m2t.args.parseStrings
                % If strings is an actual string (labels etc) we need to return a
                % cell containing the string
                c = cellstr(strings);
                return
            end

            % Make sure we have a valid interpreter set up
            if ~any(strcmpi(interpreter, {'latex', 'tex', 'none'}))
                LibAv05.M2T_userWarning(m2t, 'Don''t know interpreter ''%s''. Default handling.', interpreter);
                interpreter = 'tex';
            end

            strings = M2T_cellstrOneLinePerCell(strings);

            % Now loop over the strings and return them pretty-printed in c.
            c = cell(1, length(strings));
            for k = 1:length(strings)
                % linear indexing for independence of cell array dimensions
                s = strings{k};

                % The interpreter property of the text element defines how the string
                % is parsed
                switch lower(interpreter)
                    case 'latex' % Basic subset of the LaTeX markup language

                        % Replace $$...$$ with $...$ for groups, but otherwise leave
                        % untouched.
                        % Displaymath \[...\] seems to be unsupported by TikZ/PGF.
                        % If this changes, use '\\[$2\\]' as replacement below.
                        % Do not escape dollar in replacement string (e.g., "\$$2\$"),
                        % since this is not properly handled by octave 3.8.2.
                        string = regexprep(s, '(\$\$)(.*?)(\$\$)', '$$2$');

                    case 'tex' % Subset of plain TeX markup language

                        % Deal with UTF8 characters.
                        string = s;

                        % degree symbol following "^" or "_" needs to be escaped
                        string = regexprep(string, '([\^\_])Â°', '$1{{}^\\circ}');
                        string = strrep(string, 'Â°', '^\circ');
                        string = strrep(string, 'âˆž', '\infty');

                        % Parse string piece-wise in a separate function.
                        string = M2T_parseTexString(m2t, string);

                    case 'none' % Literal characters
                        % Make special characters TeX compatible

                        string = strrep(s, '\', '\textbackslash{}');
                        % Note: '{' and '}' can't be converted to '\{' and '\}',
                        %       respectively, via strrep(...) as this would lead to
                        %       backslashes converted to '\textbackslash\{\}' because
                        %       the backslash was converted to '\textbackslash{}' in
                        %       the previous step. Using regular expressions with
                        %       negative look-behind makes sure any braces in 'string'
                        %       were not introduced by escaped backslashes.
                        %       Also keep in mind that escaping braces before backslashes
                        %       would not remedy the issue -- in that case 'string' would
                        %       contain backslashes introduced by brace escaping that are
                        %       not supposed to be printable characters.
                        repl = M2T_switchMatOct('\\{', '\{');
                        string = regexprep(string, '(?<!\\textbackslash){', repl);
                        repl = M2T_switchMatOct('\\}', '\}');
                        string = regexprep(string, '(?<!\\textbackslash{)}', repl);
                        string = strrep(string, '$', '\$');
                        string = strrep(string, '%', '\%');
                        string = strrep(string, '_', '\_');
                        string = strrep(string, '^', '\textasciicircum{}');
                        string = strrep(string, '#', '\#');
                        string = strrep(string, '&', '\&');
                        string = strrep(string, '~', '\textasciitilde{}'); % or '\~{}'
                        % Clean up: remove superfluous '{}' if it's followed by a backslash
                        string = strrep(string, '{}\', '\');
                        % Clean up: remove superfluous '{}' at the end of 'string'
                        string = regexprep(string, '\{\}$', '');

                        % Make sure to return a string and not a cellstr.
                        if iscellstr(string)
                            string = string{1};
                        end
                    otherwise
                        error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_prettyPrint', 'Unknown interpreter');
                end
                c{k} = string;
            end
        end
        % ==============================================================================
        function strings = M2T_cellstrOneLinePerCell(strings)
            % convert to cellstr that contains only one-line strings
            if ischar(strings)
                strings = cellstr(strings);
            elseif iscellstr(strings)
                cs = cell(1, length(strings));
                for s = 1:length(strings)
                    tmp = cellstr(strings{s});
                    cs{s} = tmp;
                end
                strings = cs;
            else
                error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_cellstrOneLinePerCell', ...
                    'Data type not understood.');
            end
        end
        % ==============================================================================
        function parsed = M2T_parseTexString(m2t, string)
            if iscellstr(string)
                % Convert cell string to regular string, otherwise MATLAB complains
                string = string{:};
            end

            % Get the position of all braces
            bracesPos = regexp(string, '\{|\}');

            % Exclude braces that are part of any of these MATLAB-supported TeX commands:
            % \color{...}  \color[...]{...}  \fontname{...}  \fontsize{...}
            [sCmd, eCmd] = regexp(string, '\\(color(\[[^\]]*\])?|fontname|fontsize)\{[^}]*\}');
            for i = 1:length(sCmd)
                bracesPos(bracesPos >= sCmd(i) & bracesPos <= eCmd(i)) = [];
            end

            % Exclude braces that are preceded by an odd number of backslashes which
            % means the brace is escaped and thus to be printed, not a grouping brace
            expr = '(?<!\\)(\\\\)*\\(\{|\})';
            escaped = regexp(string, expr, 'end');
            % It's necessary to go over 'string' with the same RegEx again to catch
            % overlapping matches, e.g. string == '\{\}'. In such a case the simple
            % regexp(...) above only finds the first brace. What we have to do is look
            % only at the part of 'string' that starts with the first brace but doesn't
            % encompass its escaping backslash. Iterating over all previously found
            % matches makes sure all overlapping matches are found, too. That way even
            % cases like string == '\{\} \{\}' are handled correctly.
            % The call to unique(...) is not necessary to get the behavior described, but
            % by removing duplicates in 'escaped' it's cleaner than without.
            for i = escaped
                escaped = unique([escaped, regexp(string(i:end), expr, 'end') + i-1]);
            end
            % Now do the actual removal of escaped braces
            for i = 1:length(escaped)
                bracesPos(bracesPos == escaped(i)) = [];
            end

            parsed = '';
            % Have a virtual brace one character left of where the actual string
            % begins (remember, MATLAB strings start counting at 1, not 0). This is
            % to make sure substrings left of the first brace get parsed, too.
            prevBracePos = 0;
            % Iterate over all the brace positions in order to split up 'string'
            % at those positions and then parse the substrings. A virtual brace is
            % added right of where the actual string ends to make sure substrings
            % right of the right-most brace get parsed as well.
            for currBracePos = [bracesPos, length(string)+1]
                if (prevBracePos + 1) < currBracePos
                    % Parse the substring between (but not including) prevBracePos
                    % and currBracePos, i.e. between the previous brace and the
                    % current one (but only if there actually is a non-empty
                    % substring). Then append it to the output string.
                    substring = string(prevBracePos+1 : currBracePos-1);
                    parsed = [parsed, LibAv05.M2T_parseTexSubstring(m2t, substring)];
                end
                if currBracePos <= length(string)
                    % Append the brace itself to the output string, but only if the
                    % current brace position is within the limits of the string, i.e.
                    % don't append anything for the last, virtual brace that is only
                    % there to enable parsing of substrings beyond the right-most
                    % actual brace.
                    brace = string(currBracePos);
                    parsed = [parsed, brace];
                end
                % The current brace position will be next iteration's previous one
                prevBracePos = currBracePos;
            end

            % Enclose everything in $...$ to use math mode
            parsed = ['$' parsed '$'];
            % ...except when everything is text
            parsed = regexprep(parsed, '^\$\\text\{([^}]*)\}\$$', '$1');
            % start-> $ \text {(non-}) } $<-end
            % ...or when the parsed string is empty
            parsed = regexprep(parsed, '^\$\$$', '');

            % Ensure math mode for pipe symbol (issue #587)
            parsed = strrep(parsed, '|', '$|$');
        end
        % ==============================================================================
        function string = M2T_parseTexSubstring(m2t, string)
            origstr = string; % keep this for warning messages

            % Font families (italic, bold, etc.) get a trailing '{}' because they may be
            % followed by a letter which would produce an error in (La)TeX.
            for i = {'it', 'bf', 'rm', 'sl'}
                string = strrep(string, ['\' i{:}], ['\' i{:} '{}']);
            end

            % The same holds true for special characters like \alpha
            % The list of MATLAB-supported TeX characters was taken from
            % http://www.mathworks.com/help/techdoc/ref/text_props.html#String
            named = {'alpha', 'angle', 'ast', 'beta', 'gamma', 'delta',     ...
                'epsilon', 'zeta', 'eta', 'theta', 'vartheta', 'iota', ...
                'kappa', 'lambda', 'mu', 'nu', 'xi', 'pi', 'rho',      ...
                'sigma', 'varsigma', 'tau', 'equiv', 'Im', 'otimes',   ...
                'cap', 'int', 'rfloor', 'lfloor', 'perp', 'wedge',     ...
                'rceil', 'vee', 'langle', 'upsilon', 'phi', 'chi',     ...
                'psi', 'omega', 'Gamma', 'Delta', 'Theta', 'Lambda',   ...
                'Xi', 'Pi', 'Sigma', 'Upsilon', 'Phi', 'Psi', 'Omega', ...
                'forall', 'exists', 'ni', 'cong', 'approx', 'Re',      ...
                'oplus', 'cup', 'subseteq', 'lceil', 'cdot', 'neg',    ...
                'times', 'surd', 'varpi', 'rangle', 'sim', 'leq',      ...
                'infty', 'clubsuit', 'diamondsuit', 'heartsuit',       ...
                'spadesuit', 'leftrightarrow', 'leftarrow',            ...
                'Leftarrow', 'uparrow', 'rightarrow', 'Rightarrow',    ...
                'downarrow', 'circ', 'pm', 'geq', 'propto', 'partial', ...
                'bullet', 'div', 'neq', 'aleph', 'wp', 'oslash',       ...
                'supseteq', 'nabla', 'ldots', 'prime', '0', 'mid',     ...
                'copyright'                                            };
            for i = named
                string = strrep(string, ['\' i{:}], ['\' i{:} '{}']);
                % FIXME: Only append '{}' if there's an odd number of backslashes
                %        in front of the items from 'named'. If it's an even
                %        number instead, that means there's an escaped (printable)
                %        backslash and some text like "alpha" after that.
            end
            % Some special characters' names are subsets of others, e.g. '\o' is
            % a subset of '\omega'. This would produce undesired double-escapes.
            % For example if '\o' was converted to '\o{}' after '\omega' has been
            % converted to '\omega{}' this would result in '\o{}mega{}' instead of
            % '\omega{}'. Had '\o' been converted to '\o{}' _before_ '\omega' is
            % converted then the result would be '\o{}mega' and thus also wrong.
            % To circumvent the problem all those special character names that are
            % subsets of others are now converted using a regular expression that
            % uses negative lookahead. The special handling of the backslash is
            % required for MATLAB/Octave compatibility.
            string = regexprep(string, '(\\)o(?!mega|times|plus|slash)', '$1o{}');
            string = regexprep(string, '(\\)in(?!t|fty)', '$1in{}');
            string = regexprep(string, '(\\)subset(?!eq)', '$1subset{}');
            string = regexprep(string, '(\\)supset(?!eq)', '$1supset{}');

            % Convert '\0{}' (TeX text mode) to '\emptyset{}' (TeX math mode)
            string = strrep(string, '\0{}', '\emptyset{}');

            % Add skip to \fontsize
            % This is required for a successful LaTeX run on the output as in contrast
            % to MATLAB/Octave it requires the skip parameter (even if it's zero)
            string = regexprep(string, '(\\fontsize\{[^}]*\})', '$1{0}');

            % Put '\o{}' inside \text{...} as it is a text mode symbol that does not
            % exist in math mode (and LaTeX gives a warning if you use it in math mode)
            string = strrep(string, '\o{}', '\text{\o{}}');

            % Put everything that isn't a TeX command inside \text{...}
            expr = '(\\[a-zA-Z]+(\[[^\]]*\])?(\{[^}]*\}){1,2})';
            % |( \cmd  )( [...]?  )( {...}{1,2} )|
            % (              subset $1               )
            repl = '}$1\\text{';
            string = regexprep(string, expr, repl);
            % ...\alpha{}... -> ...}\alpha{}\text{...
            string = ['\text{' string '}'];
            % ...}\alpha{}\text{... -> \text{...}\alpha{}\text{...}

            % '_' has to be in math mode so long as it's not escaped as '\_' in which
            % case it remains as-is. Extra care has to be taken to make sure any
            % backslashes in front of the underscore are not themselves escaped and
            % thus printable backslashes. This is the case if there's an even number
            % of backslashes in a row.
            repl = '$1}_\\text{';
            string = regexprep(string, '(?<!\\)((\\\\)*)_', repl);

            % '^' has to be in math mode so long as it's not escaped as '\^' in which
            % case it is expressed as '\textasciicircum{}' for compatibility with
            % regular TeX. Same thing here regarding even/odd number of backslashes
            % as in the case of underscores above.
            repl = '$1\\textasciicircum{}';
            string = regexprep(string, '(?<!\\)((\\\\)*)\\\^', repl);
            repl = '$1}^\\text{';
            string = regexprep(string, '(?<!\\)((\\\\)*)\^', repl);

            % '<' and '>' has to be either in math mode or needs to be typeset as
            % '\textless' and '\textgreater' in textmode
            % This is handled better, if 'LibAv05.M2T_parseStringsAsMath' is activated
            if m2t.args.LibAv05.M2T_parseStringsAsMath == 0
                string = regexprep(string, '<', '\\textless{}');
                string = regexprep(string, '>', '\\textgreater{}');
            end

            % Move font styles like \bf into the \text{} command.
            expr = '(\\bf|\\it|\\rm|\\fontname)({\w*})+(\\text{)';
            while regexp(string, expr)
                string = regexprep(string, expr, '$3$1$2');
            end

            % Replace Fontnames
            [dummy, dummy, dummy, dummy, fonts, dummy, subStrings] = regexp(string, '\\fontname{(\w*)}'); %#ok
            fonts = M2T_fonts2tex(fonts);
            subStrings = [subStrings; fonts, {''}];
            string = cell2mat(subStrings(:)');

            % Merge adjacent \text fields:
            string = M2T_mergeAdjacentTexCmds(string, '\text');

            % '\\' has to be escaped to '\textbackslash{}'
            % This cannot be done with strrep(...) as it would replace e.g. 4 backslashes
            % with three times the replacement string because it finds overlapping matches
            % (see http://www.mathworks.de/help/techdoc/ref/strrep.html)
            % Note: Octave's backslash handling is broken. Even though its output does
            % not resemble MATLAB's, the same m2t code is used for either software. That
            % way MATLAB-compatible code produces the same LibAv05.M2T_matlab2tikz output no matter
            % which software it's executed in. So long as this MATLAB incompatibility
            % remains in Octave you're probably better off not using backslashes in TeX
            % text anyway.
            string = regexprep(string, '(\\)\\', '$1textbackslash{}');

            % '_', '^', '{', and '}' are already escaped properly, even in MATLAB's TeX
            % dialect (and if they're not, that's intentional)

            % Escape "$", "%", and "#" to make them compatible to true TeX while in
            % MATLAB/Octave they are not escaped
            string = strrep(string, '$', '\$');
            string = strrep(string, '%', '\%');
            string = strrep(string, '#', '\#');

            % Escape "Â§" as "\S" since it can give UTF-8 problems otherwise.
            % The TeX string 'a_Â§' in particular lead to problems in Octave 3.6.0.
            % m2t transcoded that string into '$\text{a}_\text{*}\text{#}$' with
            % * = 0xC2 and # = 0xA7 which corresponds with the two-byte UTF-8
            % encoding. Even though this looks like an Octave bug that shows
            % during the '..._\text{abc}' to '..._\text{a}\text{bc}' conversion,
            % it's best to include the workaround here.
            string = strrep(string, 'Â§', '\S{}');

            string = M2T_escapeAmpersands(m2t, string, origstr);
            string = M2T_escapeTildes(m2t, string, origstr);

            % Convert '..._\text{abc}' and '...^\text{abc}' to '..._\text{a}\text{bc}'
            % and '...^\text{a}\text{bc}', respectively.
            % Things get a little more complicated if instead of 'a' it's e.g. '$'. The
            % latter has been converted to '\$' by now and simply extracting the first
            % character from '\text{\$bc}' would result in '\text{$}\text{$bc}' which
            % is syntactically wrong. Instead the whole command '\$' has to be moved in
            % front of the \text{...} block, e.g. '..._\text{\$bc}' -> '..._\$\text{bc}'.
            % Note that the problem does not occur for the majority of special characters
            % like '\alpha' because they use math mode and therefore are never inside a
            % \text{...} block to begin with. This means that the number of special
            % characters affected by this issue is actually quite small:
            %   $ # % & _ { } \o Â§ ~ \ ^
            expr = ['(_|\^)(\\text)\{([^}\\]|\\\$|\\#|\\%|\\&|\\_|\\\{|\\\}|', ...
                ... %   (_/^)(\text) {(non-}\| \$ | \#| \%| \&| \_| \{ | \} |
                ... %   ($1)( $2 )  (                 $3                      ->
                '\\o\{\}|\\S\{\}|\\textasciitilde\{\}|\\textbackslash\{\}|', ...
                ... %    \o{}  | \S{}  | \textasciitilde{}  | \textbackslash{}  |
                ... %  <-                         $3                                 ->
                '\\textasciicircum\{\})'];
            %    \textasciicircum{} )
            %  <-      $3           )
            string = regexprep(string, expr, '$1$2{$3}$2{');

            string = M2T_parseStringsAsMath(m2t, string);

            % Clean up: remove empty \text{}
            string = strrep(string, '\text{}', '');
            % \text{}\alpha{}\text{...} -> \alpha{}\text{...}

            % Clean up: convert '{}\' to '\' unless it's prefixed by a backslash which
            % means the opening brace is escaped and thus a printable character instead
            % of a grouping brace.
            string = regexprep(string, '(?<!\\)\{\}(\\)', '$1');
            % \alpha{}\text{...} -> \alpha\text{...}

            % Clean up: convert '{}}' to '}' unless it's prefixed by a backslash
            string = regexprep(string, '(?<!\\)\{\}\}', '}');

            % Clean up: remove '{}' at the end of 'string' unless it's prefixed by a
            % backslash
            string = regexprep(string, '(?<!\\)\{\}$', '');
        end
        % ==============================================================================
        function string = M2T_escapeTildes(m2t, string, origstr)
            % Escape plain "~" in MATLAB and replace escaped "\~" in Octave with a proper
            % escape sequence. An un-escaped "~" produces weird output in Octave, thus
            % give a warning in that case
            switch getEnvironment
                case 'MATLAB'
                    string = strrep(string, '~', '\textasciitilde{}'); % or '\~{}'
                case 'Octave'
                    string = strrep(string, '\~', '\textasciitilde{}'); % ditto
                    if regexp(string, '(?<!\\)~')
                        LibAv05.M2T_userWarning(m2t,                                     ...
                            ['TeX string ''%s'' contains un-escaped ''~''. ' ...
                            'For proper display in Octave you probably '     ...
                            'want to escape it even though that''s '         ...
                            'incompatible with MATLAB. '                     ...
                            'In the LibAv05.M2T_matlab2tikz output it will have its '    ...
                            'usual TeX function as a non-breaking space.'],  ...
                            origstr)
                    end
                otherwise
                    errorUnknownEnvironment();
            end
        end
        % ==============================================================================
        function string = M2T_escapeAmpersands(m2t, string, origstr)
            % Escape plain "&" in MATLAB and replace it and the following character with
            % a space in Octave unless the "&" is already escaped
            switch getEnvironment
                case 'MATLAB'
                    string = strrep(string, '&', '\&');
                case 'Octave'
                    % Ampersands should already be escaped in Octave.
                    % Octave (LibAv05.M2T_tested with 3.6.0) handles un-escaped ampersands a little
                    % funny in that it removes the following character, if there is one:
                    % 'abc&def'      -> 'abc ef'
                    % 'abc&\deltaef' -> 'abc ef'
                    % 'abc&$ef'      -> 'abc ef'
                    % 'abcdef&'      -> 'abcdef'
                    % Don't remove closing brace after '&' as this would result in
                    % unbalanced braces
                    string = regexprep(string, '(?<!\\)&(?!})', ' ');
                    string = regexprep(string, '(?<!\\)&}', '}');
                    if regexp(string, '(?<!\\)&\\')
                        % If there's a backslash after the ampersand, that means not only
                        % the backslash should be removed but the whole escape sequence,
                        % e.g. '\delta' or '\$'. Actually the '\delta' case is the
                        % trickier one since by now 'string' would have been turned from
                        % 'abc&\deltaef' into '\text{abc&}\delta{}\text{ef}', i.e. after
                        % the ampersand first comes a closing brace and then '\delta';
                        % the latter as well as the ampersand itself should be removed
                        % while the brace must remain in place to avoid unbalanced braces.
                        LibAv05.M2T_userWarning(m2t,                                                ...
                            ['TeX string ''%s'' contains a special character '  ...
                            'after an un-escaped ''&''. The output generated ' ...
                            'by LibAv05.M2T_matlab2tikz will not precisely match that '    ...
                            'which you see in Octave itself in that the '      ...
                            'special character and the preceding ''&'' is '    ...
                            'not replaced with a space.'], origstr)
                    end
                otherwise
                    errorUnknownEnvironment();
            end
        end
        % ==============================================================================
        function [string] = M2T_parseStringsAsMath(m2t, string)
            % Some further processing makes the output behave more like TeX math mode,
            % but only if the LibAv05.M2T_matlab2tikz parameter LibAv05.M2T_parseStringsAsMath=true.
            if m2t.args.LibAv05.M2T_parseStringsAsMath

                % Some characters should be in math mode: =-+/,.()<>0-9
                expr = '(\\text)\{([^}=\-+/,.()<>0-9]*)([=\-+/,.()<>0-9]+)([^}]*)\}';
                %    \text  {(any non-"x"/'}'char)( any "x" char  )(non-}) }
                %  ( $1 )  (       $2        )(      $3       )( $4)
                while regexp(string, expr)
                    % Iterating is necessary to catch all occurrences. See above.
                    string = regexprep(string, expr, '$1{$2}$3$1{$4}');
                end

                % \text{ } should be a math-mode space
                string = regexprep(string, '\\text\{(\s+)}', '$1');

                % '<<' probably means 'much smaller than', i.e. '\ll'
                repl = M2T_switchMatOct('$1\\ll{}$2', '$1\ll{}$2');
                string = regexprep(string, '([^<])<<([^<])', repl);

                % '>>' probably means 'much greater than', i.e. '\gg'
                repl = M2T_switchMatOct('$1\\gg{}$2', '$1\gg{}$2');
                string = regexprep(string, '([^>])>>([^>])', repl);

                % Single letters are most likely variables and thus should be in math mode
                string = regexprep(string, '\\text\{([a-zA-Z])\}', '$1');

            end
        end
        % ==============================================================================
        function tex = M2T_fonts2tex(fonts)
            % Returns a tex command for each fontname in the cell array fonts.
            if ~iscell(fonts)
                error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_fonts2tex', ...
                         'Expecting a cell array as input.');
            end
            tex = cell(size(fonts));

            for ii = 1:numel(fonts)
                font = fonts{ii}{1};

                % List of known fonts.
                switch lower(font)
                    case 'courier'
                        tex{ii} = '\ttfamily{}';
                    case 'times'
                        tex{ii} = '\rmfamily{}';
                    case {'arial', 'helvetica'}
                        tex{ii} = '\sffamily{}';
                    otherwise
                        warning('LibAv05.M2T_matlab2tikz:LibAv05.M2T_fonts2tex', ...
                            'Unknown font ''%s''. Using tex default font.',font);
                        % Unknown font -> Switch to standard font.
                        tex{ii} = '\rm{}';
                end
            end
        end
        % ==============================================================================
        function string = M2T_mergeAdjacentTexCmds(string, cmd)
            % Merges adjacent tex commands like \text into one command
            % If necessary, add a backslash
            if cmd(1) ~= '\'
                cmd = ['\' cmd];
            end
            % Link each bracket to the corresponding bracket
            link = zeros(size(string));
            pos = [regexp([' ' string], '([^\\]{)'), ...
                regexp([' ' string], '([^\\]})')];
            pos = sort(pos);
            ii = 1;
            while ii <= numel(pos)
                if string(pos(ii)) == '}'
                    link(pos(ii-1)) = pos(ii);
                    link(pos(ii)) = pos(ii - 1);
                    pos([ii-1, ii]) = [];
                    ii = ii - 1;
                else
                    ii = ii + 1;
                end
            end
            % Find dispensable commands
            pos = regexp(string, ['}\' cmd '{']);
            delete = zeros(0,1);
            len = numel(cmd);
            for p = pos
                l = link(p);
                if l > len && isequal(string(l-len:l-1), cmd)
                    delete(end+1,1) = p;
                end
            end
            %   3. Remove these commands (starting from the back
            delete = repmat(delete, 1, len+2) + repmat(0:len+1,numel(delete), 1);
            string(delete(:)) = [];
        end
        function dims = M2T_pos2dims(pos)
            % Position quadruplet [left, bottom, width, height] to dimension structure
            dims = struct('left' , pos(1), 'bottom', pos(2));
            if numel(pos) == 4
                dims.width  = pos(3);
                dims.height = pos(4);
                dims.right  = dims.left   + dims.width;
                dims.top    = dims.bottom + dims.height;
            end
        end
        % OPTION ARRAYS ================================================================
        function opts = M2T_opts_new()
            % create a new options array
            opts = cell(0,2);
        end
        function opts = M2T_opts_add(opts, key, value)
            % add a key-value pair to an options array (with duplication check)
            if ~exist('value','var')
                value = '';
            end
            value = char(value);

            % Check if the key already exists.
            if LibAv05.M2T_opts_has(opts, key)
                oldValue = M2T_opts_get(opts, key);
                if isequal(value, oldValue)
                    return; % no action needed: value already present
                else
                    error('LibAv05.M2T_matlab2tikz:LibAv05.M2T_opts_add', ...
                         ['Trying to add (%s, %s) to options, but it already ' ...
                          'contains the conflicting key-value pair (%s, %s).'], ...
                          key, value, key, oldValue);
                end
            end
            opts = M2T_opts_append(opts, key, value);
        end
        function opts = M2T_opts_addSubOpts(opts, key, subOpts)
            % add a key={Opts} pair to an options array
            formatted = ['{' LibAv05.M2T_opts_print(subOpts) '}'];
            opts      = M2T_opts_add(opts, key, formatted);
        end
        function bool = M2T_opts_has(opts, key)
            % returns true if the options array contains the key
            bool = ~isempty(opts) && ismember(key, opts(:,1));
        end
        function value = M2T_opts_get(opts, key)
            % returns the value(s) stored for a key in an options array
            idx = find(ismember(opts(:,1), key));
            switch numel(idx)
                case 1
                    value = opts{idx,2}; % just the value
                otherwise
                    value = opts(idx,2); % as cell array
            end
        end
        function opts = M2T_opts_append(opts, key, value)
            % append a key-value pair to an options array (duplicate keys allowed)
            if ~exist('value','var')
                value = '';
            end
            value = char(value);
            if ~(LibAv05.M2T_opts_has(opts, key) && isequal(LibAv05.M2T_opts_get(opts, key), value))
                opts = cat(1, opts, {key, value});
            end
        end
        function opts = M2T_opts_append_userdefined(opts, userDefined)
            % appends user-defined options to an options array
            % the userDefined options can come either as a single string or a cellstr that
            % is already TikZ-formatted. The internal 2D cell format is NOT supported.
            if ~isempty(userDefined)
                if ischar(userDefined)
                    userDefined = {userDefined};
                end
                for k = 1:length(userDefined)
                    opts = M2T_opts_append(opts, userDefined{k});
                end
            end
        end
        function opts = M2T_opts_copy(opts_from, name_from, opts, name_to)
            % copies an option (if it exists) from one option array to another one
            if ~exist('name_to', 'var') || isempty(name_to)
                name_to = name_from;
            end
            if LibAv05.M2T_opts_has(opts_from, name_from)
                value = M2T_opts_get(opts_from, name_from);
                opts = M2T_opts_append(opts, name_to, value);
            end
        end
        function opts = M2T_opts_remove(opts, varargin)
            % remove some key-value pairs from an options array
            keysToDelete = varargin;
            idxToDelete = ismember(opts(:,1), keysToDelete);
            opts(idxToDelete, :) = [];
        end
        function opts = M2T_opts_merge(opts, varargin)
            % merge multiple options arrays
            for jArg = 1:numel(varargin)
                opts2 = varargin{jArg};
                for k = 1:size(opts2, 1)
                    opts = M2T_opts_append(opts, opts2{k,1}, opts2{k,2});
                end
            end
        end
        function str = M2T_opts_print(opts, sep)
            % pretty print an options array
            if ~exist('sep','var') || ~ischar(sep)
                sep = ', ';
            end
            nOpts = size(opts,1);
            c = cell(1,nOpts);
            for k = 1:nOpts
                if isempty(opts{k,2})
                    c{k} = sprintf('%s', opts{k,1});
                else
                    c{k} = sprintf('%s=%s', opts{k,1}, opts{k,2});
                end
            end
            str = m2tstrLibAv05.M2T_join(c, sep);
        end
        % ==============================================================================
        function m2t = M2T_m2t_addAxisOption(m2t, key, value)
            % Adds an option to the last axesContainer
            if ~exist('value','var')
                value = '';
            end
            m2t.axes{end}.options = M2T_opts_add(m2t.axes{end}.options, key, value);
        end
        % ==============================================================================
        function bool = M2T_isHG2()
            % Checks if graphics system is HG2 (true) or HG1 (false).
            % HG1 : MATLAB up to R2014a and currently all OCTAVE versions
            % HG2 : MATLAB starting from R2014b (version 8.4)
            [env, envVersion] = getEnvironment();
            bool = strcmpi(env,'MATLAB') && ~isVersionBelow(envVersion, [8,4]);
        end
        % ==============================================================================
        function str = M2T_formatAspectRatio(m2t, values)
            % format the aspect ratio. Behind the scenes, LibAv05.M2T_formatDim is used
            strs = arrayfun(@LibAv05.M2T_formatDim, values, 'UniformOutput', false);
            str = M2T_join(m2t, strs, ' ');
        end
        % ==============================================================================
        function str = M2T_formatDim(value, unit)
            % format the value for use as a TeX dimension
            if ~exist('unit','var') || isempty(unit)
                unit = '';
            end
            tolerance = 1e-7;
            value  = round(value/tolerance)*tolerance;
            if value == 1 && ~isempty(unit) && unit(1) == '\'
                str = unit; % just use the unit
            else
                % LaTeX has support for single precision (about 6.5 decimal places),
                % but such accuracy is overkill for positioning. We clip to three
                % decimals to overcome numerical rounding issues that tend to be very
                % platform and version dependent. See also #604.
                str = sprintf('%.3f', value);
                str = regexprep(str, '(\d*\.\d*?)0+$', '$1'); % remove trailing zeros
                str = regexprep(str, '\.$', ''); % remove trailing period
                str = [str unit];
            end
        end
        % ==============================================================================
        function [retval] = M2T_switchMatOct(matlabValue, octaveValue)
            % Returns a different value for MATLAB and Octave
            switch getEnvironment
                case 'MATLAB'
                    retval = matlabValue;
                case 'Octave'
                    retval = octaveValue;
                otherwise
                    errorUnknownEnvironment();
            end
        end
        % ==============================================================================
        function M2T_checkDeprecatedEnvironment(minimalVersions)
            [env, envVersion] = getEnvironment();
            if isfield(minimalVersions, env)
                minVersion = minimalVersions.(env);
                envWithVersion = sprintf('%s %s', env, minVersion.name);

                if isVersionBelow(envVersion, minVersion.num)
                    ID = 'LibAv05.M2T_matlab2tikz:deprecatedEnvironment';

                    warningMessage = ['\n', repmat('=',1,80), '\n\n', ...
                        '  LibAv05.M2T_matlab2tikz is LibAv05.M2T_tested and developed on   %s   and newer.\n', ...
                        '  This script may still be able to handle your plots, but if you\n', ...
                        '  hit a bug, please consider upgrading your environment first.\n', ...
                        '  Type "warning off %s" to suppress this warning.\n', ...
                        '\n', repmat('=',1,80), ];
                    warning(ID, warningMessage, envWithVersion, ID);

                end
            else
                errorUnknownEnvironment();
            end
        end
        % ==============================================================================
        function m2t = M2T_needsPgfplotsVersion(m2t, minVersion)
            if isVersionBelow(m2t.pgfplotsVersion, minVersion)
                m2t.pgfplotsVersion = minVersion;
            end
        end
        % ==============================================================================
        function str = M2T_formatPgfplotsVersion(version)
            version = versionArray(version);
            if all(isfinite(version))
                str = sprintf('%d.',version);
                str = str(1:end-1); % remove the last period
            else
                str = 'newest';
            end
        end
        % ==============================================================================
        function [formatted,treeish] = M2T_VersionControlIdentifier()
            % This function gives the (git) commit ID of LibAv05.M2T_matlab2tikz
            %
            % This assumes the standard directory structure as used by Nico's master branch:
            %     SOMEPATH/src/LibAv05.M2T_matlab2tikz.m with a .git directory in SOMEPATH.
            %
            % The HEAD of that repository is determined from file system information only
            % by following dynamic references (e.g. ref:refs/heds/master) in branch files
            % until an absolute commit hash (e.g. 1a3c9d1...) is found.
            % NOTE: Packed branch references are NOT supported by this approach
            MAXITER     = 10; % stop following dynamic references after a while
            formatted   = '';
            REFPREFIX   = 'ref:';
            isReference = @(treeish)(any(strfind(treeish, REFPREFIX)));
            treeish     = [REFPREFIX 'HEAD'];
            try
                % get the LibAv05.M2T_matlab2tikz directory
                m2tDir = fileparts(mfilename('fullpath'));
                gitDir = fullfile(m2tDir,'..','.git');

                nIter = 1;
                while isReference(treeish)
                    refName    = treeish(numel(REFPREFIX)+1:end);
                    branchFile = fullfile(gitDir, refName);

                    if exist(branchFile, 'file') && nIter < MAXITER
                        % The FID is reused in every iteration, so `onCleanup` cannot
                        % be used to `fclose(fid)`. But since there is very little that
                        % can go wrong in a single `fscanf`, it's probably best to leave
                        % this part as it is for the time being.
                        fid     = fopen(branchFile,'r');
                        treeish = fscanf(fid,'%s');
                        fclose(fid);
                        nIter   = nIter + 1;
                    else % no branch file or iteration limit reached
                        treeish = '';
                        return;
                    end
                end
            catch
                treeish = '';
            end
            if ~isempty(treeish)
                formatted = sprintf('(commit %s)',treeish);
            end
        end
        % ==============================================================================




        
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %-------------------------------------------------------------
        % Rutinas de Fiteo
        %-------------------------------------------------------------
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function Peakfit_______________________l() 
        end  % end function   
        %-------------------------------------------------------------
        % Rutinas peakfit
        %-------------------------------------------------------------        
        function [FitResults,GOF,baseline,coeff,residual,xi,yi,BootResults]=peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA,clipheight)
        % A command-line peak fitting program for time-series signals, written as a
        % self-contained Matlab function in a single m-file. Uses a non-linear
        % optimization algorithm to decompose a complex, overlapping-peak signal
        % into its component parts. The objective is to determine whether your
        % signal can be represented as the sum of fundamental underlying peaks
        % shapes. Accepts signals of any length, including those with non-integer
        % and non-uniform x-values. Fits any number of peaks of any of 45 curve
        % shapes. This is a command line version, usable from a remote terminal. It
        % is capable of making multiple trial fits with sightly different starting
        % values and taking the one with the lowest mean fit error (example 6). It
        % can estimate the standard deviation of peak parameters from a single
        % signal using the bootstrap method (example 10).
        %
        % If you are unsure what input arguments to use (number of peaks, shape,
        % etc.), try fitting your signal using the interactive peak fitter ipf.m
        % (Matlab only), which uses single keystrokes to pan and zoom, select
        % number of peaks, peak shape, baseline mode, etc. Then press W to print
        % out the peakfit command line for that fit.
        %
        % Important: the data matrix "signal" must be a 2xn or nx2 matrix; the x
        % and y vectors must be separate rows or columns and not concatenated into
        % one long vector.
        %
        % Version 8.1: July 2016, adds  shape 46 (quadslope), for fitting
        % curved baselines. See example 33.
        %
        % For more details, see
        % http://terpconnect.umd.edu/~toh/spectrum/CurveFittingC.html and
        % http://terpconnect.umd.edu/~toh/spectrum/InteractivePeakFitter.htm
        %
        % peakfit(signal);       
        % Performs an iterative least-squares fit of a single Gaussian peak to the
        % data matrix "signal", which has x values in column 1 and Y values in
        % column 2 (e.g. [x y]). 
        %
        % peakfit(signal,center,window);
        % Fits a single Gaussian peak to a portion of the matrix "signal". The
        % portion is centered on the x-value "center" and has width "window" (in x
        % units).
        % 
        % peakfit(signal,center,window,NumPeaks..);
        % "NumPeaks" = number of peaks in the model (default is 1 if not
        % specified). No limit to maximum number of peaks in version 3.1
        % 
        % peakfit(signal,center,window,NumPeaks,peakshape); 
        % "peakshape" specifies the peak shape of the model: (1=Gaussian (default),
        % 2=Lorentzian, 3=logistic distribution, 4=Pearson, 5=exponentionally
        % broadened Gaussian; 6=equal-width Gaussians; 7=Equal-width Lorentzians;
        % 8=exponentionally broadened equal-width Gaussian, 9=exponential pulse,
        % 10=up-sigmoid (logistic function), 11=Fixed-width Gaussian, 12=Fixed-
        % width Lorentzian; 13=Gaussian/Lorentzian blend; 14=Bifurcated Gaussian,
        % 15=Breit-Wigner-Fano, 16=Fixed-position Gaussians; 17=Fixed-position 
        % Lorentzians; 18=exponentionally broadened Lorentzian; 19=alpha function; 
        % 20=Voigt profile; 21=triangular; 22=multiple shapes; 23=down-sigmoid; 
        % 25=lognormal; 26=slope; 27=Gaussian first derivative; 28=polynomial; 
        % 29=piecewise linear; 30=variable-alpha Voigt; 31=variable time constant 
        % ExpGaussian; 32=variable shape Pearson; 33=variable Gaussian/Lorentzian
        % blend, 34=fixed-width Voigt; 35=fixed-width Gaussian/Lorentzian blend; 
        % 36=fixed-width exponentionally-broadened Gaussian; 37=fixed-width 
        % Pearson; 38=variable time constant ExpLorentzian; 40=sine wave;
        % 41=rectangle; 42=flattened Gaussian; 43:Gompertz function (3 variables);
        % 44=1-exp(-k*x); 45: Four-parameter logistic; 46=quadratic/parabola.
        %
        % peakfit(signal,center,window,NumPeaks,peakshape,extra) 
        % 'extra' specifies the value of 'extra', used only in the Voigt, Pearson,
        % exponentionally broadened Gaussian, Gaussian/Lorentzian blend, and
        % bifurcated Gaussian and Lorentzian shapes to fine-tune the peak shape.
        % 
        % peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials);
        % Performs "NumTrials" trial fits and selects the best one (with lowest
        % fitting error). NumTrials can be any positive integer (default is 1).
        %
        % peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start)
        % Specifies the first guesses vector "firstguess" for the peak positions
        % and widths. Must be expressed as a vector , in square brackets, e.g.
        % start=[position1 width1 position2 width2 ...]
        % 
        % peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero) 
        % 'autozero' sets the baseline correction mode:
        % autozero=0 (default) does not subtract baseline from data segment; 
        % autozero=1 interpolates a linear baseline from the edges of the data 
        %            segment and subtracts it from the signal (assumes that the 
        %            peak returns to the baseline at the edges of the signal); 
        % autozero=2 is like mode 1 except that it computes a quadratic curved baseline; 
        % autozero=3 compensates for a flat baseline without reference to the 
        %            signal itself (best if the peak does not return to the
        %            baseline at the edges of the signal).
        %
        % peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,...
        % autozero,fixedparameters)
        % 'fixedparameters' specifies fixed values for widths (shapes 10, 11) or
        % positions (shapes 16, 17)
        % 
        % peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,...
        % autozero,fixedparameters,plots)
        % 'plots' controls graphic plotting: 0=no plot; 1=plots draw as usual (default)
        % 
        % peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,...
        % autozero,fixedparameters,plots,bipolar)
        % 'bipolar' = 0 constrain peaks heights to be positions; 'bipolar' = 1
        % allows positive ands negative peak heights.
        %
        % peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,...
        % autozero,fixedparameters,plots,bipolar,minwidth)
        % 'minwidth' sets the minmimum allowed peak width. The default if not
        % specified is equal to the x-axis interval. Must be a vector of minimum
        % widths, one value for each peak, if the multiple peak shape is chosen, 
        % as in example 17 and 18.
        %
        % peakfit(signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,...
        % autozero,fixedparameters,plots,bipolar,minwidth)
        %  'DELTA' (14th input argument) controls the restart variance when
        %  NumTrials>1. Default value is 1.0. Larger values give more variance.
        %  Version 5.8 and later only. 
        % [FitResults,FitError]=peakfit(signal,center,window...) Returns the
        % FitResults vector in the order peak number, peak position, peak height,
        % peak width, and peak area), and the FitError (the percent RMS
        % difference between the data and the model in the selected segment of that
        % data) of the best fit.
        %
        % [FitResults,LowestError,BestStart,xi,yi,BootResults]=peakfit(signal,...)
        % Prints out parameter error estimates for each peak (bootstrap method).
        %
        % Optional output parameters 
        % 1. FitResults: a table of model peak parameters, one row for each peak,
        %    listing Peak number, Peak position, Height, Width, and Peak area.
        % 2. GOF: Goodness of Fit, a vector containing the rms fitting error of the
        %    best trial fit and the R-squared (coefficient of determination).
        % 3. Baseline, the polynomial coefficients of  the baseline in linear
        %    and quadratic baseline modes (1 and 2) or the value of the constant
        %    baseline in flat baseline mode.
        % 3. coeff: Coefficients for the polynomial fit (shape 28 only; for other 
        %    shapes, coeff=0)
        % 5. residual: the difference between the data and the best fit.
        % 6. xi: vector containing 600 interploated x-values for the model peaks. 
        % 7. yi: matrix containing the y values of each model peak at each xi. 
        %    Type plot(xi,yi(1,:)) to plot peak 1 or plot(xi,yi) to plot all peaks
        % 8. BootResults: a table of bootstrap precision results for a each peak
        %    and peak parameter.
        %  
        % Example 1a: Fits exp(-x)^2 with a single Gaussian peak model.
        % >> x=[0:.1:10]';y=exp(-(x-5).^2);peakfit([x y])
        %    Peak number  Peak position   Height     Width      Peak area
        %         1            5            1        1.665       1.7725
        %
        % Example 1b: Fits histogram of 10000 random numbers with a single Gaussian
        % >> [N,X]=hist(randn(size(1:10000)));peakfit([X;N])
        %   Peak number position  Height     Width     Area
        %     1.0000    0.0659   234.4104    2.4592   607.3360
        %
        % Example 1c: Fits small set of manually entered y data to a single Gaussian peak model.
        % >> y=[0 1 2 4 6 7 6 4 2 1 0 ];x=1:length(y);
        % >> peakfit([x;y],length(y)/2,length(y),0,0,0,0,0,0)
        %
        % Example 2:
        % x=[0:.01:10];y=exp(-(x-5).^2)+randn(size(x));peakfit([x;y])
        % Measurement of very noisy peak with signal-to-noise ratio = 1.
        % ans =
        %         1         5.0279       0.9272      1.7948      1.7716
        %
        % Example 3:
        % >> x=[0:.1:10];y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));
        % >> peakfit([x' y'],0,0,2)
        % Fits a noisy two-peak signal with a double Gaussian model (NumPeaks=2).
        % ans =
        %         1       3.0001      0.49489        1.642      0.86504
        %         2       4.9927       1.0016       1.6597       1.7696
        %
        % Example 4:
        % >> x=1:10;y=ones(size(x))./(1+(x-5).^2);peakfit(y,0,0,1,2)
        % Fit Lorentzian (peakshape=2) located at x=50, height=1, width=2.
        % ans =
        %        1           50      0.99974       1.9971       3.1079
        %
        % Example 5: 
        % >> x=[0:.005:1];y=humps(x);peakfit([x' y'],.3,.7,1,4,3);
        % Fits a portion of the humps function, 0.7 units wide and centered on 
        % x=0.3, with a single (NumPeaks=1) Pearson function (peakshape=4) with
        % extra=3 (controls shape of Pearson function).
        %
        % Example 6: 
        % >> x=[0:.005:1];y=(humps(x)+humps(x-.13)).^3;smatrix=[x' y'];
        % >> [FitResults,FitError]=peakfit(smatrix,.4,.7,2,1,0,10)
        % Creates a data matrix 'smatrix', fits a portion to a two-peak Gaussian
        % model, takes the best of 10 trials.  Returns FitResults and FitError.
        % FitResults =
        %          1      0.31056  2.0125e+006      0.11057  2.3689e+005
        %          2      0.41529  2.2403e+006      0.12033  2.8696e+005
        % FitError =
        %         1.1899
        %
        % Example 7:
        % >> peakfit([x' y'],.4,.7,2,1,0,10,[.3 .1 .5 .1]);
        % As above, but specifies the first-guess position and width of the two
        % peaks, in the order [position1 width1 position2 width2]
        %
        % Example 8: (Version 4 only)
        % Demonstration of the four autozero modes, for a single Gaussian on flat
        %  baseline, with position=10, height=1, and width=1.66. Autozero mode
        %  is specified by the 9th input argument (0,1,2, or 3).
        % >> x=8:.05:12;y=1+exp(-(x-10).^2);
        % >> [FitResults,FitError]=peakfit([x;y],0,0,1,1,0,1,0,0)
        % Autozero=0 means to ignore the baseline (default mode if not specified)
        % FitResults =
        %         1           10       1.8561        3.612       5.7641
        % FitError =
        %         5.387
        % >> [FitResults,FitError]=peakfit([x;y],0,0,1,1,0,1,0,1)
        % Autozero=1 subtracts linear baseline from edge to edge.
        % Does not work well because signal does not return to baseline at edges.
        % FitResults =
        %         1       9.9984      0.96153        1.559       1.5916
        % FitError =
        %        1.9801
        % >> [FitResults,FitError]=peakfit([x;y],0,0,1,1,0,1,0,2)
        % Autozero=1 subtracts quadratic baseline from edge to edge.
        % Does not work well because signal does not return to baseline at edges.
        % FitResults =
        %         1       9.9996      0.81749       1.4384       1.2503
        % FitError =
        %        1.8204
        % Autozero=3: Flat baseline mode, measures baseline by regression
        % >> [FitResults,FitError,Baseline]=peakfit([x;y],0,0,1,1,0,1,0,3)
        % FitResults =
        %         1           10       1.0001       1.6653       1.7645
        % Baseline =
        %     0.0037056
        % FitError =
        %       0.99985
        %
        % Example 9:
        % x=[0:.1:10];y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));
        % [FitResults,FitError]=peakfit([x' y'],0,0,2,11,0,0,0,0,[1.666 1.666])
        % Same as example 3, fit with fixed-width Gaussian (shape 11), width=1.666
        % 
        % Example 10: (Version 3 or later; Prints out parameter error estimates)
        % x=0:.05:9;y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.01*randn(1,length(x));
        % [FitResults,GOF,baseline,coeff,residual,xi,yi,BootResults]=peakfit([x;y],0,0,2,6,0,1,0,0,0);
        %
        % Example 11: (Version 3.2 or later)
        % x=[0:.005:1];y=humps(x);[FitResults,FitError]=peakfit([x' y'],0.54,0.93,2,13,15,10,0,0,0) 
        %
        % FitResults =
        %         1      0.30078       190.41      0.19131       23.064
        %         2      0.89788       39.552      0.33448       6.1999
        % FitError = 
        %       0.34502
        % Fits both peaks of the Humps function with a Gaussian/Lorentzian blend
        % (shape 13) that is 15% Gaussian (Extra=15).
        % 
        % Example 12:  (Version 3.2 or later)
        % >> x=[0:.1:10];y=exp(-(x-4).^2)+.5*exp(-(x-5).^2)+.01*randn(size(x));
        % >> [FitResults,FitError]=peakfit([x' y'],0,0,1,14,45,10,0,0,0) 
        % FitResults =
        %         1       4.2028       1.2315        4.077       2.6723
        % FitError =
        %       0.84461
        % Fit a slightly asymmetrical peak with a bifurcated Gaussian (shape 14)
        % 
        % Example 13:  (Version 3.3 or later)
        % >> x=[0:.1:10]';y=exp(-(x-5).^2);peakfit([x y],0,0,1,1,0,0,0,0,0,0)
        % Example 1 without plotting (11th input argument = 0, default is 1)
        % 
        % Example 14:  (Version 3.9 or later)
        % Exponentially broadened Lorentzian with position=9, height=1.
        % x=[0:.1:20]; 
        % L=lorentzian(x,9,1);
        % L1=ExpBroaden(L',-10)+0.02.*randn(size(x))';
        % [FitResults,FitError]=peakfit([x;L1'],0,0,1,18,10)
        %
        % Example 15: Fitting humps function with two Voigts (version 7.45)
        % FitResults,FitError]=peakfit(humps(0:.01:2),60,120,2,20,1.7,5,0)
        % FitResults =
        %         1       31.099       94.768       19.432       2375.3
        %         2       90.128       20.052       32.973       783.34
        % FitError =
        %       0.72829      0.99915
        %
        % Example 16: (Version 4.3 or later) Set +/- mode to 1 (bipolar)
        % x=[0:.1:10];y=exp(-(x-5).^2)-.5*exp(-(x-3).^2)+.1*randn(size(x));
        % peakfit([x' y'],0,0,2,1,0,1,0,0,0,1,1)
        % FitResults =
        %         1       3.1636      -0.5433         1.62      -0.9369
        %         2       4.9487      0.96859       1.8456       1.9029
        % FitError =
        %        8.2757
        %
        % Example 17: Version 5 or later. Fits humps function to a model consisting 
        % of one Pearson (shape=4, extra=3) and one Gaussian (shape=1), flat
        % baseline mode=3, NumTrials=10.
        % x=[0:.005:1.2];y=humps(x);[FitResults,FitError]=peakfit([x' y'],0,0,2,[2 1],[0 0])
        % FitResults =
        %         1      0.30154       84.671      0.27892       17.085
        %         2      0.88522       11.545      0.20825       2.5399
        % Baseline =
        %         0.901
        % FitError =
        %        10.457
        %
        % Example 18: 5 peaks, 5 different shapes, all heights = 1, widths = 3.
        % x=0:.1:60;
        % y=modelpeaks2(x,[1 2 3 4 5],[1 1 1 1 1],[10 20 30 40 50],...
        % [3 3 3 3 3],[0 0 0 2 -20])+.01*randn(size(x));
        % peakfit([x' y'],0,0,5,[1 2 3 4 5],[0 0 0 2 -20])
        %
        % Example 19: Minimum width constraint (13th input argument)
        % x=1:30;y=gaussian(x,15,8)+.05*randn(size(x));
        % No constraint:
        % peakfit([x;y],0,0,5,1,0,10,0,0,0,1,0,0);
        % Widths constrained to values above 7:
        % peakfit([x;y],0,0,5,1,0,10,0,0,0,1,0,7);
        %
        % Example 20: Noise test with peak height = RMS noise = 1.
        % x=[-5:.02:5];y=exp(-(x).^2)+randn(size(x));P=peakfit([x;y],0,0,1,1,0,10,0,0,0,1,1);
        %
        % Example 21a: Gaussian peak on strong sloped straight-line baseline, 2-peak
        % fit with variable-slope straight line ("third peak" is shape 26, peakfit version 6 only).
        % x=8:.05:12;y=x+exp(-(x-10).^2);
        % [FitResults,FitError]=peakfit([x;y],0,0,2,[1 26],[1 1],1,0)
        % FitResults =           
        %         1           10            1       1.6651       1.7642
        %         2        4.485      0.22297         0.05       40.045
        % FitError =0.093
        %
        % Example 21b: As above, but with TWO signal peaks, peakshape=[1 1 26],
        % helped by adding rough first guesses {'start') using the 'polyfit'
        % function to generate first guesses for the sloping baseline:
        %  x=7:.05:13;y=x/2+exp(-(x-9).^2)+exp(-(x-11).^2);
        % start=[8 1 10 1 polyfit(x,y,1)];
        % [FitResults,FitError]=peakfit([x;y],0,0,3,[1 1 26],[1 1 1],1,start)
        %
        % Example 22: Segmented linear fit (Shape 29, peakfit version 6 only):
        % x=[0.9:.005:1.7];y=humps(x);
        % peakfit([x' y'],0,0,9,29,0,10,0,0,0,1,1)
        %
        % Example 23: Sixth degree polynomial fit (Shape 28, peakfit version 6 only)
        %   x=[0.3:.005:1.7];y=humps(x);y=y+cumsum(y);
        %   peakfit([x' y'],0,0,1,28,6,10,0,0,0,1,1)
        %
        % Example 24: Effect of quantization of independent (x) and dependent (y) variables.
        % x=.5+round([2:.02:7.5]);y=exp(-(x-5).^2)+randn(size(x))/10;peakfit([x;y])
        % x=[2:.01:8];y=exp(-(x-5).^2)+.1.*randn(size(x));y=.1.*round(10.*y);peakfit([x;y])
        %
        % Example 25: Variable-alpha Voigt functions, shape 30. (Version 7.45 and
        % above only). FitResults has an added 6th column for the measured alphas
        % of each peak.
        % x=[0:.005:1.3];y=humps(x);[FitResults,FitError]=peakfit([x' y'],60,120,2,30,15)
        % FitResults =
        %      1      0.30147       95.797      0.19391       24.486       2.2834
        %      2      0.89186       19.474      0.33404       7.3552      0.39824
        % FitError =
        %       0.84006      0.99887
        %
        % Example 26: Variable time constant exponentially-broadened Gaussian
        % function, shape 31. (Version 7 and above only). FitResults has an added
        % 6th column for the measured  time constant of each peak.
        % [FitResults,FitError]=peakfit(DataMatrix3,1860.5,364,2,31,32.9731,5,[1810 60 30 1910 60 30])
        % FitResults =
        %     1       1800.6       1.8581       60.169       119.01       32.781
        %     2       32.781      0.48491       1900.4        30.79       33.443
        % FitError =
        %      0.076651      0.99999
        %
        % Example 27 Pearson variable shape, PeakShape=32,(Version 7 and above
        % only). Requires modelpeaks2 function in path.
        % x=1:.1:30;y=modelpeaks2(x,[4 4],[1 1],[10 20],[5 5],[1 10]);
        % [FitResults,FitError]=peakfit([x;y],0,0,2,32,10,5)
        %
        % Example 28 Gaussian/Lorentzian blend variable shape, PeakShape=33
        % (Version 7 and above only). Requires modelpeaks2 function in path.
        % x=1:.1:30;y=modelpeaks2(x,[13 13],[1 1],[10 20],[3 3],[20 80]);
        % [FitResults,FitError]=peakfit([x;y],0,0,2,33,0,5)
        %
        % Example 29:  Fixed-position Gaussian (shape 16), positions=[3 5]. 
        % x=0:.1:10;y=exp(-(x-5).^2)+.5*exp(-(x-3).^2)+.1*randn(size(x));
        % [FitResults,FitError]=peakfit([x' y'],0,0,2,16,0,0,0,0,[3 5])
        %
        % Example 30: Fixed-width Gaussian/Lorentzian blend (shape 35), 
        % width=[3 3], compared to variable-width fit (shape 13). Fitting
        % error is larger but nevertheless results are more accurate.
        % x=0:.1:10;y=GL(x,4,3,50)+.5*GL(x,6,3,50)+.1*randn(size(x));
        % [FitResults,FitError]=peakfit([x;y],0,0,2,13,50,1)
        %          1       4.0632       1.0545       3.2182       3.9242
        %          2       6.2736      0.41234       2.8114       1.3585
        %  FitError =  6.4311      0.95211 
        % [FitResults,FitError]=peakfit([x;y],0,0,2,35,50,1,0,0,[3 3])
        % FitResults =
        %          1       3.9527       1.0048            3       3.5138
        %          2       6.1007       0.5008            3       1.7502
        % 
        %  FitError =  6.4783      0.95141
        %
        % Example 31: Variable time constant exponentially-broadened Lorentzian
        % function, shape 38. (Version 7.7 and above only). FitResults has an added
        % 6th column for the measured time constant.
        % x=[1:100]';
        % y=explorentzian(x',40,5,-10)+.01*randn(size(x));
        % peakfit([x y],0,0,1,38,0,10)
        %
        % Example 32: 3-parameter logistic (Gompertz), shape 43. (Version 7.9 and
        % above only). Parameters labeled Bo, Kh, and L. FitResults extended to 6
        % columns. 
        % t=0:.1:10; Bo=6;Kh=3;L=4;
        % y=Bo*exp(-exp((Kh*exp(1)/Bo)*(L-t)+1))+.1.*randn(size(t));
        % [FitResults,GOF]=peakfit([t;y],0,0,1,43)
        %
        % Example 33: Shape 46, 'quadslope' (version 8 and above). Two overlapping
        % Gaussians (positions=9,11; heights=1; widths=1.666) on a curved baseline,
        % using a 3-peak fit with peakshape=[1 1 46], default NumTrials and start.
        % x=7:.05:13;
        % y=x.^2/50+exp(-(x-9).^2)+exp(-(x-11).^2)+.01.*randn(size(x));
        % [FitResults,FitError]=peakfit([x;y],0,0,3,[1 1 46],[1 1 1])

        % Copyright (c) 2015, Thomas C. O'Haver

        % Permission is hereby granted, free of charge, to any person obtaining a copy
        % of this software and associated documentation files (the "Software"), to deal
        % in the Software without restriction, including without limitation the rights
        % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        % copies of the Software, and to permit persons to whom the Software is
        % furnished to do so, subject to the following conditions:
        % 
        % The above copyright notice and this permission notice shall be included in
        % all copies or substantial portions of the Software.
        % 
        % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
        % THE SOFTWARE.
        % 
        global AA xxx PEAKHEIGHTS FIXEDPARAMETERS AUTOZERO delta BIPOLAR CLIPHEIGHT
        % format short g
        format compact
        warning off all
        NumArgOut=nargout;
        datasize=size(signal);
        if datasize(1)<datasize(2),signal=signal';end
        datasize=size(signal);
        if datasize(2)==1, %  'signal' is a vector; Must be peakfit(Y-vector)
            X=1:length(signal); % Create an independent variable vector
            Y=signal;
        else
            % 'signal' is a matrix. Must be peakfit(DataMatrix)
            X=signal(:,1); % Split matrix argument 
            Y=signal(:,2);
        end
        X=reshape(X,1,length(X)); % Adjust X and Y vector shape to 1 x n (rather than n x 1)
        Y=reshape(Y,1,length(Y));
        % If necessary, flip the data vectors so that X increases
        if X(1)>X(length(X)),
            disp('X-axis flipped.')
            X=fliplr(X);
            Y=fliplr(Y);
        end

        % Isolate desired segment from data set for curve fitting
        if nargin==1 || nargin==2,center=(max(X)-min(X))/2;window=max(X)-min(X);end
        % Y=Y-min(Y);
        xoffset=0;
        n1=LibAv05.val2ind(X,center-window/2);
        n2=LibAv05.val2ind(X,center+window/2);
        if window==0,n1=1;n2=length(X);end
        xx=X(n1:n2)-xoffset;
        yy=Y(n1:n2);
        ShapeString='Gaussian';
        coeff=0;
        CLIPHEIGHT=max(Y);
        LOGPLOT=0;
        % Define values of any missing arguments
        % (signal,center,window,NumPeaks,peakshape,extra,NumTrials,start,autozero,fixedparameters,plots,bipolar,minwidth,DELTA)
        switch nargin
            case 1
                NumPeaks=1;
                peakshape=1;
                extra=1;
                NumTrials=1;
                xx=X;yy=Y;
                start=LibAv05.calcstart(xx,NumPeaks,xoffset);
                AUTOZERO=0;
                plots=1;
                BIPOLAR=0;
                MINWIDTH=xx(2)-xx(1);
                delta=1;
                CLIPHEIGHT=max(Y);
            case 2
                NumPeaks=1;
                peakshape=1;
                extra=1;
                NumTrials=1;
                xx=signal;yy=center;
                start=LibAv05.calcstart(xx,NumPeaks,xoffset);
                AUTOZERO=0;
                plots=1;
                BIPOLAR=0;
                MINWIDTH=xx(2)-xx(1);
                delta=1;
                CLIPHEIGHT=max(Y);
            case 3
                NumPeaks=1;
                peakshape=1;
                extra=1;
                NumTrials=1;
                start=LibAv05.calcstart(xx,NumPeaks,xoffset);
                AUTOZERO=0;
                FIXEDPARAMETERS=0;
                plots=1;
                BIPOLAR=0;
                MINWIDTH=xx(2)-xx(1);
                delta=1;
                CLIPHEIGHT=max(Y);
            case 4 % Numpeaks specified in arguments
                peakshape=1;
                extra=1;
                NumTrials=1;
                start=LibAv05.calcstart(xx,NumPeaks,xoffset);
                AUTOZERO=0;
                FIXEDPARAMETERS=0;
                plots=1;
                BIPOLAR=0;
                MINWIDTH=xx(2)-xx(1);
                delta=1;
                CLIPHEIGHT=max(Y);
            case 5 % Numpeaks, peakshape specified in arguments
                extra=ones(1,NumPeaks);
                NumTrials=1;
                start=LibAv05.calcstart(xx,NumPeaks,xoffset);
                AUTOZERO=0;
                FIXEDPARAMETERS=0;
                plots=1;
                BIPOLAR=0;
                MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
                delta=1;
                CLIPHEIGHT=max(Y);
            case 6
                NumTrials=1;
                start=LibAv05.calcstart(xx,NumPeaks,xoffset);
                AUTOZERO=0;
                FIXEDPARAMETERS=0;
                plots=1;
                BIPOLAR=0;
                MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
                delta=1;
            case 7
                start=LibAv05.calcstart(xx,NumPeaks,xoffset);
                AUTOZERO=0;
                FIXEDPARAMETERS=0;
                plots=1;
                BIPOLAR=0;
                MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
                delta=1;
                CLIPHEIGHT=max(Y);
            case 8
                AUTOZERO=0;
                FIXEDPARAMETERS=0;
                plots=1;
                BIPOLAR=0;
                MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
                delta=1;
                CLIPHEIGHT=max(Y);
            case 9
                % initialstart=start % testing
                AUTOZERO=autozero;
                FIXEDPARAMETERS=0;
                plots=1;
                BIPOLAR=0;
                MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
                delta=1;
            case 10
                AUTOZERO=autozero;
                FIXEDPARAMETERS=fixedparameters;
                plots=1;
                BIPOLAR=0;
                MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
                delta=1;
            case 11
                AUTOZERO=autozero;
                FIXEDPARAMETERS=fixedparameters;
                BIPOLAR=0;
                MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
                delta=1;
                CLIPHEIGHT=max(Y);
            case 12
                AUTOZERO=autozero;
                FIXEDPARAMETERS=fixedparameters;
                BIPOLAR=bipolar;
                MINWIDTH=zeros(size(peakshape))+(xx(2)-xx(1));
                delta=1;
                CLIPHEIGHT=max(Y);
            case 13
                AUTOZERO=autozero;
                FIXEDPARAMETERS=fixedparameters;
                BIPOLAR=bipolar;
                MINWIDTH=minwidth;
                delta=1;
            case 14
                AUTOZERO=autozero;
                FIXEDPARAMETERS=fixedparameters;
                BIPOLAR=bipolar;
                MINWIDTH=minwidth;
                delta=DELTA;
                CLIPHEIGHT=max(Y);
            case 15
                AUTOZERO=autozero;
                FIXEDPARAMETERS=fixedparameters;
                BIPOLAR=bipolar;
                MINWIDTH=minwidth;
                delta=DELTA;
                CLIPHEIGHT=clipheight;
            otherwise
        end % switch nargin
        % Saturation Code, skips points greater than set maximum
        if CLIPHEIGHT<max(Y),
            apnt=1;
            for pnt=1:length(xx),
                if yy(pnt)<CLIPHEIGHT,
                    axx(apnt)=xx(pnt);
                    ayy(apnt)=yy(pnt);
                    apnt=apnt+1;
                end
            end
            xx=axx;yy=ayy;
        end
        % Default values for placeholder zeros1
        if NumTrials==0;NumTrials=1;end
        shapesvector=peakshape;
        if isscalar(peakshape),
        else
            % disp('peakshape is vector');
            shapesvector=peakshape;
            NumPeaks=length(peakshape);
            peakshape=22;
        end
        if peakshape==0;peakshape=1;end
        if NumPeaks==0;NumPeaks=1;end
        if start==0;start=LibAv05.calcstart(xx,NumPeaks,xoffset);end
        firststart=start; % <<<<<<<<<<<
        if FIXEDPARAMETERS==0, FIXEDPARAMETERS=length(xx)/10;end
        if peakshape==16;FIXEDPOSITIONS=fixedparameters;end
        if peakshape==17;FIXEDPOSITIONS=fixedparameters;end
        if AUTOZERO>3,AUTOZERO=3;disp('AUTPOZERO must be between 0 and 3');end
        if AUTOZERO<0,AUTOZERO=0;disp('AUTPOZERO must be between 0 and 3');end
        Heights=zeros(1,NumPeaks);
        FitResults=zeros(NumPeaks,6);

        % % Remove linear baseline from data segment if AUTOZERO==1
        baseline=0;
        bkgcoef=0;
        bkgsize=round(length(xx)/10);
        if bkgsize<2,bkgsize=2;end
        lxx=length(xx);
        if AUTOZERO==1, % linear autozero operation  
            XX1=xx(1:round(lxx/bkgsize));
            XX2=xx((lxx-round(lxx/bkgsize)):lxx);
            Y1=yy(1:(round(length(xx)/bkgsize)));
            Y2=yy((lxx-round(lxx/bkgsize)):lxx);
            bkgcoef=polyfit([XX1,XX2],[Y1,Y2],1);  % Fit straight line to sub-group of points
            bkg=polyval(bkgcoef,xx);
            yy=yy-bkg;
        end % if
        if AUTOZERO==2, % Quadratic autozero operation  
            XX1=xx(1:round(lxx/bkgsize));
            XX2=xx((lxx-round(lxx/bkgsize)):lxx);
            Y1=yy(1:round(length(xx)/bkgsize));
            Y2=yy((lxx-round(lxx/bkgsize)):lxx);
            bkgcoef=polyfit([XX1,XX2],[Y1,Y2],2);  % Fit parabola to sub-group of points
            bkg=polyval(bkgcoef,xx);
            yy=yy-bkg;
        end % if autozero

        PEAKHEIGHTS=zeros(1,NumPeaks);
        n=length(xx);
        newstart=start;
        % Assign ShapStrings
        switch peakshape(1)
            case 1
                ShapeString='Gaussian';
            case 2
                ShapeString='Lorentzian';
            case 3
                ShapeString='Logistic';
            case 4
                ShapeString='Pearson';
            case 5
                ShapeString='ExpGaussian';
            case 6
                ShapeString='Equal width Gaussians';
            case 7
                ShapeString='Equal width Lorentzians';
            case 8
                ShapeString='Exp. equal width Gaussians';
            case 9
                ShapeString='Exponential Pulse';
            case 10
                ShapeString='Up Sigmoid (logistic function)';
            case 23
                ShapeString='Down Sigmoid (logistic function)';  
            case 11
                ShapeString='Fixed-width Gaussian';
            case 12
                ShapeString='Fixed-width Lorentzian';
            case 13
                ShapeString='Gaussian/Lorentzian blend';
            case 14
                ShapeString='BiGaussian';    
            case 15
                ShapeString='Breit-Wigner-Fano';   
            case 16
                ShapeString='Fixed-position Gaussians';
            case 17
                ShapeString='Fixed-position Lorentzians';
            case 18
                ShapeString='Exp. Lorentzian';
            case 19
                ShapeString='Alpha function';
            case 20
                ShapeString='Voigt (equal alphas)';
            case 21
                ShapeString='triangular';
            case 22
                ShapeString=num2str(shapesvector);
            case 24
                ShapeString='Negative Binomial Distribution';
            case 25
                ShapeString='Lognormal Distribution';
            case 26
                ShapeString='slope';
            case 27
                ShapeString='First derivative';
            case 28
                ShapeString='Polynomial';
            case 29
                ShapeString='Segmented linear';
            case 30
                ShapeString='Voigt (variable alphas)';
            case 31
                ShapeString='ExpGaussian (var. time constant)';
            case 32
                ShapeString='Pearson (var. shape constant)';
            case 33
                ShapeString='Variable Gaussian/Lorentzian';
            case 34
                ShapeString='Fixed-width Voigt';
            case 35
                ShapeString='Fixed-width G/L blend';
            case 36
                ShapeString='Fixed-width ExpGaussian';
            case 37
                ShapeString='Fixed-width Pearson';
            case 38
                ShapeString='ExpLorentzian (var. time constant)'; 
            case 40
                ShapeString='Sine wave';
            case 41
                ShapeString='Rectangular pulse';
            case 42
                ShapeString='Flattened Gaussian';  
            case 43
                ShapeString='3-parameter Gompertz.';  
            case 44
                ShapeString='1-exp(-k*t)';  
            case 45
                ShapeString='4-parameter logistic'; 
            case 46
                ShapeString='Quadratic Baseline'; 
            case 47
                ShapeString='LibAv05.blackbody'; 
            otherwise
        end % switch peakshape

        % Perform peak fitting for selected peak shape using fminsearch function
        options = optimset('TolX',.0000001,'Display','off','MaxFunEvals',1000 );
        LowestError=1000; % or any big number greater than largest error expected
        FitParameters=zeros(1,NumPeaks.*3); 
        BestStart=zeros(1,NumPeaks.*2); 
        height=zeros(1,NumPeaks); 
        bestmodel=zeros(size(yy));
        TrialParameters=zeros(1,NumPeaks.*3);
        for k=1:NumTrials,
            % StartMatrix(k,:)=newstart;
            % disp(['Trial number ' num2str(k) ] ) % optionally prints the current trial number as progress indicator
            switch peakshape(1)
                case 1
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitgaussian(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 2
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitlorentzian(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 3
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitlogistic(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 4
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitpearson(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 5
                    zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
                    zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 6
                    cwnewstart(1)=newstart(1);
                    for pc=2:NumPeaks,
                        cwnewstart(pc)=newstart(2.*pc-1);
                    end
                    cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitewgaussian(lambda,xx,yy)),cwnewstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(NumPeaks+1)<MINWIDTH,
                            TrialParameters(NumPeaks+1)=MINWIDTH;
                        end
                    end
                case 7
                    cwnewstart(1)=newstart(1);
                    for pc=2:NumPeaks,
                        cwnewstart(pc)=newstart(2.*pc-1);
                    end
                    cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitewlorentzian(lambda,xx,yy)),cwnewstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(NumPeaks+1)<MINWIDTH,
                            TrialParameters(NumPeaks+1)=MINWIDTH;
                        end
                    end
                case 8
                    cwnewstart(1)=newstart(1);
                    for pc=2:NumPeaks,
                        cwnewstart(pc)=newstart(2.*pc-1);
                    end
                    cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(NumPeaks+1)<MINWIDTH,
                            TrialParameters(NumPeaks+1)=MINWIDTH;
                        end
                    end
                case 9
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexppulse(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 10
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitupsigmoid(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 23
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitdownsigmoid(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 11
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end
                    % fixedstart=fixedstart
                    TrialParameters=fminsearch(@(lambda)(LibAv05.FitFWGaussian(lambda,xx,yy)),fixedstart,options);
                case 12
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.FitFWLorentzian(lambda,xx,yy)),fixedstart,options);
                case 13
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitGL(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 14
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitBiGaussian(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 15
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitBWF(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 16
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
                        fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.FitFPGaussian(lambda,xx,yy)),fixedstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(Peak)<MINWIDTH,
                            TrialParameters(Peak)=MINWIDTH;
                        end
                    end
                case 17
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
                        fixedstart(pc)=fixedstart(pc)+.1*(rand-.5).*fixedstart(pc);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(Peak)<MINWIDTH,
                            TrialParameters(Peak)=MINWIDTH;
                        end
                    end
                case 18
                    zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
                    zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 19
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitalphafunction(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 20
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitvoigt(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 21
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fittriangular(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 22
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitmultiple(lambda,xx,yy,NumPeaks,shapesvector,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH(Peak),
                            TrialParameters(2*Peak)=MINWIDTH(Peak);
                        end
                    end
                case 24
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitnbinpdf(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 25
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitlognpdf(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 26
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitlinslope(lambda,xx,yy)),polyfit(xx,yy,1),options);
                     coeff=TrialParameters;
                case 27
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitd1gauss(lambda,xx,yy)),newstart,options);
                case 28
                    coeff=LibAv05.fitpolynomial(xx,yy,extra);
                    TrialParameters=coeff;
                case 29
                    cnewstart(1)=newstart(1);
                    for pc=2:NumPeaks,
                        cnewstart(pc)=newstart(2.*pc-1)+(delta*(rand-.5)/50);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitsegmented(lambda,xx,yy)),cnewstart,options);
                case 30
                    % newstart=newstart % testing
                    nn=max(xx)-min(xx);
                    start=[];
                    startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
                    for marker=1:NumPeaks,
                        markx=startpos(marker)+ xoffset;
                        start=[start markx nn/5 extra];
                    end % for marker
                     newstart=start;
                    for parameter=1:3:3*NumPeaks,
                        newstart(parameter)=newstart(parameter)*(1+randn/100);
                        newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                        newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
                    end
        %             case30newstart=newstart % uncomment for testing
        %             sizestartcase30=size(start) % uncomment for testing

                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitvoigtv(lambda,xx,yy)),start);
                 case 31
                    nn=max(xx)-min(xx);
                    start=[];
                    startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
                    for marker=1:NumPeaks,
                        markx=startpos(marker)+ xoffset;
                        start=[start markx nn/5 extra];
                    end % for marker
                     newstart=start;
                    for parameter=1:3:3*NumPeaks,
                        newstart(parameter)=newstart(parameter)*(1+randn/100);
                        newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                        newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
                    end
                    % case31newstart=newstart % uncomment for testing
                    % sizestartcase31=size(start) % uncomment for testing
                    zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
                    zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexpgaussianv(lambda,zxx,zyy)),newstart);
                case 32
                    nn=max(xx)-min(xx);
                    start=[];
                    startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
                    for marker=1:NumPeaks,
                        markx=startpos(marker)+ xoffset;
                        start=[start markx nn/5 extra];
                    end % for marker
                     newstart=start;
                    for parameter=1:3:3*NumPeaks,
                        newstart(parameter)=newstart(parameter)*(1+randn/100);
                        newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                        newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
                    end
                    case32newstart=newstart % uncomment for testing
                    sizestartcase32=size(start) % uncomment for testing
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitpearsonv(lambda,xx,yy)),newstart);
                case 33
                     nn=max(xx)-min(xx);
                    start=[];
                    startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
                    for marker=1:NumPeaks,
                        markx=startpos(marker)+ xoffset;
                        start=[start markx nn/5 extra];
                    end % for marker
                     newstart=start;
                    for parameter=1:3:3*NumPeaks,
                        newstart(parameter)=newstart(parameter)*(1+randn/100);
                        newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                        newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
                    end
                    % newstart=newstart % uncomment for testing
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitGLv(lambda,xx,yy)),newstart);
                case 34
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end  
                    % fixedstart=fixedstart % uncomment for testing
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitFWVoigt(lambda,xx,yy,extra)),fixedstart,options);
                case 35
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end            
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitFWGL(lambda,xx,yy,extra)),fixedstart,options);
                case 36
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end            
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitFWExpGaussian(lambda,xx,yy,extra)),fixedstart,options);
                case 37
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end            
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitFWPearson(lambda,xx,yy,extra)),fixedstart,options);
                case 38
                    nn=max(xx)-min(xx);
                    start=[];
                    startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
                    for marker=1:NumPeaks,
                        markx=startpos(marker)+ xoffset;
                        start=[start markx nn/5 extra];
                    end % for marker
                     newstart=start;
                    for parameter=1:3:3*NumPeaks,
                        newstart(parameter)=newstart(parameter)*(1+randn/100);
                        newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                        newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
                    end
                    % newstart=newstart
                    zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
                    zyy=[ones(size(yy)).*yy(1) yy zeros(size(yy)).*yy(length(yy)) ];
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexplorentzianv(lambda,zxx,zyy)),newstart);
                case 40
                     TrialParameters=fminsearch(@(lambda)(LibAv05.fitsine(lambda,xx,yy)),newstart,options); 
                case 41
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitrectangle(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 42
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitngaussian(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 43
                     nn=max(xx)-min(xx);
                    start=[];
                    startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
                    for marker=1:NumPeaks,
                        markx=startpos(marker)+ xoffset;
                        start=[start markx nn/5 extra];
                    end % for marker
                     newstart=start;
                    for parameter=1:3:3*NumPeaks,
                        newstart(parameter)=newstart(parameter)*(1+randn/100);
                        newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                        newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitGompertz(lambda,xx,yy)),newstart);

                case 44
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitOneMinusExp(lambda,xx,yy)),newstart,options);
                 case 45
                     nn=max(xx)-min(xx);
                    start=[];
                    startpos=[nn/(NumPeaks+1):nn/(NumPeaks+1):nn-(nn/(NumPeaks+1))]+min(xx);
                    for marker=1:NumPeaks,
                        markx=startpos(marker)+ xoffset;
                        start=[start markx nn/5 extra];
                    end % for marker
                     newstart=start;
                    for parameter=1:3:3*NumPeaks,
                        newstart(parameter)=newstart(parameter)*(1+randn/100);
                        newstart(parameter+1)=newstart(parameter+1)*(1+randn/20);
                        newstart(parameter+2)=newstart(parameter+1)*(1+randn/20);
                    end
                    % newstart=newstart % uncomment for testing
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitFourPL(lambda,xx,yy)),newstart);
                case 46
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitquadslope(lambda,xx,yy)),newstart,options);
                case 47
                    bbstart=3000;
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitblackbody(lambda,xx,yy)),bbstart,options);
                otherwise
            end % switch peakshape

            % Check variables
            % sizeNewstart=size(newstart) % uncomment for testing

        % Construct model from Trial parameters
        A=zeros(NumPeaks,n);
        for m=1:NumPeaks,
            switch peakshape(1)
                case 1
                    A(m,:)=LibAv05.gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 2
                    A(m,:)=LibAv05.lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 3
                    A(m,:)=LibAv05.logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 4
                    A(m,:)=LibAv05.pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
                case 5
                    A(m,:)=LibAv05.expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
                case 6
                    A(m,:)=LibAv05.gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
                case 7
                    A(m,:)=LibAv05.lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
                case 8
                    A(m,:)=LibAv05.expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';
                case 9
                    A(m,:)=LibAv05.exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 10
                    A(m,:)=LibAv05.upsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 11
                    A(m,:)=LibAv05.gaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
                case 12
                    A(m,:)=LibAv05.lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
                case 13
                    A(m,:)=LibAv05.GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
                case 14
                    A(m,:)=LibAv05.BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
                case 15
                    A(m,:)=LibAv05.BWF(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
                case 16
                    A(m,:)=LibAv05.gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
                case 17
                    A(m,:)=LibAv05.lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
                case 18
                    A(m,:)=LibAv05.explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
                case 19
                    A(m,:)=LibAv05.alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 20
                    A(m,:)=LibAv05.voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);        
                case 21
                    A(m,:)=LibAv05.triangular(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 22
                    A(m,:)=LibAv05.peakfunction(shapesvector(m),xx,TrialParameters(2*m-1),TrialParameters(2*m),extra(m));        
                case 23
                    A(m,:)=LibAv05.downsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));        
                case 24
                    A(m,:)=LibAv05.nbinpdf(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 25
                    A(m,:)=LibAv05.lognormal(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 26
                    A(m,:)=LibAv05.linslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 27
                    A(m,:)=LibAv05.d1gauss(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 28
                    A(m,:)=LibAv05.polynomial(xx,coeff);
                case 29
                    A(m,:)=LibAv05.segmented(xx,yy,PEAKHEIGHTS);
                case 30
                    A(m,:)=LibAv05.voigt(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
                case 31
                    A(m,:)=LibAv05.expgaussian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),-TrialParameters(3*m));        
                case 32
                    A(m,:)=LibAv05.pearson(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
                case 33
                    A(m,:)=LibAv05.GL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));
                case 34
                     width(m)=abs(FIXEDPARAMETERS(m));
        %                 gD(m)=width(m);
        %                 gL(m)=extra.*gD(m);
        %                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2))
                    A(m,:)=LibAv05.voigt(xx,TrialParameters(m), width(m),extra);
                case 35
                    A(m,:)=LibAv05.GL(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
                case 36
                    A(m,:)=LibAv05.expgaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
                case 37
                    A(m,:)=LibAv05.pearson(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);    
                case 38
                    A(m,:)=LibAv05.explorentzian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),-TrialParameters(3*m));        
                case 40
                    A(m,:)=LibAv05.sine(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 41
                    A(m,:)=LibAv05.rectangle(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 42
                    A(m,:)=LibAv05.ngaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
                case 43
                    A(m,:)=LibAv05.Gompertz(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
                case 44
                    A(m,:)=LibAv05.OneMinusExp(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 45
                    A(m,:)=LibAv05.FourPL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
                case 46
                    A(m,:)=LibAv05.quadslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                case 47
                    A(m,:)=LibAv05.blackbody(xx,TrialParameters(m));
            end % switch
            xxrange=max(xx)-min(xx);
            for parameter=1:2:2*NumPeaks,
                newstart(parameter)=newstart(parameter)+(xxrange.*(delta*randn)./(NumPeaks+1));
                newstart(parameter+1)=newstart(parameter+1)*(1+delta*(rand-.5)/100);
            end
        end % for NumPeaks
        % newstart=newstart; % <<<<<<<<<< % error check
        % Multiplies each row by the corresponding amplitude and adds them up
        if peakshape(1)==29, % Segmented linear
            model=LibAv05.segmented(xx,yy,PEAKHEIGHTS);
            TrialParameters=PEAKHEIGHTS;
            Heights=ones(size(PEAKHEIGHTS));
        else
            if AUTOZERO==3,
                baseline=PEAKHEIGHTS(1);
                Heights=PEAKHEIGHTS(2:1+NumPeaks);
                model=Heights'*A+baseline;
            else
        %          sizePeakHeights=size(PEAKHEIGHTS) %  % uncomment for testing
        %          SizeA=size(A)  % uncomment for testing
                model=PEAKHEIGHTS'*A;
                Heights=PEAKHEIGHTS;
                baseline=0;
            end
        end
        if peakshape(1)==28, % polynomial;
            model=LibAv05.polynomial(xx,coeff);
            TrialParameters=PEAKHEIGHTS;
            Heights=ones(size(PEAKHEIGHTS));
        end
        % Compare trial model to data segment and compute the fit error
            MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
          % Take only the single fit that has the lowest MeanFitError
          if MeanFitError<LowestError, 
              if min(Heights)>=-BIPOLAR*10^100,  % Consider only fits with positive peak heights
                LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
                FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
                BestStart=newstart; % Assign BestStart to the start with the lowest MeanFitError
                height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
                bestmodel=model; % Assign bestmodel to the model with the lowest MeanFitError
              end % if min(PEAKHEIGHTS)>0
          end % if MeanFitError<LowestError
        %  ErrorVector(k)=MeanFitError; %  % uncomment for testing
        end % for k (NumTrials)
            Rsquared=1-(norm(yy-bestmodel)./norm(yy-mean(yy)));
            SStot=sum((yy-mean(yy)).^2);
            SSres=sum((yy-bestmodel).^2);
            Rsquared=1-(SSres./SStot);
            GOF=[LowestError Rsquared];
        % Uncomment following 4 lines to monitor trail fit starts and errors.
        % StartMatrix=StartMatrix;
        % ErrorVector=ErrorVector;
        % matrix=[StartMatrix ErrorVector']
        % std(StartMatrix)
        % Construct model from best-fit parameters
        AA=zeros(NumPeaks,600);
        xxx=linspace(min(xx),max(xx),600);
        % xxx=linspace(min(xx)-length(xx),max(xx)+length(xx),200);
        for m=1:NumPeaks,
           switch peakshape(1)
            case 1
                AA(m,:)=LibAv05.gaussian(xxx,FitParameters(2*m-1),FitParameters(2*m));
            case 2
                AA(m,:)=LibAv05.lorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m));
            case 3
                AA(m,:)=LibAv05.logistic(xxx,FitParameters(2*m-1),FitParameters(2*m));
            case 4
                AA(m,:)=LibAv05.pearson(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
            case 5
                AA(m,:)=LibAv05.expgaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
            case 6
                AA(m,:)=LibAv05.gaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
            case 7
                AA(m,:)=LibAv05.lorentzian(xxx,FitParameters(m),FitParameters(NumPeaks+1));
            case 8
                AA(m,:)=LibAv05.expgaussian(xxx,FitParameters(m),FitParameters(NumPeaks+1),-extra*length(xxx)./length(xx))';
            case 9
                AA(m,:)=LibAv05.exppulse(xxx,FitParameters(2*m-1),FitParameters(2*m));  
            case 10
                AA(m,:)=LibAv05.upsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));   
            case 11
                AA(m,:)=LibAv05.gaussian(xxx,FitParameters(m),FIXEDPARAMETERS(m));
            case 12
                AA(m,:)=LibAv05.lorentzian(xxx,FitParameters(m),FIXEDPARAMETERS(m));
            case 13
                AA(m,:)=LibAv05.GL(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
            case 14
                AA(m,:)=LibAv05.BiGaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
            case 15
                AA(m,:)=LibAv05.BWF(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
            case 16
                AA(m,:)=LibAv05.gaussian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
            case 17
                AA(m,:)=LibAv05.lorentzian(xxx,FIXEDPOSITIONS(m),FitParameters(m));
            case 18
                AA(m,:)=LibAv05.explorentzian(xxx,FitParameters(2*m-1),FitParameters(2*m),-extra*length(xxx)./length(xx))';
            case 19
                AA(m,:)=LibAv05.alphafunction(xxx,FitParameters(2*m-1),FitParameters(2*m));
            case 20
                AA(m,:)=LibAv05.voigt(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);       
            case 21
                AA(m,:)=LibAv05.triangular(xxx,FitParameters(2*m-1),FitParameters(2*m));
            case 22
                AA(m,:)=LibAv05.peakfunction(shapesvector(m),xxx,FitParameters(2*m-1),FitParameters(2*m),extra(m));        
            case 23
                AA(m,:)=LibAv05.downsigmoid(xxx,FitParameters(2*m-1),FitParameters(2*m));  
            case 24
                AA(m,:)=LibAv05.nbinpdf(xxx,FitParameters(2*m-1),FitParameters(2*m));    
            case 25
                AA(m,:)=LibAv05.lognormal(xxx,FitParameters(2*m-1),FitParameters(2*m));    
            case 26
                AA(m,:)=LibAv05.linslope(xxx,FitParameters(2*m-1),FitParameters(2*m));   
            case 27
                AA(m,:)=LibAv05.d1gauss(xxx,FitParameters(2*m-1),FitParameters(2*m));  
            case 28
                AA(m,:)=LibAv05.polynomial(xxx,coeff);
            case 29
            case 30
                AA(m,:)=LibAv05.voigt(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));        
            case 31
                AA(m,:)=LibAv05.expgaussian(xxx,FitParameters(3*m-2),FitParameters(3*m-1),-FitParameters(3*m)*length(xxx)./length(xx));        
            case 32
                AA(m,:)=LibAv05.pearson(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));        
            case 33
                AA(m,:)=LibAv05.GL(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m)); 
            case 34
                          width(m)=abs(FIXEDPARAMETERS(m));
        %                 gD(m)=width(m);
        %                 gL(m)=extra.*gD(m);
        %                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 +
        %                 gD(m).^2));
                AA(m,:)=LibAv05.voigt(xxx,FitParameters(m),width(m),extra);
            case 35
                AA(m,:)=LibAv05.GL(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);
            case 36
                AA(m,:)=LibAv05.expgaussian(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);    
            case 37
                AA(m,:)=LibAv05.pearson(xxx,FitParameters(m),FIXEDPARAMETERS(m),extra);
            case 38
                AA(m,:)=LibAv05.explorentzian(xxx,FitParameters(3*m-2),FitParameters(3*m-1),-FitParameters(3*m)*length(xxx)./length(xx));        
            case 40
                AA(m,:)=LibAv05.sine(xx,FitParameters(2*m-1),FitParameters(2*m));
            case 41
                AA(m,:)=LibAv05.rectangle(xxx,FitParameters(2*m-1),FitParameters(2*m));
            case 42
                AA(m,:)=LibAv05.ngaussian(xxx,FitParameters(2*m-1),FitParameters(2*m),extra);
            case 43
                AA(m,:)=LibAv05.Gompertz(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));  
            case 44
                AA(m,:)=LibAv05.OneMinusExp(xxx,FitParameters(2*m-1),FitParameters(2*m));
            case 45
                AA(m,:)=LibAv05.FourPL(xxx,FitParameters(3*m-2),FitParameters(3*m-1),FitParameters(3*m));    
            case 46
                AA(m,:)=LibAv05.quadslope(xxx,FitParameters(2*m-1),FitParameters(2*m));
            case 47
                AA(m,:)=LibAv05.blackbody(xxx,FitParameters(m));

               otherwise
          end % switch
        end % for NumPeaks

        % Multiplies each row by the corresponding amplitude and adds them up
        if peakshape(1)==29, % Segmented linear
            mmodel=LibAv05.segmented(xx,yy,PEAKHEIGHTS);
            baseline=0;
        else
            heightsize=size(height');
            AAsize=size(AA);
            if heightsize(2)==AAsize(1),
                mmodel=height'*AA+baseline;
            else
                mmodel=height*AA+baseline;
            end
        end
        % Top half of the figure shows original signal and the fitted model.
        if plots,
            subplot(2,1,1);plot(xx+xoffset,yy,'b.'); % Plot the original signal in blue dots
            hold on
        end
        if peakshape(1)==28, % Polynomial
             yi=LibAv05.polynomial(xxx,coeff);
        else
            for m=1:NumPeaks,
                if plots, plot(xxx+xoffset,height(m)*AA(m,:)+baseline,'g'),end  % Plot the individual component peaks in green lines
                area(m)=trapz(xxx+xoffset,height(m)*AA(m,:)); % Compute the area of each component peak using trapezoidal method
                yi(m,:)=height(m)*AA(m,:); % Place y values of individual model peaks into matrix yi
            end
        end
        xi=xxx+xoffset; % Place the x-values of the individual model peaks into xi
        residual=yy-bestmodel;

        if plots,
            % Mark starting peak positions with vertical dashed magenta lines
            if peakshape(1)==16||peakshape(1)==17
            else
                if peakshape(1)==29, % Segmented linear
                    subplot(2,1,1);plot([PEAKHEIGHTS' PEAKHEIGHTS'],[0 max(yy)],'m--')
                else
                    for marker=1:NumPeaks,
                        markx=BestStart((2*marker)-1);
                        subplot(2,1,1);plot([markx+xoffset markx+xoffset],[0 max(yy)],'m--')
                    end % for
                end
            end % if peakshape

            % Plot the total model (sum of component peaks) in red lines
            if peakshape(1)==29, % Segmented linear
                mmodel=LibAv05.segmented(xx,yy,PEAKHEIGHTS);
               plot(xx+xoffset,mmodel,'r');  
            else
               plot(xxx+xoffset,mmodel,'r');  
            end
            hold off;
            lyy=min(yy);
            uyy=max(yy)+(max(yy)-min(yy))/10;
            if BIPOLAR,
                axis([min(xx) max(xx) lyy uyy]);
                ylabel('+ - mode')
            else
                axis([min(xx) max(xx) 0 uyy]);
                ylabel('+ mode')
            end
            switch AUTOZERO,
                case 0
                    title(['peakfit.m Version 8   No baseline correction'])
                case 1
                    title(['peakfit.m Version 8   Linear baseline subtraction'])
                case 2
                    title(['peakfit.m Version 8   Quadratic subtraction baseline'])
                case 3
                    title(['peakfit.m Version 8   Flat baseline correction'])
            end

            switch peakshape(1)
                case {4,20,34,37,42}
                    xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Shape Constant = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
                case {5,8,18,36}
                    xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      Time Constant = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
                case 13
                    xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      % Gaussian = ' num2str(extra)   '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000)  ] )
                case {14,15,22,35}
                    xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH) '      extra = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
                case 28
                    xlabel(['Shape = ' ShapeString '      Order = ' num2str(extra)  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(1000*LowestError)/1000) ] )
                case 43
                    xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '        Error = ' num2str(round(1000*LowestError)/1000) '%   R2 = ' num2str(round(100000*Rsquared)/100000) ] )
                 otherwise
                    if peakshape(1)==29, % Segmented linear
                        xlabel(['Breakpoints = ' num2str(NumPeaks) '     Shape = ' ShapeString  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
                    else
                        xlabel(['Peaks = ' num2str(NumPeaks) '     Shape = ' ShapeString '     Min. Width = ' num2str(MINWIDTH)  '     Error = ' num2str(round(1000*LowestError)/1000) '%  R2 = ' num2str(round(100000*Rsquared)/100000) ] )
                    end % if peakshape(1)==29
            end % switch peakshape(1)

            % Bottom half of the figure shows the residuals and displays RMS error
            % between original signal and model
            % residual=yy-bestmodel;
            subplot(2,1,2);plot(xx+xoffset,residual,'m.')
            axis([min(xx)+xoffset max(xx)+xoffset min(residual) max(residual)]);
            xlabel('Residual Plot')
            if NumTrials>1,
               title(['Best of ' num2str(NumTrials) ' fits'])
            else
               title(['Single fit'])
            end
        end % if plots

        % Put results into a matrix FitResults, one row for each peak, showing peak index number,
        % position, amplitude, and width.
        FitResults=zeros(NumPeaks,6);
        %  FitParameters=FitParameters
        switch peakshape(1),
            case {6,7,8}, % equal-width peak models only
                for m=1:NumPeaks,
                    if m==1,
                        FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
                    else
                        FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
                    end
                end
            case {11,12,34,35,36,37}, % Fixed-width shapes only
                for m=1:NumPeaks,
                    width(m)=abs(FitParameters(m));
                    if peakshape==34,
                        gD(m)=width(m);
                        gL(m)=extra.*gD(m);
                        width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
                    end
                    if m==1,
                        FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
                    else
                        FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
                    end
                end
            case {16,17}, % Fixed-position shapes only
                for m=1:NumPeaks,
                    if m==1,
                        FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
                    else
                        FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
                    end
                end
            case 28,   % Simple polynomial fit
                FitResults=PEAKHEIGHTS;
            case 29, % Segmented linear fit
                FitResults=PEAKHEIGHTS;
            case {30,31,32,33,38,43} % Special case of shapes with 3 iterated variables
                for m=1:NumPeaks,
                    width(m)=abs(FitParameters(3*m-1));
                    if peakshape==30,
                        gD(m)=width(m);
                        gL(m)=FitParameters(3*m).*gD(m);
                        width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
                    end
                    if m==1,
                        FitResults=[round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)];
                    else
                        FitResults=[FitResults ; [round(m) FitParameters(3*m-2) height(m) width(m) area(m)] FitParameters(3*m)];
                    end
                end
            case 47 % Shapes with 1 iterated variable
                FitParameters=FitParameters
            otherwise % Normal shapes with 2 iterated variables
                for m=1:NumPeaks,
                    width(m)=abs(FitParameters(2*m));
                    if peakshape==20,
                        gD=width(m);
                        gL=extra.*gD;
                        width(m) = 2.*(0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2));
                    end
                    if m==1,
                        FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)];
                    else
                        FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)]];
                    end % if m==1

                end % for m=1:NumPeaks,
        end % switch peakshape(1)

        % Rearrange fit results for Gompertz to Bo, Kh, and L
        if peakshape(1)==43;
            for m=1:NumPeaks,
                FitResults(m,2)=FitResults(m,2).*FitResults(m,3);
                FitResults(m,4)=FitResults(m,3).*FitResults(m,4);
                FitResults(m,3)=1;
            end
        end

        %Sort FitResults
        FitResults=sortrows(FitResults,2);

        % Display Fit Results on lower graph
        if plots,
            % Display Fit Results on lower  graph
            subplot(2,1,2);
            startx=min(xx)+(max(xx)-min(xx))./20;
            dxx=(max(xx)-min(xx))./10;
            dyy=((max(residual)-min(residual))./10);
            starty=max(residual)-dyy;
            FigureSize=get(gcf,'Position');
            switch peakshape(1)
                case {9,19,10,23,40}  % Pulse and sigmoid shapes only
                    text(startx,starty+dyy/2,['Peak #          tau1           Height           tau2             Area'] );
                case 28, % Polynomial
                    text(startx,starty+dyy/2,['Polynomial coefficients'] );
                case 29 % Segmented linear
                     text(startx,starty+dyy/2,['x-axis breakpoints'] );
                case {30,31,32,33,38} % Special case of shapes with 3 iterated variables
                    text(startx,starty+dyy/2,['Peak #          Position        Height         Width             Area       Shape factor'] );            
                case 43 % 3 parameter Gompertz
                    text(startx,starty+dyy/2,['Peak #           Bo             Height            Kh                Area                 L'] );            
                otherwise
                    text(startx,starty+dyy/2,['Peak #          Position         Height         Width             Area   '] );
            end
            % Display FitResults using sprintf
            if peakshape(1)==28||peakshape(1)==29, % Polynomial or segmented linear
                for number=1:length(FitResults),
                    column=1;
                    itemstring=sprintf('%0.4g',FitResults(number));
                    xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
                    yposition=starty-number.*dyy.*(400./FigureSize(4));
                    text(xposition,yposition,['                ' itemstring]);
                end
            else
                for peaknumber=1:NumPeaks,
                    for column=1:5,
                        itemstring=sprintf('%0.4g',FitResults(peaknumber,column));
                        xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
                        yposition=starty-peaknumber.*dyy.*(400./FigureSize(4));
                        text(xposition,yposition,itemstring);
                    end
                end
                xposition=startx;
                yposition=starty-(peaknumber+1).*dyy.*(400./FigureSize(4));
                if AUTOZERO==3,
                    text(xposition,yposition,[ 'Baseline= ' num2str(baseline) ]);
                end % if AUTOZERO
            end % if peakshape(1)
            if peakshape(1)==30 || peakshape(1)==31 || peakshape(1)==32 || peakshape(1)==33 || peakshape(1)==38 || peakshape(1)==43,
                for peaknumber=1:NumPeaks,
                    column=6;
                    itemstring=sprintf('%0.4g',FitParameters(3*peaknumber));
                    xposition=startx+(1.7.*dxx.*(column-1).*(600./FigureSize(3)));
                    yposition=starty-peaknumber.*dyy.*(400./FigureSize(4));
                    text(xposition,yposition,itemstring);
                end
            end
        end % if plots

        if NumArgOut==8,
            if plots,disp('Computing bootstrap sampling statistics.....'),end
            BootstrapResultsMatrix=zeros(6,100,NumPeaks);
            BootstrapErrorMatrix=zeros(1,100,NumPeaks);
            clear bx by
            tic;
            for trial=1:100,
                n=1;
                bx=xx;
                by=yy;
                while n<length(xx)-1,
                    if rand>.5,
                        bx(n)=xx(n+1);
                        by(n)=yy(n+1);
                    end
                    n=n+1;
                end
                bx=bx+xoffset;
                [FitResults,BootFitError]=LibAv05.fitpeaks(bx,by,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,FIXEDPARAMETERS,shapesvector);
                for peak=1:NumPeaks,
                    switch peakshape(1)
                        case {30,31,32,33,38,43}
                            BootstrapResultsMatrix(1:6,trial,peak)=FitResults(peak,1:6);
                        otherwise
                            BootstrapResultsMatrix(1:5,trial,peak)=FitResults(peak,1:5);
                    end
                    BootstrapErrorMatrix(:,trial,peak)=BootFitError;
                end
            end
            if plots,toc;end
            for peak=1:NumPeaks,
                if plots,
                    disp(' ')
                    % Label columns for bootstrap results
                    switch peakshape(1)
                        case {9,19,10,23,40}  % Pulse and sigmoid shapes only
                            disp(['Peak #',num2str(peak) '         tau1           Height         tau2           Area'] );
                        case {30,31,32,33,38} % Special case of shapes with 3 iterated variables
                            disp(['Peak #',num2str(peak) '         Position        Height         Width             Area       Shape factor'] );
                        case 43 % 3 parameter Gompertz
                            disp(['Peak #',num2str(peak) '         Bo             Height         Kh                       L'] );
                        otherwise
                            disp(['Peak #',num2str(peak) '         Position         Height         Width             Area'] );
                    end

                end % if plots
                BootstrapMean=mean(real(BootstrapResultsMatrix(:,:,peak)'));
                BootstrapSTD=std(BootstrapResultsMatrix(:,:,peak)');
                BootstrapLibAv05.iqr=LibAv05.iqr(BootstrapResultsMatrix(:,:,peak)');
                PercentRSD=100.*BootstrapSTD./BootstrapMean;
                PercentLibAv05.iqr=100.*BootstrapLibAv05.iqr./BootstrapMean;
                BootstrapMean=BootstrapMean(2:6);
                BootstrapSTD=BootstrapSTD(2:6);
                BootstrapLibAv05.iqr=BootstrapLibAv05.iqr(2:6);
                PercentRSD=PercentRSD(2:6);
                PercentLibAv05.iqr=PercentLibAv05.iqr(2:6);
                format short g
                if plots,
                    disp(['Bootstrap Mean: ', num2str(BootstrapMean)])
                    disp(['Bootstrap STD:  ', num2str(BootstrapSTD)])
                    disp(['Bootstrap LibAv05.iqr:  ', num2str(BootstrapLibAv05.iqr)])
                    disp(['Percent RSD:    ', num2str(PercentRSD)])
                    disp(['Percent LibAv05.iqr:    ', num2str(PercentLibAv05.iqr)])
                end % if plots
                BootResults(peak,:)=[BootstrapMean BootstrapSTD PercentRSD BootstrapLibAv05.iqr PercentLibAv05.iqr];
            end % peak=1:NumPeaks,
        end % if NumArgOut==8,
        if AUTOZERO==3;
        else
            baseline=bkgcoef;
        end
        
        end
        
        % ----------------------------------------------------------------------
        function [FitResults,LowestError]=fitpeaks(xx,yy,NumPeaks,peakshape,extra,NumTrials,start,AUTOZERO,fixedparameters,shapesvector)
        % Based on peakfit Version 3: June, 2012. 
        global PEAKHEIGHTS FIXEDPARAMETERS BIPOLAR MINWIDTH coeff
        format short g
        format compact
        warning off all
        FIXEDPARAMETERS=fixedparameters;
        xoffset=0;
        if start==0;start=LibAv05.calcstart(xx,NumPeaks,xoffset);end
        PEAKHEIGHTS=zeros(1,NumPeaks);
        n=length(xx);
        newstart=start;
        coeff=0;
        LOGPLOT=0;

        % Perform peak fitting for selected peak shape using fminsearch function
        options = optimset('TolX',.000001,'Display','off','MaxFunEvals',1000 );
        LowestError=1000; % or any big number greater than largest error expected
        FitParameters=zeros(1,NumPeaks.*2); 
        BestStart=zeros(1,NumPeaks.*2); 
        height=zeros(1,NumPeaks); 
        bestmodel=zeros(size(yy));
        for k=1:NumTrials,
            % StartVector=newstart
            switch peakshape(1)
                case 1
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitgaussian(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 2
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitlorentzian(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 3
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitlogistic(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 4
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitpearson(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 5
                    zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
                    zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexpgaussian(lambda,zxx,zyy,-extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 6
                    cwnewstart(1)=newstart(1);
                    for pc=2:NumPeaks,
                        cwnewstart(pc)=newstart(2.*pc-1);
                    end
                    cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitewgaussian(lambda,xx,yy)),cwnewstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(NumPeaks+1)<MINWIDTH,
                            TrialParameters(NumPeaks+1)=MINWIDTH;
                        end
                    end
                case 7
                    cwnewstart(1)=newstart(1);
                    for pc=2:NumPeaks,
                        cwnewstart(pc)=newstart(2.*pc-1);
                    end
                    cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitewlorentzian(lambda,xx,yy)),cwnewstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(NumPeaks+1)<MINWIDTH,
                            TrialParameters(NumPeaks+1)=MINWIDTH;
                        end
                    end
                case 8
                    cwnewstart(1)=newstart(1);
                    for pc=2:NumPeaks,
                        cwnewstart(pc)=newstart(2.*pc-1);
                    end
                    cwnewstart(NumPeaks+1)=(max(xx)-min(xx))/5;
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexpewgaussian(lambda,xx,yy,-extra)),cwnewstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(NumPeaks+1)<MINWIDTH,
                            TrialParameters(NumPeaks+1)=MINWIDTH;
                        end
                    end
                case 9
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexppulse(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 10
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitupsigmoid(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 11
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.FitFWGaussian(lambda,xx,yy)),fixedstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 12
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.FitFWLorentzian(lambda,xx,yy)),fixedstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 13
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitGL(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 14
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitBiGaussian(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 15
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitBWF(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 16
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.FitFPGaussian(lambda,xx,yy)),fixedstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(Peak)<MINWIDTH,
                            TrialParameters(Peak)=MINWIDTH;
                        end
                    end
                case 17
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=(max(xx)-min(xx))./(NumPeaks+1);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.FitFPLorentzian(lambda,xx,yy)),fixedstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(Peak)<MINWIDTH,
                            TrialParameters(Peak)=MINWIDTH;
                        end
                    end
                case 18
                    zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
                    zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexplorentzian(lambda,zxx,zyy,-extra)),newstart,options);
                case 19
                    TrialParameters=fminsearch(@(lambda)(LibAv05.alphafunction(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 20
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitvoigt(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 21
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fittriangular(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 22
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitmultiple(lambda,xx,yy,NumPeaks,shapesvector,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 23
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitdownsigmoid(lambda,xx,yy)),newstart,optionst);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 24
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitnbinpdf(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 25
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitlognpdf(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 26
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitlinslope(lambda,xx,yy)),polyfit(xx,yy,1),options);
                coeff=TrialParameters;
                case 27
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitd1gauss(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 28
                    TrialParameters=LibAv05.fitpolynomial(xx,yy,extra);
                case 29
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitsegmented(lambda,xx,yy)),newstart,options);
                case 30
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitvoigtv(lambda,xx,yy)),newstart);
                case 31
                    zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
                    zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexpgaussianv(lambda,zxx,zyy)),newstart);
                case 32
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitpearsonv(lambda,xx,yy)),newstart);
                case 33
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitGLv(lambda,xx,yy)),newstart);
                case 34
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitFWVoigt(lambda,xx,yy,extra)),fixedstart,options);
                case 35
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitFWGL(lambda,xx,yy,extra)),fixedstart,options);
                case 36
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitFWExpGaussian(lambda,xx,yy,extra)),fixedstart,options);
                case 37
                    fixedstart=[];
                    for pc=1:NumPeaks,
                        fixedstart(pc)=min(xx)+pc.*(max(xx)-min(xx))./(NumPeaks+1);
                    end
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitFWPearson(lambda,xx,yy,extra)),fixedstart,options);
                case 38
                    zxx=[zeros(size(xx)) xx zeros(size(xx)) ];
                    zyy=[zeros(size(yy)) yy zeros(size(yy)) ];
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitexplorentzianv(lambda,zxx,zyy)),newstart);         
                case 40
                     TrialParameters=fminsearch(@(lambda)(LibAv05.fitsine(lambda,xx,yy)),newstart,options);
                case 41
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitrectangle(lambda,xx,yy)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 42
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitngaussian(lambda,xx,yy,extra)),newstart,options);
                    for Peak=1:NumPeaks;
                        if TrialParameters(2*Peak)<MINWIDTH,
                            TrialParameters(2*Peak)=MINWIDTH;
                        end
                    end
                case 43
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitGompertz(lambda,xx,yy)),newstart);
                case 44
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitOneMinusExp(lambda,xx,yy)),newstart);           
                case 45
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitFourPL(lambda,xx,yy)),newstart);  
                case 46
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitquadslope(lambda,xx,yy)),newstart,options);
                case 47
                    bbstart=3000;
                    TrialParameters=fminsearch(@(lambda)(LibAv05.fitblackbody(lambda,xx,yy)),bbstart,options);
                otherwise
            end % switch peakshape

        for peaks=1:NumPeaks,
             peakindex=2*peaks-1;
             newstart(peakindex)=start(peakindex)-xoffset;
        end

            % Construct model from Trial parameters
            A=zeros(NumPeaks,n);
            for m=1:NumPeaks,
                switch peakshape(1)
                    case 1
                        A(m,:)=LibAv05.gaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 2
                        A(m,:)=LibAv05.lorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 3
                        A(m,:)=LibAv05.logistic(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 4
                        A(m,:)=LibAv05.pearson(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
                    case 5
                        A(m,:)=LibAv05.expgaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
                    case 6
                        A(m,:)=LibAv05.gaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
                    case 7
                        A(m,:)=LibAv05.lorentzian(xx,TrialParameters(m),TrialParameters(NumPeaks+1));
                    case 8
                        A(m,:)=LibAv05.expgaussian(xx,TrialParameters(m),TrialParameters(NumPeaks+1),-extra)';
                    case 9
                        A(m,:)=LibAv05.exppulse(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 10
                        A(m,:)=LibAv05.upsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 11
                        A(m,:)=LibAv05.gaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
                    case 12
                        A(m,:)=LibAv05.lorentzian(xx,TrialParameters(m),FIXEDPARAMETERS(m));
                    case 13
                        A(m,:)=LibAv05.GL(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
                    case 14
                        A(m,:)=LibAv05.BiGaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
                    case 15
                        A(m,:)=LibAv05.BWF(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
                    case 16
                        A(m,:)=LibAv05.gaussian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
                    case 17
                        A(m,:)=LibAv05.lorentzian(xx,FIXEDPOSITIONS(m),TrialParameters(m));
                    case 18
                        A(m,:)=LibAv05.explorentzian(xx,TrialParameters(2*m-1),TrialParameters(2*m),-extra)';
                    case 19
                        A(m,:)=LibAv05.alphafunction(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 20
                        A(m,:)=LibAv05.voigt(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
                    case 21
                        A(m,:)=LibAv05.triangular(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 22
                        A(m,:)=LibAv05.peakfunction(shapesvector(m),xx,TrialParameters(2*m-1),TrialParameters(2*m),extra(m));
                    case 23
                        A(m,:)=LibAv05.downsigmoid(xx,TrialParameters(2*m-1),TrialParameters(2*m));      
                    case 24
                        A(m,:)=LibAv05.nbinpdf(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 25
                        A(m,:)=LibAv05.lognormal(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 26
                        A(m,:)=LibAv05.linslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 27
                        A(m,:)=LibAv05.d1gauss(xx,TrialParameters(2*m-1),TrialParameters(2*m));       
                    case 28
                        A(m,:)=LibAv05.polynomial(xx,TrialParameters(2*m-1),TrialParameters(2*m));       
                    case 29
                        A(m,:)=LibAv05.segmented(xx,yy,PEAKHEIGHTS);
                    case 30
                        A(m,:)=LibAv05.voigt(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
                    case 31
                        A(m,:)=LibAv05.expgaussian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
                    case 32
                        A(m,:)=LibAv05.pearson(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
                    case 33
                        A(m,:)=LibAv05.GL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));
                    case 34
                        width(m)=abs(FIXEDPARAMETERS(m));

        %                 gD(m)=width(m);
        %                 gL(m)=extra.*gD(m);
        %                 width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2))

                        A(m,:)=LibAv05.voigt(xx,TrialParameters(m), width(m),extra);
                    case 35
                        A(m,:)=LibAv05.GL(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
                    case 36
                        A(m,:)=LibAv05.expgaussian(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
                    case 37
                        A(m,:)=LibAv05.pearson(xx,TrialParameters(m),FIXEDPARAMETERS(m),extra);
                    case 38
                        A(m,:)=LibAv05.explorentzian(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));                    
                    case 40
                        A(m,:)=LibAv05.sine(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 41
                        A(m,:)=LibAv05.rectangle(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 42
                        A(m,:)=LibAv05.ngaussian(xx,TrialParameters(2*m-1),TrialParameters(2*m),extra);
                    case 43
                        A(m,:)=LibAv05.Gompertz(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));
                    case 44
                        A(m,:)=LibAv05.OneMinusExp(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 45
                        A(m,:)=LibAv05.FourPL(xx,TrialParameters(3*m-2),TrialParameters(3*m-1),TrialParameters(3*m));        
                    case 46
                        A(m,:)=LibAv05.quadslope(xx,TrialParameters(2*m-1),TrialParameters(2*m));
                    case 47
                        A(m,:)=LibAv05.blackbody(xx,TrialParameters(m));

                end % switch
            end % for

            % Multiplies each row by the corresponding amplitude and adds them up
            if peakshape(1)==29, % Segmented linear
                model=LibAv05.segmented(xx,yy,PEAKHEIGHTS);
                TrialParameters=coeff;
                Heights=ones(size(coeff));
            else
                if AUTOZERO==3,
                    baseline=PEAKHEIGHTS(1);
                    Heights=PEAKHEIGHTS(2:1+NumPeaks);
                    model=Heights'*A+baseline;
                else
                    model=PEAKHEIGHTS'*A;
                    Heights=PEAKHEIGHTS;
                    baseline=0;
                end
            end

            % Compare trial model to data segment and compute the fit error
            MeanFitError=100*norm(yy-model)./(sqrt(n)*max(yy));
            % Take only the single fit that has the lowest MeanFitError
            if MeanFitError<LowestError,
                if min(Heights)>=-BIPOLAR*10^100,  % Consider only fits with positive peak heights
                    LowestError=MeanFitError;  % Assign LowestError to the lowest MeanFitError
                    FitParameters=TrialParameters;  % Assign FitParameters to the fit with the lowest MeanFitError
                    height=Heights; % Assign height to the PEAKHEIGHTS with the lowest MeanFitError
                end % if min(PEAKHEIGHTS)>0
            end % if MeanFitError<LowestError
        end % for k (NumTrials)
            Rsquared=1-(norm(yy-bestmodel)./norm(yy-mean(yy)));
            SStot=sum((yy-mean(yy)).^2);
            SSres=sum((yy-bestmodel).^2);
            Rsquared=1-(SSres./SStot);
            GOF=[LowestError Rsquared];
        for m=1:NumPeaks,
            area(m)=trapz(xx+xoffset,height(m)*A(m,:)); % Compute the area of each component peak using trapezoidal method
        end

        % Put results into a matrix FitResults, one row for each peak, showing peak index number,
        % position, amplitude, and width.
        FitResults=zeros(NumPeaks,6);
        switch peakshape(1),
            case {6,7,8}, % equal-width peak models only
                for m=1:NumPeaks,
                    if m==1,
                        FitResults=[[round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
                    else
                        FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) abs(FitParameters(NumPeaks+1)) area(m)]];
                    end
                end
            case {11,12,35,36,37}, % Fixed-width shapes only
                for m=1:NumPeaks,
                    if m==1,
                        FitResults=[[round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
                    else
                        FitResults=[FitResults ; [round(m) FitParameters(m)+xoffset height(m) FIXEDPARAMETERS(m) area(m)]];
                    end
                end

            case {16,17}, % Fixed-position shapes only
                for m=1:NumPeaks,
                    if m==1,
                        FitResults=[round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)];
                    else
                        FitResults=[FitResults ; [round(m) FIXEDPOSITIONS(m) height(m) FitParameters(m) area(m)]];
                    end
                end
            case 28,   % Simple polynomial fit
                FitResults=PEAKHEIGHTS;
            case 29, % Segmented linear fit
                FitResults=PEAKHEIGHTS;
            case {30,31,32,33,38,43} % Special case of shapes with 3 iterated variables
                for m=1:NumPeaks,
                    width(m)=abs(FitParameters(3*m-1));
                    if peakshape==30,
                        gD(m)=width(m);
                        gL(m)=FitParameters(3*m).*gD(m);
                        width(m) = 2.*(0.5346*gL(m) + sqrt(0.2166*gL(m).^2 + gD(m).^2));
                    end
                    if m==1,
                        FitResults=[round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)];
                    else
                        FitResults=[FitResults ; [round(m) FitParameters(3*m-2) height(m) width(m) area(m) FitParameters(3*m)]];
                    end
                end
            otherwise % Normal shapes with 2 iterated variables
                for m=1:NumPeaks,
                    width(m)=abs(FitParameters(2*m));
                    if peakshape==20,
                        gD=width(m);
                        gL=extra.*gD;
                        width(m) = 2.*(0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2));
                    end
                    if m==1,
                        FitResults=[round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)];
                    else
                        FitResults=[FitResults ; [round(m) FitParameters(2*m-1)+xoffset height(m) width(m) area(m)]];
                    end % if m==1
                end % for m=1:NumPeaks,
        end % switch peakshape(1)
        if peakshape==34,
                DW=2*(0.5346*a*1.2772 + sqrt(0.2166*a*1.2772.^2 + 1.2772.^2))
        end
        end
        
        % ----------------------------------------------------------------------
        function start=calcstart(xx,NumPeaks,xoffset)
          n=max(xx)-min(xx);
          start=[];
          startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
          for marker=1:NumPeaks,
              markx=startpos(marker)+ xoffset;
              start=[start markx n/ (3.*NumPeaks)];
          end % for marker
        end
        % ----------------------------------------------------------------------
        function [index,closestval]=val2ind(x,val)
        % Returns the index and the value of the element of vector x that is closest to val
        % If more than one element is equally close, returns vectors of indicies and values
        % Tom O'Haver (toh@umd.edu) October 2006
        % Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
        % [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
        dif=abs(x-val);
        index=find((dif-min(dif))==0);
        closestval=x(index);
        end
        % ----------------------------------------------------------------------
        function err = fitgaussian(lambda,t,y)
        % Fitting function for multiple Gaussian peaks.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        numpeaks=round(length(lambda)/2);
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
        %    if lambda(2*j)<MINWIDTH,lambda(2*j)=MINWIDTH;end
            A(:,j) = LibAv05.gaussian(t,lambda(2*j-1),lambda(2*j))';
        end 
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitngaussian(lambda,t,y,shapeconstant)
        %   Fitting functions for multiple flattened Gaussian peaks.
        % T. C. O'Haver (toh@umd.edu),   Version 1.3, October 23, 2006.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = LibAv05.ngaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        % ----------------------------------------------------------------------
        function err = fitewgaussian(lambda,t,y)
        % Fitting function for multiple Gaussian peaks with equal peak widths.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        numpeaks=round(length(lambda)-1);
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
            A(:,j) = LibAv05.gaussian(t,lambda(j),lambda(numpeaks+1))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        % ----------------------------------------------------------------------
        function err = FitFWGaussian(lambda,t,y)
        %	Fitting function for multiple fixed width Gaussians
        global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
        numpeaks=round(length(lambda));
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
            A(:,j) = LibAv05.gaussian(t,lambda(j),FIXEDPARAMETERS(j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        % ----------------------------------------------------------------------
        function err = FitFPGaussian(lambda,t,y)
        %	Fitting function for multiple fixed-position Gaussians
        global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
        numpeaks=round(length(lambda));
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
            A(:,j) = LibAv05.gaussian(t,FIXEDPARAMETERS(j), lambda(j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        % ----------------------------------------------------------------------
        function err = FitFPLorentzian(lambda,t,y)
        %	Fitting function for multiple fixed-position Lorentzians
        global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR
        numpeaks=round(length(lambda));
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
            A(:,j) = LibAv05.lorentzian(t,FIXEDPARAMETERS(j), lambda(j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        err = norm(z-y');
        end
        % ----------------------------------------------------------------------
        function err = FitFWLorentzian(lambda,t,y)
        %	Fitting function for multiple fixed width Lorentzians
        global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
        numpeaks=round(length(lambda));
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
            A(:,j) = LibAv05.lorentzian(t,lambda(j),FIXEDPARAMETERS(j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        % ----------------------------------------------------------------------
        function err = fitewlorentzian(lambda,t,y)
        % Fitting function for multiple  Lorentzian band signal with equal peak widths.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        numpeaks=round(length(lambda)-1);
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
            A(:,j) = LibAv05.lorentzian(t,lambda(j),lambda(numpeaks+1))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        % ----------------------------------------------------------------------
        function g = gaussian(x,pos,wid)
        %  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
        %  X may be scalar, vector, or matrix, pos and wid both scalar
        % Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
        % plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
        g = exp(-((x-pos)./(0.60056120439323.*wid)).^2);
        end
        % ----------------------------------------------------------------------
        function g = ngaussian(x,pos,wid,n)
        %  ngaussian(x,pos,wid) = flattened Gaussian centered on x=pos, half-width=wid
        %  x may be scalar, vector, or matrix, pos and wid both scalar
        % Shape is Gaussian when n=1. Becomes more rectangular as n increases.
        %  T. C. O'Haver, 1988, revised 2014
        % Example: ngaussian([1 2 3],1,2,1) gives result [1.0000    0.5000    0.0625]
        if n>0,
            g = 1-(10.^-(n.*LibAv05.gaussian(x,pos,wid)));
            g=g./max(g);
        else
            g = LibAv05.gaussian(x,pos,wid);
        end
        end
        % ----------------------------------------------------------------------
        function err = fitlorentzian(lambda,t,y)
        %	Fitting function for multiple lorentzians, lambda(1)=position, lambda(2)=width
        %	Fitgauss assumes a lorentzian function 
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = LibAv05.lorentzian(t,lambda(2*j-1),lambda(2*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        % [lambda PEAKHEIGHTS err]
        % ----------------------------------------------------------------------
        function g = lorentzian(x,position,width)
        % lorentzian(x,position,width) Lorentzian function.
        % where x may be scalar, vector, or matrix
        % position and width scalar
        % T. C. O'Haver, 1988
        % Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
        g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);
        end
        % ----------------------------------------------------------------------
        function err = fitlogistic(lambda,t,y)
        %	Fitting function for multiple logistic peaks, lambda(1)=position, lambda(2)=width
        %	between the data and the values computed by the current
        %	function of lambda.  Fitlogistic assumes a logistic function 
        %  T. C. O'Haver, May 2006
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = LibAv05.logistic(t,lambda(2*j-1),lambda(2*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g = logistic(x,pos,wid)
        % logistic function.  pos=position; wid=half-width (both scalar)
        % logistic(x,pos,wid), where x may be scalar, vector, or matrix
        % pos=position; wid=half-width (both scalar)
        % T. C. O'Haver, 1991 
        ni = exp(-((x-pos)/(.477.*wid)) .^2);
        g = (2.*ni)./(1+ni);
        end
        % ----------------------------------------------------------------------
        function err = fittriangular(lambda,t,y)
        %	Fitting function for multiple triangular, lambda(1)=position, lambda(2)=width
        %	between the data and the values computed by the current
        %	function of lambda.  Fittriangular assumes a triangular function 
        %  T. C. O'Haver, May 2006
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = LibAv05.triangular(t,lambda(2*j-1),lambda(2*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g = triangular(x,pos,wid)
        %Triangle function.  pos=position; wid=half-width (both scalar)
        %trianglar(x,pos,wid), where x may be scalar or vector,
        %pos=position; wid=half-width (both scalar)
        % T. C. O'Haver, 1991
        % Example
        % x=[0:.1:10];plot(x,trianglar(x,5.5,2.3),'.')
        g=1-(1./wid) .*abs(x-pos);
        for i=1:length(x),  
        if g(i)<0,g(i)=0;end
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitrectangle(lambda,t,y)
        %	Fitting function for multiple rectangle, lambda(1)=position, lambda(2)=width
        %	between the data and the values computed by the current
        %	function of lambda.  Fitrectangle assumes a rectangle function 
        %  T. C. O'Haver, May 2016
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = LibAv05.rectangle(t,lambda(2*j-1),lambda(2*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g = rectangle(x,pos,wid)
        %rectangle function.  pos=position; wid=half-width (both scalar)
        %rectangle(x,pos,wid), where x may be scalar or vector,
        %pos=position; wid=half-width (both scalar)
        % T. C. O'Haver, 2016
        % Example
        % x=[0:.1:10];plot(x,rectangle(x,5.5,2.3),'.')
        g=zeros(size(x));
        hw=wid./2;
        for i=1:length(x),  
        if x(i)<pos-hw,g(i)=0;end
        if x(i)>pos-hw,g(i)=1;end
        if x(i)>pos+hw,g(i)=0;end
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitpearson(lambda,t,y,shapeconstant)
        %   Fitting functions for multiple Pearson 7 bands.
        % T. C. O'Haver (toh@umd.edu),   Version 1.3, October 23, 2006.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = pearson(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitpearsonv(lambda,t,y)
        % Fitting functions for multiple pearson functions with independently variable
        % percent Gaussian
        % T. C. O'Haver (toh@umd.edu), 2015.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/3));
        for j = 1:length(lambda)/3,
            A(:,j) = LibAv05.pearson(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitFWPearson(lambda,t,y,shapeconstant)
        %	Fitting function for multiple fixed width Pearson7
        global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
        numpeaks=round(length(lambda));
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
            A(:,j) = LibAv05.pearson(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g = pearson(x,pos,wid,m)
        % Pearson VII function. 
        % g = pearson(x,pos,wid,m) where x may be scalar, vector, or matrix
        % pos=position; wid=half-width (both scalar)
        % m=some number
        %  T. C. O'Haver, 1990  
        g=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m;
        end
        
        % ----------------------------------------------------------------------
        function err = fitexpgaussian(lambda,t,y,timeconstant)
        %   Fitting functions for multiple exponentially-broadened Gaussian bands signal.
        %  T. C. O'Haver, October 23, 2006.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = LibAv05.expgaussian(t,lambda(2*j-1),lambda(2*j),timeconstant);
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitexplorentzian(lambda,t,y,timeconstant)
        %   Fitting functions for multiple exponentially-broadened lorentzian band signal.
        %  T. C. O'Haver, 2013.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = LibAv05.explorentzian(t,lambda(2*j-1),lambda(2*j),timeconstant);
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitexplorentzianv(lambda,t,y)
        % Fitting functions for multiple exponentially-broadened Gaussians with
        % independently variable time constants
        % T. C. O'Haver (toh@umd.edu), 2015.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/3));
        for j = 1:length(lambda)/3,
            A(:,j) = LibAv05.explorentzian(t,lambda(3*j-2),lambda(3*j-1),-lambda(3*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitexpewgaussian(lambda,t,y,timeconstant)
        % Fitting function for multiple exponentially-broadened Gaussian bands with equal peak widths.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        numpeaks=round(length(lambda)-1);
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
            A(:,j) = LibAv05.expgaussian(t,lambda(j),lambda(numpeaks+1),timeconstant);
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        
        end
        
        % ----------------------------------------------------------------------
        function err = fitexpgaussianv(lambda,t,y)
        % Fitting functions for multiple exponentially-broadened Gaussians with
        % independently variable time constants
        % T. C. O'Haver (toh@umd.edu), 2015.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/3));
        for j = 1:length(lambda)/3,
            A(:,j) = LibAv05.expgaussian(t,lambda(3*j-2),lambda(3*j-1),-lambda(3*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        
        end
        
        
        % ----------------------------------------------------------------------
        function err = fitFWExpGaussian(lambda,t,y,shapeconstant)
        %	Fitting function for multiple fixed width exponentially-broadened gaussian
        global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
        numpeaks=round(length(lambda));
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
            A(:,j) = LibAv05.expgaussian(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        
        end
        
        
        % ----------------------------------------------------------------------
        function g = expgaussian(x,pos,wid,timeconstant)
        %  Exponentially-broadened gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
        %  x may be scalar, vector, or matrix, pos and wid both scalar
        %  T. C. O'Haver, 2006
        g = exp(-((x-pos)./(0.60056120439323.*wid)) .^2);
        g = LibAv05.ExpBroaden(g',timeconstant);
        end
        
        % ----------------------------------------------------------------------
        function g = explorentzian(x,pos,wid,timeconstant)
        %  Exponentially-broadened lorentzian(x,pos,wid) = lorentzian peak centered on pos, half-width=wid
        %  x may be scalar, vector, or matrix, pos and wid both scalar
        %  T. C. O'Haver, 2013
        g = ones(size(x))./(1+((x-pos)./(0.5.*wid)).^2);
        g = LibAv05.ExpBroaden(g',timeconstant);
        end
        
        % ----------------------------------------------------------------------
        function yb = ExpBroaden(y,t)
        % ExpBroaden(y,t) zero pads y and convolutes result by an exponential decay
        % of time constant t by multiplying Fourier transforms and inverse
        % transforming the result.
        hly=round(length(y)./2);
        ey=[y(1).*ones(1,hly)';y;y(length(y)).*ones(1,hly)'];
        fy=fft(ey);
        a=exp(-(1:length(fy))./t);
        fa=fft(a);
        fy1=fy.*fa';
        ybz=real(ifft(fy1))./sum(a);
        yb=ybz(hly+2:length(ybz)-hly+1);
        end
        
        % ----------------------------------------------------------------------
        function err = fitexppulse(tau,x,y)
        % Iterative fit of the sum of exponential pulses
        % of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(x),round(length(tau)/2));
        for j = 1:length(tau)/2,
            A(:,j) = exppulse(x,tau(2*j-1),tau(2*j));
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g = exppulse(x,t1,t2)
        % Exponential pulse of the form 
        % g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
        e=(x-t1)./t2;
        p = 4*exp(-e).*(1-exp(-e));
        p=p .* (p>0);
        g = p';
        end
        
        % ----------------------------------------------------------------------
        function err = fitalphafunction(tau,x,y)
        % Iterative fit of the sum of alpha funciton
        % of the form Height.*exp(-tau1.*x).*(1-exp(-tau2.*x)))
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(x),round(length(tau)/2));
        for j = 1:length(tau)/2,
            A(:,j) = alphafunction(x,tau(2*j-1),tau(2*j));
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g = alphafunction(x,pos,spoint)
        % alpha function.  pos=position; wid=half-width (both scalar)
        % alphafunction(x,pos,wid), where x may be scalar, vector, or matrix
        % pos=position; wid=half-width (both scalar)
        % Taekyung Kwon, July 2013  
        g = (x-spoint)./pos.*exp(1-(x-spoint)./pos);
        for m=1:length(x);if g(m)<0;g(m)=0;end;end
        end
        
        % ----------------------------------------------------------------------
        function err = fitdownsigmoid(tau,x,y)
        % Fitting function for iterative fit to the sum of multiple
        % downward moving sigmiods 
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(x),round(length(tau)/2));
        for j = 1:length(tau)/2,
            A(:,j) = downsigmoid(x,tau(2*j-1),tau(2*j));
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitupsigmoid(tau,x,y)
        % Fitting function for iterative fit to the sum of multiple
        % upwards moving sigmiods
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(x),round(length(tau)/2));
        for j = 1:length(tau)/2,
            A(:,j) = upsigmoid(x,tau(2*j-1),tau(2*j));
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g=downsigmoid(x,t1,t2)
         % down step sigmoid
        g=.5-.5*erf(real((x-t1)/sqrt(2*t2)));
        end
        
        % ----------------------------------------------------------------------
        function g=upsigmoid(x,t1,t2)
        % up step sigmoid
        g=1/2 + 1/2* erf(real((x-t1)/sqrt(2*t2))); 
        end
        
        % ----------------------------------------------------------------------
        function err = fitGL(lambda,t,y,shapeconstant)
        %   Fitting functions for multiple Gaussian/Lorentzian blend peaks
        % T. C. O'Haver (toh@umd.edu), 2012.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = GL(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitFWGL(lambda,t,y,shapeconstant)
        %	Fitting function for a multiple fixed width Gaussian/Lorentzian blend
        global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
        numpeaks=round(length(lambda));
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
            A(:,j) = GL(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        
        % ----------------------------------------------------------------------
        function err = fitGLv(lambda,t,y)
        % Fitting functions for multiple Gaussian/Lorentzian blend functions with
        % independently variable percent Gaussian
        % T. C. O'Haver (toh@umd.edu), 2015.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/3));
        for j = 1:length(lambda)/3,
            A(:,j) = GL(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g = GL(x,pos,wid,m)
        % Gaussian/Lorentzian blend. m = percent Gaussian character
        % pos=position; wid=half-width
        % m = percent Gaussian character.
        %  T. C. O'Haver, 2012
        % sizex=size(x)
        % sizepos=size(pos)
        % sizewid=size(wid)
        % sizem=size(m)
        g=2.*((m/100).*LibAv05.gaussian(x,pos,wid)+(1-(m(1)/100)).*LibAv05.lorentzian(x,pos,wid))/2;
        end
        
        % ----------------------------------------------------------------------
        function err = fitvoigt(lambda,t,y,shapeconstant)
        % Fitting functions for multiple Voigt profile function
        % T. C. O'Haver (toh@umd.edu), 2013.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = LibAv05.voigt(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitFWVoigt(lambda,t,y,shapeconstant)
        %	Fitting function for multiple fixed width Voigt
        global PEAKHEIGHTS AUTOZERO FIXEDPARAMETERS BIPOLAR LOGPLOT
        numpeaks=round(length(lambda));
        % numpeaksfitFWVoigt=numpeaks
        A = zeros(length(t),numpeaks);
        for j = 1:numpeaks,
            A(:,j) = voigt(t,lambda(j),FIXEDPARAMETERS(j),shapeconstant)';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitvoigtv(lambda,t,y)
        % Fitting functions for multiple Voigt profile function with independently variable
        % alphas
        % T. C. O'Haver (toh@umd.edu), 2015.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/3));
        for j = 1:length(lambda)/3,
            A(:,j) = voigt(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g=voigt(xx,pos,gD,alpha)
        % Voigt profile function. xx is the independent variable (energy,
        % wavelength, etc), gD is the Doppler (Gaussian) width, and alpha is the
        % shape constant (ratio of the Lorentzian width gL to the Doppler width gD.
        % Based on Chong Tao's "Voigt lineshape spectrum simulation", 
        % File ID: #26707
        % alpha=alpha
        gL=alpha.*gD;
        gV = 0.5346*gL + sqrt(0.2166*gL.^2 + gD.^2);
        x = gL/gV;
        % sizeabs=size(abs(xx-pos))
        % sizegV=size(gV)
        y = abs(xx-pos)./gV;
        g = 1/(2*gV*(1.065 + 0.447*x + 0.058*x^2))*((1-x)*exp(-0.693.*y.^2) + (x./(1+y.^2)) + 0.016*(1-x)*x*(exp(-0.0841.*y.^2.25)-1./(1 + 0.021.*y.^2.25)));
        g=g./max(g);
        end
        
        % ----------------------------------------------------------------------
        function err = fitBiGaussian(lambda,t,y,shapeconstant)
        %   Fitting functions for multiple BiGaussian.
        % T. C. O'Haver (toh@umd.edu),  2012.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = LibAv05.BiGaussian(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g = BiGaussian(x,pos,wid,m)
        % BiGaussian (different widths on leading edge and trailing edge).
        % pos=position; wid=width 
        % m determines shape; symmetrical if m=50.
        %  T. C. O'Haver, 2012
        lx=length(x);
        hx=LibAv05.val2ind(x,pos);
        g(1:hx)=LibAv05.gaussian(x(1:hx),pos,wid*(m/100));
        g(hx+1:lx)=LibAv05.gaussian(x(hx+1:lx),pos,wid*(1-m/100));
        end
        
        % ----------------------------------------------------------------------
        function err = fitBWF(lambda,t,y,shapeconstant)
        %   Fitting function for multiple Breit-Wigner-Fano.
        % T. C. O'Haver (toh@umd.edu),  2014.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = BWF(t,lambda(2*j-1),lambda(2*j),shapeconstant)';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g = BWF(x,pos,wid,m)
        % BWF (Breit-Wigner-Fano) http://en.wikipedia.org/wiki/Fano_resonance
        % pos=position; wid=width; m=Fano factor
        %  T. C. O'Haver, 2014
        y=((m*wid/2+x-pos).^2)./(((wid/2).^2)+(x-pos).^2);
        % y=((1+(x-pos./(m.*wid))).^2)./(1+((x-pos)./wid).^2);
        g=y./max(y);
        end
        
        % ----------------------------------------------------------------------
        function err = fitnbinpdf(tau,x,y)
        % Fitting function for iterative fit to the sum of multiple
        % Negative Binomial Distributions
        % (http://www.mathworks.com/help/stats/negative-binomial-distribution.html)
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(x),round(length(tau)/2));
        for j = 1:length(tau)/2,
            A(:,j) = nbinpdf(x,tau(2*j-1),tau(2*j));
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function err = fitlognpdf(tau,x,y)
        % Fitting function for iterative fit to the sum of multiple
        % Lognormal Distributions
        % (http://www.mathworks.com/help/stats/lognormal-distribution.html)
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(x),round(length(tau)/2));
        for j = 1:length(tau)/2,
            A(:,j) = LibAv05.lognormal(x,tau(2*j-1),tau(2*j));
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g = lognormal(x,pos,wid)
        % lognormal function.  pos=position; wid=half-width (both scalar)
        % lognormal(x,pos,wid), where x may be scalar, vector, or matrix
        % pos=position; wid=half-width (both scalar)
        % T. C. O'Haver, 1991  
        g = exp(-(log(x/pos)/(0.01.*wid)) .^2);
        end
        
        % ----------------------------------------------------------------------
        function err = fitsine(tau,x,y)
        % Fitting function for iterative fit to the sum of multiple
        % sine waves (alpha test, NRFPT)
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(x),round(length(tau)/2));
        for j = 1:length(tau)/2,
            A(:,j) = sine(x,tau(2*j-1),tau(2*j));
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        
        % ----------------------------------------------------------------------
        function g=sine(t,f,phase) 
        % Sine wave (alpha test)
        g=sin(2*pi*f*(t+phase));
        
        end
        % ----------------------------------------------------------------------
        function err = fitd1gauss(lambda,t,y)
        %   Fitting functions for multiple first derivative of a Gaussian
        %  T. C. O'Haver, 2014
        global PEAKHEIGHTS AUTOZERO BIPOLAR
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = d1gauss(t,lambda(2*j-1),lambda(2*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        err = norm(z-y');
        end
        % ----------------------------------------------------------------------
        function y=d1gauss(x,p,w)
        % First derivative of Gaussian (alpha test)
        y=-(5.54518.*(x-p).*exp(-(2.77259.*(p-x).^2)./w^2))./w.^2;
        y=y./max(y);
        end
        % ----------------------------------------------------------------------
        function coeff = fitpolynomial(t,y,order)
        coeff=polyfit(t,y,order);
        end
        % order=order
        % coeff=coeff
        % ----------------------------------------------------------------------
        function y=polynomial(t,coeff)
        y=polyval(coeff,t);
        end
        % ----------------------------------------------------------------------
        function err = fitsegmented(lambda,t,y)
        %   Fitting functions for articulated segmented linear
        %  T. C. O'Haver, 2014
        global LOGPLOT
        breakpoints=[t(1) lambda max(t)];
        z = segmented(t,y,breakpoints);
        % lengthz=length(z);
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y);
        end
        end
        % ----------------------------------------------------------------------
        function yi=segmented(x,y,segs)
        global PEAKHEIGHTS
        clear yy
        for n=1:length(segs)
          yind=LibAv05.val2ind(x,segs(n));
          yy(n)=y(yind(1));
        end
        yi=interp1(segs,yy,x);
        PEAKHEIGHTS=segs;
        end
        % ----------------------------------------------------------------------
        function err = fitlinslope(tau,x,y)
        % Fitting function for iterative fit to linear function
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(x),round(length(tau)/2));
        for j = 1:length(tau)/2,
            z = (x.*tau(2*j-1)+tau(2*j))';
            A(:,j) = z./max(z);
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        % ----------------------------------------------------------------------
        function y=linslope(x,slope,intercept)
        y=x.*slope+intercept;
        % y=y./max(y);
        end
        % ----------------------------------------------------------------------
        function err = fitquadslope(lambda,x,y)
        % Fitting function for iterative fit to linear function
        global PEAKHEIGHTS LOGPLOT
        A = zeros(length(x),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
             A(:,j)=quadslope(x,lambda(1),lambda(2));
        end
        PEAKHEIGHTS=A\y';
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        % ----------------------------------------------------------------------
        function y=quadslope(x,boa,coa) % normalized quadratic
         y=(x.^2+(boa).*x+coa);
        y=y./max(y);
        end
        % ----------------------------------------------------------------------
        function err = fitGompertz(lambda,t,y)
        % Fitting functions for multiple Gompertz function
        % T. C. O'Haver (toh@umd.edu), 2016.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/3));
        % Sizelambda=size(lambda)
        for j = 1:length(lambda)/3,
            A(:,j) = Gompertz(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j));
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        % PEAKHEIGHTS=1;
        % SizeA=size(A)
        % sizwPEAKHEIGHTS=size(PEAKHEIGHTS)
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        % ----------------------------------------------------------------------
        function y=Gompertz(t,Bo,Kh,L)
        % A LibAv05.Gompertz curve or Gompertz function, named after Benjamin Gompertz, is
        % a sigmoid function. It is a type of mathematical model for a time series,
        % where growth is slowest at the start and end of a time period. The
        % right-hand or future value asymptote of the function is approached much
        % more gradually by the curve than the left-hand or lower valued asymptote,
        % in contrast to the simple logistic function in which both asymptotes are
        % approached by the curve symmetrically. It is a special case of the
        % generalized logistic function.
        % 
        % Example:
        % x=1:.1:10;y=LibAv05.gompertz(x,6,3,4);plot(x,y)
        y=Bo*exp(-exp((Kh*exp(1)/Bo)*(L-t) +1));
        end
        % ----------------------------------------------------------------------
        function err = fitFourPL(lambda,t,y)
        % Fitting functions for multiple Gompertz function
        % T. C. O'Haver (toh@umd.edu), 2016.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT
        A = zeros(length(t),round(length(lambda)/3));
        for j = 1:length(lambda)/3,
            A(:,j) = FourPL(t,lambda(3*j-2),lambda(3*j-1),lambda(3*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        PEAKHEIGHTS=1;
        z = A*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        end
        % ----------------------------------------------------------------------
        function y=FourPL(x,miny,slope,ip)
        % Normalized four parameters logistic regression  
        % https://psg.hitachi-solutions.com/masterplex/blog/the-4-parameter-logisti
        % c-4pl-nonlinear-regression-model
        % miny = minimum asymptote. In an ELISA assay where you have a standard
        % curve, this can be thought of as the response value at 0 standard
        % concentration. slope = Hill slope. The Hill Slope or slope factor refers
        % to the steepness of the curve. It could either be positive or negative.
        % As the absolute value of the Hill slope increases, so does the steepness
        % of the curve. ip = inflection point: The inflection point is defined as the
        % point on the curve where the curvature changes direction or signs. This
        % can be better explained if you can imagine the concavity of a sigmoidal
        % curve. The inflection point is where the curve changes from being concave
        % upwards to concave downwards. maxy = maximum asymptote. In an ELISA assay
        % where you have a standard curve, this can be thought of as the response
        % value for infinite standard concentration.
        % Example:
        % x=0:20;
        % miny=0;slope=5;ip=10;d=0;maxy=10;
        % y=LibAv05.FourPL(x,miny,slope,ip,maxy);plot(x,y)
        %
        y = 1+(miny-1)./(1+(x./ip).^slope);
        end
        % ----------------------------------------------------------------------
        function err = fitOneMinusExp(lambda,t,y)
        %   Fitting functions for multiple first derivative of a Gaussian
        %  T. C. O'Haver, 2014
        global PEAKHEIGHTS AUTOZERO BIPOLAR
        A = zeros(length(t),round(length(lambda)/2));
        for j = 1:length(lambda)/2,
            A(:,j) = LibAv05.OneMinusExp(t,lambda(2*j-1),lambda(2*j))';
        end
        if AUTOZERO==3,A=[ones(size(y))' A];end
        if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        z = A*PEAKHEIGHTS;
        err = norm(z-y');
        end
        % ----------------------------------------------------------------------
        function g = OneMinusExp(x,pos,wid)
        % LibAv05.OneMinusExp(x,pos,wid) = 1-exp(-wid.*(x-pos));
        % Example:  x=0:10;y=1-exp(-(.5.*x));plot(x,y)
        g = 1-exp(-wid.*(x-pos));
        end
        % ----------------------------------------------------------------------
        
        function err = fitblackbody(lambda,wavelength,y)
        %  Fitting function for a LibAv05.blackbody spectrum.
        %  T. C. O'Haver, May 2008
        global PEAKHEIGHTS AUTOZERO BIPOLAR
        % sizelambda=size(lambda)
        radiance = LibAv05.blackbody(wavelength,lambda(1));
        % if AUTOZERO==3,A=[ones(size(y))' A];end
        % if BIPOLAR,PEAKHEIGHTS=A\y';else PEAKHEIGHTS=abs(A\y');end
        % sizePEAKHEIGHTS=size(PEAKHEIGHTS)
        % sizey=size(y)
        PEAKHEIGHTS = radiance/y;
        z = radiance*PEAKHEIGHTS;
        err = norm(z-y);
        end
        
        % ----------------------------------------------------------------------
        function radiance=blackbody(wavelength,temperature)
        radiance = 1.19111E+16*wavelength.^(-5)./(exp(14380000./(wavelength*temperature))-1);
        end
        % ----------------------------------------------------------------------
        function b=iqr(a)
        % b = LibAv05.iqr(a)  returns the interquartile range of the values in a.  For
        %  vector input, b is the difference between the 75th and 25th percentiles
        %  of a.  For matrix input, b is a row vector containing the interquartile
        %  range of each column of a.
        %  T. C. O'Haver, 2012
        mina=min(a);
        sizea=size(a);
        NumCols=sizea(2);
        for n=1:NumCols,b(:,n)=a(:,n)-mina(n);end
        Sorteda=sort(b);
        lx=length(Sorteda);
        SecondQuartile=round(lx/4);
        FourthQuartile=3*round(lx/4);
        b=abs(Sorteda(FourthQuartile,:)-Sorteda(SecondQuartile,:));
        
        end
        %-------------------------------------------------------------
        % LibAv05.fitmultiple
        %------------------------------------------------------------- 
        % ----------------------------------------------------------------------
        function err = fitmultiple(lambda,xx,y,numpeaks,shapesvector,extra)
        % Fitting function for a multiple-shape band signal.
        % The sequence of peak shapes are defined by the vector "shape".
        % The vector "extra" determines the shape of variable-shape peaks.
        global PEAKHEIGHTS AUTOZERO BIPOLAR LOGPLOT coeff FIXEDPARAMETERS
        % FIXEDPARAMETERS=FIXEDPARAMETERS % testing
        % PEAKHEIGHTS=PEAKHEIGHTS % testing
        % lambda=lambda
        % numpeaks=round(length(lambda)/2)
        % sizeshapesvector=size(shapesvector)
        A=zeros(numpeaks,length(xx));
        for j = 1:numpeaks,
            if shapesvector(j)==28,
                coeff=polyfit(xx,y,extra(j));
                A(j,:) = polyval(coeff,xx);
            else
                % sizeA=size(A)
                switch shapesvector(j)
                    case 1
                        A(j,:)=LibAv05.gaussian(xx,lambda(2*j-1),lambda(2*j));
                    case 2
                        A(j,:)=LibAv05.lorentzian(xx,lambda(2*j-1),lambda(2*j));
                    case 3
                        A(j,:)=LibAv05.logistic(xx,lambda(2*j-1),lambda(2*j));
                    case 4
                        A(j,:)=LibAv05.pearson(xx,lambda(2*j-1),lambda(2*j),extra(j));
                    case 5
                        A(j,:)=LibAv05.expgaussian(xx,lambda(2*j-1),lambda(2*j),-extra(j))';
                    case 6
                        A(j,:)=LibAv05.gaussian(xx,lambda(j),lambda(NumPeaks+1));
                    case 7
                        A(j,:)=LibAv05.lorentzian(xx,lambda(j),lambda(NumPeaks+1));
                    case 8
                        A(j,:)=LibAv05.expgaussian(xx,lambda(j),lambda(NumPeaks+1),-extra(j))';
                    case 9
                        A(j,:)=LibAv05.exppulse(xx,lambda(2*j-1),lambda(2*j));
                    case 10
                        A(j,:)=LibAv05.upsigmoid(xx,lambda(2*j-1),lambda(2*j));
                    case 11
                        A(j,:)=LibAv05.gaussian(xx,lambda(j),FIXEDPARAMETERS(j));
                    case 12
                        A(j,:)=LibAv05.lorentzian(xx,lambda(j),FIXEDPARAMETERS(j));
                    case 13
                        A(j,:)=LibAv05.GL(xx,lambda(2*j-1),lambda(2*j),extra(j));
                    case 14
                        A(j,:)=LibAv05.BiGaussian(xx,lambda(2*j-1),lambda(2*j),extra(j));
                    case 15
                        A(j,:)=LibAv05.BWF(xx,lambda(2*j-1),lambda(2*j),extra(j));
                    case 16
                        A(j,:)=LibAv05.gaussian(xx,FIXEDPOSITIONS(j),lambda(j));
                    case 17
                        A(j,:)=LibAv05.lorentzian(xx,FIXEDPOSITIONS(j),lambda(j));
                    case 18
                        A(j,:)=LibAv05.explorentzian(xx,lambda(2*j-1),lambda(2*j),-extra(j))';
                    case 19
                        A(j,:)=LibAv05.alphafunction(xx,lambda(2*j-1),lambda(2*j));
                    case 20
                        A(j,:)=LibAv05.voigt(xx,lambda(2*j-1),lambda(2*j),extra(j));
                    case 21
                        A(j,:)=LibAv05.triangular(xx,lambda(2*j-1),lambda(2*j));
                    case 22
                        A(j,:)=LibAv05.peakfunction(shapesvector(j),xx,lambda(2*j-1),lambda(2*j),extra(j));
                    case 23
                        A(j,:)=LibAv05.downsigmoid(xx,lambda(2*j-1),lambda(2*j));      
                    case 24
                        A(j,:)=LibAv05.nbinpdf(xx,lambda(2*j-1),lambda(2*j));
                    case 25
                        A(j,:)=LibAv05.lognormal(xx,lambda(2*j-1),lambda(2*j));
                    case 26
                        A(j,:)=LibAv05.linslope(xx,lambda(2*j-1),lambda(2*j));
                    case 27
                        A(j,:)=LibAv05.d1gauss(xx,lambda(2*j-1),lambda(2*j));       
                    case 28
                        A(j,:)=LibAv05.polynomial(xx,lambda(2*j-1),lambda(2*j));       
                    case 29
                        A(j,:)=LibAv05.segmented(xx,yy,PEAKHEIGHTS);
                    case 30
                        A(j,:)=LibAv05.voigt(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
                    case 31
                        A(j,:)=LibAv05.expgaussian(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
                    case 32
                        A(j,:)=LibAv05.pearson(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
                    case 33
                        A(j,:)=LibAv05.GL(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));
                    case 34
                        width(j)=abs(FIXEDPARAMETERS(j));
                         % [j lambda(j) width(j) extra(j)]
        %                 gD(j)=width(j);
        %                 gL(j)=extra.*gD(j);
        %                 width(j) = 2.*(0.5346*gL(j) + sqrt(0.2166*gL(j).^2 +
        %                 gD(j).^2))
                        A(j,:)=LibAv05.voigt(xx,lambda(j),width(j),extra(j));
                        % figure(5);plot(A(j,:));figure(1)
                    case 35
                        A(j,:)=LibAv05.GL(xx,lambda(j),FIXEDPARAMETERS(j),extra(j));
                    case 36
                        A(j,:)=LibAv05.expgaussian(xx,lambda(j),FIXEDPARAMETERS(j),extra(j));
                    case 37
                        A(j,:)=LibAv05.pearson(xx,lambda(j),FIXEDPARAMETERS(j),extra(j));
                    case 38
                        A(j,:)=LibAv05.explorentzian(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));                    
                    case 40
                        A(j,:)=LibAv05.sine(xx,lambda(2*j-1),lambda(2*j));
                    case 41
                        A(j,:)=LibAv05.rectangle(xx,lambda(2*    j-1),lambda(2*j));
                    case 42
                        A(j,:)=LibAv05.ngaussian(xx,lambda(2*j-1),lambda(2*j),extra(j));
                    case 43
                        A(j,:)=LibAv05.Gompertz(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));
                    case 44
                        A(j,:)=LibAv05.OneMinusExp(xx,lambda(2*j-1),lambda(2*j));
                    case 45
                        A(j,:)=LibAv05.FourPL(xx,lambda(3*j-2),lambda(3*j-1),lambda(3*j));        
                    case 46
                        A(j,:)=LibAv05.quadslope(xx,lambda(2*j-1),lambda(2*j));
                    case 47
                        A(j,:)=LibAv05.blackbody(xx,lambda(j));
                end % switch
            end % if shapesvector
        end % for j=1:numpeaks,
        if AUTOZERO==3,A=[ones(size(y))' A];end
        % sizeA=size(A)
        % max(A')
        % sizeyp=size(y')
        if BIPOLAR,PEAKHEIGHTS=A'\y';else PEAKHEIGHTS=abs(A'\y');end
        % PEAKHEIGHTS=PEAKHEIGHTS
        z = A'*PEAKHEIGHTS;
        if LOGPLOT,
            err = norm(log10(z)-log10(y)');
        else
            err = norm(z-y');
        end
        
        end
        % ----------------------------------------------------------------------
        %-------------------------------------------------------------
        % Rutinas LibAv05.peakfunction
        %------------------------------------------------------------- 
        function p=peakfunction(shape,x,pos,wid,extra,coeff)
        global FIXEDPARAMETERS
        % function that generates any of 20 peak types specified by number. 'shape'
        % specifies the shape type of each peak in the signal: "peakshape" = 1-20.
        % 1=Gaussian 2=Lorentzian, 3=logistic, 4=Pearson, 5=exponentionally
        % broadened Gaussian; 9=exponential pulse, 10=up sigmoid,
        % 13=Gaussian/Lorentzian blend; 14=BiGaussian, 15=Breit-Wigner-Fano (BWF) ,
        % 18=exponentionally broadened Lorentzian; 19=alpha function; 20=Voigt
        % profile; 21=triangular; 23=down sigmoid; 25=lognormal. "extra" is required
        % for variable-shape peaks only.
        switch shape,
            case 1
                p=LibAv05.gaussian(x,pos,wid);
            case 2
                p=LibAv05.lorentzian(x,pos,wid);
            case 3
                p=LibAv05.logistic(x,pos,wid);
            case 4
                p=LibAv05.pearson(x,pos,wid,extra);
            case 5
                p=LibAv05.expgaussian(x,pos,wid,extra);
            case 6
                p=LibAv05.gaussian(x,pos,wid);
            case 7
                p=LibAv05.lorentzian(x,pos,wid);
            case 8
                p=LibAv05.expgaussian(x,pos,wid,extra)';
            case 9
                p=LibAv05.exppulse(x,pos,wid);
            case 10
                p=LibAv05.upsigmoid(x,pos,wid);
            case 11
                p=LibAv05.gaussian(x,pos,wid);
            case 12
                p=LibAv05.lorentzian(x,pos,wid);
            case 13
                p=LibAv05.GL(x,pos,wid,extra);
            case 14
                p=LibAv05.BiGaussian(x,pos,wid,extra);
            case 15
                p=LibAv05.BWF(x,pos,wid,extra);
            case 16
                p=LibAv05.gaussian(x,pos,wid);
            case 17
                p=LibAv05.lorentzian(x,pos,wid);
            case 18
                p=LibAv05.explorentzian(x,pos,wid,extra)';
            case 19
                p=LibAv05.alphafunction(x,pos,wid);
            case 20
                p=LibAv05.voigt(x,pos,wid,extra);
            case 21
                p=LibAv05.triangular(x,pos,wid);
            case 23
                p=LibAv05.downsigmoid(x,pos,wid);
            case 25
                p=LibAv05.lognormal(x,pos,wid);
            case 26
                p=LibAv05.linslope(x,pos,wid);
            case 27
                p=LibAv05.d1gauss(x,pos,wid);
            case 28
                p=LibAv05.polynomial(x,coeff);
            case 30
                p=LibAv05.voigt(x,pos,pos,wid,extra);
            case 31
                p=LibAv05.expgaussian(x,pos,wid,-extra);
            case 32
                p=LibAv05.pearson(x,pos,wid,extra);
            case 33
                p=LibAv05.GL(x,pos,wid,extra);
            case 34
                % width(extra)=abs(FIXEDPARAMETERS(extra));
                p=LibAv05.voigt(x,pos,wid,extra);
            case 35
                p=LibAv05.GL(x,pos,wid,extra);
            case 36
                p=LibAv05.expgaussian(x,pos,wid,-extra);
            case 37
                p=LibAv05.pearson(x,pos,wid,extra);
            case 38
                p=LibAv05.explorentzian(x,pos,wid,-extra);
            case 40
                p=LibAv05.sine(x,pos,wid);
            case 41
                p=LibAv05.rectangle(x,pos,wid);
            case 42
                p=LibAv05.ngaussian(x,pos,wid,extra);
            case 43
                p=LibAv05.Gompertz(x,pos,wid,extra);
            case 44
                p=LibAv05.OneMinusExp(x,pos,wid);
            case 45
                p=LibAv05.FourPL(x,pos,wid,extra);
            case 46
                p=LibAv05.quadslope(x,pos,wid);
            case 47
                p=LibAv05.blackbody(x,pos);
            otherwise
        end % switch
        
        end
        
        %-------------------------------------------------------------
        % FIN Rutinas peakfit
        %-------------------------------------------------------------
        
        
    %----------------------    
    end % fin de methods
    %---------------------- 
%---------------------- 
end % fin de Class def
%---------------------- 