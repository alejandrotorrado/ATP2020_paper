% 21-Jul-2014 22:38:46 Hago que funcione en mac y linux (/ en vez de \)
% 15-Jul-2014 18:09:12 Separo las trayectorias en trajectories y
% trajectories_nogaps
% 29-Apr-2014 08:59:03 Elimino la encriptaci�n
% 17-Mar-2014 11:19:29 Limpieza y reorganizaci�n general. Cambio
% getComputerName por leer directamente el hostname
% 13-Mar-2014 10:26:49 Hago que guarde m�s a menudo datosegm (para que se
% guarden los tiempos)
% 15-Feb-2014 13:54:51 A�ado la resegmentaci�n
% 07-Feb-2014 09:34:22 Hago que guarde los tiempos
% 04-Feb-2014 15:13:59 Hago que pueda quedarse justo despu�s de calcular
% referencias
% 02-Feb-2014 20:02:11 A�ado el c�lculo de ramificaciones para tener en
% cuenta los trozos a los que ha podido ir un pez en cada momento. REVIERTO
% PARCIALMENTE: Lo dejo preparado, pero desactivado. Porque tampoco
% funciona demasiado bien, y me da miedo que sea m�s impredecible.
% 08-Jan-2014 11:15:14 Hago que abra parpool aunque s�lo se quiera 1
% procesador. Adem�s no abre parpool si ya est� abierto con el n�mero de
% procesadores adecuado.
% 04-Jan-2014 21:32:13 A�ado el rec�lculo de tam_mapas si es necesario
% 22-Dec-2013 02:23:38 Hago que vaya m�s r�pido si hay un solo bicho
% 21-Dec-2013 12:48:19 Cambio raizarchivos por nombrearchivo, para poder
% incluir la extensi�n
% 21-Dec-2013 11:51:49 Cambios est�ticos sobre todo
% 20-Dec-2013 20:07:17 Hago que por defecto identifique un m�ximo de 1000
% frames por trozo. Adem�s, hago que por defecto rechaze frames en los que
% haya m�s de 2*n_peces manchas.
% 12-Dec-2013 20:01:21 Hago que muestre el n�mero de procesadores que est�
% usando realmente
% 12-Dec-2013 15:37:42 Hago que guarde la versi�n de Matlab para elegir
% entre mmreader y VideoReader
% 02-Dec-2013 19:10:14 Hago que encripte autom�ticamente si no hay
% argumentos de entrada
% 01-Dec-2013 10:25:32 Hago que las trayectorias se guarden tambi�n en el
% directorio del v�deo, y mejoro el mensaje de despedida.
% 29-Nov-2013 10:04:14 Hago que por defecto sean 5 individuos.
% 29-Nov-2013 08:20:35 Hago que los errores salgan m�s elegantemente.
% Adem�s cambio la fecha en datosegm.version (esto deber�a hacerlo m�s a
% menudo)
% 28-Nov-2013 20:33:03 Hago que muestre un aviso al terminar
% 25-Nov-2013 17:46:17 Paso de trayectorias a trajectories, en incoroporo
% la transformaci�n a txt.
% 22-Nov-2013 18:18:55
% 05-Nov-2013 14:54:18 Revierto al sistema de referencias antiguo
% 25-Sep-2013 21:26:12 Meto el nuevo sistema de referencias
% 19-Jul-2013 16:08:12 Hago que se pueda forzar que coja referencias
% incluso cuando hay un solo pez
% 04-Jul-2013 20:50:02 Hago que se pueda controlar el n�mero de procesadores
% 18-Jun-2013 19:22:53 A�ado mancha2centro
% 10-Jun-2013 14:41:27
% 08-Jun-2013 22:07:07 A�ado el fichero de log
% 31-May-2013 20:49:35 A�ado la posibilidad de encriptar los archivos
% 31-May-2013 16:22:41 Hago que pueda funcionar sin entradas (para la versi�n compilada)
% 28-May-2013 19:35:30 A�ado la opci�n de s�lo calcular el datosegm b�sico. Adem�s cambio la salida, de trayectorias a datosegm
% 09-May-2013 09:56:22 Hago que cuando est� activo empezarsinmas reutilice todo por defecto (a menos que empezarsinmas sea 2)
% 30-Apr-2013 11:37:33 Hago que s�lo reutilice las detecci�n de manchas individuales si encuentra indiv
% Adem�s, hago que no abra parpool si est� en Trueno
% 26-Apr-2013 17:40:24 Hago que calcule trozos antes de buscar manchas
% individuales (lo hago para poder borrar segm.pixels cuando subo datos a
% Trueno)
% 25-Apr-2013 17:33:02 A�ado Trueno
% 17-Apr-2013 11:58:53 A�ado la opci�n "save & exit"
% 10-Apr-2013 16:40:53 A�ado la posibilidad de traquear s�lo un intervalo
% 28-Feb-2013 12:01:04 A�ado el detector de blanco y negro o color
% 18-Feb-2013 10:51:00 Hago que pueda reutilizar trozos
% 11-Feb-2013 18:09:13 Hago que pixelsmierda pueda calcularse con el m�todo de contar frames
% 24-Jan-2013 12:55:36 Actualizo a la nueva versi�n de solapamiento2trozos
% 27-Nov-2012 18:43:13 Cambio el formato de reutiliza
% 23-Nov-2012 14:36:50 Integro el panel. Adem�s, hago que pueda meterse datosegm como �nico argumento de entrada.
% Adem�s, reorganizo para que cosas como npixels se calculen al hacer la segmentaci�n, y reordeno para que est�n
% m�s claros los pasos que pueden reutilizarse.
% 23-Nov-2012 08:47:37 Quito cambiacontraste de la llamada a datosegm2segms
% 21-Nov-2012 15:15:36 Meto el panel. Adem�s reorganizo varias cosas. Quito el c�digo que correspond�a al caso de meter
% el v�deo ya segmentado (de todos modos ese c�digo era antiguo y probablemente ya no funcionar�a)
% 16-Nov-2012 14:28:25 Hago que se pueda elegir entre adquisici�n de referencias antigua o nueva.
% 18-Oct-2012 14:44:33 A�ado el panel de control
% 03-Oct-2012 18:49:01 Cambio de excluyezona_intensmed a mascara_intensmed
% 19-Sep-2012 11:31:58 Hago que funcione mejor para v�deos de un solo pez.
% Salta la clasificaci�n de manchas individuales, y no hace la
% identificaci�n usando los mapas.
% 31-Jul-2012 09:28:59 A�ado la posibilidad de excluir una zona para el
% c�lculo de intensmed
% 26-Jul-2012 20:50:24 Hago que vuelva a usar las versiones antiguas de
% datosegm2referencias y datosegm2intervalosbuenos
% 26-Jul-2012 20:39:39 A�ado el tama�o m�ximo de manchas
% 23-Jul-2012 20:42:13 A�ado el limpiador de mierda
% 27-Jun-2012 19:35:36 Hago que pueda funcionar con referencias externas
% 01-Jun-2012 18:51:26 Intento mejorar la eficiencia en el uso de memoria
% 19-May-2012 11:20:21 Hago que pueda reutilizar las refs. individuales.
% 08-May-2012 19:26:00 Modificaciones menores.
% 13-Mar-2012 20:10:03 Incluyo la comprobaci�n de n�mero de peces en cada
% mancha
% 08-Mar-2012 21:07:11 Actualizo, metiendo todos los cambios que he metido
% en identitracking_refsexternas.
% 22-Feb-2012 19:22:25 A�ado la posibilidad de invertir el contraste
% 26-Jan-2012 17:50:16 De momento, nada.
% 06-Dec-2011 11:23:43 A�ado matriznoson
% 18-Nov-2011 14:36:36 A�ado roi
% 10-Nov-2011 17:59:50 Cambio a la nueva segmentaci�n en la que s�lo tiene
% en cuenta la diferencia con el videomedio
% 14-Oct-2011 17:55:14 Cambio el c�lculo de las probabilidades de error.
% APE 11 oct 11 Viene de identitracking_masdedos

% (C) 2014 Alfonso P�rez Escudero, Gonzalo G. de Polavieja, Consejo Superior de Investigaciones Cient�ficas

% umbral=-1 significa que se ajusta manualmente.
% raizarchivos=[] significa que ya est� hecha la segmentaci�n, y buscar�
% datosegm.
%
% mascara_intensmed debe ser una matriz l�gica del mismo tama�o que los
% frames, con unos en la regi�n de la que se quiere sacar la intensidad
% media. Si se deja vac�a, al final se coger� la roi.


function datosegm=idTracker_ATP(directorio,nombrearchivo,directorio_destino,...
    n_peces,umbral_LIGHT,umbral_DARK,reutiliza,roi,cambiacontraste,referencias,...
    mascara_intensmed,solodatosegm,nRefFrames,minBlobPixelsLight,minBlobPixelsDark,...
    maxBlobPixelsLight,maxBlobPixelsDark,DL_thresh,reduceResolution,BkgFrames)

encriptar=false;
existedatosegm=false;

if ispc
    barra='\';
else
    barra='/';
end

try
    
    if nargin==0
        directorio=ultimodir;
        [nombrearchivo,directorio]=uigetfile('*.*','Select video file',directorio);
        ultimodir(directorio);
        if isequal(nombrearchivo,0)
            error('idTracker:WindowClosed','No file selected')
        end
    end
    
    if nargin~=1 % Si s�lo hay un argumento de entrada, ser� datosegm
        if nombrearchivo(end-3)=='.'
            extension=nombrearchivo(end-2:end);
            nombrearchivo=nombrearchivo(1:end-4);
        else
            extension='';
        end
        if nombrearchivo(end)=='1'
            raizarchivos=nombrearchivo(1:end-1);
        else
            raizarchivos=nombrearchivo;
        end
        
        if directorio(end)~=barra
            directorio(end+1)=barra;
        end
        
        if nargin<3 || isempty(directorio_destino)
            directorio_destino=[directorio 'segm\'];
        end
        if isempty(dir(directorio_destino))
            mkdir(directorio_destino)
        end
        if directorio_destino(end)~=barra
            directorio_destino(end+1)=barra;
        end
        
        datosegm=directorio2datosegm(directorio,raizarchivos,directorio_destino,extension);
        existedatosegm=true;
        datosegm.raizarchivo='segm';
        
        if nargin<4 || isempty(n_peces)
            n_peces=5;
        end
        if nargin<5 || isempty(umbral_LIGHT)
            umbral=.85;
        end
        if nargin<6 || isempty(reutiliza)
            reutiliza=false;
        end
        if nargin<7
            roi=[];
        end
        if nargin<8 || isempty(cambiacontraste)
            cambiacontraste=false;
        end
        if nargin<9
            referencias=[];
        end
        if nargin<10
            mascara_intensmed=[];
        end
        
        if length(reutiliza)==1
            datosegm.reutiliza.datosegm=reutiliza;
            datosegm.reutiliza.Background=reutiliza;
            datosegm.reutiliza.Segmentation=reutiliza;
            datosegm.reutiliza.Trozos=reutiliza;
            datosegm.reutiliza.Individualization=reutiliza;
            datosegm.reutiliza.Resegmentation=reutiliza;
            datosegm.reutiliza.References=reutiliza;
            datosegm.reutiliza.Identification=reutiliza;
            datosegm.reutiliza.Trajectories=reutiliza;
            datosegm.reutiliza.FillGaps=reutiliza;
        end
        
        datosegm.n_peces=n_peces;
        datosegm.umbral=umbral_LIGHT;
        datosegm.umbral_LIGHT=umbral_LIGHT;
        datosegm.umbral_DARK=umbral_DARK;
        datosegm.roi=roi;
        datosegm.cambiacontraste=cambiacontraste;
        
        % primerframe_intervalosbuenos=5000; % No se consideran los primeros 5000 frames para las referencias, porque la pared podr�a afectar. Esto es para v�deos de agresi�n.
        datosegm.primerframe_intervalosbuenos=1;
        datosegm.interval=[1 size(datosegm.frame2archivo,1)];
        % primerframe_intervalosbuenos=24*500; % Para el v�deo con mano y techo de Juli�n
        % primerframe_intervalosbuenos=20*500; % Para mi v�deo con mano y techo
        % datosegm.nframes_refs=3000; % ******************************************************** ATP CHANGED THIS (see line below)
        datosegm.nframes_refs=nRefFrames;
        datosegm.ratio_bwdist=2;
        % disp('Guarning reduceresol!')
        % datosegm.reduceresol=1; % Para moscas % ******************************************************** ATP CHANGED THIS (see line below)
        datosegm.reduceresol = reduceResolution;
        % reduceresol=3; % Para ratones
        % reduceresol=2; % Para ratones, c�mara m�s lejos
        
        datosegm.n_procesadores=Inf;
        
        % datosegm.umbral_npixels=250; % ******************************************************** ATP CHANGED THIS (see line below)
        datosegm.umbral_npixels = minBlobPixelsLight;
        datosegm.umbral_npixels_LIGHT=minBlobPixelsLight;
        datosegm.umbral_npixels_DARK=minBlobPixelsDark;
        datosegm.umbral_npixelsmax = maxBlobPixelsDark;
        datosegm.umbral_npixelsmax_LIGHT = maxBlobPixelsLight;
        datosegm.umbral_npixelsmax_DARK = maxBlobPixelsDark;
        datosegm.int_thresh = DL_thresh;
        datosegm.nBkgFrames = BkgFrames;
        datosegm.limpiamierda=true;
        datosegm.refsantiguas=false;
        
        % Cosas que todav�a falta integrar en el panel:
        datosegm.mascara_intensmed=mascara_intensmed;
    else
        datosegm=directorio;
        datosegm.reutiliza.datosegm=true;
        referencias=[];
    end
    
    datosegm.version='20140805T102158';
    datosegm.version_numero='2.1';
    versioninfo.version=datosegm.version;
    versioninfo.version_numero=datosegm.version_numero;    
    save versioninfo versioninfo
    datosegm.MatlabVersion=version;
    [a,datosegm.ordenata]=system('hostname');   
%     datosegm.ordenata=getComputerName;
    datosegm.guarda_menorescell=false;
    datosegm.max_framesportrozo=1000;
    datosegm.max_manchas.relativo=2;
    datosegm.max_manchas.absoluto=10;
    
    if encriptar || ~isfield(datosegm,'encriptar')
        datosegm.encriptar=encriptar;
    end
    
    if nargin<11 || ~solodatosegm
        
        % Cosas que creo que est�n obsoletas
        %         datosegm.umbral_dif=-.4; % Cambiado el 23 jul 12 para el tracking del v�deo guarro de 12 peces de Juli�n.
        % datosegm.usaresta=false; % Para ratones
        
        datosegm=datosegm2progreso(datosegm);
        
        h_panel=[];
        datosegm.estilopixelsmierda=1;
        
        
        if ~isfield(datosegm,'trueno') || isempty(datosegm.trueno)
            datosegm.trueno=false;
        end
%         disp('Guarning!')
        poolobj = gcp('nocreate');
        
        if isempty(poolobj)
            parpool('local');
            clear poolobj
            poolobj = gcp('nocreate');
        end        
                
        if ~isfield(datosegm,'muestrapanel') || datosegm.muestrapanel
            [datosegm,h_panel]=panel_identitracking_ATP(datosegm);
        end
        if ~isfield(datosegm,'ratio_bwdist')
            error('It seems that the video was tracked with an older version of idTracker, not compatible with the new version. Please, erase or rename the folder called "segm", and try to re-track the video from the beginning')
        end
        datosegm.tiempos.total(1)=now;
        datosegm.tiempos.tiempoguardando=0;
        if ~isfield(datosegm,'blancoynegro') || isempty(datosegm.blancoynegro)
            datosegm.tiempos.colorbn(1)=now;
            datosegm=datosegm2colorbn(datosegm);
            datosegm.tiempos.colorbn(2)=now;
        end
        if ~isfield(datosegm,'umbral_npixelsmax');
            datosegm.umbral_npixelsmax=10000*datosegm.reduceresol^2;
        end
%         datosegm.id_log=fopen([datosegm.directorio 'Log' datestr(now,30) '.txt'],'w');
        tic
        variable=datosegm;
        save([datosegm.directorio 'datosegm.mat'],'variable')
        datosegm.tiempos.tiempoguardando=datosegm.tiempos.tiempoguardando+toc;
        
%         fprintf(datosegm.id_log,'Version %s\n\n',datosegm.version)
%         fprintf(datosegm.id_log,'%s - Start\n',datestr(now,30));
        
        if ~isfield(datosegm,'saltatodo') || ~datosegm.saltatodo
            %     if isfield(datosegm,'trueno') && datosegm.trueno==1
            %         unidad=pwd;
            %         unidad=unidad(1);
            %         directorio_subemierdas=[unidad ':\ColaIdentitrackingClusterTrueno'];
            %         if isempty(dir(directorio_subemierdas))
            %             mkdir(directorio_subemierdas)
            %         end
            %         copyfile([datosegm.directorio 'datosegm.mat'],[directorio_subemierdas '\datosegm' datestr(now,30) '.mat'])
            %     end
            %     if ispc && parpool('size')==0 % S�lo abre parpool si est� en un pc, porque si es linux indica que estamos en Trueno.
            % %                 disp('Guarning! No se abre parpool!')
            %         parpool open % Lo abre con la configuraci�n por defecto par el ordenador
            %         %         disp('Guarning! parpool se abre s�lo con 6 procesadores')
            % %                 parpool open local 4
            %         %     parpool open local 8
            %     end
            
            nprocesadores_abiertos=poolobj.NumWorkers;
%             if nprocesadores_abiertos~=datosegm.n_procesadores && (datosegm.n_procesadores~=Inf || nprocesadores_abiertos~=feature('numCores'))
%                 if nprocesadores_abiertos~=0
%                     delete(gcp)
%                 end
%                 if datosegm.n_procesadores==Inf
%                     parpool open % Configuraci�n por defecto del ordenador
%                 else
%                     parpool('open','local',datosegm.n_procesadores)
%                 end
%             end
            datosegm.n_procesadores_real=poolobj.NumWorkers;
            if isfield(h_panel,'n_procesadores') && ishandle(h_panel.n_procesadores)
                set(h_panel.n_procesadores,'String',num2str(datosegm.n_procesadores_real))
                drawnow
            end
            
%             fprintf(datosegm.id_log,'%s - parpool is open\n',datestr(now,30));
            if isfield(datosegm,'empezarsinmas') && datosegm.empezarsinmas==1
                datosegm.reutiliza.Background=1;
                datosegm.reutiliza.Segmentation=1;
                datosegm.reutiliza.Individualization=1;
                datosegm.reutiliza.Trozos=1;
                datosegm.reutiliza.Resegmentation=1;
                datosegm.reutiliza.References=1;
                datosegm.reutiliza.Identification=1;
                datosegm.reutiliza.Trajectories=1;
                datosegm.reutiliza.FillGaps=1;
                variable=datosegm;
                save([datosegm.directorio 'datosegm.mat'],'variable')
            end
            
%             % Por si datosegm viene de una versi�n antigua
%             pasos={'Background' 'Segmentation' 'Individualization' 'Trozos' 'Resegmentation' 'References' 'Identification' 'Trajectories' 'FillGaps'};
%             for c_pasos=1:length(pasos)
%                 if ~isfield(datosegm.reutiliza,pasos{c_pasos})
%                     datosegm.reutiliza.(pasos{c_pasos})=1;
%                 end
%                 if ~isfield(datosegm.progreso,pasos{c_pasos})
%                     datosegm.reutiliza.(pasos{c_pasos})=0;
%                 end
%             end
            
            % V�deo medio
            n_frames=size(datosegm.frame2archivo,1);
            if datosegm.limpiamierda && (~isfield(datosegm,'videomedio') || isempty(datosegm.videomedio) || ~datosegm.reutiliza.Background)
                datosegm.tiempos.videomedio(1)=now;
                if isempty(h_panel) || ~ishandle(h_panel.fig)
                    fprintf('V�deo medio\n')
                end
                datosegm=datosegm2videomedio(datosegm,100,h_panel);
                if datosegm.limpiamierda
                    datosegm=datosegm2datosegm_pixelsmierda(datosegm);
                end
                datosegm.tiempos.videomedio(2)=now;
                tic
                variable=datosegm;
                save([datosegm.directorio 'datosegm.mat'],'variable')
                datosegm.tiempos.tiempoguardando=datosegm.tiempos.tiempoguardando+toc;
            end
%             fprintf(datosegm.id_log,'%s - Videomedio done\n',datestr(now,30));
            % % Selecci�n manual del umbral
            % if umbral==-1
            %     umbral=datosegm2umbral_manual(datosegm);
            %     datosegm.umbral=umbral;
            %     save([datosegm.directorio 'datosegm.mat'],'datosegm')
            % end
            
            % Segmentaci�n, solapamiento y mapas
            if isempty(h_panel) || ~ishandle(h_panel.fig)
                fprintf('Segmentaci�n, solapamiento y mapas\n')
            end
            if ~datosegm.reutiliza.Segmentation || isempty(dir([datosegm.directorio 'segm_' num2str(size(datosegm.archivo2frame,1)) '.mat']))
%                 fprintf(datosegm.id_log,'%s - Segmentation...\n',datestr(now,30));
                datosegm.tiempos.segm(1)=now;
                [datosegm,solapamiento,npixels,segmbuena,borde,mancha2centro,max_bwdist,bwdist_centro,max_distacentro]=datosegm2segm(datosegm,h_panel);
                datosegm.tiempos.segm(2)=now;
                tic
                disp('h1')
                variable=datosegm;
                save([datosegm.directorio 'datosegm.mat'],'variable')
                variable=solapamiento;
                save([datosegm.directorio 'solapamiento.mat'],'variable')
                npixelsyotros.npixels=npixels;
                npixelsyotros.segmbuena=segmbuena;
                npixelsyotros.borde=borde;
                npixelsyotros.mancha2centro=mancha2centro;
                npixelsyotros.max_bwdist=max_bwdist;
                npixelsyotros.bwdist_centro=bwdist_centro;
                npixelsyotros.max_distacentro=max_distacentro;
                variable=npixelsyotros;
                save([datosegm.directorio 'npixelsyotros'],'variable')
                datosegm.tiempos.tiempoguardando=datosegm.tiempos.tiempoguardando+toc;
            else
                disp('h3')
                if isempty(h_panel)
                    fprintf('Reutiliza segmentaci�n anterior.\n')
                end
                load([datosegm.directorio 'datosegm.mat'])
                datosegm=variable;
                load([datosegm.directorio 'solapamiento.mat'])
                solapamiento=variable;
                load([datosegm.directorio 'npixelsyotros'])
                npixelsyotros=variable;
                npixels=npixelsyotros.npixels;
                segmbuena=npixelsyotros.segmbuena;
                borde=npixelsyotros.borde;
                mancha2centro=npixelsyotros.mancha2centro;
            end
            clear npixelsyotros
%             fprintf(datosegm.id_log,'%s - Segmentation done\n',datestr(now,30));
            if (~isfield(datosegm,'trueno') || ~(datosegm.trueno==1))
                if (~datosegm.reutiliza.Trozos || isempty(dir([datosegm.directorio 'trozos.mat'])))
                    if datosegm.n_peces>1 || (isfield(datosegm,'siemprerefs') && datosegm.siemprerefs)
                        datosegm.tiempos.trozos(1)=now;
                        [trozos,solapos]=solapamiento2trozos(solapamiento,npixels,datosegm,mancha2centro);                        
                    else
                        trozos=solapamiento2trozos(solapamiento,npixels,datosegm,mancha2centro);
                        solapos=[];
                    end
                    [conectan,conviven,solapan]=trozos2conectatrozos(trozos,solapamiento);
                    datosegm.tiempos.trozos(2)=now;
                    tic
                    save([datosegm.directorio 'conectanconviven.mat'],'conectan','conviven','solapan')
                    variable=datosegm;
                    save([datosegm.directorio 'datosegm.mat'],'variable')
                    trozosolapos.trozos=trozos;
                    trozosolapos.solapos=solapos;
                    variable=trozosolapos;
                    save([datosegm.directorio 'trozos'],'variable')
                    datosegm.tiempos.tiempoguardando=datosegm.tiempos.tiempoguardando+toc;
                else
                    load([datosegm.directorio 'trozos.mat']);
                    trozosolapos=variable;
                    trozos=trozosolapos.trozos;
                    solapos=trozosolapos.solapos;
                    if ~isempty(dir([datosegm.directorio 'conectanconviven.mat']))
                        load([datosegm.directorio 'conectanconviven.mat'])
                    end
                end
                clear trozosolapos
            end
            
%             fprintf(datosegm.id_log,'%s - Trozos done\n',datestr(now,30));
            
            
            
            if ~isfield(datosegm,'trueno') || ~(datosegm.trueno==1) % trueno ser� 2 cuando est� preparado para continuar en el cluster.
                % Si es necesario, recalcula tam_mapas (esto no deber�a hacer falta
                % pr�cticamente nunca, s�lo cuando se ha cancelado el tracking en
                % un momento muy concreto de la segmentaci�n de algunos v�deos)
                if ~isfield(datosegm,'tam_mapas')
                    load([datosegm.directorio 'npixelsyotros.mat']);
                    npixelsyotros=variable;
                    ind_mapa=find(npixelsyotros.segmbuena & ~npixelsyotros.borde,1);
                    [ind_frame,ind_mancha]=ind2sub(size(npixelsyotros.segmbuena),ind_mapa);
                    load([datosegm.directorio 'segm_' num2str(datosegm.frame2archivo(ind_frame,1)) '.mat']);
                    segm=variable;
                    datosegm.tam_mapas=size(segm(datosegm.frame2archivo(ind_frame,2)).mapas{ind_mancha});
                    variable=datosegm;
                    save([datosegm.directorio 'datosegm.mat'],'variable')
                end
                
                % Crea las referencias para distinguir peces individuales de manchas
                % m�ltiples.
                nframes_refindiv=5;
                refs_indiv=NaN([datosegm.tam_mapas nframes_refindiv*datosegm.n_peces]);
                if datosegm.n_peces>1 || (isfield(datosegm,'siemprerefs') && datosegm.siemprerefs)
                    if isempty(referencias)
                        % Asume que los frames con igual n�mero de manchas que de peces
                        % tienen todas las manchas individuales
                        if ~datosegm.reutiliza.Individualization || isempty(dir([datosegm.directorio 'refs_indiv.mat']))
                            framesbuenos=find(datosegm.n_manchas==datosegm.n_peces);
                            % disp('Guarning, guarning, guarning, guarning!!!!! Esto s�lo funciona para los v�deos de learning de Sara!')
                            % if isempty(framesbuenos) && n_peces==2
                            %     datosegm.n_peces=1;
                            %     n_peces=1;
                            %     framesbuenos=find(datosegm.n_manchas==n_peces);
                            % end
                            c_refs=0;
                            for c_frames=1:nframes_refindiv
                                frame=framesbuenos(ceil(rand*length(framesbuenos))); % Coge uno al azar
                                archivo_act=datosegm.frame2archivo(frame,1);
                                frame_arch=datosegm.frame2archivo(frame,2);
                                load([datosegm.directorio datosegm.raizarchivo '_' num2str(archivo_act)])
                                segm=variable;
                                for c_peces=1:datosegm.n_peces
                                    if ~isempty(segm(frame_arch).mapas{c_peces})
                                        c_refs=c_refs+1;
                                        refs_indiv(:,:,:,c_refs)=segm(frame_arch).mapas{c_peces};
                                    end
                                end % c_peces
                            end % c_frames
                            refs_indiv=refs_indiv(:,:,:,1:c_refs);
                            mat_validos=datosegm.indvalidos;
                            indvalidos{1}=find(mat_validos(:,:,1));
                            mat_validos(:,:,1)=false;
                            indvalidos{2}=find(mat_validos);
                            [menores,errores_pezindiv]=comparamapas(refs_indiv,{refs_indiv},indvalidos);
                            errores_pezindiv=errores_pezindiv{1};
                            errores_pezindiv(errores_pezindiv==0)=NaN;
                            errores_pezindiv=squeeze(min(errores_pezindiv,[],1));
                            errores_pezindiv=sum(errores_pezindiv,2);
                            refsindiv.refs_indiv=refs_indiv;
                            refsindiv.errores_pezindiv=errores_pezindiv;
                            variable=refsindiv;
                            save([datosegm.directorio 'refs_indiv.mat'],'variable')
                            datosegm.umbral_errorindiv=mean(errores_pezindiv)+std(errores_pezindiv)*3; % Damos 3 sd's.
                            variable=datosegm;
                            save([datosegm.directorio 'datosegm.mat'],'variable')
                            clear segm
                        else
                            fprintf('Reutiliza referencias individuales.\n')
                            load([datosegm.directorio 'refs_indiv.mat'])
                            refsindiv=variable;
                            refs_indiv=refsindiv.refs_indiv;
                            errores_pezindiv=refsindiv.errores_pezindiv;
                        end
                        clear refsindiv
                    else
                        disp('Guarning! Usa 10 sd para las manchas individuales, y cuando son refs. internas usa 3 sd')
                        nframes_refindiv=ceil(100/datosegm.n_peces);
                        n_refs=length(referencias);
                        refs_indiv=NaN([datosegm.tam_mapas nframes_refindiv*n_refs]);
                        cframes_tot=0;
                        for c_refs=1:n_refs
                            n_frames=size(referencias{c_refs},4);
                            indices=equiespaciados(nframes_refindiv,n_frames);
                            for c_frames=indices
                                cframes_tot=cframes_tot+1;
                                refs_indiv(:,:,:,cframes_tot)=referencias{c_refs}(:,:,:,c_frames);
                            end % c_frames
                        end % c_refs
                        % c_refs=0;
                        %
                        % for c_frames=1:nframes_refindiv
                        %     frame=framesbuenos(ceil(rand*length(framesbuenos))); % Coge uno al azar
                        %     archivo_act=datosegm.frame2archivo(frame,1);
                        %     frame_arch=datosegm.frame2archivo(frame,2);
                        %     load([datosegm.directorio datosegm.raizarchivo '_' num2str(archivo_act)],'segm')
                        %     for c_peces=1:n_peces
                        %         if ~isempty(segm(frame_arch).mapas{c_peces})
                        %             c_refs=c_refs+1;
                        %             refs_indiv(:,:,:,c_refs)=segm(frame_arch).mapas{c_peces};
                        %         end
                        %     end % c_peces
                        % end % c_frames
                        % refs_indiv=refs_indiv(:,:,:,1:c_refs);
                        mat_validos=datosegm.indvalidos;
                        indvalidos{1}=find(mat_validos(:,:,1));
                        mat_validos(:,:,1)=false;
                        indvalidos{2}=find(mat_validos);
                        [menores,errores_pezindiv]=comparamapas(refs_indiv,{refs_indiv},indvalidos);
                        errores_pezindiv=errores_pezindiv{1};
                        errores_pezindiv(errores_pezindiv==0)=NaN;
                        errores_pezindiv=squeeze(min(errores_pezindiv,[],1));
                        errores_pezindiv=sum(errores_pezindiv,2);
                        refsindiv.refs_indiv=refs_indiv;
                        refsindiv.errores_pezindiv=errores_pezindiv;
                        variable=refsindiv;
                        save([datosegm.directorio 'refs_indiv.mat'],'variable')
                        %                 save([datosegm.directorio 'refs_indiv.mat'],'refs_indiv','errores_pezindiv')
                        datosegm.umbral_errorindiv=mean(errores_pezindiv)+std(errores_pezindiv)*10; % Damos 10 sd's. Cambio de 3 a 10 el 24 de feb del 12. Quiz� 3 estaba bien para los v�deos de cerca, pero no para los de lejos.
                        variable=datosegm;
                        save([datosegm.directorio 'datosegm.mat'],'variable')
                        clear refsindiv
                    end
                else
                    refs_indiv=[];
                end
                
%                 fprintf(datosegm.id_log,'%s - Refs_indiv done\n',datestr(now,30));
                
                % Distingue manchas individuales de manchas colectivas. Siempre asume que
                % cuando hay tantas manchas como peces, son todas individuales
                n_archivos=size(datosegm.archivo2frame,1);
                load([datosegm.directorio datosegm.raizarchivo '_' num2str(n_archivos)]);
                segm=variable;
                if datosegm.n_peces>1 || (isfield(datosegm,'siemprerefs') && datosegm.siemprerefs)
                    if ~datosegm.reutiliza.Individualization || ~isfield(segm,'indiv') || isempty(dir([datosegm.directorio 'indiv.mat']))
                        datosegm.tiempos.indiv(1)=now;
                        indiv=datosegm2segm_indiv(datosegm,refs_indiv,h_panel); % Indiv se va guardando desde dentro                        
                        datosegm.tiempos.indiv(2)=now;
                        variable=datosegm;
                        save([datosegm.directorio 'datosegm.mat'],'variable')
                    else
                        disp('Reutiliza identificaci�n de manchas individuales')
                        load([datosegm.directorio 'indiv.mat'])
                        indiv=variable;
                    end
                end
                
%                 fprintf(datosegm.id_log,'%s - Individualization done\n',datestr(now,30));
%                 
                % Resegmentaci�n               
                if datosegm.n_peces>1 && (~isfield(datosegm,'resegmentar') || datosegm.resegmentar)
                    n_archivos=size(datosegm.archivo2frame,1);
                    load([datosegm.directorio datosegm.raizarchivo '_' num2str(n_archivos)])
                    segm=variable;
                    if ~datosegm.reutiliza.Resegmentation || ~isfield(segm,'resegmentado')
                        % Resegmentaci�n (todo se guarda dentro)
                        datosegm.tiempos.resegmentacion(1)=now;
                        [datosegm,npixelsyotros,solapamiento,trozos,solapos,conectan,conviven,solapan,indiv]=datosegm2resegmentacion(datosegm,h_panel);
                        campos=fieldnames(npixelsyotros);
                        for c_campos=1:length(campos)
                            eval([campos{c_campos} '=npixelsyotros.' campos{c_campos} ';'])
                        end
                        % Se asegura de que el �ltimo segm tenga el campo
                        % resegmentado, para saber luego si este paso ha
                        % terminado
                        n_archivos=size(datosegm.archivo2frame,1);
                        load([datosegm.directorio datosegm.raizarchivo '_' num2str(n_archivos)])
                        segm=variable;
                        if ~isfield(segm,'resegmentado')
                            segm(1).resegmentado=[];
                        end
                        variable=segm;
                        save([datosegm.directorio datosegm.raizarchivo '_' num2str(n_archivos)],'variable')
                        datosegm.tiempos.resegmentacion(2)=now;
                        variable=datosegm;
                        save([datosegm.directorio 'datosegm.mat'],'variable')
                    else
                        disp('Reutiliza resegmentaci�n')
                    end
                end
                
                
                if isempty(referencias)
                    if datosegm.n_peces>1 || (isfield(datosegm,'siemprerefs') && datosegm.siemprerefs)
                        % Intervalos v�lidos para las referencias
                        fprintf('Intervalos buenos\n')
                        if ~datosegm.reutiliza.References || isempty(dir([datosegm.directorio 'intervalosbuenos.mat']))
                            %         intervalosbuenos=datosegm2intervalosbuenos(datosegm,primerframe_i
                            %         ntervalosbuenos);
                            if datosegm.refsantiguas
                                intervalosbuenos=datosegm2intervalosbuenos_20120724T110320(datosegm,datosegm.primerframe_intervalosbuenos);
                            else
                                datosegm.tiempos.intervalosbuenos(1)=now;                                
                                intervalosbuenos=datosegm2intervalosbuenos(datosegm,trozos,solapos,indiv,segmbuena,borde,datosegm.primerframe_intervalosbuenos,1);
                                datosegm.tiempos.intervalosbuenos(2)=now;
                                variable=datosegm;
                                save([datosegm.directorio 'datosegm.mat'],'variable')
                            end
                            tic
                            variable=intervalosbuenos;
                            save([datosegm.directorio 'intervalosbuenos.mat'],'variable')
                            datosegm.tiempos.tiempoguardando=datosegm.tiempos.tiempoguardando+toc;
                        else
                            fprintf('Reutiliza intervalosbuenos.\n')
                            load([datosegm.directorio 'intervalosbuenos.mat'])
                            intervalosbuenos=variable;
                        end
                        
                        % Referencias
                        fprintf('Referencias\n')
                        if ~datosegm.reutiliza.References || isempty(dir([datosegm.directorio 'referencias.mat']))
                            %         [referencias,framesescogidos]=datosegm2referencias(datosegm,intervalosbuenos,nframes_ref);
                            if datosegm.refsantiguas
                                [referencias,framesescogidos]=datosegm2referencias_20120724T172537(datosegm,intervalosbuenos,datosegm.nframes_refs);
                                listamapas=[];
                            else
                                datosegm.tiempos.referencias(1)=now;
                                %                         [referencias,listamapas]=datosegm2referencias_renuevo(datosegm,intervalosbuenos,trozos,conviven,datosegm.nframes_refs,h_panel);
                                [referencias,listamapas]=datosegm2referencias(datosegm,intervalosbuenos,trozos,datosegm.nframes_refs,h_panel);
                                datosegm.tiempos.referencias(2)=now;
                                variable=datosegm;
                                save([datosegm.directorio 'datosegm.mat'],'variable')
                            end
                            tic
                            refs.referencias=referencias;
                            refs.listamapas=listamapas;
                            variable=refs;
                            save([datosegm.directorio 'referencias.mat'],'variable')
                            datosegm.tiempos.tiempoguardando=datosegm.tiempos.tiempoguardando+toc;
                        else
                            fprintf('Reutiliza referencias.\n')
                            load([datosegm.directorio 'referencias.mat'])
                            refs=variable;
                            referencias=refs.referencias;
                            listamapas=refs.listamapas;
                        end
                        clear refs
                    else
                        refs.referencias=[];
                        refs.listamapas=[];
                        variable=refs;
                        save([datosegm.directorio 'referencias'],'variable');
                        clear refs
                    end
                else
                    refs.referencias=referencias;
                    refs.listamapas=[];
                    variable=refs;
                    save([datosegm.directorio 'referencias'],'variable')
                    clear refs
                end
                
%                 fprintf(datosegm.id_log,'%s - References done\n',datestr(now,30));
                
                % % An�lisis de las referencias para buscar probabilidades de error
                % fprintf('An�lisis de probabilidades de error\n')
                % if ~reutiliza(7) || isempty(whos('maxprob'))
                %     mat_validos=datosegm.indvalidos;
                %     indvalidos{1}=find(mat_validos(:,:,1));
                %     mat_validos(:,:,1)=false;
                %     indvalidos{2}=find(mat_validos);
                %     [probasignacion,asignacion,n_validos,menores]=referencias2probasignacion(referencias,indvalidos,200);
                %     save([datosegm.directorio 'probasignacion.mat'],'probasignacion','asignacion','n_validos','menores')
                % else
                %     fprintf('Reutiliza an�lisis de errores.\n')
                % end
                
                % Identificaci�n
                % prob_id=datosegm2segm_id(datosegm,referencias,maxprob,binsc,indvalidos,reutiliza);
                % save([datosegm.directorio 'prob_id.mat'],'prob_id')
                
                % Trozos, probabilidades y trayectorias
                % solapamiento=datosegm2solapamiento(datosegm);
                % disp('Guarning, guarning!!')
                
                if ~isfield(datosegm,'solohastareferencias') || ~datosegm.solohastareferencias
                    
                    % Limpio la memoria con la esperanza de que ayude algo
                    a=who;
                    for c=1:length(a)
                        if ~strcmp(a{c},'datosegm') && ~strcmp(a{c},'h_panel') && ~strcmp(a{c},'a')
                            clear(a{c})
                        end
                    end
                    clear a
                    encriptar=datosegm.encriptar;
                    
                    if datosegm.n_peces>1 || (isfield(datosegm,'siemprerefs') && datosegm.siemprerefs)
                        load([datosegm.directorio 'referencias'])
                        refs=variable;
                        referencias=refs.referencias;
                        clear refs
                    end
                    load([datosegm.directorio 'trozos'])
                    trozosolapos=variable;
                    trozos=trozosolapos.trozos;
                    solapos=trozosolapos.solapos;
                    clear trozosolapos
                    
                    if datosegm.progreso.Identification<1 || ~datosegm.reutiliza.Identification
                        if datosegm.n_peces>1 || (isfield(datosegm,'siemprerefs') && datosegm.siemprerefs)
                            mat_validos=datosegm.indvalidos;
                            indvalidos{1}=find(mat_validos(:,:,1));
                            mat_validos(:,:,1)=false;
                            indvalidos{2}=find(mat_validos);
                            datosegm.tiempos.identificacion(1)=now;
                            mancha2id=trozos2mancha2id(datosegm,trozos,solapos,indvalidos,referencias,[],1,h_panel); % mancha2id se guarda autom�ticamente desde dentro del programa
                            datosegm.tiempos.identificacion(2)=now;
                            variable=datosegm;
                            save([datosegm.directorio 'datosegm.mat'],'variable')
                        else
                            load([datosegm.directorio 'npixelsyotros'])
                            npixelsyotros=variable;
                            load([datosegm.directorio 'conectanconviven'])
                            identificaciones.mancha2id=trozos2mancha2id_unbicho(trozos,npixelsyotros.mancha2centro,conviven);
                            mancha2id=identificaciones.mancha2id;
                            variable=identificaciones;
                            save([datosegm.directorio 'mancha2id'],'variable')
                        end
                    else
                        load([datosegm.directorio 'mancha2id.mat']);
                        man2id=variable;
                        mancha2id=man2id.mancha2id;
                        disp('Reutiliza identificaci�n')
                    end
                    
%                     fprintf(datosegm.id_log,'%s - Identification done\n',datestr(now,30));
                    
                    if datosegm.n_peces>1 || (isfield(datosegm,'siemprerefs') && datosegm.siemprerefs)
                        clear referencias
                        set(h_panel.waitTrajectories,'XData',[0 0 .1 .1])
                        set(h_panel.textowaitTrajectories,'String',[num2str(round(sum(.1)*100)) ' %'])
                        load([datosegm.directorio 'solapamiento']);     
                        solapamiento=variable;
                        set(h_panel.waitTrajectories,'XData',[0 0 .25 .25])
                        set(h_panel.textowaitTrajectories,'String',[num2str(round(sum(.25)*100)) ' %'])
                        idtrozos=mancha2id2idtrozos(datosegm,trozos,solapos,mancha2id);
                        probtrozos=idtrozos2probtrozos(idtrozos);
                        % idtrozos=trozos2id_trozos(datosegm,trozos,solapos,indvalidos,referencias,[]);
                        idprobtrozos.idtrozos=idtrozos;
                        idprobtrozos.probtrozos=probtrozos;
                        variable=idprobtrozos;
                        save([datosegm.directorio 'idtrozos.mat'],'variable')
                        clear idprobtrozos
                        
                        set(h_panel.waitTrajectories,'XData',[0 0 .5 .5])
                        set(h_panel.textowaitTrajectories,'String',[num2str(round(sum(.5)*100)) ' %'])
                        
%                         fprintf(datosegm.id_log,'%s - idtrozos done\n',datestr(now,30));
                        
                        % logprob_idtrozos=trozosprobid2idtrozos(trozos,prob_id);
                        % [mancha2pez,trozo2pez,prob_error,logprob_rel]=prob_idtrozos2idtrozos(trozos,logprob_idtrozos);
                        % matriznoson=trozos2matriznoson(trozos);
                        % save([datosegm.directorio 'matriznoson.mat'],'matriznoson')
                        
                        
%                         fprintf(datosegm.id_log,'%s - conectanconviven done\n',datestr(now,30));
% % %                         
                        % C�lculo del n�mero de peces en cada mancha
                        % [mancha2indiv,mancha2borde]=datosegm2mancha2indiv(datosegm);
                        % save([datosegm.directorio 'mancha2indiv.mat'],'mancha2indiv','mancha2borde')
                        load([datosegm.directorio 'conectanconviven.mat'])
                        [mancha2pez,trozo2pez,probtrozos_relac]=probtrozos2identidades(trozos,probtrozos,conviven);
                        % Para tener en cuenta la din�mica, usar esto:
                        %             [probramificaciones,relevantes]=trozos2probramificaciones(trozos,solapan);
                        %             matriznoson=(probramificaciones==0 & relevantes)*.9; % Asumo que hay una probabilidad de 0.1 de que haya un error en matriznoson
                        %             matriznoson(conviven)=1;
                        %             [mancha2pez,trozo2pez,probtrozos_relac]=probtrozos2identidades(trozos,probtrozos,matriznoson);
                        
                        % [mancha2pez,trozo2pez,framesafavor,id_relac]=idtrozos2identidades(trozos,idtrozos,conviven);
                        man2pez.mancha2pez=mancha2pez;
                        man2pez.trozo2pez=trozo2pez;
                        man2pez.probtrozos_relac=probtrozos_relac;
                        
                        
                        %             man2pez.probramificaciones=probramificaciones;
                        %             man2pez.relevantes=relevantes;
                        %             man2pez.matriznoson=matriznoson;
                    else
                        mancha2pez=mancha2id;
                        mancha2pez(mancha2pez==0)=NaN;
%                         trozo2pez=NaN(1,max(trozos(:)));
%                         trozo2pez(trozos(mancha2pez==1))=1;
                        man2pez.mancha2pez=mancha2pez;
%                         man2pez.trozo2pez=trozo2pez;
                        probtrozos_relac=[];
                    end
                    variable=man2pez;
                    save([datosegm.directorio 'mancha2pez.mat'],'variable')
                    clear man2pez
%                     fprintf(datosegm.id_log,'%s - mancha2pez done\n',datestr(now,30));

                    set(h_panel.waitTrajectories,'XData',[0 0 .75 .75])
                    set(h_panel.textowaitTrajectories,'String',[num2str(round(sum(.75)*100)) ' %'])
                        
%                     variable=load([datosegm.directorio 'npixelsyotros.mat']);
%                     mancha2centro=variable.mancha2centro;
%                     [trajectories,probtrajectories]=mancha2pez2trayectorias(datosegm,mancha2pez,trozos,probtrozos_relac,mancha2centro);
%                     save([datosegm.directorio 'trajectories.mat'],'trajectories','probtrajectories')
%                     save([datosegm.directorio_videos 'trajectories.mat'],'trajectories','probtrajectories')
% %                     trajectories2txt(trajectories,[datosegm.directorio 'trajectories.txt'])
%                     trajectories2txt(trajectories,[datosegm.directorio_videos 'trajectories.txt'])
                    
                    set(h_panel.waitTrajectories,'XData',[0 0 1 1])
                    set(h_panel.textowaitTrajectories,'String',[num2str(round(sum(1)*100)) ' %'])
                    
                    load([datosegm.directorio 'mancha2pez.mat'])
                    man2pez=variable;
                    load([datosegm.directorio 'npixelsyotros.mat'])
                    npixelsyotros=variable;
                    [trajectories,probtrajectories]=mancha2pez2trayectorias(datosegm,man2pez.mancha2pez,trozos,[],npixelsyotros.mancha2centro);
                    save([datosegm.directorio 'trajectories.mat'],'trajectories','probtrajectories')
                    save([datosegm.directorio_videos 'trajectories.mat'],'trajectories','probtrajectories')
                    trajectories2txt(trajectories,probtrajectories,[datosegm.directorio_videos 'trajectories.txt'])
                    
                    
                    datosegm.tiempos.fillgaps(1)=now;
                    if datosegm.n_peces>1
                        datosegm2smartinterp2(datosegm,h_panel);
                        load([datosegm.directorio 'mancha2pez_nogaps.mat'])          
                        man2pez=variable;
                        [trajectories,probtrajectories]=mancha2pez2trayectorias(datosegm,man2pez.mancha2pez,trozos,man2pez.probtrozos_relac,man2pez.mancha2centro);
                        probtrajectories(man2pez.tiporefit==1)=-1;
                        probtrajectories(man2pez.tiporefit>=2)=-2;
                        save([datosegm.directorio 'trajectories_nogaps.mat'],'trajectories','probtrajectories')
                        save([datosegm.directorio_videos 'trajectories_nogaps.mat'],'trajectories','probtrajectories')
                        trajectories2txt(trajectories,probtrajectories,[datosegm.directorio_videos 'trajectories_nogaps.txt'])
                    end
                    
                    
                    
                    
                    progreso=1;
                    set(h_panel.waitFillGaps,'XData',[0 0 progreso progreso])
                    set(h_panel.textowaitFillGaps,'String',[num2str(round(progreso*100)) ' %'])

                    datosegm.tiempos.fillgaps(2)=now;
                    
                    
                    datosegm.tiempos.total(2)=now;
                    variable=datosegm;
                    save([datosegm.directorio 'datosegm.mat'],'variable')
%                     fprintf(datosegm.id_log,'%s - Fin.\n',datestr(now,30));
                    %         msgbox(sprintf('Tracking finished! :-)\n\nThe results are in the files named ''trajectories''\nin folder %s',datosegm.directorio),'Job done')
                    poolobj = gcp('nocreate');
                    if ~isempty(poolobj)
                        delete(poolobj);
                    end
                    if ~isfield(datosegm,'muestrapanel') || datosegm.muestrapanel
                        close(h_panel.fig)
                        despedida(datosegm,trajectories,probtrajectories)
                    end
                end % if not soloreferencias
            else
                datosegm.trueno=2; % Para que s� contin�e la pr�xima vez que se ejecute
                variable=datosegm;
                save([datosegm.directorio 'datosegm.mat'],'variable')
            end % if not trueno
        else
            close(h_panel.fig)
        end % if no saltatodo
%         fclose(datosegm.id_log);
    end % if no saltatodo (por solodatosegm)
catch me
    try
        delete(gcp);
    catch
    end
    if nargin==0
        boton='Exit & create error log file';
        if strcmpi(me.identifier,'idTracker:WindowClosed')
            boton = questdlg('It seems you have closed idTracker. The tracking has stoped.','idTracker has closed','Exit','Exit & create error log file','Exit');
        end
        if strcmpi(boton,'Exit & create error log file')
            erroraco.me.identifier=me.identifier; % Voy campo por campo para que quede en un struct en vez de en un error object
            erroraco.me.message=me.message;
            erroraco.me.cause=me.cause;
            erroraco.me.stack=me.stack;
            try
                if exist('datosegm','var')
                    erroraco.datosegm=datosegm;
                    variable=erroraco;
                    save([datosegm.directorio 'idTrackerError' datestr(now,30)],'variable')
                    errorfile=[datosegm.directorio 'idTrackerError' datestr(now,30)];
                else
%                     unidad=pwd;
%                     unidad=unidad(1);
                    variable=erroraco;
                    save(['.\idTrackerError' datestr(now,30)],'variable')
                    errorfile=['.\idTrackerError' datestr(now,30)];
                end
            catch
                error(sprintf('An error has occured, and idTracker must close.\n\nError log file could not be created (maybe your disk is full?)\n\nPlease, try again and if the error persist report it to bugs.idtracker@gmail.com\n\nSorry for the inconvenience!'))
            end
            error(sprintf('An error has occured, and idTracker must close.\nThe error was\n\n%s\n\nAn error log file has been created at\n%s\n\nPlease, try again and if the error persists send the error log file to bugs.idtracker@gmail.com\n\nSorry for the inconvenience!',me.message,errorfile))
        end % if hay error log
    else
        getReport(me)
    end
end
