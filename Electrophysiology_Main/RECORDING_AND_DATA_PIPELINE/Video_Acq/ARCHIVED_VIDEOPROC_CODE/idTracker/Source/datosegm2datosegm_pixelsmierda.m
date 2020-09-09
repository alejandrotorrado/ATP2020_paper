% APE 11 feb 12

% (C) 2014 Alfonso P�rez Escudero, Gonzalo G. de Polavieja, Consejo Superior de Investigaciones Cient�ficas

function datosegm=datosegm2datosegm_pixelsmierda(datosegm)

if isfield(datosegm,'estilopixelsmierda') && datosegm.estilopixelsmierda==2 % Usando n�mero de frames
    datosegm.pixelsmierda=datosegm.videomedio_cuentaframes>.5;
else % M�todo tradicional
    
     % ATP ADDED THIS
    umbral = datosegm.umbral_LIGHT(1);
    %--------------------------------------
    
    datosegm.pixelsmierda=datosegm.videomedio<umbral;
end