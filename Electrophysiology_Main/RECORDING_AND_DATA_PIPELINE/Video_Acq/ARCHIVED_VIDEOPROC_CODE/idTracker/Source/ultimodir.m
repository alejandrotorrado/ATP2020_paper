% APE 4 dic 13

% (C) 2014 Alfonso P�rez Escudero, Gonzalo G. de Polavieja, Consejo Superior de Investigaciones Cient�ficas

function directorio=ultimodir(directorio)

if nargin==1 && isstr(directorio) % Si tiene input, es que hay que guardar el directorio
    save UltimoDir directorio    
else % Si no tiene input, es que hay que recuperar el directorio
    if ~isempty(dir('.\UltimoDir.mat'))
        load UltimoDir
    else
        directorio='';
    end
end