function klustakwik_matlab_hpc(command)

%   KLUSTAKWIK_MATLAB_HPC - Call KlustaKwik from Matlab
%
%   KLUSTAKWIK_MATLAB(COMMAND_STRING)
%
%  This function calls KlustaKwik from Matlab with the string
%  argument COMMAND_STRING.  It assumes there is a
%  KlustaKwik directory called 'KlustaKwik' located one level up
%  from the directory where KlustaKwik_matlab is located (type
%  which KlustaKwik_matlab to learn this).
%  If COMMAND_STRING is not given, KlustaKwik is called with no
%  arguments.
%
% This version is setup to run commands on Linux and specifically using the
% copy of KlustaKwik compiled on ATP's home folder on the HPC cluster.


if nargin==0, command = ''; end;

% THIS IS WHERE YOU SET THE LOCATION OF THE KLUSTAKWIK EXECUTABLE FILE
kexecute= '/home/atorrpac/bin/KlustaKwik';
        
command_string = [kexecute ' ' command ' < /dev/null'];
tic
[status, ~ ]  = system(command_string);
% disp(status)
toc
