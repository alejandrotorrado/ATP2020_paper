function klustakwik_matlab(command)

%   KLUSTAKWIK_MATLAB - Call KlustaKwik from Matlab
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
%  If we are running on a PC, KLUSTAKWIK_MATLAB attempts to call
%  [ONELEVELUP]/KlustaKwik/windows/KlustaKwik.exe
%
%  If we are running on a Mac, KLUSTAKWIK_MATLAB attempts to call
%  [ONELEVELUP]/KlustaKwik/macosx-i386/KlustaKwik
%
%  If we are running on Linux, KLUSTAKWIK_MATLAB attempts to call
%  [ONELEVELUP]/KlustaKwik/linux-pentium/KlustaKwik
%
%  The KlustaKwik executable is launched through the '!' command
%  line function (see 'help !'). This assumes the user does not need
%  to interact with the KlustaKwik process, although any text printed
%  to the console will be displayed.


global myname;

if isempty(myname);
    startup;
end

disp(myname)



if nargin==0, command = ''; end;

klusta_filename = which('klustakwik_matlab');

[klustapath_root,filename]  = fileparts(klusta_filename);
[klustapath_onelevelup]     = fileparts(klustapath_root);
klustapath = [klustapath_onelevelup filesep 'KlustaKwik'];

arch = computer;

if ~isempty(strfind(computer,'MAC')),
    kexecute = [klustapath filesep 'macosx-i386' filesep 'KlustaKwik']; %updated to KKWIK v2.0 on 1/4/14 KBH
elseif ~isempty(strfind(computer,'PC')),
    
    
    switch (myname);
        case ('iDell');
            kexecute= [klustapath filesep 'windows' filesep 'KlustaKwik_Latest.exe'];
        case ('coypu');
            kexecute= [klustapath filesep 'windows' filesep 'KlustaKwik_Latest.exe'];
        case ('badger');
            kexecute= [klustapath filesep 'windows' filesep 'KlustaKwik.exe'];
        case ('marmoset');
            kexecute= [klustapath filesep 'windows' filesep 'KlustaKwik_Latest.exe'];
    end
elseif ~isempty(strfind(computer,'LNX')),
    kexecute= [klustapath filesep 'linux-pentium' filesep 'KlustaKwik'];
else
    error(['Sorry, platform ' computer ' is not supported.']);
end;
command_string = ['!' kexecute ' ' command];
tic
eval(command_string);
toc
