
%%                                                          SAVE DATA  2.02
%
%                                                Fabio Crameri, 21.06.2019

function [overwriteAll] = f_saveData( SAVE )

% SAVE.Directory              = '~/Documents/PLOTS/';
% SAVE.DataName               = 'test';
% SAVE.data                   = [1 2 3; 4 5 6];
% SAVE.dat                    = logical(0);
% SAVE.txt                    = logical(0);
% SAVE.mat                    = logical(1);
% SAVE.write2existing         = logical(0);
% SAVE.overwrite              = logical(0);    %this is a global variable!
% SAVE.Pause                  = logical(0);

% [SAVE.overwriteAll] = f_saveData( SAVE );

%**************************************************************************
%% DEFAULTS
SAVE.dummy = 0;
if ~isfield(SAVE,'Directory');          SAVE.Directory      = ['~',filesep,'Desktop',filesep,'SLData',filesep];	end
if ~isfield(SAVE,'DataName');           SAVE.DataName               = 'test';                   end
if ~isfield(SAVE,'DataTitle');          SAVE.DataTitle              = {'none'};                 end
if ~isfield(SAVE,'data');               SAVE.data                   = [1 2 3; 4 5 6];           end
if ~isfield(SAVE,'dat');                SAVE.dat                    = logical(0);               end
if ~isfield(SAVE,'txt');                SAVE.txt                    = logical(0);               end
if ~isfield(SAVE,'mat');                SAVE.mat                	= logical(1);               end
if ~isfield(SAVE,'write2existing');  	SAVE.write2existing      	= logical(0);               end
if ~isfield(SAVE,'overwrite');          SAVE.overwrite              = logical(0);               end
if ~isfield(SAVE,'Pause');              SAVE.Pause                  = logical(0);               end
if ~isfield(SAVE,'overwriteAll');       SAVE.overwriteAll           = SAVE.overwrite;           end
if ~isfield(SAVE,'BulletString');    	SAVE.BulletString       	= char(9862);               end
if ~isfield(SAVE,'PointingString');   	SAVE.PointingString       	= char(10551);              end
overwriteAll    = SAVE.overwriteAll;

%% SAVING DATA
if SAVE.Pause; pause; end
currentDir = pwd;

% CHECK IF DIRECTORY ALREADY EXISTS
if ~exist(SAVE.Directory,'dir')
    c=questdlg(['Create directory "',SAVE.Directory,'"?'],'Directory does not exist!', ...
        'Yes','No','Yes');
    switch c
        case 'Yes'
            %create directory
            mkdir(SAVE.Directory);
        case 'No'
            error('So, change save-directory and re-run again.')
    end
end

if ~isempty(SAVE.Directory)
    cd(SAVE.Directory)
end

%check for existing +data folder
saveDirectory = [SAVE.Directory,'+data/'];
if ~exist(saveDirectory,'dir'); mkdir(saveDirectory); end
cd(saveDirectory)

%check for existing data files
skipSaving = false;
if SAVE.write2existing  %only works with .mat files
    if exist([saveDirectory,SAVE.DataName,'.mat'],'file')
        load([saveDirectory,SAVE.DataName,'.mat']); %load previously saved 'saveData' from .mat file
        SAVE.data = [saveData; SAVE.data];  %Adds new data as a new row  %NEEDS TESTING........
    end
    
elseif ~SAVE.overwrite && ~overwriteAll
    bull(1) = exist([saveDirectory,SAVE.DataName,'.dat'],'file'); %compare
    bull(2) = exist([saveDirectory,SAVE.DataName,'.mat'],'file'); %compare
    if max(bull)
        b=questdlg(['Overwrite "',SAVE.DataName,'"?'],'Data file allready exists!', ...
            'Yes','Yes to all','No','Yes');
        switch b
            case 'Yes'
                %continue
            case 'Yes to all'
                overwriteAll = true;
                %             case 'Keep both'
                %                 %change filename and keep both
                %                 for i_fnumber=1:999
                %                     if i_fnumber==999; skipSaving = true; disp('file could not be saved: numbering only goes to "_99".'); end
                %                     num_string = num2str(100+i_fnumber);
                %                     num_string(1)='';
                %                     DataName2 = [SAVE.DataName,'_',num_string]; %update filename
                %                     if ~exist([DataName2,'.dat'],'file') && ...
                %                             ~exist([DataName2,'.mat'],'file')
                %                         SAVE.DataName = DataName2;
                %                         break
                %                     end
                %                 end
            case 'No'
                skipSaving = true;
        end
    end
end

if skipSaving
    %skip saving this data
    
else
    %% SAVE DATA TO A FILE
    % save the data matrix 'SAVE.data' as a text file
    % to the current directory
    saveData        = SAVE.data;
    saveDataTitle   = SAVE.DataTitle;
    if isfield(SAVE,'DataTitle') && ~strcmp(SAVE.DataTitle(1),'none') %with title line(s)
        if SAVE.dat || SAVE.txt
            formatSpecTitle = ''; formatSpec = '';
            for iEntry=1:size(SAVE.DataTitle,1)
                if iEntry==1;	formatSpecTitle = [formatSpecTitle,'%s'];
                else;           formatSpecTitle = [formatSpecTitle,' %s'];
                end
            end
            formatSpecTitle   	= [formatSpecTitle,'\n'];
            for iEntry=1:size(saveData,2)
                if iEntry==1;   formatSpec = [formatSpec,'%f'];
                else;           formatSpec = [formatSpec,' %f'];
                end
            end
            formatSpec       	= [formatSpec,'\n'];
            if SAVE.dat
                fid             = fopen([SAVE.DataName,'.dat'],'w');
            elseif SAVE.txt
                fid             = fopen([SAVE.DataName,'.txt'],'w');
            end
            fprintf(fid,formatSpecTitle, SAVE.DataTitle');
            fprintf(fid,formatSpec, saveData');
            fclose(fid);
        end
        if SAVE.mat
            saveDataTitle       = strsplit(saveDataTitle,', ');
            save([SAVE.DataName,'.mat'],'saveData','saveDataTitle')
        end
        
    else %without title line
        if SAVE.dat
            save([SAVE.DataName,'.dat'],'saveData','-ascii')
            
        end
        if SAVE.txt
            save([SAVE.DataName,'.txt'],'saveData','-ascii')
            
        end
        if SAVE.mat
            save([SAVE.DataName,'.mat'],'saveData')
            
        end
    end
    
    if ~isempty(SAVE.Directory)
        cd ..
    end
    
    %% DISPLAY INFORMATION
    disp(' ')
    disp(['   ',SAVE.BulletString,' Data Saved'])
    disp(['     ',SAVE.PointingString,SAVE.Directory,'+data',filesep,SAVE.DataName])
    
end
cd(currentDir)



















