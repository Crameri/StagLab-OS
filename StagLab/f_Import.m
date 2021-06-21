
%%                                                         IMPORT DATA 1.46
%
%                                                Fabio Crameri, 11.05.2021

function [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH)

%FUNCTION SETUP
DESVARIA.Task    = 'set special characters';
[~,~,STYLE,~] = f_DesignVaria(DESVARIA,[],[],[],[],[],[]);

%TASKS
if strcmp(DATA.Task,'ImportFieldData')
    [DATA,PLOT] = f_importFieldData(DATA,FILE,GRID,PLOT,SWITCH);
    
elseif strcmp(DATA.Task,'ImportTimedat')
    attemptedFileFix  	= false;
    if isfield(SWITCH,'enforceFixingTimedatFile') && SWITCH.enforceFixingTimedatFile
        [DATA] = f_fixTimedatFile(DATA); %try fixing time.dat file
        attemptedFileFix  	= true;
    end
    try
        [DATA] = f_importTimedat(DATA);
    catch
        [DATA] = f_fixTimedatFile(DATA); %try fixing time.dat file
        attemptedFileFix  	= true;
        [DATA] = f_importTimedat(DATA);
    end
    if attemptedFileFix 
        warning off backtrace; warning('Corrupt time.dat file has been fixed successfully.'); warning on backtrace; 
        disp(['                   ',STYLE.SCHAR.downrightArrow,' ',DATA.filename])
    end
    PLOT    = NaN;  %default output

elseif strcmp(DATA.Task,'ImportRefstat')
    [DATA] = f_importRefstat(FILE,DATA);
    PLOT    = NaN;  %default output

end

end



function [DATA,PLOT] = f_importFieldData(DATA,FILE,GRID,PLOT,SWITCH)
%%                                             IMPORT FIELD-DATA FILES 1.1
% Input:
%    DATA.Field2Import           = 'Density';
%    DATA.FieldAbbreviation      = 'RHO';
%                                                Fabio Crameri, 13.07.2017

%Defaults
DATA.NotFound           = false;
if ~isfield(DATA,'StopExecutionIfNotFound');    DATA.StopExecutionIfNotFound = false; end

%Check for and read data
if strcmp(DATA.Field2Import,'Velocity')  %VECTOR DATA
    
    if strcmp(GRID.Type,'yinyang')
        if ~strcmp(DATA.FieldAbbreviation,'VAR') && ...
                isfield(PLOT,[DATA.FieldAbbreviation,'X_3Dyin']) && isfield(PLOT,[DATA.FieldAbbreviation,'X_3Dyang']) && ...
                isfield(PLOT,[DATA.FieldAbbreviation,'Y_3Dyin']) && isfield(PLOT,[DATA.FieldAbbreviation,'Y_3Dyang']) && ...
                isfield(PLOT,[DATA.FieldAbbreviation,'Z_3Dyin']) && isfield(PLOT,[DATA.FieldAbbreviation,'Z_3Dyang'])
            %nothing to be done
        else
            [~,~,~,~,VARX_3D,VARY_3D,VARZ_3D,P_3D,VARX_3Dyang,VARY_3Dyang,VARZ_3Dyang,P_3Dyang,~,~] ...
                = f_readStagYY(FILE.directory,FILE.name,FILE.number,DATA.Field2Import,SWITCH);
            if ischar(VARX_3D)
                if DATA.StopExecutionIfNotFound
                   error(VARX_3D);
                end
                if SWITCH.Verbose
                    warning off backtrace
                    disp(' '); warning(VARX_3D); disp(' ');
                    warning on backtrace
                    DATA.NotFound           = true;
                end
                VARX_3D = ones(GRID.nx,GRID.ny,GRID.nz)*NaN; VARX_3Dyang = VARX_3D;
                VARY_3D = VARX_3D;      VARY_3Dyang = VARX_3D;  
                VARZ_3D = VARX_3D;      VARZ_3Dyang = VARX_3D;
                P_3D = VARX_3D;         P_3Dyang = VARX_3D;
            end
            eval(['PLOT.',DATA.FieldAbbreviation,'X_3Dyin   	= VARX_3D;'])
            eval(['PLOT.',DATA.FieldAbbreviation,'Y_3Dyin   	= VARY_3D;'])
            eval(['PLOT.',DATA.FieldAbbreviation,'Z_3Dyin   	= VARZ_3D;'])
            PLOT.P_3Dyin        = P_3D;
            eval(['PLOT.',DATA.FieldAbbreviation,'X_3Dyang   	= VARX_3Dyang;'])
            eval(['PLOT.',DATA.FieldAbbreviation,'Y_3Dyang   	= VARY_3Dyang;'])
            eval(['PLOT.',DATA.FieldAbbreviation,'Z_3Dyang   	= VARZ_3Dyang;'])
            PLOT.P_3Dyang   	= P_3Dyang;
        end
        eval(['PLOT.',DATA.FieldAbbreviation,'H_3Dyin 	= sqrt(PLOT.',DATA.FieldAbbreviation,'X_3Dyin.^2 +PLOT.',DATA.FieldAbbreviation,'Y_3Dyin.^2);'])   %absolute horizontal velocity
        eval(['PLOT.',DATA.FieldAbbreviation,'H_3Dyang = sqrt(PLOT.',DATA.FieldAbbreviation,'X_3Dyang.^2 +PLOT.',DATA.FieldAbbreviation,'Y_3Dyang.^2);'])   %absolute horizontal velocity
        eval(['PLOT.',DATA.FieldAbbreviation,'_3Dyin 	= sqrt(PLOT.',DATA.FieldAbbreviation,'X_3Dyin.^2 +PLOT.',DATA.FieldAbbreviation,'Y_3Dyin.^2 +PLOT.',DATA.FieldAbbreviation,'Z_3Dyin.^2);'])   %absolute velocity
        eval(['PLOT.',DATA.FieldAbbreviation,'_3Dyang 	= sqrt(PLOT.',DATA.FieldAbbreviation,'X_3Dyang.^2 +PLOT.',DATA.FieldAbbreviation,'Y_3Dyang.^2 +PLOT.',DATA.FieldAbbreviation,'Z_3Dyang.^2);'])   %absolute velocity
        %         VAR_3D          	= sqrt(VARX_3D.^2 +VARY_3D.^2 +VARZ_3D.^2);     %absolute velocity
%         VAR_3Dyang        	= sqrt(VARX_3Dyang.^2 +VARY_3Dyang.^2 +VARZ_3Dyang.^2);     %absolute velocity

    else %all other grid types
        if ~strcmp(DATA.FieldAbbreviation,'VAR') && ...
                isfield(PLOT,[DATA.FieldAbbreviation,'X_3D']) && isfield(PLOT,[DATA.FieldAbbreviation,'Y_3D']) && isfield(PLOT,[DATA.FieldAbbreviation,'Z_3D'])
            %nothing to be done
        else
            [~,~,~,~,VARX_3D,VARY_3D,VARZ_3D,P_3D,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,DATA.Field2Import,SWITCH);
            if ischar(VARX_3D)
                if DATA.StopExecutionIfNotFound
                   error(VARX_3D);
                end
                if SWITCH.Verbose
                    warning off backtrace
                    disp(' '); warning(VARX_3D); disp(' ');
                    warning on backtrace
                    DATA.NotFound           = true;
                end
                VARX_3D = ones(GRID.nx,GRID.ny,GRID.nz)*NaN;
                VARY_3D = VARX_3D; 
                VARZ_3D = VARX_3D;
                P_3D    = VARX_3D;
            end
            if size(VARX_3D,1)==1 %exchange x and y
                dummy_3D = zeros(size(VARX_3D,2),size(VARX_3D,1),size(VARX_3D,3));
                dummy_3D(:,1,:) = VARX_3D(1,:,:); 	VARX_3D = dummy_3D;
                dummy_3D(:,1,:) = VARY_3D(1,:,:);  	VARY_3D = dummy_3D;
                dummy_3D(:,1,:) = VARZ_3D(1,:,:);  	VARZ_3D = dummy_3D;
                dummy_3D = VARX_3D; VARX_3D = VARY_3D; VARY_3D = dummy_3D;
                dummy_3D(:,1,:) = P_3D(1,:,:);     	P_3D = dummy_3D;
            end
            if strcmp(GRID.Type,'spherical2D')
                VARX_3Ds     = -VARZ_3D.*sin(PLOT.X_3D)-VARX_3D.*cos(PLOT.X_3D);
                VARZ_3Ds     = VARX_3D.*sin(PLOT.X_3D)-VARZ_3D.*cos(PLOT.X_3D);
                eval(['PLOT.',DATA.FieldAbbreviation,'X_3Ds   	= VARX_3Ds;'])
                eval(['PLOT.',DATA.FieldAbbreviation,'Z_3Ds   	= VARZ_3Ds;'])
            end
            eval(['PLOT.',DATA.FieldAbbreviation,'X_3D   	= VARX_3D;'])
            eval(['PLOT.',DATA.FieldAbbreviation,'Y_3D   	= VARY_3D;'])
            eval(['PLOT.',DATA.FieldAbbreviation,'Z_3D   	= VARZ_3D;'])
            PLOT.P_3D           = P_3D;
        end
        eval(['PLOT.',DATA.FieldAbbreviation,'H_3D	= sqrt(PLOT.',DATA.FieldAbbreviation,'X_3D.^2 +PLOT.',DATA.FieldAbbreviation,'Y_3D.^2);'])   %absolute horizontal velocity
        eval(['PLOT.',DATA.FieldAbbreviation,'_3D	= sqrt(PLOT.',DATA.FieldAbbreviation,'X_3D.^2 +PLOT.',DATA.FieldAbbreviation,'Y_3D.^2 +PLOT.',DATA.FieldAbbreviation,'Z_3D.^2);'])   %absolute velocity
%         VAR_3D                = sqrt(VARX_3D.^2 + VARY_3D.^2 + VARZ_3D.^2);     %absolute velocity
    end
    
else %SCALAR DATA
    
    if strcmp(GRID.Type,'yinyang')
        if ~strcmp(DATA.FieldAbbreviation,'VAR') && ...
                isfield(PLOT,[DATA.FieldAbbreviation,'_3Dyin']) && isfield(PLOT,[DATA.FieldAbbreviation,'_3Dyang'])
            %nothing to be done
        else
            [~,~,~,~,VAR_3D,VAR_3Dyang,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,DATA.Field2Import,SWITCH);
            if ischar(VAR_3D)
                if DATA.StopExecutionIfNotFound
                   error(VAR_3D);
                end
                if SWITCH.Verbose
                    warning off backtrace
                    disp(' '); warning(VAR_3D); disp(' ');
                    warning on backtrace
                    DATA.NotFound           = true;
                end
                VAR_3D = ones(GRID.nx,GRID.ny,GRID.nz)*NaN; VAR_3Dyang = VAR_3D;
            end
            eval(['PLOT.',DATA.FieldAbbreviation,'_3Dyin   	= VAR_3D;'])
            eval(['PLOT.',DATA.FieldAbbreviation,'_3Dyang  	= VAR_3Dyang;'])
        end
        
    else %all other grid types
        if ~strcmp(DATA.FieldAbbreviation,'VAR') && ...
                isfield(PLOT,[DATA.FieldAbbreviation,'_3D'])
            %nothing to be done
        else
            [~,~,~,~,VAR_3D,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,DATA.Field2Import,SWITCH);
            if ischar(VAR_3D)
                if DATA.StopExecutionIfNotFound
                   error(VAR_3D);
                end
                if SWITCH.Verbose
                    warning off backtrace
                    disp(' '); warning(VAR_3D); disp(' ');
                    warning on backtrace
                    DATA.NotFound           = true;
                end
                VAR_3D = ones(GRID.nx,GRID.ny,GRID.nz)*NaN;
            end
            if size(VAR_3D,1)==1 %exchange x and y
                dummy_3D = zeros(size(VAR_3D,2),size(VAR_3D,1),size(VAR_3D,3));
                dummy_3D(:,1,:) = VAR_3D(1,:,:);    VAR_3D = dummy_3D;
            end
            eval(['PLOT.',DATA.FieldAbbreviation,'_3D      = VAR_3D;'])
        end
    end
end

%Remove fields
DATA = rmfield(DATA,'StopExecutionIfNotFound');

end



function [DATA] = f_importTimedat(DATA)
%%                                       IMPORT STAGYY'S TIMEDAT FILES 4.0
%Input:
% DATA.filename
% DATA.startRow  (optional)
% DATA.endRow    (optional)
% DATA.titleRow  (optional)
%                                                Fabio Crameri, 10.09.2018

%IMPORTFILE Import numeric data from a text file as column vectors.
%  Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.

%% Initialize variables.
if ~isfield(DATA,'titleRow');           DATA.titleRow     = true;    end    %data file includes header
delimiter           = ' ';
% DATA.startRow       = 2;
DATA.endRow         = inf;

%% Format string for each line of text:
% For more information, see the TEXTSCAN documentation.
formatSpecTitle = ''; formatSpec = '';
for iEntry=1:DATA.numberEntries
    formatSpecTitle     = [formatSpecTitle,'%s'];
    formatSpec          = [formatSpec,'%f'];
end
formatSpec              = [formatSpec,'%[^\n\r]'];

%% Open the text file.
fileID = fopen(DATA.filename,'r');

if DATA.titleRow
    DATA.TitleStartRow  = 1;
    DATA.TitleEndRow    = 1;
    %% Read columns of title according to format string.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    titleArray = textscan(fileID,formatSpecTitle,DATA.TitleEndRow(1)-DATA.TitleStartRow(1)+1,'Delimiter',delimiter,'MultipleDelimsAsOne',true,'ReturnOnError',false);
    
    % Title comparison check
    if strcmp(DATA.title(1,end),titleArray{1,end}) && ...
            size(titleArray,2)==size(DATA.title,2)
    else
        warning('!!! check new format of data file !!!');
    end
end

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
DATA.Array = textscan(fileID,formatSpec,DATA.endRow(1)-DATA.startRow(1)+1,'Delimiter',delimiter,'MultipleDelimsAsOne',true,'HeaderLines',DATA.startRow(1)-1,'ReturnOnError',false);
for block=2:length(DATA.startRow)
    frewind(fileID);
    DATA.ArrayBlock = textscan(fileID,formatSpec,DATA.endRow(block)-DATA.startRow(block)+1,'Delimiter',delimiter,'MultipleDelimsAsOne',true,'HeaderLines',DATA.startRow(block)-1,'ReturnOnError',false);
    for col=1:length(DATA.Array)
        DATA.Array{col} = [DATA.Array{col};DATA.ArrayBlock{col}];
    end
end

%remove last obsolete entry
if iscell(DATA.Array{end})
    DATA.Array(end)     = [];
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
% istep = DATA.Array{:, 1};
% time = DATA.Array{:, 2};
% F_top = DATA.Array{:, 3};
% F_bot = DATA.Array{:, 4};
% Tmin = DATA.Array{:, 5};
% Tmean = DATA.Array{:, 6};
% Tmax = DATA.Array{:, 7};
% Vmin = DATA.Array{:, 8};
% Vrms = DATA.Array{:, 9};
% Vmax = DATA.Array{:, 10};
% eta_min = DATA.Array{:, 11};
% eta_mean = DATA.Array{:, 12};
% eta_max = DATA.Array{:, 13};
% ra_eff = DATA.Array{:, 14};
% Nu_top = DATA.Array{:, 15};
% Nu_bot = DATA.Array{:, 16};
% C_min = DATA.Array{:, 17};
% C_mean = DATA.Array{:, 18};
% C_max = DATA.Array{:, 19};
% F_mean = DATA.Array{:, 20};
% F_max = DATA.Array{:, 21};
% erupt_rate = DATA.Array{:, 22};
% erupta = DATA.Array{:, 23};
% erupt_heatflux = DATA.Array{:, 24};
% entrainment = DATA.Array{:, 25};
% Cmass_error = DATA.Array{:, 26};
% H_int = DATA.Array{:, 27};

% r_innercore = DATA.Array{:, 28}; 
% Tsurf = DATA.Array{:, 29};  
% Tcmb = DATA.Array{:, 30};

%% Extract data needed
% checkIndivRange = [0 DATA.individual_range(1,:) 0];
% % maybe more efficient using this structure & using intercept: DATA.Array{strcmp(DATA.title,'time')}(2,1);
% for i=1:size(DATA.title,2)
%     if ~max(strcmp(DATA.title(1,i),DATA.tag))==1 && ~strcmp(DATA.title(1,i),'time') && ...
%             ~checkIndivRange(1,i+2) && ~checkIndivRange(1,i)  %check if not min or max data
%
%         %delete unused rows
%         DATA.Array{:,i} = [];
%     end
% end

%adjust DATA.Array length  (not sure why this is needed...)
% DATA.Array(28) = [];
end


function [DATA] = f_fixTimedatFile(DATA)
%%                                        FIXING CORRUPT TIMEDAT FILES 1.0

%                                                Fabio Crameri, 10.09.2018
corruptFile     = DATA.filename;
cleanFile       = strrep(corruptFile,'.dat','Fixed.dat');

lines           = strsplit(fileread(corruptFile),'\n');  %read file and split into lines
lines           = unique(lines,'stable');  %remove duplicate lines

fid             = fopen(cleanFile,'w');    %open new file for writing
timestep = zeros(size(lines,2)-2,1)*NaN;
ProblemFound = false;
ilnew = 1;
linesNew(ilnew)             = lines(1); %copy first line
for il=2:size(lines,2)-1
    lineCell                = lines(il);
    lineStr                 = lines{il};
    if ~ProblemFound && il>2
        timestepPrevious    = timestep(il-2,1);
    end
    timestepCurrent         = sscanf(lineStr,'%f',1);
    if  il>2 && timestepCurrent~=timestepPrevious+1
        ProblemFound        = true;
        timestep(il-1,1)    = NaN;
    else
        ilnew               = ilnew+1;
        ProblemFound        = false;
        timestep(il-1,1)    = timestepCurrent;
        
        linesNew(ilnew)    	= lines(il);
    end
end
%write to new file
fwrite(fid, strjoin(linesNew,'\n'),'char'); %merge lines and write
fclose(fid);
%function output
DATA.filename   = cleanFile;
end



function [PARA] = f_importRefstat(FILE,DATA)
%%                                       IMPORT STAGYY'S REFSTAT FILES 1.2

%                                                Fabio Crameri, 21.05.2020
if isfield(DATA,'nz')
    onlyFirstEntries = false;
    nz = DATA.nz;
    nz = nz*2; %double as many points as vertical grid points
else
    onlyFirstEntries = true;
    nz              = 0;    %reads only one line of data
end
filestem        = [FILE.directory,FILE.name,'_refstat.dat'];

%% Check title format
fileID          = fopen(filestem,'r');
entry         	= fscanf(fileID,'%s[^\n\r]',[1 1]);
fclose(fileID);
if strcmp(entry,'SYSTEM') %new StagYY version
    startRow0   = 4;
else %old StagYY version
    startRow0   = 3;
end

%% Find "COMBINED ADIABAT" line
fileID = fopen(filestem,'r');
numberLines       	= 0;
cc                	= 0;
tline           	= fgetl(fileID);
while cc<4 && ~strcmp(num2str(tline),'-1')  %tline becomes -1 when file ended during the loop
    tline           = fgetl(fileID);
    numberLines     = numberLines+1;
    if ~isempty(strfind(tline,'COMBINED ADIABAT'))
        %if contains(tline,'COMBINED ADIABAT') %not compatible with MatLab versions older than 2016
        cc  = cc+1;
        startRow_CA = numberLines+2;
        endRow_CA   = startRow_CA+nz;
    end
end
fclose(fileID);
if onlyFirstEntries
    numPhases   = 1;
else
    %% Find number of phases (i.e., column packages)
    numPhases   = round((startRow_CA-startRow0)/nz);
    numPhases   = numPhases+1; %adding the combined adiabate package
end

%% Read data
numPhases=2; %<<<<<<<<<<<<<<<< prevent reading all phases - needs further implementation
for icol=1:numPhases
    %% Initialize variables (find start and end row)
    delimiter       = ' ';
    if icol==1 
        startRow    = startRow0;
        endRow      = startRow+nz;
    elseif icol<numPhases
        startRow    = startRow+nz+2;  %%% NEEDS CHECK WHETHER THERE IS INDEED ONLY ONE HEADER ABOVE EACH SUBSEQUENT PACKAGE...
        endRow      = endRow+nz+2;
    elseif icol==numPhases %combined adiabat
        startRow    = startRow_CA;
        endRow      = endRow_CA;
    end

    %% Format string for each line of text:
    formatSpec      = '%f%f%f%f%f%f%f%*s%[^\n\r]';
    
    %% Open the text file.
    fileID = fopen(filestem,'r');
    
    %% Read columns of data according to format string.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    for block=2:length(startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end
    
    %% Close the text file.
    fclose(fileID);
    
    %% Allocate imported array to column variable names
    if icol>1 && icol==numPhases %last (combined adiabat) columns
        PARA.CA_z(:,1)          = dataArray{:, 1};
        PARA.CA_Tref(:,1)   	= dataArray{:, 2};
        PARA.CA_rho(:,1)        = dataArray{:, 3};
        PARA.CA_expan(:,1)      = dataArray{:, 4};
        PARA.CA_Tcond(:,1)  	= dataArray{:, 5};
        PARA.CA_P(:,1)          = dataArray{:, 6};
    else
        PARA.z(:,icol)          = dataArray{:, 1};
        PARA.Tref(:,icol)       = dataArray{:, 2};
        PARA.rho(:,icol)        = dataArray{:, 3};
        PARA.expan(:,icol)      = dataArray{:, 4};
        PARA.Cp(:,icol)         = dataArray{:, 5};
        PARA.Tkappa(:,icol)     = dataArray{:, 6};
        PARA.Tcond(:,icol)      = dataArray{:, 7};
    end
end

end






