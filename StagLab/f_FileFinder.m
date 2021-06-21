
%%                                                         FILE FINDER 5.67
% input:
%  FILE.name
%  FILE.number (optional)
%  FILE.stemRead
%
% can check automatically for alternative folder structures as e.g.,:
%  /folder/+op/<filename>...
%                                                Fabio Crameri, 31.05.2021

%% NOTES
% searches only for a few fields/files (in SL_FieldPlot mode).

function [FILE] = f_FileFinder(FILE,SWITCH,SAVE,STYLE)
%% DEFAULTS
if ~exist('SWITCH','var') || ~isfield(SWITCH,'GeodynamicCode'); SWITCH.GeodynamicCode   = 'StagYY';     end
if ~exist('SAVE','var') || ~isfield(SAVE,'app');                SAVE.app                = 'na';         end
if ~exist('STYLE','var') || ~isfield(STYLE,'SCHAR') || ~isfield(STYLE.SCHAR,'doubleSidedArrow'); STYLE.SCHAR.doubleSidedArrow = '<->'; end
if ~isfield(SWITCH,'StagYYoutput');         SWITCH.StagYYoutput     = 'Binary'; 	end
if ~isfield(SWITCH,'Verbose');              SWITCH.Verbose          = false;      	end
if ~isfield(SWITCH,'checkForFolder');       SWITCH.checkForFolder   = true;     	end
if ~isfield(FILE,'number');                 FILE.number             = NaN;       	end
FILE.NotFound               = false;
FILE.NotFoundError          = '';
noFileInFolder              = false;

%% FUNCTION SETUP
if strcmp(SAVE.app,'SL_FieldPlot')
    if strcmpi(SWITCH.GeodynamicCode,'StagYY') && strcmpi(SWITCH.StagYYoutput,'HDF5') %hdf5
        fileNameExtensions  	= {'.h5','.h5','.h5','.h5'};
        fieldNames            	= {'Temperature','Viscosity','Stress','Density'};
        fieldName              	= fieldNames{1,1};
        fieldNatures         	= {'field','field','field','field'};
        fieldNature          	= fieldNatures{1,1};
        
    elseif strcmpi(SWITCH.GeodynamicCode,'StagYY') %binary
        fileNameExtensions  	= {'_t0*','_t1*','_t2*','_t3*','_t4*','_t5*','_t6*','_t7*','_t8*','_t9*',...
            '_eta*','_ed*','_str*','_rho*','_cs*','_sage*'}; %using only graph data, diagnostics like grid spacing will be messed up!
        fieldNames            	= {'Temperature','Temperature','Temperature','Temperature','Temperature',...
            'Temperature','Temperature','Temperature','Temperature','Temperature',...
            'Viscosity','Strain rate','Stress','Density','Topography','Surface age'};
        fieldName              	= fieldNames{1,1};
        fieldNatures          	= {'field','field','field','field','field','field','field','field','field','field',...
            'field','field','field','field','graph','graph'};
        fieldNature          	= fieldNatures{1,1};
        
    elseif strcmpi(SWITCH.GeodynamicCode,'Fluidity')
        fileNameExtensions   	= {'.csv'};
        fieldNames            	= {'Temperature'};
        fieldName              	= fieldNames{1,1};
        fieldNatures          	= {'field'};
        fieldNature          	= fieldNatures{1,1};
%         SWITCH.checkForFolder   = false;

    elseif strcmpi(SWITCH.GeodynamicCode,'Aspect') && strcmpi(SWITCH.ASPECToutput,'HDF5') %hdf5
        fileNameExtensions  	= {'.xdmf'};
        fieldNames            	= {'Temperature'};
        fieldName              	= fieldNames{1,1};
        fieldNatures         	= {'field','field','field'};
        fieldNature          	= fieldNatures{1,1};
        
    else
        error('Adjust here to use StagLab with output files from another geodynamic code.')
    end
elseif strcmp(SAVE.app,'SL_RadialProfile')
    fileNameExtensions    	= {'_rprof.dat'};
elseif strcmp(SAVE.app,'SL_TimeGraph')
    fileNameExtensions    	= {'_time.dat'};
else
    error('Specify which app is used here to set fileFinder mode!')
end
fileNameExtension           = fileNameExtensions{1,1};

%% PROBLEM CHECKS
%check for directory separator at the end of string
if ~strcmp(FILE.stemRead(1,end),filesep)
    FILE.stemRead           = [FILE.stemRead,filesep];
end
if ~strcmp(FILE.stemSave(1,end),filesep)
    FILE.stemSave           = [FILE.stemSave,filesep];
end

%% FOLDER STRUCTURE
if SWITCH.checkForFolder
    %% CHECK CURRENT FOLDER
    stringToFind                    = [FILE.name,fileNameExtensions{1,1}];  %check for string of the file
    filesFound                    	= dir([FILE.stemRead,stringToFind]);
    filesFound([filesFound.isdir])	= [];  %remove directories
    
    %% CHECK FOR ALTERNATIVE FIELDS
    if isempty(filesFound) && strcmp(SAVE.app,'SL_FieldPlot') %file not found in current directory
        if strcmpi(SWITCH.GeodynamicCode,'StagYY') && strcmpi(SWITCH.StagYYoutput,'HDF5') %hdf5
            for iFields=1:size(fileNameExtensions,2)
                stringToFind                    = [fieldNames{1,iFields},'_*',fileNameExtensions{1,iFields}];  %check for string of the file
                filesFound                    	= dir([FILE.stemRead,stringToFind]);
                filesFound([filesFound.isdir])	= [];  %remove directories
                if ~isempty(filesFound) %file not found in current directory
                    fileNameExtension           = fileNameExtensions{1,iFields};
                    fieldName                   = fieldNames{1,iFields};
                    fieldNature                 = fieldNatures{1,iFields};
                    break
                end
            end
            
        elseif strcmpi(SWITCH.GeodynamicCode,'StagYY') %binary
            for iFields=1:size(fileNameExtensions,2)
                stringToFind                    = [FILE.name,fileNameExtensions{1,iFields}];  %check for string of the file
                filesFound                    	= dir([FILE.stemRead,stringToFind]);
                filesFound([filesFound.isdir])	= [];  %remove directories
                if ~isempty(filesFound) %file not found in current directory
                    fileNameExtension           = fileNameExtensions{1,iFields};
                    fieldName                   = fieldNames{1,iFields};
                    fieldNature                 = fieldNatures{1,iFields};
                    break
                end
            end
            
        elseif strcmpi(SWITCH.GeodynamicCode,'Fluidity')
            %nothing to do!
            
        elseif strcmpi(SWITCH.GeodynamicCode,'Aspect') && strcmpi(SWITCH.ASPECToutput,'HDF5') %hdf5
            for iFields=1:size(fileNameExtensions,2)
                stringToFind                    = 'solution.xdmf';  %check for string of the file
                filesFound                    	= dir([FILE.stemRead,stringToFind]);
                filesFound([filesFound.isdir])	= [];  %remove directories
                if ~isempty(filesFound) %file not found in current directory
                    fileNameExtension           = fileNameExtensions{1,iFields};
                    fieldName                   = fieldNames{1,iFields};
                    fieldNature                 = fieldNatures{1,iFields};
                    break
                end
            end
            
            if SWITCH.Verbose; warning('could implement this better: Check here!'); end
            
        else
            error('Adjust here to use StagLab with output files from another geodynamic code.')
        end
    end
    
    %% CHECK FOR ALTERNATIVE FOLDERS
    if isempty(filesFound) %file not found in current directory
        % INPUT-FILE FOLDERS TO CHECK
        ReadingDir2check        	= { strcat(FILE.stemRead,'+op',filesep,FILE.name,filesep); %list deeper folders first
                                        strcat(FILE.stemRead,FILE.name,filesep) };
        % OUTPUT-FILE FOLDERS TO CHECK
        WritingDir2check        	= { strcat(FILE.stemRead,'+im',filesep,FILE.name,filesep) }; %list deeper folders first
        
        % CHECK FOR INPUT-FILES FOLDER
        for idir=1:size(ReadingDir2check,1)
            filesFound            	= false;
            dir2check               = ReadingDir2check{idir};
            for iFields=1:size(fileNameExtensions,2)
                stringToFind      	= [FILE.name,fileNameExtensions{1,iFields}];  %check for string of the file
                if exist(dir2check,'dir') && ~isempty(dir([dir2check,stringToFind]))
                    FILE.stemRead 	= dir2check;
                    filesFound    	= true;
                end
                if filesFound
                    break
                end
            end
            if filesFound
                fileNameExtension   = fileNameExtensions{1,iFields};
                if strcmp(SAVE.app,'SL_FieldPlot')
                    fieldName         	= fieldNames{1,iFields};
                    fieldNature     	= fieldNatures{1,iFields};
                end
                break
            end
        end
        % CHECK FOR OUTPUT-FILES FOLDER (+im)
        for idir=1:size(WritingDir2check,1)
            directoryFound          = false;            
            dir2check               = WritingDir2check{idir};
            if exist(dir2check,'dir')
                FILE.stemSave       = dir2check;
                directoryFound    	= true;
            end
            if directoryFound
                break
            end
        end
    else %file found
        %nothing to do: continue
    end
end

%% REMOVE ADDED NUMBER TO FILE-NAME-EXTENSION
if strcmpi(SWITCH.GeodynamicCode,'StagYY') && ~strcmpi(SWITCH.StagYYoutput,'HDF5') %StagYY binary only
    if strcmp(fileNameExtension(1,end-1),'1') || strcmp(fileNameExtension(1,end-1),'2') || ...
            strcmp(fileNameExtension(1,end-1),'3') || strcmp(fileNameExtension(1,end-1),'4') ||...
            strcmp(fileNameExtension(1,end-1),'5') || strcmp(fileNameExtension(1,end-1),'6') ||...
            strcmp(fileNameExtension(1,end-1),'7') || strcmp(fileNameExtension(1,end-1),'8') ||...
            strcmp(fileNameExtension(1,end-1),'9') || strcmp(fileNameExtension(1,end-1),'0')
        fileNameExtension = [fileNameExtension(1,1:end-2),fileNameExtension(1,end)];
    end
end

%% CHECK FOR NUMBERED INPUT FILE
if strcmp(SAVE.app,'SL_FieldPlot')
    FILE.LastFile = false;
    if strcmpi(SWITCH.GeodynamicCode,'Fluidity')
        numberingFormat     = '';
        fnumberString       = '';
        filetofind          = [FILE.stemRead,FILE.name,fileNameExtension];
        if ~exist(filetofind,'file')
            noFileInFolder 	= true;
        end
        
    elseif strcmpi(SWITCH.GeodynamicCode,'StagYY') && strcmpi(SWITCH.StagYYoutput,'HDF5') %hdf5
        numberingFormat     = '';
        fnumberString       = '';
        filetofind          = [FILE.stemRead,fieldName,'_*',fileNameExtension];
        if isempty(dir(filetofind))
            noFileInFolder 	= true;
        end
        
    elseif strcmpi(SWITCH.GeodynamicCode,'Aspect') && strcmpi(SWITCH.ASPECToutput,'HDF5') %hdf5
        numberingFormat     = '';
        fnumberString       = '';
        filetofind          = [FILE.stemRead,'solution/solution-00000.h5'];
        if isempty(dir(filetofind))
            noFileInFolder 	= true;
        end
        
    else
        %check for current file
        if strcmpi(SWITCH.GeodynamicCode,'StagYY') %binary
            numberingFormat = '00000';
            filetofindOrig 	= [FILE.stemRead,FILE.name,fileNameExtension(1,1:end-1),numberingFormat]; %full file name for zero numbered Temperature file
            % elseif strcmpi(SWITCH.GeodynamicCode,'->CODENAME<-')
            %filetofindOrig  = [FILE.stemRead,FILE.name,'/',FILE.name,'_t00000']; %full file name for zero numbered Temperature file
        else
            error('Adjust here to use StagLab with output files from another geodynamic code.')
        end
        filetofind          = filetofindOrig;
        filetofind2         = filetofindOrig;
        nextfiletofind      = filetofindOrig;
        
        fnumberString       = num2str(FILE.number);
        if ~strcmp(filetofind(1,end-size(fnumberString,2)+1),'0'); error('Filenumber is too big!'); end
        filetofind(1,end-size(fnumberString,2)+1:end) = fnumberString;
        
        %already check for next file
        fnumberString2      = num2str(FILE.number+1);
        if ~strcmp(nextfiletofind(1,end-size(fnumberString2,2)+1),'0'); error('Filenumber is too big!'); end
        nextfiletofind(1,end-size(fnumberString2,2)+1:end) = fnumberString2;
        
        %% CONVERT FOR WINDOWS PC
        if ispc %just to be sure
            filetofindOrig   = strrep(filetofind,'/','\');
            filetofind2      = strrep(filetofind,'/','\');
            filetofind       = strrep(filetofind,'/','\');
            nextfiletofind   = strrep(filetofind,'/','\');
        end
        
        %% IF FILE NOT FOUND CHECK FOR LATEST FILE
        %list all files
        list1 = [];
        for inum=0:9 %loop largest number to make sure it takes the temperature file
            list0           = dir([FILE.stemRead,FILE.name,fileNameExtension(1,1:end-1),num2str(inum),'*']);
            list1           = [list1;list0];
        end
        nameList            = {list1.name};
        if isempty(nameList)
            noFileInFolder  = true;
        else
            nameList(strcmp(nameList(1,:),[FILE.name,'_time.dat'])) = [];    %remove time dat file to not be confused with '_t00...'
            %get first output number
            sortedList   	= sort(nameList);
            firstName   	= sortedList{1};
            dummy           = regexp(firstName,fileNameExtension(1,1:end-1),'split');
            firstNumberStr  = dummy{end};
            firstNumber     = str2double(firstNumberStr);
            %get last output number
            latestName   	= sortedList{length(sortedList)};
            dummy           = regexp(latestName,fileNameExtension(1,1:end-1),'split');
            lastNumberStr   = dummy{end};
            lastNumber      = str2double(lastNumberStr);
            
            %FILE.LastFile = false; %may be removed if not necessary
            if ~exist(filetofind,'file')
                if FILE.number<firstNumber %check for first output file number
                    filetofind2(1,end-size(firstNumberStr,2)+1:end) = firstNumberStr;
                    if exist(filetofind2,'file')   % THIS IS SLOW
                        FILE.number     = firstNumber;
                    end
                    warning off backtrace
                    warning(['Specified file number is too low: First available file number is ',firstNumberStr,'!'])
                    warning on backtrace
                else %check for last output file number
                    filetofind2(1,end-size(lastNumberStr,2)+1:end) = lastNumberStr;
                    if exist(filetofind2,'file')   % THIS IS SLOW
                        FILE.number     = lastNumber;
                        FILE.LastFile   = true; %last file that exists reached
                    end
                end
            elseif ~exist(nextfiletofind,'file')
                FILE.LastFile = true; %last file that exists reached
            end
        end
        fnumberString       = num2str(FILE.number); %update file number
    end
    
elseif strcmp(SAVE.app,'SL_RadialProfile') || strcmp(SAVE.app,'SL_TimeGraph')
    filetofind          = [FILE.stemRead,FILE.name,fileNameExtension];
    if ~exist(filetofind,'file')
        noFileInFolder 	= true;
    end
end

%% ERROR MESSAGE IF NO FILE AT ALL IS NOT FOUND
errorMessage = '';
if FILE.number==-1 || noFileInFolder
    folderParts     = regexp(filetofind,filesep,'split');
    if ~strcmp(folderParts{1},'.') && ~strcmp(folderParts{1},'..') && ~strcmp(folderParts{1},'~') && ~ispc
        prefix      = '/';
    else
        prefix      = '';
    end
    for idir=1:size(folderParts,2)
        if idir==1 && isempty(folderParts{1})
            continue
        end
        if idir==size(folderParts,2)
            %errorMessage    = [newline STYLE.SCHAR.problemWarning,' FILE NOT FOUND!' newline '   ',STYLE.SCHAR.downrightArrow,' ',folderParts{end},'  in ',[prefix,fullfile(folderParts{1:end-1})] newline]; %newline is only available since 2016b
            errorMessage    = [char(13) STYLE.SCHAR.problemWarning,' FILE NOT FOUND!' char(13) '   ',STYLE.SCHAR.downrightArrow,' ',folderParts{end},'  in ',[prefix,fullfile(folderParts{1:end-1})] char(13)];
            break
        end
        if ~exist([prefix,fullfile(folderParts{1:idir})],'dir')
            currentFolder   = [prefix,fullfile(folderParts{1:idir})];
            %errorMessage    = [newline STYLE.SCHAR.problemWarning,' DIRECTORY NOT FOUND!' newline '   ',STYLE.SCHAR.downrightArrow,' ',currentFolder newline]; %newline is only available since 2016b
            errorMessage    = [char(13) STYLE.SCHAR.problemWarning,' DIRECTORY NOT FOUND!' char(13) '   ',STYLE.SCHAR.downrightArrow,' ',currentFolder char(13)];
            break
        end
    end
    FILE.NotFound           = true;
    FILE.NotFoundError      = errorMessage;
    return
    
end

%% DISPLAY INFO
disp(['Name:              ',FILE.name])
if strcmp(SAVE.app,'SL_FieldPlot')
    disp(['Number:            ',num2str(FILE.number)])
    if SWITCH.Verbose && exist('firstNumber','var')
        numberAvailableRangeStr     = [num2str(firstNumber),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(lastNumber)];
        disp(['Numbers available: ',numberAvailableRangeStr])
    end
end
if SWITCH.Verbose
    disp(['Folder:            ',num2str(FILE.stemRead)])
end

%% ADDITIONAL OUTPUT
if strcmp(SAVE.app,'SL_FieldPlot')
    FILE.found              = [FILE.stemRead,FILE.name,fileNameExtension(1,1:end-1),numberingFormat];
    FILE.found(end-length(fnumberString)+1:end)	= fnumberString;
    FILE.foundFieldName    	= fieldName;
    FILE.foundFieldNature   = fieldNature;
elseif strcmp(SAVE.app,'SL_RadialProfile') || strcmp(SAVE.app,'SL_TimeGraph')
    FILE.found              = filetofind;
end
FILE.directory              = FILE.stemRead;
if ispc %just to be sure
    FILE.found              = strrep(FILE.found,'/','\');
    FILE.directory          = strrep(FILE.directory,'/','\');
end
