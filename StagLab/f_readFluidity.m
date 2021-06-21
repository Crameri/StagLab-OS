
%%                                                 READ FLUIDITY FILES 1.11
% 
% NOTES
% This needs to read an .csv output file (e.g., exported from Paraview) and to transform it into a MATLAB file
%
% Things to adjust:
% - Calculate second invariant of the stress tensor instead of just xx stress component
% - Adjust Deformation Mechanism annotation in colorbar
% - currently uses interpolated Topography (on field grid)
% 
% SWITCH needs to be a field containing SWITCH.Precision
% if file couldn't be found DataIn is a string containing an error message.
%  -> hence use:   if ischar(X_3D); error(X_3D); end     to check if file is found.
%
%                                 Fanny Garel, Original script, 21.09.2017
%                                   Fabio Crameri, Adjustments, 23.01.2018

function [DataIn] = f_readFluidity(directory,fnameInput,fnameNumber,Type,SWITCH,ioMode)

startDir = pwd;
cd(directory);

%% INPUT: DETAILS FOR SPECIFIC FILES TO READ (change here!)
factorGridRefinement    = 2;
MaterialName            = 'Normal';
% fnameInput              = 'example_csv_file_with_free_surface';
numberString            = fnameNumber;
fname                   = [fnameInput,'.csv'];

%% OTHER INPUT
FileFormat = 'n';           % native - default
% FileFormat = 'l';           % Little Endian
% FileFormat = 'b';           % Big    Endian

%% ERROR CHECKS
if strcmp(Type,'Tracers') || strcmp(Type,'Tracer Info')
    cd(startDir);
    for i=1:20
        DataIn{i}    = ['Reading tracers not implemented for Fluidity-Mode!'];
    end
    return
    
end

%% FIELD TO READ
field = [];
dataType    = 'field';
% time
field       = char(field,['"',MaterialName,'::Time"']);
% coordinates
field    	= char(field,'"Points:0"'); %x-component
field       = char(field,'"Points:1"'); %z-component
field(1,:) 	= [];
switch Type
    case 'Velocity'
        field       = char(field,['"',MaterialName,'::Velocity:0"']); %x-component
        field       = char(field,['"',MaterialName,'::Velocity:1"']); %y-component
        field       = char(field,['"',MaterialName,'::Velocity:2"']); %z-component
        field       = char(field,['"',MaterialName,'::Pressure"']); %z-component
        scalardata  = false;
    case 'Temperature'
        field       = char(field,['"',MaterialName,'::Temperature"']);
        scalardata  = true;
    case 'Viscosity'
        field       = char(field,['"',MaterialName,'::Viscosity:0"']); %x-component
        scalardata  = true;
    case 'Composition'
        field       = char(field,['"',MaterialName,'::MaterialVolumeFraction"']);
        scalardata  = true;
%     case 'Melt Fraction'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
    case 'Stress'
        field       = char(field,['"',MaterialName,'::Stress:0"']);
        scalardata  = true;
    case 'Strain rate'
        field       = char(field,['"',MaterialName,'::SecondInvariant_StrainRate"']);
        scalardata  = true;
    case 'Density'
        field       = char(field,['"',MaterialName,'::Density"']);
        scalardata  = true;
%     case 'Air'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Primordial'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Basalt'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Harzburgite'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Cont. crust'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Metal'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Water'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Dyn pressure'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Topography self-grav'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Crustal thickness'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Heat flux'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Age'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Geoid'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Toroidal'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = false;
%     case 'Poloidal'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = false;
    case 'Deformation mechanism'
        field       = char(field,['"',MaterialName,'::DeformationMechanism"']);
        scalardata  = true;
%     case 'zz-Stress component'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
%     case 'Topo normal stress'
%         field       = char(field,['"',MaterialName,'::xxx"']);
%         scalardata  = true;
    case 'Topography'
        field       = char(field,['"',MaterialName,'::FreeSurface"']);
        scalardata  = true;
        dataType    = 'graph';
    otherwise
        % The file does not exist and we should stop processing data
        file_stem_now = pwd;
        %warning(['The file - ',file_stem_now,filesep,fname,' - does not exist!'])
        cd(startDir);
        for i=1:10
            DataIn{i}    = ['Unknown property: ',Type];
        end
        return
        
end
if ~exist(fname,'file')
    % The file does not exist and we should stop processing data
    file_stem_now = pwd;
    %warning(['The file - ',file_stem_now,filesep,fname,' - does not exist!'])
    cd(startDir);
    for i=1:10
        DataIn{i}    = ['The file - ',file_stem_now,filesep,fname,' - does not exist!'];
    end
    return

end

if strcmp(SWITCH.Precision,'single')
    PrecisionString = {'single' 'int32'};
elseif strcmp(SWITCH.Precision,'double')
    PrecisionString = {'double' 'int64'};
else
    error('Precision unknown!')
end

if scalardata
    nval    =   1;      % temperature has only one value
else
    nval    =   4;      % assumed that we have a velocity-pressure file
end

%==========================================================================
%% READ DATA
fid         = fopen(fname,'r',FileFormat);        % Open File
header0    	= textscan(fid, [repmat('%s', 1, 100) '%*[^\n]'], 1,'delimiter', ',');
fclose(fid);

header      = string(header0);

%specify indices for different fields of interest
dummy      	= strfind(header,['"',MaterialName,'::Time"']);
idxTime    	= find(not(cellfun('isempty',dummy))); %time
% Index = find(contains(header,['"',MaterialName,'::Time"')); %only for Matlab 2016b and later
dummy      	= strfind(header,'"Points:0"');
idxXcoord  	= find(not(cellfun('isempty',dummy))); %x
dummy      	= strfind(header,'"Points:1"');
idxZcoord  	= find(not(cellfun('isempty',dummy))); %y

if scalardata
    dummy           = strfind(header,strtrim(field(end,:)));
    idxField        = find(not(cellfun('isempty',dummy))); %field
else
    dummy           = strfind(header,strtrim(field(end-3,:)));
    idxField(1,1)  	= find(not(cellfun('isempty',dummy))); %field x
    dummy           = strfind(header,strtrim(field(end-2,:)));
    idxField(2,1)  	= find(not(cellfun('isempty',dummy))); %field y
    dummy           = strfind(header,strtrim(field(end-1,:)));
    idxField(3,1)  	= find(not(cellfun('isempty',dummy))); %field z
    dummy           = strfind(header,strtrim(field(end,:)));
    idxField(4,1)  	= find(not(cellfun('isempty',dummy))); %field pressure
end

%setup data-format vector
Fields_names = [];
format_reading = [];
% for istr=1:length(header)
%     if sum(strcmp(header{istr},field))==1
%         clearvars temp
%         temp            = field(find(strcmp(header{istr},field)==1),:);
%         Fields_names    = char(Fields_names,temp(10:find(temp==['"',1,'last')-1));
%         format_reading  = [format_reading '%f '];
%     else
%         format_reading  = [format_reading '%*f '];
%     end
% end
for istr=1:length(header)
    format_reading  = [format_reading '%f '];
end

% read data
fid                 = fopen(fname,'r');
C                   = textscan(fid,format_reading,'delimiter',',','HeaderLines',1);
fclose(fid);  
% assign data
istep               = NaN;              %time step
time                = C{idxTime}(1); 	%total time
x0                  = C{idxXcoord};     % x-coordinates
y0                  = 1;                % y-coordinates
z0                  = C{idxZcoord};    	% z-coordinates

depthDomain        = max(z0)-min(z0);
widthDomain        = max(x0)-min(x0);
xBegin             = 0;
xEnd               = widthDomain;
%calculate automatic grid size
areaDomain          = depthDomain*widthDomain;
numGridPoints       = length(x0);
averageCellArea   	= areaDomain/numGridPoints;
averageCellWidth    = sqrt(averageCellArea);
nx                  = widthDomain/averageCellWidth;  % resolution for remeshing - change it for fancier figures
ny                  = 1;
nz                  = depthDomain/averageCellWidth;  % resolution for remeshing - change it for fancier figures
nx                  = round(nx*factorGridRefinement);
nz                  = round(nz*factorGridRefinement);

nb                  = 1;    %number blocks (for YinYang 2, otherwise 1)
Aspect              = widthDomain/depthDomain;    %aspect ratio
rcmb                = 1;    %radius cmb

xVec                = linspace(xBegin,xEnd,nx);
yVec                = 1;
zVec                = linspace(0,depthDomain,nz);
[z,x]               = meshgrid(zVec,xVec);

% process data field
[Xtemp,order,~]     = unique(x0(z0==depthDomain));
clearvars temp
warning off
if scalardata
    dummy             	= C{idxField}(z0==depthDomain);
    DATA_Surface        = dummy(order);
    DATA_Surface2       = interp1(Xtemp,DATA_Surface,xVec);  %interpolate onto field grid points
    DATA             	= griddata(z0,x0,C{idxField},z,x);
else
    dummy             	= C{idxField(1,1)}(z0==depthDomain);
    DATA_Surface        = dummy(order);
    DATA_Surface2       = interp1(Xtemp,DATA_Surface,xVec);  %interpolate onto field grid points
    VX                  = griddata(z0,x0,C{idxField(1,1)},z,x);
    VZ                  = griddata(z0,x0,C{idxField(2,1)},z,x); %% CURRENTLY ONLY IMPLEMENTED FOR 2-D
    VY                  = griddata(z0,x0,C{idxField(3,1)},z,x);
    P                   = griddata(z0,x0,C{idxField(4,1)},z,x);
end
warning on

if strcmp(dataType,'graph')
    zVec            = 0;
elseif strcmp(dataType,'field')
      
else
    error('unknown data type!')
end  
    
[Y_3D, X_3D, Z_3D]  = meshgrid(yVec,xVec,zVec);   	%create standard field mesh
if scalardata
    DATA_3D           	= zeros(size(X_3D));
    if strcmp(dataType,'graph')
        DATA_3D(:,1,1) 	= DATA_Surface2;
        
    elseif strcmp(dataType,'field')
        DATA_3D(:,1,:)	= DATA;
        
    else
        error('unknown data type!')
    end
else
    VX_3D           	= zeros(size(X_3D));
    VX_3D(:,1,:)        = VX;
    VY_3D           	= zeros(size(X_3D));
    VY_3D(:,1,:)        = VY;
    VZ_3D           	= zeros(size(X_3D));
    VZ_3D(:,1,:)        = VZ;
    P_3D                = zeros(size(X_3D));
    P_3D(:,1,:)         = P;
end
%==========================================================================


%% CHECKS
%create current coordinates (that might be different from the standard field coordinates)
if ~exist('X_3Dc','var'); X_3Dc = NaN; end
if ~exist('Y_3Dc','var'); Y_3Dc = NaN; end
if ~exist('Z_3Dc','var'); Z_3Dc = NaN; end


%% DECREASE SIZE OF INPUT DATA (do not change!)
decrease_filesize = false;
nr_saveX = 2;   %saves every 'nr_saveX' number in x-direction
nr_saveY = 2;   %                              in y-direction
nr_saveZ = 2;   %                              in z-direction

if decrease_filesize
    if nb==1
        % no ying-yang
        numx = size(X_3D,1); numy = size(X_3D,2); numz = size(X_3D,3);
        dummy = X_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz); X_3D = dummy;
        dummy = Y_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz); Y_3D = dummy;
        dummy = Z_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz); Z_3D = dummy;
        switch Type
            case 'Velocity'
                dummy = VX_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz); VX_3D = dummy;
                dummy = VY_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz); VY_3D = dummy;
                dummy = VZ_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz); VZ_3D = dummy;
                dummy = P_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz); P_3D = dummy;
            otherwise
                dummy = DATA_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz); DATA_3D = dummy;
        end
    else
        % ying-yang grid
        numx = size(X_3D,1); numy = size(X_3D,2); numz = size(X_3D,3);
        dummy = X_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz,:); X_3D = dummy;
        dummy = Y_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz,:); Y_3D = dummy;
        dummy = Z_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz,:); Z_3D = dummy;
        switch Type
            case 'Velocity'
                dummy = VX_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz,:); VX_3D = dummy;
                dummy = VY_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz,:); VY_3D = dummy;
                dummy = VZ_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz,:); VZ_3D = dummy;
                dummy = P_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz,:); P_3D = dummy;
            otherwise
                dummy = DATA_3D(1:nr_saveX:numx,1:nr_saveY:numy,1:nr_saveZ:numz,:); DATA_3D = dummy;
        end
    end
    disp(['...saved every ',num2str(nr_saveX),'. x-grid point, ',num2str(nr_saveY),'. y-grid point, ',num2str(nr_saveZ),'. z-grid point'])
end


%% OUTPUT
if ioMode==1
    % prepare output data
    DataIn{1} = nb;
    DataIn{2} = X_3D;
    DataIn{3} = Y_3D;
    DataIn{4} = Z_3D;
    if nb==1
        % no yin-yang
        if ~scalardata
            DataIn{5}    = VX_3D;
            DataIn{6}    = VY_3D;
            DataIn{7}    = VZ_3D;
            DataIn{8}    = P_3D;
            DataIn{9}    = time;
            DataIn{10}   = rcmb;
            DataIn{11}   = fname;
            DataIn{12}   = Aspect;
        else
            DataIn{5}    = DATA_3D;
            DataIn{6}    = time;
            DataIn{7}    = rcmb;
            DataIn{8}    = fname;
            DataIn{9}    = Aspect;
        end
    else
        % yin-yang grid
        if ~scalardata
            DataIn{5}    = VX_3D(:,:,:,1);
            DataIn{6}    = VY_3D(:,:,:,1);
            DataIn{7}    = VZ_3D(:,:,:,1);
            DataIn{8}    = P_3D (:,:,:,1);
            
            DataIn{9}    = VX_3D(:,:,:,2);
            DataIn{10}   = VY_3D(:,:,:,2);
            DataIn{11}   = VZ_3D(:,:,:,2);
            DataIn{12}   = P_3D (:,:,:,2);
            
            DataIn{13}   = time;
            DataIn{14}   = rcmb;
            DataIn{15}   = fname;
            DataIn{16}   = Aspect;
        else
            DataIn{5}    = DATA_3D(:,:,:,1);
            DataIn{6}    = DATA_3D(:,:,:,2);
            DataIn{7}    = time;
            DataIn{8}    = rcmb;
            DataIn{9}    = fname;
            DataIn{10}   = Aspect;
        end
    end
    
elseif ioMode==2 % * THIS VERSION NEEDS LESS DATA TRANSFER *
    % prepare output data
    READ.nb                 = nb;
    READ.X_3D               = X_3D;
    READ.Y_3D               = Y_3D;
    READ.Z_3D               = Z_3D;
    READ.X_3Dc             	= X_3Dc;
    READ.Y_3Dc            	= Y_3Dc;
    READ.Z_3Dc            	= Z_3Dc;
    READ.time               = time;
    READ.rcmb               = rcmb;
    READ.fname              = fname;
    READ.Aspect             = Aspect;
    if nb==1
        % no yin-yang
        if ~scalardata
            READ.VX_3D      = VX_3D;
            READ.VY_3D      = VY_3D;
            READ.VZ_3D   	= VZ_3D;
            READ.P_3D       = P_3D;
        else
            READ.VAR_3D   	= DATA_3D;
        end
    else
        % yin-yang grid
        if ~scalardata
            READ.VX_3Dyin   = VX_3D(:,:,:,1);
            READ.VY_3Dyin   = VY_3D(:,:,:,1);
            READ.VZ_3Dyin   = VZ_3D(:,:,:,1);
            READ.P_3Dyin    = P_3D(:,:,:,1);
            
            READ.VX_3Dyang	= VX_3D(:,:,:,2);
            READ.VY_3Dyang	= VY_3D(:,:,:,2);
            READ.VZ_3Dyang 	= VZ_3D(:,:,:,2);
            READ.P_3Dyang  	= P_3D(:,:,:,2);
        else
            READ.VAR_3Dyin	= DATA_3D(:,:,:,1);
            READ.VAR_3Dyang	= DATA_3D(:,:,:,2);
        end
    end
    DataIn{1}           	= READ;
end

cd(startDir);



