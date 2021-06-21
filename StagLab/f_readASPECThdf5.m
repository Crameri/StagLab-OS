
%%                                                  READ ASPECT HDF5 0.10
% 
% NOTES
% SWITCH needs to be a field containing SWITCH.Precision
% add to parfile: SWITCH.ASPECToutput = 'HDF5';
% if file couldn't be found DataIn is a string containing an error message.
%  -> hence use:   if ischar(X_3D); error(X_3D); end     to check if file is found.
% 
%                                                Fabio Crameri, 30.05.2021


function [DataIn] = f_readASPECThdf5(directory,fnameInput,fnameNumber,Type,SWITCH,ioMode)

%input to automatise in the future
nb                  = 1;            % number blocks (for YinYang 2, otherwise 1)




%% READ FILE
startDir            = pwd;
cd(directory);

%% ERROR CHECKS
if strcmp(Type,'Tracers') || strcmp(Type,'Tracer Info')
    cd(startDir);
    for i=1:20
        DataIn{i}    = ['Reading tracers not implemented for HDF5-Mode!'];
    end
    return
    
end

% Field name,  end-string for binary output,  flag for scalar data,  string for HDF5 output
AvailableFields     = {
%     'Tracers',                  false,      'Tracers'
%     'Tracer Info',              false,      'xxx'
    'Velocity',                 false,      'velocity'
    'Horizontal velocity',      false,      'velocity'
    'Radial velocity',          false,      'velocity'
    'Pressure',                 true,       'p'
    'Streamfunction',           false,      'velocity'
    'Temperature',              true,       'T'
%     'Viscosity',                true,       'xxx'
%     'Composition',              true,       'xxx'
%     'Melt Fraction',            true,       'xxx'
%     'Stress',                   true,       'xxx'
%     'Horizontal princ. stress', false,      'xxx'
%     'Radial princ. stress',     false,      'xxx'
%     'Strain rate',              true,       'xxx'
%     'Density',              	true,       'xxx'
%     'Air',                      true,       'xxx'
%     'Primordial',               true,       'xxx'
%     'Basalt',                   true,       'xxx'
%     'Harzburgite',              true,       'xxx'
%     'Cont. crust',              true,       'xxx'
%     'Metal',                    true,       'xxx'
%     'Water',                    true,       'xxx'
%     'Melt',                     true,       'xxx'
%     'Dyn pressure',             true,       'xxx'
%     'Topography self-grav',     true,       'xxx'
%    'Crustal thickness',        true,       'CrustThickness' % crustal thickness doesn't work 
%     'Heat flux',                true,       'xxx'
%     'Surface age',              true,       'xxx'
%     'Age since melted',         true,       'xxx'
%     'Geoid',                    true,       'xxx'
%     'Toroidal',                 false,      'xxx'
%     'Poloidal',                 false,      'xxx'
%     'Deformation mechanism',    true,       'xxx'
%     'zz-Stress component',      true,       'xxx'
%     'Topo normal stress',       true,       'xxx'
%     'Topography',               true,       'xxx'
};
if max(strcmp(AvailableFields(:,1),Type))==0
    % The file does not exist and we should stop processing data
    cd(startDir);
    for i=1:10
        DataIn{i}    = ['Unknown property: ',Type];
    end
    return
    
end
fname       = Type;
VARname     = char(AvailableFields(strcmp(AvailableFields(:,1),Type),3));
scalardata  = cell2mat(AvailableFields(strcmp(AvailableFields(:,1),Type),2));

%% CHECK GRID & COORD DETAILS
%% Find setup information
cd('solution/')

% Find grid information
GRIDfilename        = ['mesh-',num2str(0,'%05d'),'.h5'];  %****** adjust if it is not always zero
info                = h5info(GRIDfilename);
dummy               = length(info.Datasets);
for i=1:dummy
    GridParameterNames{i}  	= info.Datasets(i).Name;
end
clearvars dummy
%GridParameterNames
cellSize          	= info.Datasets(1).Dataspace.Size;
nodeSize          	= info.Datasets(2).Dataspace.Size;

numCells            = cellSize(2);
numNodes            = nodeSize(2);


%% Read grid data
CellNumbers       	= h5read(GRIDfilename,'/cells'); 
NodePositions      	= h5read(GRIDfilename,'/nodes');


%% FIND TIME and model setup
% target timestep
iStep               = fnameNumber;    
TIMEfilename        = ['../solution.xdmf'];

%find first time value
fid                 = fopen(TIMEfilename,'r');
find1=0; lineCount=0; timeFound=0; 
secondTimeStringReached=false; lastLineReached=false;
GeometryFound=false;
while ~secondTimeStringReached && ~lastLineReached
    find1 = find1+1; 
    lineCount = lineCount+1;
    xmlstrs{find1,1} = fgetl(fid);

    %Find first time value
    String2Find     = '<Time Value="';
    EndString2Find1 = '"/>';
    strNum1 = strfind(xmlstrs{find1,1},String2Find);
    if ~isempty(strNum1)
        timeFound = timeFound+1;
        strEndNum = strfind(xmlstrs{find1,1},EndString2Find1);
        dummy{1,1}  = xmlstrs{find1,1}(strNum1+length(String2Find):strEndNum-1);
        if timeFound==1
            StartTime   = str2double(dummy);
            lineCount   = 0;
            
        elseif timeFound==2
            SecondTime 	= str2double(dummy);
            linesPerTimeStep = lineCount;   %****** use this to make reading times more efficient by jumping to relevant line (after below information was extracted)
            secondTimeStringReached = true;
        end
    end

    %check for GeometryType
    String2Find     = '<Geometry GeometryType="';
    EndString2Find1 = '">';
    strNum1 = strfind(xmlstrs{find1,1},String2Find);
    if ~GeometryFound && ~isempty(strNum1)
        strEndNum = strfind(xmlstrs{find1,1},EndString2Find1);
        GeometryType  = char( xmlstrs{find1,1}(strNum1+length(String2Find):strEndNum-1) );
        GeometryFound = true;
    end
    
    %check for last line
    String2Find     = '</Xdmf>';
    strNum1 = strfind(xmlstrs{find1,1},String2Find);
    if ~isempty(strNum1)
        lastLineReached = true;
    end
end
fclose(fid);
clearvars dummy


%find current time value
%find first time value
fid                 = fopen(TIMEfilename,'r');
find1=0; lineCount=0; timeStepCount=-1;
currentTimeStringReached=false; lastLineReached=false;
while ~currentTimeStringReached && ~lastLineReached
    find1 = find1+1; 
    lineCount = lineCount+1;
    xmlstrs{find1,1} = fgetl(fid);

    %Find first time value
    String2Find     = '<Time Value="';
    EndString2Find1 = '"/>';
    strNum1 = strfind(xmlstrs{find1,1},String2Find);
    if ~isempty(strNum1)
        timeStepCount   = timeStepCount+1;
        strEndNum = strfind(xmlstrs{find1,1},EndString2Find1);
        dummy{1,1}  = xmlstrs{find1,1}(strNum1+length(String2Find):strEndNum-1);
        if timeStepCount==iStep
            CurrentTime   = str2double(dummy);
            currentTimeStringReached = true;
        end
    end

    %check for last line
    String2Find     = '</Xdmf>';
    strNum1 = strfind(xmlstrs{find1,1},String2Find);
    if ~isempty(strNum1)
        lastLineReached = true;
    end
end
fclose(fid);



time                = CurrentTime;


%% Setup interpolated grid
if strcmpi(GeometryType,'XY')
    xdiff           = max(NodePositions(1,:))-min(NodePositions(1,:));
    zdiff           = max(NodePositions(2,:))-min(NodePositions(2,:));
    AspectRatio     = round(xdiff/zdiff);
    width           = AspectRatio;
    height          = 1;
    nx              = sqrt( (width/height)*numNodes+((width-height)^2)/(4*height^2) )-(width-height)/(2*height);
    ny              = 1;
    nz              = numNodes/nx;
    dx              = xdiff/(nx-1);
    dy              = 1;
    dz              = zdiff/(nz-1);
    
else
    error(['Model geometry ',GeometryType,' not yet supported. Check here!'])
end


%% FIND FILE 
% Read parameter field file
VARfilename         = ['solution-',num2str(fnameNumber,'%05d'),'.h5'];
% Error check
if ~exist(VARfilename,'file')
    file_stem_now 	= pwd;
    cd(startDir);
    for i=1:20
        DataIn{i}  	= ['The file - ',directory,'solution',filesep,VARfilename,' - does not exist!'];
    end
    return
    
end
info                = h5info(VARfilename);
NumberParameters    = length(info.Datasets);
for iParameter=1:NumberParameters
    ParameterNames{iParameter}  	= info.Datasets(iParameter).Name;
end
% ParameterNames



%% Setup and interpolate grid
%X_3D = NaN(nx,ny,nz); Y_3D = NaN(nx,ny,nz); Z_3D = NaN(nx,ny,nz);
DATA_3D = NaN(nx,ny,nz); DATAX_3D = NaN(nx,ny,nz); DATAY_3D = NaN(nx,ny,nz); DATAZ_3D = NaN(nx,ny,nz);

if strcmpi(GeometryType,'XY')
    xGrid           = min(NodePositions(1,:)):dx:max(NodePositions(1,:));
    yGrid           = 0;
    zGrid           = min(NodePositions(2,:)):dz:max(NodePositions(2,:));
    [Y_3D,X_3D,Z_3D]= meshgrid(yGrid,xGrid,zGrid);
    [z_2D,x_2D]   	= meshgrid(zGrid,xGrid); %2-D grid
    
    % Read parameter field data
    DATAfull        = h5read(VARfilename,['/',VARname]);
    if scalardata
        if ~SWITCH.Verbose; warning off; end
        DATA          	= griddata(NodePositions(2,:),NodePositions(1,:),DATAfull,z_2D,x_2D);
        if ~SWITCH.Verbose; warning on; end
        
        % Assign data
        DATA_3D(:,:,:) 	= DATA;
    else %vector data
        
        DATAfullx   = DATAfull(1,:);
        DATAfully   = DATAfull(2,:);
        %DATAfullz   = DATAfull(3,:);
        if ~SWITCH.Verbose; warning off; end
        DATAx          	= griddata(NodePositions(2,:),NodePositions(1,:),DATAfullx,z_2D,x_2D);
        DATAy          	= griddata(NodePositions(2,:),NodePositions(1,:),DATAfully,z_2D,x_2D);
        if ~SWITCH.Verbose; warning on; end
        
        % Assign data
        DATAX_3D(:,:,:) 	= DATAx;
        DATAY_3D(:,:,:) 	= DATAx.*0;
        DATAZ_3D(:,:,:) 	= DATAy;
    end
else
    error(['Model geometry ',GeometryType,' not yet supported. Check here!'])
end








%% CHECK GEOMETRY
if min(Z_3D(:))<0   %MIGHT NEED ADJUSTMENT (probably not optimal)******************************
    GType           = 'spherical2D';
else
    GType           = 'Cartesian';
end

%% DATA ADJUSTMENTS
if strcmp(GType,'spherical2D')
    %Transformation to spherical2D coordinates (x-data needs to be in radians)
    if nx>1
        hDummy      = atan(abs(Z_3D./X_3D));
        zDummy      = sqrt(X_3D.^2+Z_3D.^2);
        hDummy(Z_3D>0 & X_3D<0) = pi-hDummy(Z_3D>0 & X_3D<0);
        hDummy(Z_3D<0 & X_3D<0) = pi+hDummy(Z_3D<0 & X_3D<0);
        hDummy(Z_3D<0 & X_3D>0) = 2*pi-hDummy(Z_3D<0 & X_3D>0);
    else
        hDummy      = atan(abs(Z_3D./Y_3D));
        zDummy      = sqrt(Y_3D.^2+Z_3D.^2);
        hDummy(Z_3D>0 & Y_3D<0) = pi-hDummy(Z_3D>0 & Y_3D<0);
        hDummy(Z_3D<0 & Y_3D<0) = pi+hDummy(Z_3D<0 & Y_3D<0);
        hDummy(Z_3D<0 & Y_3D>0) = 2*pi-hDummy(Z_3D<0 & Y_3D>0);
    end
    rcmb            = min(zDummy(:));    %radius cmb
    zDummy          = zDummy-min(zDummy(:));
    Y_3D            = hDummy;
    Z_3D            = zDummy;
else
    rcmb            = -1; %StagYY's flag for Cartesian geometry
end

if nx>1
    Aspect        	= nx/nz;    %aspect ratio
else
    Aspect        	= ny/nz;    %aspect ratio
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
            DataIn{5}    = DATAX_3D;
            DataIn{6}    = DATAY_3D;
            DataIn{7}    = DATAZ_3D;
            DataIn{8}    = NaN;
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
            DataIn{5}    = DATAX_3D(:,:,:,1);
            DataIn{6}    = DATAY_3D(:,:,:,1);
            DataIn{7}    = DATAZ_3D(:,:,:,1);
            DataIn{8}    = NaN;
            
            DataIn{9}    = DATAX_3D(:,:,:,2);
            DataIn{10}   = DATAY_3D(:,:,:,2);
            DataIn{11}   = DATAZ_3D(:,:,:,2);
            DataIn{12}   = NaN;
            
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
    READ.time               = time;
    READ.rcmb               = rcmb;
    READ.fname              = fname;
    READ.Aspect             = Aspect;
    if nb==1
        % no yin-yang
        if ~scalardata
            READ.VX_3D      = DATAX_3D;
            READ.VY_3D      = DATAY_3D;
            READ.VZ_3D   	= DATAZ_3D;
            READ.P_3D       = NaN;
        else
            READ.VAR_3D   	= DATA_3D;
        end
    else
        % yin-yang grid
        if ~scalardata
            READ.VX_3Dyin   = DATAX_3D(:,:,:,1);
            READ.VY_3Dyin   = DATAY_3D(:,:,:,1);
            READ.VZ_3Dyin   = DATAZ_3D(:,:,:,1);
            READ.P_3Dyin    = NaN;
            
            READ.VX_3Dyang	= DATAX_3D(:,:,:,2);
            READ.VY_3Dyang	= DATAY_3D(:,:,:,2);
            READ.VZ_3Dyang 	= DATAZ_3D(:,:,:,2);
            READ.P_3Dyang  	= NaN;
        else
            READ.VAR_3Dyin	= DATA_3D(:,:,:,1);
            READ.VAR_3Dyang	= DATA_3D(:,:,:,2);
        end
    end
    DataIn{1}           	= READ;
end

cd(startDir);

end



