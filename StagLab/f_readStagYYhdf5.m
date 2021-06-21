
%%                                                  READ STAG-YY HDF5 0.20
%    . contains vec2mat
%    . contains col2row
% 
% NOTES
% SWITCH needs to be a field containing SWITCH.Precision
% add to parfile: SWITCH.StagYYoutput = 'HDF5';
% if file couldn't be found DataIn is a string containing an error message.
%  -> hence use:   if ischar(X_3D); error(X_3D); end     to check if file is found.
% 
%                            Kiran Chotalia, Original script, 22.01.2018
%                            Fabio Crameri, Further adjustments, 24.01.2018
%                            Anna Guelcher, New improvements, 24.09.2020
%                                   - check grid & coord details (still beta)
%                                   - improvements to find correct file:
%                                     - find TIME 
%                                     - FIND NUMFILES and NUMSTEPS
%                                     - FIND FILE
%                                   - load data 
%                                   - reshape data if vertical splitting
%                                   into nodes
%                                   - adjustments for vector data (Fabio)

function [DataIn] = f_readStagYYhdf5(directory,fnameInput,fnameNumber,Type,SWITCH,ioMode)

warning('Processing HDF5 Output from StagYY is still in beta, and not fully functional!')

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
    'Velocity',                 false,      'Velocity'
    'Horizontal velocity',      false,      'Velocity'
    'Radial velocity',          false,      'Velocity'
    'Pressure',                 true,       'Pressure'
    'Streamfunction',           false,      'Velocity'
    'Temperature',              true,       'Temperature'
    'Viscosity',                true,       'Viscosity'
    'Composition',              true,       'Composition'
    'Melt Fraction',            true,       'MeltFrac'
    'Stress',                   true,       'Stress'
%     'Horizontal princ. stress', false,      'xxx'
%     'Radial princ. stress',     false,      'xxx'
    'Strain rate',              true,       'StrainRate'
    'Density',              	true,       'Density'
%     'Air',                      true,       'xxx'
%     'Primordial',               true,       'xxx'
    'Basalt',                   true,       'Basalt'
    'Harzburgite',              true,       'Harzburgite'
%     'Cont. crust',              true,       'xxx'
%     'Metal',                    true,       'xxx'
    'Water',                    true,       'water'
%     'Melt',                     true,       'xxx'
    'Dyn pressure',             true,       'DynamicPressure'
%     'Topography self-grav',     true,       'xxx'
%    'Crustal thickness',        true,       'CrustThickness' % crustal thickness doesn't work 
%     'Heat flux',                true,       'xxx'
%     'Surface age',              true,       'xxx'
    'Age since melted',         true,       'Age'
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

%% CHECK GRID & COORD DETAILS
fileCoordInfo1   	= dir('*_NodeCoordinates.h5');  %check for the coordinates files
fileCoordInfo1([fileCoordInfo1.isdir])	= [];                   %remove directories
fileCoordInfo   	= fileCoordInfo1(end).name;         %only check first one found
numNodes            = str2double(fileCoordInfo(19:23));
%subdomain grid:
nxNode           	= str2double(fileCoordInfo(1:5));
nyNode            	= str2double(fileCoordInfo(7:11));
nzNode            	= str2double(fileCoordInfo(13:17));
%global grid (for now input, has to be given!)
nx                  = 1;
ny                  = 512;
nz                  = 96;
% Converting the subdomain grid numbers to global grid numbers does not work yet, 
% as it depends on the geometrical division of the main domain (# cores, etc.)
% Also, vertical splitting is not yet implemented in the reading of the 3D data.
NumNxNode           = nx/nxNode;
NumNyNode           = ny/nyNode;
NumNzNode           = nz/nzNode;

%check whether sum(numNodes*(subdomain grid size)) fits into global grid size (nx*ny*nz) 
if (numNodes/(NumNxNode * NumNyNode * NumNzNode)~=1)
    error('global grid and subdomain grid dimensions do not match')  
end

%% FIND TIME 
% target timestep
iStep               = fnameNumber;    
TIMEfilename        = ['time_botT.h5'];
TIMEname            = ['time_botT_',num2str(iStep,'%05d')];
time                = h5read(TIMEfilename,['/',TIMEname]);
time                = time(1);

%% FIND NUMFILES and NUMSTEPS (AG)
i0File              = 0;    % start number of file 
i0Node              = 1;    % start number of node
idxType             = find(strcmp(AvailableFields(:,1),Type));
TypeLength          = strlength(AvailableFields{idxType,3});

% Find # h5 files (varies per model, max. file size and set-up) 
fileFieldStringInfo1    = dir([AvailableFields{idxType,3},'_*_',num2str(i0Node,'%05d'),'.h5']);      %check for the coordinates files
fileFieldStringInfo1([fileFieldStringInfo1.isdir])	= [];   % remove directories
fileFieldStringInfo     = fileFieldStringInfo1(end).name;  	% only check last one 
numFiles            	= str2double(fileFieldStringInfo((TypeLength+2):(TypeLength+6)));

% Find # steps written in each file
VARfilename         = [AvailableFields{idxType,3},'_',num2str(i0File,'%05d'),'_',num2str(i0Node,'%05d'),'.h5'];
info                = h5info(VARfilename);
numSteps            = size(info.Datasets,1); %This is how many timesteps are written in this file

%% FIND FILE 
nb                  = 1;            % number blocks (for YinYang 2, otherwise 1)
iNode               = 1;            % start number of node [0,numNodes]
iFile = floor(iStep/numSteps);  % this is the correct file number, still beta for vector data files
fname               = [AvailableFields{idxType,3},'_',num2str(iFile,'%05d'),'_',num2str(iNode,'%05d')];
scalardata          = AvailableFields{idxType,2};

%% ERROR CHECK
if ~exist([fname,'.h5'],'file')
    file_stem_now 	= pwd;
    cd(startDir);
    for i=1:20
        DataIn{i}  	= ['The file - ',file_stem_now,filesep,fname,'.h5 - does not exist!'];
    end
    return
    
end

%% OPEN FILE

X_3D = NaN(nx,ny,nz); Y_3D = NaN(nx,ny,nz); Z_3D = NaN(nx,ny,nz);
DATA_3D = NaN(nx,ny,nz); DATAX_3D = NaN(nx,ny,nz); DATAY_3D = NaN(nx,ny,nz); DATAZ_3D = NaN(nx,ny,nz);
Xfull = []; Yfull = []; Zfull = []; DATAfull = []; DATAfullX = []; DATAfullY = []; DATAfullZ = [];
for iNode=1:numNodes %loop through nodes
    %% Load Grid Data
    GRIDfilename    = [num2str(nxNode,'%05d'),'_',num2str(nyNode,'%05d'),'_',num2str(nzNode,'%05d'),'_',num2str(iNode,'%05d'),'_NodeCoordinates.h5'];
    X1(:,:)       	= h5read(GRIDfilename,'/X');
    Y1(:,:)      	= h5read(GRIDfilename,'/Y');
    Z1(:,:)     	= h5read(GRIDfilename,'/Z');
    
    %% Load Field Data
    VARfilename     = [AvailableFields{idxType,3},'_',num2str(iFile,'%05d'),'_',num2str(iNode,'%05d'),'.h5'];
    VARname         = [AvailableFields{idxType,3},'_',num2str(iNode,'%05d'),'_',num2str(iStep,'%05d')];
    if scalardata
        VAR1(:,:)  	= h5read(VARfilename,['/',VARname]);
    else %vector data
        dummy(:,:,:)= h5read(VARfilename,['/',VARname]);
        if min([nx ny nz])>1
            VAR1X(:,:) 	= dummy(1,:,:); 
            VAR1Y(:,:)	= dummy(2,:,:);
            VAR1Z(:,:) 	= dummy(3,:,:);
        elseif nx==1              
            VAR1X(:,:) 	= dummy(3,:,:); %zeroes, since nx = 1
            VAR1Y(:,:)	= dummy(1,:,:);
            VAR1Z(:,:) 	= dummy(2,:,:);
        elseif ny==1 
            VAR1X(:,:) 	= dummy(1,:,:); 
            VAR1Y(:,:)	= dummy(3,:,:); %zeroes, since ny = 1
            VAR1Z(:,:) 	= dummy(2,:,:);
        elseif nz==1 
        end
    end
    
    %% Flip rows and columns   
    X1              = col2row(X1);
    Y1              = col2row(Y1);
    Z1              = col2row(Z1);
    if scalardata
        VAR1        = col2row(VAR1);
    else %vector data
        VAR1X    	= col2row(VAR1X);
        VAR1Y      	= col2row(VAR1Y);
        VAR1Z    	= col2row(VAR1Z);
    end
    
    %% Trim Grid Data
    X1              = X1(1:end-1,1:end-1);
    Y1              = Y1(1:end-1,1:end-1);
    Z1              = Z1(1:end-1,1:end-1);
    if ~scalardata %velocity is on nodes
        VAR1X      	= VAR1X(1:end-1,1:end-1);
        VAR1Y      	= VAR1Y(1:end-1,1:end-1);
        VAR1Z      	= VAR1Z(1:end-1,1:end-1);
    end
    
    %% Load Data
    Xfull           = [Xfull, X1];
    Yfull           = [Yfull, Y1];
    Zfull           = [Zfull, Z1];
    if scalardata
        DATAfull 	= [DATAfull, VAR1];
    else %vector data
        DATAfullX 	= [DATAfullX, VAR1X];
        DATAfullY  	= [DATAfullY, VAR1Y];
        DATAfullZ 	= [DATAfullZ, VAR1Z];
    end
    clearvars VAR1 VAR1X VAR1Y VAR1Z X1 Y1 Z1
end

%% RESHAPE DATA 
% AG Needed when vertical splitting is done
% only works for 2D y-z geometry for now
if NumNzNode==2
    if (nx==1)
        Xfull_new = zeros(nz,ny);
        Yfull_new = zeros(nz,ny);
        Zfull_new = zeros(nz,ny);
        
        Xfull_new(1:nzNode,:) = Xfull(:,1:ny);
        Yfull_new(1:nzNode,:) = Yfull(:,1:ny);
        Zfull_new(1:nzNode,:) = Zfull(:,1:ny);
        Xfull_new(nzNode+1:end,:) = Xfull(:,ny+1:end);
        Yfull_new(nzNode+1:end,:) = Yfull(:,ny+1:end);
        Zfull_new(nzNode+1:end,:) = Zfull(:,ny+1:end);
        if scalardata
            DATAfull_new = zeros(nz,ny);
            
            DATAfull_new(1:nzNode,:)     = DATAfull(:,1:ny);
            DATAfull_new(nzNode+1:end,:) = DATAfull(:,ny+1:end);

        else
            DATAfullX_new = zeros(nz,ny);
            DATAfullY_new = zeros(nz,ny);
            DATAfullZ_new = zeros(nz,ny);
            
            DATAfullX_new(1:nzNode,:)     = DATAfullX(:,1:ny);
            DATAfullY_new(1:nzNode,:)     = DATAfullY(:,1:ny);
            DATAfullZ_new(1:nzNode,:)     = DATAfullZ(:,1:ny);
            DATAfullX_new(nzNode+1:end,:) = DATAfullX(:,ny+1:end);
            DATAfullY_new(nzNode+1:end,:) = DATAfullY(:,ny+1:end);
            DATAfullZ_new(nzNode+1:end,:) = DATAfullZ(:,ny+1:end);
        end    
    elseif(ny==1)
        % TO DO for x-z geometry
    end
elseif (NumNzNode==1)
    Xfull_new = Xfull;
    Yfull_new = Yfull;
    Zfull_new = Zfull;
    if scalardata
        DATAfull_new = Datafull;
    else
        DATAfullX_new = DatafullX;
        DATAfullY_new = DatafullY;
        DATAfullZ_new = DatafullZ;
    end    
end

%% ASSIGN DATA
X_3D(:,:,:)     	= Xfull_new';
Y_3D(:,:,:)      	= Yfull_new';
Z_3D(:,:,:)       	= Zfull_new';
if scalardata
    DATA_3D(:,:,:) 	= DATAfull_new';
else
    DATAX_3D(:,:,:)	= DATAfullX_new';
    DATAY_3D(:,:,:)	= DATAfullY_new';
    DATAZ_3D(:,:,:)	= DATAfullZ_new';
end

%% CHECK GEOMETRY
if min(Z_3D(:))<0   %MIGHT NEED ADJUSTMENT (probably not optimal)
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


%% NEW ADUSTMENT FOR VECTOR DATA (VX VY ETC)
% Only works for ny-nz field, TO DO: nx-nz. 
if strcmp(GType,'spherical2D')
    %Transformation to match StagYY binary output
    DATAX_3D  	= -DATAX_3D;
    DATAY_3D  	= -DATAY_3D;
    DATAZ_3D   	= -DATAZ_3D;
    
    dummyZ      = (-DATAZ_3D.*sin(Y_3D)-DATAY_3D.*cos(Y_3D)); %convert back to Cartesian form (HDF5 vector data is already in spherical form, for some reason)
    dummyY      = (DATAY_3D.*sin(Y_3D)-DATAZ_3D.*cos(Y_3D));
    DATAZ_3D   	= dummyZ;
    DATAY_3D  	= dummyY;
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



function OUT = col2row(IN)
%% Changes to Grid
%  from
%      1     5     9    13    17
%      2     6    10    14    18
%      3     7    11    15    19
%      4     8    12    16    20
%  to
%      1     2     3     4     5
%      6     7     8     9    10
%     11    12    13    14    15
%     16    17    18    19    20
% 
% for changing HDF5 output from V84_StagYY so it's easier to understand and plot
% 
%                                               Kiran Chotalia, 22.01.2018

% Matrix where numbers read down columns
x = IN;
% size of matrix
sz = size(x);
% Takes input and transforms it into a vector
x = x(:)';
% Vector changed back to matrix so numbers read across rows
OUT = vec2mat(x,sz(2));
end



function [mat, padded] = vec2mat(vec, matCol, padding)
%VEC2MAT Convert a vector into a matrix.
%   MAT = VEC2MAT(VEC, MATCOL) converts the vector VEC into a matrix with
%   MATCOL columns, creating one row at a time. If the length of VEC is not
%   a multiple of MATCOL, then the function places extra entries of 0 in the
%   last row of MAT.
%
%   MAT = VEC2MAT(VEC, MATCOL, PADDING) is the same as the first syntax,
%   except that the extra entries are taken from the matrix PADDING, in order.
%   If the number of elements in PADDING is insufficient, the last element is
%   used for the remaining entries.
%
%   [MAT, PADDED] = VEC2MAT(...) returns an integer PADDED that indicates
%   how many extra entries were placed in the last row of MAT.
%
%   See also RESHAPE.

%   Copyright 1996-2014 The MathWorks, Inc.

narginchk(2,3);	% 2 or 3 inputs required
if ndims(vec) > 2
    error(message('comm:vec2mat:InvalidVec'));
elseif (length(matCol) ~= 1 || ~isfinite(matCol) || ~isreal(matCol)...
        || floor(matCol) ~= matCol || matCol < 1)
    error(message('comm:vec2mat:InvalidMatcol'));
end

[vecRow, vecCol] = size(vec);
vecLen = vecRow*vecCol;
if vecCol == matCol
    mat = vec;
    padded = 0;
    return;			% nothing to do
elseif vecRow > 1
    vec = reshape(vec.', 1, vecLen);
end

try
    if nargin < 3 || isempty(padding)
        padding = cast(0, class(vec));	% default padding
    else
        padding = cast(padding(:).', class(vec));
    end
catch exception
    throw(exception)
end
paddingLen = length(padding);

matRow = ceil(vecLen/matCol);
padded = matRow*matCol - vecLen;	% number of elements to be padded
vec = [vec, padding(1:min(padded, paddingLen)),...
       repmat(padding(paddingLen), 1, padded-paddingLen)];	% padding
mat = reshape(vec, matCol, matRow).';
end
