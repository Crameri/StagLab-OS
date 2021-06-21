
%%                                                        READ STAG-YY 4.62
%    . calls f_readOther
%    . contains f_readTracers
%
%% ABOUT
% Reads a STAGYY output file and transforms it into a MATLAB file
% Syntax:
%     For scalar fields:
%       [X_3D, Y_3D, Z_3D, DATA_3D] = f_readStagYY(directory, fname_input, fname_number, {'viscosity','temperature'}, SWITCH, ioMode);
%     For vector fields:
%       [X_3D, Y_3D, Z_3D, VX_3D, VY_3D, VZ_3D, P_3D] = f_readStagYY(directory, fname_input, fname_number, {'viscosity','temperature'}, SWITCH, ioMode);
%     or simply:
%       [READ] = ReadStag3D(...);
%
% To check if file is found, use:
%       if ischar(X_3D); error(X_3D); end
%
%                                   initially written by Boris Kaus, <2010
%                                    modified by Paul Tackley, August 2010
%                               last modified by Fabio Crameri, 31.05.2021

function [varargout] = f_readStagYY(directory,fnameInput,fnameNumber,Type,SWITCH,ioMode)

%% DEFAULTS
if ~exist('ioMode','var') || isempty(ioMode); ioMode = 1; end

%% CHECK WHICH GEODYNAMIC CODE TO READ
if ~isfield(SWITCH,'GeodynamicCode') || (strcmpi(SWITCH.GeodynamicCode,'StagYY') && ~strcmpi(SWITCH.StagYYoutput,'HDF5')) %StagYY binary output
    % use the function below
  
elseif strcmpi(SWITCH.GeodynamicCode,'StagYY') && strcmpi(SWITCH.StagYYoutput,'HDF5') %StagYY hdf5 output
    [varargout] = f_readStagYYhdf5(directory,fnameInput,fnameNumber,Type,SWITCH,ioMode);
    return
    
elseif strcmpi(SWITCH.GeodynamicCode,'Fluidity')
    [varargout] = f_readFluidity(directory,fnameInput,fnameNumber,Type,SWITCH,ioMode);
    return
    
elseif strcmpi(SWITCH.GeodynamicCode,'Aspect') && strcmpi(SWITCH.ASPECToutput,'HDF5') %Aspect hdf5 output
    [varargout] = f_readASPECThdf5(directory,fnameInput,fnameNumber,Type,SWITCH,ioMode);
    return
    
else %if strcmpi(SWITCH.GeodynamicCode,'->CODENAME<-')
    [varargout] = f_readOther(directory,fnameInput,fnameNumber,Type,SWITCH,ioMode);
    return
end

%% READ FILE
start_dir           = pwd;
cd(directory);

FileFormat          = 'n';           % native - default
% FileFormat          = 'l';           % Little Endian
% FileFormat          = 'b';           % Big    Endian

if fnameNumber<10000
    numberString    = num2str(10000+fnameNumber);
    numberString(1) ='0';
else
    numberString    = num2str(fnameNumber);
end

% Field name,  end-string for binary output,  flag for scalar data,  string for HDF5 output
AvailableFields     = {
    'Tracers',                  '_tra',     false,      'Tracers'
    'Tracer Info',              '_tra',     false,      'xxx'
    'Velocity',                 '_vp',      false,      'Velocity'
    'Horizontal velocity',      '_vp',      false,      'Velocity'
    'Radial velocity',          '_vp',      false,      'Velocity'
    'Pressure',                 '_vp',      false,      'Pressure'
    'Streamfunction',           '_vp',      false,      'Velocity'
    'Temperature',              '_t',       true,       'Temperature'
    'Viscosity',                '_eta',     true,       'Viscosity'
    'Composition',              '_c',       true,       'Composition'
    'Melt Fraction',            '_f',       true,       'MeltFrac'
    'Stress',                   '_str',     true,       'Stress'
    'Horizontal princ. stress', '_sx',      false,      'xxx'
    'Radial princ. stress',     '_sx',      false,      'xxx'
    'Strain rate',              '_ed',      true,       'StrainRate'
    'Density',              	'_rho',     true,       'Density'
    'Phase',                   	'_ph',      true,       'xxx'
    'Air',                      '_air',     true,       'xxx'
    'Primordial',               '_prm',     true,       'xxx'
    'Basalt',                   '_bs',      true,       'Basalt'
    'Harzburgite',              '_hz',      true,       'Harzburgite'
    'Cont. crust',              '_cc',      true,       'xxx'
    'Metal',                    '_mtl',     true,       'xxx'
    'Water',                    '_wtr',     true,       'water'
    'Melt',                     '_mlt',     true,       'xxx'
    'Dyn pressure',             '_pd',      true,       'DynamicPressure'
    'Topography (self-grav)',  	'_csg',     true,       'xxx'
    'Crustal thickness',        '_cr',      true,       'CrustThickness'
    'Heat flux',                '_hf',      true,       'xxx'
    'Surface age',              '_sage',    true,       'xxx'
    'Age since melted',         '_age',     true,       'Age'
    'Geoid',                    '_g',       true,       'xxx'
    'Toroidal',                 '_to',      false,      'xxx'
    'Poloidal',                 '_po',      false,      'xxx'
    'Deformation mechanism',    '_defm',    true,       'xxx'
    'zz-Stress component',      '_nstr',    true,       'xxx'
    'Topo normal stress',       '_nst',     true,       'xxx'
    'Topography',               '_cs',   	true,       'xxx'
};
if max(strcmp(AvailableFields(:,1),Type))==0 %&& max(strcmp(AvailableFields(:,1),Type))==0 %NEEDS TO CHECK FOR FIELD LABEL INSTEAD!!! (IN CASE OF RESIDUAL TEMPERATURE ETC)
    error(['Unknown field: ',Type])
end

%% VARIABLE SETUP
idxType             = find(strcmp(AvailableFields(:,1),Type));
fname               = [fnameInput,AvailableFields{idxType,2},numberString];
scalardata          = AvailableFields{idxType,3};
%check for alternative file endings
if strcmp(Type,'Topography')
    if ~exist(fname,'file')
        fname       = [fnameInput,'_sc',numberString];
        if exist(fname,'file'); disp('     -> _sc used for topography'); end
    end
    if ~exist(fname,'file')
        fname       = [fnameInput,'_fstopo',numberString];
        if exist(fname,'file'); disp('     -> _fstopo used for topography'); end
    end
    if ~exist(fname,'file')
        fname       = [fnameInput,'_cs',numberString]; %just to have the original name for error message
    end
end

%% ERROR CHECK
if ~exist(fname,'file')
    file_stem_now       = pwd;
    cd(start_dir);
    for i=1:20
        varargout{i}    = ['The file - ',file_stem_now,filesep,fname,' - does not exist!'];
    end
    return
    
end

%% SET DATA FORMAT
% Precision
if ~isfield(SWITCH,'Precision') || strcmp(SWITCH.Precision,'single')
    PrecisionString     = {'single','int32'};
elseif strcmp(SWITCH.Precision,'double')
    PrecisionString     = {'double','int64'};
else
    error('Precision unknown!')
end
% Scalar or vector data
if scalardata
    nval                = 1;      % temperature has only one value
else
    nval                = 4;      % assumed that we have a velocity-pressure file
end

%% IMPORT TRACERS
if strcmp(Type,'Tracers') || strcmp(Type,'Tracer Info')
    FILE.Folder         = directory;
    FILE.Name           = fnameInput;
    FILE.Number         = fnameNumber;
    [TRACER] = f_readTracers(Type,FILE,SWITCH);
    
    if isempty(TRACER.NumberTracers)  %if tracer file is empty
        file_stem_now       = pwd;
        cd(start_dir);
        for i=1:20
            varargout{i}    = ['The file - ',file_stem_now,filesep,fname,' - is empty!'];
        end
        return
        
    end
    % prepare output data
    varargout{1}      	= TRACER.NumberTracers;
    varargout{2}        = TRACER.NumberVariables;
    varargout{3}      	= TRACER.VariableNames;
    varargout{4}        = TRACER.data;
    varargout{5}        = TRACER.IdealMass;
    varargout{6}        = TRACER.NumberTraceElements;
    varargout{7}        = TRACER.OutgassedAmount;
    varargout{8}        = TRACER.Time;
    varargout{9}        = TRACER.Rcmb;
    varargout{10}       = TRACER.TimeStep;
    varargout{11}       = TRACER.AspectRatio;
    varargout{12}       = TRACER.nb;

    cd(start_dir);
    return
    
end

%% OPEN FILE
fid             = fopen(fname,'r',FileFormat);              % Open File

%==========================================================================
%% READ HEADER
%==========================================================================
magic       = fread(fid,1,char(PrecisionString(1,2)));  	% Version
if isempty(magic)
    %warning off backtrace; warning(['File ',fname,' seems to be empty! - File read aborted.']); warning on backtrace;
    cd(start_dir);
    for i=1:20
        varargout{i}    = ['File ',fname,' seems to be empty! - File read aborted.'];
    end
    return
    
end
if magic>8000 %64-bit
    PrecisionString     = {'double','int64'};
    fclose(fid);
    fid                 = fopen(fname,'r',FileFormat);              % Re-Open File
    magic           	= fread(fid,1,char(PrecisionString(1,2)));
    magic           	= magic-8000;
end
if (magic<100 && nval>1) || (magic>300 && nval==1)          % check #components
     error('wrong number of components in field')
end
magic       = mod(magic,100);
if magic>=9 && nval==4
    xyp     = 1;     % extra ghost point in x & y direction
else
    xyp     = 0;
end

nxtot       = fread(fid,1,char(PrecisionString(1,2)));     	% nx total
nytot       = fread(fid,1,char(PrecisionString(1,2)));  	% ny total
nztot       = fread(fid,1,char(PrecisionString(1,2)));      % nz total
nblocks     = fread(fid,1,char(PrecisionString(1,2)));    	% # of blocks
Aspect      = fread(fid,2,char(PrecisionString(1,1)));     	% Aspect ratio
nnx         = fread(fid,1,char(PrecisionString(1,2)));  	% Number of parallel subdomains
nny         = fread(fid,1,char(PrecisionString(1,2)));  	%  in the x,y,z and b directions
nnz         = fread(fid,1,char(PrecisionString(1,2)));      %
nnb         = fread(fid,1,char(PrecisionString(1,2)));  	%

nz2         = nztot*2 + 1;
zg          = fread(fid,nz2,char(PrecisionString(1,1)));  	% z-coordinates

% compute nx, ny, nz and nb PER CPU
nx          =   nxtot/nnx;
ny          =   nytot/nny;
nz          =   nztot/nnz;
nb          =   nblocks/nnb;
npi         =   (nx+xyp)*(ny+xyp)*nz*nb*nval;      % the number of values per 'read' block

rcmb      	= fread(fid,1,char(PrecisionString(1,1)));
istep      	= fread(fid,1,char(PrecisionString(1,2)));
time       	= fread(fid,1,char(PrecisionString(1,1)));
erupta_total= fread(fid,1,char(PrecisionString(1,1)));
botT_val 	= fread(fid,1,char(PrecisionString(1,1)));

x           = fread(fid,nxtot,char(PrecisionString(1,1)));      % x-coordinates
y           = fread(fid,nytot,char(PrecisionString(1,1)));      % y-coordinates
z           = fread(fid,nztot,char(PrecisionString(1,1)));      % z-coordinates

% read the parallel blocks
if scalardata
    DATA_3D = zeros(nxtot,nytot,nztot);
else
    scalefac= fread(fid,1,char(PrecisionString(1,1)));     	% scale factor
    VX_3D   = zeros(nxtot,nytot,nztot);                     %   Vx
    VY_3D   = zeros(nxtot,nytot,nztot);                     %   Vy
    VZ_3D   = zeros(nxtot,nytot,nztot);                     %   Vz
    P_3D    = zeros(nxtot,nytot,nztot);                     %   Pressure
end

for ibc=1:nnb       % loop over parallel subdomains
    for izc=1:nnz
        for iyc=1:nny
            for ixc=1:nnx
                data_CPU            = fread(fid,npi,char(PrecisionString(1,1)));      % read the data for this CPU
                
                % Create a 3D matrix from these data
                if scalardata
                    data_CPU_3D     = reshape(data_CPU, [nx ny nz nb]);
                else
                    data_CPU_3D     = reshape(data_CPU*scalefac, [nval nx+xyp ny+xyp nz nb]);
                end
                
                % Add local 3D matrix to global matrix
                if scalardata
                    % Scalar data
                    DATA_3D((ixc-1)*nx+(1:nx), (iyc-1)*ny+(1:ny), (izc-1)*nz+(1:nz),(ibc-1)*nb+(1:nb))  = data_CPU_3D;
                else
                    % velocity-pressure data
                    VX_3D((ixc-1)*nx+(1:nx), (iyc-1)*ny+(1:ny), (izc-1)*nz+(1:nz), (ibc-1)*nb+(1:nb))   = squeeze(data_CPU_3D(1,1:nx,1:ny,:,:));
                    VY_3D((ixc-1)*nx+(1:nx), (iyc-1)*ny+(1:ny), (izc-1)*nz+(1:nz), (ibc-1)*nb+(1:nb))   = squeeze(data_CPU_3D(2,1:nx,1:ny,:,:));
                    VZ_3D((ixc-1)*nx+(1:nx), (iyc-1)*ny+(1:ny), (izc-1)*nz+(1:nz), (ibc-1)*nb+(1:nb))   = squeeze(data_CPU_3D(3,1:nx,1:ny,:,:));
                    P_3D( (ixc-1)*nx+(1:nx), (iyc-1)*ny+(1:ny), (izc-1)*nz+(1:nz), (ibc-1)*nb+(1:nb))   = squeeze(data_CPU_3D(4,1:nx,1:ny,:,:));
                end
            end
        end
    end
end
fclose(fid);                                % close file

[Y_3D, X_3D, Z_3D]      = meshgrid(y,x,z);

%======================================
%% DECREASE DATA SIZE
decreaseFilesize   = false;
nrSaveX             = 2;   %saves every 'nr_saveX' number in x-direction
nrSaveY             = 2;   %                              in y-direction
nrSaveZ             = 2;   %                              in z-direction

if decreaseFilesize
    if nblocks==1
        % no ying-yang
        numx    = size(X_3D,1); numy = size(X_3D,2); numz = size(X_3D,3);
        dummy   = X_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz); X_3D = dummy;
        dummy   = Y_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz); Y_3D = dummy;
        dummy   = Z_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz); Z_3D = dummy;
        if ~scalardata
            dummy   = VX_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz); VX_3D = dummy;
            dummy   = VY_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz); VY_3D = dummy;
            dummy   = VZ_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz); VZ_3D = dummy;
            dummy   = P_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz); P_3D = dummy;
        else
            dummy   = DATA_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz); DATA_3D = dummy;
        end
    else
        % ying-yang grid
        numx    = size(X_3D,1); numy = size(X_3D,2); numz = size(X_3D,3);
        dummy   = X_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz,:); X_3D = dummy;
        dummy   = Y_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz,:); Y_3D = dummy;
        dummy   = Z_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz,:); Z_3D = dummy;
        if ~scalardata
            dummy   = VX_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz,:); VX_3D = dummy;
            dummy   = VY_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz,:); VY_3D = dummy;
            dummy   = VZ_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz,:); VZ_3D = dummy;
            dummy   = P_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz,:); P_3D = dummy;
        else
            dummy   = DATA_3D(1:nrSaveX:numx,1:nrSaveY:numy,1:nrSaveZ:numz,:); DATA_3D = dummy;
        end
    end
    disp(['...saved every ',num2str(nrSaveX),'. x-grid point, ',num2str(nrSaveY),'. y-grid point, ',num2str(nrSaveZ),'. z-grid point'])
end
%======================================

%======================================
%% MIRROR DATA HORIZONTALLY
% SWITCH.flipDataHorizontally = true;
if isfield(SWITCH,'flipDataHorizontally') && SWITCH.flipDataHorizontally
    if nxtot==1 || nytot==1 %2-D geometry
        if logical(0) %flip coordinates
            %         X_3D        = max(X_3D(:))+min(X_3D(:))-X_3D;
            %         Y_3D        = max(X_3D(:))+min(X_3D(:))-X_3D;
            X_3D        = flipud(X_3D);
            Y_3D        = flipud(Y_3D);
        else %flip data
            if scalardata
                DATA_3D = flipud(DATA_3D);
            else
                VX_3D 	= -flipud(VX_3D);
                VY_3D  	= flipud(VY_3D);
                VZ_3D   = flipud(VZ_3D);
                P_3D    = flipud(P_3D);
            end
        end
    else %3-D
        if nb==1 %Cartesian 3-D
            if scalardata
                DATA_3D = flipud(DATA_3D);
            else
                VX_3D 	= -flipud(VX_3D);
                VY_3D  	= flipud(VY_3D);
                VZ_3D   = flipud(VZ_3D);
                P_3D    = flipud(P_3D);
            end
        end
    end
end
%======================================

%======================================
%% SHIFT DATA HORIZONTALLY
if isfield(SWITCH,'shiftDataX') && SWITCH.shiftDataX>0 && SWITCH.shiftDataX<1
    nx2shift    = floor(size(X_3D,1) *SWITCH.shiftDataX);  %shift to the left by half the domain width
    if scalardata
        DATA_3D = [DATA_3D;DATA_3D(1:nx2shift,:,:)];
        DATA_3D(1:nx2shift,:,:) = [];
    else
        VX_3D = [VX_3D;VX_3D(1:nx2shift,:,:)];
        VX_3D(1:nx2shift,:,:) = [];
        VY_3D = [VY_3D;VY_3D(1:nx2shift,:,:)];
        VY_3D(1:nx2shift,:,:) = [];
        VZ_3D = [VZ_3D;VZ_3D(1:nx2shift,:,:)];
        VZ_3D(1:nx2shift,:,:) = [];
        P_3D = [P_3D;P_3D(1:nx2shift,:,:)];
        P_3D(1:nx2shift,:,:) = [];
    end
end
if isfield(SWITCH,'shiftDataY') && SWITCH.shiftDataY>0 && SWITCH.shiftDataY<1
    ny2shift    = floor(size(X_3D,2) *SWITCH.shiftDataY);  %shift to the left by half the domain width
    if scalardata
        DATA_3D = [DATA_3D,DATA_3D(:,1:ny2shift,:)];
        DATA_3D(:,1:ny2shift,:,:) = [];
    else
        VX_3D = [VX_3D,VX_3D(:,1:ny2shift,:)];
        VX_3D(:,1:ny2shift,:,:) = [];
        VY_3D = [VY_3D,VY_3D(:,1:ny2shift,:,:)];
        VY_3D(:,1:ny2shift,:,:) = [];
        VZ_3D = [VZ_3D,VZ_3D(:,1:ny2shift,:)];
        VZ_3D(:,1:ny2shift,:,:) = [];
        P_3D = [P_3D,P_3D(:,1:ny2shift,:)];
        P_3D(:,1:ny2shift,:,:) = [];
    end
end
%======================================

%======================================
%% EXTRACT 2-D SLICE OUT OF 3-D FIELD
if nxtot~=1 && nytot~=1 && nblocks==1 && SWITCH.ThreeToTwoD     %3-D Cartesian geometry 
    %find y-gridpoint corresponding to slice
    sliceYidx               = round(size(Y_3D,2)/2); %middle of domain
    %take slice of data
    X_3D                    = X_3D(:,sliceYidx,:);
    Y_3D                    = Y_3D(:,sliceYidx,:);
    Z_3D                    = Z_3D(:,sliceYidx,:);
    if ~scalardata
        VX_3D               = VX_3D(:,sliceYidx,:);
        VY_3D               = VY_3D(:,sliceYidx,:);
        VZ_3D               = VZ_3D(:,sliceYidx,:);
        P_3D                = P_3D(:,sliceYidx,:);
    else
        DATA_3D             = DATA_3D(:,sliceYidx,:);
    end
    %adjust other variables
end
%======================================

%% OUTPUT
if ioMode==1
    % prepare output data
    varargout{1}            = nblocks;
    varargout{2}            = X_3D;
    varargout{3}            = Y_3D;
    varargout{4}            = Z_3D;
    if nblocks==1
        % no yin-yang
        if ~scalardata
            varargout{5}    = VX_3D;
            varargout{6}    = VY_3D;
            varargout{7}    = VZ_3D;
            varargout{8}    = P_3D;
            varargout{9}    = time;
            varargout{10}   = rcmb;
            varargout{11}   = fname;
            varargout{12}   = Aspect;
        else
            varargout{5}    = DATA_3D;
            varargout{6}    = time;
            varargout{7}    = rcmb;
            varargout{8}    = fname;
            varargout{9}    = Aspect;
        end
    else
        % yin-yang grid
        if ~scalardata
            varargout{5}    = VX_3D(:,:,:,1);
            varargout{6}    = VY_3D(:,:,:,1);
            varargout{7}    = VZ_3D(:,:,:,1);
            varargout{8}    = P_3D (:,:,:,1);
            
            varargout{9}    = VX_3D(:,:,:,2);
            varargout{10}   = VY_3D(:,:,:,2);
            varargout{11}   = VZ_3D(:,:,:,2);
            varargout{12}   = P_3D (:,:,:,2);
            
            varargout{13}   = time;
            varargout{14}   = rcmb;
            varargout{15}   = fname;
            varargout{16}   = Aspect;
        else
            varargout{5}    = DATA_3D(:,:,:,1);
            varargout{6}    = DATA_3D(:,:,:,2);
            varargout{7}    = time;
            varargout{8}    = rcmb;
            varargout{9}    = fname;
            varargout{10}   = Aspect;
        end
    end
    
elseif ioMode==2 % * THIS VERSION NEEDS LESS DATA TRANSFER *
    % prepare output data
    READ.nb                 = nblocks;
    READ.X_3D               = X_3D;
    READ.Y_3D               = Y_3D;
    READ.Z_3D               = Z_3D;
    READ.time               = time;
    READ.rcmb               = rcmb;
    READ.fname              = fname;
    READ.Aspect             = Aspect;
    if nblocks==1
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
    varargout{1}            = READ;
end

cd(start_dir);




%%                                                        READ TRACERS 1.0

%                                                Fabio Crameri, 27.06.2017

function [TRACER] = f_readTracers(Type,FILE,SWITCH)

%% INPUT
if ~exist('FILE','var')
    FILE.Name           = 'n4';
    FILE.Number         = 0;
    FILE.Folder         = '/work/StagYY/+op/n4/';
    FILE.Folder         = '/work/StagYYF_170201/+op/n4/';
end

if ~exist('FILE','var') || ~isfield(FILE,'Name');         	FILE.Name           = 'test';       end
if ~isfield(FILE,'Number');                                	FILE.Number         = 0;            end
if ~isfield(FILE,'Folder');                               	FILE.Folder         = '/work/';  	end
if ~exist('SWITCH','var') || ~isfield(SWITCH,'Precision');	SWITCH.Precision  	= 'single';    	end  %'single' or 'double'

FileFormat              = 'n';           % native - default
% FileFormat              = 'l';           % Little Endian
% FileFormat              = 'b';           % Big    Endian
plotTracerPositions     = logical(0);

%% VARIABLES ADJUSTMENTS
% File name
if FILE.Number<10000
    numberString        = num2str(10000+FILE.Number);
    numberString(1)     ='0';
else
    numberString        = num2str(10000+FILE.Number);
end
FileNameFull            = [FILE.Folder,FILE.Name];
FileNameFullTracer     	= [FileNameFull,'_tra',numberString];
% Precision
if ~isfield(SWITCH,'Precision') || strcmp(SWITCH.Precision,'single')
    PrecisionString     = {'single','int32'};
elseif strcmp(SWITCH.Precision,'double')
    PrecisionString     = {'double','int64'};
else
    error('Precision unknown!')
end

%% OPENING FILE
fid                     = fopen(FileNameFullTracer,'r',FileFormat);              % Open File

%% READING TRACER HEADING
nbFull               	= fread(fid,1,char(PrecisionString(1,2)));
if nbFull>8000 %64-bit
    PrecisionString     = {'double','int64'};
    fclose(fid);
    fid             	= fopen(FileNameFullTracer,'r',FileFormat);              % Open File
    nbFull           	= fread(fid,1,char(PrecisionString(1,2)));
    nbFull           	= nbFull-8000;
end
TRACER.nb            	= mod(nbFull,100);
%      nbFull = nbtot+200   ! 1=compatible with old tracer files
%                           ! +100 adds labels for tracer vars, r_cmb and shape information
%                           ! +200 adds trace element outgassing information

TRACER.AspectRatio      = fread(fid,2,char(PrecisionString(1,1)));
TRACER.TimeStep       	= fread(fid,1,char(PrecisionString(1,2)));
TRACER.Time          	= fread(fid,1,char(PrecisionString(1,1)));
TRACER.NumberVariables 	= fread(fid,1,char(PrecisionString(1,2)));
TRACER.NumberTracers 	= fread(fid,1,char(PrecisionString(1,2)));
if TRACER.nb==2
    TRACER.NumberTracersYang 	= fread(fid,1,char(PrecisionString(1,2)));
end
TRACER.IdealMass      	= fread(fid,1,char(PrecisionString(1,1)));
TRACER.Rcmb             = [];
TRACER.VariableNames    = [];
cart0_cyl1_sph2         = [];
if nbFull>=100
    cart0_cyl1_sph2     = fread(fid,1,char(PrecisionString(1,2)));
    if cart0_cyl1_sph2==2 || cart0_cyl1_sph2==3
        TRACER.Rcmb    	= fread(fid,1,char(PrecisionString(1,1)));
    end
    tracerVarName     	= fread(fid,TRACER.NumberVariables*16,'uint8=>char');
    TRACER.VariableNames = cell(1,TRACER.NumberVariables);
    for iTracerName=1:TRACER.NumberVariables
        TRACER.VariableNames(1,iTracerName) = {strcat(tracerVarName((iTracerName-1)*16+1:iTracerName*16))'};
    end
end
if isempty(cart0_cyl1_sph2); TRACER.Geometry = [];
elseif cart0_cyl1_sph2==0; TRACER.Geometry = 'Cartesian';
elseif cart0_cyl1_sph2==1; TRACER.Geometry = 'spherical2D';
elseif cart0_cyl1_sph2==2; TRACER.Geometry = 'spherical';
end
for iVar=1:TRACER.NumberVariables
    if iVar==1
        TRACER.VariableNamesString	= TRACER.VariableNames{iVar};
    else
        TRACER.VariableNamesString	= [TRACER.VariableNamesString,',',TRACER.VariableNames{iVar}];
    end
end

TRACER.NumberTraceElements      = [];
TRACER.OutgassedAmount          = [];
if nbFull>=200
    TRACER.NumberTraceElements	= fread(fid,1,char(PrecisionString(1,2)));
    TRACER.OutgassedAmount     	= fread(fid,TRACER.NumberTraceElements,char(PrecisionString(1,1)));
end

%% DISPLAY FILE INFORMATIONS
if logical(0)
    disp(' ')
    disp(['Time Step:                ',num2str(TRACER.TimeStep)]);
    disp(['Time:                     ',num2str(TRACER.Time)]);
    disp(' ')
    disp(['#Tracers:                 ',num2str(TRACER.NumberTracers)]);
    disp(['#Tracer Variables:        ',num2str(TRACER.NumberVariables),' (',TRACER.VariableNamesString,')']);
    disp(' ')
    disp(['Tracer Ideal Mass:        ',num2str(TRACER.IdealMass)]);
    disp(['Domain Aspect Ratio:      ',num2str(TRACER.AspectRatio(1),3),' x ',num2str(TRACER.AspectRatio(2),3)]);
    if cart0_cyl1_sph2==0;      disp('Geometry:                 Cartesian');
    elseif cart0_cyl1_sph2==1;  disp(['Geometry:                 Cylindrical with r_cmb = ',num2str(TRACER.Rcmb)]);
    elseif cart0_cyl1_sph2==2;  disp(['Geometry:                 Spherical with r_cmb = ',num2str(TRACER.Rcmb)]);
    end
end

%% READING ACTUAL TRACERS 
TRACER.data = [];
if strcmp(Type,'Tracers')
    if isempty(TRACER.NumberTracers)
        warning off backtrace
        warning('Could not read tracer data: Tracer file is empty!')
        warning on backtrace
    else
        TresholdNumberTracers       = 30e3;   %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        if SWITCH.EfficientTracerReading && TRACER.NumberTracers>TresholdNumberTracers
            ReduceNumTracers        = true;
        else
            ReduceNumTracers        = false;
        end
        if ReduceNumTracers
            %find number of tracers to delete
            if ReduceNumTracers
                NumTracersToRead  	= TRACER.NumberTracers;
                numBatch2skip = 0;
                sizeBatch = TRACER.NumberVariables;
                while NumTracersToRead>TresholdNumberTracers
                    numBatch2skip       = numBatch2skip+1;
                    sizeBatch           = numBatch2skip*TRACER.NumberVariables;
                    %calculate revised total number of tracers to read
                    NumTracersToRead    = ceil(TRACER.NumberTracers*TRACER.NumberVariables / sizeBatch);
                end
            end
        else %read all
            NumTracersToRead        = TRACER.NumberTracers;
            sizeBatch               = TRACER.NumberVariables;
        end
        
        if SWITCH.closeOldWaitbars; wbOld = findall(0,'tag','TMWWaitbar'); delete(wbOld); end
        wb = waitbar(0,'Reading Tracers...');
        TRACER.data = zeros(NumTracersToRead,TRACER.NumberVariables);
        for itracer=1:NumTracersToRead
            waitbar(itracer/NumTracersToRead,wb)
            
            dummy                   = fread(fid,sizeBatch,char(PrecisionString(1,1)))';
            TRACER.data(itracer,:)	= dummy(1,1:TRACER.NumberVariables); %keep reading the correct amount of tracer variables
            
        end
        close(wb)
    end
end
fclose(fid);

%% TEST FIGURE
if plotTracerPositions && strcmp(Type,'Tracers')
    figure(22)
    scatter(TRACER.data(:,1),TRACER.data(:,3))
    axis equal
    axis tight
end






