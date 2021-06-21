
%%                                               READ OTHER CODE FILES 1.23
%
%                                                Fabio Crameri, 27.06.2018
%
%% NOTES
% This needs to read an output file and to transform it into a MATLAB file
%
% Syntax:
%     For scalar fields:
%
%       [X_3D, Y_3D, Z_3D, DATA_3D] = f_readOther(fname_input, fname_number, {'viscosity','temperature'}, SWITCH, ioMode);
%
%   For vector fields:
%
%       [X_3D, Y_3D, Z_3D, VX_3D, VY_3D, VZ_3D, P_3D] = f_readOther(fname_input, fname_number, {'viscosity','temperature'}, SWITCH, ioMode);
%
%   or simply:
%
%       [READ] = ReadStag3D(...);
%
% SWITCH needs to be a field containing SWITCH.Precision
% if file couldn't be found DataIn is a string containing an error message.
%  -> hence use:   if ischar(X_3D); error(X_3D); end     to check if file is found.

function [DataIn] = f_readOther(directory, fname_input, fname_number, Type, SWITCH, ioMode)

error('Adjust here to use StagLab with output files from another geodynamic code.')

start_dir = pwd;
cd(directory);

FileFormat = 'n';           % native - default
% FileFormat = 'l';           % Little Endian
% FileFormat = 'b';           % Big    Endian

%% ERROR CHECKS
if strcmp(Type,'Tracers') || strcmp(Type,'Tracer Info')
    cd(start_dir);
    for i=1:20
        DataIn{i}    = ['Reading tracers not implemented for >thisCode<-Mode!'];
    end
    return
    
end

%% CREATE COMPLETE FILENAME FOR SPECIFIC FIELDS TO READ (change here!)
if fname_number<10000
    number_string = num2str(10000+fname_number);
    number_string(1)='0';
else
    number_string = num2str(fname_number);
end

switch Type
    case 'Velocity'
        fname       = [fname_input,'_vp',number_string];
        scalardata  = false;
    case 'Temperature'
        fname       = [fname_input,'_t',number_string];
        scalardata  = true;
    case 'Viscosity'
        fname       = [fname_input,'_eta',number_string];
        scalardata  = true;
    case 'Composition'
        fname       = [fname_input,'_c',number_string];
        scalardata  = true;
    case 'Melt Fraction'
        fname       = [fname_input,'_f',number_string];
        scalardata  = true;
    case 'Stress'
        fname       = [fname_input,'_str',number_string];
        scalardata  = true;
    case 'Strain rate'
        fname       = [fname_input,'_ed',number_string];
        scalardata  = true;
    case 'Density'
        fname       = [fname_input,'_rho',number_string];
        scalardata  = true;
    case 'Phase'
        fname       = [fname_input,'_ph',number_string];
        scalardata  = true;
    case 'Air'
        fname       = [fname_input,'_air',number_string];
        scalardata  = true;
    case 'Primordial'
        fname       = [fname_input,'_prm',number_string];
        scalardata  = true;
    case 'Basalt'
        fname       = [fname_input,'_bs',number_string];
        scalardata  = true;
    case 'Harzburgite'
        fname       = [fname_input,'_hz',number_string];
        scalardata  = true;
    case 'Cont. crust'
        fname       = [fname_input,'_cc',number_string];
        scalardata  = true;
    case 'Metal'
        fname       = [fname_input,'_mtl',number_string];
        scalardata  = true;
    case 'Water'
        fname       = [fname_input,'_wtr',number_string];
        scalardata  = true;
    case 'Dyn pressure'
        fname       = [fname_input,'_pd',number_string];
        scalardata  = true;
    case 'Topography self-grav'
        fname       = [fname_input,'_csg',number_string];
        scalardata  = true;
    case 'Crustal thickness'
        fname       = [fname_input,'_cr',number_string];
        scalardata  = true;
    case 'Heat flux'
        fname       = [fname_input,'_hf',number_string];
        scalardata  = true;
    case 'Age'
        fname       = [fname_input,'_age',number_string];
        scalardata  = true;
    case 'Geoid'
        fname       = [fname_input,'_g',number_string];
        scalardata  = true;
    case 'Toroidal'
        fname       = [fname_input,'_to',number_string];
        scalardata  = false;
    case 'Poloidal'
        fname       = [fname_input,'_po',number_string];
        scalardata  = false;
    case 'Deformation mechanism'
        fname       = [fname_input,'_defm',number_string];
        scalardata  = true;
    case 'zz-Stress component'
        fname       = [fname_input,'_nstr',number_string];
        scalardata  = true;
    case 'Topo normal stress'
        fname       = [fname_input,'_nst',number_string];
        scalardata  = true;
    case 'Topography'
        fname       = [fname_input,'_cs',number_string];
        scalardata  = true;
    otherwise
        error(['Unknown property: ',Type])
end
if ~exist(fname,'file')
    % The file does not exist and we should stop processing data
    file_stem_now = pwd;
    %warning(['The file - ',file_stem_now,filesep,fname,' - does not exist!'])
    cd(start_dir);
    for i=1:10
        DataIn{i}    = ['The file - ',file_stem_now,filesep,fname,' - does not exist!'];
    end
    
    return
end

if strcmp(SWITCH.Precision,'single')
    PrecisionString={'single' 'int32'};
elseif strcmp(SWITCH.Precision,'double')
    PrecisionString={'double' 'int64'};
else
    error('Precision unknown!')
end

if scalardata
    nval    =   1;      % temperature has only one value
else
    nval    =   4;      % assumed that we have a velocity-pressure file
end

%==========================================================================
%% READ DATA (add here!)

fid         =   fopen(fname,'r',FileFormat);        % Open File

% read data
nx          = 1;
ny          = 1;
nz          = 1;
nb          = 1;    %number blocks (for YinYang 2, otherwise 1)
Aspect      = 1;    %aspect ratio

rcmb      	= 1;    %radius cmb
istep      	= 1;    %time step
time       	= 1;    %total time

x           = 1;   	% x-coordinates
y           = 1;   	% y-coordinates
z           = 1;  	% z-coordinates

[Y_3D, X_3D, Z_3D]  = meshgrid(y,x,z);      %create mesh

% read the fields
if scalardata
    DATA_3D = zeros(nx,ny,nz);
else
    VX_3D   = zeros(nx,ny,nz);                     %   Vx
    VY_3D   = zeros(nx,ny,nz);                     %   Vy
    VZ_3D   = zeros(nx,ny,nz);                     %   Vz
    P_3D    = zeros(nx,ny,nz);                     %   Pressure
end

fclose(fid);                                % close file
%==========================================================================



%% DECREASE SIZE OF INPUT DATA (do not change!)
decrease_filesize = false;
nr_saveX = 2;   %saves every 'nr_saveX' number in x-direction
nr_saveY = 2;   %                              in y-direction
nr_saveZ = 2;   %                              in z-direction

if decrease_filesize
    if nblocks==1
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

%% PREPARE OUTPUT DATA (do not change!)
DataIn{1}            = nb;
DataIn{2}            = X_3D;
DataIn{3}            = Y_3D;
DataIn{4}            = Z_3D;
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

cd(start_dir);



