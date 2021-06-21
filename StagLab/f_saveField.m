
%%                                                         SAVE FIELD  2.12
%
%                                                Fabio Crameri, 21.06.2019

function f_saveField(x2d,y2d,z2d,var2d,FIELD,FILE,SAVE)

if ~isfield(SAVE,'dat');                SAVE.dat                    = logical(0);           end
if ~isfield(SAVE,'mat');                SAVE.mat                    = logical(1);           end
if ~isfield(SAVE,'Directory');          SAVE.Directory      = ['~',filesep,'Desktop',filesep,FILE.name,filesep];	end
if ~isfield(SAVE,'BulletString');     	SAVE.BulletString       	= char(9862);           end
if ~isfield(SAVE,'PointingString');   	SAVE.PointingString       	= char(10551);        	end
% if ~isfield(SAVE,'overwrite');          SAVE.overwrite              = logical(0);           end

%% INPUT
exitAfterSaving     = logical(0);
calculatePHI        = logical(0);

%% DEFAULTS
successfulSaving    = false;
x2d; %these are used later on!
y2d;
z2d;

%% CHECK GRID DIMENSION
if size(y2d,2)==1
    gridDim = '2-D';
elseif size(y2d,2)>1
    gridDim = '3-D';
else
    error('could not figure out grid dimension')
end

%% OPTIONAL FIELD CONVERSIONS
if calculatePHI && strcmp(FIELD.name,'Temperature') %calculate phi from T-field
    T_adi = 0.5;
    T_liq = 2.*T_adi;
    T_sol = T_liq*(1.-log(0.79))^(-1.);
    %     T_adi = 0.5;
    %     T_liq = 2.35*T_adi;
    %     T_sol = T_liq*(1.-log(0.79))^(-1.);
    T_adi = T_adi *2500; %dimensionalize
    T_liq = T_liq *2500; %dimensionalize
    T_sol = T_sol *2500; %dimensionalize
    
    phi2d = (var2d-T_sol)./(T_liq-T_sol);
else
    phi2d = 1; %just to be defined
end

%% ADJUST 3-D DATA to save as .dat
if strcmp(gridDim,'3-D')
    %make array data
    [x3d,y3d,z3d] = meshgrid(x2d,y2d,z2d);
    x_array     = x3d(:);
    y_array     = y3d(:);
    z_array     = z3d(:);
    var2d_array = var2d(:);
    phi_array   = phi2d(:);
end

%% ADJUST SAVE DIRECTORY WITH +field FOLDER
saveDirectory = [SAVE.Directory,'+field/'];

%% CHECK IF FILE DIRECTORY ALREADY EXISTS
if ~exist(saveDirectory,'dir')
    c=questdlg(['Create directory "',saveDirectory,'"?'],'Directory does not exist!', ...
        'Yes','No','Yes');
    switch c
        case 'Yes'
            %create directory
            mkdir(saveDirectory);
        case 'No'
            error('So, change save-directory and re-run again.')
    end
end

%% SAVE DATA
if SAVE.dat
    if strcmp(gridDim,'3-D') %change to array data
        x2d0 = x2d; y2d0 = y2d; z2d0 = z2d; var2d0 = var2d;
        phi2d0 = phi2d;
        x2d = x_array; y2d = y_array; z2d = z_array; var2d = var2d_array;
        phi2d = phi_array;
    end
    save([saveDirectory,filesep,'x2d.dat'],'-ascii', 'x2d');
    if strcmp(gridDim,'3-D'); save([saveDirectory,filesep,'y2d.dat'],'-ascii', 'y2d'); end
    save([saveDirectory,filesep,'z2d.dat'],'-ascii', 'z2d');
    save([saveDirectory,filesep,FIELD.symbol,'_',num2str(FILE.number),'.dat'],'-ascii', 'var2d');
    %save additional data
    if calculatePHI; save([saveDirectory,filesep,'phi_',num2str(FILE.number),'.dat'],'-ascii', 'phi2d'); end
    if strcmp(gridDim,'3-D') %reset to matrix data
        x2d = x2d0; y2d = y2d0; z2d = z2d0; var2d = var2d0;
        phi2d = phi2d0;
    end
    successfulSaving    = true;
end
if SAVE.mat
    save([saveDirectory,filesep,'x2d.mat'], 'x2d');
    if strcmp(gridDim,'3-D'); save([saveDirectory,filesep,'y2d.mat'], 'y2d'); end
    save([saveDirectory,filesep,'z2d.mat'], 'z2d');
    save([saveDirectory,filesep,FIELD.symbol,'_',num2str(FILE.number),'.mat'], 'var2d');
    %save additional data
    if calculatePHI; save([saveDirectory,filesep,'phi_',num2str(FILE.number),'.mat'], 'phi2d'); end
    successfulSaving    = true;
end

%% DISPLAY INFORMATION
if successfulSaving
    disp(' ')
    disp(['   ',SAVE.BulletString,' Field Saved (',FIELD.symbol,')'])
    disp(['     ',SAVE.PointingString,saveDirectory])
end

%% EXIT CHECK
if exitAfterSaving
    error('exit after successful file saving...')
end


end
