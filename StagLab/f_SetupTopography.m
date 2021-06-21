
%%                                                    SETUP TOPOGRAPHY 2.43
%
%                                                Fabio Crameri, 22.11.2018
%% NOTES
% topography component calculation not implemented for non-dimensional mode.

function [TOPO,PLOT,SAVE] = f_SetupTopography(TOPO,FIELD,GRID,SWITCH,SETUP,PLOT,FILE,STYLE,SAVE,PLATE)

%% SETUP VARIABLES
if TOPO.field && ~strcmp(FIELD.name,'Topography')
    plotCMB 	= false;
    plotSURF	= true;
    Tfieldname  = 'Topography';
elseif strcmp(FIELD.name,'Topography')
    plotCMB 	= TOPO.cmb;
    plotSURF	= TOPO.surf;
    Tfieldname  = 'Topography';
elseif strcmp(FIELD.name,'Dyn. topography') || strcmp(FIELD.name,'Iso. topography') || strcmp(FIELD.name,'Res. topography')
    plotCMB 	= false;
    plotSURF	= true;
    Tfieldname  = 'Topography';
else
    error('fc: unknown Tfieldname')
end

%% READ TOPOGRAPHY
if TOPO.ascii  %topo from ascii file (old stag versions)
    error('topography from ascii output file not supported anylonger.')
    if FILE.number<10000
        number_string = num2str(10000+FILE.number);
        number_string(1)='0';
    else
        number_string = num2str(10000+FILE.number);
    end
    filename    = strcat(FILE.directory,FILE.name,'_sc',number_string,'.dat');
    if ~exist(filename,'file')        %else try:  _fstopo00000.dat
        filename    = strcat(FILE.directory,FILE.name,'_fstopo',number_string,'.dat'); disp('     -> _fstopo.dat used for topography!');
    end
    if ~exist(filename,'file')  %else: _sc00000.dat
        filename    = strcat(FILE.directory,FILE.name,'_cs',number_string,'.dat'); disp('     -> _cs.dat used for topography!'); %try: _cs00000.dat
    end
    topo_dat = load(filename);
    
    %x = topo_dat(:,1);
    if strcmp(GRID.Dim,'2-D') %size(topo_dat,2)<4  %2-D
        topo2d(:,1) = topo_dat(:,3); %cmb topo
        topo2d(:,2) = topo_dat(:,2); %surf topo
    elseif strcmp(GRID.Dim,'3-D')
        %y = topo_dat(:,2);
        topo2d(:,1) = topo_dat(:,4); %cmb topo
        topo2d(:,2) = topo_dat(:,3); %surf topo
    end
    
else  %topo from binary outputfile (new stagyy versions)
    % Read topo information
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = Tfieldname;
    DATA.FieldAbbreviation      = 'TOPO';
    DATA.StopExecutionIfNotFound = true;
    [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        TOPO_3D = PLOT.TOPO_3Dyin; TOPO_3Dyang = PLOT.TOPO_3Dyang;
    else %all other grid types
        TOPO_3D = PLOT.TOPO_3D;
    end
    if strcmp(GRID.Dim,'2-D')
        topo2d(:,:)     = TOPO_3D(:,1,:);
    elseif strcmp(GRID.Dim,'3-D')
        topo2d(:,:,:)   = TOPO_3D(:,:,:);
    end
    %     %actually not needed - just to be sure:  %might only be correct for 2-D?
    %     topo2d(:,1) = topo2d(:,1)-mean(topo2d(:,1));    % remove mean value (CMB topo)
    %     topo2d(:,2) = topo2d(:,2)-mean(topo2d(:,2));    % remove mean value (surf topo)
end

%% HORIZONTAL GRID VARIABLES
%find index close to sealevel
% idxSealevel = find(abs(GRID.Z_3Dp(1,1,:))==min(abs(GRID.Z_3Dp(1,1,:))));
% if strcmp(GRID.Type,'Cartesian')
%     x_PLOT(:,:) 	= GRID.X_3Dp(:,1,idxSealevel);
%     y_PLOT(:,:) 	= GRID.Y_3Dp(1,:,idxSealevel);
% elseif strcmp(GRID.Type,'spherical2D')
%     x_PLOT(:,:) 	= GRID.X_3Dsp(:,1,end-idxSealevel); %account for not-flipped grid in non-Cartesian geometry
%     y_PLOT(:,:) 	= GRID.Y_3Dsp(1,:,end-idxSealevel); %account for not-flipped grid in non-Cartesian geometry
% end

%% VARIABLES FOR DIFFERENCE PLOT
nloopComp = 1;
if SWITCH.plotDifference
    if ~isfield(PLOT,'topo2d_A') %for odd cases/snapshots
        PLOT.topo2d_A = topo2d;
        return  %end execution of function here
    else
        if strcmp(TOPO.difference,'diff')
            topo2d_diff = PLOT.topo2d_A - topo2d;
            topo2d = topo2d_diff;
            clearvars topo2d_diff
            PLOT = rmfield(PLOT,'topo2d_A');
        elseif strcmp(TOPO.difference,'comp')
            nloopComp = 2;
        else
            error(['TOPO.difference = ',TOPO.difference,' not available!'])
        end
    end
end

%% ERROR CHECK/PREVENTION
if size(topo2d,2)<2  %array needs to have two columns
    topo2d  = [topo2d, topo2d];
end

for iloopComp = 1:nloopComp
    if SWITCH.plotDifference && strcmp(TOPO.difference,'comp')
        if iloopComp==1
            PLOT.topo2d_B = topo2d;
            topo2d = PLOT.topo2d_A;
        else
            topo2d = PLOT.topo2d_B;
        end
    end
    topoSmooth_2 = zeros(size(topo2d));
    topoSmooth_1 = topo2d;
    
    %% SMOOTHING
    %check if curve fitting toolbox is installed
    try
        smooth(1:5); %it is present
        SmoothingPossible       = true;
    catch me %it is not present
        warning(me.message)
        SmoothingPossible       = false;
        topoSmooth              = NaN;
    end
    if strcmp(GRID.Dim,'2-D')
        %------------------------
        % SMOOTHING 2-D
        %------------------------
        smoothingMethod = 3;
        if SmoothingPossible
            if smoothingMethod==1 %DOESN'T WORK FOR 2-D YET!!!!!!!!! TOPO BECOMES SMALL!?
                filterSize = 5; %has to be an odd number!!!!
                %F = [.05 .1 .05; .1 .4 .1; .05 .1 .05]; %filter
                %F = [1/9 1/9 1/9; 1/9 1/9 1/9; 1/9 1/9 1/9]; %like 2-D filter - with one loop!
                F = fspecial('gaussian',filterSize); %filterSize is the size of matrix
                F = F(:,ceil(end/2)); %convert for 2D
                
                %add ghost points to prevent side effects from smoothing
                dummy = topoSmooth_1; %surf and cmb
                for ii=1:filterSize
                    dummy = [dummy(1,:); dummy; dummy(end,:)];
                end
                topoSmooth_2 = zeros(size(dummy));
                %smoothing
                for jj=1:17
                    topoSmooth_2(:,1) = conv(dummy(:,1),F,'same');
                    topoSmooth_2(:,2) = conv(dummy(:,2),F,'same');
                    dummy = topoSmooth_2;
                end
                %remove ghostpoints
                topoSmooth_2 = dummy(filterSize+1:end-filterSize,:);
                
            elseif smoothingMethod==2
                for jj=1:2
                    topoSmooth_2(:,1) = smooth(topoSmooth_1(:,1),5,'loess');
                    topoSmooth_2(:,2) = smooth(topoSmooth_1(:,2),5,'loess');
                    topoSmooth_1 = topoSmooth_2;
                end
            elseif smoothingMethod==3
                for jj=1:2
                    topoSmooth_2(:,1) = smooth(topoSmooth_1(:,1),5,'moving');
                    topoSmooth_2(:,2) = smooth(topoSmooth_1(:,2),5,'moving');
                    topoSmooth_1 = topoSmooth_2;
                end
            end
        end
        topoSmooth = topoSmooth_2;
        %------------------------
        if (plotCMB && plotSURF)
            % nothing to be done
        elseif (plotCMB)
            topo2d(:,2)=[];
            topoSmooth(:,2)=[];
        elseif (plotSURF)
            topo2d(:,1)=[];
            topoSmooth(:,1)=[];
        end
        
    elseif strcmp(GRID.Dim,'3-D')
        if TOPO.ascii
            if (plotCMB && plotSURF)
                disp('only one topography plot available at once!! -fc')
                z = topo2d(:,2); %surf topo
            elseif (plotCMB)
                z = topo2d(:,1); %cmb topo
            elseif (plotSURF)
                z = topo2d(:,2); %surf topo
            end
            low = min(x); high = max(x); tix = low:0.01:high;
            low = min(y); high = max(y); tiy = low:0.01:high;
            
            [xi,yi] = meshgrid(tix,tiy);
            topoP = griddata(x,y,z,xi,yi);
        else
            if (plotCMB && plotSURF)
                disp('only one topography plot available at once!! -fc')
                topoP = topo2d(:,:,2); %surf topo
            elseif (plotCMB)
                topoP = topo2d(:,:,1); %cmb topo
            elseif (plotSURF)
                topoP = topo2d(:,:,2); %surf topo
            end
        end
        if SWITCH.plotTOPOsmooth
            %------------------------
            % SMOOTHING 3-D
            %------------------------
            smoothingMethod = 1;
            if smoothingMethod==1
                filterSize = 5;
                %F = [.05 .1 .05; .1 .4 .1; .05 .1 .05]; %filter
                %F = [1/9 1/9 1/9; 1/9 1/9 1/9; 1/9 1/9 1/9]; %like 2-D filter - with one loop!
                F = fspecial('gaussian',filterSize); %filterSize is the size of matrix
                %F = fspecial('gaussian');
                
                %add ghost points to prevent side effects from smoothing
                dummy = topoP;
                for ii=1:filterSize
                    dummy = [dummy(1,:); dummy; dummy(end,:)];
                    dummy = [dummy(:,1), dummy, dummy(:,end)];
                end
                %smoothing
                for jj=1:17
                    topoSmooth = conv2(dummy,F,'same');
                    dummy = topoSmooth;
                end
                %remove ghostpoints
                topoSmooth = dummy(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
                
            elseif smoothingMethod==2
                %according to smooth2a
                Nr = 3; %rows
                Nc = 3; %colums
                % You end up replacing element "i" by the mean of a (2*Nr+1)-by-
                % (2*Nc+1) rectangle centered on element "i".
                [row,col] = size(topoP);
                eL = spdiags(ones(row,2*Nr+1),(-Nr:Nr),row,row);
                eR = spdiags(ones(col,2*Nc+1),(-Nc:Nc),col,col);
                % Setting all "NaN" elements of "matrixIn" to zero so that these will not
                % affect the summation.  (If this isn't done, any sum that includes a NaN
                % will also become NaN.)
                A = isnan(topoP);
                topoP(A) = 0;
                % For each element, we have to count how many non-NaN elements went into
                % the sums.  This is so we can divide by that number to get a mean.  We use
                % the same matrices to do this (ie, "eL" and "eR").
                nrmlize = eL*(~A)*eR;
                nrmlize(A) = NaN;
                % Actually taking the mean.
                topoSmooth = eL*topoP*eR;
                topoSmooth = topoSmooth./nrmlize;
            end
        end
    end
    
    %% DIMENSIONALISATION
    if strcmp(GRID.Dim,'2-D')
        topoP           = topo2d        *GRID.dimFactor;
        topoSmoothP     = topoSmooth   	*GRID.dimFactor;
    else
        topoP           = topoP         *GRID.dimFactor;
        topoSmoothP    	= topoSmooth    *GRID.dimFactor;
    end
end

%% FUNCTION OUTPUT
% if strcmp(FIELD.name,'Topography')
    TOPO.topo2d               	= topoP;            %actual topo
    if SWITCH.plotTOPOsmooth
        TOPO.topo2dp          	= topoSmoothP;      %smoothed topo
    else
        TOPO.topo2dp         	= topoP;            %actual topo
    end
% else 
%     error('Field undefined: check here!')
% end





%% TOPOGRAPHY COMPONENTS (ISOSTATIC, RESIDUAL, DYNAMIC)

if ~TOPO.calculateComponents || ~SWITCH.DimensionalMode || (strcmp(SETUP.topBC,'sticky-air') && ~SWITCH.ActualDepth)
    if TOPO.calculateComponents
        %% ERROR CHECKS & WARNING MESSAGES
        if ~SWITCH.DimensionalMode
            warning('Topography-component routine only works in dimensional mode (i.e., SWITCH.DimensionalMode = true)!')
        end
        % this needs actual depth switched on
        if ~SWITCH.ActualDepth && strcmp(SETUP.topBC,'sticky-air')
            warning('Topography-component routine needs actual depth switched on (i.e., SWITCH.ActualDepth = true)!')
        end
    end
    
    %% DUMMY FUNCTION OUTPUT
    TOPO.topoIso2d      	= TOPO.topo2d*NaN;
    TOPO.topoRes2d      	= TOPO.topo2d*NaN;

else
    %% INPUT
    compensationDepth   = 100 *1e3;   %184 *1e3;  	%[m] point at which upper plate does not flow
    getRidgeTopography  = logical(0);
    dynamicTopography   = logical(1);
    DynTopoMethod       = 1;                        %1: from sigma_zz, 2: from residual pressure at plate base, 3:from integrated residual pressure at plate base
    
    %% DEFAULTS
    scaleByViscosity    = logical(0);   %scale normal (i.e., sigma_zz) by viscosity
    useActualPlateETA   = logical(0);   %use actual plate (surface/core) viscosity for scaling zz-stress values
    useActualPlateRHOa  = logical(1);   %use actual plate (integrated) density for calculating dynamic topography
    useActualPlateRHOb  = logical(0);   %use actual plate (surface/core) density for calculating dynamic topography
    
    %% READ TEMPERATURE ---------------------------------------------------
    % Read temperature field information
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Temperature';
    DATA.FieldAbbreviation      = 'T';
    DATA.StopExecutionIfNotFound = true;
    [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        T_3D        = PLOT.T_3Dyin; T_3Dyang   = PLOT.T_3Dyang;
    else %all other grid types
        T_3D        = PLOT.T_3D;
    end
    
    %% READ VISCOSITY ---------------------------------------------------
    % Read viscosity field information
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Viscosity';
    DATA.FieldAbbreviation      = 'ETA';
    DATA.StopExecutionIfNotFound = true;
    [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        ETA_3D        = PLOT.ETA_3Dyin; ETA_3Dyang   = PLOT.ETA_3Dyang;
    else %all other grid types
        ETA_3D        = PLOT.ETA_3D;
    end
    
    %% READ DENSITY ---------------------------------------------------
    % Read density field information
    DensityFieldAvailable = true;
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Density';
    DATA.FieldAbbreviation      = 'RHO';
    [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        RHO_3D        = PLOT.RHO_3Dyin; RHO_3Dyang   = PLOT.RHO_3Dyang;
    else %all other grid types
        RHO_3D        = PLOT.RHO_3D;
    end
    if DATA.NotFound; DensityFieldAvailable = false; end
    
    %% READ topography -----------------------------------------------------------
    if getRidgeTopography && ~isfield(TOPO,'topo2d') && ~isfield(TOPO,'z_surf')
        [TOPO]      = f_readTopo3D(TOPO,SWITCH,GRID,FILE);  %non-dim or dim [m]
    end
    %% READ pressure component ---------------------------------------------------
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Velocity';
    DATA.FieldAbbreviation      = 'V';
    DATA.StopExecutionIfNotFound = true;
    [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        P_3D = PLOT.P_3Dyin; P_3Dyang = PLOT.P_3Dyang;
    else %all other grid types
        P_3D = PLOT.P_3D;
    end    
    
    %% DUMMIES
    NSTRESS_3D = ones(size(P_3D)).*NaN;
    
    %% ERROR CHECKS
    if ~DensityFieldAvailable && DynTopoMethod==3
        DynTopoMethod       = 2; if SWITCH.Verbose; warning('DynTopoMethod changed to option 2. Check here!'); end
    end
    if ~DensityFieldAvailable && (useActualPlateRHOa || useActualPlateRHOb)
        useActualPlateRHOa  = false;
        useActualPlateRHOb  = false; if SWITCH.Verbose; warning('No actual density used for topography component calculation. Check here!'); end
    end
    
    if dynamicTopography && DynTopoMethod==1
        dynamicTopographyFailed = false;
        %% READ zz-stress component ---------------------------------------------------
        % Read zz-stress component information
        Stress2use  = 'sigma_zz';   %'sigma_zz' or 'principalStress'
        if strcmp(Stress2use,'sigma_zz')
            DATA.Task                   = 'ImportFieldData';
            DATA.Field2Import           = 'zz-Stress component';
            DATA.FieldAbbreviation      = 'NSTRESS';
            [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
            if strcmp(GRID.Type,'yinyang')
                NSTRESS_3D   	= PLOT.NSTRESS_3Dyin; NSTRESS_3Dyang   = PLOT.NSTRESS_3Dyang;
            else %all other grid types
                NSTRESS_3D     	= PLOT.NSTRESS_3D;
            end
            if DATA.NotFound; dynamicTopographyFailed = true; end
            
        elseif strcmp(Stress2use,'principalStress')
            DATA.Task                   = 'ImportFieldData';
            DATA.Field2Import           = 'Radial princ. stress';
            DATA.FieldAbbreviation      = 'NSTRESS';
            [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
            if strcmp(GRID.Type,'yinyang')
                NSTRESS_3D   	= PLOT.NSTRESS_3Dyin; NSTRESS_3Dyang   = PLOT.NSTRESS_3Dyang;
            else %all other grid types
                NSTRESS_3D     	= PLOT.NSTRESS_3D;
            end
            if DATA.NotFound; dynamicTopographyFailed = true; end
            
            NSTRESS_3D = -NSTRESS_3D;
        end
        %% DUMMIES
        if dynamicTopographyFailed
            NSTRESS_3D = ones(size(GRID.X_3D)).*NaN;
        end
        if ~exist('P_3D','var')
            P_3D = ones(size(GRID.X_3D)).*NaN;
        end
    end
    
    
    %% GRID VARIABLES
    %find index close to sealevel
%     idxSealevel = find(abs(GRID.Z_3Dp(1,1,:))==min(abs(GRID.Z_3Dp(1,1,:))));
%     if strcmp(GRID.Type,'Cartesian')
%         x2d(:,:) 	= GRID.X_3Dp(:,1,idxSealevel);   %in [plotting dimension]
%         y2d(:,:) 	= GRID.Y_3Dp(1,:,idxSealevel);
%     elseif strcmp(GRID.Type,'spherical2D')
%         x2d(:,:) 	= GRID.X_3Dsp(:,1,idxSealevel);   %in [plotting dimension]
%         y2d(:,:) 	= GRID.Y_3Dsp(1,:,idxSealevel);
%     end
    z2d(:,:,:) 	= GRID.Z_3Dp(:,:,:);   %in [plotting dimension]
    
    %% DERIVING PARAMETERS AND DIMENSIONALIZATION
    NSTRESS_3D  = NSTRESS_3D .*SETUP.stressscale;  %[MPa]
    P_3D        = P_3D .*SETUP.Pscale; %[GPa]
    T_3D        = T_3D .*SETUP.Tscale; %in [K]
    ETA_3D      = ETA_3D .*SETUP.etascale; %in [Pa s]
    RHO_3D      = RHO_3D .*SETUP.rhoscale; %in [kg/m^3]
    z2d         = z2d ./GRID.m2p./1e3; %depth in [km]
    
    if useActualPlateETA || useActualPlateRHOa || useActualPlateRHOb || logical(1) %use always for iso and res topography !!
        % define lithosphere (near) surface values
        deltaTsurface       = 0.08; %0.03: surface, 0.15: core  %fraction of Tmax
        T_PlateSurface      = min(T_3D(:))+deltaTsurface*max(T_3D(:));
        LITHO.Tisovalue     = T_PlateSurface;  %Temperature isovalue
        LITHO.zmax          = PLOT.lithoThickness(1,2); %maximum depth in [km]
        LITHO.zmin          = 20;                       %minimum depth in [km] (needed here to ensure no air in it)
        LITHO.gridpoints    = NaN;                      %define to switch on gridpoints output
        [LITHO] = f_PlateThickness(T_3D,z2d,LITHO);
        plateSurfaceIdxs    = LITHO.gridpoints;         %only used for this
        
        % define lithosphere core values
        T_PlateCore         = PLOT.lithoThickness(1,1)/2;
        LITHO.Tisovalue     = T_PlateCore;              %Temperature isovalue
        [LITHO] = f_PlateThickness(T_3D,z2d,LITHO);
        plateCoreIdxs       = LITHO.gridpoints;         %only used for this
    end
    
    % derive lithosphere thickness by T contour
    T_PlateBase         = PLOT.lithoThickness(1,1);
    LITHO.Tisovalue     = T_PlateBase;                  %Temperature isovalue
    LITHO.zmax          = PLOT.lithoThickness(1,2);     %depth in [km]
    LITHO.gridpoints    = NaN;                          %define to switch on gridpoints output
    [LITHO] = f_PlateThickness(T_3D,z2d,LITHO);
    
    % test plot depth levels
    if logical(0)
        figure(2)
        hold on
        a(:) = GRID.X_3Dp(:,1);
        b(:) = GRID.Z_3Dp(1,1,plateSurfaceIdxs);
        plot(a,b,'-b')
        a(:) = GRID.X_3Dp(:,1);
        b(:) = GRID.Z_3Dp(1,1,plateCoreIdxs);
        plot(a,b,'-k')
        a(:) = GRID.X_3Dp(:,1);
        b(:) = GRID.Z_3Dp(1,1,LITHO.gridpoints);
        plot(a,b,'-r')
        axis ij
        figure(1)
    end
    
    % get different field components at these depths
    nstressLithoBase    = zeros(size(LITHO.gridpoints,1),size(LITHO.gridpoints,2),1)*NaN;
    etaLithoBase        = nstressLithoBase;
    if useActualPlateETA; etaLithoTop = nstressLithoBase; end
    rhoLithoBase        = nstressLithoBase;
    if useActualPlateRHOa || useActualPlateRHOb; rhoLithoTop = nstressLithoBase; rhoLithoMid = nstressLithoBase; rhoPlateIntegrated = nstressLithoBase; end
%     x2d_nstress         = x2d;
%     y2d_nstress         = y2d;
    for ivalx=1:size(LITHO.gridpoints,1)  % THIS CAN BE OPTIMISED (CHECK DIAGNOSE PLATE ADJUSTMENTS!)
        for ivaly=1:size(LITHO.gridpoints,2)
            if isnan(LITHO.gridpoints(ivalx,ivaly))
                nstressLithoBase(ivalx,ivaly)   = NaN;
                etaLithoBase(ivalx,ivaly)       = NaN;
                rhoLithoBase(ivalx,ivaly)       = NaN;
%                 x2d_nstress(ivalx,ivaly)        = NaN; %with NaN values
%                 y2d_nstress(ivalx,ivaly)        = NaN; %with NaN values
                continue
            end
            nstressLithoBase(ivalx,ivaly)       = NSTRESS_3D(ivalx,ivaly,LITHO.gridpoints(ivalx,ivaly));
            etaLithoBase(ivalx,ivaly)           = ETA_3D(ivalx,ivaly,LITHO.gridpoints(ivalx,ivaly));
            if useActualPlateETA
                etaLithoTop(ivalx,ivaly)        = ETA_3D(ivalx,ivaly,plateSurfaceIdxs(ivalx,ivaly));
            end
            if DensityFieldAvailable
                rhoLithoBase(ivalx,ivaly)           = RHO_3D(ivalx,ivaly,LITHO.gridpoints(ivalx,ivaly));
                if useActualPlateRHOa || useActualPlateRHOb
                    rhoLithoTop(ivalx,ivaly)        = RHO_3D(ivalx,ivaly,plateSurfaceIdxs(ivalx,ivaly));
                    rhoLithoMid(ivalx,ivaly)        = RHO_3D(ivalx,ivaly,plateCoreIdxs(ivalx,ivaly));
                    rhoPlateIntegrated(ivalx,ivaly) = mean( RHO_3D(ivalx,ivaly,LITHO.gridpoints(ivalx,ivaly):plateSurfaceIdxs(ivalx,ivaly)) ); %[kg/m^3]
                end
            end
        end
    end
    
    %OTHER PARAMETERS
    g                       = SETUP.g;              %[m/s^2]
    nstressLithoBase        = nstressLithoBase*1e6; %[MPa]->[Pa]
    
    %DENSITY
    if useActualPlateRHOa || useActualPlateRHOb
        rhoPlate_base       = rhoLithoBase;     %smooth(rhoLithoBase,20,'moving')'; %[kg/m^3]
        rhoPlate_top        = rhoLithoTop;      %smooth(rhoLithoTop,20,'moving')'; %[kg/m^3]
        rhoPlate_mid        = rhoLithoMid;      %smooth(rhoLithoTop,20,'moving')'; %[kg/m^3]
        if useActualPlateRHOa
            rhoPlate      	= rhoPlateIntegrated;   %[kg/m^3] mean value between plate top and bot
        elseif useActualPlateRHOb
            rhoPlate      	= (rhoPlate_top+rhoPlate_base)/2; %[kg/m^3] mean value between plate top and bot
            %rhoPlate    	 = rhoPlate_mid; %[kg/m^3] mean value at plate core
        end
    else
        rhoPlate            = SETUP.rho0;           %[kg/m^3]
    end
    
    %VISCOSITY
    if useActualPlateETA
        %etaPlate            = etaLithoTop;
        etaPlate            = smooth(etaLithoTop,20,'moving')';
    else
        if scaleByViscosity; warning('No automatic detection of plate viscosity!'); end
        etaPlate          	= 1e23; %[Pa s] <<<<<<<<<<<<<<<<<<<<<<<<< MIGHT ADJUST PLATE VISCOSITY VALUE HERE
        % etaPlate           = mean(etaLithoBase(:));
    end
    
    %% CALCULATE TOPOGRAPHY COMPONENTS
    % calculate contribution to total P from plate's P, rho, and h at compensation depth
    
    % convert depth level to index and update actual depth value
    [~,idxCompDepth]        = min(abs(GRID.Z_3D(1,1,:)-compensationDepth)); %index of closest grid value
    compensationDepth       = GRID.Z_3D(1,1,idxCompDepth);      %[m] - update depth level to closest grid value
    
    %density values
    if isfield(PLATE,'UMantleDensity') && ~isnan(PLATE.UMantleDensity)
        rhoMantle           = PLATE.UMantleDensity;         	%[kg/m^3]
    else
        warning('No automatic detection of UM-density!')
        rhoMantle           = 3270;                             %[kg/m^3] some kind of mean <<<<<<<<<<<<<<<<<<<<<<<<< MIGHT ADJUST
    end
    rhoPlate;                                                   %[kg/m^3] some kind of mean
    rhoAir                  = 1;                                %[kg/m^3]
    
    %ridge topography level
    if getRidgeTopography
        [~,idx]             = min(abs(GRID.X_3Dp(:,1,1)-PLATE.Spreading(1,1)));
        topoRidge           = TOPO.topo2dp(idx,1).*1e3;
        P_topoRidge         = topoRidge *rhoMantle *g;
    end
    
    %other parameters
    topo2d                          = TOPO.topo2dp.*1e3;            %[m]
    thicknessPlate                  = LITHO.thicknessFull .*1e3;  	%[m] at every x-position <<<<<<<<<<<<<<<<<<<<<<<<< MIGHT ADJUST (currently without +ive topography)
    thicknessPlateUpToCompensationDepth	= min(compensationDepth,thicknessPlate);
    thicknessPlateUpToCompensationDepthWithTopo 	= thicknessPlateUpToCompensationDepth +topo2d; %[m] <<<<<<<<<<<<<<<<<<<<<<<<< MIGHT ADJUST (currently smooth topo)
    thicknessPlateWithTopo          = thicknessPlate +topo2d; %[m]
    positiveTopo                    = zeros(size(topo2d));
    negativeTopo                    = positiveTopo;
    positiveTopo(topo2d>0)          = topo2d(topo2d>0);
    negativeTopo(topo2d<0)          = topo2d(topo2d<0);
    compensationThickness           = compensationDepth+topo2d;
    
    
    
    %% ISOSTATIC & RESIDUAL TOPOGRAPHY
    %calculate plate's hydrostatic pressure contribution - COMPENSATION DEPTH SHOULD BE DEEPER THAN PLATE THICKNESS!!!!!!
    %deltaP = deltaRho*g*h
    deltaPh_Plate_cd       	= (rhoMantle-rhoPlate).*g.*(thicknessPlateUpToCompensationDepthWithTopo); %[Pa] - plate's pressure contribution at compensation depth
    %         deltaPh_Plate           = (rhoMantle-rhoPlate).*g.*(thicknessPlateWithTopo); %[Pa] - plate's pressure contribution
    deltaPh_Plate           = (rhoMantle-rhoPlate).*g.*thicknessPlate; %[Pa] - plate's pressure contribution
    
    %background hydrostatic mantle pressure (i.e., without plate)
    deltaPh_Mantle          = rhoMantle.*g.*(compensationThickness);  	%[Pa] - hydrostatic ridge pressure at compensation depth
    
    %total hydrostatic pressure at compensation depth
    %         Ph                     = rhoPlate.*g.*thicknessPlateWithTopo + rhoMantle.*g.*(compensationDepth-thicknessPlateUpToCompensationDepth); %[Pa] - isostatic (i.e., hydrostatic) pressure at compensation depth
    Ph                      = deltaPh_Plate_cd + deltaPh_Mantle;                %[Pa] - hydrostatic (i.e., hydrostatic) pressure at compensation depth
    
    %total pressure at compensation depth
    Ptotal                  = P_3D(:,:,idxCompDepth) .*1e9;           	%[Pa] - total (dynamic+hydrostatic) pressure at compensation depth
    
    %residual pressure at compensation depth
    Pres                    = Ptotal-Ph;                            %[Pa] - residual (~dynamic) pressure at compensation depth
    
    %residual pressure contribution of the plate
    %%%% THIS NEEDS BACKGROUND, HYDROSTATIC PRESSURE ADDED
    option = 4;
    if option==1
        deltaP_res          = Ptotal-deltaPh_Plate;
    elseif option==2
        deltaP_res      	= deltaPh_Plate+deltaPh_Mantle-Ptotal;                 	%[Pa] - residual (~dynamic) pressure contribution at compensation depth
        deltaP_res          = Ph-Ptotal;
    elseif option==3
        deltaP_res          = Ptotal -(rhoMantle-rhoPlate).*g.*(compensationDepth) -rhoMantle.*g.*(compensationDepth);
    end
    
    %isostatic topography component
    % P_isoPlate = rho*g*h    -->	h_0 = P_isoPlate/(rho*g) - depth between topo and compensation depth
    isoTopo2d   = zeros(size(thicknessPlate));
    if size(rhoPlate,1)==1 && size(rhoMantle,1)==1
        isoTopo2d       = ((rhoMantle-rhoPlate).*thicknessPlateWithTopo)./(rhoMantle-rhoAir); 	%[m] isostatic topography
    else
        isoTopo2d(rhoPlate<=rhoMantle) 	= ((rhoMantle-rhoPlate(rhoPlate<=rhoMantle)).*thicknessPlate(rhoPlate<=rhoMantle))./(rhoPlate(rhoPlate<=rhoMantle)-rhoAir);  	%[m] isostatic topography
        isoTopo2d(rhoPlate>rhoMantle) 	= ((rhoMantle-rhoPlate(rhoPlate>rhoMantle)).*thicknessPlateWithTopo(rhoPlate>rhoMantle))./(rhoMantle-rhoAir); 	%[m] isostatic topography
    end
    %     THIS IS THE VECTORISED VERSION OF THIS:
    %     for ix=1:size(thicknessPlate,1)
    %         for iy=1:size(thicknessPlate,2)
    %             if rhoPlate(ix,iy)<=rhoMantle
    %                 isoTopo2d(ix,iy)  	= ((rhoMantle-rhoPlate(ix,iy))*thicknessPlate(ix,iy))/(rhoPlate(ix,iy)-rhoAir);  	%[m] isostatic topography
    %             elseif rhoPlate(ix,iy)>rhoMantle
    %                 isoTopo2d(ix,iy)	= ((rhoMantle-rhoPlate(ix,iy))*thicknessPlateWithTopo(ix,iy))/(rhoMantle-rhoAir); 	%[m] isostatic topography
    %             end
    %         end
    %     end
   
    %set mean isostatic topography component equal zero (due to mass conservation of the incompressible model)
    if TOPO.normIsoTopoToZero
       isoTopo2d            = isoTopo2d-mean2(isoTopo2d);
    end

    %residual topography component
    if option==1
        resTopo2d           = deltaP_res./((rhoPlate-rhoAir).*g) -compensationDepth  ; 	%[m] residual topography
        %             resTopo2d           = deltaP_res./((rhoMantle).*g) -compensationDepth  ; 	%[m] residual topography
    elseif option==2
        %resTopo2d           = Pres./(rhoPlate.*g) -compensationThickness;  %[m] residual topography
        resTopo2d           = ( Pres./((rhoPlate-rhoAir).*g) )  ; 	%[m] residual topography
    elseif option==3
        resTopo2d           = deltaP_res./((rhoMantle).*g)  ; 	%[m] residual topography
    elseif option==4
        resTopo2d           = topo2d-isoTopo2d;
    end
    
    %conversion to km
    isoTopo2d               = isoTopo2d./1e3;       %[km]
    resTopo2d               = resTopo2d./1e3;       %[km]
    
    %check by adding two components, which should somewhat equal the effective topography
    %         isoTopo2d           = isoTopo2d+resTopo2d;
    
    
    %% DYNAMIC TOPOGRAPHY
    if dynamicTopography && ~dynamicTopographyFailed
        if DynTopoMethod==1
            % calculate dyn. topography from zz-stress component
            % sigma_zz = rho*g*h        -->	h_0 = sigma_zz/(rho*g)
            %                               h_dyn = h_0-z_measurementEFFECT
            
            % sigma_normal: stress at the base of the lithosphere (i.e., nstressLithoBase)
            % sigma_eff: stress at the base of the lithosphere (= nstressLithoBase ???)
            % rho: from density field or just characteristic value
            % h: depth of isoline
            
            % scale sigma_zz by viscosity (sigma = 2*eta*epsilon)
            % use viscosity difference to plate
            if scaleByViscosity
                eta_diff2plate      = etaPlate ./ etaLithoBase;
                nstressLithoBase    = nstressLithoBase .* eta_diff2plate;
            end
            
            % figure(2)
            % plot(x2d_nstress,nstressLithoBase)
            
            % calculate dyn. topography from zz-stress component
            % sigma_normal = rho*g*h    -->	h_0 = sigma_normal/(rho*g)
            dynTopo2d           = -nstressLithoBase./((rhoPlate-rhoAir).*g);  %[m]  %SHOULD BE DENSITY DIFFERENCE WITH AIR/WATER....
            dynTopo2d           = dynTopo2d./1e3;       %[km]
            
        elseif DynTopoMethod==2
            %total hydrostatic pressure at plate base MEAN DENSITY
            Ph_pb             	= rhoPlate.*g.*thicknessPlateWithTopo;  	%[Pa] - hydrostatic (i.e., hydrostatic) pressure at plate base
            
            %total pressure at plate base
            Ptotal_pb = zeros(size(LITHO.gridpointsFull));
            for ii=1:size(LITHO.gridpointsFull,1)
                for jj=1:size(LITHO.gridpointsFull,2)
                    Ptotal_pb(ii,jj) 	= P_3D(ii,jj,LITHO.gridpointsFull(ii,jj)) *1e9;    	%[Pa] - total (dynamic+hydrostatic) pressure at plate base
                end
            end
            %residual pressure at plate base
            Pres_pb           	= Ph_pb-Ptotal_pb;                       	%[Pa] - residual (~dynamic) pressure at plate base
            % calculate dynamic topography using residual pressure at plate base
            dynTopo2d           = ( Pres_pb./((rhoPlate_top-rhoAir).*g) )  ; 	%[m] residual topography
            dynTopo2d           = dynTopo2d./1e3;       %[km]
            
        elseif DynTopoMethod==3
            %total pressure at plate base
            Ptotal_pb = zeros(size(LITHO.gridpointsFull));
            for ii=1:size(LITHO.gridpointsFull,1)
                for jj=1:size(LITHO.gridpointsFull,2)
                    for kk=size(RHO_3D,3):-1:1
                        %check if plate thickness has been exceeded
                        if GRID.Z_3D(ii,jj,kk)>thicknessPlate(ii,jj)
                            Ph_pb(ii,jj)     	= Ph_pb(ii,jj)- RHO_3D(ii,jj,kk+1)*g*(GRID.dz(ii,jj,kk+1)/2); %remove half of the cell that was added last to account for pressure node in the middle of the cell
                            break
                        end
                        %total hydrostatic pressure at plate base INTEGRATED DENSITY
                        if kk==size(RHO_3D,3) %top grid point
                            Ph_pb(ii,jj)     	= RHO_3D(ii,jj,kk)*g*(GRID.dz(ii,jj,kk)/2);  	%pressure node in the middle of the cell
                        else
                            Ph_pb(ii,jj)     	= Ph_pb(ii,jj)+ RHO_3D(ii,jj,kk)*g*GRID.dz(ii,jj,kk);       %[Pa] - hydrostatic (i.e., hydrostatic) pressure at plate base
                        end
                    end
                    Ptotal_pb(ii,jj) 	= P_3D(ii,jj,LITHO.gridpointsFull(ii,jj)) *1e9;    	%[Pa] - total (dynamic+hydrostatic) pressure at plate base
                    % Ptotal_pb(ii,jj) 	= (P_3D(ii,jj,LITHO.gridpointsFull(ii,jj)) -(P_3D(ii,jj,LITHO.gridpointsFull(ii,jj))-P_3D(ii,jj,LITHO.gridpointsFull(ii,jj)+1))/2) *1e9;    	%[Pa] - total (dynamic+hydrostatic) pressure at plate base
                end
            end
            %residual pressure at plate base
            Pres_pb           	= Ph_pb-Ptotal_pb;                       	%[Pa] - residual (~dynamic) pressure at plate base
            % calculate dynamic topography using residual pressure at plate base
            dynTopo2d           = ( Pres_pb./((rhoPlate_top-rhoAir).*g) )  ; 	%[m] residual topography
            dynTopo2d           = dynTopo2d./1e3;       %[km]
        end
    end
    
    %% SMOOTHING
    for iComponent=1:3
        clearvars dummy topoSmooth_2 topoSmooth_1
        if iComponent==1
            dummy           = isoTopo2d;            %[km]
            NanLocations 	= find(isnan(isoTopo2d));
        elseif iComponent==2
            dummy           = resTopo2d;            %[km]
            NanLocations 	= find(isnan(resTopo2d));
        elseif iComponent==3
            if exist('dynTopo2d','var')
                dummy     	= dynTopo2d;            %[km]
                NanLocations= find(isnan(dynTopo2d));
            else
                continue
            end
        end
        if strcmp(GRID.Dim,'2-D') && TOPO.dynTOPOsmooth && SmoothingPossible
            topoSmooth_1 = dummy;
            for jj=1:2
                topoSmooth_2(:,1) = smooth(topoSmooth_1(:,1),5,'moving');
                topoSmooth_1 = topoSmooth_2;
            end
            dummy = topoSmooth_2;
        end
        %re-insert NaN values
        dummy(NanLocations) = NaN;
        
        if iComponent==1
            isoTopo2d      = dummy;            %[km]
        elseif iComponent==2
            resTopo2d      = dummy;            %[km]
        elseif iComponent==3
            dynTopo2d      = dummy;            %[km]
        end
    end
    
    
    %% FUNCTION OUTPUT
    TOPO.topoIso2d      	= isoTopo2d;      %actual or smoothed topo
    TOPO.topoRes2d      	= resTopo2d;      %actual or smoothed topo
    if exist('dynTopo2d','var')
        TOPO.topoDyn2d      	= dynTopo2d;      %actual or smoothed topo
    end
end

