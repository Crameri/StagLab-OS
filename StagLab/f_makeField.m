
%%                                          MAKE CUSTOM FIELD VARIABLE 1.42
% Regional Residual
% Horizontal Residual
% Horizontal-Band Residual
% Global Residual
% Dynamic Pressure
% Upwelling and downwelling
% Viscous Dissipation
%                                                Fabio Crameri, 15.08.2019
%
%% NOTES
% output variable VAR needs to consist of either VAR.var2d or
% VAR.var3d_yin and VAR.var3d_yang

% add YinYang calculation for residual temperature too - see f_DiagnosticsMantle commented lines....

function [VAR] = f_makeField(FIELD,FILE,GRID,SETUP,SWITCH,PLOT,MANTLE,TOPO)

%% DUMMIES
if ~isfield(FIELD,'name')
    clearvars FIELD
    FIELD.name  = 'undefined';
end

%% DEFINE TASK
if ~isfield(SWITCH,'customFieldTask')
    SWITCH.customFieldTask  = FIELD.name;
end


if strcmp(SWITCH.customFieldTask,'Regional residual')
    %% PLOT REGIONAL RESIDUAL-TEMPERATURE
    if strcmp(GRID.Type,'yinyang')
        if ~strcmp(FIELD.name,'Regional residual') %for the field, this is done later in f_plotField
            error('Regional Residual field is not implemented for yinyang grid.')
        end
        
    else %all other geometries
        %initiating
        VAR.var2d               = zeros(size(PLOT.T_3D));
        meanTregional           = VAR.var2d;
        VAR.meanTres            = VAR.var2d;
        VAR.minTres             = VAR.var2d;
        VAR.maxTres             = VAR.var2d;
        extentXgrid             = false;
        extentYgrid             = false;
        %flag periodic boundaries
        if strcmp(GRID.Type,'yinyang')
            extentYgrid         = true;
        elseif strcmp(GRID.Type,'Cartesian')
            if strcmp(GRID.Dim,'2-D')
                if size(PLOT.T_3D,1)>1 && size(PLOT.T_3D,2)==1 && strcmp(GRID.BCxx,'permeable'); extentXgrid = true; end
                if size(PLOT.T_3D,1)==1 && size(PLOT.T_3D,2)>1 && strcmp(GRID.BCyy,'permeable'); extentYgrid = true; end
            elseif strcmp(GRID.Dim,'3-D')
                if strcmp(GRID.BCxx,'permeable'); extentXgrid = true; end
                if strcmp(GRID.BCyy,'permeable'); extentYgrid = true; end
                warning(['Grid type ',GRID.Type,' in 3-D not implemented/tested yet!'])
            end
        elseif strcmp(GRID.Type,'spherical2D')
            if size(PLOT.T_3D,1)>1 && size(PLOT.T_3D,2)==1 && strcmp(GRID.BCxx,'permeable'); extentXgrid = true; end
            if size(PLOT.T_3D,1)==1 && size(PLOT.T_3D,2)>1 && strcmp(GRID.BCyy,'permeable'); extentYgrid = true; end
        else
            error(['Grid type ',GRID.Type,' not found!'])
        end
        %loop grid points
        for iz=1:size(PLOT.T_3D,3)
            if iz==1
                if strcmp(GRID.Type,'spherical2D')
                    X_3Dp              = GRID.X_3Dsp;
                    Y_3Dp              = GRID.Y_3Dsp;
                else
                    X_3Dp              = GRID.X_3Dp;
                    Y_3Dp              = GRID.Y_3Dp;
                end
            end
            %define horizontal extent via Rayleigh-Taylor Instability theory
            topLayerThickness       = max(0,max(GRID.Z_3Dp(:,:,iz)));
            bottomLayerThickness    = max(0,max(max(GRID.Z_3Dp(:))-max(GRID.Z_3Dp(:,:,iz))));
            b                       = min(topLayerThickness,bottomLayerThickness);     %top/bottom layer thickness
            lambda                  = 2.568 *b;                         %(see Turcotte & Schubert, Geodynamics 2nd ed. page 248)
            %get regional mean temperature
            for ix=1:size(PLOT.T_3D,1)
                for iy=1:size(PLOT.T_3D,2)
                    currentX                = X_3Dp(ix,iy,iz);
                    currentY                = Y_3Dp(ix,iy,iz);
                    currentZ                = GRID.Z_3Dp(ix,iy,iz);
                    lowerX                  = currentX-lambda;
                    lowerY                  = currentY-lambda;
                    lowerZ                  = currentZ-b;
                    higherX                 = currentX+lambda;
                    higherY                 = currentY+lambda;
                    higherZ                 = currentZ+b;
                    minX                    = min(X_3Dp(:,iy,iz)); %side boundaries
                    maxX                    = max(X_3Dp(:,iy,iz));
                    minY                    = min(Y_3Dp(ix,:,iz));
                    maxY                    = max(Y_3Dp(ix,:,iz));
                    Tx = []; Ty = []; Tz = [];
                    %get X-indices
                    if lowerX<minX && extentXgrid       %crossing side boundary
                        lowerX      	= maxX-(minX-lowerX);
                        [~,idxLowerX]   = min(abs(X_3Dp(:,iy,iz)-lowerX));
                        Tx            	= [Tx; PLOT.T_3D(1:ix,iy,iz); PLOT.T_3D(idxLowerX:end,iy,iz)];
                    else
                        [~,idxLowerX]   = min(abs(X_3Dp(:,iy,iz)-lowerX));
                        Tx            	= [Tx; PLOT.T_3D(idxLowerX:ix,iy,iz)];
                    end
                    if higherX>maxX && extentXgrid      %crossing side boundary
                        higherX         = minX+(higherX-maxX);
                        [~,idxHigherX]	= min(abs(X_3Dp(:,iy,iz)-higherX));
                        Tx            	= [Tx; PLOT.T_3D(ix:end,iy,iz); PLOT.T_3D(1:idxHigherX,iy,iz)];
                    else
                        [~,idxHigherX] 	= min(abs(X_3Dp(:,iy,iz)-higherX));
                        Tx            	= [Tx; PLOT.T_3D(ix:idxHigherX,iy,iz)];
                    end
                    %get Y-indices
                    if lowerY<minY && extentYgrid       %crossing side boundary
                        lowerY      	= maxY-(minY-lowerY);
                        [~,idxLowerY]   = min(abs(Y_3Dp(ix,:,iz)-lowerY));
                        Ty            	= [Ty; PLOT.T_3D(ix,1:iy,iz); PLOT.T_3D(ix,idxLowerY:end,iz)];
                    else
                        [~,idxLowerY]   = min(abs(Y_3Dp(ix,:,iz)-lowerY));
                        Ty            	= [Ty; PLOT.T_3D(ix,idxLowerY:iy,iz)];
                    end
                    if higherY>maxY && extentYgrid      %crossing side boundary
                        higherY         = minY+(higherY-maxY);
                        [~,idxHigherY]	= min(abs(Y_3Dp(ix,:,iz)-higherY));
                        Ty            	= [Ty; PLOT.T_3D(ix,iy:end,iz); PLOT.T_3D(ix,1:idxHigherY,iz)];
                    else
                        [~,idxHigherY] 	= min(abs(Y_3Dp(ix,:,iz)-higherY));
                        Ty            	= [Ty; PLOT.T_3D(ix,iy:idxHigherY,iz)];
                    end
                    %get Z-indices
                    [~,idxLowerZ]       = min(abs(GRID.Z_3Dp(ix,iy,:)-lowerZ));
                    Tz                  = [Tz; PLOT.T_3D(ix,iy,idxLowerZ:iz)];
                    [~,idxHigherZ]      = min(abs(GRID.Z_3Dp(ix,iy,:)-higherZ));
                    Tz                  = [Tz; PLOT.T_3D(ix,iy,iz:idxHigherZ)];
                    %compile
                    Tregional           = [Tx(:);Ty(:);Tz(:)];
                    meanTregional(ix,iy,iz)  	= mean(Tregional(:));   	%regional mean temperature
                    VAR.var2d(ix,iy,iz)     	= PLOT.T_3D(ix,iy,iz) - meanTregional(ix,iy,iz);     %residual temperature
                end
            end
        end
        for iz=1:size(PLOT.T_3D,3)
            %define horizontal extent via Rayleigh-Taylor Instability theory
            topLayerThickness       = max(0,max(GRID.Z_3Dp(:,:,iz)));
            bottomLayerThickness    = max(0,max(max(GRID.Z_3Dp(:))-max(GRID.Z_3Dp(:,:,iz))));
            b                       = min(topLayerThickness,bottomLayerThickness);     %top/bottom layer thickness
            lambda                  = 2.568 *b;                         %(see Turcotte & Schubert, Geodynamics 2nd ed. page 248)
            %get regional mean temperature
            for ix=1:size(PLOT.T_3D,1)
                for iy=1:size(PLOT.T_3D,2)
                    currentX                = X_3Dp(ix,iy,iz);
                    currentY                = Y_3Dp(ix,iy,iz);
                    currentZ                = GRID.Z_3Dp(ix,iy,iz);
                    lowerX                  = currentX-lambda;
                    lowerY                  = currentY-lambda;
                    lowerZ                  = currentZ-b;
                    higherX                 = currentX+lambda;
                    higherY                 = currentY+lambda;
                    higherZ                 = currentZ+b;
                    minX                    = min(X_3Dp(:,iy,iz)); %side boundaries
                    maxX                    = max(X_3Dp(:,iy,iz));
                    minY                    = min(Y_3Dp(ix,:,iz));
                    maxY                    = max(Y_3Dp(ix,:,iz));
                    Tx = []; Ty = []; Tz = [];
                    %get X-indices
                    if lowerX<minX && extentXgrid       %crossing side boundary
                        lowerX      	= maxX-(minX-lowerX);
                        [~,idxLowerX]   = min(abs(X_3Dp(:,iy,iz)-lowerX));
                        Tx            	= [Tx; VAR.var2d(1:ix,iy,iz); VAR.var2d(idxLowerX:end,iy,iz)];
                    else
                        [~,idxLowerX]   = min(abs(X_3Dp(:,iy,iz)-lowerX));
                        Tx            	= [Tx; VAR.var2d(idxLowerX:ix,iy,iz)];
                    end
                    if higherX>maxX && extentXgrid      %crossing side boundary
                        higherX         = minX+(higherX-maxX);
                        [~,idxHigherX]	= min(abs(X_3Dp(:,iy,iz)-higherX));
                        Tx            	= [Tx; VAR.var2d(ix:end,iy,iz); VAR.var2d(1:idxHigherX,iy,iz)];
                    else
                        [~,idxHigherX] 	= min(abs(X_3Dp(:,iy,iz)-higherX));
                        Tx            	= [Tx; VAR.var2d(ix:idxHigherX,iy,iz)];
                    end
                    %get Y-indices
                    if lowerY<minY && extentYgrid       %crossing side boundary
                        lowerY      	= maxY-(minY-lowerY);
                        [~,idxLowerY]   = min(abs(Y_3Dp(ix,:,iz)-lowerY));
                        Ty            	= [Ty; VAR.var2d(ix,1:iy,iz); VAR.var2d(ix,idxLowerY:end,iz)];
                    else
                        [~,idxLowerY]   = min(abs(Y_3Dp(ix,:,iz)-lowerY));
                        Ty            	= [Ty; VAR.var2d(ix,idxLowerY:iy,iz)];
                    end
                    if higherY>maxY && extentYgrid      %crossing side boundary
                        higherY         = minY+(higherY-maxY);
                        [~,idxHigherY]	= min(abs(Y_3Dp(ix,:,iz)-higherY));
                        Ty            	= [Ty; VAR.var2d(ix,iy:end,iz); VAR.var2d(ix,1:idxHigherY,iz)];
                    else
                        [~,idxHigherY] 	= min(abs(Y_3Dp(ix,:,iz)-higherY));
                        Ty            	= [Ty; VAR.var2d(ix,iy:idxHigherY,iz)];
                    end
                    %get Z-indices
                    [~,idxLowerZ]       = min(abs(GRID.Z_3Dp(ix,iy,:)-lowerZ));
                    Tz                  = [Tz; VAR.var2d(ix,iy,idxLowerZ:iz)];
                    [~,idxHigherZ]      = min(abs(GRID.Z_3Dp(ix,iy,:)-higherZ));
                    Tz                  = [Tz; VAR.var2d(ix,iy,iz:idxHigherZ)];
                    %compile
                    TresRegional            = [Tx(:);Ty(:);Tz(:)];
                    VAR.meanTres(:,:,iz) 	= mean(TresRegional(:));         %regional mean residual temperature
                    VAR.maxTres(:,:,iz)   	= max(TresRegional(:));
                    VAR.minTres(:,:,iz)  	= min(TresRegional(:));
                end
            end
        end
        %adjust array size
        if strcmp(GRID.Dim,'2-D')
            dummy(:,:)  = VAR.var2d(:,1,:); VAR.var2d = dummy;
        end
    end
    
    
elseif strcmp(SWITCH.customFieldTask,'Horizontal residual') || strcmp(SWITCH.customFieldTask,'Horizontal-band residual')
    %% PLOT HORIZONTAL RESIDUAL-TEMPERATURE
    %initiating
    VAR.var2d               = zeros(size(GRID.X_3D));
    meanThorizontal       	= VAR.var2d;
    VAR.meanTres            = VAR.var2d;
    VAR.minTres             = VAR.var2d;
    VAR.maxTres             = VAR.var2d;
    if strcmp(GRID.Type,'yinyang')
        VAR.var2d_yang      = VAR.var2d;
        if ~strcmp(FIELD.name,'Horizontal residual') && ~strcmp(FIELD.name,'Horizontal-band residual') %for the field, this is done later in f_plotField
            % CHECK FOR ONLY YINYANG IMPLEMENTATION --- WILL NEED TO BE UPDATED WITH OTHER OPTIONS
            if ~strcmp(SWITCH.customFieldTask,'Horizontal residual')
                warning('Horizontal-Band Residual field is not implemented yet for yinyang grid. - Switching to Horizontal-residual.')
                SWITCH.customFieldTask      = 'Horizontal residual';
            end
            if strcmp(SWITCH.customFieldTask,'Horizontal residual') %horizontal
                for iz=1:size(PLOT.T_3Dyin,3)
                    meanThorizontal(:,:,iz)   	= mean2([PLOT.T_3Dyin(:,:,iz);PLOT.T_3Dyang(:,:,iz)]);     %horizontal mean temperature
                    VAR.var2d(:,:,iz)           = PLOT.T_3Dyin(:,:,iz) - meanThorizontal(:,:,iz);            %horizontal temperature
                    VAR.var2d_yang(:,:,iz)      = PLOT.T_3Dyang(:,:,iz) - meanThorizontal(:,:,iz);           %horizontal temperature
                end
                for iz=1:size(PLOT.T_3Dyin,3)
                    dummy                       = [VAR.var2d(:,:,iz);VAR.var2d_yang(:,:,iz)];
                    VAR.meanTres(:,:,iz)      	= mean(dummy(:));   	%horizontal mean residual temperature
                    VAR.maxTres(:,:,iz)      	= max(dummy(:));
                    VAR.minTres(:,:,iz)        	= min(dummy(:));
                end

            elseif strcmp(SWITCH.customFieldTask,'Horizontal-band residual') %horizontal
                error('Horizontal-Band Residual field is not implemented for yinyang grid.')
            end
        end
        
    else %all other geometries
        if strcmp(SWITCH.customFieldTask,'Horizontal residual') %horizontal
            for iz=1:size(PLOT.T_3D,3)
                meanThorizontal(:,:,iz) 	= mean2(PLOT.T_3D(:,:,iz));         %horizontal mean temperature
                VAR.var2d(:,:,iz)           = PLOT.T_3D(:,:,iz) - meanThorizontal(:,:,iz);      %residual temperature
            end
            for iz=1:size(PLOT.T_3D,3)
                dummy                       = VAR.var2d(:,:,iz);
                VAR.meanTres(:,:,iz)        = mean(dummy(:));         %horizontal mean residual temperature
                VAR.maxTres(:,:,iz)         = max(dummy(:));
                VAR.minTres(:,:,iz)         = min(dummy(:));
            end
            
        elseif strcmp(SWITCH.customFieldTask,'Horizontal-band residual') %horizontal band
            bandThickness       = 10;
            for iz=1:size(PLOT.T_3D,3)
                zTop                        = min(size(PLOT.T_3D,3), iz+bandThickness);
                zBot                        = max(1, iz-bandThickness);
                meanThorizontal(:,:,iz)     = mean2(PLOT.T_3D(:,:,zBot:zTop));	%horizontal mean temperature
                VAR.var2d(:,:,iz)           = PLOT.T_3D(:,:,iz) - meanThorizontal(:,:,iz);      %residual temperature
            end
            for iz=1:size(PLOT.T_3D,3)
                zTop                        = min(size(PLOT.T_3D,3), iz+bandThickness);
                zBot                        = max(1, iz-bandThickness);
                dummy                       = VAR.var2d(:,:,zBot:zTop);
                VAR.meanTres(:,:,iz)        = mean(dummy(:));         %horizontal mean residual temperature
                VAR.maxTres(:,:,iz)         = max(dummy(:));
                VAR.minTres(:,:,iz)         = min(dummy(:));
            end
        end
        
        %adjust array size
        clearvars dummy
        if strcmp(GRID.Dim,'2-D')
            dummy(:,:)  = VAR.var2d(:,1,:); VAR.var2d = dummy;
        end
    end
    
    
elseif strcmp(SWITCH.customFieldTask,'Global residual')
    %% PLOT GLOBAL RESIDUAL-TEMPERATURE
    if strcmp(GRID.Type,'yinyang')
        if ~strcmp(FIELD.name,'Global residual') %for the field, this is done later in f_plotField
            Tmean_z                 = mean([PLOT.T_3Dyin(:);PLOT.T_3Dyang(:)]);     %global mean temperature
            VAR.var2d(:,:,:)        = PLOT.T_3Dyin(:,:,:) - Tmean_z;            %residual temperature
            VAR.var2d_yang(:,:,:) 	= PLOT.T_3Dyang(:,:,:) - Tmean_z;           %residual temperature
            
            dummy                   = [VAR.var2d(:);VAR.var2d_yang(:)];
            VAR.meanTres            = mean(dummy(:));   	%horizontal mean residual temperature
            VAR.maxTres             = max(dummy(:));
            VAR.minTres             = min(dummy(:));
        end
        
    else %all other geometries
        Tmean_z                 = mean(PLOT.T_3D(:));          %global mean temperature
        VAR.var2d(:,:,:)        = PLOT.T_3D(:,:,:) - Tmean_z;      %residual temperature
        
        VAR.meanTres            = mean(VAR.var2d(:));         %horizontal mean residual temperature
        VAR.maxTres             = max(VAR.var2d(:));
        VAR.minTres             = min(VAR.var2d(:));
        
        %adjust array size
        clearvars dummy
        if strcmp(GRID.Dim,'2-D')
            dummy(:,:)  = VAR.var2d(:,1,:); VAR.var2d = dummy;
        end
    end
    
    
elseif strcmp(SWITCH.customFieldTask,'Dynamic pressure')
    %% PLOT DYNAMIC PRESSURE
    %Total pressure minus hydrostatic pressure 
    %(including actual free-surface topography, if applied)
    
    % READ pressure -----------------------------------------------------------
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
    
    % READ density -----------------------------------------------------------
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Density';
    DATA.FieldAbbreviation      = 'RHO';
    DATA.StopExecutionIfNotFound = true;
    [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        RHO_3D = PLOT.RHO_3Dyin; RHO_3Dyang = PLOT.RHO_3Dyang;
    else %all other grid types
        RHO_3D = PLOT.RHO_3D;
    end
    
    % READ topography -----------------------------------------------------------
    %only if a free surface is present
    if strcmp(SETUP.topBC,'sticky-air') || strcmp(SETUP.topBC,'free')
        if strcmp(GRID.Type,'yinyang')
            error('Topo for yinyang not implemented yet.')
            
        else %all other geometries
            if ~isfield(TOPO,'topo2d') %&& ~isfield(TOPO,'z_surf')
                [TOPO]      = f_readTopo3D(TOPO,SWITCH,GRID,FILE);  %non-dim or dim [m]
            end
            
        end
    else
        TOPO.topo2d     = 0;
    end
    
    % CREATE NEW FIELD
    idxSealevel     	= find(abs(GRID.Z_3Dp(1,1,:))==min(abs(GRID.Z_3Dp(1,1,:)))); %find index close to sealevel
    idxSealevel         = idxSealevel-1;
    
    if strcmp(GRID.Type,'yinyang')
        PhydrostaticTopo_yin = zeros(size(RHO_3Dyin)); Phydrostatic_yin = PhydrostaticTopo_yin;
        PhydrostaticTopo_yang = zeros(size(RHO_3Dyin)); Phydrostatic_yang = PhydrostaticTopo_yin;
        RHOplate           	= mean(RHO_3Dyin(:,:,idxSealevel));
        RHOplate_yang    	= mean(RHO_3Dyang(:,:,idxSealevel));

        for iz=1:size(RHO_3Dyin,3)
            PhydrostaticTopo_yin(:,:,1) = RHOplate.*SETUP.g.*TOPO.topo2d_yin;
            Phydrostatic_yin(:,:,iz)  	= RHO_3Dyin(:,:,iz).*SETUP.g.*GRID.Z_3D(:,:,iz);
            
            PhydrostaticTopo_yang(:,:,1)= RHOplate_yang.*SETUP.g.*TOPO.topo2d_yang;
            Phydrostatic_yang(:,:,iz)  	= RHO_3Dyang(:,:,iz).*SETUP.g.*GRID.Z_3D(:,:,iz);
        end
        
        Phydrostatic_yin    = Phydrostatic+PhydrostaticTopo; %topography contribution (if free surface)
        
        Pdynamic            = P_3Dyin-Phydrostatic_yin;
        VAR.var3d_yin       = Pdynamic;
        
        %ADJUST THIS ACCORDING TO THE ONE BELOW.........
        
        
        Phydrostatic        = RHO_3D.*SETUP.g.*(GRID.Z_3D+TOPO.topo2d); %hydrostatic from flat surface
        Phydrostatic        = Phydrostatic+PhydrostaticTopo; %topography contribution (if free surface)
        
        Phydrostatic        = RHO_3Dintegrated_yang.*SETUP.g.*(GRID.Z_3D+TOPO.topo2d_yang);
        Pdynamic            = P_3Dyang-Phydrostatic;
        VAR.var3d_yang      = Pdynamic;
        
    else %all other geometries
        PhydrostaticTopo = zeros(size(RHO_3D)); Phydrostatic = PhydrostaticTopo;
        %HYDROSTATIC PRESSURE
        %idx=end is top and ind=1 is bottom
        %topography contribution
        RHOplate       	= mean(RHO_3D(:,:,idxSealevel));
        
        if logical(1) %summing method
            for iz=size(RHO_3D,3):-1:1
                %fill array with topo value
                PhydrostaticTopo(:,:,iz) = RHOplate.*SETUP.g.*(TOPO.topo2d/GRID.m2p);
                %pressure contribution of current depth layer
                if iz==size(RHO_3D,3) %top grid point
                    Phydrostatic(:,:,iz)  	= RHO_3D(:,:,iz).*SETUP.g.*GRID.dz(:,:,iz)/2; %Pressure node is in the middle of the cell
                else
                    Phydrostatic(:,:,iz)  	= RHO_3D(:,:,iz).*SETUP.g.*GRID.dz(:,:,iz);
                end
                if iz<size(RHO_3D,3)
                    %integrate all overlying layers and add to current one (i.e., sum current and above layer)
                    Phydrostatic(:,:,iz)  	= Phydrostatic(:,:,iz)+Phydrostatic(:,:,iz+1);
                end
            end
            
        else %mean method
            %mean of density needs to weighted for variable grid spacing though
            for iz=size(RHO_3D,3):-1:1
                %fill array with topo value
                PhydrostaticTopo(:,:,iz) = RHOplate.*SETUP.g.*(TOPO.topo2d/GRID.m2p);
                %pressure contribution of current depth layer
                Phydrostatic(:,:,iz)  	= mean(RHO_3D(:,:,iz:end),3).*SETUP.g.*GRID.Z_3D(:,:,iz);
            end
        end
        
%         Phydrostatic    = Phydrostatic+PhydrostaticTopo; %topography
%         contribution (if free surface) - not needed as density is already
%         adjusted for it
        
        Pdynamic      	= P_3D-Phydrostatic;
        VAR.var2d       = Pdynamic;
        %adjust array size
        clearvars dummy
        if strcmp(GRID.Dim,'2-D')
            dummy(:,:)  = VAR.var2d(:,1,:); VAR.var2d = dummy;
        end
    end
    


    
elseif strcmp(SWITCH.customFieldTask,'Upwelling and downwelling')
    %% PLOT UP- AND DOWN-WELLING
    %field consists of five numbers:
    %-2: active downwelling, -1: passive downwelling,
    %1: passive upwelling, 2: active upwelling, and
    %0: nothing
    
    % CREATE NEW FIELD
    if strcmp(GRID.Type,'yinyang')
        VAR.var3d_yin   = zeros(size(MANTLE.upWelling));
        VAR.var3d_yin 	= VAR.var3d_yin +MANTLE.passiveUpwelling;
        VAR.var3d_yin 	= VAR.var3d_yin +MANTLE.activeUpwelling;
        VAR.var3d_yin 	= VAR.var3d_yin -MANTLE.passiveDownwelling;
        VAR.var3d_yin 	= VAR.var3d_yin -MANTLE.activeDownwelling;
        
        error('check for the correct yang variables here!') %finish implementing yang variables here...........
        VAR.var3d_yang  = zeros(size(MANTLE.upWelling));
        VAR.var3d_yang 	= VAR.var3d_yang +MANTLE.passiveUpwelling;
        VAR.var3d_yang 	= VAR.var3d_yang +MANTLE.activeUpwelling;
        VAR.var3d_yang 	= VAR.var3d_yang -MANTLE.passiveDownwelling;
        VAR.var3d_yang 	= VAR.var3d_yang -MANTLE.activeDownwelling;
        
    else %all other geometries
        VAR.var2d       = zeros(size(MANTLE.upWelling));
        VAR.var2d       = VAR.var2d +1.*MANTLE.passiveUpwelling;
        VAR.var2d       = VAR.var2d +2.*MANTLE.activeUpwelling;
        VAR.var2d       = VAR.var2d -1.*MANTLE.passiveDownwelling;
        VAR.var2d       = VAR.var2d -2.*MANTLE.activeDownwelling;
        %adjust array size
        clearvars dummy
        if strcmp(GRID.Dim,'2-D')
            dummy(:,:)  = VAR.var2d(:,1,:); VAR.var2d = dummy;
        end
    end
    
    %TEST FIGURE
    % figure(2),clf
    % for i=1:4
    %     subplot(2,2,i)
    %     if i==1
    %         dummy(:,:) = MANTLE.passiveUpwelling(:,1,:);
    %     elseif i==2
    %         dummy(:,:) = MANTLE.activeUpwelling(:,1,:);
    %     elseif i==3
    %         dummy(:,:) = MANTLE.passiveDownwelling(:,1,:);
    %     elseif i==4
    %         dummy(:,:) = MANTLE.activeDownwelling(:,1,:);
    %     end
    %     pcolor(dummy')
    %     axis equal
    %     colorbar
    % end
    
elseif strcmp(SWITCH.customFieldTask,'Viscous dissipation')
    %STRESS
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Stress';
    DATA.FieldAbbreviation      = 'STR';
    [~,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        STR_3D = PLOT.STR_3Dyin; STR_3Dyang = PLOT.STR_3Dyang;
    else %all other grid types
        STR_3D = PLOT.STR_3D;
    end
    
    %STRAINRATE
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Strain rate';
    DATA.FieldAbbreviation      = 'EDOT';
    [~,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        EDOT_3D = PLOT.EDOT_3Dyin; EDOT_3Dyang = PLOT.EDOT_3Dyang;
    else %all other grid types
        EDOT_3D = PLOT.EDOT_3D;
    end

    %CALCULATE VISCOUS DISSIPATION (phi=sigma*edot)
    if strcmp(GRID.Type,'yinyang')
        VAR.var3d_yin   = STR_3D.*EDOT_3D;
        VAR.var3d_yang  = STR_3Dyang.*EDOT_3Dyang;
    else
        %VAR.var2d       = ETA_3D.*EDOT_3D.^2;   %different expression for viscous dissipation rate per unit volume [W/m^3]
        VAR.var2d       = STR_3D.*EDOT_3D;      %viscous dissipation rate per unit volume [W/m^3]
        
        %multiply by grid cell area
        %VAR.var2d       = VAR.var2d.*GRID.cellVolume;  %viscous dissipation rate [W/m]
        
        %adjust array size
        clearvars dummy
        if strcmp(GRID.Dim,'2-D')
            dummy(:,:)  = VAR.var2d(:,1,:); VAR.var2d = dummy;
        end
    end

end


% REMOVE LOCALLY USED FIELDS
SWITCH = rmfield(SWITCH,'customFieldTask');











