
%%                                   DIAGNOSE SURFACE FIELD VARIATION 1.41
%
% after an original skript by Tobias Rolf
%                                                Fabio Crameri, 21.11.2019
%
%% NOTES
% This currently assumes constant grid spacing/area!

function [PLOT,SAVE] = f_DiagnosticsSurface(FIELD_M,FIELD,FILE,GRID,SWITCH,PLOT,STYLE,SAVE,PLATE)

%% INPUT
removeContinentalPart           = logical(0);           %Removes continental areas from surface age data
relateToOceanicSurfaceAreaOnly  = logical(1);           %Sets total area to oceanic area instead of total model surface

numberBins                  = PLOT.sfvNumberBins;
depthLevel                  = PLOT.sfvDepthLevel;     	%set depth to diagnose

%% DEFAULTS
normalisation               = 'probability';         	%'count', 'probability', 'pdf'

%% READ INPUT DATA  
if ~strcmp(PLOT.sfvField,'Velocity') &&...  %Scalar Field Data
        ~strcmp(PLOT.sfvField,'Horizontal velocity') &&...
        ~strcmp(PLOT.sfvField,'Radial velocity') &&...
        ~strcmp(PLOT.sfvField,'Pressure') &&...
        ~strcmp(PLOT.sfvField,'Horizontal princ. stress') &&...
        ~strcmp(PLOT.sfvField,'Radial princ. stress') &&...
        ~strcmp(PLOT.sfvField,'Toroidal')  &&...
        ~strcmp(PLOT.sfvField,'Poloidal')  &&...
        ~strcmp(PLOT.sfvField,'Streamfunction')
    
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = PLOT.sfvField;
    DATA.FieldAbbreviation      = 'VAR';
    [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        VAR_3D = PLOT.VAR_3Dyin; VAR_3Dyang = PLOT.VAR_3Dyang;
    else %all other grid types
        VAR_3D = PLOT.VAR_3D;
    end
    if DATA.NotFound
        PLOT.SFVunsuccessful = true; warning(VAR_3D);
        return
        
    end
    
else %Vector Field Data
    if ~strcmp(PLOT.sfvField,'Velocity') &&...  %Scalar Field Data
            ~strcmp(PLOT.sfvField,'Horizontal velocity') &&...
            ~strcmp(PLOT.sfvField,'Radial velocity')
        PLOT.SFVunsuccessful	= true;
        warning('not implemented yet: Check here!')
        return
        
    end
    
    %velocity
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Velocity';
    DATA.FieldAbbreviation      = 'V';
    [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        VX_3D = PLOT.VX_3Dyin; VX_3Dyang = PLOT.VX_3Dyang;
        VY_3D = PLOT.VY_3Dyin; VY_3Dyang = PLOT.VY_3Dyang;
        VZ_3D = PLOT.VZ_3Dyin; VZ_3Dyang = PLOT.VZ_3Dyang;
        V_3D = PLOT.V_3Dyin; V_3Dyang = PLOT.V_3Dyang;
    else %all other grid types
        VX_3D = PLOT.VX_3D; VY_3D = PLOT.VY_3D; VZ_3D = PLOT.VZ_3D; V_3D = PLOT.V_3D;
    end
    if DATA.NotFound
        PLOT.SFVunsuccessful = true; warning(VX_3D);
        return
        
    end
    if strcmp(PLOT.sfvField,'Velocity')
        VAR_3D      = V_3D;     %absolute velocity
        if strcmp(GRID.Type,'yinyang'); VAR_3Dyang = V_3Dyang; end
    elseif strcmp(PLOT.sfvField,'Horizontal velocity')
        VAR_3D      = sqrt(VX_3D.^2 +VY_3D.^2);
        if strcmp(GRID.Type,'yinyang'); VAR_3Dyang = sqrt(VX_3Dyang.^2 +VY_3Dyang.^2); end
    elseif strcmp(PLOT.sfvField,'Radial velocity')
        VAR_3D      = VZ_3D;
        if strcmp(GRID.Type,'yinyang'); VAR_3Dyang = VZ_3Dyang; end
    end
end

%% REMOVE CONTINENTAL PART IF NEEDED
if removeContinentalPart && strcmp(PLOT.sfvField,'Surface age')
    %load continental crust
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Cont. crust';
    DATA.FieldAbbreviation      = 'CC';
    [DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        CC_3D = PLOT.CC_3Dyin; CC_3Dyang = PLOT.CC_3Dyang;
    else %all other grid types
        CC_3D = PLOT.CC_3D;
    end
    if ~isnan(CC_3D(1))
        CCvalue         = CC_3D(:,:,end);       %get horizontal profile of continental crust
        CClocation      = zeros(size(CCvalue)); %apply threshold
        CClocation(CCvalue>0.2)     = 1;
        if strcmp(GRID.Type,'yinyang')
            CCvalue         = CC_3Dyang(:,:,end);       %get horizontal profile of continental crust
            CClocation_yang	= zeros(size(CCvalue)); %apply threshold
            CClocation_yang(CCvalue>0.2)	= 1;
        end
        clearvars CCvalue
    end
    
    %process surface age
    if strcmp(GRID.Type,'yinyang')
        SAGE_3D     = PLOT.SAGE_3Dyin; SAGE_3Dyang = PLOT.SAGE_3Dyang;
        if exist('CClocation','var')
            error('Surface age splitting into oceanic areas only is not implemented yet for yinyang geometry!')
        end
    else %all other grid types
        SAGE_3D     = PLOT.SAGE_3D;
        PLOT.SAGEcont_3D 	= zeros(size(SAGE_3D))*NaN;
        PLOT.SAGEocean_3D	= PLOT.SAGEcont_3D;
        if exist('CClocation','var')
            if relateToOceanicSurfaceAreaOnly
                PLOT.SAGEcont_3D   	= SAGE_3D(CClocation==1); %this would fully remove unwanted entries
                PLOT.SAGEocean_3D 	= SAGE_3D(CClocation==0);
            else
                PLOT.SAGEcont_3D(CClocation==1)   	= SAGE_3D(CClocation==1); %this leaves NaNs in unwanted areas
                PLOT.SAGEocean_3D(CClocation==0) 	= SAGE_3D(CClocation==0);
            end
        end
    end
    
    %Update field to process with wanted array
    VAR_3D      = PLOT.SAGEocean_3D;
end

%% ADJUST DEPTH LEVEL
depthLevel                  = depthLevel *GRID.nd2dim;
dummy                       = abs((GRID.Z_3Dp(1,1,:))-depthLevel);
[~,depthLevelIdx]           = min(dummy); %index of closest value
depthLevel                  = GRID.Z_3Dp(1,1,depthLevelIdx);
clearvars dummy

%% EXTRACT NECESSARY DATA
plateInteriorArea       = logical(1);

PLOT.KindOfDataDummy            = 'na';
if (size(VAR_3D,2)<10 && size(VAR_3D,3)<10) || (size(VAR_3D,1)<10 && size(VAR_3D,3)<10) %Graph data
    PLOT.KindOfDataDummy        = 'graph';
    if size(VAR_3D,3)==2
       VAR_3D                   = VAR_3D(:,:,2);           %choose surface data (e.g., for topography)
    end
    DataP                       = VAR_3D;

elseif plateInteriorArea && exist('PLATE','var') && isfield(PLATE,'coreLevelIdx') %Graph data from field data along plate core
    %PLATE.coreStrainrate etc. are already saved (remove deletion in plate diagnostics if want to use)
    coreLevelDummy              = reshape(1:numel(PLATE.coreLevelIdx),size(PLATE.coreLevelIdx))+numel(PLATE.coreLevelIdx)*(PLATE.coreLevelIdx-1);
    plateCoreValues             = VAR_3D(coreLevelDummy);
    dummy                       = GRID.Z_3Dp(coreLevelDummy);
    depthLevel                  = median(dummy(:));
    dummy                       = abs((GRID.Z_3Dp(1,1,:))-depthLevel);
    [~,depthLevelIdx]           = min(dummy(:)); %index of closest value
    depthLevel                  = GRID.Z_3Dp(1,1,depthLevelIdx);
    clearvars dummy
    PLOT.KindOfDataDummy        = 'plateCoreGraph';
    DataP                       = plateCoreValues;
    
else %Field data
    if SWITCH.Verbose && plateInteriorArea && (~exist('PLATE','var') || ~isfield(PLATE,'coreLevelIdx'))
        warning('Plate Diagnostics was not successful!')
    end
    PLOT.KindOfDataDummy        = 'field';
    % take horizontally-flat slice
    DataP(:,:)                  = VAR_3D(:,:,depthLevelIdx);
    
end

%% DIMENSIONALISATION
fieldDim                = FIELD_M{strcmp(FIELD_M(:,2),PLOT.sfvField),4};
varScale                = FIELD_M{strcmp(FIELD_M(:,2),PLOT.sfvField),5};

DataP                   = DataP .*varScale;

%% COMPUTE TOTAL SURFACE AREA
if strcmp(GRID.Type,'spherical2D')
    xExtent                 = abs(max(GRID.X_3Dsp(:,:,depthLevelIdx)) - min(GRID.X_3Dsp(:,:,depthLevelIdx)));
    yExtent                 = abs(max(GRID.Y_3Dsp(:,:,depthLevelIdx)) - min(GRID.Y_3Dsp(:,:,depthLevelIdx)));
else
    xExtent                 = abs(max(GRID.X_3Dp(:,:,depthLevelIdx)) - min(GRID.X_3Dp(:,:,depthLevelIdx)));
    yExtent                 = abs(max(GRID.Y_3Dp(:,:,depthLevelIdx)) - min(GRID.Y_3Dp(:,:,depthLevelIdx)));
end
if xExtent==0
    totalSurfaceArea    = yExtent;
elseif yExtent==0
    totalSurfaceArea    = xExtent;
else
    totalSurfaceArea    = xExtent*yExtent;
end

%% DIAGNOSE DATA
DataMinValue            = min(DataP(:));
DataMaxValue            = max(DataP(:));
DataMeanValue          	= mean(DataP(:));
DataMedianValue     	= median(DataP(:));
DataStandardDeviation  	= std(DataP(:));

%% PLOTTING
plotAdditionLegend      = logical(1);
faceColor               = [0.2 0.2 0.2];
edgeColor               = [1 1 1];
LineMeanColor          	= [1 1 1];
AreaColor               = [0.88 0.88 0.88];
AreaAlpha               = 0.8;
if strcmp(STYLE.ColorMode,'dark')
    faceColor           = 1-faceColor;
    edgeColor           = 1-edgeColor;
    LineMeanColor    	= 1-LineMeanColor;
    AreaColor           = 1-AreaColor;
end

if isnumeric(PLOT.sfvXmin) && isnumeric(PLOT.sfvXmax)
    xMin    = PLOT.sfvXmin;
    xMax    = PLOT.sfvXmax;
elseif isnumeric(PLOT.sfvXmin)
    xMin    = PLOT.sfvXmin;
    xMax    = max(DataP(:));
elseif isnumeric(PLOT.sfvXmax)
    xMin    = min(DataP(:));
    xMax    = PLOT.sfvXmax;
else
    xMin    = min(DataP(:));
    xMax    = max(DataP(:));
end
xSpan       = xMax-xMin;
dBin        = xSpan/(numberBins);

if PLOT.sfvShowOutliers
    upperLim = max(DataP(:));
else
    upperLim = xMax;
end

hHist = histogram(DataP,[xMin:dBin:xMax-dBin,upperLim],...
    'FaceColor',faceColor,'FaceAlpha',0.8,'EdgeColor',edgeColor,'EdgeAlpha',0.8,'Normalization',normalisation);

%constant axis
if SWITCH.Colorbar(1,2)
    ylim(gca,[FIELD.cColorbarMIN,FIELD.cColorbarMAX])
end
if isnumeric(PLOT.sfvXmin) || isnumeric(PLOT.sfvXmax)
    xlim(gca,[xMin,xMax])
end

%plot additions
if PLOT.sfvAdditions
    hold on
    yAxisLim    = ylim;
    %standard deviation
    hStdDev = fill([DataMeanValue-DataStandardDeviation,DataMeanValue-DataStandardDeviation,...
        DataMeanValue+DataStandardDeviation,DataMeanValue+DataStandardDeviation],...
        [yAxisLim,fliplr(yAxisLim)],AreaColor,'EdgeColor','none','FaceAlpha',AreaAlpha);
    %mean value
    try
        LineMeanColor = [LineMeanColor,0.8]; %transparency
        hMean = plot([DataMeanValue,DataMeanValue],yAxisLim,'Color',LineMeanColor,'LineWidth',1.7);
    catch
        LineMeanColor(1,4) = [];
        hMean = plot([DataMeanValue,DataMeanValue],yAxisLim,'Color',LineMeanColor,'LineWidth',1.7);
    end
    %median value
    hMedian = plot([DataMedianValue,DataMedianValue],yAxisLim,'Color',LineMeanColor,'LineWidth',1.7,'LineStyle',':');
    %put histogram to the top
    uistack(hHist,'top');
    if plotAdditionLegend
        %legend
        legend([hStdDev,hMean,hMedian],{'Std. Deviation','Mean','Median'});
    end
end

%% ANNOTATIONS
DESVARIA.xlabelName         = PLOT.sfvField;
DESVARIA.xlabelDim          = fieldDim;
if strcmp(normalisation,'probability') || strcmp(normalisation,'pdf')
    DESVARIA.zlabelName   	= 'Probability';
    DESVARIA.zlabelDim    	= '';
elseif strcmp(normalisation,'counts')
    DESVARIA.zlabelName    	= 'Frequency';
    DESVARIA.zlabelDim     	= '';
else
    DESVARIA.zlabelName   	= '';
    DESVARIA.zlabelDim    	= '';
    warning('Normalisation not recognised: Check here!')
end
if strcmp(GRID.Dim,'3-D')    
    DESVARIA.ylabelName3D   = PLOT.sfvField;
    DESVARIA.ylabelDim      = fieldDim;
    DESVARIA.xlabelName3D  	= DESVARIA.zlabelName;
    DESVARIA.xlabelDim  	= DESVARIA.zlabelDim;
end
DESVARIA.Task	= 'create annotation strings';
[PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);

title(PLOT.titleStringCurrent)
xlabel(PLOT.xlabel)
ylabel(PLOT.ylabel)
set(gca,'YGrid','on')
set(gca,'TickDir',SWITCH.TickDirection)

DESVARIA.Task 	= 'make title';
[PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);

%% FUNCTION OUTPUT
PLOT.sfvTotalSurfaceArea        = totalSurfaceArea;
PLOT.sfvDataMinValue            = DataMinValue;
PLOT.sfvDataMaxValue            = DataMaxValue;
PLOT.sfvDataMeanValue           = DataMeanValue;
PLOT.sfvDataMedianValue         = DataMedianValue;
PLOT.sfvDataStandardDeviation 	= DataStandardDeviation;
PLOT.sfvDepthLevel              = depthLevel;       %update depth level in plotting dimension
