
%%                                                  PLOT CUSTOM GRAPH 4.21
%    . calls f_DesignColourmap
%                                                Fabio Crameri, 17.06.2021
%
%% NOTES
% only implemented for 2-D geometry, yet.
% absolute values not perfect yet, as it might inhibit normal negative values....

function [PLOT,SAVE] = f_plotCustomGraph(iPlotLocal,FIELD,GRID,SWITCH,SETUP,PLOT,FILE,STYLE,SAVE)

%% DEFAULTS
data2import                 = 'Geodynamics';
plotYmean                   = false;
indicateCurrentTime         = false;
indicateCurrentTimeAsDot    = false;        %else the lines will have different shading
normaliseYdata2median       = false;
normaliseYdata2min          = false;
normaliseYdata2max          = false;
absoluteValues              = false;
plotAcceleration            = false;
plotSinkingOverTrenchVel    = false;
graphLineStyle              = '-';          %'o' or '-'
colouredDots              	= false;

%% INPUT
graphNameFound              = true;
if strcmp(PLOT.CustomGraphName(iPlotLocal),'custom')
    dataNameX                   = PLOT.CustomDataNameX;
    dataNameY                   = PLOT.CustomDataNameY;
    numMin                      = PLOT.CustomDataNumMin;
    numMax                      = PLOT.CustomDataNumMax;
    indicateCurrentTime         = PLOT.CustomIndicateCurrentTime;
    DESVARIA.fieldName         	= PLOT.CustomFieldName;
    absoluteValues              = PLOT.CustomAbsoluteValues;
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	= PLOT.CustomYAxisMin;
    graphYAxisMax             	= PLOT.CustomYAxisMax;
    graphXAxisMin           	= PLOT.CustomXAxisMin;
    graphXAxisMax            	= PLOT.CustomXAxisMax;    
 
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'SlabSinking vs. TrenchRetreat')
    dataNameX                   = {'Trench velocity'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Slab sinking velocity'};      %see SL_FieldPlot for available fields
    plotAcceleration            = logical(0);
    numMin                      = 1;
    numMax                      = 20;
    indicateCurrentTime         = logical(0);
    DESVARIA.fieldName       	= 'Trench dynamics';   	%title
    absoluteValues              = logical(0);               %plots only absolute values
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	=   -10.0;
    graphYAxisMax             	=   10.0;
%     graphXAxisMin           	=   0;
%     graphXAxisMax            	=   40;
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'PlateVelocities vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Trench velocity', 'Theoretic trench velocity', 'Lower-plate velocity', 'Convergence velocity', 'Slab sinking velocity'};      %see SL_FieldPlot for available fields
%     dataNameY                   = {'Trench velocity'};      %see SL_FieldPlot for available fields
    plotAcceleration            = logical(0);
    numMin                      = 1;
    numMax                      = 999;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Plate velocities';   	%title
    absoluteValues              = logical(0);               %plots only absolute values
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	=   -10.0;
    graphYAxisMax             	=   10.0;
    graphXAxisMin           	=   0;
    graphXAxisMax            	=   55;


elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'PlateVelocities vs. SlabDepth')
    dataNameX                   = {'Slab-tip depth'};       %e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Trench velocity', 'Theoretic trench velocity', 'Lower-plate velocity', 'Convergence velocity', 'Slab sinking velocity'};      %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 999;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Slab-UP interaction';   	%title
    interpolateRoughGridData    = logical(1);               %applied only for some parameters
    graphXAxisMin            	=   300;
    graphXAxisMax            	=   1000;
    graphYAxisMin           	=   -10;
    graphYAxisMax            	=   10;
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'UPtilt vs. SlabDepth')
    dataNameX                   = {'Slab-tip depth'};       %e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Upper-plate tilt'};      %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 120;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Slab-UP interaction';   	%title
    interpolateRoughGridData    = logical(1);               %applied only for some parameters
    graphXAxisMin            	=   300;
    graphXAxisMax            	=   1000;
    graphYAxisMin           	=   0;
    graphYAxisMax            	=   0.1;
    plotYmean                   = logical(0);
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'Inundation vs. SlabDepth')
    dataNameX                   = {'Slab-tip depth'};       %e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Inundation'};      %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 120;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName       	= 'Slab-UP interaction';   	%title
    interpolateRoughGridData    = logical(1);               %applied only for some parameters
    graphXAxisMin            	=   300;
    graphXAxisMax            	=   1000;
    graphYAxisMin           	=   0;
    graphYAxisMax            	=   1500;
    plotYmean                   = logical(1);
    normaliseYdata2min          = true;
 
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'UPtilt vs. Time')
    dataNameX                   = {'Time'};       %e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Upper-plate tilt'};      %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 999;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Slab-UP interaction';   	%title
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphXAxisMin            	=   0;
    graphXAxisMax            	=   55;
    graphYAxisMin           	=   0;
    graphYAxisMax            	=   0.1;
    plotYmean                   = logical(0);
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'Inundation vs. Time')
    dataNameX                   = {'Time'};       %e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Inundation'};      %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 120;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Slab-UP interaction';   	%title
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphXAxisMin            	=   0;
    graphXAxisMax            	=   55;
    graphYAxisMin           	=   0;
    graphYAxisMax            	=   1500;
    plotYmean                   = logical(0);
    normaliseYdata2min          = true;
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'ShallowSlabAngle vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Shallow-slab angle'};      %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 999;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Slab-dip angle';   	%title
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin           	=   20;
    graphYAxisMax            	=   90;
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'ViscosityContrast vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Slab-mantle visc. contrast'};      %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 20;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName         	= 'Slab-mantle visc. contrast';   	%title
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphXAxisMin            	=   0;
    graphXAxisMax            	=   10;
    graphYAxisMin           	=   0;
    graphYAxisMax            	=   200;
   
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'Sinking/Trench Vel. vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Slab sinking velocity','Trench velocity'};      %see SL_FieldPlot for available fields
    plotSinkingOverTrenchVel    = true;
    numMin                      = 1;
    numMax                      = 20;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Sinking/trench vel.';   	%title
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphXAxisMin            	=   0;
    graphXAxisMax            	=   10;
    graphYAxisMin           	=   0;
    graphYAxisMax            	=   3;
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'ViscContrast vs. Sinking/Trench Vel.')
    dataNameX                   = {'Slab sinking velocity','Trench velocity'};       %e.g., 'Time', see SL_FieldPlot for available fields
    plotSinkingOverTrenchVel    = true;
    dataNameY                   = {'Slab-mantle visc. contrast'};       %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 20;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName         	= 'Slab-mantle visc. contrast';   	%title
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphXAxisMin            	=   -3;
    graphXAxisMax            	=   5;
    graphYAxisMin           	=   0;
    graphYAxisMax            	=   200;
    graphLineStyle              = 'o';          %'o' or '-'
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'ShallowSlabAngle vs. Sinking/Trench Vel.')
    dataNameX                   = {'Slab sinking velocity','Trench velocity'};       %e.g., 'Time', see SL_FieldPlot for available fields
    plotSinkingOverTrenchVel    = true;
    dataNameY                   = {'Shallow-slab angle'};       %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 20;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName         	= 'Slab-dip angle';   	%title
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphXAxisMin            	=   -3;
    graphXAxisMax            	=   5;
    graphYAxisMin           	=   0;
    graphYAxisMax            	=   90;
    graphLineStyle              = 'o';          %'o' or '-'

elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'TrenchDepth vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Trench depth','Back-arc basin depth','Island-arc height'};      %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 20;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Trench depth';   	%title
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
%     graphXAxisMin            	=   0;
%     graphXAxisMax            	=   10;
    graphYAxisMin           	=   -6;
    graphYAxisMax            	=   3;
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'TrenchDepth vs. Sinking/Trench Vel.')
    dataNameX                   = {'Slab sinking velocity','Trench velocity'};       %e.g., 'Time', see SL_FieldPlot for available fields
    plotSinkingOverTrenchVel    = true;
    dataNameY                   = {'Trench depth'};       %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 20;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Trench depth';   	%title
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphXAxisMin            	=   0;
    graphXAxisMax            	=   3;
    graphYAxisMin           	=   -10;
    graphYAxisMax            	=   1;
    graphLineStyle              = 'o';          %'o' or '-'
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'MaxPlateStress vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Max. plate stress'};      %see SL_FieldPlot for available fields
    numMin                      = 0;
    numMax                      = 86;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Max. plate stress';   	%title
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
%     graphXAxisMin            	=   0;
%     graphXAxisMax            	=   10;
%     graphYAxisMin           	=   0;
%     graphYAxisMax            	=   200;

elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'MaxYieldDepth vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Max. yield depth fraction'};      %see SL_FieldPlot for available fields
    dataNameY                   = {'Max. yield depth'};      %see SL_FieldPlot for available fields
    numMin                      = 1;
    numMax                      = 52;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Max. yield depth';   	%title
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
%     graphXAxisMin            	=   0;
%     graphXAxisMax            	=   10;
%     graphYAxisMin           	=   0;
%     graphYAxisMax            	=   200;

elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'LLSVP Locations vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'LLSVP horiz. locations'};      %see SL_FieldPlot for available fields
    plotAcceleration            = logical(0);
    numMin                      = 10;
    numMax                      = 999;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'LLSVP locations';   	%title
    absoluteValues              = logical(0);               %plots only absolute values
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	= min(GRID.X_3Dp(:));
    graphYAxisMax             	= max(GRID.X_3Dp(:));
%     graphXAxisMin           	=   0;
%     graphXAxisMax            	=   55;
    graphLineStyle              = 'o';          %'o' or '-'
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'Continent Locations vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Continent horiz. locations'};      %see SL_FieldPlot for available fields
    plotAcceleration            = logical(0);
    numMin                      = 1;
    numMax                      = 250;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Continent locations';   	%title
    absoluteValues              = logical(0);               %plots only absolute values
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	= 4700;
    graphYAxisMax             	= 5300;
%     graphYAxisMin              	= min(GRID.X_3Dp(:));
%     graphYAxisMax             	= max(GRID.X_3Dp(:));
%     graphXAxisMin           	=   0;
%     graphXAxisMax            	=   55;
    graphLineStyle              = '-';          %'o' or '-'
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'Continent Motion vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'v_{h,mean} Continent'};      %see SL_FieldPlot for available fields
    plotAcceleration            = logical(0);
    numMin                      = 1;
    numMax                      = 250;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= 'Continent velocity';   	%title
    absoluteValues              = logical(0);               %plots only absolute values
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	= -1;
    graphYAxisMax             	= 4;
%     graphXAxisMin           	=   0;
%     graphXAxisMax            	=   55;
    graphLineStyle              = '-';          %'o' or '-'
    
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'Number plumes vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Number plumes'};      %see SL_FieldPlot for available fields
    plotAcceleration            = logical(0);
    numMin                      = 1;
    numMax                      = 999;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= '';   	%title
    absoluteValues              = logical(0);               %plots only absolute values
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	= 0;
%     graphYAxisMax             	= 5;
%     graphXAxisMin           	=   0;
%     graphXAxisMax            	=   55;
    graphLineStyle              = '-';          %'o' or '-'

elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'Mid-mantle plumes vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Mid-mantle plumes'};      %see SL_FieldPlot for available fields
    plotAcceleration            = logical(0);
    numMin                      = 1;
    numMax                      = 999;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= '';   	%title
    absoluteValues              = logical(0);               %plots only absolute values
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	= 0;
%     graphYAxisMax             	= 5;
%     graphXAxisMin           	=   0;
%     graphXAxisMax            	=   55;
    graphLineStyle              = '-';          %'o' or '-'
 
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'Horiz. velocity vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'v_{h,mean} UM','v_{h,mean} MM','v_{h,mean} LM'};      %see SL_FieldPlot for available fields
    plotAcceleration            = logical(0);
    numMin                      = 1;
    numMax                      = 999;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= '';   	%title
    absoluteValues              = logical(0);               %plots only absolute values
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	= 0;
%     graphYAxisMax             	= 5;
%     graphXAxisMin           	=   0;
%     graphXAxisMax            	=   55;
    graphLineStyle              = '-';          %'o' or '-'
  
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'Surface age vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Mean surface age','Median surface age'};      %see SL_FieldPlot for available fields
    plotAcceleration            = logical(0);
    numMin                      = 1;
    numMax                      = 999;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= '';   	%title
    absoluteValues              = logical(0);               %plots only absolute values
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	= 0;
%     graphYAxisMax             	= 5;
%     graphXAxisMin           	=   0;
%     graphXAxisMax            	=   55;
    graphLineStyle              = '-';          %'o' or '-'
  
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'Crustal thickness vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Mean crustal thickness','Median crustal thickness'};      %see SL_FieldPlot for available fields
    plotAcceleration            = logical(0);
    numMin                      = 1;
    numMax                      = 999;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= '';   	%title
    absoluteValues              = logical(0);               %plots only absolute values
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	= 0;
%     graphYAxisMax             	= 5;
%     graphXAxisMin           	=   0;
%     graphXAxisMax            	=   55;
    graphLineStyle              = '-';          %'o' or '-'
   
elseif strcmp(PLOT.CustomGraphName(iPlotLocal),'Plate mobility vs. Time')
    dataNameX                   = {'Time'}; 	%e.g., 'Time', see SL_FieldPlot for available fields
    dataNameY                   = {'Plate mobility','Plate drivety'};      %see SL_FieldPlot for available fields
    plotAcceleration            = logical(0);
    numMin                      = 1;
    numMax                      = 999;
    indicateCurrentTime         = logical(1);
    DESVARIA.fieldName        	= '';   	%title
    absoluteValues              = logical(0);               %plots only absolute values
    interpolateRoughGridData    = logical(0);               %applied only for some parameters
    graphYAxisMin              	= 0;
%     graphYAxisMax             	= 2;
%     graphXAxisMin           	=   0;
%     graphXAxisMax            	=   55;
    graphLineStyle              = '-';          %'o' or '-'
    
else
    graphNameFound              = false;
    warning(['PLOT.CustomGraphName = ',PLOT.CustomGraphName{iPlotLocal},' not recognised!'])
end

%% SWITCHES CHECK
plotLegend = true;
if PLOT.loopCase>1
    plotLegend              = false;
end
if plotAcceleration && absoluteValues
    absoluteValues          = false;
end
if strcmp(graphLineStyle,'o')
    indicateCurrentTime2	= false;
else
    indicateCurrentTime2    = indicateCurrentTime;
end

%% LOOK FOR DATA FILE
if graphNameFound
    orgDir = pwd;
    fileNotFound = false;
    if ~exist([FILE.directory,'+data'],'dir')
        fileNotFound = true;
    else
        cd([FILE.directory,'+data'])
        filetofind      = [FILE.name,'_',data2import,num2str(numMax),'.mat'];
        % IF FILE NOT FOUND CHECK FOR LATEST FILE
        if ~exist(filetofind,'file')
            % check for last output file number
            for io=numMin:numMax
                filetofind	= [FILE.name,'_',data2import,num2str(io),'.mat'];
                if ispc
                    filetofind       = strrep(filetofind,'/','\');
                end
                if ~exist(filetofind,'file') %&& io~=0 %sometimes 0 has no data
                    numMax=io-1; %update max number
                    break  %exit loop
                end
            end
        end
        % ERROR MESSAGE IF NO FILE AT ALL IS NOT FOUND
        if numMax==-1 || numMax<numMin
            fileNotFound = true;
            warning(['file not found - check file directory: ',FILE.directory,'+data',filesep,filetofind])
        end
        cd(orgDir);
    end
    arraySize = abs(numMax-numMin)+1;
else
    fileNotFound = false;
end

if fileNotFound || ~graphNameFound %|| strcmp(GRID.Dim,'3-D')
    %% EMPTY PLOT IF NO FILE FOUND
    [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
    
else
    %% READ DATA
    % Read graph data produced by StagLab
    Time2check  = zeros(arraySize,1);
    xName       = cell(1,size(dataNameX,2)); xName(:) = {'NaN'};
    xGroupName 	= xName;
    xDim        = xName;
    xData       = zeros(arraySize,size(dataNameX,2));
    yName       = cell(1,size(dataNameY,2)); yName(:) = {'NaN'};
    yGroupName  = yName;
    yDim        = yName;
    yData       = zeros(arraySize,size(dataNameY,2));
    itimestep   = 0;
    for number=numMin:numMax %loop time steps
        itimestep = itimestep+1;
        % GEODYNAMIC DATA:
        % [name, dim, data]
        dummy = [FILE.directory,'+data',filesep,FILE.name,'_',data2import,num2str(number),'.mat'];
        if ispc; dummy = strrep(dummy,'/','\'); end
        load(dummy);
        
        %always read time
        Time2check(itimestep,1)  	= saveData{1,4}.*PLOT.timeConvert;
        
        %read necessary data
        try
            for idata=1:size(dataNameX,2) %read necessary data rows
                if itimestep==1 %just once
                    conversionFactorX   	= 1;
                    xName(1,idata)          = saveData(strcmp(saveData(:,1),dataNameX(1,idata)),1);
                    xGroupName(1,idata)     = saveData(strcmp(saveData(:,1),dataNameX(1,idata)),2);
                    if strcmp(xName(1,idata),'Time')
                        conversionFactorX 	= PLOT.timeConvert;      %dimensionalise time
                        xDim(1,idata)      	= {PLOT.time2plotDim};
                    else
                        xDim(1,idata)    	= saveData(strcmp(saveData(:,1),dataNameX(1,idata)),3);
                    end
                end
                dummy = saveData{strcmp(saveData(:,1),dataNameX(1,idata)),4} *conversionFactorX;
                xData(itimestep,idata)      = dummy(1,1);
                if size(dummy,1)>1
                    for idatapoints=1:size(dummy,1)
                        xDataMulti(itimestep,idata,idatapoints)   	= dummy(idatapoints,1);  %multiple data points
                    end
                end
            end
        catch errme
            warning(['The data for "',dataNameX{1,idata},'" (i.e., timestep ',num2str(number),') could not be found. You might need to recalculate and overwrite your data files!'])
            %% EMPTY PLOT IF NO FILE FOUND
            [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
            
            return
        end
        
        %read necessary data
        try
            for idata=1:size(dataNameY,2) %read necessary data rows (for plotting multiple graphs on top of each other)
                if itimestep==1 %just once
                    conversionFactorY    	= 1;
                    yName(1,idata)          = saveData(strcmp(saveData(:,1),dataNameY(1,idata)),1);
                    yGroupName(1,idata)     = saveData(strcmp(saveData(:,1),dataNameY(1,idata)),2);
                    if strcmp(yName(1,idata),'Time')
                        conversionFactorY 	= PLOT.timeConvert;     %dimensionalise time
                        yDim(1,idata)      	= {PLOT.time2plotDim};
                    else
                        yDim(1,idata)     	= saveData(strcmp(saveData(:,1),dataNameY(1,idata)),3);
                    end
                end
                dummy = saveData{strcmp(saveData(:,1),dataNameY(1,idata)),4} *conversionFactorY;
                yData(itimestep,idata)   	= dummy(1,1);  %only the data of one tracking point
                if size(dummy,1)>1
                    for idatapoints=1:size(dummy,1)
                        yDataMulti(itimestep,idata,idatapoints)   	= dummy(idatapoints,1);  %multiple data points
                    end
                end
            end
        catch errme
            warning(['The data for "',dataNameY{1,idata},'" (i.e., timestep ',num2str(number),') could not be found. You might need to recalculate and overwrite your data files!'])
            %% EMPTY PLOT IF NO FILE FOUND
            [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
            
            return
        end
    end
    
    %% DATA IMPROVEMENT
    %intelligently replace NaN's
    if exist('fit','file')
        fitFileExists   = true;
    else
        fitFileExists   = false;
        try
            fit(x,y)
        catch errme
            warning(errme.message)
        end
    end
    %smooth discrete grid data
    if interpolateRoughGridData && fitFileExists
        if strcmp(dataNameX(1,1),'Slab-tip depth') || strcmp(dataNameX(1,1),'Slab-tip horiz. position')
            %check if x-and y-data are monotonically increasing (not strictly, i.e., can be stagnant)
            if sum(all(diff(xData)>=0))<=0; xMonotonic = false; else; xMonotonic = true; end
            if sum(all(diff(yData)>=0))<=0; yMonotonic = false; else; yMonotonic = true; end
            if ~xMonotonic && ~yMonotonic %if non-monotonic on both axes
                %spline fitting over the whole graph cannot be applied!
                %do it step-wise (from point to point) instead.
                if ~all(diff(xData)) %if repeating numbers on y-axis
                    xdummy               	= zeros(size(xData));
                    xdummy(diff(xData)~=0) 	= xData(diff(xData)~=0);
                    xdummy(1)               = xData(1); %make sure first and last entries are there
                    xdummy(end)           	= xData(end); %make sure first and last entries are there
                    dummy                   = xdummy;
                    xdummy(dummy==0)       	= [];
                    xdummyFull             	= dummy;
                    ydummyFull            	= 1:size(xData,1);  ydummyFull = ydummyFull';
                    ydummy                  = ydummyFull;
                    ydummy(dummy==0)        = [];
                    % average values at repeated entries:
                    for idata=1:size(xdummyFull,1)
                        if dummy(idata,1)==0
                            %get group average
                            groupLeftIdx = idata;
                            groupRightIdx = idata; averageValue = 0;
                            ilocal = 1; groupValues = 0;
                            while groupRightIdx<size(yData,1) && dummy(groupRightIdx,1)==0
                                groupValues(ilocal,1) = yData(groupRightIdx,1);
                                ilocal = ilocal+1;
                                groupRightIdx = groupRightIdx+1;
                            end
                            if groupValues~=0
                                yData(groupLeftIdx:groupRightIdx,1) = mean(groupValues);
                            end
                        end
                    end
                end
                if ~all(diff(yData)) %if repeating numbers on y-axis
                    xdummy               	= zeros(size(yData));
                    xdummy(diff(yData)~=0) 	= yData(diff(yData)~=0);
                    xdummy(1)               = yData(1); %make sure first and last entries are there
                    xdummy(end)           	= yData(end); %make sure first and last entries are there
                    dummy                   = xdummy;
                    xdummy(dummy==0)       	= [];
                    xdummyFull             	= dummy;
                    ydummyFull            	= 1:size(yData,1);  ydummyFull = ydummyFull';
                    ydummy                  = ydummyFull;
                    ydummy(dummy==0)        = [];
                    % average values at repeated entries:
                    for idata=1:size(xdummyFull,1)
                        if dummy(idata,1)==0
                            %get group average
                            groupLeftIdx = idata;
                            groupRightIdx = idata; averageValue = 0;
                            ilocal = 1; groupValues = 0;
                            while groupRightIdx<size(xData,1) && dummy(groupRightIdx,1)==0
                                groupValues(ilocal,1) = xData(groupRightIdx,1);
                                ilocal = ilocal+1;
                                groupRightIdx = groupRightIdx+1;
                            end
                            if groupValues~=0
                                xData(groupLeftIdx:groupRightIdx,1) = mean(groupValues);
                            end
                        end
                    end
                end
                
            else
                if ~all(diff(xData)) %if repeating numbers on x-axis
                    xdummy               	= zeros(size(xData));
                    xdummy(diff(xData)~=0) 	= xData(diff(xData)~=0);
                    xdummy(1)               = xData(1); %make sure first and last entries are there
                    xdummy(end)           	= xData(end); %make sure first and last entries are there
                    dummy                   = xdummy;
                    xdummy(dummy==0)       	= [];
                    xdummyFull             	= dummy;
                    xdummyFull(dummy==0)  	= NaN;
                    ydummyFull            	= 1:size(xData,1);  ydummyFull = ydummyFull';
                    ydummy                  = ydummyFull;
                    ydummy(dummy==0)        = [];
                    
                    %remove NaN values
                    xdummy(isnan(ydummy))   = [];
                    ydummy(isnan(ydummy))   = [];
                    ydummy(isnan(xdummy))   = [];
                    xdummy(isnan(xdummy))   = [];
                    
                    % add spline to the available graph points:
                    smoothingParam      	= 0.5;         %2e-7 Smoothing of fitted line: 0: total smooth (i.e., linear), 1: no smoothing
                    fo                      = fit(ydummy(~isnan(xdummy)),xdummy(~isnan(xdummy)),'smoothingspline','smoothingparam',smoothingParam);
                    for idata=1:size(xdummyFull,1)
                        if isnan(xdummyFull(idata,1))
                            xdummyFull(idata,1) = fo(ydummyFull(idata,1));
                        end
                    end
                    xData                   = xdummyFull;
                end
                if ~all(diff(yData)) %if repeating numbers on y-axis
                    xdummy               	= zeros(size(yData));
                    xdummy(diff(yData)~=0) 	= yData(diff(yData)~=0);
                    xdummy(1)               = yData(1); %make sure first and last entries are there
                    xdummy(end)           	= yData(end); %make sure first and last entries are there
                    dummy                   = xdummy;
                    xdummy(dummy==0)       	= [];
                    xdummyFull             	= dummy;
                    xdummyFull(dummy==0)  	= NaN;
                    ydummyFull            	= 1:size(yData,1);  ydummyFull = ydummyFull';
                    ydummy                  = ydummyFull;
                    ydummy(dummy==0)        = [];
                    
                    %remove NaN values
                    xdummy(isnan(ydummy))   = [];
                    ydummy(isnan(ydummy))   = [];
                    ydummy(isnan(xdummy))   = [];
                    xdummy(isnan(xdummy))   = [];
                    
                    % add spline to the available graph points:
                    smoothingParam      	= 0.5;         %2e-7 Smoothing of fitted line: 0: total smooth (i.e., linear), 1: no smoothing
                    fo                      = fit(ydummy(~isnan(xdummy)),xdummy(~isnan(xdummy)),'smoothingspline','smoothingparam',smoothingParam);
                    for idata=1:size(xdummyFull,1)
                        if isnan(xdummyFull(idata,1))
                            xdummyFull(idata,1) = fo(ydummyFull(idata,1));
                        end
                    end
                    yData                   = xdummyFull;
                end
            end
        end
    end
    
    %% ERROR CHECKS
%     if strcmp(GRID.Type,'spherical2D')
%         % RE-FLIP DEPTH VECTOR AND ACCOUNT FOR STICKY-AIR LAYER
%         z2d = z2d +SETUP.d_air*GRID.dimFactor;  %set z=0 to real surface
%         z2d = 1*GRID.dimFactor -z2d; % reverse z-axis  (1.)
%     end
    
    %% DATA ALTERATIONS
    if absoluteValues
        yData                     	= abs(yData);
        yDataMulti                  = abs(yDataMulti);
    end
    if plotAcceleration
        dv                       	= yData(2:end,:) - yData(1:end-1,:);
        dt                        	= xData(2:end,:) - xData(1:end-1,:);
        yData                   	= dv./dt;
        if ~isempty(yData)
            yData                	= [yData;yData(1,:)*NaN]; %fill empty last row with NaNs
        end
        
        dataNameY                 	= strrep(dataNameY,'Velocity','Acceleration');
        yGroupName              	= strrep(yGroupName,'Velocity','Acceleration');
        DESVARIA.fieldName           	= strrep(DESVARIA.fieldName,'Velocities','Accelerations');
        yDim(:,:)                 	= {'cm/a^2'};
    end
    if plotSinkingOverTrenchVel
        if strcmp(dataNameX(1,1),'Slab sinking velocity')
            xData                   = xData(:,1)./xData(:,2);
            xDim                 	= {'nd'};
            xGroupName            	= {'v_{SlabSinking}/v_{Trench}'};
            xName                 	= {'v_{SlabSinking}/v_{Trench}'};
            %DESVARIA.fieldName      = 'v_{SlabSinking}/v_{Trench}';
            dataNameX             	= xGroupName;
        elseif strcmp(dataNameY(1,1),'Slab sinking velocity')
            yData                	= yData(:,1)./yData(:,2);
            yDim                   	= {'nd'};
            yGroupName            	= {'v_{SlabSinking}/v_{Trench}'};
            yName               	= {'v_{SlabSinking}/v_{Trench}'};
            DESVARIA.fieldName   	= 'v_{SlabSinking}/v_{Trench}';
            dataNameY             	= yGroupName;
        end
    end
    if normaliseYdata2median
        yData                    	= yData-median(yData(:),'omitnan');
    end
    if normaliseYdata2min
        yData                    	= yData-min(yData(:));
    end
    if normaliseYdata2max
        yData                    	= yData-max(yData(:));
    end
    
    %% DEFINE AXIS SETUP
    multipleXgraphs             = false;
    multipleYgraphs             = false;
    if size(dataNameX,2)==1 && size(dataNameY,2)>1  %single data on x-axis & multiple data on y-axis
        multipleYgraphs             = true;
        numGraphs                   = size(dataNameY,2);
        xDataPlot                   = zeros(size(yData));
        for ii=1:numGraphs; xDataPlot(:,ii) = xData; end %needs to be same size as yData
        yDataPlot                   = yData;
        DESVARIA.xlabelName      	= xName{1,1};
        DESVARIA.zlabelName      	= yGroupName{1,1};
        
    elseif size(dataNameX,2)>1 && size(dataNameY,2)==1  %multiple data on x-axis & single data on y-axis
        multipleXgraphs             = true;
        numGraphs                   = size(dataNameY,2);
        xDataPlot                   = xData;
        yDataPlot                   = zeros(size(xData));
        for ii=1:numGraphs; yDataPlot(:,ii) = yData; end %needs to be same size as xData
        DESVARIA.xlabelName       	= xGroupName{1,1};
        DESVARIA.zlabelName       	= yName{1,1};
        
    elseif size(dataNameX,2)==1 && size(dataNameY,2)==1  %single data on both axes
        numGraphs                   = 1;
        xDataPlot                   = xData;
        yDataPlot                   = yData;
        DESVARIA.xlabelName        	= xName{1,1};
        DESVARIA.zlabelName       	= yName{1,1};
        
    else
        error('Graph plot cannot have multiple data sets on both axes!')
    end
    DESVARIA.xlabelDim        	= xDim{1,1};
    DESVARIA.zlabelDim         	= yDim{1,1};

    %% GRAPH DIAGNOSTICS
    if plotYmean
        FitRangeX       = [400, 600; 700, 1000];        %x-range(s) for fit
        yMean = zeros(size(FitRangeX,1),1); yMax = yMean; yMin = yMean; yStdDev = yMean;
        if ( min(xDataPlot(:,1))<=FitRangeX(1,2) && max(xDataPlot(:,1))>=FitRangeX(2,1) ) %||... %only if range is covered
              %min(xDataPlot(:,1))<=FitRangeX(1,2)  
            for i=1:size(FitRangeX,1)
                if max(xDataPlot(:,1))<FitRangeX(2,1) %data doesn't reach 2nd range -> take last data point
                    yMax(i,1)       = yDataPlot(end,1);
                    yMin(i,1)       = yDataPlot(end,1);
                    yMean(i,1)      = yDataPlot(end,1);
                    yStdDev(i,1)	= 0;
                else
                    if logical(1) %include all points temporally after the lower limit index
                        [~,idx3]   	= min( abs(xDataPlot(xDataPlot(:,1)<=FitRangeX(i,2)) -FitRangeX(i,1)) ); %index of lower bound
                        dummy     	= yDataPlot(idx3:end,1); %values to take mean of
                        dummyX      = xDataPlot(idx3:end,1);
                        dummy(dummyX(:,1)>FitRangeX(i,2)) = NaN;
                        
                    else %points only within the limits
                        dummy     	= yDataPlot(xDataPlot(:,1)>=FitRangeX(i,1) & xDataPlot(:,1)<=FitRangeX(i,2)); %values to take mean of
                    end
                    yMax(i,1)       = max( dummy(~isnan(dummy)) );
                    yMin(i,1)       = min( dummy(~isnan(dummy)) );
                    yMean(i,1)      = mean( dummy(~isnan(dummy)) );
                    yStdDev(i,1)	= std( dummy(~isnan(dummy)) );
                end
            end
        else %if range is not covered
            yMean(:,:)      = NaN;
            yMax(:,:)       = NaN;
            yMin(:,:)       = NaN;
            yMean(:,:)      = NaN;
            yStdDev(:,:)	= NaN;
        end
    end
    
    %% SETUP FIGURE
    figure(1)
    
    %% COUNT SUBPLOTS
    PLOT.nrSubplot = PLOT.nrSubplot+1;
    
    %% SETUP SUBPLOTS
    SPOS.task   = 'createSubplot';
%     if isfield(SAVE,'PlotInPlotHax') && length(SAVE.PlotInPlotHax)>=PLOT.nrSubplot
%         SPOS.haxPlotInPlot = SAVE.PlotInPlotHax(PLOT.nrSubplot);
%     end
    [SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE);
    
    %% CURRENT AXIS
    AXcurrent = gca; %save current axes handle
    
    %% PLOTTING GRAPH
    if strcmp(dataNameY,'LLSVP Locations') %____________SPECIAL GRAPH PLOTS
        SpecialGraphPlot1        	= true;
    else
        SpecialGraphPlot1       	= false;
    end
    
    if SpecialGraphPlot1
        %data preparation
        yDataMulti(yDataMulti==0)	= NaN;   %empty data entries are set to zero, so remove them here
        nData                       = size(yDataMulti,3);
        dummy                       = zeros(size(yDataMulti,1),nData);
        for ind1=1:nData
            dummy(:,ind1)           = xDataPlot;
        end
        xDataPlot                   = dummy;
        clearvars dummy yDataPlot
        yDataPlot(:,:)              = yDataMulti(:,1,:);
        %switch adjustments
        numGraphs                   = nData;
        indicateCurrentTime2     	= indicateCurrentTime;
        colouredDots                = true;
    end
    
    %plotting
    hold on
    for iline=1:numGraphs
        if SpecialGraphPlot1 %SPECIAL GRAPH PLOTS
            legendString(1,1)       = dataNameY(1,1);
            
            GRAPH.hLine(iline,1) = plot(xDataPlot(:,iline),yDataPlot(:,iline),graphLineStyle,...
                'MarkerFaceColor',STYLE.keyLineColor,'MarkerEdgeColor','none','LineWidth',STYLE.keyLineWidth);
        
        else %STANDARD GRAPH PLOTS
            if multipleXgraphs
                legendString(iline,1)	= dataNameX(1,iline);
            elseif multipleYgraphs
                legendString(iline,1)	= dataNameY(1,iline);
            else
                legendString(iline,1)	= dataNameY(1,iline);
            end
            
            GRAPH.hLine(iline,1) = plot(xDataPlot(:,iline),yDataPlot(:,iline),graphLineStyle,...
                'MarkerFaceColor',STYLE.keyLineColor,'MarkerEdgeColor','none','LineWidth',STYLE.keyLineWidth);
        end
       
        if ( size(xDataPlot,2)==1 && size(yDataPlot,2)==1 ) || SpecialGraphPlot1 %only one line
            opColor         = STYLE.keyLineColor;
        else %multiple lines
            %Setup colour map
            CMAPPING.Task   = 'Categorical';
            CMAPPING.NumberColours = numGraphs;
            PLOT.ScientificColourMapName = STYLE.CategoricalColours;
            [SWITCH,CMAPPING] = f_DesignColourmap(SWITCH,PLOT,[],STYLE,CMAPPING);
            PLOT.ColorVector    = CMAPPING.ColourVector;
            opColor         = PLOT.ColorVector(iline,:);
        end
        opColor2            = opColor;
        if indicateCurrentTime2 && ~indicateCurrentTimeAsDot %reduce line brightness
            colourContrastFactor    = 2;
            if strcmp(STYLE.ColorMode,'light')
                opColor2          	= 1 -(1-opColor)/colourContrastFactor;
            else
                opColor2        	= opColor/colourContrastFactor;
            end
        end
        set(GRAPH.hLine(iline,1),'Color',opColor2);
        if strcmp(graphLineStyle,'o') && colouredDots
            set(GRAPH.hLine(iline,1),'MarkerFaceColor',opColor2,'MarkerEdgeColor','none');
        end
        hold on
        
        %MARKER FOR CURRENT TIME POSITION
        if indicateCurrentTime2
            %currentLineColor    = get(GRAPH.hLine(iline,1),'Color');
            dummy               = round(Time2check(:,1),4)==round(PLOT.time2plot,4);
            idxCurrent          = find(dummy==1);
            xCurrentPoint       = xDataPlot(dummy,iline);
            yCurrentPoint       = yDataPlot(dummy,iline);
            if isempty(xCurrentPoint) && strcmp(dataNameX(1,1),'Time') %if time point is outside data range
                xCurrentPoint   = PLOT.time2plot;
                yCurrentPoint   = yDataPlot(end,iline);
            elseif isempty(xCurrentPoint) && strcmp(dataNameY(1,1),'Time') %if time point is outside data range
                xCurrentPoint   = yDataPlot(end,iline);
                yCurrentPoint   = PLOT.time2plot;
            else
                %                 xCurrentPoint   = NaN; %no time graph
                %                 yCurrentPoint   = NaN;
            end
            
            if indicateCurrentTimeAsDot
                plot(xCurrentPoint,yCurrentPoint,'o',...
                    'MarkerEdgeColor','none','MarkerFaceColor',opColor);
            else
                if SpecialGraphPlot1 %SPECIAL GRAPH PLOTS
                    plot(xDataPlot(1:idxCurrent,iline),yDataPlot(1:idxCurrent,iline),graphLineStyle,...
                        'MarkerEdgeColor','none','MarkerFaceColor',opColor,'LineWidth',STYLE.keyLineWidth+0.5);
                else
                    plot(xDataPlot(1:idxCurrent,iline),yDataPlot(1:idxCurrent,iline),graphLineStyle,...
                        'Color',opColor,'LineWidth',STYLE.keyLineWidth+0.5);
                    plot(xCurrentPoint,yCurrentPoint,'o',...
                        'MarkerEdgeColor','none','MarkerFaceColor',opColor,'MarkerSize',5);
                end
            end
        end
        
        %LINE FITTINGS
        if plotYmean
            for i=1:size(yMean,1)
                hold on
                plot(FitRangeX(i,:),[yMean(i,1),yMean(i,1)],...
                    'LineWidth',STYLE.keyLineWidth,'Color',[0.5 0.5 0.5]);
                plot(FitRangeX(i,:),[yMean(i,1),yMean(i,1)],...
                    '--','LineWidth',STYLE.keyLineWidth,'Color',[0.8 0.2 0.2]);
            end
        end
    end
    
    %     xpoints = x2d_plate;
    %     upper = vh_plate;
    %     if ~strcmp(TOPO.areaFill,'none')
    %         if strcmp(TOPO.areaFill,'down')
    %             color = STYLE.keyLineColor;
    %             lower = zeros(size(upper))+(min(upper)-1/10*(max(upper)-min(upper)));
    %             if TOPO.indicateComponents && PLOT.topoLOWval<lower(1,1);
    %                 lower = zeros(size(upper))+PLOT.topoLOWval;
    %             end
    %             if SWITCH.Colorbar(1,2); lower = zeros(size(upper))+FIELD.cColorbarMIN; end %constant colorbar
    %             if SWITCH.plotDifference && strcmp(TOPO.difference,'diff'); lower = lower*0; end %in relation to 0 topography
    %
    %         elseif strcmp(TOPO.areaFill,'updown')
    %             color = STYLE.keyLineColor;
    %             lower = zeros(size(upper));
    %         else
    %             error(['TOPO.areaFill mode ',TOPO.areaFill,'not recognised!'])
    %         end
    %     end
    
    if plotLegend && ~isempty(legendString) && (size(legendString,1)>1 || ~strcmp(legendString{1,1},yGroupName{1,1}))
        plotLegend = true;
    end
    if plotLegend
        legend(GRAPH.hLine,legendString{:},'Location','SouthEast');
%         legend(GRAPH.hLine,legendString{:},'Location','NorthEast'); 
    end
    hold off
    
    DESVARIA.Task    = 'create annotation strings';
    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
    
    DESVARIA.Task    = 'make title';
    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
    
    xlabel(PLOT.xlabel)
    ylabel(PLOT.ylabel)
    
    % axes...
    xLimit              = xlim;
    yLimit              = ylim;
    if SWITCH.Colorbar(1,2)
        if ~exist('graphYAxisMin','var'); graphYAxisMin = FIELD.cColorbarMIN; end
        if ~exist('graphYAxisMax','var'); graphYAxisMax = FIELD.cColorbarMAX; end
        yLimit          = [graphYAxisMin graphYAxisMax];
        if exist('graphXAxisMin','var') && exist('graphXAxisMax','var')
            xLimit      = [graphXAxisMin graphXAxisMax];
        end
    end
    if SWITCH.AxesLimit
        if size(SWITCH.AxesLimitValues,1)>1  %more than one axis-zoom value
            xLimit      = SWITCH.AxesLimitValues(PLOT.loopCase,1:2);
        else
            xLimit      = SWITCH.AxesLimitValues(1,1:2);
        end
    end
    DESVARIA.axesLimits     = [xLimit, yLimit];
    DESVARIA.NotEqualAxes  	= true;
    DESVARIA.Task  	= 'setup axes';
    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
    
    %set(AXcurrent,'YGrid','on')
    grid on
    box on
    
    if SWITCH.Texting
        try
            text(min(xlim),max(ylim),['mean: ',num2str(yMean(1,1),3),'(',num2str(yStdDev(1,1),3),')'],'HorizontalAlignment','left','VerticalAlignment','bottom','Color',STYLE.keyColor,'FontSize',8)
            text(max(xlim),max(ylim),['mean: ',num2str(yMean(2,1),3),'(',num2str(yStdDev(2,1),3),')'],'HorizontalAlignment','right','VerticalAlignment','bottom','Color',STYLE.keyColor,'FontSize',8)
            disp(['   mean: ',num2str(yMean(1,1)),' / std. deviation: ',num2str(yStdDev(1,1))])
            disp(['   mean: ',num2str(yMean(2,1)),' / std. deviation: ',num2str(yStdDev(2,1))])
            %make data array
            if ~isfield(SAVE,'graphData')
                SAVE.graphData  = []; 
                dummy   	= 2;
            else
                dummy     	= size(SAVE.graphData,1)+2;
            end
            if dummy>11
                dummy   	= dummy-10;
                suite       = 3;
            elseif dummy>6
                dummy     	= dummy-5;
                suite       = 2;
            else
                suite       = 1;
            end
            angle           = dummy*10;
            SAVE.graphData  = [SAVE.graphData; suite, angle, yMean(1,1), yStdDev(1,1), yMean(2,1), yStdDev(2,1)];
            
            if size(SAVE.graphData,1)==15 %last one
                SAVE.Directory              = FILE.stemSave;
                SAVE.DataName               = '+graphData';
                SAVE.data                   = SAVE.graphData;
                SAVE.dat                    = logical(0);
                SAVE.txt                    = logical(0);
                SAVE.mat                    = logical(1);
                SAVE.write2existing         = logical(0);
                SAVE.overwrite              = logical(0);     %this is a global variable!
                SAVE.Pause                  = logical(0);
                
                [SAVE.overwriteAll] = f_saveData( SAVE );
            end
            
        catch 
            %nothing to do
        end
    end
    
    %% COLORBAR
    colorbar;   %define dummy colorbar
    PLOT.cb = NaN; 	%this is a flag if no colorbar is used
    
    %% SIMPLIFY FIGURE
    if SWITCH.SimplifyPlots
        %SET AXES ASPECT RATIO ----
        XaxisAR = GRID.aspectRatio(1);
        YaxisAR = 1;
        ZaxisAR = 1;
        set(gca,'PlotBoxAspectRatio',[min(PLOT.MaxAspectX,XaxisAR) YaxisAR ZaxisAR])
        %--------------------------
        
        %GATHER INFORMATION ABOUT SUBPLOT LAYOUT
        PLOT.spcurrent          = SP.current;
        PLOT.sp_chax            = AXcurrent;
        if exist('haxZ','var'); PLOT.sp_chaxz = haxZ; else; PLOT.sp_chaxz = NaN; end
        PLOT.sp_ccb             = PLOT.cb;   %colorbar handle
        PLOT.title              = PLOT.titleStringCurrent;  	%title string
        
        [PLOT] = f_DesignLayout(SWITCH,PLOT,GRID,FIELD);
    end
    
    %.........
    drawnow  %actually only needed if topo or other graph plot is included (axis change otherwise)
    %.........
    
    if isfield(PLOT,'hax'); nr_hax = size(PLOT.hax,1)+1; else; nr_hax = 1; end
    
    %% FUNCTION OUTPUT
    PLOT.hax(nr_hax,1)  = AXcurrent;
    
    %indicate axis handle as graph subplot
    SAVE.GraphPlotHandles = [SAVE.GraphPlotHandles; AXcurrent]; %indicate axis handle as graph subplot
    SAVE.PlotHandles = [SAVE.PlotHandles; AXcurrent];
    
end




