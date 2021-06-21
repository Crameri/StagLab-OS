
%%                                                          PLOT GRAPH 2.61
%
% plots individual graphs 
% (or multiple ones superposed on each other)
%
%    . calls f_DesignColourmap
%    . calls f_DesignLayout
%    . calls f_Varia
%
%                                                Fabio Crameri, 14.11.2020

function [TOPO,PLOT,SAVE] = f_plotGraph(TOPO,PLATE,FIELD,GRID,SWITCH,SETUP,PLOT,FILE,STYLE,SAVE)

%% DEFAULTS
areaFill                = 'none';
ShadedGraphFill         = true;
addWaterFill            = false;
titleStartix            = '';
titleEndix              = '';
addComponentGraphs      = false;
removeContinentalPart  	= false;
Components2Add          = {};
depthLevelofGraph       = 'surface';
    
%% SETUP GRAPH-SPECIFIC VARIABLES
if strcmp(FIELD.name,'Topography')
    Graphs                  = TOPO.topo2dp;
    legendStringGraphs      = {'Topography'};
    DESVARIA.zlabelName     = 'Topography';
    areaFill                = TOPO.areaFill;
    addWaterFill            = TOPO.water;
    addComponentGraphs      = TOPO.indicateComponents;
    Components2Add          = {'Topography', 'Iso. topography', 'Res. topography'};

elseif strcmp(FIELD.name,'Iso. topography')
    Graphs                  = TOPO.topoIso2d;
    legendStringGraphs      = {'Isostatic'};
%     DESVARIA.fieldName     = [startix,FIELD.name];
    DESVARIA.zlabelName     = 'Iso. topo';
    areaFill                = TOPO.areaFill;
    addWaterFill            = TOPO.water;
    addComponentGraphs      = TOPO.indicateComponents;
    Components2Add          = {'Topography', 'Iso. topography', 'Res. topography'};
    
elseif strcmp(FIELD.name,'Res. topography')
    Graphs                  = TOPO.topoRes2d;
    legendStringGraphs      = {'Residual'};
%     DESVARIA.fieldName     = [startix,FIELD.name];
    DESVARIA.zlabelName     = 'Res. topo';
    areaFill                = TOPO.areaFill;
    addWaterFill            = TOPO.water;
    addComponentGraphs      = TOPO.indicateComponents;
    Components2Add          = {'Topography', 'Iso. topography', 'Res. topography'};
    
elseif strcmp(FIELD.name,'Dyn. topography')
    Graphs                  = TOPO.topoDyn2d;
    legendStringGraphs      = {'Dynamic'};
%     DESVARIA.fieldName     = [startix,FIELD.name];
    DESVARIA.zlabelName     = 'Dyn. topo';
    areaFill                = TOPO.areaFill;
    addWaterFill            = TOPO.water;
    addComponentGraphs      = TOPO.indicateComponents;
    Components2Add          = {'Dyn. topography'};
    
elseif strcmp(FIELD.name,'Geoid')
    Graphs                  = FIELD.field2d ; %TOPO.geoid2dp;
    Components2Add          = {};
    if PLOT.GeoidCMB && PLOT.GeoidSurf
        legendStringGraphs 	= {'Surface', 'CMB'};
        Components2Add     	= {'Surface', 'CMB'};
    elseif PLOT.GeoidCMB
        legendStringGraphs 	= {'CMB'};
    elseif PLOT.GeoidSurf
        legendStringGraphs 	= {'Surface'};
    end
%     DESVARIA.fieldName     = [startix,FIELD.name];
    DESVARIA.zlabelName     = 'Geoid';
    DESVARIA.zlabelDim      = FIELD.dim;  %SETUP.vDim
    areaFill                = 'updown';
    addComponentGraphs      = logical(1);
    
elseif strcmp(FIELD.name,'Plate velocity')
    Graphs                  = PLATE.coreXVelocity;    %velocity in X-direction!
    Components2Add          = {};
    legendStringGraphs      = {'Plate'};
%     DESVARIA.fieldName     = [startix,FIELD.name];
    DESVARIA.zlabelName     = 'Velocity';
    DESVARIA.zlabelDim      = FIELD.dim;  %SETUP.vDim
    areaFill                = 'updown';
    addWaterFill            = false;
    
elseif strcmp(FIELD.name,'Heat flux')
    Graphs                  = FIELD.field2d;
    Components2Add          = {'Surface', 'CMB'};
    legendStringGraphs      = {'Surface', 'CMB'};
%     DESVARIA.fieldName     = [startix,FIELD.name];
    DESVARIA.zlabelName     = 'Heat flux';
    DESVARIA.zlabelDim      = FIELD.dim;  %SETUP.vDim
    areaFill                = 'updown';
    addComponentGraphs      = logical(1);
    
elseif strcmp(FIELD.name,'Surface age')
    Graphs                  = FIELD.field2d;
    Components2Add          = {};
    legendStringGraphs      = {'Surface'};
%     DESVARIA.fieldName     = [startix,FIELD.name];
    DESVARIA.zlabelName     = 'Age';
    DESVARIA.zlabelDim      = FIELD.dim;  %SETUP.vDim
    areaFill                = 'updown';
    removeContinentalPart  	= logical(0); 	%Removes continental areas from surface age data

elseif strcmp(FIELD.name,'Crustal thickness')
    Graphs                  = FIELD.field2d;
    Components2Add          = {};
    legendStringGraphs      = {'Crust'};
%     DESVARIA.fieldName     = [startix,FIELD.name];
    DESVARIA.zlabelName     = 'Thickness';
    DESVARIA.zlabelDim      = FIELD.dim;  %SETUP.vDim
    areaFill                = 'updown';

elseif strcmp(FIELD.name,'Plate-base topography')
    VARIA.Task              = 'derive plate-base topography';
    [PLOT,VARIA] = f_Varia(VARIA,FILE,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE,SETUP);
    Graphs                  = VARIA.PlateBaseTopo;
    Components2Add          = {};
    legendStringGraphs      = {'Plate Base'};
%     DESVARIA.fieldName     = [startix,FIELD.name];
    DESVARIA.zlabelName     = 'Topography';
    DESVARIA.zlabelDim      = FIELD.dim;  %SETUP.vDim
    areaFill                = 'updown';
end

%% SWITCH ADJUSTMENTS
if strcmp(GRID.Dim,'3-D')
    Components2Add     	= {};
end

%% PROBLEM CHECKS
if ~exist('Graphs','var') || isempty(Graphs) || isnan(Graphs(1,1))
    warning('The Graph cannot be plotted!')
    
    %% CREATE AN EMPTY PLOT
    [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
    return
    
end

%% SETUP COMPONENT DATA AND STRINGS
if addComponentGraphs
    for iGraphAdditions=1:size(Components2Add,2)
        if ~strcmp(FIELD.name,Components2Add{iGraphAdditions})
            if strcmp(Components2Add{iGraphAdditions},'Topography')
                Graphs              = [Graphs,TOPO.topo2dp];
                legendStringGraphs  = [legendStringGraphs,{'Actual'}];
                
            elseif strcmp(Components2Add{iGraphAdditions},'Iso. topography')
                Graphs              = [Graphs,TOPO.topoIso2d];
                legendStringGraphs  = [legendStringGraphs,{'Isostatic'}];
                
            elseif strcmp(Components2Add{iGraphAdditions},'Res. topography')
                Graphs              = [Graphs,TOPO.topoRes2d];
                legendStringGraphs  = [legendStringGraphs,{'Residual'}];
                
            elseif strcmp(Components2Add{iGraphAdditions},'Dyn. topography') && isfield(TOPO,'topoDyn2d')
                Graphs              = [Graphs,TOPO.topoDyn2d];
                legendStringGraphs  = [legendStringGraphs,{'Dynamic'}];
            end
        end
    end
end

%% SET DEFAULT VARIABLES
if ~exist('TOPO','var') || ~isfield(TOPO,'UPtiltAngleXnear')
    TOPO.UPtiltAngleXnear   = NaN;
    TOPO.UPtiltAngleXfar    = NaN;
end

%% HORIZONTAL GRID VARIABLES
%find index close to sealevel
if strcmp(depthLevelofGraph,'surface')
    idxSeaLevel     = find(abs(GRID.Z_3Dp(1,1,:))==min(abs(GRID.Z_3Dp(1,1,:))));
    idxLevel        = idxSeaLevel;
else
    error('Need to specify depth level here first!')
end
if strcmp(GRID.Type,'Cartesian')
    x_PLOT(:,:) 	= GRID.X_3Dp(:,1,idxLevel);
    y_PLOT(:,:) 	= GRID.Y_3Dp(1,:,idxLevel);
    if isfield(PLATE,'TrenchDepthX')
        UPtiltXnear    	= TOPO.UPtiltAngleXnear;
        UPtiltXfar    	= TOPO.UPtiltAngleXfar;
        TopoTrenchX  	= PLATE.TrenchDepthX;
        TopoBABX      	= PLATE.BackArcBasinX;
        TopoBABextentX  = PLATE.BackArcBasinExtentX;
        TopoInundationX = PLATE.InundationX;
        TopoIArcX       = PLATE.IslandArcX;
        TopoFBulgeX     = PLATE.ForeBulgeX;
    end
elseif strcmp(GRID.Type,'spherical2D')
    dummy = zeros(size(GRID.X_3Dsp));
    try
        dummy           = GRID.X_3Dsp -min(GRID.X_3Dsp(:,1,:)); %remove min value of each depth level
    catch %the above seems to cause problems on some machines (Windows)
        for izz=1:size(GRID.X_3Dsp,3) %loop through depth levels
            dummy(:,1,izz)	= GRID.X_3Dsp(:,1,izz) -min(GRID.X_3Dsp(:,1,izz)); %remove min value of each depth level
        end
    end
    x_PLOT(:,:) 	= dummy(:,1,idxLevel);
    y_PLOT(:,:) 	= GRID.Y_3Dsp(1,:,idxLevel);
    if isfield(PLATE,'TrenchDepthX')
        UPtiltXnear    	= TOPO.UPtiltAngleXnear -min(GRID.X_3Dsp(:,1,idxSeaLevel)); %remove min value of each depth level
        UPtiltXfar    	= TOPO.UPtiltAngleXfar -min(GRID.X_3Dsp(:,1,idxSeaLevel)); %remove min value of each depth level
        TopoTrenchX  	= PLATE.TrenchDepthX -min(GRID.X_3Dsp(:,1,idxSeaLevel));
        TopoBABX      	= PLATE.BackArcBasinX -min(GRID.X_3Dsp(:,1,idxSeaLevel));
        TopoBABextentX	= PLATE.BackArcBasinExtentX -min(GRID.X_3Dsp(:,1,idxSeaLevel));
        TopoInundationX = PLATE.InundationX -min(GRID.X_3Dsp(:,1,idxSeaLevel));
        TopoIArcX       = PLATE.IslandArcX -min(GRID.X_3Dsp(:,1,idxSeaLevel));
        TopoFBulgeX     = PLATE.ForeBulgeX -min(GRID.X_3Dsp(:,1,idxSeaLevel));
    end
else
    error('Needs to be adjusted here for other geometries!')
end

%% DATA ADJUSTMENTS
% REMOVE CONTINENTAL PART IF NEEDED
if removeContinentalPart && strcmp(FIELD.name,'Surface age')
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
        if exist('CClocation','var')
            error('Surface age splitting into oceanic areas only is not implemented yet for yinyang geometry!')
        end
    else %all other grid types
        GraphsCont 	= zeros(size(Graphs))*NaN;
        GraphsOcean	= GraphsCont;
        if exist('CClocation','var')
            GraphsCont(CClocation==1)	= Graphs(CClocation==1); %this leaves NaNs in unwanted areas
            GraphsOcean(CClocation==0) 	= Graphs(CClocation==0);
        end
    end
    
    %Update field to process with wanted array
    Graphs      = GraphsOcean;
end

%% SUBPLOT LAYOUT...
PLOT.nrSubplot = PLOT.nrSubplot+1; %count subplots

%subplot layout...
SPOS.task   = 'createSubplot';
if isfield(SAVE,'PlotInPlotHax') && length(SAVE.PlotInPlotHax)>=PLOT.nrSubplot
    SPOS.haxPlotInPlot = SAVE.PlotInPlotHax(PLOT.nrSubplot);
end
[SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE);

%% CURRENT AXIS
AXcurrent = gca; %save current axes handle

%% ANNOTATION STRINGS
if SWITCH.plotGraphVsTime
    DESVARIA.xlabelName = 'Time';
end
if strcmp(GRID.Dim,'3-D')
    DESVARIA    = rmfield(DESVARIA,'zlabelName');
    if isfield(DESVARIA,'zlabelDim')
        DESVARIA 	= rmfield(DESVARIA,'zlabelDim');
    end
end
DESVARIA.Task 	= 'create annotation strings';
[PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
PLOT.titleStringCurrent = [titleStartix,PLOT.titleStringCurrent,titleEndix]; %adjust title


if ~isfield(SAVE,'AllFilesTitleString') %is needed for legend in plot in plot mode (might be improved)
    SAVE.AllFilesTitleString    = PLOT.titleStringCurrent;
else
    SAVE.AllFilesTitleString    = [SAVE.AllFilesTitleString, {PLOT.titleStringCurrent}];
end
if ~isfield(SAVE,'AllFilesLegendString') %is needed for legend in plot in plot mode (might be improved)
    SAVE.AllFilesLegendString   = PLOT.legendStringCurrent;
else
    SAVE.AllFilesLegendString   = [SAVE.AllFilesLegendString, {PLOT.legendStringCurrent}];
end
    
if strcmp(GRID.Dim,'2-D')
    numGraphs   = size(Graphs,2);
elseif strcmp(GRID.Dim,'3-D')
    numGraphs   = size(Graphs,3);
end

%% FIGURE...
figure(1)

%% LEGEND ADJUSTMENT
if ~verLessThan('matlab','9.1.0') %MATLAB 2016b and later
    set(gcf,'defaultLegendAutoUpdate','off') %prevent autoupdating legends with added white grid lines
end

%% LOOP GRAPHS TO PLOT
for iGraph=numGraphs:-1:1
    if strcmp(GRID.Dim,'2-D')
        %% GRAPH
        graph2plot              = Graphs(:,iGraph);
        
        %% STYLE SETUP
        if iGraph==1
            color2plot          = STYLE.keyLineColor;
            areaOpacity         = 0.2;
            lineWidth           = STYLE.keyLineWidth;
            lineStyle           = '-';
            if TOPO.indicateComponents
                lineWidth      	= STYLE.keyLineWidth+0.25;
                areaOpacity   	= 0.1;
            end
        else
            if iGraph==2
                color2plot   	= [0.5 0.5 0.5];
                lineStyle      	= '--';
            elseif iGraph==3
                color2plot  	= [0.5 0.5 0.5];
                lineStyle     	= ':';
            elseif iGraph==4
                color2plot   	= [0.5 0.5 0.5];
                lineStyle    	= '-.';
            else
                color2plot   	= [0.4 0.5 0.6];
                lineStyle     	= '-';
            end
            areaOpacity         = 0.0;
            lineWidth           = max(0.25,STYLE.keyLineWidth-0.25);
            addWaterFill        = false;
        end
        if SWITCH.PlotInPlot
            color2plot = PLOT.gradColor4File;
        end
        
        %% SAVING TOPOGRAPHY
        if TOPO.saveData && strcmp(FIELD.name,'Topography')
            SAVE.Directory              = FILE.directory;
            SAVE.DataName               = [FILE.name,'_topo',num2str(FILE.number)];
            SAVE.data                   = [x_PLOT, graph2plot]; %[x, topo]
%             SAVE.mat                    = logical(1);
%             SAVE.dat                    = logical(0);
%             SAVE.txt                    = logical(0);
%             if SAVE.dat
%                 SAVE.mat                    = logical(0);
%                 SAVE.dat                    = logical(1);
%             end
            SAVE.write2existing         = logical(0);
            [SAVE.overwriteAll] = f_saveData( SAVE );
        end
        
        %% PLOTTING
        if addWaterFill && ~strcmp(areaFill,'updown') %add water layer
            colorWater = [0.4 0.6 0.7]; alphaWater = 0.2;
            if strcmp(STYLE.ColorMode,'dark')
                colorWater = [0.4 0.6 0.7]; alphaWater = 0.7;
            end
            graphNegative = graph2plot';
            graphNegative(graph2plot>=0) = NaN; %remove positive topo values
            %separate the two lines
            NanLocations = find(isnan(graphNegative));
            NanLocations = [0 NanLocations size(graphNegative,2)+1]; %add start and end point
            nLineSegment = size(NanLocations,2)-1;
            for iLineSegment=1:nLineSegment
                startLoc = NanLocations(1,iLineSegment) +1;
                endLoc = NanLocations(1,iLineSegment+1) -1;
                if endLoc<=startLoc+1; continue; end %skip this one: no line in between
                
                S_xpoints = x_PLOT(startLoc:endLoc,1)';
                S_upper = graphNegative(1,startLoc:endLoc);
                S_lower = zeros(size(S_upper));
                
                fill([min(S_xpoints); S_xpoints'; max(S_xpoints)],[min(S_lower); S_upper'; min(S_lower)],...
                    colorWater,'EdgeColor','none','FaceAlpha',alphaWater);
                hold on
            end
        end
        
        if ~SWITCH.plotGraphVsTime %PLOT TOPOGRAPHY GRAPH
            PLOT.graphHandle = plot(x_PLOT,graph2plot,'Color',color2plot,'LineWidth',lineWidth,'LineStyle',lineStyle);
            hold on
            
            graphHandlesAll(1,iGraph)     = PLOT.graphHandle;
            if numGraphs>1 && iGraph==1 && ~strcmp(SAVE.legendPlotGraphExists,FIELD.name) %if multiple graphs and they are repeated in one figure
                if strcmp(FIELD.name,'Topography') %don't list actual topography in legend
                    legend(graphHandlesAll(2:end),legendStringGraphs(2:end),'FontSize',STYLE.keyFontSize-2,'Location','SouthEast','Interpreter','none') %area
                else
                    legend(graphHandlesAll,legendStringGraphs,'FontSize',STYLE.keyFontSize-2,'Location','SouthEast','Interpreter','none') %area
                end
                SAVE.legendPlotGraphExists = FIELD.name;
            end
        
        elseif SWITCH.plotGraphVsTime %remember x and topo for graph vs. time plot
            if SAVE.LastFile; SAVE.graphX = x_PLOT'; end
            if ~isfield(SAVE,'graphTopo')
                SAVE.graphTopo = graph2plot';
                SAVE.graphTime	= PLOT.time;
            else
                SAVE.graphTopo  = [SAVE.graphTopo; graph2plot'];
                SAVE.graphTime	= [SAVE.graphTime; PLOT.time]; %non-dimensional time
            end
            if SAVE.LastFile
                SAVE.graphTime = SAVE.graphTime.*PLOT.timeConvert; %make dimensional according to last automatic time dimension
                contourf(SAVE.graphX,SAVE.graphTime,SAVE.graphTopo,PLOT.numContours,'LineColor','none');
                PLOT.cb=colorbar;
                if SWITCH.Colorbar(1,2); set(gca,'clim',[-max(abs(SAVE.graphTopo(:))),max(abs(SAVE.graphTopo(:)))]); end %constant cb
                if SWITCH.DimensionalMode || SWITCH.DimensionalInput; title(PLOT.cb,[FIELD.name,' [',FIELD.dim,']']); else; title(PLOT.cb,[FIELD.name,' [',FIELD.dim,']']); end
            end
        end
        
        if ~strcmp(areaFill,'none')
            xpoints = x_PLOT;
            upper   = graph2plot;
            if strcmp(areaFill,'down')
                lower = zeros(size(upper))+(min(upper)-1/10*(max(upper)-min(upper)));
                if TOPO.indicateComponents; PLOT.topoLOWval = lower(1,1); lower = zeros(size(upper))+min(Graphs(:))+1/10*min(Graphs(:)); end %set very low to be deeper than dynamic topo low
                if SWITCH.Colorbar(1,2); lower = zeros(size(upper))+FIELD.cColorbarMIN; end %constant axis min/max
                if SWITCH.plotDifference && strcmp(TOPO.difference,'diff'); lower = lower*0; end %in relation to 0 topography
                
            elseif strcmp(areaFill,'updown')
                lower = zeros(size(upper));
            else
                error(['areaFill mode ',areaFill,'not recognised!'])
            end
            
            if ShadedGraphFill && areaOpacity~=0
                shadedGraph     = true;
            else
                shadedGraph     = false;
            end
            
            if ~shadedGraph || (shadedGraph && ~strcmp(areaFill,'updown'))
                %             PLOT.TOPOhArea(iGraph,1) = fill([min(xpoints); xpoints; max(xpoints)],[min(lower); upper; min(lower)],color,'EdgeColor','none','FaceAlpha',areaOpacity);
                PLOT.TOPOhArea(iGraph,1) = fill([xpoints(1,1); xpoints; xpoints(end,1)],[min(lower); upper; min(lower)],color2plot,'EdgeColor','none','FaceAlpha',areaOpacity);
            end
            
            if ShadedGraphFill && areaOpacity~=0
                if areaOpacity~=0; areaOpacity = 0.2; end                
                if strcmp(areaFill,'updown')
                    DESVARIA.graph2plot     = graph2plot';
                    DESVARIA.x2plot         = x_PLOT';
                    DESVARIA.areaOpacity    = areaOpacity;
                    DESVARIA.color2plot     = color2plot;
                    DESVARIA.Task           = 'make shaded area';
                    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
                else
                    cdata           = [min(lower); upper; min(lower)];
                    cdata           = (cdata-min(cdata))/(max(cdata)-min(cdata)); %normalise
                    if strcmp(STYLE.ColorMode,'light')
                        cdata       = 1-cdata;         %flip colors
                    end
                    set(PLOT.TOPOhArea(iGraph,1),'CData',cdata,'FaceColor','interp','FaceAlpha',areaOpacity)
                    colormap(gca,'gray')
                end
            end
        end
        
        %% ANNOTATIONS
        if iGraph==numGraphs
            xlabel(PLOT.xlabel);
            ylabel(PLOT.ylabel);
            
            DESVARIA.Task    = 'make title';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            % axes
            DESVARIA.NotEqualAxes  	= true;
            DESVARIA.Task  	= 'setup axes';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            
            if SWITCH.Colorbar(1,2) && ~SWITCH.plotGraphVsTime; ylim([FIELD.cColorbarMIN FIELD.cColorbarMAX]); end
            if SWITCH.AxesLimit
                if size(SWITCH.AxesLimitValues,1)>1  %more than one axis-zoom value
                    xlim(SWITCH.AxesLimitValues(PLOT.loopCase,1:2));
                else
                    xlim(SWITCH.AxesLimitValues(1,1:2));
                end
            end
            %                 set(AXcurrent,'YGrid','on') %only horizontal grid
            if ~SWITCH.plotGraphVsTime; grid on; end
            box on
            
            if SWITCH.Texting
                text(min(xlim),max(ylim),['min val: ',num2str(min(graph2plot(:)),2)],'HorizontalAlignment','left','VerticalAlignment','bottom','Color',STYLE.keyColor,'FontName',STYLE.AllFontName,'FontSize',8)
                text(max(xlim),max(ylim),['max val: ',num2str(max(graph2plot(:)),2)],'HorizontalAlignment','right','VerticalAlignment','bottom','Color',STYLE.keyColor,'FontName',STYLE.AllFontName,'FontSize',8)
                disp(['   max.graph: ',num2str(max(graph2plot)),' / min.graph: ',num2str(min(graph2plot))])
            end
            
            if SWITCH.plotDifference && iGraph==1 %only for difference plot
                if isfield(PLOT,'titleString')
                    legend(PLOT.TOPOhArea,PLOT.titleString{1,:},'Location','southeast','Interpreter','none') %area
                else
                    legend(PLOT.TOPOhArea,PLOT.fileName,'Location','southeast','Interpreter','none') %area
                end
                
            elseif SWITCH.PlotInPlot && ~SWITCH.plotGraphVsTime %only for plot in plot
                if ~isfield(SAVE,'PlotInPlotHax') || length(SAVE.PlotInPlotHax)<PLOT.nrSubplot; SAVE.PlotInPlotHax(PLOT.nrSubplot) = gca; end %needed for next file plot to find the subplot
                if ~isfield(SAVE,'AllFilesTOPOhLine')
                    SAVE.AllFilesTOPOhLine = PLOT.graphHandle; %Remember topo handles for all individual files
                else
                    SAVE.AllFilesTOPOhLine = [SAVE.AllFilesTOPOhLine, PLOT.graphHandle]; %Remember topo handles for all individual files
                end
                if SAVE.LastFile
                    if isfield(PLOT,'legendString') %only for last one
                        idxLow      = (SP.current(1,3)-1)*PLOT.nrTimeSnapshots+1;
                        idxHigh     = SP.current(1,3)*PLOT.nrTimeSnapshots;
                        legend(SAVE.AllFilesTOPOhLine,PLOT.legendString{1,idxLow:idxHigh},'Location','southeast')
                    else
                        if isfield(SAVE,'FirstSuiteLegendString') && size(SAVE.FirstSuiteLegendString,2)>4
                            legString = [SAVE.FirstSuiteLegendString, SAVE.AllFilesLegendString(1,end)];
                            legHandle = [SAVE.FirstSuiteTOPOhLine(1,1), SAVE.AllFilesTOPOhLine(1,end)];
                        else
                            legString = SAVE.AllFilesLegendString; legHandle = SAVE.AllFilesTOPOhLine;
                        end
                        warning('off','MATLAB:legend:IgnoringExtraEntries') %suppress warning
                        legend(legHandle,legString,'Location','southeast')
                        warning('on','MATLAB:legend:IgnoringExtraEntries')
                    end
                end
            end
            if strcmp(FIELD.name,'Topography') && SWITCH.PlateDiagnostics
                %switches
                indicateUPtilt              = false;
                indicateTrench              = false;
                indicateForeBulge           = false;
                indicateIslandArc           = false;
                indicateBackArcBasin        = false;
                indicateBackArcBasinExtent  = false;
                indicateInundation          = false;
                if strcmp(FIELD.name,'Topography')
                    if PLOT.indicateUPtilt;                 indicateUPtilt = true;              end
                    if PLOT.indicateTrenchDepth;            indicateTrench = true;              end
                    if PLOT.indicateForeBulgeHeight;    	indicateForeBulge = true;           end
                    if PLOT.indicateIslandArcHeight;    	indicateIslandArc = true;           end
                    if PLOT.indicateBackArcBasinDepth;  	indicateBackArcBasin = true;        end
                    if PLOT.indicateBackArcBasinExtent;  	indicateBackArcBasinExtent = true; 	end
                    if PLOT.indicateInundation;             indicateInundation = true;          end
                end
                %plot additions
                if indicateUPtilt && isfield(TOPO,'UPtiltAngleXnear') && ~isnan(UPtiltXnear)
                    hold on
                    plot([UPtiltXnear,UPtiltXfar],[TOPO.UPtiltAngleZnear,TOPO.UPtiltAngleZfar],'-','Color',[0.75 0.1 0.1],'LineWidth',1.5)
                end
                if indicateTrench && ~isnan(PLATE.TrenchDepth)
                    hold on
                    plot(TopoTrenchX,PLATE.TrenchDepth,'or')
                end
                if indicateForeBulge && ~isnan(PLATE.ForeBulgeZ)
                    hold on
                    plot(TopoFBulgeX,PLATE.ForeBulgeZ,'or')
                end
                if indicateIslandArc && ~isnan(PLATE.IslandArcZ)
                    hold on
                    plot(TopoIArcX,PLATE.IslandArcZ,'or')
                end
                if indicateBackArcBasin && ~isnan(PLATE.BackArcBasinZ)
                    hold on
                    plot(TopoBABX,PLATE.BackArcBasinZ,'or')
                end
                if indicateBackArcBasinExtent && ~isnan(PLATE.BackArcBasinExtentZ)
                    hold on
                    plot(TopoBABextentX,PLATE.BackArcBasinExtentZ,'or')
                end
                if indicateInundation && ~isnan(PLATE.InundationZ)
                    hold on
                    plot(TopoInundationX,PLATE.InundationZ,'o','MarkerEdgeColor',[0 0.55 0.75],'MarkerFaceColor',[0 0.55 0.75]) %dark cyan
                end
            end
        end
        
    elseif strcmp(GRID.Dim,'3-D')
        %% 3-D
        
        %% GRAPH
        graph2plot              = Graphs(:,:,iGraph);
        if SWITCH.plot3Dtopo2D %2-D plot
            contourf(x_PLOT,y_PLOT,graph2plot',PLOT.numContours,'EdgeColor','none')
            axis tight
            axis equal
            
        else %3-D topography plot
            surf(y_PLOT,x_PLOT,graph2plot,'EdgeColor','none')
            % contour(xi,yi,zi)
            set(gca,'ztick',[]); %remove z ticks
            axis tight
            axis equal %they need to be in front of campos commands
            if TOPO.exagFactor>1
                set(gca,'DataAspectRatio',[1 1 1/TOPO.exagFactor])
            end
            %set camera position
            campos(PLOT.cameraPosition.*[1 1 -1/TOPO.exagFactor]) %z-axis is not flipped
            %set perspective mode
            camproj(PLOT.camProjection)
        end
        box on
        DESVARIA.Task    = 'make title';
        [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
        
        xlabel(PLOT.xlabel);
        ylabel(PLOT.ylabel);
        
        %INPUT here:
        plotTopoProfiles = logical(0);
        if strcmp(FILE.name,'fs6_30') %also used later on
            xprofLoc = []; yprofLoc = [140]; %[50, 206] %in [#gridpoints]
        else
            xprofLoc = []; yprofLoc = [140, 20]; %[50, 206] %in [#gridpoints]
        end
        %-----------
        if plotTopoProfiles
            colors4Profiles = [[0 0 0]; [0.5 0.5 0.5]; [0.3 0.3 1.0]; [1 0.3 0.3]; [0.6 0.6 0.3]; [0.3 1 0.3]; [0.5 0.5 0.5]];
            % add profile lines to previous plots
            hold on
            if ~isempty(xprofLoc)
                numProfx = size(xprofLoc,2);
                profx = y_PLOT';
                for iprof=1:numProfx
                    profx(1+iprof,:) = graph2plot(xprofLoc(1,iprof),:);
                    if SWITCH.plot3Dtopo2D %2-D plot
                        plot(y_PLOT*0+x_PLOT(xprofLoc(1,iprof),1),y_PLOT,'Color',colors4Profiles(iprof,:))
                    else
                        plot3(y_PLOT,y_PLOT.*0+x_PLOT(xprofLoc(1,iprof),1),graph2plot(xprofLoc(1,iprof),:),'Color',colors4Profiles(iprof,:))
                    end
                end
            end
            if ~isempty(yprofLoc)
                numProfy = size(yprofLoc,2);
                profy = x_PLOT';
                for iprof=1:numProfy
                    profy(1+iprof,:) = graph2plot(:,yprofLoc(1,iprof));
                    if SWITCH.plot3Dtopo2D %2-D plot
                        plot(x_PLOT,x_PLOT*0+y_PLOT(1,yprofLoc(1,iprof)),'Color',colors4Profiles(iprof,:))
                    else
                        plot3(x_PLOT.*0+y_PLOT(1,yprofLoc(1,iprof)),x_PLOT,graph2plot(:,yprofLoc(1,iprof)),'Color',colors4Profiles(iprof,:))
                    end
                end
            end
            %new figure with actual profiles
            figure(2)
            SwitchPlotSeparately = logical(1); %plot separately for each file
            if SwitchPlotSeparately; clf; end
            set(gcf,'Position',[47 1 685 300]);
            if exist('profx','var')
                if exist('profy','var'); subplot(1,2,1); end
                for iprof=2:size(profx,1)
                    hpr(iprof-1) = plot(profx(1,:),profx(iprof,:),'-','Color',colors4Profiles(iprof-1,:));
                    if iprof-1==1; dummy=1.25; else; dummy=0.75; end
                    set(hpr(iprof-1),'LineWidth',dummy)
                    hold on
                end
            end
            if exist('profy','var')
                if exist('profx','var'); subplot(1,2,2); end
                for iprof=2:size(profy,1)
                    hpr(iprof-1) = plot(profy(1,:),profy(iprof,:),'-','Color',colors4Profiles(iprof-1,:));
                    if iprof-1==1; dummy=1.25; else; dummy=0.75; end
                    set(hpr(iprof-1),'LineWidth',dummy)
                    hold on
                end
                % %                         to save profile data:
                %                         TOPO1.x = profy(1,:)
                %                         TOPO1.y = profy(iprof,:)
                %                         TOPO1.yLoc = y_PLOT(yprofLoc(1))
                %                         save('~/Desktop/TOPO1.mat','TOPO1')
            end
            if SwitchPlotSeparately || SAVE.LastFile
                axis tight
                ylim([-7 4])
                grid on
                box on
                xlabel(PLOT.xlabel)
                ylabel(PLOT.zlabel)
                %legend(hpr,'thin UP-plate part','thick UP-plate part')
                if strcmp(FILE.name,'fs6_30')
                    %legend(hpr,['y = ',num2str(y_PLOT(yprofLoc(1)),4),' km'])
                    legend(hpr,['3-D, across slab (y = ',num2str(y_PLOT(yprofLoc(1)),4),' km)'])
                else
                    %legend(hpr,['y = ',num2str(y_PLOT(yprofLoc(1)),4),' km'],['y = ',num2str(y_PLOT(yprofLoc(2)),4),' km'])
                    legend(hpr,['3-D, across slab gap (y = ',num2str(y_PLOT(yprofLoc(1)),4),' km)'],['3-D, across slab (y = ',num2str(y_PLOT(yprofLoc(2)),3),' km)'])
                end
                %pause(0.1) %let the figure time to adjust
                f_DesignFigure(STYLE);
                if logical(1) %SWITCH.fBackground
                    BACK.BulletString 	= STYLE.SCHAR.hugeBulletDark;
                    BACK.ApplyTo        = 'current';          %'current', 'all'
                    f_DesignBackground(BACK) %set plot background of each subplot
                end
                %SAVE FIGURE
                if SAVE.Figure; SAVE.FigureNr=2; f_saveFigure(SAVE); end
            end
            figure(1) %back to fig. 1
        end
        DESVARIA.Task    = 'make colorbar';
        [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
        
        clearvars low high tix tiy
    end %3-D
end



% %% FIGURE SIZE & POSITION
% figure(1)
% hold on

%% COLORBAR
if strcmp(GRID.Dim,'3-D') || SWITCH.plotGraphVsTime && SAVE.LastFile
    %% SPECIFY COLORMAP FOR CURRENT AXES
    AXhandles(1) = AXcurrent;  %big plot axis handle
    if exist('haxZ','var')
        AXhandles(2) = haxZ; %zoom plot axis handle
        numLoopsZoom = 2;
    else
        numLoopsZoom = 1;
    end
    if ~isfield(PLOT,'cb'); PLOT.cb = colorbar; end
    for iZoom=1:numLoopsZoom
        if ~SWITCH.multipleColormaps && iZoom>1; continue; end
        PLOT.hax = AXhandles(iZoom);
        [SWITCH,~] = f_DesignColourmap(SWITCH,PLOT,FIELD,[],[]);
    end
else
    colorbar;   %define dummy colorbar
    PLOT.cb = NaN; 	%this is a flag if no colorbar is used
end

%% SIMPLIFY FIGURE
if SWITCH.SimplifyPlots
    if strcmp(GRID.Dim,'2-D')
        %SET AXES ASPECT RATIO
        XaxisAR     = GRID.aspectRatio(1);
        XaxisARmax  = PLOT.MaxAspectX; if strcmp(GRID.Type,'Cartesian'); XaxisARmax = 999; end
        YaxisAR     = GRID.aspectRatio(2);
        ZaxisAR     = 1;
        set(gca,'PlotBoxAspectRatio',[min(XaxisARmax,XaxisAR) YaxisAR ZaxisAR])
    end
    
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

%% FUNCTION OUTPUT
if isfield(PLOT,'hax'); nr_hax = size(PLOT.hax,1)+1; else; nr_hax = 1; end
PLOT.hax(nr_hax,1)  = AXcurrent;
%indicate axis handle as graph subplot
SAVE.GraphPlotHandles   = [SAVE.GraphPlotHandles; AXcurrent];
SAVE.PlotHandles        = [SAVE.PlotHandles; AXcurrent];


