
%%                                                 PLOT SPECIAL FIELDS 2.0
%       . includes Plot Grid
%       . includes Plate Sketch (i.e., Tectonics)
%       . includes Parameter Table
%       . includes Surface Field Variation
%       . includes Plot Tracers
% 
%       . calls f_DiagnosticsSurface
%                                                Fabio Crameri, 08.08.2018

function [PLOT,SAVE] = f_plotFieldSpecial(iPlotLocal,FIELD_M,FIELD,FILE,GRID,SETUP,SWITCH,PLOT,STYLE,SAVE,PLATE)


%%                                                           PLOT GRID 1.1

%                                                Fabio Crameri, 12.07.2016
%% NOTES
% not implemented for 3-D yet...........

if strcmp(FIELD.name,'Grid')
    if strcmp(GRID.Dim,'3-D') % 3-D
        if strcmp(GRID.Dim,'3-D')
            warning('Grid plot not implemented for 3-D geometry!')
        end
        %% CREATE AN EMPTY PLOT
        [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
        return
        
    end
    
    %% COUNT SUBPLOTS
    PLOT.nrSubplot     = PLOT.nrSubplot+1;
    
    %% SETUP PLOT
    % figure...
    figure(1)
    
    %subplot layout...
    SPOS.task   = 'createSubplot';
    [SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE);
    
    % CURRENT AXIS
    AXcurrent = gca; %save current axes handle
    
    %% ADJUSTING DATA
    if PLOT.nCoarserGrid==1
        if strcmp(GRID.Type,'spherical2D')
            xPoints     = GRID.x2ds;
            zPoints     = GRID.z2ds;
        else
            xPoints     = GRID.x2dp;
            zPoints     = GRID.z2dp;
        end
    elseif PLOT.nCoarserGrid>1 %actual resolution
        DESVARIA.fieldName = [FIELD.name,' (',num2str(PLOT.nCoarserGrid),'\timescoarser)'];
        if strcmp(GRID.Type,'spherical2D')
            xPoints     = GRID.x2ds(1:PLOT.nCoarserGrid:end,1:PLOT.nCoarserGrid:end,1:PLOT.nCoarserGrid:end);
            zPoints     = GRID.z2ds(1:PLOT.nCoarserGrid:end,1:PLOT.nCoarserGrid:end,1:PLOT.nCoarserGrid:end);
        else
            xPoints     = GRID.x2dp(1:PLOT.nCoarserGrid:end,1:PLOT.nCoarserGrid:end,1:PLOT.nCoarserGrid:end);
            zPoints     = GRID.z2dp(1:PLOT.nCoarserGrid:end,1:PLOT.nCoarserGrid:end,1:PLOT.nCoarserGrid:end);
        end
    else
        error(['PLOT.nCoarserGrid has to be >= 0 and not ',num2str(PLOT.nCoarserGrid)])
    end
    
    %% PLOTTING
    plot(xPoints,zPoints,xPoints',zPoints','Color',1-STYLE.BlackOrWhiteColor,'LineWidth',0.25)
    
    %% TITLE STRING & TIME
    DESVARIA.Task    = 'create annotation strings';
    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
    
    DESVARIA.Task    = 'make title';
    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
    
    if ~isfield(SAVE,'AllFilesTitleString')
        SAVE.AllFilesTitleString = PLOT.titleStringCurrent;
    else
        SAVE.AllFilesTitleString = [SAVE.AllFilesTitleString, {PLOT.titleStringCurrent}];
    end
    
    %% ANNOTATION
    xlabel(PLOT.xlabel)
    ylabel(PLOT.ylabel)
    if strcmp(GRID.Type,'spherical2D') && SWITCH.spherical2DCenterText
        text(0.0,0.0,[num2str(GRID.nx),'\times',num2str(GRID.nz)],'Color',STYLE.keyColor,'FontName',STYLE.AllFontName,...
            'HorizontalAlignment','center')
    end
    
    %% AXES SETUP
    %limit axes
    if strcmp(GRID.Type,'spherical2D')
        BIndicationHeight       = 0; %adjust for plate boundary indications
        if PLOT.indicateTrench || PLOT.indicateRidge
            BIndicationHeight   = abs(max(GRID.Z_3Dp(:))-min(GRID.Z_3Dp(:))) /10;
        end
        DESVARIA.axesLimits     = [min(GRID.x2ds(:))-BIndicationHeight,max(GRID.x2ds(:))+BIndicationHeight,min(GRID.z2ds(:))-BIndicationHeight,max(GRID.z2ds(:))+BIndicationHeight];
    end
    if ~strcmp(GRID.Type,'yinyang')
        if SWITCH.AxesLimit
            if size(SWITCH.AxesLimitValues,1)>1  %individual axes values
                DESVARIA.axesLimits     = SWITCH.AxesLimitValues(PLOT.loopCase,:);
            else %all plots the same
                DESVARIA.axesLimits     = SWITCH.AxesLimitValues;
            end
        end
    end
    if strcmp(GRID.Type,'Cartesian'); DESVARIA.flipAxes = true; end
    DESVARIA.axisColour = 'none';
    DESVARIA.Task    = 'setup axes';
    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);

    PLOT.cb = NaN; 	%this is a flag if no colorbar is used
end
    



%%                               PLATE SKETCH WITH TECTONIC PARAMETERS 2.8

%                                                Fabio Crameri, 08.08.2018
%% NOTES
% not implemented for 3-D yet...........

if strcmp(FIELD.name,'Tectonics')
    %% INPUT
    indDimAtBot             = logical(1);
    fontSize                = 11;
    fontName                = STYLE.AllFontName;
    backGcolor              = 'none';
    backGcolorPlate         = 'none';
    fontColorPlate          = [1 1 1];
    % arrowL                  = '<<';
    % arrowL                  = [char(9664),char(9664)];
    % arrowL                  = [char(10094),char(10094)];
    arrowL                  = [char(10216),char(10216)];
    % arrowR                  = '>>';
    % arrowR                  = [char(9654),char(9654)];
    % arrowR                  = [char(10095),char(10095)];
    arrowR                  = [char(10217),char(10217)];
    
    showAbsoluteNumbers     = logical(1);
    showSpreadingData     	= logical(0);
    % PLOT.indicateTrench
    % PLOT.indicateRidge
    
    %% DESIGN INPUT
    plateColor              = [0.5 0.5 0.5];
    slabColor               = plateColor;
    subfaultColor           = [1 1 1];
    foreGcolor              = [0 0 0];
    spreadColor             = subfaultColor;
    precisionNum            = 2;
    BmarkerColor            = [0 0 0];
    crossColor              = [0.9 0.9 0.9];
    
    %% COLOUR ADJUSTMENTS
    if strcmp(STYLE.ColorMode,'dark')
        fontColorPlate  = 1-fontColorPlate;
        subfaultColor	= 1-subfaultColor;
        foreGcolor      = 1-foreGcolor;
        spreadColor  	= 1-spreadColor;
        crossColor      = 1-crossColor;
        BmarkerColor    = 1-BmarkerColor;
    end
    
    %% SETUP PARAMETERS
    if SWITCH.AxesLimit
        minX            = SWITCH.AxesLimitValues(1,1);
        maxX            = SWITCH.AxesLimitValues(1,2);
        minZ            = SWITCH.AxesLimitValues(1,3);
        maxZ            = SWITCH.AxesLimitValues(1,4);
    else
        minX            = min(GRID.X_3Dp(:));
        maxX            = max(GRID.X_3Dp(:));
        minZ            = min(GRID.Z_3Dp(:));
        maxZ            = max(GRID.Z_3Dp(:));
    end
    dX              = maxX-minX;
    dZ              = maxZ-minZ;
    plateTop        = dZ *1/2;
    plateThick      = dZ *1/7;
    plateBot        = plateTop-plateThick;
    plateMid        = plateTop-plateThick/2;
    widthSpread     = plateThick;
    widthSub        = plateThick;
    widthFault      = widthSub/10;
    depthSlab       = plateThick*3/4; %measured from plate bottom
    markerHeight1   = plateThick/4;
    SlabTipPosition = plateBot-depthSlab-plateThick; %this is reset below
    slabAngle       = atan((plateTop-SlabTipPosition)/widthSub);
    dxSlab          = sin(slabAngle)*plateThick;
    dzSlab          = cos(slabAngle)*plateThick;
    SlabTipPosition = plateBot-depthSlab-dzSlab;
    
    %% PLOTTING
    if (strcmp(PLATE.Subduction,'noTracking') && ~PLATE.StagnantLid) || ... %if plate boundaries could not be tracked
            (~PLOT.indicateTrench && ~PLOT.indicateRidge) || ... %no indication chosen
            strcmp(GRID.Type,'spherical2D')
        if strcmp(GRID.Type,'spherical2D')
            warning('plate parameters sketch not yet implemented for spherical2D geometry!')
        end
        if ~PLOT.indicateTrench && ~PLOT.indicateRidge
            warning('either PLOT.indicateTrench or PLOT.indicateRidge should be switched on to use Plate Sketch.')
        end
        %% CREATE AN EMPTY PLOT
        [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
        return
        
    else %if at least one plate boundary exists OR if there is a stagnant lid detected
        %% COUNT SUBPLOTS
        PLOT.nrSubplot     = PLOT.nrSubplot+1;
        
        %% SETUP PLOT
        % figure...
        figure(1)
        
        %subplot layout...
        SPOS.task   = 'createSubplot';
        [SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE);
        
        % CURRENT AXIS
        AXcurrent = gca; %save current axes handle
        
        if strcmp(GRID.Dim,'2-D')   % 2-D
            % DIMENSIONALIZATION
            % VH_3D       = VH_3D * SETUP.v0;  %[nd]->[cm/a]
            % T_3D        = T_3D * SETUP.deltaT; %[K]
            
            %% CREATE PLATE SKETCH
            minVelThreshold     = 0.001;    %[cm/a] below values are set to zero
            if PLATE.StagnantLid
                subX = [];
                spreadX = [];
                
            else
                subX            = PLATE.Subduction;
                spolarity       = PLATE.SubPolarity;
                subAngle        = PLATE.ShallowSlabAngle;
                bendingDiss     = PLATE.ViscDissBending;
                relBendingDiss  = PLATE.ViscDissBendingRel;
                deltaEtaSlabUM  = PLATE.SlabMantleViscDiff;
                
                vConv           = PLATE.PlateConvergence;
                %             vLP             = PLATE.LPvelocity;
                %             vUP             = PLATE.UPvelocity;
                vTrench         = PLATE.TrenchVelocity; %assuming perfectly single-sided subduction!
                vleftP          = PLATE.leftPvelocity;
                vrightP         = PLATE.rightPvelocity;
                
                spreadX         = PLATE.Spreading;
                vDiv            = PLATE.PlateDivergence;
                vridgeLP        = PLATE.ridgeLPvelocity;
                vridgeRP        = PLATE.ridgeRPvelocity;
                vRidge          = PLATE.RidgeVelocity; %assuming perfectly symmetric spreading velocities!
            end
            
            %adjust numbers
            if showAbsoluteNumbers
                if ~isempty(subX)
                    vrightP_plot  	= abs(vrightP);
                    vleftP_plot    	= abs(vleftP);
                    vTrench_plot    = abs(vTrench);
                end
                if ~isempty(spreadX)
                    vridgeLP_plot   = abs(vridgeLP);
                    vridgeRP_plot   = abs(vridgeRP);
                    vRidge_plot     = abs(vRidge);
                end
            else
                if ~isempty(subX)
                    vrightP_plot   	= vrightP;
                    vleftP_plot   	= vleftP;
                    vTrench_plot    = vTrench;
                end
                if ~isempty(spreadX)
                    vridgeLP_plot   = vridgeLP;
                    vridgeRP_plot   = vridgeRP;
                    vRidge_plot     = vRidge;
                end
            end
            %neglect very small values
            if ~isempty(subX)
                if abs(vrightP_plot)<minVelThreshold;   vrightP_plot = 0; end
                if abs(vleftP_plot)<minVelThreshold;    vleftP_plot = 0; end
                if abs(vTrench_plot)<minVelThreshold;   vTrench_plot = 0; end
            end
            if ~isempty(spreadX)
                if abs(vridgeLP_plot)<minVelThreshold;  vridgeLP_plot = 0; end
                if abs(vridgeRP_plot)<minVelThreshold;  vridgeRP_plot = 0; end
                if abs(vRidge_plot)<minVelThreshold;    vRidge_plot = 0; end
            end
            
            %dimensions
            if SWITCH.DimensionalMode || SWITCH.DimensionalInput
                dim_v = ' cm/a';
            else
                dim_v = ' ';
            end
            if indDimAtBot
                dim_v = '';
                if SWITCH.DimensionalMode || SWITCH.DimensionalInput
                    dim_vBot = '[cm/a]';
                else
                    dim_vBot = '[nd]';
                end
            end
            
            %background plate
            patch([minX maxX maxX minX],[plateTop plateTop plateBot plateBot],plateColor,'FaceAlpha',1,'EdgeColor','none')
            
            if ~isempty(spreadX) && PLOT.indicateRidge
                for ii=1:size(spreadX,1)
                    %spreading boundaries
                    patch([spreadX(ii,1)-widthSpread/2 spreadX(ii,1) spreadX(ii,1)+widthSpread/2],...
                        [plateBot plateTop plateBot],spreadColor,'FaceAlpha',1,'EdgeColor','none'); %bottom-up triangle
                    
                    %indicate ridges
                    patch([spreadX(ii,1)-markerHeight1/2 spreadX(ii,1) spreadX(ii,1)+markerHeight1/2],...
                        [plateTop plateTop+markerHeight1 plateTop],BmarkerColor,'FaceAlpha',1,'EdgeColor','none'); %triangle
                    
                    %indicate divergence velocities
                    string1 = ['{\color[rgb]{',num2str(plateColor),'}',arrowL,' ',num2str(vDiv(ii,1),precisionNum),dim_v,' ',arrowR,'}']; %divergence velocity
                    %indicate ridge velocities
                    if vRidge(ii,1)>minVelThreshold %->
                        dirInd_L=''; dirInd_R=arrowR;
                    elseif vRidge(ii,1)<-minVelThreshold %<-
                        dirInd_L=arrowL; dirInd_R='';
                    else %vTrench==0
                        dirInd_L=''; dirInd_R='';
                    end
                    string2 = ['{\color[rgb]{',num2str(foreGcolor),'}',dirInd_L,' ',num2str(vRidge_plot(ii,1),precisionNum),dim_v,' ',dirInd_R,'}']; %ridge velocity
                    string = {string2;string1};
                    if showSpreadingData
                        ht = text(spreadX(ii,1),plateTop+markerHeight1,string,...
                            'FontSize',fontSize,'FontName',fontName,'BackgroundColor',backGcolor,...
                            'VerticalAlignment','bottom','HorizontalAlignment','center');
                    end
                end
            end
            
            if ~isempty(subX) && PLOT.indicateTrench
                for ii=1:size(subX,1)
                    if strcmp(spolarity(ii,1),'right')
                        xval    = [subX(ii,1)-dxSlab subX(ii,1) subX(ii,1)+widthSub subX(ii,1)+widthSub-dxSlab]; %slab
                        zval    = [plateTop-dzSlab plateTop plateBot-depthSlab plateBot-depthSlab-dzSlab ]; %slab
                        xval2   = [subX(ii,1) subX(ii,1)+widthFault subX(ii,1)+widthFault+widthSub subX(ii,1)+widthSub]; %fault
                        zval2   = [plateTop plateTop plateBot-depthSlab plateBot-depthSlab]; %fault
                    elseif strcmp(spolarity(ii,1),'left')
                        xval    = [subX(ii,1)+dxSlab subX(ii,1) subX(ii,1)-widthSub subX(ii,1)-widthSub+dxSlab]; %slab
                        zval    = [plateTop-dzSlab plateTop plateBot-depthSlab plateBot-depthSlab-dzSlab ]; %slab
                        xval2   = [subX(ii,1) subX(ii,1)-widthFault subX(ii,1)-widthFault-widthSub subX(ii,1)-widthSub]; %fault
                        zval2   = [plateTop plateTop plateBot-depthSlab plateBot-depthSlab]; %fault
                    else %unknown polarity DOUBLE-SIDED
                        xval    = [subX(ii,1)+widthSub/2 subX(ii,1)+widthSub/2 subX(ii,1)-widthSub/2 subX(ii,1)-widthSub/2]; %slab
                        zval    = [plateBot-depthSlab-dzSlab plateTop plateTop plateBot-depthSlab-dzSlab]; %slab
                        xval2   = [subX(ii,1)+widthFault/2 subX(ii,1)-widthFault/2 subX(ii,1)-widthFault/2 subX(ii,1)+widthFault/2]; %fault
                        zval2   = [plateTop plateTop plateBot-depthSlab-dzSlab plateBot-depthSlab-dzSlab]; %fault
                    end
                    %slabs
                    patch(xval,zval,slabColor,'FaceAlpha',1,'EdgeColor','none')
                    %subduction faults
                    patch(xval2,zval2,subfaultColor,'FaceAlpha',1,'EdgeColor','none')
                    
                    %indicate trenches
                    patch([subX(ii,1)-markerHeight1/2 subX(ii,1) subX(ii,1)+markerHeight1/2],...
                        [plateTop+markerHeight1 plateTop plateTop+markerHeight1],BmarkerColor,'FaceAlpha',1,'EdgeColor','none'); %triangle
                    
                    %indicate convergence velocities & ...
                    string1 = ['{\color[rgb]{',num2str(plateColor),'}',arrowR,' ',num2str(vConv(ii,1),precisionNum),dim_v,' ',arrowL,'}']; %convergence velocity
                    %... trench velocities
                    if vTrench(ii,1)>minVelThreshold %->
                        dirInd_L=''; dirInd_R=arrowR;
                    elseif vTrench(ii,1)<-minVelThreshold %<-
                        dirInd_L=arrowL; dirInd_R='';
                    else %vTrench==0
                        dirInd_L=''; dirInd_R='';
                    end
                    string2 = ['{\color[rgb]{',num2str(foreGcolor),'}',dirInd_L,' ',num2str(vTrench_plot(ii,1),precisionNum),dim_v,' ',dirInd_R,'}']; %trench velocity
                    string = {string2;string1};
                    ht = text(subX(ii,1),plateTop+markerHeight1,string,...
                        'FontSize',fontSize,'FontName',fontName,'BackgroundColor',backGcolor,...
                        'VerticalAlignment','bottom','HorizontalAlignment','center');
                    % %               TO ADD A PATCH BEHIND THE TEXT USE:
                    %                 textPos = get(ht,'Extent');
                    %                 patch([textPos(1,1) textPos(1,1)+textPos(1,3) textPos(1,1)+textPos(1,3) textPos(1,1)],...
                    %                     [textPos(1,2) textPos(1,2) textPos(1,2)+textPos(1,4) textPos(1,2)+textPos(1,4)],...
                    %                     'w','FaceAlpha',1,'EdgeColor','none');
                    %                 uistack(ht,'top');
                    
                    %indicate upper-plate velocities
                    if logical(1)
                        if strcmp(spolarity,'left')
                            dummy   = vleftP;
                        elseif strcmp(spolarity,'right')
                            dummy   = vrightP;
                        else %UNKNOWN POLARITY
                            dummy   = vleftP;
                        end
                        if dummy>minVelThreshold %->
                            dirInd_L=''; dirInd_R=arrowR;
                        elseif dummy<-minVelThreshold %<-
                            dirInd_L=arrowL; dirInd_R='';
                        else %vTrench==0
                            dirInd_L=''; dirInd_R='';
                        end
                        horizontalShift = markerHeight1; %only absolute value
                        if strcmp(spolarity,'left')
                            string = [dirInd_L,' ',num2str(vleftP_plot(ii,1),precisionNum),dim_v,' ',dirInd_R];
                            dummy='right'; horizontalShift = -horizontalShift; xfac = 3;
                        elseif strcmp(spolarity,'right')
                            string = [dirInd_L,' ',num2str(vrightP_plot(ii,1),precisionNum),dim_v,' ',dirInd_R];
                            dummy='left'; xfac = 3;
                        else %UNKNOWN POLARITY
                            string = [dirInd_L,' ',num2str(vleftP_plot(ii,1),precisionNum),dim_v,' ',dirInd_R];
                            dummy='right'; horizontalShift = -horizontalShift/2; xfac = 1;
                        end %depending on subduction polarity
                        text(subX(ii,1)+xfac*horizontalShift,plateMid,string,...
                            'Color',plateColor,'FontWeight','bold','FontSize',fontSize,'FontName',fontName,'BackgroundColor',backGcolorPlate,...
                            'VerticalAlignment','middle','HorizontalAlignment',dummy); %bold background text
                        ht = text(subX(ii,1)+xfac*horizontalShift,plateMid,string,...
                            'Color',fontColorPlate,'FontSize',fontSize,'FontName',fontName,'BackgroundColor',backGcolorPlate,...
                            'VerticalAlignment','middle','HorizontalAlignment',dummy);
                    end
                    
                    %indicate lower-plate velocities
                    if strcmp(spolarity,'left')
                        dummy   = vrightP;
                    elseif strcmp(spolarity,'right')
                        dummy   = vleftP;
                    else %UNKNOWN POLARITY
                        dummy   = vrightP;
                    end
                    if dummy>minVelThreshold %->
                        dirInd_L=''; dirInd_R=arrowR;
                    elseif dummy<-minVelThreshold %<-
                        dirInd_L=arrowL; dirInd_R='';
                    else %vTrench==0
                        dirInd_L=''; dirInd_R='';
                    end
                    horizontalShift = markerHeight1; %only absolute value
                    if strcmp(spolarity,'left')
                        string = [dirInd_L,' ',num2str(vrightP_plot(ii,1),precisionNum),dim_v,' ',dirInd_R];
                        dummy='left'; xfac = 1;
                    elseif strcmp(spolarity,'right')
                        string = [dirInd_L,' ',num2str(vleftP_plot(ii,1),precisionNum),dim_v,' ',dirInd_R];
                        dummy='right'; horizontalShift = -horizontalShift; xfac = 1;
                    else %UNKNOWN POLARITY
                        string = [dirInd_L,' ',num2str(vrightP_plot(ii,1),precisionNum),dim_v,' ',dirInd_R];
                        dummy='left'; horizontalShift = horizontalShift/2; xfac = 1;
                    end %depending on subduction polarity
                    text(subX(ii,1)+xfac*horizontalShift,plateMid,string,...
                        'Color',plateColor,'FontWeight','bold','FontSize',fontSize,'FontName',fontName,'BackgroundColor',backGcolorPlate,...
                        'VerticalAlignment','middle','HorizontalAlignment',dummy); %bold background text
                    ht = text(subX(ii,1)+xfac*horizontalShift,plateMid,string,...
                        'Color',fontColorPlate,'FontSize',fontSize,'FontName',fontName,'BackgroundColor',backGcolorPlate,...
                        'VerticalAlignment','middle','HorizontalAlignment',dummy);
                    
                    %indicate shallow-slab angles & lithospheric bending dissipation
                    % %                 string = [char(10667),num2str(subAngle(ii,1),2),'',char(10666)];
                    % %                 string = [char(10658),num2str(subAngle(ii,1),2),'',char(10666)];
                    %                 string = [char(8738),' ',num2str(subAngle(ii,1),2),''];
                    %                 string = ['\theta = ',num2str(subAngle(ii,1),2),''];
                    string1 = ['\theta = ',num2str(subAngle(ii,1),2),''];
                    string2 = ['\phi^{vd}_{L,norm} = ',num2str(relBendingDiss(ii,1),2)];
                    string4 = ['\phi^{vd}_{L} = ',num2str(bendingDiss(ii,1),2),' W/m'];
                    string3 = ['\Delta\eta_{LA} = ',num2str(deltaEtaSlabUM,2)];
                    
                    string  = {string1;string2};
                    string  = string1;
                    string  = string2;
%                     string  = string4;
                    %                 string  = string3;
                    
                    ht = text(xval(1,4),zval(1,4),string,...
                        'Color',plateColor,'FontSize',fontSize,'FontName',fontName,'BackgroundColor',backGcolorPlate,...
                        'VerticalAlignment','top','HorizontalAlignment','center');
                end
            end
            
            %% TITLE STRING & TIME
            DESVARIA.Task    = 'create annotation strings';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            
            DESVARIA.Task    = 'make title';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            
            if ~isfield(SAVE,'AllFilesTitleString')
                SAVE.AllFilesTitleString = PLOT.titleStringCurrent;
            else
                SAVE.AllFilesTitleString = [SAVE.AllFilesTitleString, {PLOT.titleStringCurrent}];
            end
            
            %% SETUP AXIS
            DESVARIA.equalAxes  = false;
            DESVARIA.axesLimits = [minX maxX minZ maxZ]; %fix axis
            DESVARIA.removeAxes = true;
            DESVARIA.Task    = 'setup axes';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            
            %% ADD LABELS & INDICATIONS
            if indDimAtBot && ~PLATE.StagnantLid && ~SAVE.plateSketchIndication
                xLimits = get(gca,'XLim');  %Get the range of the x axis
                yLimits = get(gca,'YLim');  %Get the range of the y axis
                text(xLimits(1,2),SlabTipPosition,dim_vBot,...
                    'FontSize',fontSize,'FontName',fontName,'Color',plateColor,'BackgroundColor',backGcolor,...
                    'VerticalAlignment','bottom','HorizontalAlignment','right');
                SAVE.plateSketchIndication  = true;
            end
            if PLATE.StagnantLid
                string = 'stagnant lid';
                ht = text(minX+(maxX-minX)/2,plateMid,string,...
                    'Color',fontColorPlate,'FontSize',fontSize,'FontName',fontName,'BackgroundColor',backGcolorPlate,...
                    'VerticalAlignment','middle','HorizontalAlignment','center');
            end
            
            xlabel(PLOT.xlabel)
            ylabel(PLOT.ylabel)
            %set(AXcurrent,'YGrid','on')
            
            %     if SWITCH.Texting
            %         if SWITCH.DimensionalMode || SWITCH.DimensionalInput
            %             text(min(xlim),max(ylim),['min val: ',num2str(min(vh_plate(:)),2)],'HorizontalAlignment','left','VerticalAlignment','bottom','Color','FontName',STYLE.AllFontName,STYLE.keyColor,'FontSize',8)
            %             text(max(xlim),max(ylim),['max val: ',num2str(max(vh_plate(:)),2)],'HorizontalAlignment','right','VerticalAlignment','bottom','Color','FontName',STYLE.AllFontName,STYLE.keyColor,'FontSize',8)
            %             disp(['   max.vplate: ',num2str(max(vh_plate)),' / min.vplate: ',num2str(min(vh_plate))])
            %         else
            %         end
            %     end
            
        elseif strcmp(GRID.Dim,'3-D') %3-D
            error('plate sketch not yet implemented for 3-D!')
        end
    end
    PLOT.cb = NaN; 	%this is a flag if no colorbar is used
end
    
 


%%                       PLOT TABLE WITH VALUES OF TECTONIC PARAMETERS 1.1

%                                       	     Fabio Crameri, 09.06.2016
%% NOTES
% not implemented for 3-D yet...........

if strcmp(FIELD.name,'Parameter table')
    %% INPUT
    indDimAtBot             = logical(1);
    fontSize                = 11;
    fontName                = STYLE.AllFontName;
    backGcolor              = 'none';
    backGcolorPlate         = 'none';
    fontColorPlate          = [1 1 1];
    
    showAbsoluteNumbers     = logical(0);
    
    %% DESIGN INPUT
    precisionNum            = 2;
    
    %% COLOUR ADJUSTMENTS
    if strcmp(STYLE.ColorMode,'dark')
        
    end
    
    %% SETUP PARAMETERS
    
    
    %% PLOTTING
    if (strcmp(PLATE.Subduction,'noTracking') && ~PLATE.StagnantLid) || ... %if plate boundaries could not be tracked
            (~PLOT.indicateTrench && ~PLOT.indicateRidge) %no indication chosen
        if strcmp(GRID.Type,'spherical2D')
            warning('parameter table not yet implemented for spherical2D geometry!')
        end
        if ~PLOT.indicateTrench && ~PLOT.indicateRidge
            warning('either PLOT.indicateTrench or PLOT.indicateRidge should be switched on to use parameter table.')
        end
        %% CREATE AN EMPTY PLOT
        [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
        return
        
    else %if at least one plate boundary exists OR if there is a stagnant lid detected
        %% COUNT SUBPLOTS
        PLOT.nrSubplot = PLOT.nrSubplot+1;
        
        %% SETUP PLOT
        % figure...
        figure(1)
        % subplot...
        SPOS.task   = 'createSubplot';
        %         if isfield(SAVE,'PlotInPlot_FIELD2D_hax')
        %             SPOS.haxPlotInPlot = SAVE.PlotInPlot_FIELD2D_hax;
        %         end
        [SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE);
        
        % CURRENT AXIS
        AXcurrent = gca; %save current axes handle
        
        if strcmp(GRID.Dim,'2-D')   % 2-D
            %% CREATE PARAMETER TABLE
            if PLATE.StagnantLid
                subX = [];
                spreadX = [];
                
            else
                subX            = PLATE.Subduction;
                spolarity       = PLATE.SubPolarity;
                subAngle        = PLATE.ShallowSlabAngle;
                vConv           = PLATE.PlateConvergence;
                vLP             = PLATE.LPvelocity;
                vUP             = PLATE.UPvelocity;
                vTrench         = PLATE.TrenchVelocity; %assuming perfectly single-sided subduction!
                vleftP          = PLATE.leftPvelocity;
                vrightP         = PLATE.rightPvelocity;
                
                spreadX         = PLATE.Spreading;
                vDiv            = PLATE.PlateDivergence;
                vridgeLP        = PLATE.ridgeLPvelocity;
                vridgeRP        = PLATE.ridgeRPvelocity;
                vRidge          = PLATE.RidgeVelocity; %assuming perfectly symmetric spreading velocities!
            end
            
            %adjust numbers
            if showAbsoluteNumbers
                if ~isempty(subX)
                    vrightP_plot  	= abs(vrightP);
                    vleftP_plot    	= abs(vleftP);
                    vLP_plot    	= abs(vLP);
                    vUP_plot    	= abs(vUP);
                    vTrench_plot    = abs(vTrench);
                end
                if ~isempty(spreadX)
                    vridgeLP_plot   = abs(vridgeLP);
                    vridgeRP_plot   = abs(vridgeRP);
                    vRidge_plot     = abs(vRidge);
                end
            else
                if ~isempty(subX)
                    vrightP_plot   	= vrightP;
                    vleftP_plot   	= vleftP;
                    vLP_plot    	= vLP;
                    vUP_plot    	= vUP;
                    vTrench_plot    = vTrench;
                end
                if ~isempty(spreadX)
                    vridgeLP_plot   = vridgeLP;
                    vridgeRP_plot   = vridgeRP;
                    vRidge_plot     = vRidge;
                end
            end
            
            if SWITCH.DimensionalMode || SWITCH.DimensionalInput
                dim_v = ' cm/a';
            else
                dim_v = ' ';
            end
            if indDimAtBot
                dim_v = '';
                if SWITCH.DimensionalMode || SWITCH.DimensionalInput
                    dim_vBot = '[cm/a]';
                else
                    dim_vBot = '[nd]';
                end
            end
            
            %% SETUP TABLE DATA MATRIX
            dat1  	= []; %initialise data matrix
            dat2    = []; %initialise data matrix
            dat3    = []; %initialise data matrix
            ntables = 0;
            if ~isempty(subX) && PLOT.indicateTrench
                ntables = ntables+1;
                for ii=1:size(subX,1)
                    dat1 = [dat1;{['  SUBDUCTION ZONE ',num2str(ii)],[],''}];
                    dat1 = [dat1;{'Trench location',num2str(subX(ii,1),4),GRID.Xdim}];
                    dat1 = [dat1;{'Polarity',spolarity{ii,1},'-'}];
                    dat1 = [dat1;{['Shallow-slab dip (at ',num2str(PLATE.ShallowSlabAngleDepth,3),' ',GRID.Zdim,')'],num2str(subAngle(ii,1),2),''}];
                    dat1 = [dat1;{'Convergence',num2str(vConv(ii,1),2),'cm/a'}];
                    dat1 = [dat1;{'Trench/UP velocity',num2str(vTrench_plot(ii,1),2),'cm/a'}];
                    %                 dat1 = [dat1;{'UP velocity',num2str(vUP_plot(ii,1),2),'cm/a'}];
                    %dat = [dat;{'left-plate velocity',num2str(vleftP(ii,1),2),'cm/a'}];
                    dat1 = [dat1;{'LP velocity',num2str(vLP_plot(ii,1),2),'cm/a'}];
                    dat1 = [dat1;{'LP thickness',num2str(PLATE.PlateThicknessLP(ii,1),3),GRID.Zdim}];
                    dat1 = [dat1;{'UP thickness',num2str(PLATE.PlateThicknessUP(ii,1),3),GRID.Zdim}];
                    %dat1 = [dat1;{'Bending radius',num2str(PLATE.ViscDissBending(ii,1),3),GRID.Xdim}];
                    dat1 = [dat1;{'Bending dissipation',num2str(PLATE.ViscDissBending(ii,1),2),['m',STYLE.SCHAR.superscript2,'/s',STYLE.SCHAR.superscript2]}];
                    dat1 = [dat1;{'Rel. bend. dissipation',num2str(PLATE.ViscDissBendingRel(ii,1),2),'-'}];
                    dat1 = [dat1;{'Slab-mantle visc. ratio',num2str(PLATE.SlabMantleViscDiff(ii,1),2),'-'}];
                end
            else
                dat2 = [dat2;{['NO SUBDUCTION ZONE',num2str(ii)],[],''}];
            end
            if ~isempty(spreadX) && PLOT.indicateRidge
                ntables = ntables+1;
                for ii=1:size(spreadX,1)
                    %spreading boundaries
                    dat2 = [dat2;{['  SPREADING ZONE ',num2str(ii)],[],''}];
                    dat2 = [dat2;{'Ridge location',num2str(spreadX(ii,1),4),GRID.Xdim}];
                    dat2 = [dat2;{'Divergence',num2str(vDiv(ii,1),2),'cm/a'}];
                    dat2 = [dat2;{'Ridge velocity',num2str(vRidge_plot(ii,1),2),'cm/a'}];
                end
            else
                dat2 = [dat2;{['NO SPREADING ZONE',num2str(ii)],[],''}];
            end
            dat3 = [dat3;{['  OTHER'],[],''}];
            %         dat3 = [dat3;{'Slab-mantle visc. ratio',num2str(PLATE.SlabMantleViscDiff(ii,1),4),'-'}];
            if size(dat3,1)>1
                ntables = ntables+1;
            end
            
            %% CREATE TABLES
            pos     = get(gca,'position'); %save current subplot position
            %pos = plotboxpos(gca);
            for itable=1:ntables
                % columnname      = {'Parameter', 'Value', 'Units'};
                columnname      = [];
                % columnformat    = {'char', 'numeric', 'char'};
                columnformat    = {'char', 'char', 'char'};
                if itable==1
                    datCurrent	= dat1;
                elseif itable==2
                    datCurrent	= dat2;
                elseif itable==3
                    datCurrent	= dat3;
                else
                    error('check here')
                end
                
                htab(itable)  = uitable('Units','normalized','Data',datCurrent,...
                    'ColumnName',columnname,'ColumnFormat',columnformat,'RowName',[]);
                
                htab(itable).ColumnWidth = {120,45,40};
                % set(htab,'units','normalized') %adjust table extent to current subplot position
                % set(htab,'position',pos)
                %             tableHeight = min(htab(itable).Extent(4)+0.0001, pos(4));  %make sure it doesn't exceed axis height
                %tableHeight = htab(itable).Extent(4); %full table is shown, but might overlab other subplots
                if itable==1
                    htab(itable).Position(1) = pos(1);
                    htab(itable).Position(2) = pos(2)+pos(4)-htab(itable).Extent(4);
                elseif itable==2
                    htab(itable).Position(1) = pos(1)+pos(3)/2;
                    htab(itable).Position(2) = pos(2)+pos(4)-htab(itable).Extent(4);
                elseif itable==3
                    htab(itable).Position(1) = pos(1)+pos(3)/2;
                    htab(itable).Position(2) = pos(2)+pos(4)-htab(2,1).Extent(4)-htab(itable).Extent(4);
                end
                htab(itable).Position(3) = min(htab(itable).Extent(3)+0.0001, pos(3)); %make sure it doesn't exceed axis width
                htab(itable).Position(4) = pos(4); %make sure it doesn't exceed axis height
                
                
                %% STYLE TABLE
                htab(itable).ForegroundColor    = [0 0 0];
                htab(itable).FontName           = STYLE.AllFontName;
            end
            
            %% TITLE STRING & TIME
            DESVARIA.Task    = 'create annotation strings';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            
            DESVARIA.Task    = 'make title';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            
            if ~isfield(SAVE,'AllFilesTitleString')
                SAVE.AllFilesTitleString = PLOT.titleStringCurrent;
            else
                SAVE.AllFilesTitleString = [SAVE.AllFilesTitleString, {PLOT.titleStringCurrent}];
            end
            
            %% SETUP AXIS
            DESVARIA.equalAxes  = false;
            DESVARIA.removeAxes = true;
            DESVARIA.Task    = 'setup axes';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);

        elseif strcmp(GRID.Dim,'3-D') %3-D
            error('Parameter Table not yet implemented for 3-D!')
        end
    end
    PLOT.cb = NaN; 	%this is a flag if no colorbar is used
end
 
    


%%                                        SURFACE FIELD-VARIATION PLOT 1.02

%                                       	     Fabio Crameri, 13.07.2017
%% NOTES

if strcmp(FIELD.name,'Surface field variation')
    %% DEFAULTS
    PLOT.SFVunsuccessful            = false;      	%flag for failed function
    PLOT.sfvTotalSurfaceArea        = NaN;
    PLOT.sfvDataMinValue            = NaN;
    PLOT.sfvDataMaxValue            = NaN;
    PLOT.sfvDataMeanValue           = NaN;
    PLOT.sfvDataStandardDeviation 	= NaN;

    %% COUNT SUBPLOTS
    PLOT.nrSubplot          = PLOT.nrSubplot+1;
    
    %% SETUP PLOT
    % figure...
    figure(1)
    
    %subplot layout...
    SPOS.task = 'createSubplot';
    [SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE);
    
    % CURRENT AXIS
    AXcurrent               = gca; %save current axes handle
    
    %% DIAGNOSE AND PLOT SURFACE VARIATION
    [PLOT,SAVE] = f_DiagnosticsSurface(FIELD_M,FIELD,FILE,GRID,SWITCH,PLOT,STYLE,SAVE,PLATE);
    
    %% PROBLEM CHECK
    if PLOT.SFVunsuccessful
        
        %% CREATE AN EMPTY PLOT
        PLOT.nrSubplot = PLOT.nrSubplot-1; %undo previous, unsuccessful subplot count
        [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
        return
        
    end
    
    %% DISPLAY INFORMATION
    if strcmp(PLOT.KindOfDataDummy,'graph') %graph data
        disp('   Surface Field Variation:');
    elseif strcmp(PLOT.KindOfDataDummy,'plateCoreGraph') %graph data
        disp('   Surface Field Variation (from plate core):');
    else %field data
        disp(['   Surface Field Variation (',num2str(PLOT.sfvDepthLevel,3),' ',GRID.Zdim,' depth):']);
    end
    disp(['     ',STYLE.SCHAR.smallBullet,' min/max               = ',num2str(PLOT.sfvDataMinValue,4),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLOT.sfvDataMaxValue,4),' ',FIELD_M{strcmp(FIELD_M(:,2),PLOT.sfvField),4}]);
    disp(['     ',STYLE.SCHAR.smallBullet,' mean                  = ',num2str(PLOT.sfvDataMeanValue,4),' ',FIELD_M{strcmp(FIELD_M(:,2),PLOT.sfvField),4}]);
    disp(['     ',STYLE.SCHAR.smallBullet,' median                = ',num2str(PLOT.sfvDataMedianValue,4),' ',FIELD_M{strcmp(FIELD_M(:,2),PLOT.sfvField),4}]);
    disp(['     ',STYLE.SCHAR.smallBullet,' standard deviation    = ',num2str(PLOT.sfvDataStandardDeviation,3),' ',FIELD_M{strcmp(FIELD_M(:,2),PLOT.sfvField),4}]);
    PLOT = rmfield(PLOT,'KindOfDataDummy');
    
    %% FUNCTION OUTPUT
    %indicate axis handle as graph subplot
    SAVE.GraphPlotHandles  	= [SAVE.GraphPlotHandles; AXcurrent];
    PLOT.cb = NaN; 	%this is a flag if no colorbar is used
end



%%                                                  STREAM-GRAPH PLOT 0.00

%                                       	     Fabio Crameri, 11.07.2017
%% NOTES

if strcmp(FIELD.name,'Stream graph')
    %% DEFAULTS
    StreamGraphUnsuccessful     = false;

    %% COUNT SUBPLOTS
    PLOT.nrSubplot          = PLOT.nrSubplot+1;
    
    %% SETUP PLOT
    % figure...
    figure(1)
    
    %subplot layout...
    SPOS.task = 'createSubplot';
    [SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE);
    
    % CURRENT AXIS
    AXcurrent               = gca; %save current axes handle
    
    %% DIAGNOSE AND PLOT STREAM-GRAPH
    VARIA.Task    = 'create stream line';
    [PLOT,VARIA] = f_Varia(VARIA,FILE,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE,SETUP);
    
    %% PROBLEM CHECK
    if StreamGraphUnsuccessful
        
        %% CREATE AN EMPTY PLOT
        PLOT.nrSubplot = PLOT.nrSubplot-1; %undo previous, unsuccessful subplot count
        [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
        return
        
    end
    
    %% DISPLAY INFORMATION
%     if strcmp(PLOT.KindOfDataDummy,'graph') %graph data
%         disp('   Surface Field Variation:');
%     elseif strcmp(PLOT.KindOfDataDummy,'plateCoreGraph') %graph data
%         disp('   Surface Field Variation (from plate core):');
%     else %field data
%         disp(['   Surface Field Variation (',num2str(PLOT.sfvDepthLevel,3),' ',GRID.Zdim,' depth):']);
%     end
%     disp(['     ',STYLE.SCHAR.smallBullet,' min/max               = ',num2str(PLOT.sfvDataMinValue,4),' ',char(10231),' ',num2str(PLOT.sfvDataMaxValue,4),' ',FIELD_M{strcmp(FIELD_M(:,2),PLOT.sfvField),4}]);
%     disp(['     ',STYLE.SCHAR.smallBullet,' mean                  = ',num2str(PLOT.sfvDataMeanValue,4),' ',FIELD_M{strcmp(FIELD_M(:,2),PLOT.sfvField),4}]);
%     disp(['     ',STYLE.SCHAR.smallBullet,' median                = ',num2str(PLOT.sfvDataMedianValue,4),' ',FIELD_M{strcmp(FIELD_M(:,2),PLOT.sfvField),4}]);
%     disp(['     ',STYLE.SCHAR.smallBullet,' standard deviation    = ',num2str(PLOT.sfvDataStandardDeviation,3),' ',FIELD_M{strcmp(FIELD_M(:,2),PLOT.sfvField),4}]);
%     PLOT = rmfield(PLOT,'KindOfDataDummy');
    
    %% FUNCTION OUTPUT
    %indicate axis handle as graph subplot
    SAVE.GraphPlotHandles  	= [SAVE.GraphPlotHandles; AXcurrent];
    PLOT.cb = NaN; 	%this is a flag if no colorbar is used
end



%%                                                        PLOT TRACERS 1.1

%                                                Fabio Crameri, 22.10.2017
%% NOTES
% not implemented for 3-D yet...........

if strcmp(FIELD.name,'Tracers')
    if strcmp(GRID.Dim,'3-D') % 3-D
        if strcmp(GRID.Dim,'3-D')
            warning('Tracers plot not implemented for 3-D geometry!')
        end
        %% CREATE AN EMPTY PLOT
        [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
        return
        
    end
    
    %% READ INPUT DATA
    %tracers
    if strcmp(GRID.Type,'yinyang')
        if isfield(PLOT,'TraDataYin') && isfield(PLOT,'TraDataYang')
            TraData    	= PLOT.TraDataYin; TraDataYang 	= PLOT.TraDataYang;
            [~,~,TraVarName,~,~,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Tracer Info',SWITCH);
        else
            [~,~,TraVarName,TraData,~,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Tracers',SWITCH);
            if ischar(TraData(1,1))
                warning off backtrace; warning(TraData); warning on backtrace;
                %% CREATE AN EMPTY PLOT
                [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
                return
            end
        end
    else %all other grid types
        if isfield(PLOT,'TraData')
            TraData    	= PLOT.TraData;
            [~,~,TraVarName,~,~,~,~,TraTime,TraRcmb,~,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Tracer Info',SWITCH);
        else
            %[NumTracers,NumVariables,TraVarName,TraData,TraIdealMass,NumTraceElements,OutgassedAmount,TraTime,TraRcmb,TimeStep,Aspect,nb] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Tracers',SWITCH);
            [~,~,TraVarName,TraData,~,~,~,TraTime,TraRcmb,~,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Tracers',SWITCH);
            if ischar(TraData(1,1))
                warning off backtrace; warning(TraData); warning on backtrace;
                %% CREATE AN EMPTY PLOT
                [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
                return
            end
            if TraData(1,1)==0 && TraData(2,1)==0 %exchange x and y
                dummy = TraData(:,1); TraData(:,1) = TraData(:,2); TraData(:,2) = dummy;
            end
            PLOT.TraData 	= TraData;
        end
    end
    clearvars dummy
    
    %% SETUP PARAMETERS
    minX            = min(GRID.X_3Dp(:));
    maxX            = max(GRID.X_3Dp(:));
    minZ            = min(GRID.Z_3Dp(:));
    maxZ            = max(GRID.Z_3Dp(:));
    
    %% FLIP DEPTH VECTOR AND ACCOUNT FOR STICKY-AIR LAYER
    if SWITCH.DimensionalInput
        TraData(:,3) 	= SETUP.D_dim-TraData(:,3);      %reverse z-axis  (1.)
    else
        TraData(:,3) 	= 1.0-TraData(:,3);              %reverse z-axis  (1.)
    end
    if SWITCH.ActualDepth
        TraData(:,3)  	= TraData(:,3)-SETUP.d_air;      %set z=0 to real surface
    end
    
    %% DIMENSIONALISE DATA
    TraTime             = TraTime*PLOT.timeConvert;
    TraRcmb             = TraRcmb*GRID.dimFactor;
    TraData(:,1:3)      = TraData(:,1:3).*GRID.dimFactor;
    
    %% COUNT SUBPLOTS
    PLOT.nrSubplot   	= PLOT.nrSubplot+1;
    
    %% SETUP PLOT
    % figure...
    figure(1)
    
    %subplot layout...
    SPOS.task   = 'createSubplot';
    [SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE);
    
    % CURRENT AXIS
    AXcurrent = gca; %save current axes handle
    
    %% ADJUSTING DATA
    if strcmp(GRID.Type,'spherical2D')
        xPoints     = GRID.x2ds;
        zPoints     = GRID.z2ds;
    else
        xPoints     = GRID.x2dp;
        zPoints     = GRID.z2dp;
    end
    
    %% PLOTTING
    if verLessThan('matlab','9.1.0') %MATLAB 2016a and earlier
        IndexC      = strfind(TraVarName,PLOT.TracerVariable{iPlotLocal});
        idxType   	= find(not(cellfun('isempty',IndexC)));
    else
        idxType     = find(contains(TraVarName,PLOT.TracerVariable{iPlotLocal}));
    end
    %setup for zoom plot
    numLoopsZoom = 1;
    if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
        numLoopsZoom = 2;
    end
    %loop for zoom plot
    AXorig  = gca;
    for iZoom=1:numLoopsZoom
        if iZoom==2  %ZOOM PLOT
            DESVARIA.Task   = 'create magnifier';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            haxZ            = DESVARIA.haxZ;
        end
        
        %plotting
        if strcmp(PLOT.TracerVariable{iPlotLocal},'Position')
            ShowColorbar          	= false;
            %plotting
            plot(TraData(:,1),TraData(:,3),'.','Color',STYLE.keyLineColor)
            %scatter(TraData(:,1),TraData(:,3))
            
        elseif strcmp(PLOT.TracerVariable{iPlotLocal},'Type')
            DESVARIA.fieldName   	= ['Tracer ',PLOT.TracerVariable{iPlotLocal}];
            FIELD.symbol            = 'T_T';
            FIELD.dim               = '';
            FIELD.varScale          = 1;
            MarkerSize              = 5;
            ShowColorbar          	= true;
            SWITCH.ColormapAll      = {'davos', 'noflip'};
            FIELD.discreteCB        = true;
            %dimensionalisation
            TraData(:,idxType)  	= TraData(:,idxType).*FIELD.varScale;
            %plotting
            scatter(TraData(:,1),TraData(:,3),MarkerSize,TraData(:,idxType),'filled')
            
        elseif strcmp(PLOT.TracerVariable{iPlotLocal},'Mass')
            DESVARIA.fieldName         = ['Tracer ',PLOT.TracerVariable{iPlotLocal}];
            FIELD.symbol            = 'T_M';
            FIELD.dim               = '?';
            FIELD.varScale          = 1;
            MarkerSize              = 5;
            ShowColorbar          	= true;
            SWITCH.ColormapAll{1,1} = {'lajolla', 'noflip'};
            FIELD.discreteCB        = false;
            %dimensionalisation
            TraData(:,idxType)  	= TraData(:,idxType).*FIELD.varScale;
            %plotting
            scatter(TraData(:,1),TraData(:,3),MarkerSize,TraData(:,idxType),'filled')
            
        else
            warning(['PLOT.TracerVariable: "',PLOT.TracerVariable{iPlotLocal},'" not recognised!'])
            
            %% CREATE AN EMPTY PLOT
            PLOT.nrSubplot     = PLOT.nrSubplot-1;
            [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
            return
        end
    
        if iZoom==1 %main plot
            %% SAVING TRACER DATA TO FILE
            if ( SAVE.Tracers && ~SAVE.TracerSaved ) || logical(0) %save current tracers
                if ~isfield(SAVE,'Directory'); SAVE.Directory = FILE.directory; end
                SAVE.Directory              = FILE.directory;
                SAVE.DataName               = [FILE.name,'_Tracers',num2str(FILE.number)];
                SAVE.data                   = TraData;
                SAVE.dat                    = false;
                SAVE.txt                    = false;
                SAVE.mat                    = true;
                SAVE.write2existing         = false;
                [SAVE.overwriteAll] = f_saveData( SAVE );
                SAVE = rmfield(SAVE,'data');
                SAVE.TracerSaved            = true;
            end
            
            %% TITLE STRING & TIME
            DESVARIA.Task    = 'create annotation strings';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            
            DESVARIA.Task    = 'make title';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            
            if ~isfield(SAVE,'AllFilesTitleString')
                SAVE.AllFilesTitleString = PLOT.titleStringCurrent;
            else
                SAVE.AllFilesTitleString = [SAVE.AllFilesTitleString, {PLOT.titleStringCurrent}];
            end
            
            %% ANNOTATION
            xlabel(PLOT.xlabel)
            ylabel(PLOT.ylabel)
            if strcmp(GRID.Type,'spherical2D') && SWITCH.spherical2DCenterText
                text(0.0,0.0,[num2str(GRID.nx),'\times',num2str(GRID.nz)],'Color',STYLE.keyColor,'FontName',STYLE.AllFontName,...
                    'HorizontalAlignment','center')
            end
        end
        
        %% AXES SETUP
        if iZoom==1 %main plot
            % colorbar...
            if ShowColorbar
                DESVARIA.Task      = 'make colorbar';
                [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            else
                PLOT.cb         = NaN;  %this is the flag for no colorbar
            end
            
            if strcmp(GRID.Type,'spherical2D')
                BIndicationHeight       = 0; %adjust for plate boundary indications
                if PLOT.indicateTrench || PLOT.indicateRidge
                    BIndicationHeight   = abs(max(GRID.Z_3Dp(:))-min(GRID.Z_3Dp(:))) /10;
                end
                DESVARIA.axesLimits     = [min(GRID.x2ds(:))-BIndicationHeight,max(GRID.x2ds(:))+BIndicationHeight,min(GRID.z2ds(:))-BIndicationHeight,max(GRID.z2ds(:))+BIndicationHeight];
            end
            
            %% AXIS LIMIT BIG PLOT
            if ~strcmp(GRID.Type,'yinyang')
                if SWITCH.AxesLimit
                    if size(SWITCH.AxesLimitValues,1)>1  %individual axes values
                        DESVARIA.axesLimits     = SWITCH.AxesLimitValues(PLOT.loopCase,:);
                    else %all plots the same
                        DESVARIA.axesLimits     = SWITCH.AxesLimitValues;
                    end
                else
                    DESVARIA.axesLimits 	= [minX maxX minZ maxZ]; %fix axis to match model domain
                end
            end
        else %zoom plot
            if ishandle(PLOT.cb) && strcmp(get(PLOT.cb,'type'),'colorbar')
                if SWITCH.Colorbar(1,2)
                    set(gca,'clim',[FIELD.cColorbarMIN,FIELD.cColorbarMAX])
                else
                    cb = PLOT.cb;
                    set(gca,'clim',cb.Limits)
                end
            end
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            %lineColor = [0.9 0.9 0.9]; %white-ish
            %set(gca,'XColor',lineColor,'YColor',lineColor,'LineWidth',1)
            DESVARIA.axesLimits 	= PLOT.magnifierExtent(PLOT.loopCase,:);
        end
        if strcmp(GRID.Type,'Cartesian'); DESVARIA.flipAxes = true; end
        DESVARIA.Task    = 'setup axes';
        [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
        
        if iZoom==2; axes(AXorig); end %activate current big-plot axis
    end %zoom loop
end



%% COLORBAR
if ~ishandle(PLOT.cb) || ~strcmp(get(PLOT.cb,'type'),'colorbar')
    colorbar;   %define dummy colorbar
end

%% COLORMAP
if ~exist('AXorig','var'); AXorig  = gca; end
AXhandles(1)        = AXorig;  %big plot axis handle
if exist('haxZ','var')
    AXhandles(2)    = haxZ; %zoom plot axis handle
    numLoopsZoom    = 2;
else
    numLoopsZoom    = 1;
end
for iZoom=1:numLoopsZoom
    if ~SWITCH.multipleColormaps && iZoom>1; continue; end
    PLOT.hax        = AXhandles(iZoom);
    
    [SWITCH,~] = f_DesignColourmap(SWITCH,PLOT,FIELD,[],[]);
end

%% UNIVERSAL SWITCHES
if SWITCH.GridAlwaysOn
    grid on
    box on
end
set(gca,'Layer','top') % make sure axes ticks are visible

%% INSERT RECTANGLE INDICATING ZOOM PLOT
if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier && exist('haxZ','var')
    for i_rect=1:2
        hold on
        %rectangle('Position',[x,y,w,h]) %white bold in the back
        h_rectangle = rectangle('Position',[PLOT.magnifierExtent(PLOT.loopCase,1),...
            PLOT.magnifierExtent(PLOT.loopCase,3),...
            PLOT.magnifierExtent(PLOT.loopCase,2)-PLOT.magnifierExtent(PLOT.loopCase,1),...
            PLOT.magnifierExtent(PLOT.loopCase,4)-PLOT.magnifierExtent(PLOT.loopCase,3)]);
        h_rectangle.LineWidth = 0.5;
        h_rectangle.EdgeColor = [1 1 1];
        %rectangle('Position',[x,y,w,h])
        h_rectangle = rectangle('Position',[PLOT.magnifierExtent(PLOT.loopCase,1),...
            PLOT.magnifierExtent(PLOT.loopCase,3),...
            PLOT.magnifierExtent(PLOT.loopCase,2)-PLOT.magnifierExtent(PLOT.loopCase,1),...
            PLOT.magnifierExtent(PLOT.loopCase,4)-PLOT.magnifierExtent(PLOT.loopCase,3)]);
        h_rectangle.LineWidth = 0.25;
        h_rectangle.EdgeColor = [0 0 0];
    end
    
    axes(haxZ) %set zoom plot back on top
end

%% SIMPLIFY FIGURE
if SWITCH.SimplifyPlots
    %SET AXES ASPECT RATIO
    XaxisAR     = GRID.aspectRatio(1);
    XaxisARmax  = PLOT.MaxAspectX; 
    if strcmp(GRID.Type,'Cartesian') && ~strcmp(FIELD.name,'Surface field variation'); XaxisARmax = 999; end
    if strcmp(GRID.Type,'spherical2D') && strcmp(FIELD.name,'Grid'); XaxisAR = 1; end
    YaxisAR     = 1;
    ZaxisAR     = 1;
    set(gca,'PlotBoxAspectRatio',[min(XaxisARmax,XaxisAR) YaxisAR ZaxisAR])

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
if isfield(PLOT,'hax') && exist('haxZ','var')
    nrHax = size(PLOT.hax,1)+1; %new magnifier added
elseif isfield(PLOT,'hax')
    nrHax = size(PLOT.hax,1); %no new magnifier added
else
    nrHax = 1; %no magnifier at all
end
PLOT.hax(nrHax,1)       = AXcurrent;
SAVE.PlotHandles        = [SAVE.PlotHandles; AXcurrent];

if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier && exist('haxZ','var')
    SAVE.PlotHandles    = [SAVE.PlotHandles; haxZ]; %remember graph handles (zoom subplot)
    if isfield(PLOT,'haxZ'); nrHaxZ = size(PLOT.haxZ,1)+1; PLOT.haxZ(nrHaxZ,1) = haxZ;
    else; nrHaxZ = 1; PLOT.haxZ(nrHaxZ,1) = haxZ; end
end

%% FUNCTION SPECIFIC OUTPUT
% sFieldParameters       = 1; %just a dummy parameter so far



