
%%                                                          EMPTY PLOT 2.00
%   .calls f_DesignLayout
%   .calls f_DesignLayoutPosition
%   .calls f_DesignVaria
%                                                Fabio Crameri, 07.08.2018

function [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE)

%% VARIABLE DEFINITIONS
if strcmp(SAVE.app,'SL_FieldPlot')
    if strcmp(GRID.Type,'spherical2D')
        minX     	= min(GRID.x2ds(:));
        minZ      	= min(GRID.z2ds(:));
        maxX       	= max(GRID.x2ds(:));
        maxZ      	= max(GRID.z2ds(:));
    elseif strcmp(GRID.Type,'yinyang')
        minX     	= 0;
        minZ     	= 0;
        maxX      	= 1;
        maxZ      	= 1;
    else
        minX      	= min(GRID.X_3Dp(:));
        minZ    	= min(GRID.Z_3Dp(:));
        maxX    	= max(GRID.X_3Dp(:));
        maxZ      	= max(GRID.Z_3Dp(:));
    end
else
    minX            = 0;
    minZ            = 0;
    maxX            = 1;
    maxZ            = 1;
end
dZ                  = maxZ-minZ;
crossThick          = dZ *1/5;

%% DESIGN INPUT
crossColor          = [0.96 0.96 0.96];
annotationColor     = [0.75 0.75 0.75];
if strcmp(STYLE.ColorMode,'dark')
    crossColor      = 1-crossColor;
    annotationColor	= 1-annotationColor;
end

if strcmp(SAVE.app,'SL_FieldPlot')
    %% COUNT SUBPLOTS
    PLOT.nrSubplot  = PLOT.nrSubplot+1;
    
    %% SETUP FIGURE
    figure(1)
    
    %% SUBPLOT LAYOUT
    SPOS.task       = 'createSubplot';
    [SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE);
    
    %% CURRENT AXIS
    AXcurrent       = gca; %save current axes handle
end

%% DRAW CROSS PATCH
xleft               = minX; % + dX/4;
xright              = maxX; % - dX/4;
% dX2     = xright-xleft;
% mid     = minZ+dZ/2;
%     patch([xleft, xleft+plateThick, xleft+dX2/2, xright-plateThick, xright, xright-dX2/2+plateThick/2, xright, xright-plateThick, xleft+dX2/2, xleft+plateThick, xleft, xleft+dX2/2-plateThick/2],...
%           [minZ,  minZ,             mid-plateThick/2,     minZ,            minZ,     mid,                maxZ,   maxZ,              mid+plateThick/2,     maxZ,            maxZ, mid],...
%         plateColor,'FaceAlpha',0.1,'EdgeColor','none')
patch([xleft, xleft+crossThick, xright, xright-crossThick],...
    [maxZ,  maxZ,             minZ,            minZ],...
    crossColor,'FaceAlpha',1.0,'EdgeColor','none')
patch([xleft, xleft+crossThick, xright, xright-crossThick],...
    [minZ,  minZ,             maxZ,            maxZ],...
    crossColor,'FaceAlpha',1.0,'EdgeColor','none')

if strcmp(SAVE.app,'SL_FieldPlot')
    %% ANNOTATIONS
    %title
    DESVARIA.fieldName	= [FIELD.name,' n.a.'];
    DESVARIA.Task	= 'create annotation strings';
    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
    
    DESVARIA.Task	= 'make title';
    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
else
    fieldName   	= [FIELD.name,' n.a.'];
    text((maxX-minX)/2,(maxZ-minZ)/2,fieldName,'FontName',STYLE.AllFontName,'FontWeight','bold',...
        'FontSize',STYLE.keyFontSize+2,'Color',annotationColor,'HorizontalAlignment','center','VerticalAlignment','middle');
end

%adjust axis
axis([minX maxX minZ maxZ]) %fix axis
axis off
if strcmp(SAVE.app,'SL_FieldPlot') && strcmp(GRID.Type,'yinyang')
    axis equal
end

if strcmp(SAVE.app,'SL_FieldPlot')
    %% COLORBAR
    DESVARIA.Task    = 'make colorbar'; %make dummy colorbar
    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
    % colorbar;   %define dummy colorbar
    PLOT.cb = NaN; 	%this is a flag if no colorbar is used

    %% SIMPLIFY FIGURE ==================================
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
        PLOT.title              = 'Empty plot';  	%title string
        
        [PLOT] = f_DesignLayout(SWITCH,PLOT,GRID,FIELD);
    end
end

%.........
drawnow  %actually only needed if topo or other graph plot is included (axis change otherwise)
%.........

%% FUNCTION OUTPUT
if strcmp(SAVE.app,'SL_FieldPlot')
    if isfield(PLOT,'hax'); nr_hax = size(PLOT.hax,1)+1; else; nr_hax = 1; end
    PLOT.hax(nr_hax,1)  = AXcurrent;
    %indicate axis handle as graph subplot
    SAVE.GraphPlotHandles   = [SAVE.GraphPlotHandles; AXcurrent];
    SAVE.PlotHandles        = [SAVE.PlotHandles; AXcurrent];
end



