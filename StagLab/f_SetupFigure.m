
%%                                              FIGURE SIZE & POSITION 2.61
% Input:
% FPOS.subplotLayout    = [2 3];
% [FPOS] = f_SetupFigure(FPOS,GRID,PLOT);
%                                                Fabio Crameri, 25.04.2020

function [FPOS] = f_SetupFigure(FPOS,SAVE,GRID,SWITCH)
%% defaults
if ~isfield(FPOS,'nrGraphs');        	FPOS.nrGraphs    	= 0;          	end
if ~isfield(FPOS,'GraphAspectX');    	FPOS.GraphAspectX 	= 2.0;       	end
if ~isfield(FPOS,'Default');            FPOS.Default        = 'TopRight'; 	end
if ~exist('SAVE','var') || ~isfield(SAVE,'app');	SAVE.app = 'SL_FieldPlot'; 	end

%% get current figure position
figPosCurr  = get(gcf,'Position');
FPOS.size   = figPosCurr(1,3:4);

if strcmp(SAVE.app,'SL_FieldPlot')
    if SWITCH.UsePanel
        
        %NEED ADJUSTMENTS HERE TOO, TO ACCOUNT FOR CYLINDTRICAL, 3-D AND YY ADJUSTMENTS DONE BELOW!
        
        standardUnitWidth       = 390;  %y = increasingFactor^(x)
        standardUnitHeight      = 350;
        marginWidth             = 80; %Approximately (in reference 1:1 plot)
        marginHeight            = 30; %Approximately (in reference 1:1 plot)
        maxFigureSize           = 1700;
        maxFigureWidth          = 1700;
        maxFigureHeight         = 1700;
        increasingFactor        = 0.8;
        %PlotBoxAspectRatio   	= pbaspect;
        
        %% get plot nx:nz aspect ratio
        %get subplot aspect ratio
        if isfield(FPOS,'setAxisLimit') && FPOS.setAxisLimit && SWITCH.AxesEqual
            spAspectRatio       = (FPOS.axisLimit(1,2)-FPOS.axisLimit(1,1))/(FPOS.axisLimit(1,4)-FPOS.axisLimit(1,3));
        else
            if strcmp(GRID.Dim,'2-D')
                spAspectRatio   = GRID.nx/GRID.nz;
                if strcmp(GRID.Type,'spherical2D')
                    if GRID.partialAnnulus
                        if GRID.aspectDomain(1,1)<=2.5
                            %nothing to be done
                        elseif GRID.aspectDomain(1,1)<=4.0
                            spAspectRatio 	= 1.5;
                        else
                            spAspectRatio 	= 1.0;
                        end
                    else
                        spAspectRatio 	= 1.0;
                    end
                end
            elseif strcmp(GRID.Dim,'3-D')
                if strcmp(GRID.Type,'yinyang')
                    if strcmp(SWITCH.yyPlotMode,'map')
                        spAspectRatio   = 2.2; %dummy value - maybe check!
                    else %mapOnSphere and isosurface
                        spAspectRatio   = 1.2; %dummy value - maybe check!
                    end
                else
                    %spAspectRatio   = GRID.ny/GRID.nz;
                    spAspectRatio   = 2; %due to the 3-D view, the panel aspect is always around an aspect ratio of 2
                end
            end
        end
        %% Widen standard extent corresponding to its aspect ratio
        if isnan(spAspectRatio) && ~strcmp(GRID.Type,'yinyang') %GRID.nz can be NaN
            if SWITCH.Verbose; warning('Plot aspect ratio could not be determined: Figure size might not be appropriate.'); end
        else
            %standardUnitWidth   = standardUnitWidth*spAspectRatio; %%%%%%%%% -marginWidth*spAspectRatio-1;
            standardUnitWidth   = standardUnitWidth*spAspectRatio -marginWidth*(spAspectRatio-1);
        end
    end
    
    %% set figure width and height
    % FigureWidth = NumSubPlot*SubPlotWidth + margin*(NumSubPlot+1)
    %'standard'
    if SWITCH.UsePanel
        %% Multiply extent with number subplots 
%         allPlotsFigureWidth   	= standardUnitWidth +maxFigureWidth-(increasingFactor^(FPOS.subplotLayout(1,2)-1)*maxFigureWidth); %WidthAddition = maxX-(0.6^(numX-1)*maxX)
%         allPlotsFigureHeight  	= standardUnitHeight +maxFigureHeight-(increasingFactor^(FPOS.subplotLayout(1,1)-1)*maxFigureHeight); %HeightAddition = maxZ-(0.6^(numZ-1)*maxZ)
        allPlotsFigureWidth   	= standardUnitWidth *FPOS.subplotLayout(1,2);
        allPlotsFigureHeight  	= standardUnitHeight *FPOS.subplotLayout(1,1);
        %MULTIPLE PANELS: account for less margin between multiple subplots (for two subplots, one margin is saved)
        if FPOS.subplotLayout(1,2)>1; allPlotsFigureWidth = allPlotsFigureWidth -marginWidth*(FPOS.subplotLayout(1,2)-1); end
        if FPOS.subplotLayout(1,1)>1; allPlotsFigureHeight = allPlotsFigureHeight -marginHeight*(FPOS.subplotLayout(1,1)-1); end
        %ASPECT RATIO: account for multiplying panel margins as well
        if isnan(spAspectRatio) && ~strcmp(GRID.Type,'yinyang') %GRID.nz can be NaN
            %nothing to be done.
        else
            allPlotsFigureWidth   	= allPlotsFigureWidth -marginWidth*spAspectRatio;
        end
        
        %% Limit to max. figure size
        allPlotsMaxFigureExtent	= max(allPlotsFigureWidth,allPlotsFigureHeight);
        diff2maxFigureSize      = maxFigureSize-allPlotsMaxFigureExtent;
        if diff2maxFigureSize<0 %current figure is larger than maximum figure size
            MultiplePanelsAdjustment = maxFigureSize-allPlotsMaxFigureExtent;
        else
            MultiplePanelsAdjustment = diff2maxFigureSize-(increasingFactor^(max(FPOS.subplotLayout(:))-1)*(diff2maxFigureSize));
        end
        FractionSizeAdjustment  = 1/allPlotsMaxFigureExtent*MultiplePanelsAdjustment;
        
        FPOS.size               = [allPlotsFigureWidth allPlotsFigureHeight];
        FPOS.size               = FPOS.size .*(1+FractionSizeAdjustment);
        
        
    else
        if FPOS.subplotLayout(1,1)==1 && FPOS.subplotLayout(1,2)==1;        FPOS.size = [600 300];       	%[1 1]
        elseif FPOS.subplotLayout(1,1)==1 && FPOS.subplotLayout(1,2)==2;    FPOS.size = [900 250];        	%[1 2]
        elseif FPOS.subplotLayout(1,1)==1 && FPOS.subplotLayout(1,2)==3;    FPOS.size = [1200 200];       	%[1 3]
        elseif FPOS.subplotLayout(1,1)==1 && FPOS.subplotLayout(1,2)==4;    FPOS.size = [1400 200];      	%[1 4] test
        elseif FPOS.subplotLayout(1,1)==1 && FPOS.subplotLayout(1,2)>4;     FPOS.size = [1400 200];       	%[1 big]
            
        elseif FPOS.subplotLayout(1,1)==2 && FPOS.subplotLayout(1,2)==1;    FPOS.size = [630 600];       	%[2 1]
        elseif FPOS.subplotLayout(1,1)==2 && FPOS.subplotLayout(1,2)==2;    FPOS.size = [930 479];        	%[2 2]
        elseif FPOS.subplotLayout(1,1)==2 && FPOS.subplotLayout(1,2)==3; 	FPOS.size = [1435 479];      	%[2 3]
        elseif FPOS.subplotLayout(1,1)==2 && FPOS.subplotLayout(1,2)==4;    FPOS.size = [1435 400];      	%[2 4] test
        elseif FPOS.subplotLayout(1,1)==2 && FPOS.subplotLayout(1,2)>4;     FPOS.size = [1700 300];       	%[2 big] test
            
        elseif FPOS.subplotLayout(1,1)==3 && FPOS.subplotLayout(1,2)==1;    FPOS.size = [530 800];       	%[3 1]
        elseif FPOS.subplotLayout(1,1)==3 && FPOS.subplotLayout(1,2)==2;    FPOS.size = [840 600];      	%[3 2]
        elseif FPOS.subplotLayout(1,1)==3 && FPOS.subplotLayout(1,2)==3;    FPOS.size = [1050 500];      	%[3 3]
        elseif FPOS.subplotLayout(1,1)==3 && FPOS.subplotLayout(1,2)==4;    FPOS.size = [1400 450];     	%[3 4] ok
        elseif FPOS.subplotLayout(1,1)==3 && FPOS.subplotLayout(1,2)==5;    FPOS.size = [1700 450];     	%[3 5] ok
        elseif FPOS.subplotLayout(1,1)==3 && FPOS.subplotLayout(1,2)>5;     FPOS.size = [1700 400];       	%[3 big] ok
            
        elseif FPOS.subplotLayout(1,1)==4 && FPOS.subplotLayout(1,2)==1;    FPOS.size = [430 900];      	%[4 1] ok
        elseif FPOS.subplotLayout(1,1)==4 && FPOS.subplotLayout(1,2)==2;    FPOS.size = [840 700];       	%[4 2] ok
        elseif FPOS.subplotLayout(1,1)==4 && FPOS.subplotLayout(1,2)==3;    FPOS.size = [1050 650];      	%[4 3]
        elseif FPOS.subplotLayout(1,1)==4 && FPOS.subplotLayout(1,2)==4;    FPOS.size = [1253 700];     	%[4 4]
        elseif FPOS.subplotLayout(1,1)==4 && FPOS.subplotLayout(1,2)>4;     FPOS.size = [1400 700];      	%[4 big] test
            
        elseif FPOS.subplotLayout(1,1)==5 && FPOS.subplotLayout(1,2)==1;    FPOS.size = [420 820];      	%[5 1] ok
        elseif FPOS.subplotLayout(1,1)==5 && FPOS.subplotLayout(1,2)==2;    FPOS.size = [840 820];      	%[5 2] ok
        elseif FPOS.subplotLayout(1,1)==5 && FPOS.subplotLayout(1,2)==3;    FPOS.size = [1150 820];       	%[5 3] ok
        elseif FPOS.subplotLayout(1,1)==5 && FPOS.subplotLayout(1,2)==4;    FPOS.size = [1450 820];       	%[5 4] ok
        elseif FPOS.subplotLayout(1,1)==5 && FPOS.subplotLayout(1,2)>4;     FPOS.size = [1450 820];       	%[5 big] test
            
        elseif FPOS.subplotLayout(1,1)>5 && FPOS.subplotLayout(1,2)==1;     FPOS.size = [420 820];      	%[big 1] test
        elseif FPOS.subplotLayout(1,1)>5 && FPOS.subplotLayout(1,2)==2;     FPOS.size = [840 1450];      	%[big 2] test
        elseif FPOS.subplotLayout(1,1)>5 && FPOS.subplotLayout(1,2)==3;     FPOS.size = [820 1150];       	%[big 3] test
        elseif FPOS.subplotLayout(1,1)>5 && FPOS.subplotLayout(1,2)==4;     FPOS.size = [820 820];       	%[big 4] test
            
        elseif FPOS.subplotLayout(1,1)<=FPOS.subplotLayout(1,2);           	FPOS.size = [1400 1000];      	%[big big] ok
        elseif FPOS.subplotLayout(1,1)>FPOS.subplotLayout(1,2);             FPOS.size = [1000 1400];      	%[big big] ok
            
        else
            disp(['...subplot layout [',num2str(FPOS.subplotLayout),'] is not available'])
            disp('   -->> using current figure size.')
        end
    end
    
    %% account for addition of zoom plots
    if SWITCH.Magnifier
        FPOS.size(1,1)      = FPOS.size(1,1)*11/10;         %make figure wider
    end
    
    if ~SWITCH.UsePanel
        %% account for plot nx:nz aspect ratio & 3-D plots
        %get subplot aspect ratio
        if isfield(FPOS,'setAxisLimit') && FPOS.setAxisLimit && SWITCH.AxesEqual
            spAspectRatio       = (FPOS.axisLimit(1,2)-FPOS.axisLimit(1,1))/(FPOS.axisLimit(1,4)-FPOS.axisLimit(1,3));
        else
            if strcmp(GRID.Dim,'2-D')
                spAspectRatio   = GRID.nx/GRID.nz;
            elseif strcmp(GRID.Dim,'3-D')
                spAspectRatio   = GRID.ny/GRID.nz;
            end
        end
        if isnan(spAspectRatio) && ~strcmp(GRID.Type,'yinyang') %GRID.nz can be NaN
            if SWITCH.Verbose; warning('Plot aspect ratio could not be determined: Figure size might not be appropriate.'); end
        end
    end
    if SWITCH.UsePanel
        %apply to appropriate figure size
        if strcmp(GRID.Dim,'2-D')
            if strcmp(GRID.Type,'Cartesian')
                %nothing to be done.
            elseif strcmp(GRID.Type,'spherical2D')
                if GRID.partialAnnulus
                    if GRID.aspectDomain(1,1)<=2.5
                        %nothing to be done
                    elseif GRID.aspectDomain(1,1)<=4.0
                        FPOS.size       = FPOS.size *1.1;        	%decrease overall size
                    else
                        FPOS.size       = FPOS.size *1.25;        	%increase overall size
                    end
                else
                    FPOS.size           = FPOS.size *1.25;         	%increase overall size
                end
            else
                error('no geometry found! check here.')
            end
        elseif strcmp(GRID.Dim,'3-D')
            if strcmp(GRID.Type,'Cartesian')
                FPOS.size(1,2) = FPOS.size(1,2) *1.;     %larger overall height; 1.3
                FPOS.size = FPOS.size *1.0;                 %increase overall size; 1.5
                
                if spAspectRatio==1 %for 1:1 aspect ratio
                elseif spAspectRatio==4 %for 4:1 aspect ratio
                    FPOS.size = FPOS.size *0.75;            %decrease overall size
                elseif spAspectRatio==8 %for 8:1 aspect ratio
                    FPOS.size = FPOS.size *0.7;           	%decrease overall size
                end
            elseif strcmp(GRID.Type,'yinyang')
                if strcmp(SWITCH.yyPlotMode,'map')
                    if FPOS.yyNumSlices==1
%                         FPOS.size(1,1)  = FPOS.size(1,1) *0.7; 	%overall width
%                         FPOS.size(1,2)  = FPOS.size(1,2) *1.1; 	%overall height
                        FPOS.size       = FPOS.size *1.2;    	%overall size
                    elseif FPOS.yyNumSlices==2
                        FPOS.size(1,2)  = FPOS.size(1,2) *1.2; 	%larger overall height
                        FPOS.size       = FPOS.size *1.0;    	%decrease overall size
                    elseif FPOS.yyNumSlices==3
                        FPOS.size(1,2)  = FPOS.size(1,2) *1.6; 	%larger overall height
                        FPOS.size       = FPOS.size *1.0;    	%decrease overall size
                    end
                else %mapOnSphere and isosurface
                    FPOS.size(1,2)  = FPOS.size(1,2) *1.4; 	%overall height
                    FPOS.size       = FPOS.size *0.9;    	%overall size
                end
                
            else
                error('no geometry found! check here.')
            end
            
            %special 3-D
            if FPOS.subplotLayout(1,1)==3 && FPOS.subplotLayout(1,2)==2; FPOS.size = [870 850]; end  %[3 2] 3d_evolution
        end
        
    else %no panel (old subplot)
        %apply to appropriate figure size
        if strcmp(GRID.Dim,'2-D')
            if strcmp(GRID.Type,'Cartesian')
                if spAspectRatio>=7 %for 8:1 aspect ratio
                    FPOS.size(1,1)      = FPOS.size(1,1) *2.3;     	%extend width
                    FPOS.size           = FPOS.size *0.7;           %decrease overall size
                elseif spAspectRatio>=3 %for 4:1 aspect ratio
                    FPOS.size(1,1)      = FPOS.size(1,1) *1.7;   	%extend width
                    FPOS.size           = FPOS.size *0.75;          %decrease overall size
                elseif spAspectRatio>=1.5 %for 2:1 aspect ratio
                    %nothing to do
                else %for 1:1 aspect ratio
                    FPOS.size(1,1)      = FPOS.size(1,1) *0.75;   	%shorten width
                end
            elseif strcmp(GRID.Type,'spherical2D')
                if GRID.partialAnnulus
                    if GRID.aspectDomain(1,1)<=2.5
                        %nothing to be done
                    elseif GRID.aspectDomain(1,1)<=4.0
                        FPOS.size(1,2)  = FPOS.size(1,2) *1.25;  	%extend height
                        FPOS.size       = FPOS.size *1.0;           %decrease overall size
                    else
                        %subplots are in fact nearly 1:1 aspect
                        FPOS.size(1,2)  = FPOS.size(1,2) *1.5;      %extend height
                        FPOS.size       = FPOS.size *1.0;        	%decrease overall size
                    end
                else
                    %subplots are in fact nearly 1:1 aspect
                    FPOS.size(1,2)      = FPOS.size(1,2) *1.5;      %extend height
                    FPOS.size           = FPOS.size *1.0;         	%decrease overall size
                end
            else
                error('no geometry found! check here.')
            end
            
            %% account for graph plots
            %undo if graph plots are present
            if FPOS.nrGraphs>0 && FPOS.GraphAspectX<spAspectRatio && strcmp(GRID.Dim,'2-D')
                FPOS.size(1,2)      = FPOS.size(1,2) *1.5;      %extend height
            end
            
        elseif strcmp(GRID.Dim,'3-D')
            if strcmp(GRID.Type,'Cartesian')
                FPOS.size(1,2) = FPOS.size(1,2) *1.3;     %larger overall height
                FPOS.size = FPOS.size *1.5;                 %increase overall size
                
                if spAspectRatio==1 %for 1:1 aspect ratio
                    FPOS.size(1,1) = FPOS.size(1,1) *0.75;     %shorten width
                elseif spAspectRatio==4 %for 4:1 aspect ratio
                    FPOS.size(1,1) = FPOS.size(1,1) *1.1;     %extend width
                    FPOS.size = FPOS.size *0.75;            %decrease overall size
                elseif spAspectRatio==8 %for 8:1 aspect ratio
                    FPOS.size(1,1) = FPOS.size(1,1) *1.45;     %extend width
                    FPOS.size = FPOS.size *0.7;           	%decrease overall size
                end
            elseif strcmp(GRID.Type,'yinyang')
                if strcmp(SWITCH.yyPlotMode,'map')
                    if FPOS.yyNumSlices==1
                        FPOS.size(1,2)  = FPOS.size(1,2) *0.8; 	%same overall height
                        FPOS.size       = FPOS.size *1.2;    	%same overall size
                    elseif FPOS.yyNumSlices==2
                        FPOS.size(1,2)  = FPOS.size(1,2) *1.8; 	%larger overall height
                        FPOS.size       = FPOS.size *0.9;    	%decrease overall size
                    elseif FPOS.yyNumSlices==3
                        FPOS.size(1,2)  = FPOS.size(1,2) *2.5; 	%larger overall height
                        FPOS.size       = FPOS.size *0.8;    	%decrease overall size
                    end
                else %mapOnSphere and isosurface
                    FPOS.size(1,2)  = FPOS.size(1,2) *1.4; 	%same overall height
                    FPOS.size       = FPOS.size *0.9;    	%same overall size
                end
                
            else
                error('no geometry found! check here.')
            end
            
            %special 3-D
            if FPOS.subplotLayout(1,1)==3 && FPOS.subplotLayout(1,2)==2; FPOS.size = [870 850]; end  %[3 2] 3d_evolution
        end
    end

elseif strcmp(SAVE.app,'SL_RadialProfile')
    % set figure width and height
    if FPOS.subplotLayout(1,1)==1 && FPOS.subplotLayout(1,2)==1;      	FPOS.size = [250 450];	%[1 1]
    elseif FPOS.subplotLayout(1,1)==1 && FPOS.subplotLayout(1,2)==2; 	FPOS.size = [500 450]; 	%[1 2]
    elseif FPOS.subplotLayout(1,1)==1 && FPOS.subplotLayout(1,2)==3; 	FPOS.size = [750 450];	%[1 3]
    elseif FPOS.subplotLayout(1,1)==1 && FPOS.subplotLayout(1,2)==4;  	FPOS.size = [1000 450];	%[1 4]
    elseif FPOS.subplotLayout(1,1)==2 && FPOS.subplotLayout(1,2)==3; 	FPOS.size = [670 700];	%[2 3]
    elseif FPOS.subplotLayout(1,1)==2 && FPOS.subplotLayout(1,2)==4;  	FPOS.size = [893 700];	%[2 4]
    elseif FPOS.subplotLayout(1,1)==2 && FPOS.subplotLayout(1,2)==5;  	FPOS.size = [1120 700];	%[2 5]
    elseif FPOS.subplotLayout(1,1)==2 && FPOS.subplotLayout(1,2)==6;  	FPOS.size = [1340 700];	%[2 6]
    else
        FPOS.size   = [1340 700]; %dummy
    end
    
elseif strcmp(SAVE.app,'SL_TimeGraph')
    % set figure width and height
    if FPOS.subplotLayout(1,1)==1 && FPOS.subplotLayout(1,2)==1;      	FPOS.size = [717 255];	%[1 1]
    elseif FPOS.subplotLayout(1,1)==2 && FPOS.subplotLayout(1,2)==1;   	FPOS.size = [478 330];	%[2 1]
    elseif FPOS.subplotLayout(1,1)==3 && FPOS.subplotLayout(1,2)==1;  	FPOS.size = [478 490];	%[3 1]
    elseif FPOS.subplotLayout(1,1)==4 && FPOS.subplotLayout(1,2)==1; 	FPOS.size = [478 655];	%[4 1]
    elseif FPOS.subplotLayout(1,1)==1 && FPOS.subplotLayout(1,2)==2;  	FPOS.size = [956 170];	%[1 2]
    elseif FPOS.subplotLayout(1,1)==3 && FPOS.subplotLayout(1,2)==2;  	FPOS.size = [956 490];	%[3 2]
    elseif FPOS.subplotLayout(1,1)==4 && FPOS.subplotLayout(1,2)==2; 	FPOS.size = [956 655];	%[4 2]
    elseif FPOS.subplotLayout(1,1)==5 && FPOS.subplotLayout(1,2)==2;  	FPOS.size = [956 820];	%[5 2]
    elseif FPOS.subplotLayout(1,1)==6 && FPOS.subplotLayout(1,2)==2;  	FPOS.size = [956 980];	%[6 2]
    else
        FPOS.size   = [956 980]; %dummy
    end
end

%% set figure position
screenSize  = get(groot,'Screensize');
if strcmp(FPOS.Default,'BottomLeft')
    FPOS.x      = 47;
    FPOS.y      = 1;
    
elseif strcmp(FPOS.Default,'BottomRight')
    FPOS.x      = screenSize(3)-FPOS.size(1,1);
    FPOS.y      = 1;
    
elseif strcmp(FPOS.Default,'TopLeft')
    FPOS.x      = 47;
    FPOS.y      = screenSize(4)-FPOS.size(1,2);
    
elseif strcmp(FPOS.Default,'TopRight')
    FPOS.x      = screenSize(3)-FPOS.size(1,1);
    FPOS.y      = screenSize(4)-FPOS.size(1,2);
    
else %if no match: bottom left
    FPOS.x      = 47;
    FPOS.y      = 1;
end

%% apply figure position
FPOS.position           = [FPOS.x FPOS.y FPOS.size(1,1) FPOS.size(1,2)];
set(gcf,'Position',FPOS.position);

%% output
if exist('spAspectRatio','var')
    FPOS.spAspectRatio	= spAspectRatio;
end

end
