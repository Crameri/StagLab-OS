
%%                                             DEFAULT INPUT VARIABLES 3.21
%                                                         for SL_FieldPlot
% 
%                                                Fabio Crameri, 31.05.2021

function [IN,PLOT,SWITCH,TOPO,STYLE,SAVE] = f_Defaults(SAVE)
    if ~exist('SL_FieldPlot','file'); f_INSTALL; end
    %% fAIO ACTION: UPDATE OLD PARFILES
    if ~(ismcc || isdeployed) && (~isfield(SAVE,'count') || SAVE.count==0) %only first time during one execution
        fAIo.task   = 'updateOldParfilesMinimal';
        [~,SAVE] = f_AIo(fAIo,[],[],[],[],SAVE);
    end
    
    %% COLLECTING STAGLAB PERFORMANCE AND ERROR DATA **********************
    % StagLab collects user-specific performance and error data
    % that is stored in /StagLab/fAIo/cache.mat
    % The data is used to keep StagLab up-to-date, make parfiles forward-compatible,
    % and provide information in case of StagLab user-specific bugs and errors.
    SWITCH.AnonymiseStagLabData  = 	logical(0);     %Anonymise data
    % *********************************************************************
    
    %% DEBUGGING
    SWITCH.Verbose              =   logical(0);                         %Switch for extended display output
    SWITCH.sendErrorLog         =   logical(0);                         %Sends an error log file to the developer for debugging purposes
    SWITCH.spPause              =   0.0;                                %Pausing after subplot creation (might help with shifted subplots MatLab problem)
    
    %% INPUT FILE(S)
    IN.Name                     =   { 'test' 'test'};                   %Filename  (it's size controls the # of files processed)
    IN.Number                   =   [   1    1 ];                       %Filenumber (use iLoop to loop through multiple files)
    IN.Parameter                =   [  	11   11 ];                      %Parameter dimensionalisation (see f_Dimensions)
    IN.Folder                   =   {'/work/stagyy/' '/work/stagyy/'};  %File directory
%     PLOT.titleString            =   {'XXX' 'XXX' 'XXX'  };            %Adds manual titles if uncommented
%     PLOT.legendString           =   {'XXX' 'XXX' 'XXX'  };            %Adds manual legend names if uncommented (e.g., for graph plots)

    SWITCH.GeodynamicCode       =   'StagYY';                           %To use StagLab with a different code, run and adjust at locations of error messages
        SWITCH.StagYYoutput     =   'Binary';                           %'Binary' or 'HDF5' (only for use with StagYY)
        SWITCH.ASPECToutput     =   'HDF5';                             %Currently only 'HDF5' (only for use with ASPECT)
    SWITCH.checkForFolder       =   logical(1);                         %Scans locally for alternative folders
    SWITCH.Precision            =   {'single'};                         %Input data precision: 'single' or 'double'
    
    %% StagLab MODES
    SWITCH.DimensionalMode   	=   logical(0);   	%Converts to dimensional output as set in f_Dimensions.m
        SWITCH.modelSetup       =   'standard';   	%General model parameter setup; see f_Dimensions.m
    STYLE.ColorMode             =   'light';        %'light' or 'dark' background
    SWITCH.TimeEvolutionMode    =   logical(0);     %Adds subtle guide arrows between subplots indicating time evolution
        SWITCH.timeDirection    =   'topbot';     	%'topbot', 'leftright'
    SWITCH.AnalysisMode         =   logical(0); 	%Creates figures for analysis purposes - not recommended for publication!
    SWITCH.QuickMode         	=   logical(0);     %Adjust parfile switches for quick runtime to e.g. save data efficiently
    SWITCH.ModelSketchMode    	=   logical(0); 	%Creates simplified sketches to explain basic model setup

    %% FIGURE POSITIONING
    SWITCH.FigurePosition       =   'auto';         %'auto' or set manually by e.g., [47 1 1300 350]
    PLOT.FigureDefaultPosition 	=   'TopRight';     %Figure position on screen: 'BottomLeft', 'BottomRight', 'TopLeft', or 'TopRight'

    %% SUBPLOT LAYOUT
%    PLOT.Layout                 =   [2 2];     	%Manual subplot layout; only takes effect if uncommented; [z x]
    SWITCH.ReverseLayout       	=   false;          %e.g., plot 3x2 instead of 2x3
    SWITCH.BackgroundGuides  	=   logical(1);     %Adds subtle background guides to visually connect subplots
    PLOT.MaxAspectX             =   4.0;            %The limit for the subplot aspect-ratio in x-direction for graph plots
    PLOT.shiftLeft              =   0;              %shift subplot positions closer or further away from each other; e.g., 0.01
    PLOT.shiftUpward            =   0;              %shift subplot positions closer or further away from each other; e.g., 0.01
    SWITCH.UsePanel             =   logical(1);     %Uses novel subplot distribution (still under testing)
    
    %% PLOT STYLING
    SWITCH.PlotDesign           =   logical(1);     %Plot design improvements
    STYLE.Mode                  =   'Custom';       %'Custom': adjustable plot design | 'Nature','Science','EGU','AGU','Frontiers': preset plot design
    STYLE.AllFontName           =   'Archivo';  	%e.g. 'Helvetica','ETH Light','Archivo'
    STYLE.keyFontSize           =   10;
    STYLE.Brackets4Dimensions   =   'none';         %Puts dimensions into 'Parentheses','Brackets', or 'none'
    SWITCH.BackgroundDesign    	=   logical(1);  	%Plot-background design improvements for graph supplots
    SWITCH.LegendDesign        	=   logical(1);     %Legend design improvements
        STYLE.fitLegendPosition	= 	logical(1);   	%Fixes legends tightly to the axes
    SWITCH.SimplifyPlots      	=   logical(1);     %Simplify axes description
    
    PLOT.Style                  =   'pcolor';       %Either 'contour' or 'pcolor'
        PLOT.numContours        =   100;            %Number of contours in contour plots
    SWITCH.closeAnnulus         =   logical(1);     %Adds one row of data points to visually close the spherical2D plots
    %axes
    SWITCH.AxesOff              =   logical(0);
    SWITCH.onlyXaxis            = 	logical(1);     %Currently only implemented for spherical2D plots
    SWITCH.AxesEqual            =   logical(1);
    SWITCH.AxesLimit            =   logical(0);     
        SWITCH.AxesLimitValues  =   [500 1300 -5 200; 1500 2300 -5 200; 500 1300 -5 200; 1600 2400 -5 200]; %Axes limit for big plots; [xmin xmax ymin ymax] in actual plotting dimensions
    SWITCH.ActualDepth          =   logical(1);     %Accounts for sticky air - only applied for sticky air cases
    SWITCH.TickDirection    	=   'out';          %'in' or 'out' - tickmarks location
    STYLE.removeAxisRuler    	=   logical(1);     %Removes axis ruler in 2-D plots
    SWITCH.GridAlwaysOn         =   logical(0);
    %colours
    STYLE.keyColor              =   [0.5 0.5 0.5];  %Key colour used for text and axes
    STYLE.keyLineColor          =   [0.1 0.1 0.1];  %Key line colour in graph plots
    STYLE.keyLineWidth          =   1.0;            %Key line width in graph plots
    STYLE.annotationColor      	=   [0.2 0.2 0.2];  %Key colour for annotation text
    STYLE.annotationBackColor  	=   [0.98 0.98 0.98];  %Key colour for annotation backgrounds
    STYLE.CategoricalColours  	=   'batlow';     	%Name of a Scientific colour map (www.fabiocrameri.ch/colourmaps) to be used for lines or markers
    %other
    PLOT.caseColor              =   [ 0.4 0.4 0.4   ;  0.8 0.4 0.0  ;  0.8 0.0 0.0];
    PLOT.lineStyle              =   {'-'    '--'    ':'};
    PLOT.lineWidth              =   [ 1.25   1.25   1.25  ];
    PLOT.startGradColor         =   [0.7 0.7 0.7];  %Starting color for gradual specific file coloring (e.g., for PlotInPlot)
    
    %% ANNOTATIONS
    SWITCH.Annotation           =   logical(1);     %Subplot annotations with e.g., 'a', 'b', 'c'
        STYLE.annotationLocation=   'botLeft';    	%e.g., 'topOutTop','topOutCorner','topLeft','botLeft','topRight','botRight'
        PLOT.annotationString  	=   {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' 's' 't' 'u' 'v' 'w' 'x' 'y' 'z'};  %subplot annotation strings
        PLOT.annIfMoreThanOne   =   logical(1);     %Adds annotation only if more than 1 subplot
    SWITCH.Title                =   logical(1);
        PLOT.titleTimeOnlyOnce  =   logical(0);     %Only displays time in title once per figure
        PLOT.titleOnlyTime      =   logical(0);   	%For a title showing only time
        PLOT.titleNoTime        =   logical(0);     %For a title showing no time
    SWITCH.Texting              =   logical(0);     %Write out resolution and max values of the field
    SWITCH.spherical2DCenterText=   logical(0);     %Adds #gridpoints to the center of the spherical annulus
    
    %% PLOT ADDITIONS
    SWITCH.Magnifier          	=   logical(0);     %Adds an additional plot for a close up
        PLOT.magnifierPosition 	=   'default';      %'default', 'northwest', 'southwest', 'west', 'centre', 'north', 'south', 'east', 'northeast', 'southeast'
        PLOT.magnifierExtent  	=   [2500 3500 -50 500 ; 900 1900 -50 500 ; 1300 2300 -50 500]; %[xmin xmax ymin ymax] in actual plotting dimensions
    PLOT.zProfile               =   logical(0);     %Takes profile at non-dim value of x (or z if PLOT.z_prof uncommented)
        PLOT.x_prof             =   1.01;   
        %PLOT.z_prof             =   0.1;
    PLOT.indicateDepth          =   logical(0);     %Adds horizontal depth indicator lines
        PLOT.indicateDepthLevels=   [410, 660];     %Specify one or multiple depth levels

    %% COLORBAR
    SWITCH.Colorbar             =   logical(1);     
        SWITCH.ConstantColorbar	=   logical(0);     %Fix colorbar limits to preset values
      	PLOT.cbPos              =   'right';        %'right' 'bottom'
    	PLOT.cbLabelPos         =   'right';        %'right' 'top'
%     %use these to set the min/max colorbar values of each field:
%     IN.graph_min            =   -15;          
%     IN.graph_max            =   15;
%     IN.topo_min             =   -0.3;
%     IN.topo_max             =   0.3;
%     IN.isotopo_min      	  =   -0.3;
%     IN.isotopo_max      	  =   0.3;
%     IN.restopo_min      	  =   -0.3;
%     IN.restopo_max      	  =   0.3;
%     IN.dyntopo_min      	  =   -0.3;
%     IN.dyntopo_max      	  =   0.3;
%     IN.Ph_min               =   0;
%     IN.Ph_max               =   3;
%     IN.C_min                =   1;
%     IN.C_max                =   4;
%     IN.T_min                =   300;
%     IN.T_max                =   1150;
%     IN.Tres_min             =   -400;
%     IN.Tres_max         	  =   4000;
%     IN.eta_min              =   1e19; 
%     IN.eta_max              =   1e28; 
%     IN.str_min              =   0;
%     IN.str_max              =   5.8;
%     IN.nstr_min             =   -800;
%     IN.nstr_max             =   800;
%     IN.rho_min              =   3200;
%     IN.rho_max              =   3450;
%     IN.P_min                =   0;            
%     IN.P_max                =   140;
%     IN.diss_min             =   0;            
%     IN.diss_max             =   5e-6;
%     IN.sfv_min              =   0;            
%     IN.sfv_max              =   0.5;
%     IN.water_min            =   0;
%     IN.water_max            =   1;
%     % ...and each individual subplot use brackets:
%     IN.edot_max             =   [3e-20  1.05e-18];
%     IN.edot_min             =   [0      0       ];
%     IN.vel_min              =   -0.01;
%     IN.vel_max              =   0.01;
%     IN.sfun_min             =   -200;
%     IN.sfun_max             =   200;

	%% COLOUR MAP
  	SWITCH.ColormapAll          =   {'auto', 'noflip'}; %e.g. 'auto' (to adjust automatically), or 'davos' // 'flip' or 'noflip'
    SWITCH.MultiManualColourMaps=   logical(0);         %Uses multiple different manual colour-maps (defined in PLOT.MultiManualColourMaps) for each field
    PLOT.MultiManualColourMaps  =   {'davos', 'noflip'; 'batlow', 'noflip'};
    SWITCH.AlternativeColourmaps=   logical(0);         %Uses alternative colour-map suite
    SWITCH.DiscreteColormap     =   logical(0);
    PLOT.NumberColormapColors 	=   NaN;                %Give number of discrete colours or set 'NaN' to have it automatically adjusted to number of tick labels
    SWITCH.enforceMinMaxCbLabels=   logical(0);         %Enforces tick labels for min and max colorbar values
    SWITCH.cbarLimitIndicTop    =   [0]; SWITCH.cbarLimitColorTop = [0.81 0.81 0.81]; %Adds specific color to top/bot to indicate exceeding values
    SWITCH.cbarLimitIndicBot    =   [0]; SWITCH.cbarLimitColorBot = [0.81 0.81 0.81]; %(only one number sets it for all plots, otherwise specify for every subplot using [0 1 0 ...])
    
    %% 3-D CARTESIAN MODEL: ISOSURFACE  -----------------------------------
    PLOT.cameraPosition         =   [59240 35752 -37631];  	%Camera position for 3-D figures; standard is top: [59240 35752 -37631]; bot: [59240 35752 +37631]
    PLOT.cameraPosBotView       =   logical(0);       	%Clips view to bottom up (only for 3-D box plots)
    PLOT.camProjection          =   'orthographic';  	%Camera projection: either 'perspective' or 'orthographic'
    PLOT.lighting               =   'gouraud';          %Lighting options, e.g., none, flat, gouraud
    PLOT.camlight               =   'headlight';     	%Camlight options, e.g., headlight, left, right
    PLOT.opacity                =   [0.3 0];         	%Opacity of the isosurfaces   1: last one is opaque 0: last one is also transparent]
    PLOT.isoVal                 =   [0.2 0.45 0.8];     %Values of isosurfaces (ranging between 0 and 1)

    PLOT.isoValSpecial          =   [];                 %Plots nicely colored isosurface
    PLOT.isoColorSpecial        =   [1.0 0.5 0.4];      %e.g., blue: [0.5 0.6 1.0] or red: [1.0 0.5 0.4]
    PLOT.opacitySpecial         =   [1];                %Opacity of the special isosurface(s)
    
    %% 3-D CARTESIAN MODEL: PROCESS SLICE AS 2-D  -------------------------
    SWITCH.ThreeToTwoD          =   logical(0);     %Switches to 2-D mode for a specified y-slice (see f_ReadStagYY)

    %% 3-D CARTESIAN MODEL: SLICING  --------------------------------------
    PLOT.Slices                 =   logical(0);     %Switch to add slices
        PLOT.sliceX             =   [50];           %In plotting values e.g.,[0:0.01:1]
        PLOT.sliceY             =   [50];
    	PLOT.sliceZ             =   [2850];
        PLOT.SliceOpacity       =   0.5;

    PLOT.SliceSpecial           =   logical(0);
        PLOT.sliceXS            =   [];             %in plotting values
        PLOT.sliceYS            =   [];
        PLOT.sliceZS            =   [];
    
   	PLOT.sliceAddQuiver         =   logical(0);
        PLOT.sliceXV            =   [];             %in plotting values for Quiver
        PLOT.sliceYV            =   [];
        PLOT.sliceZV            =   [];
        
    PLOT.SliceAddTopo           =   logical(0);
%                                 TOPO.surf=logical(1); TOPO.cmb=logical(0);
%                                 TOPO.true  =logical(0);
%                                 TOPO.smooth=logical(1);
%                                 SWITCH.ColormapAllT = {'topo5'  'noflip'};
%                                 SWITCH.constantColorbarT = logical(1);

    PLOT.sliceXT                =   [];             %in plotting values
    PLOT.sliceYT                =   [];
    PLOT.sliceZT                =   [0.01  1.0];    %[surf cmb]
    PLOT.opacity_plotT          =   1.0;
    %% --------------------------------------------------------------------   
    
    %% YINYANG MODELS  ----------------------------------------------------
    SWITCH.yyPlotMode           =   'map';          %'isosurface' or 'map' or 'mapOnSphere'
    
    %% YINYANG MODEL: ISOSURFACE ------------------------------------------ %NOT FULLY IMPLEMENTED YET
    % see 3-D Cartesian Model plot options
    
    %% YINYANG MODEL: MAPS  -----------------------------------------------
    PLOT.yyZlevelV              =   [0.02 0.69 0.91];  	%Non-dim whole-domain depth: 0.94 (plates); ~0.85 (slabs)
    PLOT.yyNumSlices            =   1;                  %Number of horizontal slices per subplot
    
    SWITCH.yyFlipVertical       =   logical(1);     	%Flips map vertically to be according to spherical plot

    PLOT.yyOriginLong           =   -180;            	%Origin longitude in degrees
    SWITCH.yyLimitMap           =   logical(0);
      	PLOT.yyLimitN           =   45;
      	PLOT.yyLimitS           =   -45;
     	PLOT.yyLimitW           =   -45;
      	PLOT.yyLimitE         	=   45;
    
    SWITCH.yyLabelMeridian      =   'off';
    SWITCH.yyLabelParallel      =   'on';       
    PLOT.yyPlabelTicks          =   [-75 -60 -45 -30 -15 0 15 30 45 60];
    PLOT.yyPlabelLocation       =   'west';  %east | {west} | prime | scalar longitude

    PLOT.yyAnnTitle             =   logical(1);
    PLOT.yyAnnDepth             =   logical(1);
    PLOT.yyDepthGraph           =   logical(0);         %Adds a small depth graph - DEOSN'T WORK YET WITH AXIS STUFF IN STAGPLOT........

    SWITCH.yyOnlyTopQuiver      = 	logical(0);         %Puts arrows only on uppermost slice
    SWITCH.yyOnlyBotQuiver      =  	logical(0);         %Puts arrows only on lowermost slice
    %% --------------------------------------------------------------------    
    
    %% POST-PROCESSING
    SWITCH.PlateDiagnostics         =   logical(0);
    SWITCH.MantleDiagnostics        =   logical(0);
    PLOT.indicateTrench             =   logical(0);     PLOT.IndicationDepth = 0; %-0.0263; %[nd]
    PLOT.indicateRidge              =   logical(0);     PLOT.IndicationSize = 20;
    PLOT.indicateUpperPlate         =   logical(0);     PLOT.PlateIndSize = 6; %upper plate indication
    PLOT.indicateLowerPlate         =   logical(0);     %Lower-plate indication
    PLOT.indicatePlateCore          =   logical(0);     %Depth of the plate core
    PLOT.indicateShallowSlabDip   	=   logical(0);
    PLOT.indicateSlabTip            =   logical(0);
    PLOT.indicateSlabTipDip         =   logical(0);
    PLOT.indicatePlateFit           =   logical(0);     %Fit of the subducting plate
    PLOT.indicateTrenchDepth        =   logical(0);     %Regional topography feature
    PLOT.indicateForeBulgeHeight	=   logical(0);     %Regional topography feature
    PLOT.indicateIslandArcHeight  	=   logical(0);     %Regional topography feature
    PLOT.indicateBackArcBasinDepth 	=   logical(0);     %Regional topography feature
    PLOT.indicateBackArcBasinExtent	=   logical(0);     %Regional topography feature
    PLOT.indicateInundation         =   logical(0);     %Regional topography feature
    PLOT.indicateUPtilt             =   logical(0);     %Upper-plate tilt
    PLOT.indicateBending            =   logical(0);
    PLOT.indicatePlumes             =   logical(0);     %Contour outline of mantle plumes 
       	PLOT.pHot                   =   0.4;            %Fraction of difference Tmax-Tmean and Tmin-Tmean, resp. (T_threshold = Tmean+p_hot*(Tmax-Tmean) )
        PLOT.pCold                  =   0.4;            %Fraction of difference Tmax-Tmean and Tmin-Tmean, resp. (T_threshold = Tmean+p_hot*(Tmax-Tmean) )
      	PLOT.plumeColorHot          =   'r'; 
        PLOT.plumeColorCold         =   'b';
     	PLOT.onlyBotPlumes          =   logical(0);     %For YY-map: only plot in bottom slice
      	SWITCH.savePlumes           =   logical(0);     %Saves plume output
    PLOT.indicateCONTINENTposition  =   logical(0);     %Continent position
    PLOT.indicateLLSVPposition      =   logical(0);     %LLSVP position
    
    %% PERFORMANCE
    SWITCH.ReduceSize               =   logical(0);   	%Reduces matrix size for faster plotting
    SWITCH.EfficientTracerReading	=   logical(1);     %Reads only a selection of tracers
    SWITCH.enforceFixingTimedatFile =   logical(0);   	%Enforces attempt to automatically fix corrupt time.dat files
    
    %% DATA ALTERATION
    SWITCH.flipDataHorizontally =   logical(0);     %Flips data horizontally
    SWITCH.shiftDataX           =   0;              %Horizontal x-shift of data as fraction of model width (0-1)
    SWITCH.shiftDataY           =   0;              %Horizontal y-shift of data as fraction of model width (0-1)
    
    %% SPECIAL PLOTTING BEHAVIOUR
    SWITCH.PlotInPlot           =   logical(0); 	%For every field plots data from all files onto same plot  *not implemented for all fields yet*
    SWITCH.plotGraphVsTime      =   logical(0);     %Plots in 2-D, a selected graph (x) vs. time (y)  *only for topography, yet*
    SWITCH.plotDifference       =   logical(0);     %Plots difference between two fields
        
    %% FIELDS TO PLOT
    PLOT.Temperature            =   logical(1);
    PLOT.RegionalResidualT      =   logical(0);
    PLOT.HorizResidualT         =   logical(0);
    PLOT.HorizBandResidualT     =   logical(0);
    PLOT.GlobalResidualT        =   logical(0);
    PLOT.Composition            =   logical(0);  SWITCH.TickLabelsComposition = {'Crust','Mantle','Continent','Air','Primordial'};
    PLOT.Density                =   logical(0);
    PLOT.Viscosity              =   logical(0);  IN.eta_log = logical(1);
    
    PLOT.Velocity               =   logical(0);
    PLOT.VelocityX              =   logical(0);
    PLOT.VelocityZ              =   logical(0);
    PLOT.UpDownWelling          =   logical(0);  %Plots active and passive up- and down-welling
    PLOT.Pressure               =   logical(0);  %Total pressure
    PLOT.HeatFlux               =   logical(0);  SWITCH.plotFIELDtrue=logical(1); SWITCH.plotFIELDsmooth=logical(0); SWITCH.plot3Dfield2D=logical(0);
    PLOT.StrainRate             =   logical(0);  IN.edot_log=logical(0);
    PLOT.Stress                 =   logical(0);  IN.str_log=logical(0);
    PLOT.StressX                =   logical(0);  %Principle horizontal (x) stress component
    PLOT.StressZ                =   logical(0);  %Principle radial stress component
    PLOT.NormalStress           =   logical(0);  %Sigma_zz actually
    PLOT.ViscousDissipation  	=   logical(0);  %Sigma_ij*edot_ij

    PLOT.Topography           	=   logical(0);  TOPO.ascii=logical(0); TOPO.cmb=logical(0); TOPO.surf=logical(1);
      	TOPO.color              =   [0.81 0.81 0.81]; 
        TOPO.opacity            =   0.9; 
      	TOPO.areaFill           =   'down';         %Filled area in topography plot: 'none', 'down', 'updown'
      	TOPO.water              =   true;           %Fills negative levels in topography plot with "water", i.e., a blue area
      	SWITCH.plotTOPOtrue     =   logical(0);
     	SWITCH.plotTOPOsmooth   =   logical(1);
        TOPO.indicateComponents = 	logical(0); 
        TOPO.difference         =   'comp';         %'comp' or 'diff'
      	SWITCH.plot3Dtopo2D     =   logical(0);     %For 3-D plots
        TOPO.exagFactor         =   1;              %For 3-D plots
    PLOT.IsostaticTopography  	=   logical(0);
        TOPO.normIsoTopoToZero  =   logical(1);     %Set mean of isostatic component to zero due to volume conservation in an incompressible mantle
    PLOT.ResidualTopography  	=   logical(0);
    PLOT.DynamicTopography    	=   logical(0);  
        TOPO.dynTOPOsmooth      =   logical(1);
    PLOT.TopographySelfGrav   	=   logical(0);
    PLOT.Geoid                  =   logical(0);  PLOT.GeoidCMB=logical(0); PLOT.GeoidSurf=logical(1);
    
    PLOT.DeformationMechanism 	=   logical(0);
    PLOT.Toroidal               =   logical(0);
    PLOT.Poloidal               =   logical(0);
        
    PLOT.Phase                  =   logical(0);
    PLOT.Basalt                 =   logical(0);     %Individual compositions
    PLOT.Harzburgite            =   logical(0);
    PLOT.ContinentalCrust     	=   logical(0);
    PLOT.Air                    =   logical(0);
    PLOT.Primordial             =   logical(0);
    PLOT.Water                  =   logical(0);  IN.water_log = logical(0);
    PLOT.Melt                   =   logical(0);
    PLOT.AgeSinceLastMelted     =   logical(0);

    %these are not implemented yet:
    PLOT.DynamicPressure      	=   logical(0); %...
    PLOT.MeltFraction         	=   logical(0); %...
    
    %% GRAPHS TO PLOT
    PLOT.PlateVelocity         	=   logical(0);     %Graph indicating horizontal plate velocity
    PLOT.SurfaceAge             =   logical(0);
    PLOT.CrustalThickness     	=   logical(0);
    PLOT.PlateBaseTopography   	=   logical(0);     %For LAB topography; uses same parameter input as for lithoThickness
    PLOT.StreamGraph         	=   logical(0);     %Graph of any field data versus time along a streamline
        PLOT.StreamGraphParameter =   {'Temperature'};  %indicate which parameter(s) to plot (use field names)
    PLOT.CustomGraph            =   logical(0);     %Plot graph of any data versus any data (that has been saved previously with StagLab)
        PLOT.CustomGraphName            =   {'custom'};             %Indicate which graph(s) to plot, use 'custom' for a custom graph specified in f_plotCustomGraph.m or any pre-specified graph there
        PLOT.CustomDataNameX            =   {'Time'};               %See SL_FieldPlot for available fields; only applied if PLOT.CustomGraphName is set to 'custom'
        PLOT.CustomDataNameY            =   {'Trench velocity', 'Lower-plate velocity', 'Convergence velocity', 'Slab sinking velocity'};	%See SL_FieldPlot for available fields; only applied if PLOT.CustomGraphName is set to 'custom'
        PLOT.CustomDataNumMin           =   1;                      %Minimum data number to specify plotted range
        PLOT.CustomDataNumMax           =   999;                    %Maximum data number to specify plotted range
        PLOT.CustomIndicateCurrentTime 	=   logical(1);             %Indicates current time step
        PLOT.CustomFieldName        	=   'Plate velocities';   	%Custom plot title
        PLOT.CustomAbsoluteValues     	=   logical(0);          	%Converts all data to absolute values
        PLOT.CustomXAxisMin           	=   0;                      %Minimum to specify x-axis extent
        PLOT.CustomXAxisMax            	=   55;                     %Maximum to specify x-axis extent
        PLOT.CustomYAxisMin           	=   -10.0;                  %Minimum to specify y-axis extent
        PLOT.CustomYAxisMax           	=   10.0;                   %Maximum to specify y-axis extent

    %% SPECIAL PLOTS
    PLOT.ParameterTable      	=   logical(0);     %Plots a table with dynamic parameters
    PLOT.PlateSketch            =   logical(0);     %Plot sketch of plates with subduction parameters indication
    PLOT.Grid                   =   logical(0);  PLOT.nCoarserGrid 	= 1; %plots grid
    PLOT.Tracers              	=   logical(0);  PLOT.TracerVariable = {'Position'}; %plots the position of the tracers or any tracer-stored variable specified
    PLOT.Streamfunction       	=   logical(0);  PLOT.streamdata = 'psi'; %psi (streamfunction) or phi (v-potential)
    PLOT.Streamline             =   logical(0);     %Separate figure for streamline
    PLOT.Quiver                 =   logical(0);     %Separate figure for velocity
    PLOT.PrincStressDirection   =   logical(0);     %Separate figure for principle stress direction
    PLOT.SurfaceFieldVariation  =   logical(0);     %Bar plot of surface field variation
        PLOT.sfvField           =   'Surface age';  %Specify field for surface variation bar plot
        PLOT.sfvDepthLevel      =   0;              %[non-dimensional], Specify depth level for multi-dimensional field data, e.g., 0:top, 1:bottom
        PLOT.sfvNumberBins      =   24;             %Number of bins across the x-axis extent
        PLOT.sfvAdditions       =   logical(1);     %Indicates mean, median and standard deviation
        PLOT.sfvXmin            =   'auto';         %Automatic ('auto') or manual (number) x-axis limit
        PLOT.sfvXmax            =   'auto';         %Automatic ('auto') or manual (number) x-axis limit
        PLOT.sfvShowOutliers    =   logical(0);     %Adds all upper-bound outliers to the last bar

    %% FIELD PLOT ADDITIONS
    IN.grid_T                   =   logical(0);     %Add grid to fields
    IN.grid_eta                 =   logical(0);
    IN.grid_str                 =   logical(0);
    IN.grid_edot                =   logical(0);
    
    IN.streamline_T             =   logical(0);     %Add streamline plot to fields
    IN.streamline_eta           =   logical(0);  PLOT.streamNumStartPoints = 10; PLOT.streamLength = 300; %streamline length in [#cells]
    IN.streamline_str           =   logical(0);  PLOT.streamStartline = [0.4 0.4; 0.0 0.0; 0.1 0.9]; %[x1 x2; y1 y2; z1 z2]; vary <0-1>
    IN.streamline_v             =   logical(0);
    IN.streamline_vh            =   logical(0);
    IN.streamline_vr            =   logical(0);
    
    IN.streamfun_T              =   logical(0);     %Add streamfunction plot to fields
    IN.streamfun_eta            =   logical(0);  IN.streamNumContours = 17;
    IN.streamfun_str            =   logical(0);  PLOT.streamColor   = [0.81 0.81 0.81]; PLOT.streamColoured = 'false'; %'false','<colormapName>' (needs 'freezeColors.m')
    IN.streamfun_rho            =   logical(0);  PLOT.streamWidth   = 0.2;
    IN.streamfun_v              =   logical(0);
    IN.streamfun_vx           	=   logical(0);
    IN.streamfun_vr          	=   logical(0);
    IN.streamfun_edot           =   logical(0);
    IN.streamfun_udw            =   logical(0);

    IN.quiver_T                 =   logical(0);     %Add quiver plot to fields
    IN.quiver_eta               =   logical(0);  SWITCH.quiverNumArrows = 'auto'; %'auto' or 'manual': reduces number of arrows
    IN.quiver_str               =   logical(0);  IN.quiverNumDel    = 14;
    IN.quiver_psi               =   logical(0);  SWITCH.quiverScale = logical(1); IN.quiverScale = 17;
    IN.quiver_rho               =   logical(0);  PLOT.quiverColor   = [0.81 0.81 0.81];
    IN.quiver_nstr              =   logical(0);  PLOT.quiverWidth   = 0.25;
    IN.quiver_edot              =   logical(0);  PLOT.quiverArrowHead = 'on'; %'on' or 'off'
    IN.quiver_v                 =   logical(0);
    IN.quiver_vx              	=   logical(0);
    IN.quiver_vr              	=   logical(0);
    
    IN.defmech_T                =   logical(0);     %1: diffusion creep, 2: dislocation creep, 3: plasticity, 4: visc. cutoff
    IN.defmech_eta              =   logical(0);  IN.defmechContours	= [2 3 4];  %1.3: 1/3, 1.1: 1/10 of weakening from dislocation creep
    IN.defmech_str              =   logical(0);  PLOT.defmechColor	= [0.9 0.1 0.1; 0.7 0.0 0.0; 0.9 0.1 0.1; 0.1 0.9 0.1; 0.1 0.1 0.9];
    IN.defmech_edot             =   logical(0);
    
    IN.lithoThickness_T         =   logical(0);  PLOT.LTfieldContour= 'Temperature';
    IN.lithoThickness_eta       =   logical(0);  IN.lithoTisovalue 	= 1600; 
    IN.lithoThickness_nstr      =   logical(0);  IN.lithozmax     	= 500;
    
    IN.fieldContour_C           =   logical(0);
    IN.fieldContour_T           =   logical(0);  PLOT.fContourField	= 'Temperature';
    IN.fieldContour_eta         =   logical(0);  IN.fContourIsovalue= 1600;
    IN.fieldContour_nstr        =   logical(0);  PLOT.fieldContourColor = [0.81 0.81 0.81];
    IN.fieldContour_v           =   logical(0);  PLOT.fieldContourSize = 0.5;
    IN.fieldContour_vx          =   logical(0);
    IN.fieldContour_vr          =   logical(0);
    IN.fieldContour_str         =   logical(0);
    IN.fieldContour_edot        =   logical(0);
    IN.fieldContour_vdiss       =   logical(0);
    IN.fieldContour_psi         =   logical(0);
    IN.fieldContour_udw         =   logical(0);
    
    IN.princStress_T            =   logical(0);     %Add quiver plot to fields
    IN.princStress_eta          =   logical(0);  SWITCH.princStressNumArrows = 'auto'; %'auto' or 'manual': reduces number of arrows
    IN.princStress_str          =   logical(0);  IN.princStressNumDel    = 5;
    IN.princStress_psi          =   logical(0);  IN.princStressScale     = 2;
    IN.princStress_rho          =   logical(0);  PLOT.princStressColor   = [0.91 0.61 0.61];
    IN.princStress_nstr         =   logical(0);  PLOT.princStressWidth   = 0.25;
    IN.princStress_edot         =   logical(0);
    
    PLOT.Plume                  =   logical(0);     %!!! not implemented yet - single plot of plumes..................
    IN.plume_T                  =   logical(0);     %Add plume plot addition to fields
    IN.plume_eta                =   logical(0);
    IN.plume_str                =   logical(0);
    IN.plume_rho                =   logical(0);
    IN.plume_v                  =   logical(0);
    IN.plume_edot               =   logical(0);
    
    TOPO.field                  =   logical(0);     %Add topography to field plots
        TOPO.line               =   logical(1); 
        TOPO.area               =   logical(1); 
        PLOT.hatchAir           =   logical(1);
      	TOPO.lineColor          =   [0.2 0.2 0.2]; 
        TOPO.lineWidth          =   0.25;

        
    %% SAVING FIGURE
    SAVE.Figure               	=   logical(0);     %Saves plot to directory
        SAVE.png                =   logical(1);     SAVE.pngResolution = '-m4.0';	
        SAVE.jpg                =   logical(0); 
        SAVE.eps                =   logical(0); 
        SAVE.pdf                =   logical(0); 
        SAVE.FigureCropping   	=   logical(1);     %Crops figure to have tight borders around actual plot
        SAVE.TransparentBackground = logical(0);    %Removes the figure background (only available for .png format)
        SAVE.overwrite          =   logical(0);     %Always overwrites previous versions with the same file name
        SAVE.neverOverwrite     =   logical(0);     %Always keeps previous versions with the same file name
        SAVE.writeDirectory     =   'auto';         %'auto': save to standard folder (.../+im/...); otherwise use e.g., '/work/stagyy/'
        
    %% SAVING MOVIE
    SAVE.Movie                  =   logical(0);     %Saves movie to directory
        SAVE.avi                =   logical(0);
        SAVE.mj2                =   logical(0);
        SAVE.mp4                =   logical(1);
        SAVE.m4v                =   logical(0);
        SAVE.movFrameRate       =   10;            	%10 (default) | frames per second
        SAVE.movQuality         =   75;             %75 (default) | integer in the range [0,100]
        
    %% SAVING DATA
    SAVE.Field                  =   logical(0);     SAVE.mat = logical(1); SAVE.dat = logical(0);  %Saves plotted 2-D/3-D field according to f_saveField
    SAVE.GeodynamicDiagnostics 	=   logical(0);     %Saves geodynamic data (see list of data in SL_FieldPlot.m)
    SAVE.Tracers                =   logical(0);     %Saves processed tracer data
    SAVE.timeDat                =   logical(0);     %Saves time [nd, seconds, years]
    TOPO.saveData               =   logical(0);     %Saves topo: [x, topo, topoSmooth]
    SAVE.plateVelocityData      =   logical(0);     %Saves plate velocity vector data: []
    SWITCH.trackC               =   logical(0);     %Tracks and saves specified point of the composition field
    SAVE.subductionZoneData     =   logical(0);     %OLD (use SAVE.GeodynamicDiagnostics instead!): saves subduction zone data: [x-trench, subduction-polarity, trench/UP-vel, LP-vel, Convergence-vel,...]
end

