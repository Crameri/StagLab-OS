
%%                                             StagLab PARAMETER FILE: 3-D
%                              (for all available options see f_Defaults.m)
clear; clf;
SAVE.StartNumber       =   14;      %Set SAVE.StartNumber==SAVE.EndNumber to write only one file
SAVE.StepNumber        =   1;
SAVE.EndNumber         =   14;      %Set big in order to continue to last file

for iLoop=SAVE.StartNumber:SAVE.StepNumber:SAVE.EndNumber
clearvars -except iLoop SAVE; if ~exist('f_Defaults','file'); f_INSTALL; end
    [IN,PLOT,SWITCH,TOPO,STYLE,SAVE] = f_Defaults(SAVE); %get default variables

    %% DEBUGGING
    SWITCH.sendErrorLog         =   logical(0);                         %Sends an error log file to the developer for debugging purposes

    %% INPUT FILE(S)
    IN.Name                     =   { 'test3d'    	};        	%Filename (it's size controls the # of files processed)
    IN.Number                   =   [   iLoop       ];          %Filenumber (use iLoop to loop through multiple files)
    IN.Parameter                =   [  	11        	];        	%Parameter dimensionalisation (see f_Dimensions)
    IN.Folder                   =   {'/work/stagyy/'};          %File directory
    %PLOT.titleString            =   {'XXX'          'XXX'};   	%Adds manual titles if uncommented

    %% StagLab MODES
    SWITCH.DimensionalMode   	=   logical(0);   	%Converts to dimensional output as set in f_Dimensions.m
    STYLE.ColorMode             =   'light';        %'light' or 'dark' background
    SWITCH.AnalysisMode         =   logical(0); 	%Creates figures for analysis purposes - not recommended for publication!
    SWITCH.QuickMode         	=   logical(0);     %Adjust parfile switches for quick runtime to e.g. save data efficiently

    %% SUBPLOT LAYOUT
    %PLOT.Layout                 =   [2 2];     	%Manual subplot layout; only takes effect if uncommented; [z x]
    SWITCH.ReverseLayout       	=   false;          %e.g., plot 3x2 instead of 2x3
    SWITCH.TimeEvolutionMode    =   logical(0);     %Adds arrows between subplots indicating time evolution
        SWITCH.timeDirection    =   'topbot';     	%'topbot', 'leftright'

    %% PLOT STYLING
    SWITCH.AxesLimit            =   logical(0);     
        SWITCH.AxesLimitValues  =   [500 1300 -5 200]; %Axes limit for big plots
   
    %% ANNOTATIONS
    SWITCH.Annotation           =   logical(1);     %Subplot annotations with e.g., 'a', 'b', 'c'
    SWITCH.Title                =   logical(1);     
        PLOT.titleOnlyTime      =   logical(0);   	%For a title showing only time
        PLOT.titleNoTime        =   logical(0);     %For a title showing no time
    
    %% PLOT ADDITIONS
    SWITCH.Magnifier          	=   logical(0);     %Adds an additional plot for a close up
        PLOT.magnifierExtent  	=   [2500 3500 -50 500]; %[xmin xmax ymin ymax] in actual plotting dimensions
    PLOT.indicateDepth          =   logical(0);     %Adds horizontal depth indicator lines
        PLOT.indicateDepthLevels=   [410, 660];     %Specify one or multiple depth levels

    %% COLORBAR
    SWITCH.ConstantColorbar     =   logical(0);     %Fix colorbar limits to preset values
    % use these to set the min/max colorbar values of each field (see f_Defaults.m for more options):
    %IN.topo_min                 =   -0.3;
    %IN.topo_max                 =   0.3;
    %IN.C_min                    =   0;
    %IN.C_max                    =   4;
    %IN.T_min                    =   300;
    %IN.T_max                    =   1150;
    % ...and each individual subplot use brackets:
    %IN.edot_max                 =   [3e-20  1.05e-18];
    %IN.edot_min                 =   [0      0       ];

	%% COLOUR MAP
    SWITCH.DiscreteColormap     =   logical(0);
    PLOT.NumberColormapColors 	=   NaN;            %Give number of discrete colours or set 'NaN' to have it automatically adjusted to number of tick labels
    SWITCH.cbarLimitIndicTop    =   [0]; SWITCH.cbarLimitColorTop = [0.81 0.81 0.81]; %Adds specific color to top/bot to indicate exceeding values
    SWITCH.cbarLimitIndicBot    =   [0]; SWITCH.cbarLimitColorBot = [0.81 0.81 0.81]; %(only one number sets it for all plots, otherwise specify for every subplot using [0 1 0 ...])
    
    %% POST-PROCESSING
    SWITCH.PlateDiagnostics         =   logical(0);
    SWITCH.MantleDiagnostics        =   logical(0);
    PLOT.indicateTrench             =   logical(0);
    PLOT.indicateRidge              =   logical(0);
    PLOT.indicateUpperPlate         =   logical(0);
    PLOT.indicateLowerPlate         =   logical(0);
    PLOT.indicatePlateCore          =   logical(0);     %Depth of the plate core
    PLOT.indicateShallowSlabDip            =   logical(0);
    PLOT.indicateSlabTip            =   logical(0);
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
    PLOT.UpDownWelling          =   logical(0);     %Plots active and passive up- and down-welling
    PLOT.Pressure               =   logical(0);     %Total pressure
    PLOT.HeatFlux               =   logical(0);  SWITCH.plotFIELDtrue=logical(1); SWITCH.plotFIELDsmooth=logical(0);
    PLOT.StrainRate             =   logical(0);  IN.edot_log=logical(0);
    PLOT.Stress                 =   logical(0);  IN.str_log=logical(0);
    PLOT.ViscousDissipation  	=   logical(0);     %Sigma_ij*edot_ij

    PLOT.Topography           	=   logical(0);  TOPO.cmb=logical(0); TOPO.surf=logical(1);
        SWITCH.plotTOPOtrue     =   logical(0);
     	SWITCH.plotTOPOsmooth   =   logical(1);
        TOPO.indicateComponents = 	logical(0);     %Indication for isostatic and residual topography 
    PLOT.IsostaticTopography  	=   logical(0);
    PLOT.ResidualTopography  	=   logical(0);
    PLOT.Geoid                  =   logical(0);	 PLOT.GeoidCMB=logical(0); PLOT.GeoidSurf=logical(1);
    
    PLOT.DeformationMechanism 	=   logical(0);
        
    PLOT.Basalt                 =   logical(0);     %Individual compositions
    PLOT.Harzburgite            =   logical(0);
    PLOT.ContinentalCrust     	=   logical(0);
    PLOT.Air                    =   logical(0);
    PLOT.Primordial             =   logical(0);
    PLOT.Water                  =   logical(0);
    PLOT.Melt                   =   logical(0);
    PLOT.AgeSinceLastMelted     =   logical(0);
    
    %% GRAPHS TO PLOT
    PLOT.PlateVelocity         	=   logical(0);     %Graph indicating horizontal plate velocity
    PLOT.SurfaceAge             =   logical(0);
    PLOT.CrustalThickness     	=   logical(0);
    PLOT.PlateBaseTopography   	=   logical(0);     %For LAB topography; uses same parameter input as for lithoThickness
    PLOT.StreamGraph         	=   logical(0);     %Graph of any field data versus time along a streamline
        PLOT.StreamGraphParameter =   {'Temperature'};  %indicate which parameter(s) to plot (use field names)
    PLOT.CustomGraph            =   logical(0);     %Plot graph of any data versus any data (that has been saved previously with StagLab)
        PLOT.CustomGraphName    =   {'custom'};     %Indicate which graph(s) to plot, use 'custom' for a custom graph specified in f_plotCustomGraph.m or any pre-specified graph there
    
    %% SPECIAL PLOTS
    PLOT.ParameterTable      	=   logical(0);     %Plots a table with dynamic parameters
    PLOT.PlateSketch            =   logical(0);     %Plot sketch of plates with subduction parameters indication
    PLOT.Grid                   =   logical(0);  PLOT.nCoarserGrid 	= 1; %plots grid
    PLOT.Tracers              	=   logical(0);  PLOT.TracerVariable = {'Position'}; %plots the position of the tracers or any tracer-stored variable specified
    PLOT.Streamfunction       	=   logical(0);
    PLOT.SurfaceFieldVariation  =   logical(0);     %Bar plot of surface field variation
        PLOT.sfvField           =   'Surface age';  %Specify field for surface variation bar plot
        
    %% FIELD PLOT ADDITIONS
    IN.grid_T                   =   logical(0);     %Add grid to fields
    IN.grid_eta                 =   logical(0);
    IN.grid_str                 =   logical(0);
    IN.grid_edot                =   logical(0);
    
    IN.streamline_T             =   logical(0);     %Add streamline plot to fields
    IN.streamline_eta           =   logical(0);  PLOT.streamNumStartPoints = 10; PLOT.streamLength = 300; %streamline length in [#cells]
    IN.streamline_str           =   logical(0);  PLOT.streamStartline = [0.4 0.4; 0.0 0.0; 0.1 0.9]; %[x1 x2; y1 y2; z1 z2]; vary <0-1>
    
    IN.streamfun_T              =   logical(0);     %Add streamfunction plot to fields
    IN.streamfun_eta            =   logical(0);  IN.streamNumContours = 17;
    IN.streamfun_str            =   logical(0);
    IN.streamfun_rho            =   logical(0);
    IN.streamfun_v              =   logical(0);
    IN.streamfun_edot           =   logical(0);

    IN.quiver_T                 =   logical(0);     %Add quiver plot to fields
    IN.quiver_eta               =   logical(0);  SWITCH.quiverNumArrows = 'auto'; %'auto' or 'manual': reduces number of arrows
    IN.quiver_str               =   logical(0);  IN.quiverNumDel    = 14; %number of arrows to delete
    IN.quiver_psi               =   logical(0);  IN.quiverScale     = 17;
    IN.quiver_rho               =   logical(0);
    IN.quiver_nstr              =   logical(0);
    IN.quiver_edot              =   logical(0);
    
    IN.defmech_T                =   logical(0);     %1: diffusion creep, 2: dislocation creep, 3: plasticity, 4: visc. cutoff
    IN.defmech_eta              =   logical(0);
    IN.defmech_str              =   logical(0);
    IN.defmech_edot             =   logical(0);
    
    IN.lithoThickness_T         =   logical(0);  PLOT.LTfieldContour= 'Temperature';
    IN.lithoThickness_eta       =   logical(0);  IN.lithoTisovalue 	= 1600; 
    IN.lithoThickness_nstr      =   logical(0);  IN.lithozmax     	= 500;
    
    IN.fieldContour_T           =   logical(0);  PLOT.fContourField	= 'Temperature';
    IN.fieldContour_eta         =   logical(0);  IN.fContourIsovalue= 1600;
    IN.fieldContour_nstr        =   logical(0);
    IN.fieldContour_vx          =   logical(0);
    IN.fieldContour_str         =   logical(0);
    IN.fieldContour_edot        =   logical(0);
    IN.fieldContour_vdiss       =   logical(0);
    IN.fieldContour_psi         =   logical(0);
    
    IN.princStress_T            =   logical(0);     %Add quiver plot to fields
    IN.princStress_eta          =   logical(0);  SWITCH.princStressNumArrows = 'auto'; %'auto' or 'manual': reduces number of arrows
    IN.princStress_str          =   logical(0);  IN.princStressNumDel    = 5;
    IN.princStress_psi          =   logical(0);  IN.princStressScale     = 2;
    IN.princStress_rho          =   logical(0);
    IN.princStress_nstr         =   logical(0);
    IN.princStress_edot         =   logical(0);
    
    IN.plume_T                  =   logical(0);     %Add plume plot addition to fields
    IN.plume_eta                =   logical(0);
    IN.plume_str                =   logical(0);
    IN.plume_rho                =   logical(0);
    IN.plume_v                  =   logical(0);
    IN.plume_edot               =   logical(0);
    
    TOPO.field                  =   logical(0);     %Add topography to field plots
        TOPO.line               =   logical(1); 
        
    %% SAVING FIGURE
    SAVE.Figure               	=   logical(0);     %Saves plot to directory

    %% SAVING MOVIE
    SAVE.Movie                  =   logical(0);     %Saves movie to directory

    %% SAVING DATA
    SAVE.Field                  =   logical(0);     %Saves plotted 2-D/3-D field according to f_saveField
    SAVE.GeodynamicDiagnostics 	=   logical(0);     %Saves geodynamic data (see list of data in SL_FieldPlot.m)
    TOPO.saveData               =   logical(0);     %Saves topo: [x, topo, topoSmooth]
    
    %% RUN MAIN ROUTINE >>>
    [SAVE] = SL_FieldPlot(IN,PLOT,SWITCH,TOPO,STYLE,SAVE); %run main routine
    if SAVE.LastFile; break; end %last file reached
end
    
    
