
%%                                                      SL_FieldPlot 16.000
%
%                                               Plots parameter-field data
%    . calls f_AIo
%    . calls f_Connectivity
%    . calls f_Defaults
%    . calls f_DefaultsColorbar
%    . calls f_DesignAnnotations
%    . calls f_DesignBackground
%    . calls f_DesignColourmap
%    . calls f_DesignFigure
%    . calls f_DesignLayout
%    . calls f_DesignLayoutPosition
%    . calls f_DesignLegend
%    . calls f_DesignPlotRelation
%    . calls f_DesignSimplify
%    . calls f_DesignVaria
%    . calls f_DiagnosticsMantle
%    . calls f_DiagnosticsPlate
%    . calls f_DiagnosticsSurface
%    . calls f_Dimensions
%    . calls f_EmptyPlot
%    . calls f_FileFinder
%    . calls f_Import
%    . calls f_IncludedFit
%    . calls f_makeField
%    . calls f_PlateThickness
%    . calls f_plotField
%    . calls f_plotFieldSpecial
%    . calls f_plotCustomGraph
%    . calls f_plotGraph
%    . calls f_readFluidity
%    . calls f_readOther
%    . calls f_readStagYY
%    . calls f_readStagYYhdf5
%    . calls f_readTopo3D
%    . calls f_saveData
%    . calls f_saveField
%    . calls f_saveFigure
%    . calls f_saveMovie
%    . calls f_SetupFigure
%    . calls f_SetupTopography
%    . calls f_Setup2DField
%    . calls f_Startup
%    . calls f_trackPlumes   %only for YY
%    . calls f_Varia
%    . calls f_YYtoMap
%
%    . calls hatchfill2
%    . calls MinVolEllipse
%    . calls flowfun
%    . calls cumsimp
%    . calls panel
%    . calls plotboxpos
%
%                                               17.06.2021 : Fabio Crameri

%% UPDATES StagLab >6.0

%% UPDATES StagLab 6.0
% *Beta compatibility for ASPECT code*
% *Deeper integration of StagYY's HDF5 output format*
% *Added Scientific categorical colour maps*
% *fix for usage with other system time formats*
% *fix for updated StagLab repository*

%% UPDATES StagLab 5.0
% *introducing continent diagnostics*
% *introducing LLSVP diagnostics*
% *introducing slab-tip diagnostics*
% *introducing horizontal mantle flow diagnostics*
% *introducing panel*
% *introducing transparent figure background*
% *introducing automated error logging*
% *introducing categorical Scientific Colour Maps*
% *introducing journal-specific plot design*
% *additional StagLab-data output*
% *improvements and extensions to YinYang mode*
% *additional parameter field additions*
% *option to make difference plot with vector data*
% *less interruptive updating of old parfiles*
% *improved Windows compatibility*
% *special character fix for Windows*
% *improved code design*
% *stability and speed improvements for loop mode*
% *bug fixes*

%% UPDATES StagLab 4.0
% *re-introducing multi-subduction-zone tracking*
% *introducing plume-mobility diagnostics*
% *improved Windows compatibility*
% *improved compatibility with latest StagYY version*
% *improvements to yinyang mode*
% *additional parameter fields*
% *extended suite of scientific colour-maps*
% *analysis mode for SL_RadialProfile and SL_TimeGraph*
% *flexibility extensions to SL_RadialProfile and SL_TimeGraph*
% *automatic fixing of corrupt time.dat files*
% *stability improvements*
% *bug fixes*

%% UPDATES StagLab 3.0
% *introducing automated installation and testing*
% *introducing 2-D mode for 3-D*
% *introducing analysis mode*
% *introducing tracer plot*
% *introducing surface-variation histogram plot*
% *introducing topography diagnostics*
% *introducing perceptually-uniform colour schemes*
% *option to discretize colour maps*
% *option to set default figure position on screen*
% *option to shift data horizontally*
% *support for partial spherical2D geometry*
% *magnifier support for spherical2D geometry*
% *additional parameter fields*
% *additional plate diagnostics*
% *refined visual design*
% *improved file finder*
% *improved code design*
% *improved speed*
% *major improvements to SL_RadialProfile and SL_TimeGraph*
% *bug fixes*

%% UPDATES StagLab 2.0
% *introducing mantle-dynamics diagnostics*
% *introducing tectonic diagnostics*
% *introducing topography components (isostatic,residual)*
% *introducing plot for up- and downwelling*
% *introducing parameter table*
% *introducing plot-in-plot mode*
% *introducing movies*
% *introducing fAIo*
% *more parameter fields added*
% *automatic detection of side-boundary v-condition*
% *less-disruptive error handling*
% *cleaner plot design and layout*
% *improved colormaps*
% *improved file finder*
% *improved display output*
% *improved stability of design-routines*
% *improved saving and plotting of tectonic data*
% *improvements towards convertibility to other geodynamic codes*
% *bug fixes*

%% UPDATES StagLab 1.0
% *introducing StagLab*
% *now supports YinYang horizontal maps*
% *hot and cold plume tracking*
% *3-D Cartesian plate boundary tracking*
% *code speed optimisations: deriving lithosphere thickness*
% *option added to plot horizontal residual temperature*
% *option added to plot heat flux*
% *option added to plot temporal evolution of tectonic parameters*
% *option added to save figure to specific directory*
% *improved code design*
% *improved user friendliness*
% *bug fix for missing axis labelling with odd number of subplots*
% *bug fixed that led to empty plate sketch plot*
% *bug fixes and updates to the dimensionalisation*

function [SAVE] = SL_FieldPlot(IN,PLOT,SWITCH,TOPO,STYLE,SAVE)
if exist('IN','var') %if called from a parfile
    %% DEFAULT VARIABLES
    %...see f_Defaults.m
else
    clear; clf; close all;
    %% INPUT (...use parameter file instead)
    SAVE.StartNumber = 1; SAVE.EndNumber = 1;
    [IN,PLOT,SWITCH,TOPO,STYLE,SAVE] = f_Defaults(SAVE); %get default variables
    IN.Name     	=   {  'testDIM'  'test' };	% filename
    IN.Number       =   [   2           3 	];
    IN.Parameter    =   [   11          11 	];
    IN.Folder       =   {   '~/work/plotTest/' };
end

%% STARTUP PROCEDURE
SAVE.app = 'SL_FieldPlot'; SAVE.appVersion = 16.0; %only 3 digits to prevent unecessary parfile updates
if ~isfield(SAVE,'count'); SAVE.count = 0; end
DESVARIA.Task    = 'set special characters';
[~,DESVARIA,STYLE,~] = f_DesignVaria(DESVARIA,[],[],[],[],STYLE,[]);
[fAIo,SAVE,SWITCH] = f_Startup(STYLE,SAVE,SWITCH);
if ~strcmp(fAIo.Status,'all fine!'); return; end
if SAVE.count==1; close all; end %only once during entire execution

try

    %% SPECIAL SETUPS
    if SWITCH.ModelSketchMode %optimise plots for simple model-setup representation
        SWITCH.AnalysisMode         = true;
        SWITCH.AxesEqual            = false;
        SWITCH.AxesLimit            = true;
        SWITCH.AxesLimitValues      = [0 6080 -50 660]; %might be adjusted
        SWITCH.ConstantColorbar     = false;
        SWITCH.DiscreteColormap     = true;
        PLOT.NumberColormapColors 	= NaN;
    end
    if SWITCH.AnalysisMode %optimise plots for accurate analysis
        if ~SWITCH.ModelSketchMode
            SWITCH.PlateDiagnostics     = true;
            PLOT.indicateTrench         = true;
            SWITCH.MantleDiagnostics    = true;
            SWITCH.Texting              = true;
        end
        SWITCH.PlotDesign           = true;
        STYLE.AllFontName           = 'Helvetica';
        STYLE.keyFontSize           = 9;
        TOPO.area                   = false;
        SWITCH.SimplifyPlots      	= true;
        SWITCH.BackgroundDesign    	= false;
        SWITCH.Annotation           = false;
        SWITCH.spherical2DCenterText= true;
        SWITCH.onlyXaxis            = false;
        STYLE.removeAxisRuler    	= false;
        STYLE.keyColor              = [0 0 0];
        SWITCH.ConstantColorbar   	= false;
        SWITCH.GridAlwaysOn         = true;
        STYLE.MinorAxisTicks      	= true;
    end
    
    DESVARIA.Task = 'set plot design mode'; %set journal design mode
    [PLOT,DESVARIA,STYLE,~] = f_DesignVaria(DESVARIA,[],[],[],PLOT,STYLE,[]);
    
    if SWITCH.QuickMode %optimise parfile input for quick running
        fields = fieldnames(PLOT);            %switch off all plotting fields
        for fn=fields'
            if islogical(PLOT.(fn{1})) || ~isempty(strfind(fn{1},'_T'))
                PLOT.(fn{1}) = false;
                if SWITCH.Verbose; disp(fn{1}); end
            end
        end
        PLOT.Temperature            = true;     %plot only temperature
        SWITCH.PlateDiagnostics     = true;
        PLOT.indicateTrench         = true;
        SWITCH.MantleDiagnostics    = true;
        SWITCH.PlotDesign           = false;    %switch off additional things
        SWITCH.SimplifyPlots      	= false;
        SWITCH.Annotation           = false;
        SWITCH.Colorbar             = false;
        SWITCH.Magnifier          	= false;
        PLOT.numContours            = 5;
        SWITCH.BackgroundGuides   	= false;
        SAVE.Figure                	= false;
        SAVE.Movie                  = false;
        TOPO.area                   = false;
        if ~TOPO.saveData; TOPO.field = false; end
    end
    if SWITCH.sendErrorLog
        SWITCH.Verbose              = true;
    end
    
    %% VARIABLE CONVERSIONS
    SAVE.BulletString       = STYLE.SCHAR.hugeBulletLight;
    SAVE.PointingString     = STYLE.SCHAR.downrightArrow;
    PLOT.gridAddition     	= [IN.grid_T IN.grid_eta IN.grid_str IN.grid_edot];
    PLOT.quiverM            = [PLOT.Quiver IN.quiver_T IN.quiver_eta IN.quiver_str IN.quiverNumDel IN.quiverScale IN.quiver_psi IN.quiver_rho IN.quiver_nstr IN.quiver_edot ...
                                IN.quiver_v IN.quiver_vx IN.quiver_vr];
    PLOT.princStressM    	= [PLOT.PrincStressDirection IN.princStressNumDel IN.princStressScale IN.princStress_T IN.princStress_eta IN.princStress_str IN.princStress_psi ...
                                IN.princStress_rho IN.princStress_nstr IN.princStress_edot];
    PLOT.streamfun          = [PLOT.Streamfunction IN.streamfun_T IN.streamfun_eta IN.streamfun_str IN.streamNumContours IN.streamfun_rho IN.streamfun_v IN.streamfun_vx ...
                                IN.streamfun_vr IN.streamfun_edot IN.streamfun_udw];
    PLOT.streamline         = [PLOT.Streamline IN.streamline_T IN.streamline_eta IN.streamline_str IN.streamline_v IN.streamline_vh IN.streamline_vr];
    PLOT.defmech            = [IN.defmech_T IN.defmech_eta IN.defmech_str IN.defmech_edot IN.defmechContours];
    PLOT.lithoThickness     = [IN.lithoTisovalue IN.lithozmax IN.lithoThickness_T IN.lithoThickness_eta IN.lithoThickness_nstr];
    PLOT.fieldContour       = [IN.fContourIsovalue,IN.fieldContour_C,IN.fieldContour_T,IN.fieldContour_eta,IN.fieldContour_nstr,IN.fieldContour_v,IN.fieldContour_vx,...
                                IN.fieldContour_vr,IN.fieldContour_str,IN.fieldContour_edot,IN.fieldContour_vdiss,IN.fieldContour_psi,IN.fieldContour_udw];
    PLOT.plume              = [PLOT.Plume IN.plume_T IN.plume_eta IN.plume_str IN.plume_rho IN.plume_v IN.plume_edot];
    
    fieldString = '';
    nrFields=0; nrGraphs=0;
    if PLOT.Composition;            nrFields = nrFields+1; fieldString=strcat(fieldString,'-c');  end
    if PLOT.Temperature;            nrFields = nrFields+1; fieldString=strcat(fieldString,'-t');  end
    if PLOT.RegionalResidualT;      nrFields = nrFields+1; fieldString=strcat(fieldString,'-t'); end
    if PLOT.HorizResidualT;         nrFields = nrFields+1; fieldString=strcat(fieldString,'-t'); end
    if PLOT.HorizBandResidualT;     nrFields = nrFields+1; fieldString=strcat(fieldString,'-t'); end
    if PLOT.GlobalResidualT;        nrFields = nrFields+1; fieldString=strcat(fieldString,'-t'); end
    if PLOT.Viscosity;              nrFields = nrFields+1; fieldString=strcat(fieldString,'-eta');  end
    if PLOT.Stress;                 nrFields = nrFields+1; fieldString=strcat(fieldString,'-str');  end
    if PLOT.StressX;                nrFields = nrFields+1; fieldString=strcat(fieldString,'-sx');  end
    if PLOT.StressZ;                nrFields = nrFields+1; fieldString=strcat(fieldString,'-sz');  end
    if PLOT.StrainRate;             nrFields = nrFields+1; fieldString=strcat(fieldString,'-edot');  end
    if PLOT.Density;                nrFields = nrFields+1; fieldString=strcat(fieldString,'-rho');  end
    if PLOT.Velocity;               nrFields = nrFields+1; fieldString=strcat(fieldString,'-v');  end
    if PLOT.VelocityX;              nrFields = nrFields+1; fieldString=strcat(fieldString,'-vx');  end
    if PLOT.VelocityZ;              nrFields = nrFields+1; fieldString=strcat(fieldString,'-vz');  end
    if PLOT.Pressure;               nrFields = nrFields+1; fieldString=strcat(fieldString,'-p');  end
    if PLOT.DynamicPressure;      	nrFields = nrFields+1; fieldString=strcat(fieldString,'-pd');  end
    if PLOT.HeatFlux;               nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-hf'); end
    if PLOT.MeltFraction;         	nrFields = nrFields+1; fieldString=strcat(fieldString,'-f');  end
    if PLOT.CrustalThickness;     	nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-cr');  end
    if PLOT.DeformationMechanism; 	nrFields = nrFields+1; fieldString=strcat(fieldString,'-dlc'); end
    if PLOT.Phase;                  nrFields = nrFields+1; fieldString=strcat(fieldString,'-ph');  end
    if PLOT.Primordial;             nrFields = nrFields+1; fieldString=strcat(fieldString,'-prm');  end
    if PLOT.Air;                    nrFields = nrFields+1; fieldString=strcat(fieldString,'-air');  end
    if PLOT.ContinentalCrust;     	nrFields = nrFields+1; fieldString=strcat(fieldString,'-cc');  end
    if PLOT.Basalt;                 nrFields = nrFields+1; fieldString=strcat(fieldString,'-bs');  end
    if PLOT.Harzburgite;            nrFields = nrFields+1; fieldString=strcat(fieldString,'-hz');  end
    if PLOT.Water;                  nrFields = nrFields+1; fieldString=strcat(fieldString,'-wtr');  end
    if PLOT.Melt;                   nrFields = nrFields+1; fieldString=strcat(fieldString,'-mlt');  end
    if PLOT.SurfaceAge;             nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-sage');  end
    if PLOT.AgeSinceLastMelted;  	nrFields = nrFields+1; fieldString=strcat(fieldString,'-age');  end
    if PLOT.Topography;           	nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-topo');  end
    if PLOT.TopographySelfGrav;    	nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-topoSF');  end
    if PLOT.PlateVelocity;       	nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-pvh');  end
    if PLOT.PlateBaseTopography;   	nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-iltopo');  end
    if PLOT.Toroidal;               nrFields = nrFields+1; fieldString=strcat(fieldString,'-to');  end
    if PLOT.Poloidal;               nrFields = nrFields+1; fieldString=strcat(fieldString,'-po');  end
    if PLOT.Streamfunction;        	nrFields = nrFields+1; fieldString=strcat(fieldString,'-psi');  end
    if PLOT.UpDownWelling;          nrFields = nrFields+1; fieldString=strcat(fieldString,'-udwell');  end
    if PLOT.Geoid;                  nrFields = nrFields+1; fieldString=strcat(fieldString,'-g'); end
    if PLOT.DynamicTopography;   	nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-dtopo'); end
    if PLOT.IsostaticTopography;  	nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-itopo'); end
    if PLOT.ResidualTopography;    	nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-rtopo'); end
    if PLOT.NormalStress;           nrFields = nrFields+1; fieldString=strcat(fieldString,'-nstr'); end
    if PLOT.ViscousDissipation;   	nrFields = nrFields+1; fieldString=strcat(fieldString,'-vdiss'); end
    if PLOT.SurfaceFieldVariation;  nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-sfvar'); end
    if PLOT.PlateSketch;            nrFields = nrFields+1; nrGraphs = nrGraphs+1; fieldString=strcat(fieldString,'-plate'); end
    if PLOT.ParameterTable;      	nrFields = nrFields+1; fieldString=strcat(fieldString,'-partable'); end
    if PLOT.Grid;                   nrFields = nrFields+1; fieldString=strcat(fieldString,'-grid'); end
    if PLOT.Tracers;                nrFields = nrFields+length(PLOT.TracerVariable); fieldString=strcat(fieldString,'-tra'); end
    if PLOT.StreamGraph;            nrFields = nrFields+length(PLOT.StreamGraphParameter); nrGraphs = nrGraphs+length(PLOT.StreamGraphParameter); fieldString=strcat(fieldString,'-sgraph'); end
    if PLOT.CustomGraph;            nrFields = nrFields+length(PLOT.CustomGraphName); nrGraphs = nrGraphs+length(PLOT.CustomGraphName); fieldString=strcat(fieldString,'-tgraph'); end
    
    PLOT.nrSubplot              = 0;
    PLOT.loopPlot               = 0;
    nrModels                    = size(unique(IN.Name),2);
    nrFiles                  	= size(IN.Name,2);
    nrTimeSnapshots             = nrFiles/nrModels;
    nrPlots                     = nrFields*nrFiles;
    
    %% INITIALISE GLOBAL VARIABLES
    SAVE.plateSketchIndication  = false;
    SAVE.legendPlotGraphExists  = 'no';
    SAVE.princStressPlotExists  = false;
    SAVE.quiverPlotExists       = false;
    SAVE.TracerSaved            = false;
    
    %% SWITCH ADJUSTMENTS & ERROR CHECKS
    if nrFields==0; error('No field to plot: choose a field to plot!'); end
    if SWITCH.plotDifference
        if nrFiles~=2  %only 2 simulations/snapshots allowed for direct comparison
            error('To use SWITCH.plotDifference, exactly 2 plots are needed; i.e. size(IN.Name,2) =! 2')
        end
        nrPlots = nrPlots/2;
    end
    if SWITCH.plotGraphVsTime
        SWITCH.PlotInPlot       = true;
    end
    if SWITCH.PlotInPlot
        nrPlots                 = nrFields*nrModels;
        %     if ~isfield(SAVE,'PlotInPlotHax'); SAVE.PlotInPlotHax = zeros(nrPlots)*NaN; end  %initialise vector
        if nrFiles==1 && SAVE.StartNumber==SAVE.EndNumber; error('SWITCH.PlotInPlot switched on makes only sense with more than one file to plot!'); end
        if SWITCH.Texting; SWITCH.Texting = false; if SWITCH.Verbose; disp('  (Texting has been switched off)'); end; end
        if ~PLOT.titleOnlyTime; PLOT.titleOnlyTime = true; end
    end
    if isfield(PLOT,'Layout') && (PLOT.Layout(1,1)*PLOT.Layout(1,2))<nrPlots
        if SWITCH.Verbose; warning(['PLOT.Layout = [',num2str(PLOT.Layout(1,1)),' ',num2str(PLOT.Layout(1,2)),'] yields not enough subplots for the current figure: It is now automatically adjusted!']); end
        PLOT                    = rmfield(PLOT,'Layout');
    end
    if strcmp(STYLE.ColorMode,'dark') && ~SWITCH.PlotDesign
        STYLE.ColorMode = 'light'; warning('Dark colour mode not implemented without SWITCH.PlotDesign being switched on!');
    end
    if strcmp(SWITCH.ColormapAll(1,1),'default'); SWITCH.ColormapAll{1,1} = 'davos'; end
    if SWITCH.PlotDesign && ~SWITCH.SimplifyPlots && ~SWITCH.AnalysisMode
        SWITCH.SimplifyPlots = true; warning('SWITCH.SimplifyPlots has to be switched on if SWITCH.PlotDesign is switched on!');
    end
    if PLOT.indicatePlumes ||  max(PLOT.plume(:))==1 || PLOT.UpDownWelling
        SWITCH.MantleDiagnostics    = true; %will currently be switched off later for YY
    end
    if PLOT.PlateSketch || PLOT.ParameterTable || PLOT.PlateVelocity || ...
            PLOT.indicateTrench || PLOT.indicateRidge || PLOT.indicatePlateFit || ...
            PLOT.indicateShallowSlabDip || PLOT.indicateBending || PLOT.indicateSlabTip || ...
            (PLOT.Topography && TOPO.indicateComponents) || PLOT.DynamicTopography || PLOT.IsostaticTopography || PLOT.ResidualTopography %needed for the UM-density
        SWITCH.PlateDiagnostics     = true;
    end
    if (PLOT.Topography && TOPO.indicateComponents) || PLOT.IsostaticTopography || PLOT.ResidualTopography || PLOT.DynamicTopography
        TOPO.calculateComponents    = true;
    else
        TOPO.calculateComponents    = false;
    end
    if strcmpi(SWITCH.GeodynamicCode,'StagYY') && ~isempty(IN.Number(IN.Number>=99999))
        IN.Number(IN.Number>=99999)	= 99998; %largest possible file number currently
    end
    if strcmpi(SWITCH.GeodynamicCode,'Aspect')
        warning off backtrace;
        warning('Processing HDF5 Output from ASPECT is still in beta, and not fully functional!')
        warning on backtrace
    end
    if PLOT.titleOnlyTime && PLOT.titleNoTime
        warning off backtrace; warning('PLOT.titleOnlyTime and PLOT.titleNoTime cannot be switched on both at the same time!'); warning on backtrace;
        PLOT.titleOnlyTime = false; PLOT.titleNoTime = false;
    end
    % check for length of number input arrays
    if size(IN.Number,2)<nrFiles
        dummy = ones(1,nrFiles-size(IN.Number,2)); dummy(:) = IN.Number(1,1);
        IN.Number           = [IN.Number, dummy];
    end
    if size(IN.Parameter,2)<nrFiles
        dummy = ones(1,nrFiles-size(IN.Parameter,2)); dummy(:) = IN.Parameter(1,1);
        IN.Parameter        = [IN.Parameter, dummy];
    end
    if size(PLOT.magnifierExtent,1)<nrFiles
        dummy = repmat(PLOT.magnifierExtent(1,:),nrFiles-size(PLOT.magnifierExtent,1),1);
        PLOT.magnifierExtent = [PLOT.magnifierExtent; dummy];
    end
    % check for length of string input arrays
    if size(IN.Folder,2)<nrFiles
        dummy = cell(1,nrFiles-size(IN.Folder,2)); dummy(:) = IN.Folder(1,1);
        IN.Folder           = [IN.Folder, dummy];
    end
    if isfield(PLOT,'titleString') && size(PLOT.titleString,2)<nrPlots
        dummy = cell(1,nrPlots-size(PLOT.titleString,2)); dummy(:) = {''};
        PLOT.titleString    = [PLOT.titleString, dummy];
    end
    
    %% CHECK PLATFORM
    if ismac; if SWITCH.Verbose; disp('Platform:          Mac'); end
    elseif isunix; if SWITCH.Verbose; disp('Platform:          Linux'); end
    elseif ispc; if SWITCH.Verbose; disp('Platform:          Windows PC'); end
    end
    
    %% CHECK MATLAB RELEASE VERSION to USE NEWEST FEATURES
    MatlabRelease = version('-release');
    if verLessThan('matlab','8.4.0') %MATLAB 2014a and earlier
        disp(['>>>>>>>       your are currently using MatLab',MatlabRelease,'                  <<<<<<<<<'])
        disp('>>>>>>> to use all available tools update to MatLab2014b or newer! <<<<<<<<<')
        SWITCH.multipleColormaps        = false;
        SWITCH.MultiManualColourMaps    = false;
        if strcmp(SWITCH.ColormapAll(1,1),'auto'); SWITCH.ColormapAll{1,1} = 'parula'; end
        SWITCH.MversionSupported        = false;
    else %MATLAB 2014b and later
        if strcmp(SWITCH.ColormapAll(1,1),'auto') || SWITCH.MultiManualColourMaps
            SWITCH.multipleColormaps    = true;
        else
            SWITCH.multipleColormaps    = false;
        end
        SWITCH.MversionSupported        = true;
    end
    
    %% SETUP FIGURE AND COLOUR SCHEME
    %colour scheme
    PLOT.keyColorLight              = min([1 1 1],max([0 0 0],STYLE.keyColor+0.11));
    if strcmp(STYLE.ColorMode,'light')
        STYLE.ColorModeBW           = 'white';
        STYLE.BlackOrWhiteColor   	= [1 1 1];
    elseif strcmp(STYLE.ColorMode,'dark')
        STYLE.ColorModeBW           = 'black';
        STYLE.BlackOrWhiteColor   	= [0 0 0];
        STYLE.keyColor            	= 1-STYLE.keyColor; %invert key colours
        STYLE.keyLineColor          = 1-STYLE.keyLineColor;
        TOPO.color                  = 1-TOPO.color;
        TOPO.lineColor              = 1-TOPO.lineColor;
        STYLE.annotationColor       = 1-STYLE.annotationColor;
        STYLE.annotationBackColor   = 1-STYLE.annotationBackColor;
    else
        warning(['Colour Mode ''',STYLE.ColorMode,''' not recognised!'])
        STYLE.ColorMode             = 'light';
        STYLE.ColorModeBW           = 'white';
    end
    if ~SWITCH.PlotInPlot
        clf(figure(1));
        set(gcf,'color',STYLE.ColorModeBW); %old version used: colordef(figure(1),STYLE.ColorModeBW);
        %     if PLOT.CustomGraph
        %         warning('colordef can currently erase line colouring in graph plots (reported MatLab bug). Type "reset(groot)" to fix temporarily.')
        %     end
    end
    
    %% MATLAB VERSION ADJUSTMENTS
    if ~verLessThan('matlab', '9.1.0') %MATLAB 2016b and later
        set(gcf,'defaultLegendAutoUpdate','off') %prevent autoupdating legends
    end
    
    %% GET SUBPLOT LAYOUT AUTOMATICALLY
    if ~isfield(PLOT,'Layout') %if not predefined subplot layout
        if SWITCH.PlotInPlot; nrFiles = nrModels; end
        if nrFiles==1 && nrPlots>=4
            subplotLayout   = [2 ceil(nrPlots/2)];
        elseif nrPlots>=4 && nrModels==1
            subplotLayout   = [nrFiles nrPlots/nrFiles];
        else
            subplotLayout   = [nrModels nrPlots/nrModels];
        end
        if SWITCH.ReverseLayout; subplotLayout = fliplr(subplotLayout); end
        PLOT.Layout         = subplotLayout;
    end
    %set panel
    if SWITCH.UsePanel
        SAVE.P	= panel();
        nzSubplots  = PLOT.Layout(1);
        nxSubplots  = PLOT.Layout(2);
        SAVE.P.pack(nzSubplots,nxSubplots);   %set subplot layout
        % set margins
        %             SAVE.P.de.margin = 42; %descendence
        %             SAVE.P(1,1).marginbottom = 12;
        %             SAVE.P(2).marginleft = 20;
        SAVE.P.margin = [32 33 30 26];  %figure edge margins: [left bottom right top]; e.g., [23 24 21 17]
        SAVE.P.marginbottom = 32;
        if SWITCH.Verbose; warning('Panel routine still in testing: Margins might need to be adjusted here'); end
    end
    
    %% ADJUST VARIABLES
    if isfield(SAVE,'GraphPlotHandles'); SAVE.GraphPlotHandles = []; end
    if isfield(SAVE,'PlotHandles'); SAVE.PlotHandles = []; end
    
    %% LOOP FILES
    PLOT.spNum(1,3) = 1; %subplot counter
    for loopCase=1:nrFiles
        if loopCase>1 && strcmp(IN.Name{1,loopCase},IN.Name{1,loopCase-1}) %same file as before
            %skip things?
            fileRepetition  = true;
            loopTimeSnapshot= loopTimeSnapshot+1;
        else
            fileRepetition  = false;
            loopTimeSnapshot= 1;
        end
        
        %% CLEAR OLD VARIABLES
        clearvars GRID
        clearvars io i x2d y2d z2d var2d X3D Y3D Z3D VAR3D topo2d
        fieldsToDelete  = {'topo2d','topo2dp','topoIso2d','topoRes2d','topoDyn2d','isotopo2d'};
        TOPO = rmfield(TOPO,fieldsToDelete(isfield(TOPO,fieldsToDelete)));
        for ib=1:3
            fieldsToDelete  = {'T_3D','ETA_3D','RHO_3D','STR_3D','EDOT_3D','C_3D','BASALT_3D','TOPO_3D','AIR_3D','NSTRESS_3D','DEFM_3',...
                'VD_3D','V_3D','VX_3D','VY_3D','VZ_3D','P_3D'};
            if ib==1; dummy = fieldsToDelete;
            elseif ib==2; dummy = strcat(fieldsToDelete,'yin');
            elseif ib==3; dummy = strcat(fieldsToDelete,'yang');
            end
            PLOT = rmfield(PLOT,dummy(isfield(PLOT,dummy)));
        end
        
        %% SETUP CURRENT VARIABLES
        PLOT.loopCase       = loopCase;
        FILE.name           = IN.Name{1,loopCase}; PLOT.fname = FILE.name;
        FILE.number         = IN.Number(1,loopCase);
        FILE.stemRead     	= IN.Folder{1,loopCase};
        FILE.stemSave       = FILE.stemRead;
        SWITCH.parameter    = IN.Parameter(1,loopCase);
        if SWITCH.PlotInPlot
            if fileRepetition
                PLOT.loopPlot   = max(0,PLOT.loopPlot-1);
                PLOT.nrSubplot  = max(0,PLOT.nrSubplot-1); %reset to previous plot for each case
            else
                if nrTimeSnapshots==1 && isfield(SAVE,'AllFilesTOPOhLine') %only for SAVE.StartNumber:SAVE.EndNumber loops
                    if ~isfield(SAVE,'FirstSuiteTOPOhLine')
                        SAVE.FirstSuiteLegendString	= SAVE.AllFilesLegendString;
                        SAVE.FirstSuiteTOPOhLine    = SAVE.AllFilesTOPOhLine; %remember all lines of first suite
                        %                 else
                        %                     SAVE.FirstSuiteLegendString	= [SAVE.FirstSuiteLegendString, SAVE.AllFilesLegendString];
                        %                     SAVE.FirstSuiteTOPOhLine    = [SAVE.FirstSuiteTOPOhLine, SAVE.AllFilesTOPOhLine]; %remember all lines of first suite
                    end
                end
                if isfield(SAVE,'AllFilesTOPOhLine') && length(SAVE.AllFilesTOPOhLine)>4;   SAVE = rmfield(SAVE,'AllFilesTOPOhLine'); end   %forget about previous handles
                %if isfield(SAVE,'PlotInPlotHax'); SAVE = rmfield(SAVE,'PlotInPlotHax'); end %forget about previous axis
            end
            TOPO.water = false; TOPO.areaFill = 'none';
            %gradual coloring for each time sequence
            if SAVE.EndNumber-SAVE.StartNumber>0 %needs more than one time snapshot to work (using SAVE.StartNumber)
                dDarkening = PLOT.startGradColor./(ceil((SAVE.EndNumber-SAVE.StartNumber)/SAVE.StepNumber)+1-1); %incremental increase in darkness
                PLOT.gradColor4File = PLOT.startGradColor - dDarkening.*ceil((IN.Number(1,1)-SAVE.StartNumber)/SAVE.StepNumber-0.999); %last one is black [0 0 0]
                clearvars dDarkening
            elseif nrTimeSnapshots>1 %or needs more than one timesnapshot to work (using FILE.number)
                dDarkening = PLOT.startGradColor./(nrTimeSnapshots-1); %incremental increase in darkness
                PLOT.gradColor4File = PLOT.startGradColor - dDarkening.*(loopTimeSnapshot-1); %last one is black [0 0 0]
                clearvars dDarkening
            elseif nrFiles>1 %or needs more than one file to work
                dDarkening = PLOT.startGradColor./(nrFiles-1); %incremental increase in darkness
                PLOT.gradColor4File = PLOT.startGradColor - dDarkening.*(loopCase-1); %last one is black [0 0 0]
                clearvars dDarkening
            else
                error('Plot-In-Plot only works with multiple models or time snapshots!');
            end
            if strcmp(STYLE.ColorMode,'dark')
                PLOT.gradColor4File	= 1-PLOT.gradColor4File;
            end
        end
        
        %% CHECK FOR INPUT FILE
        [FILE] = f_FileFinder(FILE,SWITCH,SAVE,STYLE);
        if FILE.NotFound
            SAVE.LastFile     	= true;
            close all
            disp([FILE.NotFoundError])
            return
            
        end
        if FILE.number==SAVE.EndNumber %might also be last file as defined by user
            FILE.LastFile       = true;
        end
        SAVE.LastFile           = FILE.LastFile;
        if SAVE.StartNumber==SAVE.EndNumber; SAVE.LastFile = true; end %only one figure
        PLOT.number             = FILE.number;
        IN.Number(1,loopCase)   = FILE.number;
        
        %% READ INPUT FILE  (for input field and grid variables)
        [READ] = f_readStagYY(FILE.directory,FILE.name,FILE.number,FILE.foundFieldName,SWITCH,2);
        if ischar(READ); error(READ); end
        % GRID ADJUSTMENTS
        if  isfield(READ,'X_3D') && size(READ.X_3D,1)==1 %exchange x and y
            dummy_3D = zeros(size(READ.Y_3D,2),size(READ.Y_3D,1),size(READ.Y_3D,3));
            dummy_3Dx = zeros(size(READ.Y_3D,2),size(READ.Y_3D,1),size(READ.Y_3D,3));
            dummy_3Dx(:,1,:) = READ.X_3D(1,:,:);
            dummy_3D(:,1,:) = READ.Y_3D(1,:,:);     READ.Y_3D = dummy_3Dx; READ.X_3D = dummy_3D;
            dummy_3D(:,1,:) = READ.Z_3D(1,:,:);     READ.Z_3D = dummy_3D;
            dummy_3D(:,1,:) = READ.VAR_3D(1,:,:); 	READ.VAR_3D = dummy_3D;
            READ.Aspect = flipud(READ.Aspect);
            clearvars dummy_3D dummy_3Dx
        end
        X_3D = READ.X_3D; Y_3D = READ.Y_3D; Z_3D = READ.Z_3D;
        
        %% EXTRACT DETAILS OF INPUT FILE
        GRID.aspectDomain = READ.Aspect; %aspect ratio of model-domain (not used yet!)
        %grid size
        nx = size(Z_3D,1); ny = size(Z_3D,2); nz = size(Z_3D,3);
        if ~strcmp(FILE.foundFieldNature,'field')
            nz          = NaN; %nz can not be determined from e.g. graph data like topography fields
        end
        %grid type
        if strcmpi(SWITCH.GeodynamicCode,'StagYY')
            if READ.nb>1 %check grid type
                GRID.Type = 'yinyang';
            elseif READ.rcmb==-1 %this is stag's flag for Cartesian geometry
                GRID.Type = 'Cartesian';
            else
                GRID.Type = 'spherical2D';
            end
        elseif strcmpi(SWITCH.GeodynamicCode,'Fluidity')
            GRID.Type = 'Cartesian';  %only Cartesian geometry implemented
        elseif strcmpi(SWITCH.GeodynamicCode,'Aspect')
            GRID.Type = 'Cartesian';  %only Cartesian geometry implemented
        else
            error('adjust here for a different geodynamic code!')
        end
        if ((nx>1) && (ny>1) && (nz>1)) || strcmp(GRID.Type,'yinyang') %check dimension
            GRID.Dim    = '3-D';
        else
            GRID.Dim    = '2-D';
        end
        GRID.nx = nx; GRID.ny = ny; GRID.nz = nz; GRID.nb = READ.nb;
        %data format
        if strcmp(FILE.foundFieldNature,'field') && max(Z_3D(:))>0.5 && max(Z_3D(:))<1.5  %if 0.5<D<1.5 then non-dimensional  %%%NOT PERFECT
            SWITCH.DimensionalInput = false; dimString = 'non-dimensional';
        elseif ~strcmp(FILE.foundFieldNature,'field')
            if strcmp(GRID.Type,'Cartesian') && max(X_3D(:))>0.5 && max(X_3D(:))<9.0 %Z_3D not available
                SWITCH.DimensionalInput = false; dimString = 'non-dimensional';
            else
                warning('Data format (dimensional or non-dimensional) could not be determined!')
                SWITCH.DimensionalInput = false; dimString = 'unknown'; %unknown
                SWITCH.DimensionalMode 	= false;
            end
        else
            SWITCH.DimensionalInput = true; dimString = 'dimensional';
        end
        
        %% DISPLAY INFO
        disp(['Data:              ',SWITCH.GeodynamicCode,', ',dimString])
        gridSizeString = '';
        if nx>1; gridSizeString = num2str(nx); end
        if ny>1; if nx>1; gridSizeString = [gridSizeString,STYLE.SCHAR.timesCross,num2str(ny)]; else; gridSizeString = num2str(ny); end; end
        if nz>1; gridSizeString = [gridSizeString,STYLE.SCHAR.timesCross,num2str(nz)]; end
        if isnan(nz); gridSizeString = [gridSizeString,STYLE.SCHAR.timesCross,num2str(nz)]; end %in case nz could not be determined
        if isfield(READ,'nb') && READ.nb>1; gridSizeString = [gridSizeString,STYLE.SCHAR.timesCross,num2str(READ.nb)]; end
        disp(['Model Domain:      ',GRID.Dim,' ',GRID.Type,' (',gridSizeString,')'])
        
        %% fAIO ACTION: UPDATE CACHE
        if SAVE.count==1 %only first time during one execution
            fAIo.task   = 'handleCache';
            [fAIo,SAVE] = f_AIo(fAIo,SWITCH,FILE,GRID,STYLE,SAVE);
        end
        
        %% ADJUSTMENTS FOR GRID GEOMETRIES
        if SWITCH.Magnifier && ~strcmp(GRID.Type,'Cartesian')
            %         SWITCH.Magnifier     	= false;	warning(['Magnifier not implemented for ',GRID.Type,' geometry!']);
        end
        if strcmp(GRID.Type,'yinyang') && ~SWITCH.UsePanel
            if strcmp(SWITCH.yyPlotMode,'map')
                if ~isfield(STYLE,'cbShiftX');	STYLE.cbShiftX = 0.0; end
                STYLE.cbShiftX          = STYLE.cbShiftX -0.04; %*(PLOT.Layout(2)-1); %shift colorbar to the left
            else %mapOnSphere and isosurface
                if ~isfield(STYLE,'cbShiftX');	STYLE.cbShiftX = 0.0; end
                STYLE.cbShiftX          = STYLE.cbShiftX -0.04; %shift colorbar to the left
            end
        end
        if SWITCH.onlyXaxis && ~strcmp(GRID.Type,'spherical2D')
            SWITCH.onlyXaxis        = false;
        end
        if strcmp(GRID.Type,'spherical2D') && GRID.aspectDomain(1,1)<6.2
            GRID.partialAnnulus    = true;
            SWITCH.closeAnnulus    = false;    %switch off for partial spherical2D geometries
        else
            GRID.partialAnnulus    = false;
        end
        if SWITCH.AxesLimit || (strcmp(GRID.Type,'spherical2D') && GRID.partialAnnulus)
            SWITCH.spherical2DCenterText = false;
        end
        
        %% ADJUST CYLINDRICAL POSITION
        GRID.partialAnnulus    = false;
        if strcmp(GRID.Type,'spherical2D')
            minX                = min(READ.X_3D(:));
            maxX                = max(READ.X_3D(:));
            diffX               = maxX-minX;
            if diffX<6.2    %only for partial spherical annuluss
                GRID.partialAnnulus = true;
                X_3D            = X_3D +(pi-diffX/2);
                READ.X_3D       = READ.X_3D +(pi-diffX/2);
            end
        end
        
        %% KEEP T-FIELD FOR LATER USE
        if strcmp(FILE.foundFieldName,'Temperature')
            if strcmp(GRID.Type,'yinyang') %save T-field for later use
                PLOT.T_3Dyin    = READ.VAR_3Dyin;
                PLOT.T_3Dyang   = READ.VAR_3Dyang;
            elseif strcmp(GRID.Dim,'2-D')
                PLOT.T_3D       = READ.VAR_3D;
            elseif strcmp(GRID.Dim,'3-D')
                PLOT.T_3D       = READ.VAR_3D;
            end
        end
        
        %% READ VELOCITY/PRESSURE FILE, if vp-data exists, AND KEEP FOR LATER USE
        try
            [READ2] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Velocity',SWITCH,2);
            if  isfield(READ2,'VX_3D') && size(READ2.VX_3D,1)==1 %exchange x and y
                dummy_3D = zeros(size(READ2.VX_3D,2),size(READ2.VX_3D,1),size(READ2.VX_3D,3));
                dummy_3Dx = dummy_3D;
                dummy_3Dx(:,1,:)    = READ2.VX_3D(1,:,:);
                dummy_3D(:,1,:)     = READ2.VY_3D(1,:,:);  	READ2.VY_3D = dummy_3Dx;    READ2.VX_3D = dummy_3D;
                dummy_3D(:,1,:)     = READ2.VZ_3D(1,:,:); 	READ2.VZ_3D = dummy_3D;
                dummy_3D(:,1,:)     = READ2.P_3D(1,:,:);   	READ2.P_3D = dummy_3D;
                clearvars dummy_3D dummy_3Dx
            end
            if strcmp(GRID.Type,'yinyang')
                PLOT.VX_3Dyin = READ2.VX_3Dyin;	PLOT.VX_3Dyang = READ2.VX_3Dyang; %save velocity fields for later use
                PLOT.VY_3Dyin = READ2.VY_3Dyin;	PLOT.VY_3Dyang = READ2.VY_3Dyang;
                PLOT.VZ_3Dyin = READ2.VZ_3Dyin;	PLOT.VZ_3Dyang = READ2.VZ_3Dyang;
                PLOT.P_3Dyin = READ2.P_3Dyin; 	PLOT.P_3Dyang = READ2.P_3Dyang;
            else
                if strcmp(GRID.Type,'spherical2D')
                    PLOT.VX_3Ds     = -READ2.VZ_3D.*sin(READ.X_3D)-READ2.VX_3D.*cos(READ.X_3D);
                    PLOT.VZ_3Ds     = READ2.VX_3D.*sin(READ.X_3D)-READ2.VZ_3D.*cos(READ.X_3D);
                end
                PLOT.VX_3D = READ2.VX_3D; PLOT.VY_3D = READ2.VY_3D; PLOT.VZ_3D = READ2.VZ_3D; %save velocity fields for later use
                PLOT.P_3D   = READ2.P_3D; %save pressure field for later
            end
        catch
            if SWITCH.Verbose; warning('Velocity file not found!'); end
        end
        clearvars READ2 velocityFile dummy_3Dx dummy_3D
        
        %% READ TIME.DAT FILE, if exists
        filestemTimedat   	= strcat(FILE.directory,FILE.name,'_time.dat');
        timedatExists    	= true;
        if ~exist(filestemTimedat,'file')
            timedatExists  	= false;
            if SWITCH.Verbose; warning([filestemTimedat,' could not be found.']); end
        else
            try
                %find number of header entries in first line automatically
                data1                   = fopen(filestemTimedat,'r');
                iEntry                  = 0;
                entry                   = 1;
                while (entry~='0')
                    entry               = fscanf(data1,'%s',[1 1]);
                    iEntry              = iEntry+1;
                    TIMEDAT.title{1,iEntry} = entry;
                end
                fclose(data1);
                TIMEDAT.title(:,end)    = [];
                %find number of columns
                fid                     = fopen(filestemTimedat,'r');
                nEntries=0; iLine=0;
                while nEntries<10 %line has less than 10 entries/columns (i.e., is empty)
                    iLine               = iLine+1;
                    nEntries            = numel(regexp(fgetl(fid),'\s*([^\s]*)\s*'));
                end
                fclose(fid);
                try
                    TIMEDAT.filename        = filestemTimedat;
                    TIMEDAT.titleRow        = true;
                    TIMEDAT.numberEntries  	= nEntries;
                    TIMEDAT.startRow       	= iLine;
                    if nEntries==0
                        TIMEDAT.titleRow 	= false;
                        TIMEDAT.startRow 	= 1;
                    end
                    TIMEDAT.Task            = 'ImportTimedat';
                    [TIMEDAT,~] = f_Import(TIMEDAT,FILE,[],[],SWITCH);
                    PLOT.timedatStep        = TIMEDAT.Array{1,strcmp(TIMEDAT.title,'istep')};
                    PLOT.timedatTime        = TIMEDAT.Array{1,strcmp(TIMEDAT.title,'time')};
                    PLOT.timedatVrms        = TIMEDAT.Array{1,strcmp(TIMEDAT.title,'Vrms')};
                catch erme
                    timedatExists           = false;
                    if SWITCH.Verbose; warning(erme); end
                end
                if ~(nEntries==39 || nEntries==30 || nEntries==27 || nEntries==17 || nEntries==0)
                    warning('Number of columns in time.dat is different from expected formats! Check if this needs an update!')
                end
            catch
                timedatExists           = false;
                if SWITCH.Verbose; warning(['There was a problem reading ',filestemTimedat,'.']); end
            end
            clearvars iEntry entry nEntries iLine data1
        end
        
        %% READ TRACER INFORMATION
        if ~strcmp(GRID.Type,'yinyang')
            %tracerFileExists        = true;
            [NumTracers,NumVariables,TraVarName,~,TraIdealMass,TraNumTraceElements,TraOutgassedAmount,TraTime,TraRcmb,TraTimeStep,TraAspect,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Tracer Info',SWITCH);
            if ischar(NumTracers(1,1))
                %tracerFileExists  	= false;
                if SWITCH.Verbose; disp(['-> ',NumTracers]); end
            else
                for iVar=1:NumVariables
                    if iVar==1; TraVariableNamesString = TraVarName{iVar};
                    else; TraVariableNamesString = [TraVariableNamesString,',',TraVarName{iVar}]; end
                end
                if SWITCH.Verbose
                    disp('Tracers            ')
                    disp([' Number (Type):    ',num2str(NumTracers),' x ',num2str(NumVariables),' (',TraVariableNamesString,')'])
                    disp([' Time Step (Time): ',num2str(TraTimeStep),' (',num2str(TraTime,3),')'])
                    disp([' Ideal Mass:       ',num2str(TraIdealMass)])
                    disp([' Box Aspect Ratio: ',num2str(TraAspect(1),3),' x ',num2str(TraAspect(2),3)])
                    disp([' Radius CMB:       ',num2str(TraRcmb)])
                    disp([' #Trace Elements:  ',num2str(TraNumTraceElements)])
                    disp([' Outgassed Amount: ',num2str(TraOutgassedAmount)])
                else
                    disp(['Tracers:           ',num2str(NumTracers),' x ',num2str(NumVariables),' (',TraVariableNamesString,')'])
                end
            end
        else
            if SWITCH.Verbose; warning('Reading Tracers currently omitted for yinyang models!'); end
        end
        
        %% SET FIGURE SIZE
        if loopCase==1
            h=figure(1);
            if strcmp(SWITCH.FigurePosition,'auto')     %set figure position automatically
                FPOS.subplotLayout      = PLOT.Layout;
                FPOS.Default            = PLOT.FigureDefaultPosition;
                FPOS.setAxisLimit       = SWITCH.AxesLimit;
                FPOS.axisLimit          = SWITCH.AxesLimitValues;
                FPOS.nrGraphs           = nrGraphs;
                FPOS.GraphAspectX       = PLOT.MaxAspectX;
                if strcmp(GRID.Type,'yinyang'); FPOS.yyNumSlices = PLOT.yyNumSlices; end
                [FPOS] = f_SetupFigure(FPOS,SAVE,GRID,SWITCH);
                if nrGraphs>0 && FPOS.spAspectRatio>PLOT.MaxAspectX
                    PLOT.shiftUpward   	= PLOT.shiftUpward-0.03;
                end
            else
                set(h,'Position',SWITCH.FigurePosition);     %set figure position manually
            end
            clearvars FPOS
        end
        
        %% DIMENSIONAL PARAMETERS
        [SETUP,SWITCH] = f_Dimensions(SWITCH,FILE);
        if strcmp(GRID.Dim,'2-D')
            SETUP.volumescalePlot   = SETUP.volumescale;    %[m^2]
            SETUP.volumeDimPlot     = SETUP.areaDim;        %[m^2]
        else %3-D
            SETUP.volumescalePlot   = SETUP.volumescale;    %[m^3]
            SETUP.volumeDimPlot     = SETUP.volumeDim;      %[m^3]
        end
        if strcmp(SETUP.topBC,'sticky-air')
            SWITCH.ActualDepth  = true;
        else
            SWITCH.ActualDepth  = false;
        end
        disp(['Top BC:            ',SETUP.topBC])
        
        %% CHECK FOR SIDE-BOUNDARY CONDITION if possible
        GRID.BCxx = 'unknown'; %default
        GRID.BCyy = 'unknown'; %default
        if isfield(PLOT,'VX_3D') && isfield(PLOT,'VY_3D')
            if size(PLOT.VX_3D,1)>1 %x-boundary
                if sum(PLOT.VX_3D(1,ceil(end/2),:))==0
                    GRID.BCxx   = 'impermeable';
                else
                    GRID.BCxx   = 'permeable';
                end
            end
            if size(PLOT.VY_3D,2)>1 %y-boundary
                if sum(PLOT.VY_3D(ceil(end/2),1,:))==0
                    GRID.BCyy   = 'impermeable';
                else
                    GRID.BCyy   = 'permeable';
                end
            end
            if strcmp(GRID.Type,'yinyang')
            elseif ~strcmp(GRID.BCxx,'unknown') && ~strcmp(GRID.BCyy,'unknown') %3-D
                disp(['Side BCs:          ',GRID.BCxx,'(x), ',GRID.BCyy,'(y)'])
            elseif ~strcmp(GRID.BCxx,'unknown') %2-D x
                disp(['Side BC:           ',GRID.BCxx])
            elseif ~strcmp(GRID.BCyy,'unknown') %2-D y
                disp(['Side BC:           ',GRID.BCyy])
            end
        end
        
        %% SETUP TIME
        PLOT.time       = READ.time;     %[nd] or [dim]
        PLOT.time_dim   = READ.time*SETUP.timescale/SETUP.secyear;  % get dimensional time [yr]
        if PLOT.time_dim/1e6>100
            PLOT.time2plot = PLOT.time_dim/1e9; PLOT.time2plotDim = 'Gyr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear/1e9;
        elseif PLOT.time_dim/1e6<1e-7 %0.0000001
            PLOT.time2plot = PLOT.time_dim*SETUP.secyear; PLOT.time2plotDim = 's'; PLOT.timeConvert = SETUP.timescale;
        elseif PLOT.time_dim/1e6<1e-4 %0.0001
            PLOT.time2plot = PLOT.time_dim; PLOT.time2plotDim = 'yr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear;
        elseif PLOT.time_dim/1e6<1e-1 %0.1
            PLOT.time2plot = PLOT.time_dim/1e3; PLOT.time2plotDim = 'kyr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear/1e3;
        else
            PLOT.time2plot = PLOT.time_dim/1e6; PLOT.time2plotDim = 'Myr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear/1e6;
        end
        if ~SWITCH.DimensionalMode && ~SWITCH.DimensionalInput; PLOT.time2plot = PLOT.time; PLOT.time2plotDim = 'nd'; PLOT.timeConvert = 1; end
        
        %SET MANTLE TURN-OVER TIME: t_mt = t_dim*v_rms/D
        % with v_rms from the beginning to the current timestep
        if timedatExists %needs data from time.dat
            currentStep    	= find(PLOT.timedatStep==FILE.number);
            v_rmsCurrent    = mean(PLOT.timedatVrms(1:currentStep));
            PLOT.timeMT     = PLOT.time*v_rmsCurrent/SETUP.D;
            PLOT.dtMT       = SETUP.D/v_rmsCurrent *PLOT.timeConvert;
        else
            PLOT.timeMT     = NaN;
            PLOT.dtMT       = NaN;
        end
        clearvars currentStep
        
        % INDICATE TIME
        timeString0	= [num2str(PLOT.time,3),' | ',num2str(PLOT.time2plot,3),' ',PLOT.time2plotDim];
        if ~SWITCH.DimensionalMode && ~SWITCH.DimensionalInput; timeString0	= num2str(PLOT.time,3); end
        timeString	= timeString0;
        if ~isnan(PLOT.timeMT)
            timeString = [timeString0,' | ',num2str(PLOT.timeMT,2),' mantle turnovers (vrms = ',num2str(v_rmsCurrent,2),', dt_MT = ',num2str(PLOT.dtMT,2),' ',PLOT.time2plotDim,')'];
        end
        disp(['Time:              ',timeString])
        
        %% SETUP GRID
        if strcmp(GRID.Dim,'3-D')
            x2d         = X_3D(:,1,1);
            y2d         = Y_3D(1,:,1);
            z2d         = Z_3D(1,1,:);
        elseif strcmp(GRID.Dim,'2-D')
            x2d(:,:)    = X_3D(:,1,:);
            y2d(:,:)    = Y_3D(:,1,:);	%just in order to have it defined
            z2d(:,:)    = Z_3D(:,1,:);
        end
        
        % FLIP DEPTH VECTOR AND ACCOUNT FOR STICKY-AIR LAYER
        GRID.Z_3Dndnf   = Z_3D; %not flipped z array (non-dimensionalised later)
        if SWITCH.DimensionalInput
            z2d         = SETUP.D_dim-z2d;      %reverse z-axis  (1.)
            Z_3D        = SETUP.D_dim-Z_3D;
        else
            z2d         = 1.0-z2d;              %reverse z-axis  (1.)
            Z_3D        = 1.0-Z_3D;
        end
        if SWITCH.ActualDepth
            z2d         = z2d-SETUP.d_air;      %set z=0 to real surface
            Z_3D        = Z_3D-SETUP.d_air;     %set z=0 to real surface
        end
        
        % DIMENSIONALISATION OF GLOBAL GRID VARIABLES
        % GRID.dimFactor    - conversion factor from input dimension to plot dimension
        % GRID.ndFactor     - conversion factor from input dimension to non-dimensional
        % GRID.m2p          - conversion factor from [m] to plot dimension
        if SWITCH.DimensionalInput
            GRID.ndFactor   = 1/SETUP.D_dim;            %[nd]
            if isfield(GRID,'Z_3Dndnf'); GRID.Z_3Dndnf = GRID.Z_3Dndnf/SETUP.D_dim; end %[nd]
        else
            GRID.ndFactor   = 1;                        %[nd]
        end
        if SWITCH.DimensionalMode || SWITCH.DimensionalInput
            GRID.dimFactor  = SETUP.lengthscale;        %[m]
            %plotting variables
            if abs(max(max(Z_3D(:)*SETUP.lengthscale)))>1e3 % >1km
                GRID.dimFactor  = GRID.dimFactor /1e3; %[km]
                GRID.Xdim = 'km'; GRID.Ydim = 'km'; GRID.Zdim = 'km';
                GRID.m2p        = 1e-3;                 %[m]->[km]
            elseif abs(max(max(Z_3D(:)*SETUP.lengthscale)))>1 % >1m
                GRID.dimFactor  = GRID.dimFactor;       %[m]
                GRID.Xdim = 'm'; GRID.Ydim = 'm'; GRID.Zdim = 'm';
                GRID.m2p        = 1;                    %[m]->[m]
            else % <1m
                GRID.dimFactor  = GRID.dimFactor *1e2;  %[cm]
                GRID.Xdim = 'cm'; GRID.Ydim = 'cm'; GRID.Zdim = 'cm';
                GRID.m2p        = 1e2;                  %[m]->[cm]
            end
        else
            GRID.dimFactor  = GRID.ndFactor;            %[nd]
            GRID.Xdim = 'nd'; GRID.Ydim = 'nd'; GRID.Zdim = 'nd';
            GRID.m2p        = 1/SETUP.D_dim;            %[m]->[nd]
        end
        if strcmp(GRID.Type,'spherical2D') || strcmp(GRID.Type,'yinyang') %horizontal grid values are in degrees
            ndFactor    = 1;
            dimFactor 	= 1;
            lenghtscale = 1;
        else
            ndFactor    = GRID.ndFactor;
            dimFactor   = GRID.dimFactor;
            lenghtscale = SETUP.lengthscale;
        end
        GRID.X_3Dnd     = X_3D .*ndFactor;              %[nd]
        GRID.Y_3Dnd     = Y_3D .*ndFactor;              %[nd]
        GRID.Z_3Dnd     = Z_3D .*GRID.ndFactor;         %[nd]
        GRID.x2d_nd     = x2d .*ndFactor;               %[nd]
        GRID.y2d_nd     = y2d .*ndFactor;               %[nd]
        GRID.z2d_nd     = z2d .*GRID.ndFactor;          %[nd]
        GRID.rcmb_nd    = READ.rcmb *GRID.ndFactor;     %[nd]
        GRID.X_3Dp      = X_3D .*dimFactor;             %[nd], [km], [m], or [cm]
        GRID.Y_3Dp      = Y_3D .*dimFactor;
        GRID.Z_3Dp      = Z_3D .*GRID.dimFactor;
        GRID.x2dp       = x2d .*dimFactor;
        GRID.y2dp       = y2d .*dimFactor;
        GRID.z2dp       = z2d .*GRID.dimFactor;
        GRID.rcmb_p     = READ.rcmb *GRID.dimFactor;
        GRID.nd2dim     = 1/GRID.ndFactor*GRID.dimFactor;
        if SWITCH.DimensionalMode || SWITCH.DimensionalInput
            X_3D        = X_3D .*lenghtscale;           %[m] SI-units
            Y_3D        = Y_3D .*lenghtscale;           %[m]
            Z_3D        = Z_3D .*SETUP.lengthscale;     %[m]
            %         GRID.x2d	= x2d .*lenghtscale;	%[m]
            %         GRID.y2d	= y2d .*lenghtscale;	%[m]
            %         GRID.z2d  	= z2d .*SETUP.lengthscale;	%[m]
            READ.rcmb        = READ.rcmb *SETUP.lengthscale;	%[m]
        end
        GRID.Rnd       	= GRID.Z_3Dndnf+GRID.rcmb_nd;   %[nd] (1:cmb, end:surface)
        GRID.Rp     	= GRID.Rnd .*GRID.nd2dim;       %[plotting dimension]
        if SWITCH.DimensionalMode || SWITCH.DimensionalInput; GRID.R = GRID.Rp./GRID.m2p; else; GRID.R = GRID.Rp; end  %[m] or [nd]
        
        %SPHERICAL COORDINATES
        if strcmp(GRID.Type,'spherical2D') || strcmp(GRID.Type,'yinyang') %horizontal grid values are in degrees
            %conversion from radian to horizontal distance (L = r*alpha)
            GRID.X_3Dsp     = GRID.X_3Dnd .*GRID.Rp;    %[distance]
            GRID.Y_3Dsp     = GRID.Y_3Dnd .*GRID.Rp;    %[distance]
            if strcmp(GRID.Type,'spherical2D')
                clearvars dummy
                dummy(:,:)  = GRID.Rp(:,1,:);
                GRID.x2dsp  = GRID.x2d_nd .*dummy;      %[distance]
            end
        end
        if strcmp(GRID.Type,'spherical2D')
            if SWITCH.closeAnnulus %add a row to close spherical annulus plot
                GRID.x2d_nd(end+1,:) = GRID.x2d_nd(1,:); GRID.Z_3Dndnf(end+1,:,:) = GRID.Z_3Dndnf(1,:,:);
                GRID.Rnd(end+1,:,:) = GRID.Rnd(1,:,:); GRID.Rp(end+1,:,:) = GRID.Rp(1,:,:);
            end
            %Transformation to spherical2D coordinates (x-data needs to be in degrees)
            clearvars Rnd2d
            Rnd2d(:,:)  = GRID.Rp(:,1,:);
            GRID.x2ds   = Rnd2d.*-sin(GRID.x2d_nd); %.*cos(y2d_nd);
            GRID.z2ds   = Rnd2d.*-cos(GRID.x2d_nd);
            if SWITCH.closeAnnulus %remove row again
                GRID.x2d_nd(end,:) = []; GRID.Z_3Dndnf(end,:,:) = [];
                GRID.Rnd(end,:,:) = []; GRID.Rp(end,:,:) = [];
            end
        end
        
        %YINYANG COORDINATES
        %...
        
        %DEFINE GLOBAL GRID VARIABLES
        GRID.X_3D = X_3D;   GRID.Y_3D = Y_3D;   GRID.Z_3D = Z_3D; %[nd] or [m] SI-units
        GRID.nx   = nx;   	GRID.ny   = ny;   	GRID.nz   = nz;
        %     GRID.aspectRatio(1) = nx/nz; GRID.aspectRatio(2) = ny/nz;
        
        %SET GRID SPACING
        %dublicate last row for later consistency of array sizes
        if size(GRID.X_3D,1)>1;     GRID.dx = abs(GRID.X_3D(2:end,:,:)-GRID.X_3D(1:end-1,:,:));     GRID.dx(end+1,:,:) = GRID.dx(end,:,:);      else; GRID.dx = 1; end %[m] or [nd]
        if size(GRID.Y_3D,2)>1;     GRID.dy = abs(GRID.Y_3D(:,2:end,:)-GRID.Y_3D(:,1:end-1,:));     GRID.dy(:,end+1,:) = GRID.dy(:,end,:);      else; GRID.dy = 1; end %[m] or [nd]
        if size(GRID.Z_3D,3)>1;     GRID.dz = abs(GRID.Z_3D(:,:,2:end)-GRID.Z_3D(:,:,1:end-1));     GRID.dz(:,:,end+1) = GRID.dz(:,:,end);      else; GRID.dz = 1; end %[m] or [nd]
        if size(GRID.X_3Dp,1)>1;    GRID.dxp = abs(GRID.X_3Dp(2:end,:,:)-GRID.X_3Dp(1:end-1,:,:));  GRID.dxp(end+1,:,:) = GRID.dxp(end,:,:);    else; GRID.dxp = 1; end %[plotting dim]
        if size(GRID.Y_3Dp,2)>1;    GRID.dyp = abs(GRID.Y_3Dp(:,2:end,:)-GRID.Y_3Dp(:,1:end-1,:));  GRID.dyp(:,end+1,:) = GRID.dyp(:,end,:);    else; GRID.dyp = 1; end %[plotting dim]
        if size(GRID.Z_3Dp,3)>1;    GRID.dzp = abs(GRID.Z_3Dp(:,:,2:end)-GRID.Z_3Dp(:,:,1:end-1));  GRID.dzp(:,:,end+1) = GRID.dzp(:,:,end);    else; GRID.dzp = 1; end %[plotting dim]
        if strcmp(GRID.Type,'spherical2D') || strcmp(GRID.Type,'yinyang') %convert degrees to distance (d = theta*r)
            GRID.dxr    	= GRID.dx;  	%[rad]
            GRID.dyr     	= GRID.dy;      %[rad]
            if size(GRID.X_3D,1)>1;     GRID.dx = GRID.dx.*GRID.R;      end 	%[m] or [nd]
            if size(GRID.Y_3D,2)>1;     GRID.dy = GRID.dy.*GRID.R;      end   	%[m] or [nd]
            if size(GRID.X_3Dp,1)>1;    GRID.dxp = GRID.dxp.*GRID.Rp;   end     %[plotting dim]
            if size(GRID.Y_3Dp,2)>1;    GRID.dyp = GRID.dyp.*GRID.Rp;   end     %[plotting dim]
        end
        
        %GRID CELL VOLUME
        %last row is dublicated for subsequent consistency with array sizes of fields
        GRID.cellVolume = zeros(size(GRID.X_3D)); GRID.cellVolumeP = GRID.cellVolume;
        GRID.cellVolume         = GRID.dx .* GRID.dy .* GRID.dz;        %vector [m] or [nd]
        GRID.cellVolumeP        = GRID.dxp .* GRID.dyp .* GRID.dzp;     %vector [plotting dim]
        if strcmp(GRID.Dim,'2-D'); GRID.VolDim = [GRID.Zdim,'^2']; else; GRID.VolDim = [GRID.Zdim,'^3']; end
        if strcmp(GRID.Zdim,'nd'); GRID.VolDim = 'nd'; end
        %grid cell output
        GRID.spacing = 'uniform'; %dummy
        if range(GRID.dxp(:))<1e-10 %uniform x-spacing
            dxString        = ['dx = ',num2str(GRID.dxp(1),3),' ',GRID.Xdim];
        else %variable x-spacing
            GRID.spacing    = 'variable';
            dxString        = ['dx = ',num2str(min(GRID.dxp(:)),3),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(max(GRID.dxp(:)),3),' ',GRID.Xdim];
        end
        if range(GRID.dyp(:))<1e-10 %uniform y-spacing
            dyString        = [', dy = ',num2str(GRID.dyp(1),3),' ',GRID.Ydim];
        else %variable y-spacing
            GRID.spacing    = 'variable';
            dyString        = [', dy = ',num2str(min(GRID.dyp(:)),3),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(max(GRID.dyp(:)),3),' ',GRID.Ydim];
        end
        if range(GRID.dzp(:))<1e-10 %uniform z-spacing
            dzString        = [', dz = ',num2str(GRID.dzp(1),3),' ',GRID.Zdim];
        else %variable z-spacing
            GRID.spacing    = 'variable';
            dzString        = [', dz = ',num2str(min(GRID.dzp(:)),3),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(max(GRID.dzp(:)),3),' ',GRID.Zdim];
        end
        if strcmp(GRID.spacing,'uniform')
            gridCellVolString = [num2str(mean(GRID.cellVolumeP(:)),3),' ',GRID.VolDim];
            if round(mean(GRID.cellVolumeP(:)))~=round(GRID.cellVolumeP(1)); warning('cell volume is NOT constant!'); end
        else %variable
            gridCellVolString = [num2str(min(GRID.cellVolumeP(:)),3),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(max(GRID.cellVolumeP(:)),3),' ',GRID.VolDim];
        end
        if strcmp(GRID.Dim,'2-D')
            dyString        = '';
        end
        gridSpacingString   = [dxString,dyString,dzString];
        
        % DISPLAY INFO
        disp(['Grid Spacing:      ',GRID.spacing]);
        disp(['                   ',gridSpacingString])
        if SWITCH.Verbose; disp(['Grid Cell Volume:  ',gridCellVolString]); end
        
        %% DEFINE ASPECT RATIO OF SUBPLOTS
        if strcmp(GRID.Type,'spherical2D') || strcmp(GRID.Type,'yinyang')
            %         zidx = round(size(GRID.X_3Dsp,3)/2); %+/- mid-mantle
            zidx            = 1; %cmb
            GRID.aspectRatio(1) = abs(max(GRID.X_3Dsp(:,zidx))-min(GRID.X_3Dsp(:,zidx))) / abs(max(GRID.Z_3Dp(:))-min(GRID.Z_3Dp(:)));
            GRID.aspectRatio(2) = abs(max(GRID.Y_3Dsp(:,zidx))-min(GRID.Y_3Dsp(:,zidx))) / abs(max(GRID.Z_3Dp(:))-min(GRID.Z_3Dp(:)));
        else
            GRID.aspectRatio(1) = abs(max(GRID.X_3D(:))-min(GRID.X_3D(:))) / abs(max(GRID.Z_3D(:))-min(GRID.Z_3D(:)));
            GRID.aspectRatio(2) = abs(max(GRID.Y_3D(:))-min(GRID.Y_3D(:))) / abs(max(GRID.Z_3D(:))-min(GRID.Z_3D(:)));
        end
        if SWITCH.AxesLimit;  GRID.aspectRatio(1) = (SWITCH.AxesLimitValues(1,2)-SWITCH.AxesLimitValues(1,1))/(SWITCH.AxesLimitValues(1,4)-SWITCH.AxesLimitValues(1,3)); end %if manual axis extent
        if GRID.aspectRatio(1)<0.05; GRID.aspectRatio(1) = 1; end %has to be positive and finite
        if GRID.aspectRatio(2)<0.05; GRID.aspectRatio(2) = 1; end %has to be positive and finite
        GRID.rcmb = READ.rcmb;
        
        %% ERROR CHECKS
        if ( max(z2d(:)) > 10*max(x2d(:)) ) && ~(strcmp(GRID.Type,'spherical2D') || SWITCH.DimensionalInput)
            error('z-values is much larger than x-values: did you adjust the dimensionalisation parameters? check this!!!')
        end
        if SWITCH.PlotInPlot
            if strcmp(GRID.Dim,'3-D'); error('SWITCH.PlotInPlot cannot be used with 3-D data!'); end
            %if SAVE.LastFile && ~isfield(SAVE,'AllFilesTitleString'); error('SWITCH.PlotInPlot needs multiple time steps with SAVE.EndNumber > SAVE.StartNumber!'); end
            if SAVE.LastFile && ~isfield(SAVE,'AllFilesTitleString') && nrFiles==1; error('SWITCH.PlotInPlot needs multiple time steps with SAVE.EndNumber > SAVE.StartNumber!'); end
        end
        if SWITCH.plotGraphVsTime
            if strcmp(GRID.Dim,'3-D'); error('SWITCH.plotGraphVsTime cannot be used with 3-D data!'); end
            if SAVE.LastFile && ~isfield(SAVE,'graphTime'); error('SWITCH.plotGraphVsTime needs to plot multiple time steps!'); end
        end
        if strcmp(GRID.Dim,'3-D') && (PLOT.Streamfunction)
            error('Streamfunction plot is not implemented for 3-D geometry: Switch it off!')
        end
        if strcmp(GRID.Type,'yinyang')
            if SWITCH.plotDifference; error('SWITCH.plotDifference not implemented yet for YinYang geometry.'); end
            if SWITCH.Magnifier; SWITCH.Magnifier = false; warning('Yinyang mode: set SWITCH.Magnifier = false'); end %zoom plot not possible in yinyang mode
            if SWITCH.trackC; SWITCH.trackC = false; warning('Yinyang mode: set SWITCH.trackC = false'); end %Compositional tracking not implemented yet for yinyang mode
        end
        
        %% SET UP COLORBAR RANGE
        if  SWITCH.DimensionalMode || SWITCH.DimensionalInput
            dimensional     = true;
        else
            dimensional     = false;
        end
        [CB] = f_DefaultsColorbar(IN,dimensional,loopCase);
        
        %% VARIABLE TRANSITIONS
        PLOT.multipleCBMinMax   = CB.multipleMinMax;
        PLOT.nrPlots            = nrPlots;
        PLOT.nrFiles            = nrFiles;
        PLOT.nrModels           = nrModels;
        PLOT.nrTimeSnapshots    = nrTimeSnapshots;
        PLOT.fileName           = IN.Name;
        FIELD.save              = SAVE.Field;
        
        %% SET UP FIELD VARIABLES
        %  	1                           2                           3                 	4                   5                   6               7           8               9               10              11              12        	13       	14          15
        %  	plot                        name                        symbol              dim                 varscale            colormap        cmap_flip  	cColorbarMIN    cColorbarMAX    logarithmic     zeroCentredCB  	discreteCB	2ndSuiteColormaps       type
        FIELD_M = {
            PLOT.Grid                   'Grid'                      '-'                 GRID.Zdim           GRID.dimFactor      'grayC'         'noflip'	0             	1               false           0               0           'grayC'  	'flip'      'field'
            PLOT.ParameterTable       	'Parameter table'           '-'                 GRID.Zdim           GRID.dimFactor      'grayC'         'noflip'	0             	1               false           0               0           'grayC'     'flip'      'graph'
            PLOT.PlateSketch            'Tectonics'                 '-'                 GRID.Zdim           GRID.dimFactor      'grayC'         'noflip'	0             	1               false           0               0           'grayC'   	'flip'      'graph'
            PLOT.CustomGraph            'Custom graph'              '-'                 PLOT.time2plotDim   PLOT.timeConvert 	'grayC'         'noflip'    CB.graph_min 	CB.graph_max  	false           0               0           'grayC'   	'flip'      'graph'
            PLOT.StreamGraph            'Stream graph'              '-'                 PLOT.time2plotDim   PLOT.timeConvert 	'grayC'         'noflip'    CB.graph_min 	CB.graph_max  	false           0               0           'grayC'    	'flip'      'field'
            PLOT.Topography          	'Topography'                'z_{topo}'          GRID.Zdim           GRID.dimFactor      'broc'          'noflip'    CB.topo_min     CB.topo_max     CB.topo_log     1               0           'lisbon'  	'noflip'    'graph'
            PLOT.IsostaticTopography  	'Iso. topography'           'z_{isotopo}'       GRID.Zdim           GRID.dimFactor   	'broc'      	'noflip'    CB.isotopo_min  CB.isotopo_max	CB.isotopo_log  1               0           'lisbon'  	'noflip'    'graph'
            PLOT.ResidualTopography    	'Res. topography'           'z_{restopo}'       GRID.Zdim           GRID.dimFactor   	'broc'       	'noflip'    CB.restopo_min  CB.restopo_max	CB.restopo_log  1               0           'lisbon'  	'noflip'    'graph'
            PLOT.DynamicTopography    	'Dyn. topography'           'z_{dyntopo}'       GRID.Zdim           GRID.dimFactor   	'broc'        	'noflip'    CB.dyntopo_min  CB.dyntopo_max	CB.dyntopo_log  1               0           'lisbon'   	'noflip'    'graph'
            PLOT.PlateBaseTopography   	'Plate-base topography'   	'z_{topo}'          GRID.Zdim           GRID.dimFactor    	'broc'          'noflip'    CB.isotopo_min  CB.isotopo_max	CB.isotopo_log  1               0           'lisbon'   	'noflip'    'graph'
            PLOT.TopographySelfGrav   	'Topography (self-grav)'  	'z_{topo,sg}'     	GRID.Zdim           GRID.dimFactor      'broc'          'noflip'    CB.topo_min     CB.topo_max     CB.topo_log     1               0           'lisbon'   	'noflip'    'graph'
            PLOT.Geoid                  'Geoid'                     'g'                 GRID.Zdim        	GRID.dimFactor    	'vik'           'noflip'   	CB.geoid_min    CB.geoid_max    CB.geoid_log    1               0           'berlin'   	'noflip'    'graph'
            PLOT.PlateVelocity        	'Plate velocity'            'v_h'             	SETUP.vDim        	SETUP.Vscale       	'cork'       	'noflip'    -CB.vel_max     CB.vel_max    	CB.vel_log      1               0           'tofino'    'noflip'    'graph'
            PLOT.HeatFlux               'Heat flux'                 '\phi'              SETUP.HFDim        	SETUP.HFscale     	'grayC'         'noflip'    0               150            	false           0               0           'grayC'   	'flip'      'graph'
            PLOT.SurfaceAge         	'Surface age'             	'a_{surf}'        	PLOT.time2plotDim  	PLOT.timeConvert  	'turku'         'flip'      0               10            	false           0               0           'turku'   	'flip'      'graph'
            PLOT.CrustalThickness    	'Crustal thickness'       	'd_{crust}'     	GRID.Zdim           GRID.dimFactor      'grayC'         'noflip'    0            	10              false           0               0           'grayC'    	'flip'      'graph'
            PLOT.Tracers               	'Tracers'                  	'-'                 GRID.Zdim           GRID.dimFactor      'batlow'       	'noflip'	0             	1               false           0               0           'batlow'   	'flip'      'field'
            PLOT.Phase                	'Phase'                   	'Ph'                ''                  1                   'tokyo'       	'flip'   	CB.Ph_min       CB.Ph_max     	CB.Ph_log     	0               1           'tokyo'   	'noflip'    'field'
            PLOT.Composition            'Composition'               'C'                 ''                  1                   'davos'         'noflip'   	CB.C_min        CB.C_max     	CB.C_log    	0               0           'davos'    	'noflip'   	'field'
            PLOT.Air                    'Air'                       'C_{air}'         	SETUP.CompDim     	SETUP.Compscale    	'oslo'         	'flip'   	0               100           	false           0               0           'oslo'     	'noflip'    'field'
            PLOT.Basalt                 'Basalt'                    'C_{bs}'         	SETUP.CompDim     	SETUP.Compscale    	'oslo'       	'flip'      0               100            	false           0               0           'oslo'     	'noflip'    'field'
            PLOT.Harzburgite            'Harzburgite'               'C_{hz}'         	SETUP.CompDim    	SETUP.Compscale  	'oslo'        	'flip'      0               100          	false           0               0           'oslo'    	'noflip'    'field'
            PLOT.ContinentalCrust      	'Cont. crust'               'C_{cc}'         	SETUP.CompDim     	SETUP.Compscale    	'oslo'       	'flip'      0               100         	false           0               0           'oslo'     	'noflip'    'field'
            PLOT.Primordial             'Primordial'                'C_{prim}'         	SETUP.CompDim    	SETUP.Compscale   	'oslo'       	'flip'      0               100          	false           0               0           'oslo'     	'noflip'    'field'
            PLOT.Water                  'Water'                     'C_{wtr}'         	SETUP.CompDim     	SETUP.Compscale    	'oslo'         	'flip'      CB.water_min  	CB.water_max  	CB.water_log  	0               0           'oslo'     	'noflip'    'field'
            PLOT.Melt                   'Melt'                      'C_{mlt}'         	''                  1                   'bilbao'       	'noflip'	0               1               false           0               0           'turku'    	'noflip'    'field'
            PLOT.AgeSinceLastMelted   	'Age since melted'         	'a_{solid}'        	PLOT.time2plotDim  	PLOT.timeConvert 	'devon'         'flip'      0               10              false           0               0           'devon'   	'noflip'    'field'
            PLOT.Temperature            'Temperature'               'T'                 SETUP.TDim        	SETUP.Tscale        'lajolla'       'flip'      CB.T_min        CB.T_max        CB.T_log        0               0           'lajolla'	'flip'      'field'
            PLOT.RegionalResidualT      'Regional residual'         'T_{res,R}'       	SETUP.TDim       	SETUP.Tscale      	'vik'           'noflip'    CB.Tres_min     CB.Tres_max     CB.Tres_log     1               0           'berlin'   	'noflip'    'field'
            PLOT.HorizResidualT         'Horizontal residual'       'T_{res,H}'       	SETUP.TDim       	SETUP.Tscale      	'vik'           'noflip'    CB.Tres_min     CB.Tres_max     CB.Tres_log     1               0           'berlin'  	'noflip'    'field'
            PLOT.HorizBandResidualT     'Horizontal-band residual' 	'T_{res,H}'       	SETUP.TDim       	SETUP.Tscale      	'vik'           'noflip'    CB.Tres_min     CB.Tres_max     CB.Tres_log     1               0           'berlin'   	'noflip'    'field'
            PLOT.GlobalResidualT        'Global residual'           'T_{res,G}'       	SETUP.TDim       	SETUP.Tscale      	'vik'           'noflip'    CB.Tres_min     CB.Tres_max     CB.Tres_log     1               0           'berlin'  	'noflip'    'field'
            PLOT.Density                'Density'                   '\rho'           	SETUP.rhoDim       	SETUP.rhoscale    	'bilbao'    	'noflip'  	CB.rho_min      CB.rho_max      CB.rho_log      0               0           'turku'    	'flip'      'field'
            PLOT.Viscosity              'Viscosity'                 '\eta'              SETUP.etaDim      	SETUP.etascale     	'davos'         'noflip'  	CB.eta_min    	CB.eta_max      CB.eta_log      0               0           'davos'  	'noflip'    'field'
            PLOT.Stress                 'Stress'                    '\sigma'            SETUP.stressDim   	SETUP.stressscale 	'acton'        	'flip'  	CB.str_min     	CB.str_max    	CB.str_log      0               0           'acton'   	'noflip'    'field'
            PLOT.StressX                'Horizontal princ. stress'	'\sigma_x'         	SETUP.stressDim   	SETUP.stressscale 	'acton'     	'flip'  	CB.str_min   	CB.str_max    	CB.str_log      0               0           'acton'   	'noflip'    'field'
            PLOT.StressZ                'Radial princ. stress'    	'\sigma_z'        	SETUP.stressDim   	SETUP.stressscale 	'acton'       	'flip'  	CB.str_min   	CB.str_max    	CB.str_log      0               0           'acton'    	'noflip'    'field'
            PLOT.NormalStress           'zz-Stress component'       '\sigma_{zz}'       SETUP.stressDim   	SETUP.stressscale  	'vik'           'noflip'    CB.nstr_min     CB.nstr_max     false           1               0           'berlin'   	'noflip'    'field'
            PLOT.StrainRate             'Strain rate'               'UCvarepsilondot'   SETUP.edotDim   	SETUP.edotscale    	'bilbao'        'noflip'  	CB.edot_min     CB.edot_max     CB.edot_log     0               0           'turku'    	'noflip'    'field'
            PLOT.Velocity               'Velocity'                  'v'                 SETUP.vDim        	SETUP.Vscale       	'osloB'        	'flip'      CB.vel_min      CB.vel_max      CB.vel_log      0               0           'oslo'    	'noflip'    'field'
            PLOT.VelocityX              'Horizontal velocity'       'v_x'               SETUP.vDim         	SETUP.Vscale       	'cork'          'noflip'	-CB.vel_max     CB.vel_max      CB.vel_log      1               0           'tofino'  	'noflip'    'field'
            PLOT.VelocityZ              'Radial velocity'           'v_z'               SETUP.vDim         	SETUP.Vscale      	'cork'          'noflip'	-CB.vel_max     CB.vel_max      CB.vel_log      1               0           'tofino'   	'noflip'    'field'
            PLOT.Pressure               'Pressure'                  'P'                 SETUP.PDim        	SETUP.Pscale      	'bilbao'    	'noflip'  	CB.P_min        CB.P_max        CB.P_log        0               0           'turku'   	'noflip'    'field'
            PLOT.DynamicPressure      	'Dynamic pressure'        	'P_{dyn}'        	SETUP.PDim        	SETUP.Pscale      	'vik'           'noflip'  	CB.P_min        CB.P_max        CB.P_log        1               0           'berlin'  	'noflip'    'field'
            PLOT.Toroidal               'Toroidal'                  'v_{tor}'           ''                  1                   'vik'           'noflip'    0               1               false           1               0           'berlin'   	'noflip'    'field'
            PLOT.Poloidal               'Poloidal'                  'v_{pol}'         	''                  1                   'vik'           'noflip'    0               1               false           1               0           'berlin'   	'noflip'    'field'
            PLOT.Streamfunction       	'Streamfunction'            '\psi'              SETUP.streamFDim    SETUP.streamFscale 	'cork'          'noflip'    CB.sfun_min     CB.sfun_max   	CB.sfun_log   	1               0           'tofino'   	'noflip'    'field'
            PLOT.DeformationMechanism 	'Deformation mechanism'     'M_{def}'           'mechanism'      	1                   'batlow'       	'noflip'  	1               4               false           0               0           'batlow'   	'noflip'  	'field'
            PLOT.ViscousDissipation    	'Viscous dissipation'       '\phi_{vd}'         SETUP.dissDim       SETUP.dissscale     'bilbao'        'noflip'   	CB.diss_min   	CB.diss_max   	CB.diss_log   	0               0           'turku'   	'noflip'    'field'
            PLOT.UpDownWelling          'Upwelling and downwelling'	'-'                 ''                  1                   'cork'       	'noflip'    -2              2               false           1               1           'tofino' 	'noflip'    'field'
            PLOT.SurfaceFieldVariation 	'Surface field variation' 	'-'                 GRID.Zdim           GRID.dimFactor      'grayC'         'noflip'	CB.sfv_min      CB.sfv_max  	CB.sfv_log     	0               0           'grayC'   	'flip'      'graph'
            };
        
        %% UNICODE VARIABLES
        %use the following character sequence to add unicode variables:
        %UCvarepsilondot => char(941)
        %UCvarepsilon => char(949)
        FIELD_M(:,3) = strrep(FIELD_M(:,3),'UCvarepsilondot',{char(941)}); %does only work if string consists ONLY of unicode character
        FIELD_M(:,3) = strrep(FIELD_M(:,3),'UCvarepsilon',{char(949)}); %does only work if string consists ONLY of unicode character
        
        %% PREVENTING VISUALISATION PITFALLS
        DESVARIA.FIELD_M    = FIELD_M;
        DESVARIA.Task       = 'prevent visualisation pitfalls';
        [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
        FIELD_M             = DESVARIA.FIELD_M;
        
        %% WAITBAR SWITCH
        if SWITCH.closeOldWaitbars; wbOld = findall(0,'tag','TMWWaitbar'); delete(wbOld); end %close open ghost waitbars
        if strcmp(GRID.Type,'yinyang') && logical(0) %switch on waitbar
            SWITCH.waitbar  = true;
            disp('   ...reading')
            if isfield(PLOT,'wb'); PLOT.wb = waitbar(0,PLOT.wb,'reading...'); else; PLOT.wb = waitbar(0,'reading...'); end
        else
            SWITCH.waitbar  = false;
        end
        
        %% PLATE DIAGNOSTICS
        if SWITCH.PlateDiagnostics
            % GET PLATE TECTONIC DIAGNOSTICS
            [PLATE,PLOT] = f_DiagnosticsPlate(FILE,GRID,SETUP,SWITCH,PLOT,STYLE);
            %save for later
            PLOT.IndicationTrench   = PLATE.Subduction;
            if isfield(PLATE,'SubPolarity'); PLOT.SubPolarity = PLATE.SubPolarity; end
            if isfield(PLATE,'Spreading'); PLOT.IndicationRidge = PLATE.Spreading; end
            TOPO.UPtiltAngleXnear   = PLATE.UPtiltAngleXnear; TOPO.UPtiltAngleXfar   = PLATE.UPtiltAngleXfar;
            TOPO.UPtiltAngleZnear   = PLATE.UPtiltAngleZnear; TOPO.UPtiltAngleZfar   = PLATE.UPtiltAngleZfar;
        else
            PLATE.dummy     = NaN; %just to be defined
        end
        
        %% MANTLE DIAGNOSTICS
        if SWITCH.MantleDiagnostics
            % GET MANTLE FLOW DIAGNOSTICS
            [MANTLE] = f_DiagnosticsMantle(FILE,GRID,SETUP,SWITCH,PLOT,STYLE);
        else
            MANTLE.dummy    = NaN; %just to be defined
        end
        
        %% SAVING DATA
        % PLATE VELOCITY (VECTOR) DATA
        if SAVE.plateVelocityData
            SAVE.Directory              = FILE.directory;
            SAVE.DataName               = [FILE.name,'_platevel',num2str(FILE.number)];
            SAVE.data                   = [PLATE.PlateVelocity(:,1), GRID.x2d(:,1)]; %plate-V, x
            SAVE.dat                    = false;
            SAVE.txt                    = false;
            SAVE.mat                    = true;
            SAVE.write2existing         = false;
            [SAVE.overwriteAll] = f_saveData( SAVE );
            SAVE = rmfield(SAVE,'data');
        end
        % GEODYNAMIC DATA
        %ADD:
        %THEORETIC PLATE AGE
        %WATER RETENTION
        if SAVE.GeodynamicDiagnostics
            if strcmp(GRID.Type,'spherical2D')
                DistOrRad   = 'rad';
            else
                DistOrRad   = GRID.Xdim;
            end
            SAVE.Directory              = FILE.directory;
            SAVE.DataName               = [FILE.name,'_Geodynamics',num2str(FILE.number)];
            SAVE.data                   = {
                %name                           %group name             %dimension         	%dimensionalisation
                'Time',                         'Time',                 SETUP.timeDim,   	PLOT.time                       %can be converted by PLOT.timeConvert to PLOT.time2plotDim
                'Time seconds',             	'Time',                 's',                PLOT.time_dim*SETUP.secyear
                'Time years',                	'Time',                 'a',              	PLOT.time_dim
                'Mantle-transit time',      	'Time',                 '#MT',            	PLOT.timeMT
                };
            if exist('PLATE','var') && isfield(PLATE,'Subduction')
                if isnan(PLATE.SlabTipPosition); PLATE.SlabTipPosition(:,2) = NaN; end
                SAVE.data                   = [SAVE.data;
                    {
                    'Trench position',              'Position',             GRID.Xdim,       	PLATE.Subduction(:,1)'
                    'Subduction polarity',          'Subduction polarity', 	'',              	PLATE.SubPolarityNum(:,1)'      %polarity: '-1': left, '1': right, or '0': unknown
                    'Trench velocity',              'Velocity',             SETUP.vDim,      	PLATE.TrenchVelocity(:,1)'      %[cm/a]
                    'Theoretic trench velocity',    'Velocity',             SETUP.vDim,      	PLATE.theoreticTrenchVelocity(:,1)'	%[cm/a]
                    'Upper-plate velocity',         'Velocity',             SETUP.vDim,       	PLATE.UPvelocity(:,1)'         	%[cm/a]
                    'Lower-plate velocity',         'Velocity',             SETUP.vDim,      	PLATE.LPvelocity(:,1)'          %[cm/a]
                    'Convergence velocity',         'Velocity',             SETUP.vDim,       	PLATE.PlateConvergence(:,1)'    %[cm/a]
                    'Slab sinking velocity',       	'Velocity',             SETUP.vDim,        	PLATE.SlabSinkingVelocity(:,1)'	%[cm/a]
                    'Max. plate velocity',          'Velocity',             SETUP.vDim,        	PLATE.PlateVelocityMax          %[cm/a]
                    'Plate velocity difference',   	'Velocity',             SETUP.vDim,        	PLATE.PlateVelocityDiff        	%[cm/a]
                    'RMS plate velocity',           'Velocity',             SETUP.vDim,        	PLATE.PlateVelocityRMS          %[cm/a]
                    'Plate mobility',               'Dynamics',             '',                 PLATE.PlateMobility             %[nd]
                    'Plate drivety',                'Dynamics',             '',                 PLATE.PlateDrivety              %[nd] - PLATE.PlateVelocityRMS/PLATE.UMantleHVelocityMean  - is plate driven by itself or by the mantle?
                    'Slab-tip horiz. position',   	'Slab-tip position',   	GRID.Xdim,       	PLATE.SlabTipPosition(:,1)'     %[plotting dimension] or [degree]
                    'Slab-tip depth',               'Slab-tip position',   	GRID.Zdim,       	PLATE.SlabTipPosition(:,2)'     %[plotting dimension]
                    'Slab-tip angle',               'Slab angle',           '{\circ}',          PLATE.SlabTipAngle(:,1)'        %[degree]
                    'Slab-tip horiz. velocity',     'Velocity',             SETUP.vDim,        	PLATE.SlabTipVX(:,1)'           %[cm/a]
                    'Slab-tip sinking velocity',    'Velocity',             SETUP.vDim,        	PLATE.SlabTipVZ(:,1)'           %[cm/a]
                    'Shallow-slab angle',         	'Slab angle',           '{\circ}',          PLATE.ShallowSlabAngle(:,1)'   	%[degree]
                    'Slab viscosity',               'Viscosity',            SETUP.etaDim,    	PLATE.SlabViscosity(:,1)'
                    'Slab density',                 'Density',              SETUP.rhoDim,       PLATE.SlabDensityMax(:,1)'
                    'Upper-mantle viscosity',       'Viscosity',            SETUP.etaDim,    	PLATE.UMantleViscosity
                    'Upper-mantle density',         'Density',              SETUP.rhoDim,    	PLATE.UMantleDensity
                    'Upper-mantle max. velocity', 	'Velocity',          	SETUP.vDim,         PLATE.UMantleVelocityMax        %[cm/a]
                    'UM max. horiz. velocity',      'Velocity',          	SETUP.vDim,         PLATE.UMantleHVelocityMax     	%[cm/a]
                    'UM max. radial velocity',      'Velocity',          	SETUP.vDim,         PLATE.UMantleRVelocityMax    	%[cm/a]
                    'Slab-mantle visc. contrast',   'Viscosity contrast', 	'',              	PLATE.SlabMantleViscDiff(:,1)'
                    'Left-plate thickness',         'Thickness',            GRID.Zdim,       	PLATE.PlateThicknessL(:,1)'
                    'Right-plate thickness',        'Thickness',            GRID.Zdim,        	PLATE.PlateThicknessR(:,1)'
                    'Lower-plate thickness',        'Thickness',            GRID.Zdim,       	PLATE.PlateThicknessLP(:,1)'
                    'Upper-plate thickness',        'Thickness',            GRID.Zdim,        	PLATE.PlateThicknessUP(:,1)'    %[plotting dimension]
                    'Plate bending radius',         'Radius',               GRID.Zdim,        	PLATE.BendingRadius(:,1)'       %[plotting dimension]
                    'Bending dissipation',          'Dissipation',          'W/m',              PLATE.ViscDissBending(:,1)'     %[N/s]=[W/m]
                    'Rel. bending dissipation', 	'Rel. dissipation',     '',               	PLATE.ViscDissBendingRel(:,1)'
                    'Viscous plate dissipation',   	'Dissipation',          'W/m',              PLATE.ViscDissipationPlate    	%[W/m] %within 1300-K isovolume
                    'Viscous plate dissipation*',  	'Dissipation',          'W/m',              PLATE.ViscDissipationPlateNoCrust %[W/m] %without crustal areas
                    'Max. plate-core viscosity',  	'Viscosity',           	SETUP.etaDim,       PLATE.coreViscosityMax
                    'Min. plate-core strain rate',	'Strain rate',        	SETUP.edotDim,  	PLATE.coreStrainrateMin
                    'Max. plate-core strain rate', 	'Strain rate',        	SETUP.edotDim,  	PLATE.coreStrainrateMax
                    'Max. plate-core stress',     	'Stress',               SETUP.stressDim,  	PLATE.coreStressMax
                    'Max. plate stress',            'Stress',               SETUP.stressDim,  	PLATE.maxStress
                    'LAB depth',                    'Depth',                GRID.Zdim,       	PLATE.LABdepth
                    'Max. yield depth',             'Depth',                GRID.Zdim,       	PLATE.maxYieldDepth
                    'Max. yield depth fraction',    'Depth fraction',       '',                 PLATE.maxYieldDepthFraction     %fraction of plate thickness
                    'Trench depth',                 'Surface elevation',   	GRID.Zdim,          PLATE.TrenchDepth(:,1)'       	%[plotting dimension]
                    'Trench depth X',             	'Surface elevation',   	GRID.Xdim,          PLATE.TrenchDepthX(:,1)'       	%[plotting dimension]
                    'Back-arc basin position',    	'Surface elevation',   	GRID.Xdim,          PLATE.BackArcBasinX(:,1)'      	%[plotting dimension]
                    'Back-arc basin depth',        	'Surface elevation',   	GRID.Zdim,          PLATE.BackArcBasinZ(:,1)'     	%[plotting dimension]
                    'Back-arc basin extent',    	'Surface elevation',   	GRID.Xdim,          PLATE.BackArcBasinExtentX(:,1)'	%[plotting dimension]
                    'Inundation',                   'Surface elevation',   	GRID.Xdim,          PLATE.InundationDistance(:,1)'	%[plotting dimension]
                    'Island-arc position',         	'Surface elevation',   	GRID.Xdim,          PLATE.IslandArcX(:,1)'       	%[plotting dimension]
                    'Island-arc height',         	'Surface elevation',   	GRID.Zdim,          PLATE.IslandArcZ(:,1)'       	%[plotting dimension]
                    'Fore-bulge position',       	'Surface elevation',   	GRID.Xdim,          PLATE.ForeBulgeX(:,1)'       	%[plotting dimension]
                    'Fore-bulge height',           	'Surface elevation',   	GRID.Zdim,          PLATE.ForeBulgeZ(:,1)'       	%[plotting dimension]
                    'Back-arc basin volume',    	'Surface elevation',   	[GRID.Zdim,'^2'], 	PLATE.BackArcBasinArea(:,1)'  	%[plotting dimension]
                    'Upper-plate tilt',            	'Upper-plate tilt',    	'{\circ}',          PLATE.UPtiltAngle(:,1)'       	%[degree]
                    'Subduction-flow rate',        	'Flow rate',          	'm^2/s',            PLATE.subductionFlowRate      	%[m^2/s] or nd (volumetric flow rate)
                    'Mean surface age',             'Surface age',        	'Ma',               PLATE.SurfaceAgeMean            %[Ma]
                    'Median surface age',        	'Surface age',        	'Ma',               PLATE.SurfaceAgeMedian          %[Ma]
                    'Mean cont. surface age',    	'Surface age',        	'Ma',               PLATE.SurfaceAgeMeanCont     	%[Ma]
                    'Median cont. surface age',   	'Surface age',        	'Ma',               PLATE.SurfaceAgeMedianCont     	%[Ma]
                    'Mean ocean surface age',     	'Surface age',        	'Ma',               PLATE.SurfaceAgeMeanOcean      	%[Ma]
                    'Median ocean surface age',    	'Surface age',        	'Ma',               PLATE.SurfaceAgeMedianOcean    	%[Ma]
                    'Mean crustal thickness',      	'Crust',                GRID.Zdim,          PLATE.CrustalThicknessMean    	%[plotting dimension]
                    'Median crustal thickness',  	'Crust',                GRID.Zdim,          PLATE.CrustalThicknessMedian  	%[plotting dimension]
                    'Mean cont. crustal thickness',  	'Crust',        	GRID.Zdim,          PLATE.CrustalThicknessMeanCont    	%[plotting dimension]
                    'Median cont. crustal thickness',  	'Crust',            GRID.Zdim,          PLATE.CrustalThicknessMedianCont  	%[plotting dimension]
                    'Mean oceanic crustal thickness',  	'Crust',            GRID.Zdim,          PLATE.CrustalThicknessMeanOcean    	%[plotting dimension]
                    'Median oceanic crustal thickness',	'Crust',          	GRID.Zdim,          PLATE.CrustalThicknessMedianOcean  	%[plotting dimension]
                    }];
            end
            if exist('MANTLE','var') && ~isfield(MANTLE,'dummy')
                SAVE.data                   = [SAVE.data;
                    {
                    'v_{h,max} UM',             'Flow velocity',     	SETUP.vDim,       	MANTLE.VHmaxUM                  %[cm/a]
                    'v_{h,max} MM',             'Flow velocity',     	SETUP.vDim,       	MANTLE.VHmaxMM                  %[cm/a]
                    'v_{h,max} LM',             'Flow velocity',     	SETUP.vDim,       	MANTLE.VHmaxLM                  %[cm/a]
                    'v_{h,mean} UM',          	'Flow velocity',     	SETUP.vDim,       	MANTLE.VHmeanUM              	%[cm/a]
                    'v_{h,mean} MM',        	'Flow velocity',     	SETUP.vDim,       	MANTLE.VHmeanMM              	%[cm/a]
                    'v_{h,mean} LM',          	'Flow velocity',     	SETUP.vDim,       	MANTLE.VHmeanLM              	%[cm/a]
                    'Number slabs',             'Number slabs',         '#',                MANTLE.numColdPlumes          	%[number count]
                    'Upper-mantle slabs',       'Number slabs',     	'#',                MANTLE.plumeColdNumberUM      	%[number count]
                    'Mid-mantle slabs',         'Number slabs',     	'#',                MANTLE.plumeColdNumberMM      	%[number count]
                    'Lower-mantle slabs',       'Number slabs',     	'#',                MANTLE.plumeColdNumberLM      	%[number count]
                    'Active downwelling volume','Slab dynamics',     	SETUP.volumeDimPlot,MANTLE.ActDownwellingVolume 	%[nd] or [m^2]or[m^3]
                    'Active part. downwelling',	'Slab dynamics',     	'%',                MANTLE.ActDownwellingVolPerc 	%[% of total downwelling]
                    'v_{h,max} UM slabs',    	'Slab dynamics',     	SETUP.vDim,       	MANTLE.plumeColdVHmaxUM      	%[cm/a]
                    'v_{h,min} UM slabs',      	'Slab dynamics',     	SETUP.vDim,       	MANTLE.plumeColdVHminUM      	%[cm/a]
                    'v_{h,mean} UM slabs',    	'Slab dynamics',     	SETUP.vDim,       	MANTLE.plumeColdVHmeanUM      	%[cm/a]
                    'v_{h,max} MM slabs',      	'Slab dynamics',     	SETUP.vDim,       	MANTLE.plumeColdVHmaxMM      	%[cm/a]
                    'v_{h,min} MM slabs',     	'Slab dynamics',     	SETUP.vDim,       	MANTLE.plumeColdVHminMM      	%[cm/a]
                    'v_{h,mean} MM slabs',     	'Slab dynamics',     	SETUP.vDim,       	MANTLE.plumeColdVHmeanMM      	%[cm/a]
                    'v_{h,max} LM slabs',     	'Slab dynamics',     	SETUP.vDim,       	MANTLE.plumeColdVHmaxLM      	%[cm/a]
                    'v_{h,min} LM slabs',      	'Slab dynamics',     	SETUP.vDim,       	MANTLE.plumeColdVHminLM      	%[cm/a]
                    'v_{h,mean} LM slabs',  	'Slab dynamics',     	SETUP.vDim,       	MANTLE.plumeColdVHmeanLM      	%[cm/a]
                    'Number plumes',        	'Number plumes',     	'#',                MANTLE.numHotPlumes          	%[number count]
                    'Upper-mantle plumes',  	'Number plumes',     	'#',                MANTLE.plumeHotNumberUM      	%[number count]
                    'Mid-mantle plumes',        'Number plumes',     	'#',                MANTLE.plumeHotNumberMM      	%[number count]
                    'Lower-mantle plumes',  	'Number plumes',     	'#',                MANTLE.plumeHotNumberLM      	%[number count]
                    'Active upwelling volume',  'Plume dynamics',     	SETUP.volumeDimPlot,MANTLE.ActUpwellingVolume       %[nd] or [m^2]or[m^3]
                    'Active part. upwelling',	'Plume dynamics',     	'%',                MANTLE.ActUpwellingVolPerc      %[% of total upwelling]
                    'v_{h,max} UM plumes',    	'Plume dynamics',     	SETUP.vDim,       	MANTLE.plumeHotVHmaxUM      	%[cm/a]
                    'v_{h,min} UM plumes',      'Plume dynamics',     	SETUP.vDim,       	MANTLE.plumeHotVHminUM      	%[cm/a]
                    'v_{h,mean} UM plumes',    	'Plume dynamics',     	SETUP.vDim,       	MANTLE.plumeHotVHmeanUM      	%[cm/a]
                    'v_{h,max} MM plumes',     	'Plume dynamics',     	SETUP.vDim,       	MANTLE.plumeHotVHmaxMM      	%[cm/a]
                    'v_{h,min} MM plumes',     	'Plume dynamics',     	SETUP.vDim,       	MANTLE.plumeHotVHminMM      	%[cm/a]
                    'v_{h,mean} MM plumes',   	'Plume dynamics',     	SETUP.vDim,       	MANTLE.plumeHotVHmeanMM      	%[cm/a]
                    'v_{h,max} LM plumes',     	'Plume dynamics',     	SETUP.vDim,       	MANTLE.plumeHotVHmaxLM      	%[cm/a]
                    'v_{h,min} LM plumes',    	'Plume dynamics',     	SETUP.vDim,       	MANTLE.plumeHotVHminLM      	%[cm/a]
                    'v_{h,mean} LM plumes',  	'Plume dynamics',     	SETUP.vDim,       	MANTLE.plumeHotVHmeanLM      	%[cm/a]
                    'Continents',               'Number continents',  	'#',                MANTLE.contNumberUM             %[number count]
                    'v_{h,max} Continent',      'Continent dynamics',  	SETUP.vDim,       	MANTLE.contVHmaxUM              %[cm/a]
                    'v_{h,min} Continent',    	'Continent dynamics',  	SETUP.vDim,       	MANTLE.contVHminUM              %[cm/a]
                    'v_{h,mean} Continent',  	'Continent dynamics', 	SETUP.vDim,       	MANTLE.contVHmeanUM             %[cm/a]
                    'Continent horiz. locations','Continent location', 	DistOrRad,          MANTLE.contLocationX            %[plotting dimension]
                    'Continent vert. locations', 'Continent location', 	GRID.Zdim,          MANTLE.contLocationZ            %[plotting dimension]
                    'Lower-mantle LLSVPs',      'Number LLSVPs',     	'#',                MANTLE.llsvpNumberLM            %[number count]
                    'v_{h,max} LM LLSVPs',     	'LLSVP dynamics',     	SETUP.vDim,       	MANTLE.llsvpVHmaxLM             %[cm/a]
                    'v_{h,min} LM LLSVPs',    	'LLSVP dynamics',     	SETUP.vDim,       	MANTLE.llsvpVHminLM             %[cm/a]
                    'v_{h,mean} LM LLSVPs',  	'LLSVP dynamics',     	SETUP.vDim,       	MANTLE.llsvpVHmeanLM            %[cm/a]
                    'LLSVP horiz. locations',  	'LLSVP location',    	DistOrRad,          MANTLE.llsvpLocationX           %[plotting dimension]
                    'LLSVP vert. locations',   	'LLSVP location',     	GRID.Zdim,          MANTLE.llsvpLocationZ           %[plotting dimension]
                    }];
            end
            SAVE.dat                    = false;
            SAVE.txt                    = false;
            SAVE.mat                    = true;
            SAVE.write2existing         = false;
            [SAVE.overwriteAll] = f_saveData( SAVE );
            SAVE = rmfield(SAVE,'data');
        elseif SAVE.GeodynamicDiagnostics
            warning('No geodynamic data found - No data saved.')
        end
        % TIME DATA
        if SAVE.timeDat || SAVE.subductionZoneData || SAVE.plateVelocityData || SAVE.Tracers
            SAVE.Directory              = FILE.directory;
            SAVE.DataName               = [FILE.name,'_time',num2str(FILE.number)];
            SAVE.data                   = [PLOT.time, PLOT.time_dim*SETUP.secyear,PLOT.time_dim]; %[nd, seconds, years]
            SAVE.dat                    = false;
            SAVE.txt                    = false;
            SAVE.mat                    = true;
            SAVE.write2existing         = false;
            
            [SAVE.overwriteAll] = f_saveData( SAVE );
            SAVE = rmfield(SAVE,'data');
        end
        
        %% VARIABLE CONVERSION AND CLEARANCE
        SWITCH.Colorbar     = [ SWITCH.Colorbar    SWITCH.ConstantColorbar ];
        clearvars READ CB X_3D Y_3D Z_3D x2d y2d z2d
        
        %% INITIALISE VARIABLES
        TOPO.varScale = FIELD_M{strcmp(FIELD_M(:,2),'Topography'),5};
        
        PLOT.loopField = 0; FILE.loop_field = PLOT.loopField;
        if ~isfield(SAVE,'GraphPlotHandles'); SAVE.GraphPlotHandles = []; end
        if ~isfield(SAVE,'PlotHandles'); SAVE.PlotHandles = []; end
        
        %% LOOP PARAMETER FIELDS
        for fieldLoop=1:size(FIELD_M,1)
            if SWITCH.waitbar; waitbar(fieldLoop/size(FIELD_M,1)*0.4,PLOT.wb); end
            if FIELD_M{fieldLoop,1}  %check if plotting current field
                if PLOT.loopPlot==0 || ~SWITCH.PlotInPlot || (SWITCH.PlotInPlot && ~fileRepetition)
                    PLOT.loopPlot = PLOT.loopPlot+1; %track no. of PLOTS plots per file
                end
                %setup field parameters
                FIELD.name          = FIELD_M{fieldLoop,2};
                FIELD.symbol        = FIELD_M{fieldLoop,3};
                FIELD.dim           = FIELD_M{fieldLoop,4};
                FIELD.varScale      = FIELD_M{fieldLoop,5};
                if SWITCH.multipleColormaps
                    if SWITCH.MultiManualColourMaps
                        SWITCH.ColormapAll  = [PLOT.MultiManualColourMaps(min(PLOT.loopPlot,size(PLOT.MultiManualColourMaps,1)),1), ...
                            PLOT.MultiManualColourMaps(min(PLOT.loopPlot,size(PLOT.MultiManualColourMaps,1)),2)];     %multiple manual colormaps
                    else
                        if SWITCH.AlternativeColourmaps || strcmp(STYLE.ColorMode,'dark')
                            SWITCH.ColormapAll  = [FIELD_M(fieldLoop,13)  FIELD_M(fieldLoop,14)];	%alternative-suite colormaps
                        else
                            SWITCH.ColormapAll  = [FIELD_M(fieldLoop,6)  FIELD_M(fieldLoop,7)];     %standard-suite colormaps
                        end
                    end
                end
                FIELD.cColorbarMIN	= FIELD_M{fieldLoop,8};
                FIELD.cColorbarMAX	= FIELD_M{fieldLoop,9};
                FIELD.logarithmic	= FIELD_M{fieldLoop,10};
                FIELD.zeroCentredCB	= FIELD_M{fieldLoop,11};
                FIELD.discreteCB	= FIELD_M{fieldLoop,12};
                FIELD.plotType      = FIELD_M{fieldLoop,15};
                
                if strcmp(FIELD.name,'Tectonics') || strcmp(FIELD.name,'Grid') || strcmp(FIELD.name,'Parameter table') || ...  %check if special field
                        strcmp(FIELD.name,'Surface field variation') || strcmp(FIELD.name,'Tracers') || strcmp(FIELD.name,'Stream graph')
                    dummy = 1;
                    if strcmp(FIELD.name,'Tracers');        dummy = length(PLOT.TracerVariable); end
                    if strcmp(FIELD.name,'Stream graph');   dummy = length(PLOT.StreamGraphParameter); end
                    for iPlotLocal=1:dummy %allow for multiple plots
                        if strcmp(GRID.Type,'yinyang'); error('FieldSpecial is not yet implemented for YinYang grid.'); end
                        [PLOT,SAVE] = f_plotFieldSpecial(iPlotLocal,FIELD_M,FIELD,FILE,GRID,SETUP,SWITCH,PLOT,STYLE,SAVE,PLATE);
                        PLOT.loopPlot = PLOT.loopPlot+1; %track no. of PLOTS plots per file
                    end
                    if isfield(PLOT,'TraData'); PLOT = rmfield(PLOT,'TraData'); end
                    continue %jump to next loop
                end
                if strcmp(FIELD.name,'Custom graph') %check if temporal graph
                    for iPlotLocal=1:size(PLOT.CustomGraphName,2) %allow for multiple graphs
                        if strcmp(GRID.Type,'yinyang'); warning('Custom graph still in testing phase for 3-D grids.'); end
                        [PLOT,SAVE] = f_plotCustomGraph(iPlotLocal,FIELD,GRID,SWITCH,SETUP,PLOT,FILE,STYLE,SAVE);
                        PLOT.loopPlot = PLOT.loopPlot+1; %track no. of PLOTS plots per file
                    end
                    continue %jump to next loop
                end
                if strcmp(FIELD.name,'Topography') || (TOPO.field && ~isfield(TOPO,'topo2d')) || ...  %check if topography
                        strcmp(FIELD.name,'Dyn. topography') || strcmp(FIELD.name,'Iso. topography') || strcmp(FIELD.name,'Res. topography')
                    if ~strcmp(GRID.Type,'yinyang') %ignore that for yy-grid
                        [TOPO,PLOT,SAVE] = f_SetupTopography(TOPO,FIELD,GRID,SWITCH,SETUP,PLOT,FILE,STYLE,SAVE,PLATE);
                    end
                end
                if strcmp(FIELD.name,'Heat flux') || strcmp(FIELD.name,'Surface age') || strcmp(FIELD.name,'Crustal thickness') || ... %check if 2-D field
                        strcmp(FIELD.name,'Geoid')
                    [FIELD,PLOT,SAVE] = f_Setup2DField(FIELD,GRID,SWITCH,SETUP,PLOT,FILE,SAVE);
                end
                if strcmp(FIELD.name,'Plate velocity') || strcmp(FIELD.name,'Topography') || strcmp(FIELD.name,'Geoid') || ...  %plot graphs
                        strcmp(FIELD.name,'Dyn. topography') || strcmp(FIELD.name,'Iso. topography') || strcmp(FIELD.name,'Res. topography') || ...
                        strcmp(FIELD.name,'Heat flux') || strcmp(FIELD.name,'Surface age') || strcmp(FIELD.name,'Crustal thickness') || ...
                        strcmp(FIELD.name,'Plate-base topography') || strcmp(FIELD.name,'Topography (self-grav)')
                    if ~strcmp(GRID.Type,'yinyang') %ignore plotting topo and 2-D fields here for yy-grid
                        [TOPO,PLOT,SAVE] = f_plotGraph(TOPO,PLATE,FIELD,GRID,SWITCH,SETUP,PLOT,FILE,STYLE,SAVE);
                        continue; %jump to next loop
                    else
                        if strcmp(FIELD.name,'Plate velocity'); error('Plate Velocity Plot not yet implemented for YinYang grid.'); end
                        if strcmp(FIELD.name,'Plate-base topography'); error('Plate-Base Topography Plot not yet implemented for YinYang grid.'); end
                    end
                end
                
                PLOT.loopField = PLOT.loopField+1; FILE.loop_field = PLOT.loopField; %track no. of PARAMETER FIELDS plots per file
                if ~strcmp(FIELD.name,'Velocity') &&...
                        ~strcmp(FIELD.name,'Horizontal velocity') &&...
                        ~strcmp(FIELD.name,'Radial velocity') &&...
                        ~strcmp(FIELD.name,'Pressure') &&...
                        ~strcmp(FIELD.name,'Horizontal princ. stress') &&...
                        ~strcmp(FIELD.name,'Radial princ. stress') &&...
                        ~strcmp(FIELD.name,'Toroidal')  &&...
                        ~strcmp(FIELD.name,'Poloidal')  &&...
                        ~strcmp(FIELD.name,'Streamfunction')
                    
                    if (strcmp(FIELD.name,'Horizontal residual') && ~strcmp(GRID.Type,'yinyang')) ||... %res T for yy is done in f_plotField
                            (strcmp(FIELD.name,'Horizontal-band residual') && ~strcmp(GRID.Type,'yinyang')) ||...
                            (strcmp(FIELD.name,'Global residual') && ~strcmp(GRID.Type,'yinyang')) ||...
                            (strcmp(FIELD.name,'Regional residual') && ~strcmp(GRID.Type,'yinyang')) ||...
                            strcmp(FIELD.name,'Upwelling and downwelling') || strcmp(FIELD.name,'Dynamic pressure') ||...
                            strcmp(FIELD.name,'Viscous dissipation')
                        %% CREATE CUSTOM SCALAR FIELD -------------------------------------------------
                        [VAR] = f_makeField(FIELD,FILE,GRID,SETUP,SWITCH,PLOT,MANTLE,TOPO);
                        
                    else
                        %% READ SCALAR FIELD ---------------------------------------------------
                        if strcmp(GRID.Type,'yinyang')
                            if (strcmp(FIELD.name,'Horizontal residual')) ||... %res T for yy is done in f_plotField
                                    (strcmp(FIELD.name,'Horizontal-band residual')) ||...
                                    (strcmp(FIELD.name,'Global residual')) ||...
                                    (strcmp(FIELD.name,'Regional residual'))
                                [~,~,~,~,VAR_3Dyin,VAR_3Dyang,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Temperature',SWITCH); %should be moved into f_makeField......
                            else
                                [~,~,~,~,VAR_3Dyin,VAR_3Dyang,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,FIELD.name,SWITCH);
                            end
                            if ischar(VAR_3Dyin); [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE); warning off backtrace; warning(VAR_3Dyin); warning on backtrace; continue; end
                            VAR.var3d_yin   = VAR_3Dyin;
                            VAR.var3d_yang  = VAR_3Dyang;
                            
                        else %all other geometries
                            [~,~,~,~,VAR_3D,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,FIELD.name,SWITCH);
                            if ischar(VAR_3D); [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE); warning off backtrace; warning(VAR_3D); warning on backtrace; continue; end
                            if size(VAR_3D,1)==1 %exchange x and y
                                dummy_3D = zeros(size(VAR_3D,2),size(VAR_3D,1),size(VAR_3D,3));
                                dummy_3D(:,1,:) = VAR_3D(1,:,:); VAR_3D = dummy_3D;
                            end
                            if strcmp(GRID.Dim,'3-D')
                                var2d       = VAR_3D(:,:,:);
                            elseif strcmp(GRID.Dim,'2-D')
                                var2d(:,:)  = VAR_3D(:,1,:);
                            end
                            if strcmp(FIELD.name,'Deformation mechanism')
                                if strcmp(GRID.Dim,'3-D')
                                    disp('WARNING: def.mech. this needs to be adjusted here for 3-D plots! 3DXXX')
                                elseif strcmp(GRID.Dim,'2-D')
                                    %make sure values are the exact numbers they should be down to 0.00x:
                                    var2d(:,:)  = round(var2d*1e3)/1e3;
                                    if var2d(2,2)==1 && var2d(3,2)==4 %check for artificial values
                                        var2d(2,2) = var2d(4,2);  %and remove them
                                        var2d(3,2) = var2d(4,2);
                                    end
                                end
                            end
                            if SWITCH.plotDifference && mod(loopCase,2) %only for odd cases
                                if loopCase>2; error('more than 1 field not yet implemented!'); end
                                var2d_A     = var2d;
                                continue
                            elseif SWITCH.plotDifference && ~mod(loopCase,2) %only for even cases
                                VAR.var2d   = var2d_A-var2d;  %take difference
                                clearvars var2d_A
                            else
                                VAR.var2d   = var2d;
                            end
                        end
                    end
                    [PLOT,SAVE] =  f_plotField( VAR,FIELD,FIELD_M,GRID,SETUP,SWITCH,PLOT,TOPO,FILE,STYLE,SAVE,PLATE,MANTLE );
                    clearvars VAR var2d VAR_3D VAR_3Dyin VAR_3Dyang
                    
                else
                    %% READ VECTOR FIELD ---------------------------------------------------
                    % Read pressure & velocity information
                    if strcmp(GRID.Type,'yinyang')
                        if strcmp(FIELD.name,'Toroidal') || strcmp(FIELD.name,'Poloidal')
                            [~,~,~,~,VX_3Dyin,VY_3Dyin,VZ_3Dyin,P_3Dyin,VX_3Dyang,VY_3Dyang,VZ_3Dyang,P_3Dyang,~,~] ...
                                = f_readStagYY(FILE.directory,FILE.name,FILE.number,FIELD.name,SWITCH);
                        else
                            [~,~,~,~,VX_3Dyin,VY_3Dyin,VZ_3Dyin,P_3Dyin,VX_3Dyang,VY_3Dyang,VZ_3Dyang,P_3Dyang,~,~] ...
                                = f_readStagYY(FILE.directory,FILE.name,FILE.number,FIELD.name,SWITCH);
                        end
                        if ischar(VX_3Dyin); [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE); warning off backtrace; warning(VAR_3Dyin); warning on backtrace; continue; end
                        if strcmp(FIELD.name,'Horizontal velocity') || strcmp(FIELD.name,'Horizontal princ. stress');	VAR_3Dyin = VX_3Dyin; VAR_3Dyang = VX_3Dyang;
                        elseif strcmp(FIELD.name,'y-Velocity'); VAR_3Dyin = VY_3Dyin; VAR_3Dyang = VY_3Dyang;
                        elseif strcmp(FIELD.name,'Radial velocity') || strcmp(FIELD.name,'Radial princ. stress');      VAR_3Dyin = VZ_3Dyin; VAR_3Dyang = VZ_3Dyang;
                        elseif strcmp(FIELD.name,'Pressure');   VAR_3Dyin = P_3Dyin;  VAR_3Dyang = P_3Dyang;
                        elseif strcmp(FIELD.name,'Velocity') || ...
                                strcmp(FIELD.name,'Toroidal') || strcmp(FIELD.name,'Poloidal') %combine components
                            VAR_3Dyin  = sqrt(VX_3Dyin.^2 +VY_3Dyin.^2 +VZ_3Dyin.^2);
                            VAR_3Dyang = sqrt(VX_3Dyang.^2 +VY_3Dyang.^2 +VZ_3Dyang.^2);
                        elseif strcmp(FIELD.name,'Streamfunction') %nothing to do
                        end
                        VAR.var3d_yin   = VAR_3Dyin;
                        VAR.var3d_yang  = VAR_3Dyang;
                        
                    else
                        if (strcmp(FIELD.name,'Velocity') || strcmp(FIELD.name,'Horizontal velocity') || strcmp(FIELD.name,'Radial velocity') || strcmp(FIELD.name,'Pressure')) && ...
                                (isfield(PLOT,'VX_3D') && isfield(PLOT,'VY_3D') && isfield(PLOT,'VZ_3D'))
                            VX_3D = PLOT.VX_3D; VY_3D = PLOT.VY_3D; VZ_3D = PLOT.VZ_3D; P_3D = PLOT.P_3D;  %ORIGINAL
                            
                            



%                             warning('This needs check if it works with binary AND hdf5 output')
%                                 if strcmp(GRID.Type,'spherical2D')
%                                     VX_3D = PLOT.VX_3D; VY_3D = PLOT.VY_3D; VZ_3D = -PLOT.VZ_3D; P_3D = PLOT.P_3D;
%                                 elseif  strcmp(GRID.Type,'Cartesian')
%                                     VX_3D = PLOT.VX_3D; VY_3D = PLOT.VY_3D; VZ_3D = -PLOT.VZ_3D; P_3D = PLOT.P_3D; %adjustment for flip of depth vector
%                                 elseif  strcmp(GRID.Type,'yinyang')
%                                     %...
%                                 end
                            
                                
                                
                                
                                
                        else
                            [~,~,~,~,VX_3D,VY_3D,VZ_3D,P_3D,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,FIELD.name,SWITCH);
                            if ischar(VX_3D); [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE); warning off backtrace; warning(VX_3D); warning on backtrace; continue; end
                            if size(P_3D,1)==1 %exchange x and y
                                dummy_3D = zeros(size(VX_3D,2),size(VX_3D,1),size(VX_3D,3));
                                dummy_3Dx = zeros(size(VX_3D,2),size(VX_3D,1),size(VX_3D,3));
                                dummy_3D(:,1,:)     = VZ_3D(1,:,:); VZ_3D = dummy_3D;
                                dummy_3Dx(:,1,:)    = VX_3D(1,:,:);
                                dummy_3D(:,1,:)     = VY_3D(1,:,:); VY_3D = dummy_3Dx; VX_3D = dummy_3D;
                                dummy_3D(:,1,:)     = P_3D(1,:,:);  P_3D = dummy_3D;
                            end
                        end
                        if strcmp(GRID.Dim,'3-D')
                            vx2d        = VX_3D(:,:,:);
                            vy2d        = VY_3D(:,:,:);
                            vz2d        = VZ_3D(:,:,:);
                            p2d         = P_3D(:,:,:);
                        elseif strcmp(GRID.Dim,'2-D')
                            vx2d(:,:)   = VX_3D(:,1,:);
                            vy2d(:,:)   = VY_3D(:,1,:);
                            vz2d(:,:)   = VZ_3D(:,1,:);
                            p2d(:,:)    = P_3D(:,1,:);
                        end
                        v2d = sqrt(vx2d.^2 + vy2d.^2 + vz2d.^2); %absolute velocity
                        if strcmp(FIELD.name,'Velocity');                   VAR.var2d = v2d; end    %PLOT.Velocity
                        if strcmp(FIELD.name,'Horizontal velocity') || strcmp(FIELD.name,'Horizontal princ. stress');  VAR.var2d = vx2d; end   %PLOT.VelocityX, PLOT.StressX
                        if strcmp(FIELD.name,'Radial velocity') || strcmp(FIELD.name,'Radial princ. stress');          VAR.var2d = vz2d; end   %PLOT.VelocityZ, PLOT.StressZ
                        if strcmp(FIELD.name,'Horizontal princ. stress');	VAR.var2d = abs(vx2d); end   %MIGHT REMOVE THAT IF STAGYY OUTPUT SEEMS MORE CORRECT
                        if strcmp(FIELD.name,'Radial princ. stress');     	VAR.var2d = abs(vz2d); end   %MIGHT REMOVE THAT IF STAGYY OUTPUT SEEMS MORE CORRECT
                        if strcmp(FIELD.name,'Pressure');                	VAR.var2d = p2d; end    %PLOT.Pressure
                        if strcmp(FIELD.name,'Streamfunction')
                            vx2d = vx2d*SETUP.Vscale; vz2d = vz2d*SETUP.Vscale; %[cm/a]
                            [phi2d,psi2d] = flowfun(vz2d,vx2d);  % calculate the velocity potential (phi) and streamfunction (psi)
                            if strcmp(PLOT.streamdata,'phi'); stream2d = phi2d; elseif strcmp(PLOT.streamdata,'psi'); stream2d = psi2d; end
                            VAR.var2d   = stream2d;
                        end
                        if exist('VAR','var') && isfield(VAR,'var2d') && ...
                                SWITCH.plotDifference && mod(loopCase,2) %only for odd cases
                            if loopCase>2; error('more than 1 field not yet implemented!'); end
                            v2d_A     = VAR.var2d;
                            continue
                        elseif SWITCH.plotDifference && ~mod(loopCase,2) %only for even cases
                            VAR.var2d   = v2d_A-VAR.var2d;  %take difference
                            clearvars v2d_A
                        end
                    end
                    [PLOT,SAVE] = f_plotField( VAR,FIELD,FIELD_M,GRID,SETUP,SWITCH,PLOT,TOPO,FILE,STYLE,SAVE,PLATE,MANTLE );
                    clearvars VAR v2d vx2d vy2d vz2d p2d stream2d psi2d phi2d
                end  %velocity field
            end %field plot
        end %field loop
        if isfield(SWITCH,'dimensional_output') && SWITCH.DimensionalModeensional_output; SWITCH.DimensionalMode=true; end %set back if dimensional output
        disp(' ')
    end %loop case
    
    if ~SWITCH.PlotInPlot || SAVE.LastFile %if PlotInPlot then only for last file
        %% DESIGN DEFAULTS
        STYLE.FontSize              = STYLE.keyFontSize;
        STYLE.legendFontSize      	= STYLE.keyFontSize-2;
        STYLE.titleFontSize         = STYLE.FontSize+6;
        STYLE.titleFontWeight      	= 'normal';
        
        %% SWITCH ADJUSTEMENTS
        if strcmp(GRID.Type,'yinyang')
            SWITCH.BackgroundDesign     = false; %doesn't work yet
            %         if PLOT.yyNumSlices>1 %for multiple slices
            %             STYLE.cbFactorX = 1.0;
            %         else
            %             STYLE.cbShiftX  = 0.05;
            STYLE.cbFactorY = 0.7;
            %         end
        end
        
        %% SIMPLIFY PLOTS
        if SWITCH.SimplifyPlots
            if SWITCH.AnalysisMode %only remove redundant colorbars
                SWITCH.SimplifyTask     = 'remove unused colorbars';
            end
            [SWITCH,PLOT] = f_DesignSimplify(SWITCH,PLOT,GRID,STYLE);
        end
        
        %% STYLE PLOTS
        if SWITCH.PlotDesign
            STYLE.App                   = SAVE.app; %needed in f_DesignLegend
            STYLE.GraphPlotHandles      = SAVE.GraphPlotHandles;
            STYLE.onlyXaxis             = SWITCH.onlyXaxis;
            if isfield(PLOT,'titleColor'); STYLE.titleColor = PLOT.titleColor; end
            if isfield(PLOT,'ticksColorY'); STYLE.ticksColorY = PLOT.ticksColorY; end
            if SWITCH.LegendDesign; STYLE.Legend = true; else; STYLE.Legend = false; end
            STYLE.axisLineWidth      	= 0.25;  %1.25;
            if ~isfield(STYLE,'titleLocation'); STYLE.titleLocation = 'centered'; if strcmp(SWITCH.timeDirection,'topbot'); STYLE.titleLocation = 'centered'; end; end
            
            %pause(0.1) %let the figure time to adjust
            drawnow
            f_DesignFigure(STYLE,SWITCH);
        end
        
        %% STYLE BACKGROUND
        if SWITCH.BackgroundDesign && ~(strcmp(GRID.Dim,'3-D') && ~SWITCH.plot3Dtopo2D)
            BACK.ColorMode              = STYLE.ColorMode;
            for isp=1:size(SAVE.GraphPlotHandles,1) %only for graph plots
                axes(SAVE.GraphPlotHandles(isp,1)); %activate current graph axis
                BACK.BulletString       = STYLE.SCHAR.hugeBulletDark;
                BACK.ApplyTo            = 'current';
                if isp==1; BACK.display = true; else; BACK.display = false; end
                BACK.BackgroundMode  	= 'simple';
                f_DesignBackground(BACK)
            end
        elseif strcmp(STYLE.ColorMode,'dark') && strcmp(GRID.Dim,'3-D') %&& SWITCH.plot3Dtopo2D
            BACK.ColorMode          = STYLE.ColorMode;
            BACK.BulletString       = STYLE.SCHAR.hugeBulletDark;
            BACK.ApplyTo            = 'all';
            BACK.BackgroundMode  	= 'simple';
            f_DesignBackground(BACK)
        end
        
        %% SET AXES BACK ON TOP (if neccessary)
        if SWITCH.SimplifyPlots || SWITCH.BackgroundDesign %in case they were altered
            for i=1:size(SAVE.PlotHandles,1) %loop all subplot axes (including zoom plots)
                %pause(0.1)
                drawnow
                axes(SAVE.PlotHandles(i,1)); %set back on top
            end
        end
        
        %% INSERT BACKGROUND RECTANGLES
        if SWITCH.BackgroundGuides && nrFiles>1 && ...
                ~strcmp(GRID.Dim,'3-D') && ~strcmp(GRID.Type,'spherical2D')
            STYLE.rectangleColor    = STYLE.keyColor;
            f_DesignPlotRelation(PLOT,SWITCH,STYLE,'insertBackgroundRectangles')
        end
        
        %% INSERT TIME ARROWS
        if SWITCH.TimeEvolutionMode
            STYLE.arrowColor        = STYLE.keyColor;
            f_DesignPlotRelation(PLOT,SWITCH,STYLE,'insertTimeArrows')
        end
        
        %% ANNOTATIONS
        if SWITCH.Annotation
            f_DesignAnnotations(PLOT,STYLE);
        end
        
        %% fAIO ACTION: FINISHING UP
        fAIo.task = 'finishingUp';
        [fAIo,SAVE] = f_AIo(fAIo,SWITCH,FILE,GRID,STYLE,SAVE);
        
        %% SAVING DATA FILE
        if strcmp(STYLE.ColorMode,'light')
            SAVE.blackBackground        = 0;
        elseif strcmp(STYLE.ColorMode,'dark')
            SAVE.blackBackground        = 1;
        end
        if SAVE.Figure || (SAVE.Movie && SAVE.LastFile)
            if SWITCH.waitbar
                disp('   ...saving')
                PLOT.wb = waitbar(0.9,PLOT.wb,'saving...');
            end
            stringX                         = '';
            if strcmp(STYLE.ColorMode,'dark')
                stringX                     = 'Dark';
            end
            
            % SAVING NAME
            if SWITCH.plotDifference; string0 = 'diff'; else; string0 = ''; end
            if nrFiles==1
                stringA     = [IN.Name{1,1},'(',num2str(IN.Number(1,1)),')'];
                stringA2    = IN.Name{1,1};
            else
                stringA     = [IN.Name{1,1},'(',num2str(IN.Number(1,1)),')_vs_',IN.Name{1,2},'(',num2str(IN.Number(1,2)),')'];
                stringA2    = [IN.Name{1,1},'_vs_',IN.Name{1,2}];
            end
            if nrFiles>2
                stringB     = ['_vs_',IN.Name{1,3},'(',num2str(IN.Number(1,3)),')'];
                stringB2    = ['_vs_',IN.Name{1,3}];
            else
                stringB     = ''; stringB2 = '';
            end
            stringC         = fieldString;
            stringD         = SWITCH.ColormapAll{1,1};
            if strcmp(SWITCH.ColormapAll{1,2},'flip'); stringE = ['_',SWITCH.ColormapAll{1,2}]; else; stringE = ''; end
            
            SAVE.FigureName = ['+',string0,stringA,stringB,stringC,'_',stringD,stringE,stringX];
            SAVE.MovieName  = ['++',string0,stringA2,stringB2,stringC,stringE,stringX];
            
            % SAVING DIRECTORY
            if strcmp(SAVE.writeDirectory,'auto') %save to standard folder
                SAVE.Directory = FILE.stemSave;
            else
                SAVE.Directory = SAVE.writeDirectory;
            end
        end
        % SAVING FIGURE
        if SAVE.Figure
            SAVE.FigureNr               = 1;
            SAVE.keepAspect             = true;
            f_saveFigure( SAVE );
        end
        % SAVING MOVIE
        if SAVE.Movie
            lastFrame   = SAVE.LastFile;
            [SAVE] = f_saveMovie(SAVE,lastFrame);
        end
        if SAVE.Figure || (SAVE.Movie && SAVE.LastFile)
            if SWITCH.waitbar; waitbar(1,PLOT.wb); end
        end
    end
    
    %% KEEPING DATA
    if isfield(PLOT,'topCval'); SAVE.topCval = PLOT.topCval; end %keeping composition tracker data

catch errorMessage
    if SWITCH.sendErrorLog
        disp(' ')
        disp('ERROR MESSAGE:')
        disp(' ')
        disp(errorMessage.getReport('extended', 'hyperlinks','off'))
        %delete error log file if sent successfully
        diary OFF
        
        %% fAIO ACTION: SENDING ERROR LOG
        fAIo.task   = 'sendErrorLog';
        [fAIo,SAVE] = f_AIo(fAIo,SWITCH,[],[],STYLE,SAVE);
    else
        %Throw an error
        disp(' ');
        fprintf(2,errorMessage.getReport('extended')); pause(0.1);
        disp(' '); beep; SAVE.LastFile = true;
        return;
    end
end

%% FINISHING UP
if SWITCH.waitbar
    PLOT.wb = waitbar(1,PLOT.wb,'--- FINISHED ---');
    pause(0.5)
    close(PLOT.wb)
end

%% fAIO ACTION: END MESSAGE
if SAVE.LastFile
    fAIo.task = 'EndMessage';
    [fAIo,SAVE] = f_AIo(fAIo,SWITCH,FILE,GRID,STYLE,SAVE);
else
    disp(' ')
end




%% INTERNAL FUNCTIONS


%% RANGE
function y = range(x,dim)
    % Replacement function for MatLab's 'range' (which necessitates some extra libraries)
    % USAGE:
    % range(X)
    % y = range(X,dim)
    %
    % Description (from Matlab documentation):
    %   Y = RANGE(X) returns the range of the values in X.  For a vector input,
    %   Y is the difference between the maximum and minimum values.  For a
    %   matrix input, Y is a vector containing the range for each column.  For
    %   N-D arrays, RANGE operates along the first non-singleton dimension.
    %
    %   RANGE treats NaNs as missing values, and ignores them.
    %   Y = RANGE(X,DIM) operates along the dimension DIM.
    %
    %                                         Marcel Thielmann, 14.11.2017
    if nargin < 2
        y = max(x) - min(x);
    else
        y = max(x,[],dim) - min(x,[],dim);
    end
end



end
