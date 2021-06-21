
%%                                            DEFAULT INPUT VARIABLES 2.23
%                                                     for SL_RadialProfile
% 
%                                                Fabio Crameri, 17.06.2021

function [IN,PLOT,SWITCH,STYLE,SAVE] = f_DefaultsRprof
if ~exist('SL_RadialProfile','file'); f_INSTALL; end

%% COLLECTING STAGLAB PERFORMANCE AND ERROR DATA **********************
% StagLab collects user-specific performance and error data
% that is stored in /StagLab/fAIo/cache.mat
% The data is used to keep StagLab up-to-date, make parfiles forward-compatible,
% and provide information in case of StagLab user-specific bugs and errors.
SWITCH.AnonymiseStagLabData     = 	logical(0);     %Anonymise data
% *********************************************************************
    
%% DEBUGGING
SWITCH.Verbose                  =   logical(0);                 %Switch for extended display output
SWITCH.sendErrorLog             =   logical(0);              	%Sends an error log file to the developer for debugging purposes
SWITCH.readOption2           	=   logical(0);              	%Alternative option to read the rprof.dat file, which might circumvent some problems

%% INPUT FILE(S)
SWITCH.GeodynamicCode           =   'StagYY';                 	%To use StagLab with a different code, run and adjust at locations of error messages
IN.Name                         =   {'test'  'test'  };         %Filename  (it's size controls the # of files processed)
IN.Number                       =   [ 0    1 ];                 %Filenumber (i.e., timestep)
IN.Parameter                    =   [ 11   ];                   %Parameter dimensionalisation (see f_Dimensions)
IN.Folder                       =   {'/work/stagyy/'};          %File directory
%PLOT.titleString                =   {'free-slip surface' 'free surface'}; %Manual title if uncommented

%% StagLab MODES
SWITCH.DimensionalMode       	=   logical(0);                 %Converts to dimensional output as set in f_Dimensions.m
    SWITCH.modelSetup           =   'standard';                 %See f_Dimensions.m
SWITCH.AnalysisMode             =   logical(0);                 %Creates figures for analysis purposes - not recommended for publication!
SWITCH.QuickMode                =   logical(0);                 %Omits plotting profiles for faster post-processing of files
STYLE.ColorMode                 =   'light';                    %'light' or 'dark' figure background

%% FIGURE SETUP
SWITCH.FigurePosition           =   'auto';                     %'auto' or set manually by e.g., [47 1 1300 350]
PLOT.FigureDefaultPosition      =   'TopRight';                 %Figure position on screen: 'BottomLeft', 'BottomRight', 'TopLeft', or 'TopRight'

%% STYLING
SWITCH.PlotDesign               =   logical(1);
STYLE.Mode                      =   'Custom';                   %'Custom': adjustable plot design | 'Nature','Science','EGU','AGU','Frontiers': preset plot design
SWITCH.BackgroundDesign     	=   logical(1);
SWITCH.LegendDesign          	=   logical(1);
SWITCH.SimplifyPlots          	=   logical(1);
SWITCH.UsePanel                 =   logical(1);                 %Uses novel subplot distribution (still under development)
STYLE.Brackets4Dimensions       =   'none';                     %Puts dimensions into 'Parentheses','Brackets', or 'none'
%axes
STYLE.MinorAxisTicks            =   logical(0);
%font
STYLE.AllFontName               =   'Archivo';                %e.g. 'Helvetica','Corbel','Eurostyle','ETH Light','Archivo'
STYLE.keyFontSize               =   10;
%colours
STYLE.keyColor                  =   [0.5 0.5 0.5];              %Key colour used for text and axes
%lines
IN.LineWidthPlot                =   1.0;

%% AXES DETAILS
IN.constantAxis                 =   logical(0);
IN.reverseYaxis                 =   logical(1);
SWITCH.AxesLimit                =   logical(0);                 %Switch on manual y-axis limits
    IN.YAxisMinValue            =   NaN;                        %Minimum y-axis limit in plot dimension (set automatically when NaN)
    IN.YAxisMaxValue            =   600;                       	%Maximum y-axis limit in plot dimension (set automatically when NaN)

%% LABELLING
SWITCH.Title                    =   logical(0);
SWITCH.Annotation               =   logical(1);                 %Subplot annotations with e.g., 'a', 'b', 'c'
    STYLE.annotationLocation    =   'botLeft';                  %e.g., 'topOutTop','topOutCorner','topLeft','botLeft','topRight','botRight'
    PLOT.annotationString      	=   {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p'};  %subplot annotation strings
    PLOT.annIfMoreThanOne       =   logical(1);                 %Adds annotation only if more than 1 subplot
SWITCH.Texting                  =   logical(0);                 %Add text of min, mean and max values of each field

%% PLOT ADDITIONS
SWITCH.d_eta                    =   logical(0);                 %Displays viscosity increase with depth
PLOT.Solidus                    =   logical(0);                 %Adds a Solidus curve to the temperature profile
PLOT.parameterRange             =   logical(1);                 %Plotting mean, min & max
    PLOT.rangeAsArea            =   logical(1);
IN.showLegend                   =   logical(0);
    PLOT.LegendLocation         =   'SouthEast';
SWITCH.LegendOutside            =   logical(0);                 %Creates an extra space to put the legend outside of the plot

%% PARAMETER(S) TO PLOT
PLOT.Temperature                =   logical(1);         IN.T_logx       = 0; 
PLOT.Viscosity                  =   logical(1);         IN.eta_logx     = 1; 
PLOT.Density                    =   logical(0);         IN.rho_logx     = 0;
PLOT.Stress                     =   logical(1);         IN.str_logx     = 0;
PLOT.StrainRate                 =   logical(0);         IN.edot_logx    = 1; 
PLOT.Velocity                   =   logical(0);         IN.v_logx       = 0;
PLOT.VerticalVelocity           =   logical(0);     	IN.vz_logx      = 0;
PLOT.HorizontalVelocity         =   logical(1);         IN.vh_logx      = 0; 

PLOT.HorizVorticity             =   logical(0);         IN.w_logx       = 0;
PLOT.VertVorticity              =   logical(0);         IN.w_logx       = 0;
PLOT.Divergence                 =   logical(0);         IN.div_logx 	= 0;

PLOT.Advection                  =   logical(0);         IN.Hadv_logx    = 0;
PLOT.Diffusion                  =   logical(0);         IN.Hdiff_logx   = 0;
PLOT.InternalHeating            =   logical(0);         IN.Rh_logx      = 0; 
PLOT.ViscousDissipation       	=   logical(0);         IN.vd_logx      = 0; 
PLOT.AdiabaticHeating           =   logical(0);         IN.ah_logx      = 0;
PLOT.Heatflux                   =   logical(0);         IN.HF_logx      = 0; 

PLOT.Crust                      =   logical(0);         IN.c_logx       = 0;
PLOT.Air                        =   logical(0);         %log x-axis is switched on with IN.c_logx
PLOT.ContinentalCrust         	=   logical(0);         % "
PLOT.Primordial                 =   logical(0);         % "
PLOT.Fluid                      =   logical(0);         % "
PLOT.Metal                      =   logical(0);         % "

PLOT.Water                      =   logical(0);         % "

IN.plates_analyse               =   logical(0);         IN.z_level  = 0.01; %0.01 (surface); ~0.15 (slabs)

PLOT.Toroidal                   =   logical(0);         IN.To_logx      = 0; 
PLOT.ToroidalRMS                =   logical(0);         IN.To_logx      = 0;
PLOT.Poloidal                   =   logical(0);         IN.Po_logx      = 0; 
PLOT.PoloidalRMS                =   logical(0);         IN.Po_logx      = 0;
PLOT.Mobility                   =   logical(0);         IN.M_logx       = 0;
PLOT.f80                        =   logical(0);         IN.f80_logx     = 0;    %Plateness
PLOT.f90                        =   logical(0);         IN.f90_logx     = 0;    %Plateness

PLOT.RefstatCaTemperature     	=   logical(0);   
PLOT.RefstatCaDensity         	=   logical(0);         
PLOT.RefstatCaExpansivity     	=   logical(0);         
PLOT.RefstatCaConductivity     	=   logical(0);         
PLOT.RefstatCaPressure         	=   logical(0);  

%These still need to be fully implemented:
PLOT.ViscousDissipationLog     	=   logical(0);  
PLOT.AdvectionTotal         	=   logical(0);  
PLOT.ThermalConduction         	=   logical(0);  

%% COLOUR SET FOR MULTIPLE LINES
PLOT.ColorVector                =   [0.0 0.0 0.0; 0.0 0.45 0.74]; %Default colour for 1 or 2 lines
PLOT.ScientificColourMapName    =   'batlow';                	%'auto' or name of a Scientific colour map to be used for 3 or more lines 
                                                             	%e.g., 'batlow' for ordered and 'batlowS' for unordered colours - for all options see >StagLab >Colourmaps >+ScientificColourMaps_FabioCrameri.png
%% DEFINE Y-PARAMETER
PLOT.ParameterY                 =   1;

%% SAVING FIGURE
SAVE.Figure                  	=   logical(0);                 %Saves plot to directory 
    SAVE.png                    =   logical(1);     
        SAVE.pngResolution      =   '-m4.0';
        SAVE.TransparentBackground = logical(0);                %Removes the figure background (only available for .png format)
    SAVE.jpg                    =   logical(0);
        SAVE.jpgResolution    	=   '-r400';
    SAVE.eps                    =   logical(0);
    SAVE.pdf                    =   logical(0);
    IN.saveString               =   '+Profile';
  	SAVE.writeDirectory         =   'auto';                     %'auto': save to standard folder; otherwise use e.g., '/work/stagyy/'

%% SAVING DATA
SAVE.Profiles                	=   logical(0);                 %Saves processed plotted profiles to data files
    SAVE.mat                    =   logical(0);
    SAVE.dat                    =   logical(1);
    SAVE.txt                    =   logical(0);
    
end









