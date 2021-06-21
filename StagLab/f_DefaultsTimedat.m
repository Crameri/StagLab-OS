
%%                                             DEFAULT INPUT VARIABLES 2.41
%                                                          for SL_TimeGraph
% 
%                                                Fabio Crameri, 17.06.2021

function [IN,PLOT,SWITCH,STYLE,SAVE] = f_DefaultsTimedat
if ~exist('SL_TimeGraph','file'); f_INSTALL; end

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

%% INPUT FILE(S)
SWITCH.GeodynamicCode           =   'StagYY';                	%To use StagLab with a different code, run and adjust at locations of error messages
IN.Name                         =   { 'test'  'test' };         %Filename (it's size controls the # of files processed)
IN.Parameter                    =   [   11         ];           %Parameter dimensionalisation (see f_Dimensions)
IN.Folder                       =   {'/work/stagyy140211/'};    %File directory
% PLOT.titleString               =   {'case1' 'case2'};      	%Manual legend description if uncommented

%% StagLab MODES
SWITCH.DimensionalMode      	=   logical(0);                 %Converts to dimensional output as set in f_Dimensions.m
    SWITCH.modelSetup           =   'standard';             	%See f_Dimensions, e.g., 'standard', 'solomatov2004'
SWITCH.AnalysisMode             =   logical(0);                 %Creates figures for analysis purposes - not recommended for publication!
SWITCH.QuickMode                =   logical(0);                 %Omits plotting graphs for faster post-processing of files
STYLE.ColorMode                 =   'light';                    %'light' or 'dark' figure background

%% FIGURE SETUP
SWITCH.FigurePosition           =   'auto';                     %'auto' or set manually by e.g., [47 1 1300 350]
PLOT.FigureDefaultPosition      =   'TopRight';                 %Figure position on screen: 'BottomLeft', 'BottomRight', 'TopLeft', or 'TopRight'
    
%% STYLING
SWITCH.PlotDesign               =   logical(1);
STYLE.Mode                      =   'Custom';                   %'Custom': adjustable plot design | 'Nature','Science','EGU','AGU','Frontiers': preset plot design
SWITCH.BackgroundDesign      	=   logical(1);
SWITCH.LegendDesign         	=   logical(1);
SWITCH.SimplifyPlots         	=   logical(1);
SWITCH.UsePanel                 =   logical(1);                 %Uses novel subplot distribution (still under development)
STYLE.Brackets4Dimensions       =   'none';                     %Puts dimensions into 'Parentheses','Brackets', or 'none'
%axes
STYLE.MinorAxisTicks            =   logical(0);
%font
STYLE.AllFontName               =   'Archivo';                %e.g. 'Helvetica','Corbel','Eurostyle','ETH Light','Archivo'
STYLE.keyFontSize               =   10;
%colours
STYLE.keyColor                  =   [0.5 0.5 0.5];              %Key colour used for text and axes

%% APP PERFORMANCE
IN.numWrite                     =   1;                          %CURRENTLY REMOVED! writes out every 'num_write' time step
SWITCH.enforceFixingTimedatFile =   logical(0);                 %Enforces attempt to automatically fix corrupt time.dat files

%% PLOT ADDITIONS & LABELLING
SWITCH.parameterRange           =   logical(0);                 %Plotting mean, min & max
   	PLOT.rangeAsArea            =   logical(1);                 %min/max range as area
    SWITCH.MinMaxLegend         =   logical(0);              	%Legend switch for min and max lines
SWITCH.AxesLimit                =   logical(0);                 %Switch on manual x-axis limits
    IN.XAxisMinValue         	=   NaN;                        %Minimum x-axis limit in plot dimension (set automatically when NaN)
    IN.XAxisMaxValue          	=   5;                          %Maximum x-axis limit in plot dimension (set automatically when NaN)
    IN.YAxisMinValue         	=   [NaN,NaN,NaN];            	%Minimum y-axis limit in plot dimension for each subplot (set automatically when NaN)
    IN.YAxisMaxValue          	=   [NaN,NaN,NaN];            	%Maximum y-axis limit in plot dimension for each subplot (set automatically when NaN)

    %% LABELLING
SWITCH.Annotation               =   logical(1);                 %Subplot annotations with e.g., 'a', 'b', 'c'
    STYLE.annotationLocation    =   'botLeft';                  %'topOutTop','topOutCorner','topLeft','botLeft','topRight','botRight'
    PLOT.annotationString     	=   {'a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p'};  %subplot annotation strings
    PLOT.annIfMoreThanOne       =   logical(1);                 %Add annotation only if more than 1 subplot
SWITCH.Texting                  =   logical(0);                 %Add text of min, mean and max values of each field
SWITCH.LegendOutside            =   logical(0);                 %Creates an extra space to put the legend outside of the plot

%% PARAMETER(S) TO PLOT
PLOT.Temperature                =   logical(1);
PLOT.Viscosity                  =   logical(1);     IN.eta_logy     = 0;
PLOT.Velocity                   =   logical(1);     IN.v_logy       = 0;
PLOT.Composition                =   logical(0);

PLOT.Ra_eff                     =   logical(1);     IN.Ra_logy    	= 0;	%only 1 value
PLOT.Nu_top                     =   logical(0);   %only 1 value
PLOT.Nu_bot                     =   logical(0);   %only 1 value
PLOT.InternalHeating            =   logical(1);   %only 1 value
PLOT.HeatfluxTop                =   logical(0);   %only 1 value
PLOT.HeatfluxBot                =   logical(0);   %only 1 value
PLOT.EruptionRate               =   logical(0);   %only 1 value
PLOT.EruptHeatFlux              =   logical(0);   %only 1 value
PLOT.Erupta                     =   logical(0);   %only 1 value
PLOT.Entrainment                =   logical(0);   %only 1 value
PLOT.Cmass_error                =   logical(0);   %only 1 value
PLOT.Timesteps                  =   logical(0);   %only 1 value
PLOT.Performance                =   logical(1);   %only 1 value
 
%New options: need to be tested!
PLOT.InnerCoreRadius            =   logical(0);   %only 1 value
PLOT.TotalErupta                =   logical(0);   %only 1 value
PLOT.OutgassedWater             =   logical(0);   IN.waterOut_logy 	= 0;    %only 1 value
PLOT.SurfaceTemperature         =   logical(0);   %only 1 value
PLOT.CMB_Temperature            =   logical(0);   %only 1 value
PLOT.Psurf                      =   logical(0);   %only 1 value
PLOT.Weathering                 =   logical(0);   %only 1 value
PLOT.TotalWater                 =   logical(0);   %only 1 value
PLOT.MantleWater                =   logical(0);   %only 1 value
PLOT.CrustWater                 =   logical(0);   %only 1 value
PLOT.dmH2O                      =   logical(0);   %only 1 value
PLOT.OutgassedCarbon            =   logical(0);   %only 1 value
PLOT.OutgassedNitrogen          =   logical(0);   %only 1 value
PLOT.Psurf                      =   logical(0);   %only 1 value
PLOT.s_core                     =   logical(0);   %only 1 value
PLOT.tc_core                    =   logical(0);   %only 1 value
PLOT.ts_core                    =   logical(0);   %only 1 value
PLOT.OutgassedAmountWater     	=   logical(0);   IN.waterAmOut_logy = 0;    %only 1 value
PLOT.OutgassedWaterExtra        =   logical(0);   IN.waterExOut_logy = 0;    %only 1 value
PLOT.H2Odegassing             	=   logical(0);   IN.waterDeg_logy  = 0;     %only 1 value
PLOT.H2Oregassing              	=   logical(0);   IN.waterReg_logy  = 0;     %only 1 value

%% COLOUR SET FOR MULTIPLE LINES
PLOT.ColorVector                =   [0.0 0.0 0.0; 0.0 0.45 0.74]; %Default colour for 1 or 2 lines
PLOT.ScientificColourMapName    =   'batlow';                	%'auto' or name of a Scientific colour map to be used for 3 or more lines 
                                                             	%e.g., 'batlow' for ordered and 'batlowS' for unordered colours - for all options see >StagLab >Colourmaps >+ScientificColourMaps_FabioCrameri.png
%% SAVING
SAVE.Figure                   	=   logical(0);                 %Saves plot to directory   
    SAVE.png                    =   logical(1);
        SAVE.pngResolution      =   '-m4.0';
        SAVE.TransparentBackground = logical(0);                %Removes the figure background (only available for .png format)
    SAVE.jpg                    =   logical(0); 
        SAVE.jpgResolution    	=   '-r400';
    SAVE.eps                    =   logical(0); 
    SAVE.pdf                    =   logical(0); 
    IN.saveString               =   '+Evolution';
  	SAVE.writeDirectory         =   'auto';                     %'auto': save to standard folder (.../+im/...); otherwise use e.g., '/work/stagyy/'

%% SAVING DATA
SAVE.Graphs                     =   logical(0);                 %Saves processed plotted graphs to data files
    SAVE.mat                    =   logical(0);
    SAVE.dat                    =   logical(1);
    SAVE.txt                    =   logical(0);

end









