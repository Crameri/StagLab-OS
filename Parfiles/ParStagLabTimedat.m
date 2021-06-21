
%%                                         StagLab PARAMETER FILE: TIMEDAT
%                      (for all available options see f_DefaultsTimedat.m)

clear;
if ~exist('f_DefaultsTimedat','file'); f_INSTALL; end
[IN,PLOT,SWITCH,STYLE,SAVE] = f_DefaultsTimedat; %get default variables

%% DEBUGGING
SWITCH.sendErrorLog             =   logical(0);               	%Sends an error log file to the developer for debugging purposes

%% INPUT FILE(S)
IN.Name                         =   { 'test' };                 %Filename (it's size controls the # of files processed)
IN.Parameter                    =   [   11 	 ];                 %Parameter dimensionalisation (see f_Dimensions)
IN.Folder                       =   {'~/work/plotTest/'};        %File directory
% PLOT.titleString               =   {'case1' 'case2'};      	%Manual legend description if uncommented

%% StagLab MODES
SWITCH.DimensionalMode      	=   logical(0);                 %Converts to dimensional output as set in f_Dimensions.m
SWITCH.AnalysisMode             =   logical(0);                 %Creates figures for analysis purposes - not recommended for publication!
STYLE.ColorMode                 =   'dark';                    %'light' or 'dark' figure background 

%% LABELLING
SWITCH.Annotation               =   logical(1);                 %Subplot annotations with e.g., 'a', 'b', 'c'
SWITCH.Texting                  =   logical(0);                 %Add text of min, mean and max values of each field
SWITCH.LegendOutside            =   logical(0);                 %Creates an extra space to put the legend outside of the plot

%% PLOT ADDITIONS
SWITCH.parameterRange           =   logical(0);                 %Plotting mean, min & max
SWITCH.AxesLimit                =   logical(0);                 %Switch on manual x-axis limits
    IN.XAxisMinValue         	=   NaN;                        %Minimum x-axis limit in plot dimension (set automatically when NaN)
    IN.XAxisMaxValue          	=   5;                          %Maximum x-axis limit in plot dimension (set automatically when NaN)

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
PLOT.TotalErupta                =   logical(0);   %only 1 value
PLOT.OutgassedWater             =   logical(0);   %only 1 value
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

%% COLOUR SET FOR MULTIPLE LINES
PLOT.ColorVector                =   [0.0 0.0 0.0; 0.0 0.45 0.74]; %Default colour for 1 or 2 lines
PLOT.ScientificColourMapName    =   'batlow';                	%'auto' or name of a Scientific colour map to be used for 3 or more lines 
                                                             	%e.g., 'batlow' for ordered and 'batlowS' for unordered colours - for all options see >StagLab >Colourmaps >+ScientificColourMaps_FabioCrameri.png
%% SAVING
SAVE.Figure                   	=   logical(0);                 %Saves plot to directory   
    SAVE.png                    =   logical(1);
        SAVE.pngResolution      =   '-m4.0';
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


 %% RUN MAIN ROUTINE >>>
[SAVE] = SL_TimeGraph(IN,PLOT,SWITCH,STYLE,SAVE);




