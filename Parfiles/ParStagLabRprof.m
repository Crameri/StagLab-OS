
%%                                           StagLab PARAMETER FILE: RPROF
%                        (for all available options see f_DefaultsRprof.m)

clear;
if ~exist('f_DefaultsRprof','file'); f_INSTALL; end
[IN,PLOT,SWITCH,STYLE,SAVE] = f_DefaultsRprof; %get default variables

%% DEBUGGING
SWITCH.sendErrorLog             =   logical(0);                 %Sends an error log file to the developer for debugging purposes

%% INPUT FILE(S)
IN.Name                         =   { 'test' };                 %Filename  (it's size controls the # of files processed)
IN.Number                       =   [ 1     ];                  %Filenumber (i.e., timestep)
IN.Parameter                    =   [ 11    ];                  %Parameter dimensionalisation (see f_Dimensions)
IN.Folder                       =   {'~/work/stagyy/'};       	%File directory
%PLOT.titleString                =   {'free-slip surface' 'free surface'}; %Manual title if uncommented

%% StagLab MODES
SWITCH.DimensionalMode       	=   logical(0);                 %Converts to dimensional output as set in f_Dimensions.m
SWITCH.AnalysisMode             =   logical(0);                 %Creates figures for analysis purposes - not recommended for publication!
SWITCH.QuickMode                =   logical(0);                 %Omits plotting profiles for faster post-processing of files
STYLE.ColorMode                 =   'light';                    %'light' or 'dark' figure background

%% AXES DETAILS
IN.constantAxis                 =   logical(0);
SWITCH.AxesLimit                =   logical(0);                 %Switch on manual y-axis limits
    IN.YAxisMinValue            =   NaN;                        %Minimum y-axis limit in plot dimension (set automatically when NaN)
    IN.YAxisMaxValue            =   600;                       	%Maximum y-axis limit in plot dimension (set automatically when NaN)

%% LABELLING
SWITCH.Annotation               =   logical(1);                 %Subplot annotations with e.g., 'a', 'b', 'c'
SWITCH.Texting                  =   logical(0);                 %Add text of min, mean and max values of each field

%% PLOT ADDITIONS
PLOT.parameterRange             =   logical(1);                 %Plotting mean, min & max
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

IN.plates_analyse               =   logical(0);         IN.z_level  = 0.01; %0.01 (surface); ~0.15 (slabs)

PLOT.Toroidal                   =   logical(0);         IN.To_logx      = 0; 
PLOT.ToroidalRMS                =   logical(0);         IN.To_logx      = 0;
PLOT.Poloidal                   =   logical(0);         IN.Po_logx      = 0; 
PLOT.PoloidalRMS                =   logical(0);         IN.Po_logx      = 0;
PLOT.Mobility                   =   logical(0);         IN.M_logx       = 0;
PLOT.f80                        =   logical(0);         IN.f80_logx     = 0;
PLOT.f90                        =   logical(0);         IN.f90_logx     = 0;

%% COLOUR SET FOR MULTIPLE LINES
PLOT.ColorVector                =   [0.0 0.0 0.0; 0.0 0.45 0.74]; %Default colour for 1 or 2 lines
PLOT.ScientificColourMapName    =   'batlow';                	%'auto' or name of a Scientific colour map to be used for 3 or more lines 
                                                             	%e.g., 'batlow' for ordered and 'batlowS' for unordered colours - for all options see >StagLab >Colourmaps >+ScientificColourMaps_FabioCrameri.png
%% SAVING FIGURE
SAVE.Figure                  	=   logical(0);                 %Saves plot to directory 
  	SAVE.writeDirectory         =   'auto';                     %'auto': save to standard folder; otherwise use e.g., '/work/stagyy/'

%% SAVING DATA
SAVE.Profiles                	=   logical(0);                 %Saves processed plotted profiles to data files
 
%% RUN MAIN ROUTINE >>>
[SAVE] = SL_RadialProfile(IN,PLOT,SWITCH,STYLE,SAVE);



