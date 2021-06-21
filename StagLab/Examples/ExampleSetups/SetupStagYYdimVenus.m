
%%                                                 TEST PARAMETER FILE 2-D
%                                                         for SL_FieldPlot

function [TEST] = SetupStagYYdimVenus(TEST)

SAVE.StartNumber       =   880;      %Set SAVE.StartNumber==SAVE.EndNumber to write only one file
SAVE.StepNumber        =   1;
SAVE.EndNumber         =   880;      %Set big in order to continue to last file

for iLoop=SAVE.StartNumber:SAVE.StepNumber:SAVE.EndNumber
    clearvars -except iLoop SAVE; if ~exist('f_Defaults','file'); f_INSTALL; end
    [IN,PLOT,SWITCH,TOPO,STYLE,SAVE] = f_Defaults(SAVE); %get default variables

    %% INPUT FILE(S)
    IN.Name                     =   { 'StagYYdimVenus'};     	%filename; it's size controls the # of subplots
    IN.Number                   =   [   iLoop    ];           	%filenumber
    IN.Parameter                =   [  	1    ];             	%parameter dimensionalisation (see f_dimParameter)
    IN.Folder                   =   {'../'};                    %file directory

    %% SUBPLOT LAYOUT
    PLOT.ReverseLayout          =   true;                       %e.g., plot 3x2 instead of 2x3
    
    %% POST-PROCESSING
    PLOT.pHot                   =   0.5; %percentage of difference Tmax-Tmean and Tmin-Tmean, resp. (T_threshold = Tmean+p_hot*(Tmax-Tmean) )
    PLOT.pCold                  =   0.7; %percentage of difference Tmax-Tmean and Tmin-Tmean, resp. (T_threshold = Tmean+p_hot*(Tmax-Tmean) )

    %% SAVING FIGURE
    SAVE.Figure                 =   logical(1);                 %saves plot directly to directory
  	SAVE.writeDirectory         =   '../ExampleFigures/';       %'auto': save to standard folder (.../+im/...); otherwise use e.g., '/work/stagyy/'

    %% FIELDS TO PLOT
    PLOT.Temperature            =   logical(1); 
    PLOT.Viscosity              =   logical(1);
    PLOT.Topography             =   logical(1);
    
    %% SPECIAL PLOTS
    PLOT.SurfaceFieldVariation  =   logical(1);                 %bar plot of surface field variation
        PLOT.sfvField           =   'Temperature';              %specify field for surface variation bar plot
        PLOT.sfvDepthLevel      =   0;                          %[km], Specifies depth level, e.g., 0:surface 2890:CMB
        PLOT.sfvNumberBins      =   24;
            
    %% RUN MAIN ROUTINE >>>
    [SAVE] = SL_FieldPlot(IN,PLOT,SWITCH,TOPO,STYLE,SAVE);      %run main routine
    if SAVE.LastFile; break; end %last file reached
    
end

TEST.successful = true;
end

