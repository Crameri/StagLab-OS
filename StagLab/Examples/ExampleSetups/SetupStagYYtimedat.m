
%%                                          TEST PARAMETER FILE TIME GRAPH
%                                                         for SL_TimeGraph

function [TEST] = SetupStagYYtimedat(TEST)

if ~exist('f_defaultsTimedat','file'); f_INSTALL; end
[IN,PLOT,SWITCH,STYLE,SAVE] = f_DefaultsTimedat; %get default variables

%% INPUT FILE(S)
IN.Name                         =   { 'StagYYdim'	};          %filename
IN.Parameter                    =   [   11         ];
IN.Folder                       =   {'../'};                    %file directory

%% SAVING
SAVE.Figure                   	=   logical(1);                 % saves plot directly to directory   
  	SAVE.writeDirectory         =   '../ExampleFigures/';       %'auto': save to standard folder (.../+im/...); otherwise use e.g., '/work/stagyy/'

%% DIMENSIONALISATION
SWITCH.DimensionalMode      	=   logical(1);

%% PLOT ADDITIONS & LABELLING
SWITCH.parameterRange           =   logical(0);                 %plotting mean, min & max
   	PLOT.rangeAsArea            =   logical(1);                 %min/max range as area
    SWITCH.MinMaxLegend         =   logical(0);                 %legend switch for min and max lines

%% PARAMETER(S) TO PLOT
PLOT.Temperature                =   logical(1);
PLOT.Viscosity                  =   logical(1);
PLOT.Velocity                   =   logical(1);
PLOT.Composition                =   logical(0);

PLOT.Ra_eff                     =   logical(1);   %only 1 value
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

 %% RUN MAIN ROUTINE >>>
[SAVE] = SL_TimeGraph(IN,PLOT,SWITCH,STYLE,SAVE);

TEST.successful = true;
end

