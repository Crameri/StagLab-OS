
%%                                TEST PARAMETER FILE RADIAL PROFILE GRAPH
%                                                     for SL_RadialProfile

function [TEST] = SetupStagYYrprof(TEST)

if ~exist('SL_RadialProfile','file'); f_INSTALL; end
[IN,PLOT,SWITCH,STYLE,SAVE] = f_DefaultsRprof;                  %get default variables

%% INPUT FILE(S)
IN.Name                         =   {'StagYYnd'  'StagYYnd'  };	%filename
IN.Number                       =   [   0    1 ];
IN.Parameter                    =   [  	11   ];                 %see f_dimParameter.m for more details
IN.Folder                       =   {'../'};

%% DIMENSIONALISATION
SWITCH.DimensionalMode       	=   logical(1);

%% SAVING
SAVE.Figure                  	=   logical(1);                 %saves plot directly to directory 
  	SAVE.writeDirectory         =   '../ExampleFigures/';       %'auto': save to standard folder (.../+im/...); otherwise use e.g., '/work/stagyy/'

%% PARAMETER(S) TO PLOT
PLOT.Temperature                =   logical(1);         IN.T_logx       = 0;  %*
PLOT.Viscosity                  =   logical(1);         IN.eta_logx     = 1;  %*
PLOT.Density                    =   logical(0);         IN.rho_logx     = 0;
PLOT.Stress                     =   logical(1);         IN.str_logx     = 0;  %*
PLOT.StrainRate                 =   logical(0);         IN.edot_logx    = 1;  %*
PLOT.Velocity                   =   logical(0);         IN.v_logx       = 0;
PLOT.VerticalVelocity           =   logical(0);     	IN.vz_logx      = 0;
PLOT.HorizontalVelocity         =   logical(1);         IN.vh_logx      = 0;  %*

PLOT.HorizVorticity             =   logical(0);         IN.w_logx       = 0;
PLOT.VertVorticity              =   logical(0);         IN.w_logx       = 0;
PLOT.Divergence                 =   logical(0);         IN.div_logx 	= 0;

PLOT.Advection                  =   logical(0);         IN.Hadv_logx    = 0;
PLOT.Diffusion                  =   logical(0);         IN.Hdiff_logx   = 0;
PLOT.InternalHeating            =   logical(0);         IN.Rh_logx      = 0;  %*
PLOT.ViscousDissipation      	=   logical(0);         IN.vd_logx      = 0;  %*
PLOT.AdiabaticHeating           =   logical(0);         IN.ah_logx      = 0;  %*

PLOT.Heatflux                   =   logical(0);         IN.HF_logx      = 0;  %*

PLOT.Crust                      =   logical(0);         IN.c_logx       = 0;
PLOT.Air                        =   logical(0);         %log x-axis is switched on with IN.c_logx
PLOT.ContinentalCrust         	=   logical(0);         % "
PLOT.Primordial                 =   logical(0);         % "
PLOT.Fluid                      =   logical(0);         % "
PLOT.Metal                      =   logical(0);         % "

IN.plates_analyse               =   logical(0);         IN.z_level  = 0.01; %0.01 (surface); ~0.15 (slabs)

PLOT.Toroidal                   =   logical(0);         IN.To_logx      = 0;  %*
PLOT.ToroidalRMS                =   logical(0);         IN.To_logx      = 0;
PLOT.Poloidal                   =   logical(0);         IN.Po_logx      = 0;  %*
PLOT.PoloidalRMS                =   logical(0);         IN.Po_logx      = 0;
PLOT.Mobility                   =   logical(0);         IN.M_logx       = 0;
PLOT.f80                        =   logical(0);         IN.f80_logx     = 0;
PLOT.f90                        =   logical(0);         IN.f90_logx     = 0;

%% RUN MAIN ROUTINE >>>
[SAVE] = SL_RadialProfile(IN,PLOT,SWITCH,STYLE,SAVE);

TEST.successful = true;
end

