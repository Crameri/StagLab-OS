
%%                                       PARAMETER DIMENSIONALIZATION 1.71
%
%                                                Fabio Crameri, 24.02.2019

function [SETUP,SWITCH] = f_Dimensions(SWITCH,FILE)

%% DEFAULTS
if ~isfield(SWITCH,'GeodynamicCode');       SWITCH.GeodynamicCode       = 'StagYY'; end
if ~isfield(SWITCH,'useReadInParameters');  SWITCH.useReadInParameters  = true;     end

%% STANDARD PARAMETERS
% MANUAL PARAMETERS (adjust here)

SETUP.D             	= 1;                % non-dimensional value; is later converted to dimensional value.
SETUP.secyear           = 3600*24*365.25;   % year to seconds conversion factor
SETUP.g                 = 9.81;             % [m/s2]
SETUP.deltaT            = 2500;             % [K]
SETUP.tcond_dimensional = 3.0;              % th. conductivity [W/(m K)]
SETUP.rho0              = 3300;             % reference density [kg/m3]
SETUP.cp_dimensional    = 1200.0;           % reference heat capacity [J/(kg K)]
SETUP.alpha             = 3e-5;             % th. expansivity [1/K]

% AUTOMATICALLY DETERMINE PARAMETERS (from refstat.dat if possible)
[SETUP] = f_determineParametersFromRefstat(SETUP,FILE,SWITCH);

%% SPECIFIC PARAMETERS
if SWITCH.parameter==1          %Standard *********************************
    SETUP.topBC         = 'free-slip';
    SETUP.Ra0        	= 1e6;
    SETUP.D_dim         = 2890e3;
   
elseif SWITCH.parameter==5      %Ra=1e5 ***********************************
    SETUP.topBC         = 'free-slip';
    SETUP.Ra0        	= 1e5;
    SETUP.D_dim         = 2890e3;
    
elseif SWITCH.parameter==6      %Ra=1e6 ***********************************
    SETUP.topBC         = 'free-slip';
    SETUP.Ra0        	= 1e6;
    SETUP.D_dim         = 2890e3;
    
elseif SWITCH.parameter==7      %Ra=1e7 ***********************************
    SETUP.topBC         = 'free-slip';
    SETUP.Ra0        	= 1e7;
    SETUP.D_dim         = 2890e3;
    
elseif SWITCH.parameter==8      %Ra=1e8 ***********************************
    SETUP.topBC         = 'free-slip';
    SETUP.Ra0        	= 1e8;
    SETUP.D_dim         = 2890e3;
    
elseif SWITCH.parameter==11
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0        	= 1.1639e6;
    SETUP.D_dim         = 3040e3;
    SETUP.d_air         = 0.05;
    
elseif SWITCH.parameter==12
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0        	= 1.1639e6;
    SETUP.D_dim         = 1445e3;
    SETUP.d_air         = 0.05;
    
elseif SWITCH.parameter==110  %old
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0        	= 1e6;
    SETUP.D_dim         = 2890e3;
    SETUP.d_air         = 0.05;
    
elseif SWITCH.parameter==111  % free surface, e.g., Fluidity
    SETUP.topBC         = 'free';
    SETUP.Ra0        	= 1.1639e6;
    SETUP.D_dim         = 2900e3;
    SETUP.d_air         = 0.0;
    
elseif SWITCH.parameter==13  %thin air layer ******************************
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0         	= 8.6533e6;
    SETUP.D_dim         = 3000e3;
    SETUP.d_air         = 0.05;
    
elseif SWITCH.parameter==4  %Solomatov 2004 *******************************
    SETUP.topBC         = 'free-slip';
    SETUP.Ra0           = 0.1;
    SETUP.D_dim         = 600e3;
    
elseif SWITCH.parameter==14
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0         	= 0.11664;
    SETUP.D_dim         = 631.5789e3;
    SETUP.d_air         = 0.05;
    
elseif SWITCH.parameter==114 %old
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0         	= 0.1;
    SETUP.D_dim         = 600e3;
    SETUP.d_air         = 0.05;

elseif SWITCH.parameter==9     
    SETUP.topBC         = 'free-slip';
    SETUP.Ra0        	= 1e7;
    SETUP.D_dim         = 660e3;
    
elseif SWITCH.parameter==15  %Ra=1e5 sticky-air ***************************
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0         	= 1.e5;
    SETUP.D_dim         = 3040e3;
    SETUP.d_air         = 0.05;

elseif SWITCH.parameter==55  %Ra=5e5 sticky-air ***************************
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0         	= 5.e5;
    SETUP.D_dim         = 3040e3;
    SETUP.d_air         = 0.05;

elseif SWITCH.parameter==56  %Ra=5e6 sticky-air ***************************
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0         	= 5.e6;
    SETUP.D_dim         = 3040e3;
    SETUP.d_air         = 0.05;
    
elseif SWITCH.parameter==17  %Ra=1e7 sticky-air ***************************
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0         	= 1.e7;
    SETUP.D_dim         = 3040e3;
    SETUP.d_air         = 0.05;
    
elseif SWITCH.parameter==18  %Ra=1e8 sticky-air ***************************
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0         	= 1.e8;
    SETUP.D_dim         = 3040e3;
    SETUP.d_air         = 0.05;
    
elseif SWITCH.parameter==121  %Venus **************************************
    SETUP.topBC         = 'free-slip';
    SETUP.D_dim         = 2942e3;
    SETUP.Ra0         	= 1.e7;
    
elseif SWITCH.parameter==131  %Venus Rebecca ******************************
    SETUP.topBC         = 'free-slip';
    SETUP.D_dim         = 2866e3;
    SETUP.Ra0         	= 1.1154e8;
    SETUP.g                 = 8.87;             % [m/s2]
    SETUP.deltaT            = 2300;             % [K]
    SETUP.tcond_dimensional = 4.0;              % th. conductivity [W/(m K)]
    SETUP.rho0              = 3378;             % reference density [kg/m3]
    SETUP.cp_dimensional    = 1250.0;           % reference heat capacity [J/(kg K)]
    SETUP.alpha             = 2e-5;             % th. expansivity [1/K]

    
elseif SWITCH.parameter==19	%Special Lab **********************************
    SETUP.topBC         = 'free-slip';
    SETUP.D_dim         = 0.265; %2890e3;
    SETUP.rho0              = 1439;     	% reference density in kg/m3
    %     SETUP.rho0              = 388.7339;  %105.0132; %3.8873 %Ra1e5
    %     SETUP.rho0              = 869.2354;  %Ra5e5
    %     SETUP.rho0              = 1.2293e3;  %1050.132; %3.8873e1 %Ra1e6
    %     SETUP.rho0              = 2.7488e3;  %Ra5e6
    %     SETUP.rho0              = 3.8873e3;  %10501.32; %3.8873e2 %Ra1e7
    %     SETUP.rho0              = 8.6924e3;  %Ra5e7
    %     SETUP.rho0              = 1.2293e4;  %105013.2; %3.8873e3 %Ra1e8
    T_heater            = 80.15;
    SETUP.tcond_dimensional = 0.34;
    %     T_heater            = 25.9685; %deta=1.1
    %     SETUP.tcond_dimensional = 0.0048;
    %     T_heater            = 30.9324; %deta=1.1
    %     SETUP.tcond_dimensional = 0.0355;
    %     T_heater            = 39.0967; %deta=5
    %     SETUP.tcond_dimensional = 0.086;
    %     T_heater            = 45.8315; %deta=10
    %     SETUP.tcond_dimensional = 0.1277;
    %     T_heater            = 64.2944; %deta=50
    %     SETUP.tcond_dimensional = 0.2419;
    %     T_heater            = 74.2838; %deta=100
    %     SETUP.tcond_dimensional = 0.3037;
    %     T_heater            = 121.7815; %deta=490
    %     SETUP.tcond_dimensional = 0.5976;
    T_ambient               = 25.2;
    T0                      = 25.2;
    SETUP.deltaT            = T_heater-T_ambient; %54.95; %already dimensional; in [K]
    SETUP.alpha             = 3.1e-4;
    SETUP.cp_dimensional    = 2280.0;           % reference heat capacity in J/kg/K
    SETUP.kappa         = SETUP.tcond_dimensional/(SETUP.rho0*SETUP.cp_dimensional);
    eta0_in             = 1080*exp(-0.156*T0+6.25e-4*T0^2); %32.3;
    %     eta0_in             = 31.5135; %isoviscous
    SETUP.Ra0        	= (SETUP.rho0*SETUP.g*SETUP.alpha*SETUP.deltaT*SETUP.D_dim^3)/(SETUP.kappa*eta0_in);
    
elseif SWITCH.parameter==20	%Special Lab 2 ********************************
    SETUP.topBC         = 'free-slip';
    SETUP.D_dim         = 0.265; %2890e3;
    SETUP.rho0              = 1439;     	% reference density in kg/m3
    %     SETUP.rho0              = 388.7339;  %105.0132; %3.8873 %Ra1e5
    %     SETUP.rho0              = 869.2354;  %Ra5e5
    %     SETUP.rho0              = 1.2293e3;  %1050.132; %3.8873e1 %Ra1e6
    %     SETUP.rho0              = 2.7488e3;  %Ra5e6
    %     SETUP.rho0              = 3.8873e3;  %10501.32; %3.8873e2 %Ra1e7
    %     SETUP.rho0              = 8.6924e3;  %Ra5e7
    %     SETUP.rho0              = 1.2293e4;  %105013.2; %3.8873e3 %Ra1e8
    T_heater            = 80.15;
    SETUP.tcond_dimensional = 0.34e-6;
    %     T_heater            = 25.9685; %deta=1.1
    %     SETUP.tcond_dimensional = 0.0048;
    %     T_heater            = 30.9324; %deta=1.1
    %     SETUP.tcond_dimensional = 0.0355;
    %     T_heater            = 39.0967; %deta=5
    %     SETUP.tcond_dimensional = 0.086;
    %     T_heater            = 45.8315; %deta=10
    %     SETUP.tcond_dimensional = 0.1277;
    %     T_heater            = 64.2944; %deta=50
    %     SETUP.tcond_dimensional = 0.2419;
    %     T_heater            = 74.2838; %deta=100
    %     SETUP.tcond_dimensional = 0.3037;
    %     T_heater            = 121.7815; %deta=490
    %     SETUP.tcond_dimensional = 0.5976;
    T_ambient               = 25.2;
    T0                      = 25.2;
    SETUP.deltaT            = T_heater-T_ambient; %54.95; %already dimensional; [K]
    SETUP.alpha             = 3.1e-10;
    SETUP.cp_dimensional    = 2280.0;       % reference heat capacity in J/kg/K
    SETUP.kappa         = SETUP.tcond_dimensional/(SETUP.rho0*SETUP.cp_dimensional);
    eta0_in             = 1080*exp(-0.156*T0+6.25e-4*T0^2); %32.3;
    %     eta0_in             = 31.5135; %isoviscous
    SETUP.Ra0        	= (SETUP.rho0*SETUP.g*SETUP.alpha*SETUP.deltaT*SETUP.D_dim^3)/(SETUP.kappa*eta0_in);
   
elseif SWITCH.parameter==21
    SETUP.topBC         = 'sticky-air';
    SETUP.Ra0        	= 1.1639e7;
    SETUP.D_dim         = 3040e3;
    SETUP.d_air         = 0.05;
    
else  % *******************************************************************
    error(['SWITCH.parameter = ',num2str(SWITCH.parameter),' not found!'])
end



%% ADDITIONAL SPECIFIC PARAMETER SUITES
if strcmp(SWITCH.modelSetup,'standard')
    SETUP.model             = 'STANDARD';
    
elseif strcmp(SWITCH.modelSetup,'topo_bm')
    SETUP.model             = 'TOPO BENCHMARK';
    SETUP.D_dim             = 775e3;        % 850e3;	 % Depth in meters !!!
    SETUP.Ra0               = 1.5361;       % 2.0266
    SETUP.d_air             = 0.096774;     %just to define it...
    SETUP.g                 = 10;           % m/s^2
    SETUP.tcond_dimensional = 3.3;          % W/m/K
    SETUP.cp_dimensional    = 1000.0;       % J/K/kg   
    SETUP.alpha             = 1e-5;         % 1/K
    SETUP.deltaT            = 0.01;         % K
    
elseif strcmp(SWITCH.modelSetup,'solomatov2004')
    SETUP.model             = 'SOLOMATOV 2004';
    %     H_dummy = 6.3131e-12;
    %     SETUP.Tscale            = SETUP.rho0*H_dummy*SETUP.D_dim^2/SETUP.tcond_dimensional; 2500; [K]
    SETUP.deltaT            = 2500;
    
else
    error(['SWITCH.modelSetup = ',SWITCH.modelSetup,' not found!'])
end


%% USE READ-IN PARAMETERS
if SWITCH.useReadInParameters
    if ~isempty(SETUP.D_dim0)
        if SETUP.D_dim0~=SETUP.D_dim
            warning off backtrace
            warning(['The model''s box depth seems to be ',num2str(SETUP.D_dim0/1e3),...
                ' km deep, but SETUP.D_dim in f_Dimensions.m is set to ',num2str(SETUP.D_dim/1e3),...
                ' km: Check your IN.Parameter settings!'])
            warning on backtrace
        end
        
        SETUP.D_dim      	= SETUP.D_dim0;     % [m]
        
    end
end

%% CALCULATING ADDITIONAL PARAMETERS
if SWITCH.DimensionalMode || SWITCH.DimensionalInput
    SETUP.D = SETUP.D*SETUP.D_dim;
end
SETUP.kappa         = SETUP.tcond_dimensional/(SETUP.rho0*SETUP.cp_dimensional); % reference thermal diffusivity in m2/s
SETUP.eta0          = (SETUP.rho0*SETUP.g*SETUP.alpha*SETUP.deltaT*SETUP.D^3) / (SETUP.kappa*SETUP.Ra0);  % reference viscosity (varies with Ra0)
SETUP.v0            = SETUP.kappa/SETUP.D_dim *SETUP.secyear *100;                              % [cm/a]

%% CALCULATING PARAMETER SCALES
SETUP.lengthscale   = SETUP.D;                                      SETUP.lengthDim = 'm';     	% [m]
SETUP.areascale     = SETUP.D^2;                                    SETUP.areaDim = 'm^2';  	% [m^2]
SETUP.volumescale   = SETUP.D^3;                                 	SETUP.volumeDim = 'm^3';  	% [m^3]
SETUP.timescale     = SETUP.D_dim^2/SETUP.kappa;                    SETUP.timeDim = 's';     	% [s]
SETUP.Compscale   	= 100;                                          SETUP.CompDim = 'wt%';    	% [wt%]  composition
SETUP.stressscale	= SETUP.eta0*SETUP.kappa/SETUP.D_dim^2 /1e6;    SETUP.stressDim = 'MPa';	% [MPa]
SETUP.edotscale  	= SETUP.kappa/SETUP.D_dim^2;                    SETUP.edotDim = 's^{-1}';  	% [1/s]
                                                                    SETUP.edotDim2 = '1/s';
SETUP.dissscale     = SETUP.eta0*SETUP.kappa/SETUP.D_dim^2*SETUP.edotscale;         	SETUP.dissDim = 'W/m^3';    % [W/m^3]
SETUP.Tscale        = SETUP.deltaT;                                 SETUP.TDim = 'K';           % [K]
SETUP.etascale      = SETUP.eta0;                                   SETUP.etaDim = 'Pas';       % [Pa s]
SETUP.Pscale        = SETUP.stressscale/1e3;                        SETUP.PDim = 'GPa';      	% [GPa]
SETUP.Hscale        = (SETUP.kappa*SETUP.cp_dimensional*SETUP.Tscale)/SETUP.D_dim^2;    SETUP.HDim = 'W/kg';    % [W/kg]
SETUP.HRscale       = SETUP.deltaT/SETUP.D_dim^2/SETUP.kappa;       SETUP.HRDim = 'K/s';      	% [K/s]
SETUP.HFscale       = (SETUP.tcond_dimensional*SETUP.Tscale)/SETUP.D_dim *1e3;          SETUP.HFDim = 'mW/m^2';	% [mW/m2] this is milliwatt!
SETUP.Vscale        = SETUP.v0;                                     SETUP.vDim = 'cm/a';      	% [cm/a]
SETUP.streamFscale 	= 1;                                            SETUP.streamFDim = 'cm^2/s';% [cm^2/s]  conversion from dimensionalised velocity
SETUP.rhoscale     	= SETUP.rho0;                                   SETUP.rhoDim = 'kg/m^3';  	% [kg/m3]
SETUP.wscale        = SETUP.edotscale*SETUP.secyear*1e6;            SETUP.wDim = 'Ma^{-1}'; 	% [1/Ma]  vorticity scale
SETUP.divscale    	= SETUP.edotscale*SETUP.secyear*1e6;            SETUP.divDim = 'Ma^{-1}';	% [1/Ma]  divergence scale
SETUP.clapeyronscale= (SETUP.rho0*SETUP.g*SETUP.D_dim)/SETUP.deltaT; SETUP.clapeyronDim = 'kgm^{-1}s^{-2}K^{-1}';	% [kg/m/s^2/K]  clapeyron scale
SETUP.massscale     = SETUP.rhoscale*SETUP.D_dim^3;                 SETUP.massDim = 'kg';       % [kg]
SETUP.eruptratescale= SETUP.rhoscale*SETUP.D_dim/SETUP.timescale*SETUP.secyear*1e6;	SETUP.eruptrateDim = 'kgm^{-2}Ma^{-1}'; 	% [kg m^{-2} Ma^{-1}]


%% ADJUSTING VALUES
if strcmp(SETUP.topBC,'free-slip') %for free-slip cases
    SETUP.d_air         = 0.0;      %non-dimensional
    SETUP.d_air_dim     = 0.0;      %dimensional [m]
else %for sticky-air cases
    if ~isfield(SETUP,'d_air_dim'); SETUP.d_air_dim = SETUP.d_air*SETUP.D_dim; end %dimensional [m]
end

%NON-DIMENSIONAL CASE
if ~SWITCH.DimensionalMode %adjust for non-dimensional mode
    SETUP.etascale    	= 1;                        SETUP.etaDim = 'nd';              	% [nd]
    SETUP.lengthscale 	= 1;                        SETUP.lengthDim = 'nd';           	% [nd]
    SETUP.areascale     = 1;                     	SETUP.areaDim = 'nd';               % [nd]
    SETUP.volumescale   = 1;                     	SETUP.volumeDim = 'nd';             % [nd]
    SETUP.timescale     = 1;                        SETUP.timeDim = 'nd';           	% [nd]
    SETUP.stressscale	= 1;                        SETUP.stressDim = 'nd';           	% [nd]
    SETUP.edotscale  	= 1;                        SETUP.edotDim = 'nd';               % [nd]
                                                    SETUP.edotDim2 = 'nd';  
  	SETUP.massscale     = 1;                        SETUP.massDim = 'nd';               % [nd]
    SETUP.dissscale     = 1;                        SETUP.dissDim = 'nd';               % [nd]
    SETUP.Pscale        = 1;                        SETUP.PDim = 'nd';                	% [nd]
    SETUP.Hscale        = 1;                        SETUP.HDim = 'nd';                 	% [nd]
    SETUP.HRscale       = 1;                        SETUP.HRDim = 'nd';                	% [nd]
    SETUP.HFscale       = 1;                        SETUP.HFDim = 'nd';              	% [nd]
    SETUP.Vscale        = 1;                        SETUP.vDim = 'nd';              	% [nd]
    SETUP.streamFscale 	= 1;                        SETUP.streamFDim = 'nd';            % [nd]
    SETUP.Tscale        = 1;                        SETUP.TDim = 'nd';                	% [nd]
    SETUP.rhoscale     	= 1;                        SETUP.rhoDim = 'nd';              	% [nd]
 	SETUP.wscale        = 1;                        SETUP.wDim = 'nd';                  % [nd]  vorticity scale
    SETUP.divscale    	= 1;                        SETUP.divDim = 'nd';                % [nd]  divergence scale
    SETUP.clapeyronscale= 1;                        SETUP.clapeyronDim = 'nd';      	% [nd]  clapeyron scale
    SETUP.eruptratescale= 1;                        SETUP.eruptrateDim = 'nd';      	% [nd]  eruption-rate scale
    SETUP.Compscale  	= 1;                     	SETUP.CompDim = 'nd';              	% [nd]  composition
end

%DIMENSIONAL INPUT
if SWITCH.DimensionalInput
    SETUP.d_air         = SETUP.d_air_dim;
    SETUP.etascale     	= 1;                        SETUP.etaDim = 'Pas';             	% [Pa s]
    SETUP.lengthscale 	= 1;                        SETUP.lengthDim = 'm';           	% [m]
    SETUP.areascale     = 1;                     	SETUP.areaDim = 'm^2';              % [m^2]
    SETUP.volumescale   = 1;                      	SETUP.volumeDim = 'm^3';            % [m^3]
    SETUP.timescale     = 1;                        SETUP.timeDim = 's';                % [s]
    SETUP.stressscale	= 1/1e6;                    SETUP.stressDim = 'MPa';          	% [MPa]
    SETUP.edotscale  	= 1;                        SETUP.edotDim = 's^{-1}';         	% [1/s]
                                                    SETUP.edotDim2 = '1/s';
  	SETUP.massscale     = 1;                        SETUP.massDim = 'kg';               % [kg]
    SETUP.dissscale     = 1;                        SETUP.dissDim = 'W/m^3';          	% [W/m^3]
    SETUP.Pscale        = 1/1e9;                    SETUP.PDim = 'GPa';              	% [GPa]
    SETUP.Hscale        = 1;                        SETUP.HDim = 'W/kg';             	% [W/kg]
    SETUP.HRscale     	= 1;                        SETUP.HRDim = 'K/s';                % [K/s]
    SETUP.HFscale       = 1*1e3;                    SETUP.HFDim = 'mW/m^2';          	% [mW/m2]
    SETUP.Vscale        = SETUP.secyear*100;        SETUP.vDim = 'cm/a';                % [cm/a]
    SETUP.streamFscale 	= 1;                        SETUP.streamFDim = 'cm^2/s';        % [cm^2/s] conversion from dimensionalised velocity
    SETUP.Tscale        = 1;                        SETUP.TDim = 'K';                	% [K]
    SETUP.rhoscale     	= 1;                        SETUP.rhoDim = 'kg/m^3';          	% [kg/m3]
    SETUP.wscale        = 1*SETUP.secyear*1e6;      SETUP.wDim = 'Ma^{-1}';             %[1/s]->[1/Ma]  vorticity scale
    SETUP.divscale    	= 1*SETUP.secyear*1e6;      SETUP.divDim = 'Ma^{-1}';       	%[1/s]->[1/Ma]  divergence scale
    SETUP.clapeyronscale= 1;                        SETUP.clapeyronDim = 'kgm^{-1}s^{-2}K^{-1}';	% [kg/m/s^2/K]
    SETUP.eruptratescale= 1;                        SETUP.eruptrateDim = 'kgm^{-2}Ma^{-1}';      	% [kg m^{-2} Ma^{-1}]  eruption-rate scale
    SETUP.Compscale   	= 100;                     	SETUP.CompDim = 'wt%';             	% [wt%]  composition
end

%% PRINT INFORMATION
if SWITCH.DimensionalMode
    disp(['Dim. Parameters:   ',SETUP.model])
    if SWITCH.DimensionalMode || SWITCH.DimensionalInput
        disp(['D:                 ',num2str(SETUP.D/1e3),' km'])
    else
        disp(['D:                 ',num2str(SETUP.D)])
    end
    disp(['Ra0:               ',num2str(SETUP.Ra0,3)])
end

end




%INTERNAL FUNCTIONS

function [SETUP] = f_determineParametersFromRefstat(SETUP,FILE,SWITCH)

SETUP.D_dim0 = [];
if SWITCH.useReadInParameters && exist('FILE','var') && strcmpi(SWITCH.GeodynamicCode,'StagYY')
    try
        filestem = [FILE.directory,FILE.name,'_refstat.dat'];
        if exist(filestem,'file') && SWITCH.DimensionalInput
            % READ PARAMETERS FROM REFSTAT.DAT (if available)  THOSE ARE EQUAL 1 IN NON-DIMENSIONAL MODE
            PARA.Task       	= 'ImportRefstat';
            [PARA,~] = f_Import(PARA,FILE);
            SETUP.D_dim0            = PARA.z(1,1);      % domain depth (D_dim) [m]
            SETUP.tcond_dimensional = PARA.Tcond(1,1); 	% th. conductivity [W/(m K)]
            SETUP.rho0              = PARA.rho(1,1);   	% reference density [kg/m3]
            SETUP.cp_dimensional    = PARA.Cp(1,1);    	% reference heat capacity [J/(kg K)]
            SETUP.alpha             = PARA.expan(1,1);	% th. expansivity [1/K]
            %     SETUP.deltaT            = PARA.T_ref; % [K]
        end
    catch me
        if SWITCH.Verbose; warning(me.message); end
    end
end

end

