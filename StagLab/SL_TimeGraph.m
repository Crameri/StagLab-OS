
%%                                                      SL_TimeGraph 10.001
%
%                                                  plots StagYY's time.dat
%    . calls f_DesignBackground
%    . calls f_DefaultsTimedat
%    . calls f_Dimensions
%    . calls f_DesignFigure
%    . calls f_FileFinder
%    . calls f_saveData
%    . calls f_Import
%    . calls f_saveFigure
%    . calls f_SetupFigure
%    . calls f_DesignVaria
%                                               22.05.2020 : Fabio Crameri

%% UPDATES >10.0
% *preventing annotating separate legend panel*

%% UPDATES 10.0
% *introducing panel*
% *introducing automated error logging*
% *introducing categorical Scientific Colour Maps*
% *introducing journal-specific plot design*
% *special character fix for Windows*

%% UPDATES 9.0
% *introducing analysis mode*
% *option to write out min, mean and max values*
% *option to plot logarithmic y-axis*
% *option to place legend separately*
% *automatic conversion of mass fractions to wt% or ppm*
% *option for manual x- and y-axis limits*
% *additional field parameters*
% *automatic data-field detection*
% *compatibility improvements to new StagYY version*
% *automatic fixing of corrupt time.dat files*
% *bug fixes*

%% UPDATES 8.0
% *additional parameters added*
% *deeper fAIo integration*
% *using improved fileFinder*
% *more flexible file read*
% *automatic data size limitation*
% *more informative file naming*
% *option to name files manually*
% *option to handle time.dat files without header*
% *option to save processed graph data*
% *option to set default screen position of figure*
% *option for quick mode*
% *speed improvements*
% *bug fixes*

%% UPDATES 7.0
% *more stable design-routines*
% *cleaner plot design and layout*
% *automatic figure size*

%% UPDATES 6.0
% *improved plot layout*
% *improved code design*
% *bug fixes to dimensionalisation*

function [SAVE] = SL_TimeGraph(IN,PLOT,SWITCH,STYLE,SAVE)
if exist('IN','var') %if run from a parfile
    %% DEFAULT VARIABLES
    %...see f_DefaultsTimedat.m
else
    clear;
    %% INPUT (...use parameter file instead)
    SAVE.StartNumber = 1; SAVE.EndNumber = 1;
    [IN,PLOT,SWITCH,STYLE,SAVE] = f_DefaultsTimedat; %get default variables
    IN.Name         =   {  'testDIM'  'testDIM' };	% filename
    IN.Parameter    =   [   11        	];
    IN.Folder       =   {   '~/work/plotTest'   };
end

%% STARTUP PROCEDURE
SAVE.app = 'SL_TimeGraph'; SAVE.appVersion = 10.0; %only 3 digits to prevent unecessary parfile updates
if ~isfield(SAVE,'count'); SAVE.count = 0; end
DESVARIA.Task    = 'set special characters';
[~,DESVARIA,STYLE,~] = f_DesignVaria(DESVARIA,[],[],[],[],STYLE,[]);
[fAIo,SAVE,SWITCH] = f_Startup(STYLE,SAVE,SWITCH);
if ~strcmp(fAIo.Status,'all fine!'); disp(' '); disp(['                         ',STYLE.SCHAR.checkMark,' Finished.']); return; end

try 
    %% VARIABLE ADJUSTMENTS
    nrModels = length(IN.Name);
    SAVE.BulletString     	= STYLE.SCHAR.hugeBulletLight;
    SAVE.PointingString    	= STYLE.SCHAR.downrightArrow;
    if SWITCH.AnalysisMode %optimise plots for accurate analysis
        SWITCH.PlotDesign           = true;
        STYLE.AllFontName           = 'Helvetica';
        STYLE.keyFontSize           = 9;
        SWITCH.SimplifyPlots      	= true;
        SWITCH.BackgroundDesign    	= false;
        SWITCH.LegendDesign        	= true;
        SWITCH.Annotation           = false;
        STYLE.keyColor              = [0 0 0];
        STYLE.MinorAxisTicks        = true;
    end
    if SWITCH.QuickMode
        SAVE.Figure                  	=   false;
        SWITCH.PlotDesign               =   false;
        SWITCH.BackgroundDesign     	=   false;
        SWITCH.LegendDesign          	=   false;
        SWITCH.SimplifyPlots           	=   false;
        SWITCH.Annotation               =   false;
    end
    DESVARIA.Task = 'set plot design mode'; %set journal design mode
    [PLOT,DESVARIA,STYLE,~] = f_DesignVaria(DESVARIA,[],[],[],PLOT,STYLE,[]);
    if SWITCH.sendErrorLog
        SWITCH.Verbose              = true;
    end
    if isfield(PLOT,'titleString') && length(PLOT.titleString)>nrModels
        PLOT.titleString(:,nrModels+1:end)  = [];
    end
    % check for length of number input arrays
    if length(IN.Parameter)<nrModels
        dummy = ones(1,nrModels-length(IN.Parameter)); dummy(:) = IN.Parameter(1,1);
        IN.Parameter = [IN.Parameter, dummy];
    end
    % check for length of string input arrays
    if length(IN.Folder)<nrModels
        dummy = cell(1,nrModels-length(IN.Folder)); dummy(:) = IN.Folder(1,1);
        IN.Folder = [IN.Folder, dummy];
    end
    
    %% GET NUMBER OF SUBPLOTS
    nrPlots=0;                                                               % timedat label
    if PLOT.Temperature;            nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'Tmean'}; end
    if PLOT.Viscosity;              nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'eta_mean'}; end
    if PLOT.Velocity;               nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'Vrms'}; end
    if PLOT.Ra_eff;                 nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'ra_eff'}; end
    if PLOT.Nu_top;                 nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'Nu_top'}; end
    if PLOT.Nu_bot;                 nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'Nu_bot'}; end
    if PLOT.InternalHeating;        nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'H_int'}; end
    if PLOT.Composition;            nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'C_mean'}; end
    if PLOT.HeatfluxTop;            nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'F_top'}; end
    if PLOT.HeatfluxBot;            nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'F_bot'}; end
    if PLOT.EruptionRate;           nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'erupt_rate'}; end
    if PLOT.EruptHeatFlux;          nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'erupt_heatflux'}; end
    if PLOT.Erupta;                 nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'erupta'}; end
    if PLOT.Entrainment;            nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'entrainment'}; end
    if PLOT.Cmass_error;            nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'Cmass_error'}; end
    %new options - need to be adjusted & tested !
    if PLOT.TotalErupta;            nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'erupta_total'}; end
    if PLOT.OutgassedWater;         nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'outgassed_water'}; end
    if PLOT.InnerCoreRadius;        nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'r_innercore'}; warning([DATA.tag{1,nrPlots},' is a new field and still needs to be adjusted!']); end
    if PLOT.SurfaceTemperature;     nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'Tsurf'}; warning([DATA.tag{1,nrPlots},' is a new field and still needs to be adjusted!']); end
    if PLOT.CMB_Temperature;        nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'Tcmb'}; warning([DATA.tag{1,nrPlots},' is a new field and still needs to be adjusted!']); end
    if PLOT.Psurf;                  nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'Psurf'}; warning([DATA.tag{1,nrPlots},' is a new field and still needs to be adjusted!']); end
    if PLOT.Weathering;             nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'weathering'}; warning([DATA.tag{1,nrPlots},' is a new field and still needs to be adjusted!']); end
    if PLOT.TotalWater;             nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'mH2O_total'}; end
    if PLOT.MantleWater;            nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'mH2O_mantle'}; end
    if PLOT.CrustWater;             nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'mH2O_crust'}; end
    if PLOT.dmH2O;                  nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'dmH2O'}; warning([DATA.tag{1,nrPlots},' is a new field and still needs to be adjusted!']); end
    if PLOT.OutgassedCarbon;        nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'outgassed_carbon'}; warning([DATA.tag{1,nrPlots},' is a new field and still needs to be adjusted!']); end
    if PLOT.OutgassedNitrogen;      nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'outgassed_nitrogen'}; warning([DATA.tag{1,nrPlots},' is a new field and still needs to be adjusted!']); end
    if PLOT.s_core;                 nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'s_core'}; warning([DATA.tag{1,nrPlots},' is a new field and still needs to be adjusted!']); end
    if PLOT.tc_core;                nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'tc_core'}; warning([DATA.tag{1,nrPlots},' is a new field and still needs to be adjusted!']); end
    if PLOT.ts_core;                nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'ts_core'}; warning([DATA.tag{1,nrPlots},' is a new field and still needs to be adjusted!']); end
    if PLOT.OutgassedAmountWater;   nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'outgassed_amount_H2O'}; end
    if PLOT.OutgassedWaterExtra;	nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'outgassed_extra_amountH2O'}; end
    if PLOT.H2Odegassing;          	nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'wflux_d'}; end
    if PLOT.H2Oregassing;         	nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'wflux_r'}; end
    %these have to be last
    if PLOT.Timesteps;              nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'istep'}; end
    if PLOT.Performance;            nrPlots = nrPlots+1; DATA.tag(1,nrPlots)={'istep'}; end
    
    %extra axes for outside legend
    if SWITCH.LegendOutside
        nrPlots = nrPlots+1;
        PLOT.annotationString{nrPlots} = ''; %no annotation for extra legend panel (which is called second last)
    end
    
    %% LAYOUT SUBPLOTS
    if nrPlots<=4
        nxSubplots  = 1;
        nzSubplots  = nrPlots;
    else
        nxSubplots  = min(nrPlots,2);
        nzSubplots  = floor((nrPlots+1)/2);
    end
    
    %% SETUP COLOURS
    if strcmp(STYLE.ColorMode,'light')
        STYLE.ColorModeBW           = 'white';
        STYLE.BlackOrWhiteColor   	= [1 1 1];
    elseif strcmp(STYLE.ColorMode,'dark')
        STYLE.ColorModeBW           = 'black';
        STYLE.BlackOrWhiteColor   	= [0 0 0];
        STYLE.keyColor            	= 1-STYLE.keyColor; %invert key colours
    else
        warning(['Colour Mode ''',STYLE.ColorMode,''' not recognised!'])
        STYLE.ColorMode             = 'light';
        STYLE.ColorModeBW           = 'white';
    end
    
    %% SETUP FIGURE
    if ~SWITCH.QuickMode
        clf(figure(1));
        set(gcf,'Color',STYLE.ColorModeBW);
        % set figure position
        if strcmp(SWITCH.FigurePosition,'auto') %set figure position automatically
            FPOS.Default            = PLOT.FigureDefaultPosition;
            FPOS.subplotLayout      = [nzSubplots nxSubplots];
            [~] = f_SetupFigure(FPOS,SAVE,1,1);
        else %set figure position manually
            set(gcf,'Position',SWITCH.FigurePosition);
        end
        %set panel
        if SWITCH.UsePanel
            SAVE.P	= panel();
            SAVE.P.pack(nzSubplots,nxSubplots);   %set subplot layout
            % set margins
%             SAVE.P.de.margin = 42; %descendence
%             SAVE.P(1,1).marginbottom = 12;
%             SAVE.P(2).marginleft = 20;
            SAVE.P.margin = [26 23 16 16];  %figure edge margins: [left bottom right top]
            SAVE.P.marginbottom = 23;
            if SWITCH.Verbose; warning('Panel routine still in testing: Margins might need to be adjusted here'); end
        end
        
        %% COLOURING ADJUSTMENTS
        %Setup colour map
        CMAPPING.Task   = 'Categorical';
        CMAPPING.NumberColours = size(IN.Name,2);
        [SWITCH,CMAPPING] = f_DesignColourmap(SWITCH,PLOT,[],STYLE,CMAPPING);
        PLOT.ColorVector    = CMAPPING.ColourVector;
    end
    
    h_runs = zeros(size(IN.Name,2),1);
    for irun=1:size(IN.Name,2)
        clearvars -except IN PLOT SWITCH STYLE SAVE DATA irun h_runs nrPlots nxSubplots nzSubplots
        
        %% VARIABLE CONVERSION & SETUP
        FILE.name         	= IN.Name{irun};
        FILE.stemRead     	= IN.Folder{irun};
        FILE.stemSave     	= IN.Folder{irun};
        SWITCH.parameter    = IN.Parameter(irun);
        lineColor           = PLOT.ColorVector(irun,:);
        
        %% CHECK FOR INPUT FILE
        [FILE] = f_FileFinder(FILE,SWITCH,SAVE,STYLE);
        if FILE.NotFound
            SAVE.LastFile           = true;
            close all
            disp([FILE.NotFoundError])
            return
            
        end
        
        %% OPEN FILE AND CHECK FILE SIZE
        if SWITCH.closeOldWaitbars; wbOld = findall(0,'tag','TMWWaitbar'); delete(wbOld); end %close open ghost waitbars
        wb = waitbar(0,'Please wait...');
        %find header entries in first line automatically
        fid             = fopen(FILE.found,'r');
        iCol=0;
        entry=1;
        while (entry~='0')
            entry               = fscanf(fid,'%s',[1 1]);
            iCol                = iCol+1;
            DATA.title{1,iCol}  = entry;
        end
        fclose(fid);
        DATA.title(:,end)   = []; %remove the last entry (i.e., 0)
        waitbar(1/6,wb)
        %find number of columns
        fid             = fopen(FILE.found,'r');
        nColumns=0; iLine=0;
        while nColumns<10 %line has less than 10 entries/columns (i.e., is empty)
            iLine       = iLine+1;
            nColumns  	= numel(regexp(fgetl(fid),'\s*([^\s]*)\s*'));
        end
        fclose(fid);
        waitbar(2/6,wb)
        
        %% IMPORT & READING DATA
        DATA.filename       = FILE.found;
        DATA.titleRow    	= true;
        DATA.numberEntries  = nColumns;
        DATA.startRow       = iLine;
        if isempty(DATA.title)
            DATA.titleRow 	= false;
        end
        DATA.Task           = 'ImportTimedat';
        [DATA,~] = f_Import(DATA,[],[],[],SWITCH);
        waitbar(5/6,wb)
        nRows               = size(DATA.Array{1},1);
        timeRow             = find(strcmp(DATA.title,'time'));
        if ~DATA.titleRow
            timeRow         = 2;
        end
        waitbar(6/6,wb)
        close(wb)
        
        %% DISPLAY MODEL DETAILS I
        if SWITCH.Verbose
            % disp(['datafile is ',num2str(n),' x ',num2str(nLines)])
            disp(['Data Structure:       __',num2str(nColumns),'__>'])
            disp('                     |')
            disp(['                   ',num2str(nRows)])
            disp('                     |')
            disp('                     V')
        end
        
        %% DATA SIZE ADJUSTMENTS
        nRowsAdjusted = nRows;
        if SWITCH.parameterRange %this is especially slow for big data sets
            numDataPointsThreshold	= 5e3;
        else
            numDataPointsThreshold	= 5e4;
        end
        if nRows>numDataPointsThreshold
            while nRowsAdjusted>numDataPointsThreshold
                for icol=1:size(DATA.Array,2)
                    DATA.Array{1,icol}(1:2:end-1,:)    = []; %remove every second entry
                end
                nRowsAdjusted     	= size(DATA.Array{1,icol}(:),1);
            end
            disp(['                   ',STYLE.SCHAR.downrightArrow,' Data size reduced to ',num2str(nRowsAdjusted)])
        end
        
        %% fAIO ACTION: UPDATE CACHE
        fAIo.task   = 'handleCache';
        [fAIo,SAVE] = f_AIo(fAIo,SWITCH,FILE,[],STYLE,SAVE);
        
        %% CHECK FOR DIMENSIONAL INPUT:
        %get time after first timestep
        dummy = DATA.Array{1,timeRow}(2,1);
        if dummy>1e2  %if first timestep > 100 (years) then DIMENSIONAL INPUT
            SWITCH.DimensionalInput = true; dimString = 'dimensional';
            SWITCH.DimensionalMode = true; %set output to dimensional if input is dimensional
        else
            SWITCH.DimensionalInput = false; dimString = 'non-dimensional';
        end
        disp(['Data:              ',SWITCH.GeodynamicCode,', ',dimString])
        
        %% DIMENSIONALIZATION
        [SETUP,SWITCH] = f_Dimensions(SWITCH,FILE);
        
        % SET TIME
        PLOT.time         	= DATA.Array{1,timeRow};
        PLOT.time_dim     	= DATA.Array{1,timeRow}*SETUP.timescale/SETUP.secyear;  %get dimensional time [yr]
        
        dummy = max(PLOT.time_dim); %[yr]
        if dummy>100*1e6 %>100Myr
            PLOT.time2plot = PLOT.time_dim/1e9; PLOT.time2plotDim = 'Gyr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear/1e9;
        elseif dummy<1e-3 %<0.001yr
            PLOT.time2plot = PLOT.time_dim*SETUP.secyear; PLOT.time2plotDim = 's'; PLOT.timeConvert = SETUP.timescale;
        elseif dummy<1e2 %<100yr
            PLOT.time2plot = PLOT.time_dim; PLOT.time2plotDim = 'yr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear;
        elseif dummy<1e5 %<0.1Myr
            PLOT.time2plot = PLOT.time_dim/1e3; PLOT.time2plotDim = 'kyr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear/1e3;
        else
            PLOT.time2plot = PLOT.time_dim/1e6; PLOT.time2plotDim = 'Myr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear/1e6;
        end
        if ~SWITCH.DimensionalMode && ~SWITCH.DimensionalInput; PLOT.time2plot = PLOT.time; PLOT.time2plotDim = 'nd'; PLOT.timeConvert = 1; end
        
        DATA.Array{1,timeRow}       = PLOT.time2plot; %dimensionalise time array
        
        minTime                     = min(PLOT.time); %[nd or s]
        maxTime                     = max(PLOT.time); %[nd or s]
        minTimeDim                  = minTime*PLOT.timeConvert; %[plotting dimension]
        maxTimeDim                  = maxTime*PLOT.timeConvert; %[plotting dimension]
        
        %% SET UP FIELD VARIABLES
        %  	1                       2                   3                   4                   5                   6              	7
        %  	data name               variable label      time.dat label    	dim. scale      	data dimension      indiv. range 	log Y-axis
        FIELD_M = {
            'Time steps'            'cumul. time steps' 'istep'          	1                   '#'                 false           false
            'Time'                  'time'              'time'           	PLOT.timeConvert    PLOT.time2plotDim   false           false
            'Top heatflux'          'F_{top}'           'F_top'             SETUP.HFscale      	SETUP.HFDim       	false         	false
            'Bottom heatflux'       'F_{bot}'           'F_bot'             SETUP.HFscale    	SETUP.HFDim       	false        	false
            'Min. temperature'      'T_{min}'           'Tmin'              SETUP.Tscale        SETUP.TDim       	false          	false
            'Temperature'           'T_{mean}'          'Tmean'             SETUP.Tscale        SETUP.TDim       	true          	false
            'Max. temperature'      'T_{max}'           'Tmax'              SETUP.Tscale        SETUP.TDim       	false         	false
            'Min. velocity'         'v_{min}'           'Vmin'              SETUP.Vscale        SETUP.vDim       	false        	IN.v_logy
            'Velocity'              'v_{mean}'        	'Vrms'              SETUP.Vscale        SETUP.vDim       	true           	IN.v_logy
            'Max. velocity'         'v_{max}'           'Vmax'              SETUP.Vscale        SETUP.vDim       	false        	IN.v_logy
            'Min. viscosity'        '\eta_{min}'        'eta_min'         	SETUP.etascale    	SETUP.etaDim     	false         	IN.eta_logy
            'Viscosity'             '\eta_{mean}'       'eta_mean'        	SETUP.etascale     	SETUP.etaDim     	true         	IN.eta_logy
            'Max. viscosity'        '\eta_{max}'        'eta_max'         	SETUP.etascale    	SETUP.etaDim      	false           IN.eta_logy
            'Ra_{eff}'              'Ra_{eff}'          'ra_eff'           	1                   ''                  false        	IN.Ra_logy
            'Nu_{top}'              'Nu_{top}'          'Nu_top'          	1                   ''                  false          	false
            'Nu_{bot}'              'Nu_{bot}'          'Nu_bot'           	1                   ''                  false          	false
            'Min. composition'      'C_{min}'           'C_min'             100              	'wt%'              	false         	false
            'Composition'           'C_{mean}'        	'C_mean'          	100              	'wt%'             	true          	false
            'Max. composition'      'C_{max}'          	'C_max'             100              	'wt%'              	false        	false
            'Mean melt fraction'    'F_{mean}'       	'F_mean'           	1                   ''                  false         	false
            'Max. melt fraction'    'F_{max}'       	'F_max'             1                   ''                  false         	false
            'Eruption rate'         'UCvarepsilondot'  	'erupt_rate'      	SETUP.eruptratescale SETUP.eruptrateDim	false         	false
            'Erupta'                'erupta'            'erupta'          	SETUP.rhoscale      SETUP.rhoDim    	false         	false
            'Eruption heatflux'     'hf_{erupt}'        'erupt_heatflux'   	SETUP.HFscale    	SETUP.HFDim       	false         	false
            'Entrainment'           'en'                'entrainment'      	1                   ''                  false        	false
            'Conserv. error'        'Cmass_{err}'       'Cmass_error'     	1                   ''                  false        	false
            'Internal heating'      'H_{int}'           'H_int'          	SETUP.Hscale        SETUP.HDim       	false        	false
            ...%additions in latest StagYY: (need to be completed)
            'Inner-core radius'     'R_{InnerCore}'   	'r_innercore'   	SETUP.lengthscale 	SETUP.lengthDim   	false         	false
            'Surf. temperature'     'T_{surf}'          'Tsurf'          	SETUP.Tscale        SETUP.TDim       	false          	false
            'CMB temperature'       'T_{CMB}'           'Tcmb'          	SETUP.Tscale        SETUP.TDim       	false        	false
            's-core'                'S_{core}'         	's_core'          	1                   ''                  false         	false
            'tc-core'               'Tc_{core}'       	'tc_core'        	1                   ''                  false           false
            'ts-core'               'Ts_{core}'        	'ts_core'        	1                   ''                  false         	false
            'Total erupta'          'Erupta_{Total}'    'erupta_total'    	SETUP.massscale  	SETUP.massDim    	false        	false
            'Outgassed water'       'Out_{water}'       'outgassed_water'  	1                   ''                  false         	IN.waterOut_logy
            'Total H_2O degassing'  'Out_{water,total}'	'outgassed_amount_H2O' 1              	'wt%'              	false         	IN.waterAmOut_logy %already wt%
            'Extra H_2O degassing'  'Out_{water,extra}' 'outgassed_extra_amountH2O' 100         'wt%'             	false          	IN.waterExOut_logy %water lost by melting
            'Psurf'                 'P_{surf}'          'Psurf'          	1                   ''                  false         	false
            'Weathering'            'We'                'weathering'     	1                   ''                  false        	false
            'H_2O total'            'm_{H2O}'           'mH2O_total'     	SETUP.massscale  	SETUP.massDim    	false         	false
            'H_2O mantle'           'm_{H2O,Mantle}'    'mH2O_mantle'     	SETUP.massscale  	SETUP.massDim      	false          	false
            'H_2O crust'            'm_{H2O,Crust}'     'mH2O_crust'    	SETUP.massscale  	SETUP.massDim     	false         	false
            'H_2O dm'               'dm_{H2O}'          'dmH2O'          	1                   ''                  false         	false
            'Outgassed carbon'      'Out_{carb}'        'outgassed_carbon' 	1                   ''                  false       	false
            'Outgassed nitrogen'    'Out_{nitro}'       'outgassed_nitrogen' 1                	''                  false         	false
            'H_2O dehydration'      'wflux_d'        	'wflux_d'         	100                 'wt%'             	false         	IN.waterDeg_logy
            'H_2O regassing'        'wflux_r'          	'wflux_r'         	100                 'wt%'            	false         	IN.waterReg_logy
            };
        
        if nColumns==17 %old stag version (<2013)
            FIELD_M = {
                'Time steps'        'cumul. time steps' 'istep'          	1                   '#'                 false           false
                'Time'              'time'              'time'           	PLOT.timeConvert    PLOT.time2plotDim   false           false
                'Top heatflux'     	'F_{top}'           'F_top'             SETUP.HFscale      	SETUP.HFDim       	false           false
                'Bottom heatflux'  	'F_{bot}'           'F_bot'             SETUP.HFscale    	SETUP.HFDim       	false          	false
                'Min. temperature'  'T_{min}'           'Tmin'              SETUP.Tscale        SETUP.TDim       	false        	false
                'Temperature'       'T_{mean}'          'Tmean'             SETUP.Tscale        SETUP.TDim       	true          	false
                'Max. temperature'  'T_{max}'           'Tmax'              SETUP.Tscale        SETUP.TDim       	false           false
                'Min. velocity'     'v_{min}'           'Vmin'              SETUP.Vscale        SETUP.vDim       	false           IN.v_logy
                'Velocity'          'v_{mean}'        	'Vrms'              SETUP.Vscale        SETUP.vDim       	true          	IN.v_logy
                'Max. velocity'     'v_{max}'           'Vmax'              SETUP.Vscale        SETUP.vDim       	false         	IN.v_logy
                'Min. viscosity'    '\eta_{min}'        'eta_min'         	SETUP.etascale    	SETUP.etaDim     	false          	false
                'Viscosity'         '\eta_{mean}'       'eta_mean'        	SETUP.etascale     	SETUP.etaDim     	true         	false
                'Max. viscosity'    '\eta_{max}'        'eta_max'         	SETUP.etascale    	SETUP.etaDim      	false       	false
                'Ra_{eff}'          'Ra_{eff}'          'ra_eff'           	1                   ''                  false           IN.Ra_logy
                'Nu_{top}'          'Nu_{top}'          'Nu_top'          	1                   ''                  false        	false
                'Nu_{bot}'          'Nu_{bot}'          'Nu_bot'           	1                   ''                  false        	false
                'Internal heating'  'H_{int}'           'H_int'          	SETUP.Hscale        SETUP.HDim       	false        	false
                };
        end
        
        % UNICODE VARIABLES
        %use the following character sequence to add unicode variables:
        %UCvarepsilondot => char(941)
        %UCvarepsilon => char(949)
        FIELD_M(:,2) = strrep(FIELD_M(:,2),'UCvarepsilondot',{char(941)});  %does only work if string consists ONLY of unicode character
        FIELD_M(:,2) = strrep(FIELD_M(:,2),'UCvarepsilon',{char(949)});     %does only work if string consists ONLY of unicode character
        
        % FIELD VARIABLES ADJUSTMENTS
        if SWITCH.DimensionalInput
            FIELD_M{strcmp(FIELD_M(:,3),'ra_eff'),4}    = SETUP.eta0.*SETUP.Ra0;    %undo Ra_eff = Ra_0/eta_mean if dimensional mode - set new Ra-scale
        end
        if SWITCH.DimensionalMode || SWITCH.DimensionalInput
            %nothing to do
        else
            %adjust for non-dimensional parameters
            FIELD_M(:,4) = {1};
            FIELD_M(:,5) = {'nd'};
        end
        
        %% DISPLAY MODEL DETAILS II
        disp(['Top-BC:            ',SETUP.topBC])
        disp(['Time Range:        ',num2str(minTimeDim,3),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(maxTimeDim,3),' ',PLOT.time2plotDim])
        
        %% PLOTTING
        RemoveXLabel = ones(nrPlots,1);
        numSubplot = 0;
        numTexting = 0;
        AlreadyWarnedForMissingTitleRow = false;
        for iPlot=1:nrPlots
            [ixPlot,izPlot] = ind2sub([nxSubplots,nzSubplots],iPlot);
            lastUsedPlot    = nrPlots;
            %extra axes for outside legend
            if SWITCH.LegendOutside
                lastUsedPlot    = nrPlots-1;
                if iPlot==nrPlots
                    
                    continue
                end
            end
            %set parameter to plot
            numField            = find(strcmp(FIELD_M(:,3),DATA.tag(iPlot)));    %find row in StagLab data matrix
            if DATA.titleRow  %check in output data directly
                iColumnInDatafile	= find(strcmp(DATA.title,DATA.tag(iPlot)));  %find column in data file
            else %guess position
                if ~AlreadyWarnedForMissingTitleRow
                    warning off backtrace
                    disp(' '); warning('No title row was found in the output data file: Check if the data matches the fields!'); disp(' ')
                    warning on backtrace
                    AlreadyWarnedForMissingTitleRow = true;
                end
                iColumnInDatafile	= numField;
            end
            
            %skip current parameter if not available in timedat file
            fieldFound      = true;
            if isempty(iColumnInDatafile) || iColumnInDatafile>nColumns
                fieldFound      = false;
                if SWITCH.Verbose
                    warning off backtrace
                    disp(' '); warning(['''',DATA.tag{iPlot},''' not available in current time.dat file!']);
                    warning on backtrace
                end
            end
            
            %set current variables
            if fieldFound
                FIELD.name          = FIELD_M{numField,1};
                FIELD.symbol        = FIELD_M{numField,2};
                FIELD.dim           = FIELD_M{numField,5};
                FIELD.logarithmic   = false; %dummy needed for f_DesignVaria
                FIELD.logX          = false;
                FIELD.logY          = FIELD_M{numField,7};
            else
                try
                    FIELD.name          = FIELD_M{numField,1};
                catch
                    FIELD.name          = DATA.tag(iPlot);
                end
                FIELD.symbol        = '';
                FIELD.dim           = '';
                FIELD.logarithmic   = false; %dummy needed for f_DesignVaria
                FIELD.logX          = false;
                FIELD.logY          = false;
            end
            
            %dimensionalise data
            if fieldFound
                ShowParameterRange          = false;
                if SWITCH.parameterRange && FIELD_M{numField,6} %parameter range
                    ShowParameterRange      = true;
                end
                %dimensionalise
                DATA.Array{1,iColumnInDatafile}         = DATA.Array{1,iColumnInDatafile}*FIELD_M{numField,4};
                if ShowParameterRange
                    DATA.Array{1,iColumnInDatafile-1}   = DATA.Array{1,iColumnInDatafile-1}*FIELD_M{numField-1,4};
                    DATA.Array{1,iColumnInDatafile+1}   = DATA.Array{1,iColumnInDatafile+1}*FIELD_M{numField+1,4};
                end
                %automatically set to either [mass fraction], [wt%] or [ppm]
                if strcmp(FIELD.dim,'wt%') && max(DATA.Array{1,iColumnInDatafile})<1.0 %change to [ppm]
                    conversionFactorPPM     = 1e4;
                    FIELD_M{numField,4}     = FIELD_M{numField,4}*conversionFactorPPM; %convert [wt%] to [ppm]
                    FIELD_M{numField,5}     = 'ppm';
                    FIELD.dim               = 'ppm';
                    %re-dimensionalise
                    DATA.Array{1,iColumnInDatafile}         = DATA.Array{1,iColumnInDatafile}*conversionFactorPPM;
                    if ShowParameterRange
                        FIELD_M{numField-1,4}     = FIELD_M{numField-1,4}*conversionFactorPPM; %convert [wt%] to [ppm]
                        FIELD_M{numField+1,4}     = FIELD_M{numField+1,4}*conversionFactorPPM; %convert [wt%] to [ppm]
                        FIELD_M{numField-1,5}     = 'ppm';
                        FIELD_M{numField+1,5}     = 'ppm';
                        DATA.Array{1,iColumnInDatafile-1}   = DATA.Array{1,iColumnInDatafile-1}*conversionFactorPPM;
                        DATA.Array{1,iColumnInDatafile+1}   = DATA.Array{1,iColumnInDatafile+1}*conversionFactorPPM;
                    end
                end
            end
            
            if ~SWITCH.QuickMode
                %set subplot
                if nrPlots>1
                    numSubplot	= numSubplot+1;
                    if SWITCH.UsePanel
                        hax     = SAVE.P(izPlot,ixPlot).select();
                        set(gcf,'CurrentAxes',hax)
                    else
                        hax  	= subplot(nzSubplots,nxSubplots,numSubplot);
                    end
                end
                if iPlot==lastUsedPlot || (PLOT.Performance && iPlot==lastUsedPlot-1) ||...
                        (lastUsedPlot>4 && iPlot==lastUsedPlot-1) || (lastUsedPlot>4 && PLOT.Performance && iPlot==lastUsedPlot-2)
                    RemoveXLabel(iPlot,1) = false;
                end
                
                if ~fieldFound
                    %insert empty plot
                    [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,[],FIELD,STYLE,SAVE);
                    
                    continue
                    
                end
                
                %performance plot
                if PLOT.Performance && iPlot==lastUsedPlot
                    dt = zeros(size(DATA.Array{1,timeRow},1),1);
                    for i=1:size(DATA.Array{1,timeRow},1)
                        if i==size(DATA.Array{1,timeRow},1)
                            dt(i,1) = DATA.Array{1,timeRow}(i,:)-DATA.Array{1,timeRow}(i-1,:);
                        else
                            dt(i,1) = DATA.Array{1,timeRow}(i+1,:)-DATA.Array{1,timeRow}(i,:);
                        end
                    end
                    pp1 = plot(DATA.Array{iColumnInDatafile}, dt,'Color',lineColor);
                    hold on
                    
                    DESVARIA.xlabelName     = FIELD.name;
                    DESVARIA.xlabelDim      = FIELD.dim;
                    DESVARIA.zlabelName     = 'Time Step';
                    DESVARIA.zlabelDim      = PLOT.time2plotDim;
                    DESVARIA.Task           = 'create annotation strings';
                    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,[],SWITCH,PLOT,STYLE,SAVE);
                    xlabel(PLOT.xlabel)
                    ylabel(PLOT.ylabel)
                    grid on
                    if STYLE.MinorAxisTicks; grid minor; end
                    
                    DESVARIA.Task    = 'setup axes';
                    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,[],SWITCH,PLOT,STYLE,SAVE);
                    
                    if SWITCH.MinMaxLegend; legend('^{dt}/_{time step}'); end
                    h_runs(irun,1) = pp1;
                    
                    break
                    
                end
                
                %time plot
                if ShowParameterRange
                    colourMax       = [0.64 0.08 0.18];
                    colourMin       = [0.0 0.45 0.74];
                    if PLOT.rangeAsArea
                        if FIELD.logY     %logarithmic y-axis
                            set(gca,'YScale','log')
                        end
                        %min-max area
                        edgecolor   = 'none';
                        areaAlpha   = 0.2;
                        hold on
                        maxData     = DATA.Array{iColumnInDatafile+1};
                        minData     = DATA.Array{iColumnInDatafile-1};
                        xData       = DATA.Array{1,timeRow};
                        pp4 = fill([xData' fliplr(xData')],[maxData' fliplr(minData')],...
                            lineColor,'EdgeColor',edgecolor,'EdgeAlpha',areaAlpha,'FaceAlpha',areaAlpha);
                        hold on
                    else
                        %min/max graphs
                        if FIELD.logY     %logarithmic y-axis
                            pp2 = semilogy(DATA.Array{1,timeRow}, DATA.Array{iColumnInDatafile-1},'Color',colourMin); 	%minimum
                            hold on
                            pp3 = semilogy(DATA.Array{1,timeRow}, DATA.Array{iColumnInDatafile+1},'Color',colourMax);	%maximum
                        else              %normal axis
                            pp2 = plot(DATA.Array{1,timeRow}, DATA.Array{iColumnInDatafile-1},'Color',colourMin); 	%minimum
                            hold on
                            pp3 = plot(DATA.Array{1,timeRow}, DATA.Array{iColumnInDatafile+1},'Color',colourMax);	%maximum
                        end
                    end
                    if FIELD.logY     %logarithmic y-axis
                        pp1 = semilogy(DATA.Array{1,timeRow}, DATA.Array{iColumnInDatafile},'Color',lineColor);   %mean
                    else              %normal axis
                        pp1 = plot(DATA.Array{1,timeRow}, DATA.Array{iColumnInDatafile},'Color',lineColor);   %mean
                    end
                    hold on
                    
                    if SWITCH.MinMaxLegend
                        if PLOT.rangeAsArea
                            legend([pp1,pp4],[FIELD_M(numField,2),strcat(FIELD_M(numField-1,2),'-',FIELD_M(numField+1,2))])
                        else
                            legend([pp3,pp1,pp2],[FIELD_M(numField+1,2),FIELD_M(numField,2),FIELD_M(numField-1,2)],'Location','Best');
                        end
                    end
                else
                    if FIELD.logY     %logarithmic y-axis
                        pp1 = semilogy(DATA.Array{1,timeRow}, DATA.Array{iColumnInDatafile},'Color',lineColor);   %mean
                        hold on
                    else              %normal axis
                        pp1 = plot(DATA.Array{1,timeRow}, DATA.Array{iColumnInDatafile},'Color',lineColor);   %mean
                        hold on
                    end
                    if SWITCH.MinMaxLegend
                        legend(pp1,FIELD_M(numField,2),'Location','Best');
                    end
                end
                h_runs(irun,1) = pp1;
                
                %min/mean/max values
                meanVal     = nanmean(DATA.Array{iColumnInDatafile});
                if ShowParameterRange
                    minVal  = min(DATA.Array{iColumnInDatafile-1}); %minimum
                    maxVal  = max(DATA.Array{iColumnInDatafile+1}); %maximum
                else
                    minVal  = min(DATA.Array{iColumnInDatafile});   %min-mean value
                    maxVal  = max(DATA.Array{iColumnInDatafile});   %max-mean value
                end
                
                %limit time-axis and parameter-axis extent
                axis tight
                axisLimitX          = xlim; %current x-axis limits
                axisLimitY          = ylim; %current y-axis limits
                if SWITCH.AxesLimit
                    axisLimit       = axisLimitX; %current x-axis limits
                    if ~isnan(IN.XAxisMinValue)
                        axisLimit(1,1) 	= IN.XAxisMinValue;
                    end
                    if ~isnan(IN.XAxisMaxValue)
                        axisLimit(1,2) 	= IN.XAxisMaxValue;
                    end
                    DESVARIA.axesLimits     = [axisLimit(1,1),axisLimit(1,2),axisLimitY];
                    
                    axisLimit       = axisLimitY; %current y-axis limits
                    currentIdx      = min(length(IN.YAxisMinValue),iPlot);
                    if ~isnan(IN.YAxisMinValue(currentIdx))
                        axisLimit(1,1) 	= IN.YAxisMinValue(currentIdx);
                    end
                    currentIdx      = min(length(IN.YAxisMaxValue),iPlot);
                    if ~isnan(IN.YAxisMaxValue(currentIdx))
                        axisLimit(1,2) 	= IN.YAxisMaxValue(currentIdx);
                    end
                    DESVARIA.axesLimits     = [axisLimitX,axisLimit(1,1),axisLimit(1,2)];
                end
                
                %automatic y-axis extent for constant data
                if minVal==maxVal && maxVal>-1 && maxVal<1 %constant array
                    if maxVal<0
                        DESVARIA.axesLimits     = [axisLimitX,2*maxVal,0];
                    elseif maxVal>0
                        DESVARIA.axesLimits     = [axisLimitX,0,2*maxVal];
                    else %maxVal==0
                        %nothing to do
                    end
                    if SWITCH.AxesLimit && SWITCH.Verbose
                        warning('Manual axis limits might have been overwritten due to constant-parameter adjustment.');
                    end
                end
                
                %axes labelling
                dimStringTime           = FIELD_M{timeRow,5};
                DESVARIA.xlabelName     = FIELD_M{timeRow,1};
                DESVARIA.xlabelDim      = dimStringTime;
                DESVARIA.zlabelName     = FIELD.name;
                DESVARIA.zlabelDim      = FIELD.dim;
                DESVARIA.Task           = 'create annotation strings';
                [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,[],SWITCH,PLOT,STYLE,SAVE);
                if SWITCH.SimplifyPlots && RemoveXLabel(iPlot,1)
                    %no xlabel
                else
                    xlabel(PLOT.xlabel);
                end
                ylabel(PLOT.ylabel)
                grid on
                if STYLE.MinorAxisTicks; grid minor; end
                
                DESVARIA.Task    = 'setup axes';
                [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,[],SWITCH,PLOT,STYLE,SAVE);
                
                %texting
                if SWITCH.Texting
                    text(max(xlim),max(ylim),[num2str(minVal,3),' ',STYLE.SCHAR.leftArrow,' ',num2str(meanVal,3),' ',STYLE.SCHAR.rightArrow,' ',num2str(maxVal,3),' ',FIELD.dim],...
                        'HorizontalAlignment','right','VerticalAlignment','bottom','Color',STYLE.keyColor,'FontSize',8)
                end
                numTexting  = numTexting+1;
                if numTexting==1
                    disp(['Parameter Range:   ',num2str(minVal),' ',STYLE.SCHAR.leftArrow,' ',num2str(meanVal,3),' ',STYLE.SCHAR.rightArrow,' ',num2str(maxVal),' ',FIELD.dim])
                else
                    disp(['                   ',num2str(minVal),' ',STYLE.SCHAR.leftArrow,' ',num2str(meanVal,3),' ',STYLE.SCHAR.rightArrow,' ',num2str(maxVal),' ',FIELD.dim])
                end
            end
            
            %remember graph data to be saved
            if SAVE.Graphs && fieldFound
                %create array to save
                if ~exist('titleArray','var'); titleArray = {[FIELD_M{timeRow,1},'[',FIELD_M{timeRow,5},']']}; dataArray = DATA.Array{1,timeRow}; end
                titleArray              = [titleArray, {[FIELD.name,'[',FIELD.dim,']']}];
                dataArray               = [dataArray, DATA.Array{iColumnInDatafile}];
            end
        end %loop
        disp(' ')
    end
    
    %% SAVING GRAPH DATA
    if SAVE.Graphs
        SAVE.Directory     	= FILE.directory;
        SAVE.DataName      	= [FILE.name,'_TimeGraphs'];
        SAVE.DataTitle    	= strjoin(titleArray,', ');
        SAVE.data          	= dataArray;
        [SAVE.overwriteAll] = f_saveData( SAVE );
        SAVE                = rmfield(SAVE,{'data','DataTitle'});
    end
    
    %% LEGEND
    if ~SWITCH.MinMaxLegend && ~SWITCH.QuickMode %not applied if individual legends are switched on
        if isfield(PLOT,'titleString')
            if length(h_runs)>length(PLOT.titleString); PLOT.titleString{length(PLOT.titleString)+1:length(h_runs)} = '[undefined]'; end
            legendANN = PLOT.titleString; 
        else
            legendANN = IN.Name; 
        end
        hl = legend(h_runs,legendANN,'Location','NorthEast','interpreter','none');
        % hl = legend(h_runs,IN.Name,'Location','SouthEast');
        if SWITCH.LegendOutside
            numSubplot          = numSubplot+1;
            if SWITCH.UsePanel
                [ixPlot,izPlot] = ind2sub([nxSubplots,nzSubplots],numSubplot);
                hax             = SAVE.P(izPlot,ixPlot).select();
            else
                hax             = subplot(nzSubplots,nxSubplots,numSubplot);
            end
            EmptyLegendSpacePosition = get(hax,'Position');
            axis off
            LegendPosition      = get(hl,'Position');
            LegendPosition(1)   = EmptyLegendSpacePosition(1);
            LegendPosition(2)   = EmptyLegendSpacePosition(2)+EmptyLegendSpacePosition(4)-LegendPosition(4);
            set(hl,'Position',LegendPosition);
            if SWITCH.UsePanel
                [ixPlot,izPlot] = ind2sub([nxSubplots,nzSubplots],numSubplot-1);
                hax             = SAVE.P(izPlot,ixPlot).select();
            else
                hax             = subplot(nzSubplots,nxSubplots,numSubplot-1); %back to the last previous subplot
            end
        end
    end
    
    %% STYLE PLOT
    if SWITCH.PlotDesign
        %     STYLE.Legend = false; %done later
        if SWITCH.LegendDesign; STYLE.Legend = true; else; STYLE.Legend = false; end
        f_DesignFigure(STYLE);
    end
    if SWITCH.BackgroundDesign
        BACK.BulletString  	= STYLE.SCHAR.hugeBulletDark;
        BACK.ApplyTo        = 'all';          %'current', 'all'
        BACK.ColorMode      = STYLE.ColorMode;
        f_DesignBackground(BACK) %set plot background of each subplot
    end
    
    %% ANNOTATION
    if SWITCH.Annotation
        %     pause(0.5);
        f_DesignAnnotations(PLOT,STYLE);
    end
    
    %% fAIO ACTION: FINISHING UP
    if SWITCH.LegendOutside; axes(hax); end %activate second last axes
    fAIo.task = 'finishingUp';
    [fAIo,SAVE] = f_AIo(fAIo,SWITCH,FILE,[],STYLE,SAVE);
    
    %% SAVING PLOTS
    if SAVE.Figure
        % SAVING DIRECTORY
        if strcmp(SAVE.writeDirectory,'auto') %save to standard folder (e.g., .../+im/...)
            SAVE.Directory          = FILE.stemSave;
        else
            SAVE.Directory          = SAVE.writeDirectory;
        end
        if strcmp(STYLE.ColorMode,'light')
            SAVE.blackBackground    = false;
        elseif strcmp(STYLE.ColorMode,'dark')
            SAVE.blackBackground    = true;
        end
        endix                       = [strjoin(IN.Name,'_'),'_',strjoin(DATA.tag,'')];
        if SWITCH.AnalysisMode;             endix = [endix,'An'];   	end
        if SWITCH.DimensionalMode;          endix = [endix,'Dim'];   	end
        if strcmp(STYLE.ColorMode,'dark');  endix = [endix,'Dark'];   	end
        SAVE.FigureName             = strcat(IN.saveString,'_',endix);   %file name
        SAVE.FigureNr               = 1;
        f_saveFigure( SAVE )
    end

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
        disp(' '); beep;
        return;
    end
end
 
%% fAIO ACTION: END MESSAGE
fAIo.task = 'EndMessage';
[fAIo,SAVE] = f_AIo(fAIo,SWITCH,FILE,[],STYLE,SAVE);

end



