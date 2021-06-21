
%%                                                   SL_RadialProfile 8.000
%
%                                                 plots StagYY's rprof.dat
%    . calls f_DesignBackground
%    . calls f_Dimensions
%    . calls f_DesignFigure
%    . calls f_DesignLegend
%    . calls f_FileFinder
%    . calls f_saveData
%    . calls f_saveFigure
%    . calls f_SetupFigure
%    . calls f_Startup
%    . calls f_DesignVaria
%                                               17.06.2021 : Fabio Crameri

%% CONTRIBUTIONS
% Kiran Chotalia
%   29.06.2017, Option to save processed profiles
%   29.06.2017, More informative file naming

%% UPDATES >8.0

%% UPDATES 8.0
% *option to plot refstat data (combined adiabat)*
% *adjustments to extended rprof.dat output of latest StagYY version*
% *preventing annotating separate legend panel*
% *improved line colouring*
% *bug fix (annotation mismatch when using legend-outside option)*
% *bug fix (lines plotted on top of each other without min/max got erased)*
% *bug fix (frequent error due to 'if'-statement mismatching string length)*
% *error prevention (refstat.dat could not be found)
% *error prevention (refstat.dat finding ADIABAT line in data file)
% *prevent error with unequal number of IN.Numbers compared to IN.Name*

%% UPDATES 7.0
% *introducing panel*
% *introducing automated error logging*
% *introducing categorical Scientific Colour Maps*
% *introducing journal-specific plot design*
% *special character fix for Windows*

%% UPDATES 6.0
% *introducing analysis mode*
% *option for displaying min/max values*
% *option to place legend separately*
% *more robust file read*
% *automatic conversion of mass fractions to wt% or ppm*
% *manual y-axis limits*
% *additional field parameters*
% *additional error prevention*
% *bug fixes*

%% UPDATES 5.0
% *deeper fAIo integration*
% *using improved fileFinder*
% *more flexible file read*
% *more informative file naming*
% *option to name files manually*
% *option to plot solidus*
% *option to save processed profile data*
% *option to set default screen position of figure*
% *option for quick mode*
% *speed improvements*
% *bug fixes*

%% UPDATES 4.0
% *more stable design-routines*
% *cleaner plot design and layout*
% *automatic figure size*

%% UPDATES 3.0
% *improved plot layout*
% *improved code design*
% *option added for transparent legend background*
% *default variables added*

function [PLOT] = SL_RadialProfile(IN,PLOT,SWITCH,STYLE,SAVE)
if exist('IN','var') %if run from a parfile
    %% DEFAULT VARIABLES
    %...see f_DefaultsRprof.m
else
    clear;
    %% INPUT (...use parameter file instead)
    SAVE.StartNumber = 1; SAVE.EndNumber = 1;
    [IN,PLOT,SWITCH,STYLE,SAVE] = f_DefaultsRprof; %get default variables
    IN.Name         =   {  'test'  'test' };	% filename
    IN.Number       =   [   2           3 	];
    IN.Parameter    =   [   11      ];
    IN.Folder       =   {   '~/work/plotTest/'   };
end

%% STARTUP PROCEDURE
SAVE.app = 'SL_RadialProfile'; SAVE.appVersion = 8.00; %only 3 digits to prevent unecessary parfile updates
if ~isfield(SAVE,'count'); SAVE.count = 0; end
DESVARIA.Task    = 'set special characters';
[~,DESVARIA,STYLE,~] = f_DesignVaria(DESVARIA,[],[],[],[],STYLE,[]);
[fAIo,SAVE,SWITCH] = f_Startup(STYLE,SAVE,SWITCH);
if ~strcmp(fAIo.Status,'all fine!'); disp(' '); disp(['                         ',STYLE.SCHAR.checkMark,' Finished.']); return; end

try
    %% GET NUMBER OF SUBPLOTS
    nrPlots=0; readRefstat=0;
    if PLOT.Temperature;            nrPlots = nrPlots+1; parameterX(nrPlots)=2;     varString{nrPlots}='T';         fieldString{nrPlots}='Temperature'; end
    if PLOT.Heatflux;               nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;   varString{nrPlots}='hf';        fieldString{nrPlots}='Heatflux'; end
    if PLOT.Viscosity;              nrPlots = nrPlots+1; parameterX(nrPlots)=14;    varString{nrPlots}='eta';       fieldString{nrPlots}='Viscosity'; end
    if PLOT.Stress;                 nrPlots = nrPlots+1; parameterX(nrPlots)=20;    varString{nrPlots}='str';       fieldString{nrPlots}='Stress'; end
    if PLOT.StrainRate;             nrPlots = nrPlots+1; parameterX(nrPlots)=17;    varString{nrPlots}='edot';      fieldString{nrPlots}='Strain rate'; end
    if PLOT.Crust;              	nrPlots = nrPlots+1; parameterX(nrPlots)=37;    varString{nrPlots}='c';         fieldString{nrPlots}='Crust'; end
    if PLOT.Air;                    nrPlots = nrPlots+1; parameterX(nrPlots)=43;    varString{nrPlots}='air';    	fieldString{nrPlots}='Air'; end
    if PLOT.Primordial;             nrPlots = nrPlots+1; parameterX(nrPlots)=46;    varString{nrPlots}='prim';    	fieldString{nrPlots}='Primordial'; end
    if PLOT.ContinentalCrust;   	nrPlots = nrPlots+1; parameterX(nrPlots)=49;    varString{nrPlots}='cc';        fieldString{nrPlots}='Cont. crust'; end
    if PLOT.Fluid;                  nrPlots = nrPlots+1; parameterX(nrPlots)=52;    varString{nrPlots}='fl';        fieldString{nrPlots}='Fluid'; end
    if PLOT.Metal;                  nrPlots = nrPlots+1; parameterX(nrPlots)=55;    varString{nrPlots}='m';         fieldString{nrPlots}='Metal'; end
    if PLOT.Water;              	nrPlots = nrPlots+1; parameterX(nrPlots)=58;    varString{nrPlots}='water';   	fieldString{nrPlots}='Water'; end
    if PLOT.Velocity;               nrPlots = nrPlots+1; parameterX(nrPlots)=5;     varString{nrPlots}='v';         fieldString{nrPlots}='Velocity'; end
    if PLOT.VerticalVelocity;       nrPlots = nrPlots+1; parameterX(nrPlots)=8;     varString{nrPlots}='vz';        fieldString{nrPlots}='Vertical velocity'; end
    if PLOT.HorizontalVelocity;     nrPlots = nrPlots+1; parameterX(nrPlots)=11;    varString{nrPlots}='vh';        fieldString{nrPlots}='Horizontal velocity'; end
    if PLOT.VertVorticity;      	nrPlots = nrPlots+1; parameterX(nrPlots)=23;	varString{nrPlots}='wv';        fieldString{nrPlots}='Vert. vorticity'; end
    if PLOT.HorizVorticity;     	nrPlots = nrPlots+1; parameterX(nrPlots)=26;	varString{nrPlots}='wh';        fieldString{nrPlots}='Horiz. vorticity'; end
    if PLOT.Divergence;          	nrPlots = nrPlots+1; parameterX(nrPlots)=29;	varString{nrPlots}='div';    	fieldString{nrPlots}='Divergence'; end
    if PLOT.Density;                nrPlots = nrPlots+1; parameterX(nrPlots)=40;    varString{nrPlots}='rho';       fieldString{nrPlots}='Density'; end
    if PLOT.Advection;              nrPlots = nrPlots+1; parameterX(nrPlots)=32;    varString{nrPlots}='adv';       fieldString{nrPlots}='Thermal advection'; end
    if PLOT.Diffusion;              nrPlots = nrPlots+1; parameterX(nrPlots)=33;    varString{nrPlots}='diff';      fieldString{nrPlots}='Thermal diffusion'; end
    if PLOT.InternalHeating;        nrPlots = nrPlots+1; parameterX(nrPlots)=34;    varString{nrPlots}='rh';        fieldString{nrPlots}='Internal heating'; end
    if PLOT.ViscousDissipation;  	nrPlots = nrPlots+1; parameterX(nrPlots)=35;    varString{nrPlots}='diss';      fieldString{nrPlots}='Visc. dissipation'; end
    if PLOT.AdiabaticHeating;   	nrPlots = nrPlots+1; parameterX(nrPlots)=36;    varString{nrPlots}='ah';        fieldString{nrPlots}='Adiabatic heating'; end

    if PLOT.ViscousDissipationLog;	nrPlots = nrPlots+1; parameterX(nrPlots)=61;    varString{nrPlots}='vdiss';   	fieldString{nrPlots}='Visc. Dissipation Log'; end
    if PLOT.AdvectionTotal;       	nrPlots = nrPlots+1; parameterX(nrPlots)=64;    varString{nrPlots}='advT';   	fieldString{nrPlots}='Advection Total'; end
    if PLOT.ThermalConduction;    	nrPlots = nrPlots+1; parameterX(nrPlots)=67;    varString{nrPlots}='cond';   	fieldString{nrPlots}='Thermal conduction'; end
    
    if PLOT.Toroidal;               nrPlots = nrPlots+1; parameterX(nrPlots)=70;    varString{nrPlots}='to';        fieldString{nrPlots}='Toroidal'; end
    if PLOT.ToroidalRMS;            nrPlots = nrPlots+1; parameterX(nrPlots)=71;    varString{nrPlots}='rorms';  	fieldString{nrPlots}='Toroidal_rms'; end
    if PLOT.Poloidal;               nrPlots = nrPlots+1; parameterX(nrPlots)=72;    varString{nrPlots}='po';        fieldString{nrPlots}='Poloidal'; end
    if PLOT.PoloidalRMS;            nrPlots = nrPlots+1; parameterX(nrPlots)=73;    varString{nrPlots}='porms';   	fieldString{nrPlots}='Poloidal_rms'; end
    if PLOT.Mobility;               nrPlots = nrPlots+1; parameterX(nrPlots)=74;    varString{nrPlots}='M';         fieldString{nrPlots}='Mobility'; end
    if PLOT.f80;                    nrPlots = nrPlots+1; parameterX(nrPlots)=75;    varString{nrPlots}='f80';       fieldString{nrPlots}='f80'; end
    if PLOT.f90;                    nrPlots = nrPlots+1; parameterX(nrPlots)=76;    varString{nrPlots}='f90';       fieldString{nrPlots}='f90'; end
    %REFERENCE STATE DATA
%     if PLOT.RefstatTemperature;     readRefstat = 1; nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;    varString{nrPlots}='Tref';        fieldString{nrPlots}='RS Temperature'; end
%     if PLOT.RefstatDensity;         readRefstat = 1; nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;    varString{nrPlots}='rho';         fieldString{nrPlots}='RS Density'; end
%     if PLOT.RefstatExpansivity;       readRefstat = 1; nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;    varString{nrPlots}='expan';       fieldString{nrPlots}='RS Expansivity'; end
%     if PLOT.RefstatConductivity;    readRefstat = 1; nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;    varString{nrPlots}='Cp';          fieldString{nrPlots}='RS Heat capacity'; end
%     if PLOT.RefstatPressure;        readRefstat = 1; nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;    varString{nrPlots}='Tkappa';      fieldString{nrPlots}='RS Diffusivity'; end
%     if PLOT.RefstatPressure;        readRefstat = 1; nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;    varString{nrPlots}='Tcond';       fieldString{nrPlots}='RS Conductivity'; end
    if PLOT.RefstatCaTemperature; 	readRefstat = 1; nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;    varString{nrPlots}='CA_Tref';  	fieldString{nrPlots}='RefS Temperature'; end
    if PLOT.RefstatCaDensity;      	readRefstat = 1; nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;    varString{nrPlots}='CA_rho';  	fieldString{nrPlots}='RefS Density'; end
    if PLOT.RefstatCaExpansivity; 	readRefstat = 1; nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;    varString{nrPlots}='CA_expan';	fieldString{nrPlots}='RefS Expansivity'; end
    if PLOT.RefstatCaConductivity;  readRefstat = 1; nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;    varString{nrPlots}='CA_Tcond'; 	fieldString{nrPlots}='RefS Conductivity'; end
    if PLOT.RefstatCaPressure;    	readRefstat = 1; nrPlots = nrPlots+1; parameterX(nrPlots)=NaN;    varString{nrPlots}='CA_P';        fieldString{nrPlots}='RefS Pressure'; end

    %extra axes for outside legend
    if SWITCH.LegendOutside
        nrPlots = nrPlots+1;
        PLOT.annotationString{nrPlots} = ''; %no annotation for extra legend panel (which is called second last)
    end
    
    %% LAYOUT SUBPLOTS
    if nrPlots<=4
        nxSubplots  = nrPlots;
        nzSubplots  = 1;
    else
        nxSubplots  = floor((nrPlots+1)/2);
        nzSubplots  = min(nrPlots,2);
    end
    
    %% VARIABLE ADJUSTMENTS
    nrModels = size(IN.Name,2);
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
        SAVE.Figure                	= false;
        SWITCH.PlotDesign         	= false;
        SWITCH.BackgroundDesign   	= false;
        SWITCH.LegendDesign        	= false;
        SWITCH.SimplifyPlots    	= false;
        SWITCH.Annotation         	= false;
        IN.showLegend             	= false;
        SWITCH.Title              	= false;
        SWITCH.d_eta             	= false;
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
    if size(IN.Parameter,2)<nrModels
        dummy = ones(1,nrModels-size(IN.Parameter,2)); dummy(:) = IN.Parameter(1,1);
        IN.Parameter = [IN.Parameter, dummy];
    end
    if size(IN.Number,2)<nrModels
        dummy = ones(1,nrModels-size(IN.Number,2)); dummy(:) = IN.Number(1,1);
        IN.Number = [IN.Number, dummy];
        warning off backtrace; warning('IN.Number does not have the same number of entries as IN.Name.'); warning on backtrace
    end
    % check for length of string input arrays
    if size(IN.Folder,2)<nrModels
        dummy = cell(1,nrModels-size(IN.Folder,2)); dummy(:) = IN.Folder(1,1);
        IN.Folder = [IN.Folder, dummy];
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
    
    %% FIGURE SETUP
    if ~SWITCH.QuickMode
        clf(figure(1));
        if strcmp(STYLE.ColorMode,'light')
            STYLE.ColorModeBW   = 'white';
        elseif strcmp(STYLE.ColorMode,'dark')
            STYLE.ColorModeBW   = 'black';
        end
        %Setup colour map
        CMAPPING.Task   = 'Categorical';
        CMAPPING.NumberColours = size(IN.Name,2);
        [SWITCH,CMAPPING] = f_DesignColourmap(SWITCH,PLOT,[],STYLE,CMAPPING);
        PLOT.ColorVector    = CMAPPING.ColourVector;
        %Further figure setup
        set(gcf,'color',STYLE.ColorModeBW);
        % set figure position
        if strcmp(SWITCH.FigurePosition,'auto') %set figure position automatically
            FPOS.Default            = PLOT.FigureDefaultPosition;
            FPOS.subplotLayout      = [nzSubplots nxSubplots];
            [FPOS] = f_SetupFigure(FPOS,SAVE,1,1);
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
    end
    
    %% LOOP FILES
    yMinOld = 1;
    hRuns = zeros(size(IN.Name,2),1); hRunsArea = hRuns; fileNameString = '';
    for ifile=1:size(IN.Name,2)
        clearvars -except IN SWITCH PLOT STYLE SAVE ifile hRuns hRunsArea yMinOld ...
            nrPlots parameterX varString fieldString nxSubplots nzSubplots fileNameString ...
            readRefstat
        FILE.name           = IN.Name{ifile};
        FILE.number         = IN.Number(ifile);
        FILE.stemRead     	= IN.Folder{ifile};
        FILE.stemSave     	= IN.Folder{ifile};
        SWITCH.parameter    = IN.Parameter(ifile);
        lineColor           = PLOT.ColorVector(ifile,:);
        
        %% CHECK FOR INPUT FILE
        [FILE] = f_FileFinder(FILE,SWITCH,SAVE,STYLE);
        if FILE.NotFound
            SAVE.LastFile           = true;
            close all
            disp([FILE.NotFoundError])
            return
            
        end
        
        %% RPROF.DAT FILE STRUCTURE
        % OLD STAGYY VERSIONS
        %  1  2 3  4   5  6  7   8    9  10   11  12  13   14   15   16    17 18 19  20 21 22
        %  z,t0,t1,t2, v0,v1,v2, vz0,vz1,vz2, vh0,vh1,vh2, eta0,eta1,eta2, e0,e1,e2, s0,s1,s2,
        
        %  23 24 25 26 27 28  29 30 31  32   33  34    35        36          37 38 39 ...
        %  w0,w1,w2,w3,w4,w5, d0,d1,d2, adv,diff,Rh,visc.diss,adiab.heating, c0,c1,c2 ...
        
        % NEW STAGYY VERSIONS
        % 1 2     3    4    5    6    7    8     9     10    11    12    13    14     15     16     17   18   19   20   21   22
        % r Tmean Tmin Tmax vrms vmin vmax vzabs vzmin vzmax vhrms vhmin vhmax etalog etamin etamax elog emin emax slog smin smax 
        % 23    24    25    26    27    28    29   30   31   32    33     34     35         36       37    38   39 ...
        % whrms whmin whmax wzrms wzmin wzmax drms dmin dmax enadv endiff enradh enviscdiss enadiabh cmean cmin cmax 
        % 40      41     42     43      44     45     46       47      48      49     50    51    52        53       54
        % rhomean rhomin rhomax airmean airmin airmax primmean primmin primmax ccmean ccmin ccmax fmeltmean fmeltmin fmeltmax 
        % 55        56       57       58     59    60    61          62          63          64     65      66     67        68       69       
        % metalmean metalmin metalmax gsmean gsmin gsmax viscdisslog viscdissmin viscdissmax advtot advdesc advasc tcondmean tcondmin tcondmax 
        
        %% FILE DIAGNOSTICS
        %read data
        if SWITCH.closeOldWaitbars; wbOld = findall(0,'tag','TMWWaitbar'); delete(wbOld); end %close open ghost waitbars
        wb = waitbar(0,'Please wait...');
        %find #variable-entries in one line
        fid                 = fopen(FILE.found,'r');
        numEntries=0; iLine=0;
        while numEntries<15 %line has less than 15 characters (i.e., is empty)
            iLine           = iLine+1;
            numEntries      = numel(regexp(fgetl(fid),'\s*([^\s]*)\s*'));
        end
        frewind(fid);
        waitbar(1/6,wb)
        %find #depth-entries
        numberLines         = 0;
        cc                  = 0;
        tline               = fgetl(fid);
        while cc<4 && ~strcmp(num2str(tline),'-1')  %tline becomes -1 when file ended during the loop
            tline           = fgetl(fid);
            numberLines     = numberLines+1;
            if ~isempty(strfind(tline,'********'))
                %if contains(tline,'********') %not compatible with MatLab versions older than 2016
                cc  = cc+1;
                lineNumNewSet(cc) = numberLines;
            end
            cc_old          = cc;
        end
        if numel(lineNumNewSet)==1
            if SWITCH.Verbose; warning off backtrace; warning('The rprof.dat file includes only 1 time step!'); warning on backtrace;  end
            nz = numberLines -2;
        elseif numel(lineNumNewSet)==2
            nz = lineNumNewSet(2)-lineNumNewSet(1)  -1;  %num data in z-direction
        else
            nz = lineNumNewSet(3)-lineNumNewSet(2)  -1;  %num data in z-direction
        end
        waitbar(2/6,wb)
        %find #empty lines at the beginning
        numEmptyLinesBeginning  = lineNumNewSet(1);
        %find total #lines
        frewind(fid);
        try
            if ispc; error(' '); end %this is not available yet on pc systems
            option                  = 2;
            if option==1  %only for unix systems (slower - empty lines are not counted)
                [~,numLinesStr]     = system(['grep -c ".$" ',FILE.found]);
                numberLines        	= str2double(numLinesStr);
                
            elseif option==2  %only for unix systems (fast - empty lines are counted too)
                [~,numLinesStr]     = system(['wc -l ',FILE.found]);
                for iNum=1:length(numLinesStr)
                    if ~isnan(str2double(numLinesStr(iNum)))
                        firstNumber     = iNum;
                        lastNumber      = iNum;
                        iNum         	= iNum+1;
                        while ~isnan(str2double(numLinesStr(iNum)))
                            lastNumber  = iNum;
                            iNum      	= iNum+1;
                        end
                        break
                    else
                        firstNumber = [];
                    end
                end
                numberLines     	= str2double(numLinesStr(firstNumber:lastNumber));
                numberLines         = numberLines-numEmptyLinesBeginning; %remove empty first line
                %numberLines      	= numberLines-1;  %remove empty last line
            end
        catch
            if SWITCH.Verbose; disp('> fast method is not available - please be patient...'); end
            option                  = 1;
            if option==1
                numberLines         = 0;
                tline               = fgetl(fid);
                while ischar(tline)
                    tline           = fgetl(fid);
                    numberLines     = numberLines+1;
                end
                numberLines         = numberLines-numEmptyLinesBeginning;  %remove empty first line
                %numberLines         = numberLines-1;  %remove empty last line
            elseif option==2
                numberLines         = numel(textread(FILE.found,'%1c%*[^\n]'));
            end
        end
        waitbar(3/6,wb)
        %find #time-snapshots
        numberTimeSnapshots         = numberLines/(nz+1);
        if ~isreal(numberTimeSnapshots) || mod(numberTimeSnapshots,1)~=0  %if not an integer
            numberTimeSnapshots     = round(real(numberTimeSnapshots));
            if isreal(numberLines); numberLines = round(real(numberLines)); end
            if isreal(nz); nz = round(real(nz)); end
            warning('numberTimeSnapshots is not an integer: check the rprof.dat-file variation!')
        end
        %define plotted filenumber
        if FILE.number>numberTimeSnapshots-1
            FILE.number     = numberTimeSnapshots-1;    %max. file number
        else
            FILE.number     = FILE.number;              %chosen file number
        end
        
        %% DISPLAY FILE DETAILS
        disp(['Number:            ',num2str(FILE.number)])
        if SWITCH.Verbose
            % disp(['datafile is ',num2str(n),' x ',num2str(nLines)])
            disp(['Data Structure:       __',num2str(numEntries),'__>'])
            disp('                     |')
            disp(['                   ',num2str(nz)])
            disp(['                     |       for ',num2str(numberTimeSnapshots),' time steps'])
            disp('                     V')
        end
        
        %% READING DATA
        readInOption1failed = false;
        stillReading = true; iRead = 0;
        while stillReading
            iRead = iRead+1;
            frewind(fid);
            %lines to skip
            numSkipLines    = (FILE.number)*(nz+1) +numEmptyLinesBeginning;
            %jump to given line number
            line            = textscan(fid,'%s',1,'delimiter','\n','headerlines',numSkipLines-1);
            %read header
            line            = fscanf(fid,'%s',[1 1]);       %asteriks
            i_time_step     = fscanf(fid,'%i',[1 1]);       %time step
            line            = fscanf(fid,'%s %s %s',[1 3]); %string
            time            = fscanf(fid,'%f',[1 1]);       %time
            %read data
            if iRead==1; readInOption2 = SWITCH.readOption2; end
            mData = zeros(nz,numEntries)*NaN;
            if readInOption2 || readInOption1failed %this version gets messed up if some values are bad
                for j=1:nz
                    waitbar(3/6+(j/nz*3/6),wb)
                    line        = fscanf(fid,'%e',[1 numEntries]);
                    mData(j,:)  = line;
                end
                stillReading = false;
                
            elseif ~readInOption2 %replaces bad values with some other value
                emptyValue = NaN;
                formatSpec = [];
                for iEntry=1:numEntries
                    formatSpec          = [formatSpec,'%f'];
                end
                WarningDone = false;
                stillReading = false;
                for j=1:nz
                    waitbar(3/6+(j/nz*3/6),wb)
                    %line        = textscan(fid,formatSpec,[1 numEntries],'Delimiter','','WhiteSpace','','EmptyValue',emptyValue,'EndOfLine','\r\n');
                    line        = textscan(fid,formatSpec,[1 numEntries],'Delimiter','','EmptyValue',emptyValue,'TreatAsEmpty',{'Infinity','-Infinity'},'EndOfLine','\r\n');
                    try
                        mData(j,:)  = cell2mat(line);
                    catch
                        %try to fix
                        for iCheck=1:length(line)
                            if ~isreal(cell2mat(line(iCheck))) %Check for infinity values
                                line{iCheck} = 0;
                                warning off backtrace
                                warning(['Found a non-real (i.e., imaginary) number in column ',num2str(iCheck),' of the rprof.dat file. Check your StagYY run carefully!'])
                                warning on backtrace
%                                 %Switch to readInOption2
%                                 readInOption1failed = true;
%                                 stillReading = true;
%                                 break
                                
                            end
                            if isempty(cell2mat(line(iCheck))) %Check for empty entries
                                line{iCheck} = 0;
                                if SWITCH.Verbose && ~WarningDone
                                    warning off backtrace; warning(['Found an empty entry in column ',num2str(iCheck),' (and possibly other ones) of the rprof.dat file.']); warning on backtrace;
                                    WarningDone = true;
                                end
                            end
                        end
                        if readInOption1failed; break; end
                        mData(j,:)  = cell2mat(line);
                    end
                end
            end
        end
        close(wb)
        fclose(fid);
        
        % Read reference state profiles
        if readRefstat
            try
                RS.nz           = nz;
                RS.Task     	= 'ImportRefstat';
                [RS,~] = f_Import(RS,FILE,[],[],[]);
            catch
                warning('Refstat data could not be read. Check here.');
                readRefstat     = false;
            end
        else
            RS = []; %dummy variable
        end
        
        %% DIAGNOSE INPUT FILE
        if max(mData(:,1))>0.5 && max(mData(:,1))<1.5  %if 0.5<D<1.5 then non-dimensional
            SWITCH.DimensionalInput = false;    dimString = 'non-dimensional';
        else
            SWITCH.DimensionalInput = true;     dimString = 'dimensional';
        end
        disp(['Data:              ',SWITCH.GeodynamicCode,', ',dimString])
        
        %% fAIO ACTION: UPDATE CACHE
        fAIo.task   = 'handleCache';
        [fAIo,SAVE] = f_AIo(fAIo,SWITCH,FILE,[],STYLE,SAVE);
        
        %% GET DIMENSIONAL PARAMETERS
        [SETUP,SWITCH] = f_Dimensions(SWITCH,FILE);
        
        %% GET SUITABLE DEPTH VECTOR DIMENSIONS
        if SWITCH.DimensionalMode || SWITCH.DimensionalInput
            GRID.dimFactor      = SETUP.lengthscale;	%[m]
            %plotting variables
            if abs(max(mData(:,1)*SETUP.lengthscale))>1e3 % >1km
                GRID.dimFactor  = GRID.dimFactor /1e3;	%[km]
                GRID.Zdim       = 'km';
            elseif abs(max(mData(:,1)*SETUP.lengthscale))>1 % >1m
                GRID.dimFactor  = GRID.dimFactor;       %[m]
                GRID.Zdim       = 'm';
            else % <1m
                GRID.dimFactor = GRID.dimFactor *1e2;	%[cm]
                GRID.Zdim       = 'cm';
            end
        else
            GRID.ndFactor       = 1;                    %[nd]
            GRID.dimFactor      = GRID.ndFactor;        %[nd]
            GRID.Zdim           = 'nd';
        end
        
        %% SETUP DEPTH VECTOR
        %setup z-vector for grid cell tops or bottoms (e.g. vertical velocity)
        %  (is one value too large for rprof data, not sure if top or bottom.........)
        % mData(1,1)	 => bottom   => z=0
        % mData(end,1)   => surface  => z=1
        dzProf              = mData(2:end,1) - mData(1:end-1,1);
        dzBot               = (mData(1,1) - 0)*2;
        dzTop               = (1 - mData(end,1))*2;
        dzProf              = [dzBot; dzProf; dzTop]; %put additional dz for top and bottom boundary
        zEdgeProfFull       = mData(:,1) - dzProf(1:end-1,1)/2; %
        zEdgeProfFull       = [zEdgeProfFull; (zEdgeProfFull(end,1)+dzProf(end,1))];  %z from 0 to 1
        %     zEdgeProf          = zEdgeProfFull(2:end,1);   %z from ~0 to 1 %grid cell interior values, like T and p
        
        %% SET MIN/MAX VARIABLE VALUES
        if SWITCH.DimensionalMode || SWITCH.DimensionalInput
            zmin    = -500;     	zmax    = 3000;
            Tmin    = 0;            Tmax    = 2500;
            etamin  = 1e18;         etamax  = 1e29;
            edotmin = 1e-18;        edotmax = 1e-13;
            rhomin  = 2800;         rhomax  = 3400;
            strmin  = 0;            strmax  = 500;
            vmin    = 0;            vmax    = 3;
            Rhmin   = 0;            Rhmax   = 6e-12;
            hfmin   = 0;            hfmax   = 1000;
        else %non-dim
            zmin    = -0.05;     	zmax    = 0.95;
            Tmin    = 0;            Tmax    = 1;
            etamin  = 1e-4;         etamax  = 1e5;
            edotmin = 0;            edotmax = 7e4;
            rhomin  = -6e4;         rhomax  = 6e4;
            strmin  = 0;            strmax  = 7e4;
            vmin    = 0;            vmax    = 5000;
            Rhmin   = 0;            Rhmax   = 25;
            hfmin   = 0;            hfmax   = 0.5;
        end
        
        %% SETUP DATA DETAILS
        %   1                       2                       3                       4               5           6           7
        %  	name                    variable                dim. scale              dimension       min.value   max.value   logarithmic
        FIELD_M      = {
            'Depth'                 'z'                     GRID.dimFactor          GRID.Zdim       zmin        zmax        false        	%1
            'Temperature'           'T_{mean}'           	SETUP.Tscale            SETUP.TDim      Tmin      	Tmax        IN.T_logx    	%2
            'T1'                    'T_{min}'               SETUP.Tscale            SETUP.TDim      Tmin      	Tmax        IN.T_logx   	%3
            'T2'                    'T_{max}'               SETUP.Tscale            SETUP.TDim      Tmin      	Tmax        IN.T_logx     	%4
            'Velocity'              'v_{mean}'           	SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.v_logx       %5
            'V1'                    'v_{min}'               SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.v_logx       %6
            'V2'                    'v_{max}'               SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.v_logx       %7
            'Vertical velocity'     'vz_{mean}'           	SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.vz_logx      %8
            'Vz1'                   'vz_{min}'            	SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.vz_logx      %9
            'Vz2'                   'vz_{max}'           	SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.vz_logx      %10
            'Horizontal velocity'   'vh_{mean}'          	SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.vh_logx      %11
            'Vh1'                   'vh_{min}'            	SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.vh_logx      %12
            'Vh2'                   'vh_{max}'              SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.vh_logx  	%13
            'Viscosity'             '\eta_{mean}'           SETUP.etascale        	SETUP.etaDim 	etamin      etamax      IN.eta_logx  	%14
            'Eta1'                  '\eta_{min}'            SETUP.etascale        	SETUP.etaDim   	etamin      etamax      IN.eta_logx   	%15
            'Eta2'                  '\eta_{max}'            SETUP.etascale        	SETUP.etaDim   	etamin      etamax      IN.eta_logx   	%16
            'Strain rate'           'UCvarepsilondot_{mean}' SETUP.edotscale       	SETUP.edotDim  	edotmin     edotmax  	IN.edot_logx  	%17
            'E1'                    'UCvarepsilondot_{min}' SETUP.edotscale         SETUP.edotDim  	edotmin     edotmax  	IN.edot_logx  	%18
            'E2'                    'UCvarepsilondot_{max}' SETUP.edotscale         SETUP.edotDim  	edotmin     edotmax  	IN.edot_logx  	%19
            'Stress'                '\sigma_{mean}'         SETUP.stressscale       SETUP.stressDim	strmin      strmax      IN.str_logx   	%20
            'S1'                    '\sigma_{min}'          SETUP.stressscale       SETUP.stressDim	strmin      strmax      IN.str_logx  	%21
            'S2'                    '\sigma_{max}'          SETUP.stressscale       SETUP.stressDim	strmin      strmax      IN.str_logx   	%22
            'Vert. vorticity'    	'w_{h,mean}'           	SETUP.wscale        	SETUP.wDim      0           0.5         IN.w_logx       %23
            'W1'                   	'w_{h,min}'           	SETUP.wscale        	SETUP.wDim      0           0.5         IN.w_logx       %24
            'W2'                    'w_{h,max}'            	SETUP.wscale        	SETUP.wDim      0           0.5         IN.w_logx       %25
            'Horiz. vorticity'  	'w_{v,mean}'         	SETUP.wscale            SETUP.wDim      0           0.5         IN.w_logx       %26
            'W4'                    'w_{v,min}'            	SETUP.wscale            SETUP.wDim      0           0.5         IN.w_logx       %27
            'W5'                    'w_{v,max}'         	SETUP.wscale            SETUP.wDim      0           0.5         IN.w_logx       %28
            'Divergence'          	'div_{mean}'         	SETUP.divscale      	SETUP.divDim	0           1e-2       	IN.div_logx     %29
            'D1'                    'div_{min}'           	SETUP.divscale          SETUP.divDim	0           1e-2        IN.div_logx     %30
            'D2'                    'div_{max}'           	SETUP.divscale       	SETUP.divDim	0           1e-2        IN.div_logx     %31
            'Thermal advection'    	'H_{adv}'               SETUP.HRscale        	SETUP.HRDim   	0           1           IN.Hadv_logx    %32
            'Thermal diffusion'    	'H_{diff}'              SETUP.HRscale         	SETUP.HRDim  	0           1           IN.Hdiff_logx   %33
            'Internal heating'      'H_{int}'             	SETUP.Hscale         	SETUP.HDim      Rhmin       Rhmax       IN.Rh_logx      %34
            'Viscous dissipation'  	'H_{diss}'           	SETUP.HRscale          	SETUP.HRDim   	0           1           IN.vd_logx      %35
            'Adiabatic heating'    	'H_{adi}'             	SETUP.HRscale         	SETUP.HRDim   	0           1           IN.ah_logx      %36
            'Crust'                 'crust_{mean}'        	100                   	'wt%'        	0           100         IN.c_logx    	%37
            'C1'                    'c1'                    100                   	'wt%'          	0           100         IN.c_logx    	%38
            'C2'                    'c2'                    100                   	'wt%'         	0           100         IN.c_logx    	%39
            'Density'            	'rho_{mean}'         	SETUP.rhoscale      	SETUP.rhoDim  	rhomin      rhomax      IN.rho_logx     %40
            'Rho1'               	'rho_{min}'             SETUP.rhoscale          SETUP.rhoDim   	rhomin      rhomax      IN.rho_logx     %41
            'Rho2'                 	'rho_{max}'             SETUP.rhoscale        	SETUP.rhoDim   	rhomin      rhomax      IN.rho_logx     %42
            'Air'                 	'air_{mean}'         	100                   	'wt%'          	0           100         IN.c_logx    	%43
            'Air1'                  'air_{min}'            	100                   	'wt%'         	0           100         IN.c_logx      	%44
            'Air2'                 	'air_{max}'           	100                   	'wt%'        	0           100         IN.c_logx      	%45
            'Primordial'           	'prim_{mean}'       	100                   	'wt%'         	0           100         IN.c_logx    	%46
            'Prim1'                 'prim_{min}'        	100                   	'wt%'          	0           100         IN.c_logx     	%47
            'Prim2'                	'prim_{max}'           	100                   	'wt%'          	0           100         IN.c_logx     	%48
            'Cont. crust'        	'cc_{mean}'         	100                   	'wt%'          	0           100         IN.c_logx     	%49
            'Cc1'                   'cc_{min}'            	100                   	'wt%'          	0           100         IN.c_logx      	%50
            'Cc2'                 	'cc_{max}'           	100                   	'wt%'          	0           100         IN.c_logx    	%51
            'Fluids'               	'f_{mean}'              100                   	'wt%'          	0           100         IN.c_logx     	%52
            'F1'                    'f_{min}'               100                   	'wt%'          	0           100         IN.c_logx   	%53
            'F2'                	'f_{max}'           	100                   	'wt%'          	0           100         IN.c_logx     	%54
            'Metal'               	'metal_{mean}'       	100                   	'wt%'          	0           100         IN.c_logx     	%55
            'Metal1'              	'metal_{min}'        	100                   	'wt%'          	0           100         IN.c_logx    	%56
            'Metal2'            	'metal_{max}'        	100                   	'wt%'          	0           100         IN.c_logx    	%57
            'Water'               	'water_{mean}'       	100                   	'wt%'          	0           100         IN.c_logx     	%58
            'Water1'              	'water_{min}'        	100                   	'wt%'          	0           100         IN.c_logx    	%59
            'Water2'            	'water_{max}'        	100                   	'wt%'          	0           100         IN.c_logx    	%60
            ... %THESE STILL NEED TO BE FULLY IMPLEMENTED..........................
            'Visc. Dissipation Log'	'theta_{log}'        	1                   	'?'          	0           1           false           %61
            'Visc. Dissipation 1' 	'theta_{min}'        	1                   	'?'          	0           1           false           %62
            'Visc. Dissipation 2'  	'theta_{max}'        	1                   	'?'          	0           1           false           %63
            'Advection Total'      	'adv_{tot}'             1                   	'?'          	0           1           false           %64
            'Advection 1'         	'adv_{descending}'    	1                   	'?'          	0           1           false           %65
            'Advection 2'          	'adv_{ascending}'     	1                   	'?'          	0           1           false           %66
            'Thermal conduction'  	'tcond_{max}'        	1                   	'?'          	0           1           false           %67
            'Thermal conduction 1' 	'tcond_{max}'        	1                   	'?'          	0           1           false           %68
            'Thermal conduction 2' 	'tcond_{max}'        	1                   	'?'          	0           1           false           %69
            ...%
            ... %THIS IS OPTIONAL IN STAGYY
            'Toroidal prof'         'to'                    SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.To_logx      %70
            'Toroidal rms'          'to_{rms}'              SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.To_logx      %71
            'Poloidal prof'         'po'                    SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.Po_logx      %72
            'Poloidal rms'        	'po_{rms}'              SETUP.Vscale            SETUP.vDim      vmin       	vmax        IN.Po_logx      %73
            'Mobility (vrms)'       'm'                     1                       '-'             0           1           IN.M_logx       %74
            'Plateness (f80)'       'f80'                   1                       '-'             0           1           IN.f80_logx   	%75
            'Plateness (f90)'     	'f90'                   1                       '-'             0           1           IN.f90_logx   	%76
            ... %THESE ARE READ FROM REFSTAT.DAT
            'Ref. state temperature' 'T_{refstate}'         SETUP.Tscale            SETUP.TDim      Tmin      	Tmax        IN.T_logx    	%automatically found
            'Ref. state density'    'rho_{refstate}'        SETUP.rhoscale      	SETUP.rhoDim  	rhomin      rhomax      IN.rho_logx    	%automatically found
            'Ref. state expansivity' 'alpha_{refstate}'    	SETUP.alpha           	'K^{-1}'      	0           1           false           %automatically found
            'Ref. state conductivity' 'Tcond_{refstate}'    SETUP.tcond_dimensional 'Wm^{-1}K^{-1}' 0           1           false           %automatically found
            'Ref. state pressure'  	'P_{refstate}'          SETUP.Pscale          	SETUP.PDim    	0           1           false           %automatically found
            ... %THIS IS CALCULATED LOCALLY IN SL_RadialProfile
            'Heat flux'             'hf_{mean}'             1                       'W/m^2'      	hfmin       hfmax       IN.HF_logx    	%end-2
            'Hf1'                   'hf_{min}'              1                       'W/m^2'     	hfmin       hfmax       IN.HF_logx     	%end-1
            'Hf2'                   'hf_{max}'              1                       'W/m^2'       	hfmin       hfmax       IN.HF_logx    	%end
            };

        %% UNICODE VARIABLES
        %use the following character sequence to add unicode variables:
        %UCvarepsilondot => char(941)
        %UCvarepsilon => char(949)
        FIELD_M(:,2) = strrep(FIELD_M(:,2),'UCvarepsilondot',{char(941)}); %does only work if string consists ONLY of unicode character
        FIELD_M(:,2) = strrep(FIELD_M(:,2),'UCvarepsilon',{char(949)}); %does only work if string consists ONLY of unicode character
        
        %     zg(1,iz),(tprof(j,iz),j=1,3),(vprof(j,iz),j=1,3),(vzprof(j,iz),j=1,3),&
        %                       (vhprof(j,iz),j=1,3),(etaprof(j,iz),j=1,3),(eprof(j,iz),j=1,3)
        %                  write(1,'(38(1pe12.4))') (sprof(j,iz),j=1,3),((wprof(j,iz,k),j=1,3),k=1,2), &
        %                       (dprof(j,iz),j=1,3),(enprof(j,iz),j=1,5),(cprof(j,iz),j=1,3),(denprof(j,iz),j=1,3), &
        %                       (airprof(j,iz),j=1,3),(primprof(j,iz),j=1,3),(ccprof(j,iz),j=1,3),(fprof(j,iz),j=1,3), &
        %                       (metalprof(j,iz),j=1,3)
        
        %% FIND REFSTAT INDICES
        parameterX(1,ismember(fieldString(:),'RefS Temperature')) = find(strcmpi(FIELD_M(:,1),'Ref. state temperature'));
        parameterX(1,ismember(fieldString(:),'RefS Density')) = find(strcmpi(FIELD_M(:,1),'Ref. state density'));
        parameterX(1,ismember(fieldString(:),'RefS Expansivity')) = find(strcmpi(FIELD_M(:,1),'Ref. state expansivity'));
        parameterX(1,ismember(fieldString(:),'RefS Conductivity')) = find(strcmpi(FIELD_M(:,1),'Ref. state conductivity'));
        parameterX(1,ismember(fieldString(:),'RefS Pressure')) = find(strcmpi(FIELD_M(:,1),'Ref. state pressure'));

        %% SET TIME
        PLOT.time           = time;     %[nd] or [dim]
        PLOT.time_dim       = time*SETUP.timescale/SETUP.secyear;  %get dimensional time [yr]
        if PLOT.time_dim/1e6>100
            PLOT.time2plot  = PLOT.time_dim/1e9; PLOT.time2plotDim = 'Gyr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear/1e9;
        elseif PLOT.time_dim/1e6<1e-7 %0.0000001
            PLOT.time2plot  = PLOT.time_dim*SETUP.secyear; PLOT.time2plotDim = 's'; PLOT.timeConvert = SETUP.timescale;
        elseif PLOT.time_dim/1e6<1e-4 %0.0001
            PLOT.time2plot  = PLOT.time_dim; PLOT.time2plotDim = 'yr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear;
        elseif PLOT.time_dim/1e6<1e-1 %0.1
            PLOT.time2plot  = PLOT.time_dim/1e3; PLOT.time2plotDim = 'kyr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear/1e3;
        else
            PLOT.time2plot  = PLOT.time_dim/1e6; PLOT.time2plotDim = 'Myr'; PLOT.timeConvert = SETUP.timescale/SETUP.secyear/1e6;
        end
        if ~SWITCH.DimensionalMode && ~SWITCH.DimensionalInput; PLOT.time2plot = PLOT.time; PLOT.time2plotDim = 'nd'; PLOT.timeConvert = 1; end
        
        %% INDICATE TIME
        timeString0	= [num2str(PLOT.time,3),' | ',num2str(PLOT.time2plot,3),' ',PLOT.time2plotDim];
        if ~SWITCH.DimensionalMode && ~SWITCH.DimensionalInput; timeString0	= num2str(PLOT.time,3); end
        timeString	= timeString0;
        %     if ~isnan(PLOT.timeMT)
        %         timeString = [timeString0,' | ',num2str(PLOT.timeMT,2),' mantle turnovers (vrms = ',num2str(v_rmsCurrent,2),', dt_MT = ',num2str(PLOT.dtMT,2),' ',PLOT.time2plotDim,')'];
        %     end
        disp(['Time:              ',timeString])
        
        %% CALCULATE HEATFLUX FROM T-PROFILE
        if PLOT.Heatflux
            idxHeatflux                 = size(FIELD_M,1)-2;
            if ~SWITCH.DimensionalMode && ~SWITCH.DimensionalInput
                FIELD_M(idxHeatflux:idxHeatflux+2,4) = {'nd'};
                warning('Heatflux needs dimensional mode switched on!');
            end
            
            %adjust parameters calculated locally (due to optional rprof variables)
            parameterX(1,ismember(fieldString(:),'Heatflux')) = size(FIELD_M,1)-2;   %locate heatflux parameters as 3rd last entry
            
            zVector                     = mData(:,1);  %*dim_scale(PLOT.ParameterY);
            dz                          = zVector(2:end) - zVector(1:end-1);
            zVectorGradient             = zVector(2:end) - dz(1:end)./2;  % z_new = z_old + dz/2
            zVectorGradient(end+1)      = NaN; %max(zVectorGradient);
            
            T_prof                      = mData(:,2)*FIELD_M{2,3};
            T_prof_min                  = mData(:,3)*FIELD_M{2,3};
            T_prof_max                  = mData(:,4)*FIELD_M{2,3};
            
            HF_prof                     = SETUP.tcond_dimensional .*( T_prof(1:end-1)-T_prof(2:end) )./dz(1:end); %calculate heat flux from temperature: Q = k*[T(z+1)-T(z)]/dz  [W/m^2]  (upward positive)
            HF_prof_min                 = SETUP.tcond_dimensional .*( T_prof_min(1:end-1)-T_prof_min(2:end) )./dz(1:end);
            HF_prof_max                 = SETUP.tcond_dimensional .*( T_prof_max(1:end-1)-T_prof_max(2:end) )./dz(1:end);
            
            mData(1:end-1,parameterX(strcmp(fieldString,'Heatflux'))) = HF_prof;          mData(end,parameterX(strcmp(fieldString,'Heatflux'))) = NaN;
            mData(1:end-1,parameterX(strcmp(fieldString,'Heatflux'))+1) = HF_prof_min;    mData(end,parameterX(strcmp(fieldString,'Heatflux'))+1) = NaN;
            mData(1:end-1,parameterX(strcmp(fieldString,'Heatflux'))+2) = HF_prof_max;    mData(end,parameterX(strcmp(fieldString,'Heatflux'))+2) = NaN;
        end
        
        %% REVERSE Z-AXIS
        if IN.reverseYaxis
            if SWITCH.DimensionalInput
                maxZ        = SETUP.D;
            else
                maxZ        = 1.0;
            end
            mData(:,1)      = maxZ-mData(:,1);
            % zEdgeProf     = max_z-zEdgeProf; %need to take max_z for all z-data! %grid cell interior values, like T and p
            if PLOT.Heatflux; zVectorGradient = maxZ-zVectorGradient; end
            if readRefstat; RS.z = maxZ-RS.z; end
            if readRefstat; RS.CA_z = maxZ-RS.CA_z; end

            mData(:,1)      = mData(:,1)-SETUP.d_air;  %set actual depth, i.e., z=0 to rock-air interface
            % zEdgeProf = zEdgeProf-SETUP.d_air;   %set z=0 to real surface %grid cell interior values, like T and p
            if PLOT.Heatflux; zVectorGradient(:,1) = zVectorGradient(:,1)-SETUP.d_air; end
            if readRefstat; RS.z = RS.z-SETUP.d_air; end
            if readRefstat; RS.CA_z = RS.CA_z-SETUP.d_air; end
        end
        
        %% PLOTTING
        numSubplot = 0;
        nrUsedPlots = nrPlots;
        hax = zeros(nrPlots,1);
        LeftmostSubplot = zeros(nrPlots,1);
        for iPlot=1:nrPlots
            [ixPlot,izPlot] = ind2sub([nxSubplots,nzSubplots],iPlot);
            %extra axes for outside legend
            if SWITCH.LegendOutside
                nrUsedPlots    = nrPlots-1;
                if iPlot==nrPlots
                    
                    continue
                end
            end
            fieldFound          = true;
            plotParameterX      = parameterX(iPlot);
            FIELD.logX          = FIELD_M{plotParameterX,7};
            FIELD.logY          = false;
            %check for separate min/max variables
            if plotParameterX<=31 || (plotParameterX>=37 && plotParameterX<=57) %separate min/max values available
                minmaxAvailable = true;
            else
                minmaxAvailable = false;
            end
            %choose appropriate zVector:
            currentFieldString = fieldString(iPlot);
            if strcmp(currentFieldString,'Heatflux')
                zData           = zVectorGradient;  %for heat flux
            elseif length(currentFieldString{1})>=4 && strcmp(currentFieldString{1}(1:4),'RefS')
                if ~readRefstat
                    fieldFound      = false;
                    warning off backtrace; warning('refstat.dat could not be read.'); warning on backtrace
                else
                    zData           = RS.z(:,1); % .*SETUP.lengthscale;
                end
            else
                zData           = mData(:,PLOT.ParameterY);
            end
            
            %check if variable exists
            if plotParameterX>numEntries && ~(length(currentFieldString{1})==4 && strcmp(currentFieldString{1}(1:4),'RefS')) %variable number exceeds total number of entries 
                fieldFound      = false;
                warning off backtrace
                warning(['Variable "',FIELD_M{plotParameterX,1},'" could not be found in the current rprof.dat file!'])
                warning on backtrace
            end
            
            if ~SWITCH.QuickMode
                %% SETUP FIGURE
%                 figure(1); hold on
                
                %% SETUP SUBPLOT
                if nrPlots>1
                    numSubplot      = numSubplot+1;
                    if SWITCH.UsePanel
                        hax(iPlot,1) 	= SAVE.P(izPlot,ixPlot).select();
                        set(gcf,'CurrentAxes',hax(iPlot,1))
                    else
                        hax(iPlot,1)    = subplot(nzSubplots,nxSubplots,numSubplot);
                    end
                end
                if ixPlot==1
                    LeftmostSubplot(iPlot,1) = true;
                end
                
                %% SETUP PLOT APPEARANCE
                if iPlot==1; IN.LineWidthPlot = IN.LineWidthPlot+0.1; end   %wider first line
                
                if ~fieldFound
                    %insert empty plot
                    FIELD.name              = FIELD_M{plotParameterX,1}; %just to be defined
                    [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
                    
                    continue
                    
                end
            end
            
            %% DIMENSIONALISATION
            zData       = zData .*FIELD_M{PLOT.ParameterY,3};
            if length(currentFieldString{1})>=4 && strcmp(currentFieldString{1}(1:4),'RefS') %refstat
                %Combined adiabat profiles
                rmsData     = RS.(varString{iPlot})(:,1) .*FIELD_M{plotParameterX,3};
            else %rprof
                rmsData     = mData(:,plotParameterX) .*FIELD_M{plotParameterX,3};
            end
            %automatically adjust mass fractions
            if strcmp(FIELD_M{plotParameterX,4},'wt%') && max(rmsData)<1.0 %change to [ppm]
                FIELD_M{plotParameterX,3}   = FIELD_M{plotParameterX,3}*1e4; %convert [wt%] to [ppm]
                FIELD_M{plotParameterX,4}   = 'ppm';
                %re-dimensionalise
                rmsData     = mData(:,plotParameterX) .*FIELD_M{plotParameterX,3};
            end
            
            if ~SWITCH.QuickMode
                if PLOT.Solidus && strcmp(FIELD_M{plotParameterX,1},'Temperature')
                    %% CALCULATE SOLIDUS
                    % based on simple fits of the Herzberg upper mantle solidus and Zerr et al lower
                    % mantle solidus. Solidii were first plotted in (depth,T) space before fitting.
                    % A simple linear + erf fit works very well and extrapolates sensibly
                    if SWITCH.DimensionalMode
                        suppressAsth = false;
                        dummy       = zData<=660; %upper mantle indices
                        %lower mantle
                        %********************************* % calculates erfnr(x);
                        x           = zData(~dummy)./1000;
                        z           = abs(x);
                        t           = 1./(1+0.5.*z);
                        erfnr       = 1-t.*exp(-z.*z-1.26551223+t.*(1.00002368+t.*(0.37409196+  ...
                            t.*(0.09678418+t.*(-0.18628806+t.*(0.27886807+t.*(-1.13520398+ ...
                            t.*(1.48851587+t.*(-0.82215223+t.*0.17087277)))))))));
                        erfnr(x<0)      = -erfnr(x<0);
                        %*********************************
                        simpleSolidus           = zeros(size(zData));
                        simpleSolidus(~dummy)   = (2760 +0.45.*zData(~dummy) +1700.*(erfnr-1)) ./SETUP.Tscale;
                        %upper mantle
                        if suppressAsth
                            simpleSolidus(dummy) = (2050 +0.62.*zData(dummy)) ./SETUP.Tscale;
                        else
                            %********************************* % calculates erfnr(x);
                            x   = zData(dummy)./220;
                            z   = abs(x);
                            t   = 1./(1+0.5.*z);
                            erfnr   = 1-t.*exp(-z.*z-1.26551223+t.*(1.00002368+t.*(0.37409196+  ...
                                t.*(0.09678418+t.*(-0.18628806+t.*(0.27886807+t.*(-1.13520398+ ...
                                t.*(1.48851587+t.*(-0.82215223+t.*0.17087277)))))))));
                            erfnr(x<0) = -erfnr(x<0);
                            %*********************************
                            simpleSolidus(dummy) = (2050 +0.62.*zData(dummy) +660.*(erfnr-1)) ./SETUP.Tscale;
                        end
                        simpleSolidus   = simpleSolidus .*SETUP.Tscale; %dimensionalisation
                        clearvars x z t erfnr
                    else
                        warning('Solidus not implemented for non-dimensional mode!')
                    end
                    
                end
                
                %% PLOT RMS LINE (and min/max lines)
                colourMax       = [0.81 0.33 0.2];
                colourMin       = [0.4 0.64 0.18];
                colourSolidus   = [0.9 0.6 0.6];
                if FIELD.logX %LOGARITHMIC x_axis
                    semilogx(rmsData,zData,'Color',lineColor,'LineWidth',IN.LineWidthPlot); %needed to change x axis to semilogx
                    if PLOT.Solidus && strcmp(FIELD_M{plotParameterX,1},'Temperature')
                        hold on
                        hSolidus = semilogx(simpleSolidus,zData,'Color',colourSolidus,'LineWidth',IN.LineWidthPlot);
                    end
                    if PLOT.parameterRange && minmaxAvailable
                        minData = mData(:,plotParameterX+1)*FIELD_M{plotParameterX,3};
                        maxData = mData(:,plotParameterX+2)*FIELD_M{plotParameterX,3};
                        if PLOT.rangeAsArea
                            %logarithmic min/max area
                            edgecolor = 'none';
                            areaAlpha = 0.2;
                            hold on
                            pp4 = fill([maxData' fliplr(minData')],[zData' fliplr(zData')],...
                                lineColor,'EdgeColor',edgecolor,'EdgeAlpha',areaAlpha,'FaceAlpha',areaAlpha);
                        else
                            %logarithmic min/max lines
                            hold on
                            pp2 = semilogx(minData,zData,'Color',colourMin);
                            pp3 = semilogx(maxData,zData,'Color',colourMax);
                        end
                    end
                    hold on
                    pp1 = semilogx(rmsData,zData,'Color',lineColor,'LineWidth',IN.LineWidthPlot);
                    
                else %NORMAL x_axis
                    if PLOT.Solidus && strcmp(FIELD_M{plotParameterX,1},'Temperature')
                        hSolidus = plot(simpleSolidus,zData,'Color',colourSolidus,'LineWidth',IN.LineWidthPlot);
                        hold on;
                    end
                    if PLOT.parameterRange && minmaxAvailable
                        minData = mData(:,plotParameterX+1)*FIELD_M{plotParameterX,3};
                        maxData = mData(:,plotParameterX+2)*FIELD_M{plotParameterX,3};
                        if PLOT.rangeAsArea
                            %min-max area
                            edgecolor = 'none';
                            areaAlpha = 0.2;
                            hold on
                            pp4 = fill([maxData' fliplr(minData')],[zData' fliplr(zData')],...
                                lineColor,'EdgeColor',edgecolor,'EdgeAlpha',areaAlpha,'FaceAlpha',areaAlpha);
                        else
                            %min/max lines
                            hold on
                            pp2 = plot(minData,zData,'Color',colourMin);
                            pp3 = plot(maxData,zData,'Color',colourMax);
                        end
                    end
                    hold on
                    pp1 = plot(rmsData,zData,'Color',lineColor,'LineWidth',IN.LineWidthPlot);
                end
                
                %% INDIVIDUAL LEGENDS
                if PLOT.parameterRange && minmaxAvailable
                    if PLOT.rangeAsArea
                        hRunsArea(ifile,1) = pp4;
                        if IN.showLegend
                            legend([pp1,pp4],[FIELD_M(plotParameterX,2),'Min-Max Area'],'Location','Best');
                        end
                    else %plot min/max as lines
                        if IN.showLegend
                            legend([pp1,pp2,pp3],[FIELD_M(plotParameterX,2),FIELD_M(plotParameterX+1,2),FIELD_M(plotParameterX+2,2)],'Location','Best');
                        end
                    end
                else
                    if IN.showLegend
                        legend(pp1,FIELD_M(plotParameterX,2),'Location','Best');
                    end
                end
                hRuns(ifile,1)     = pp1;
                
                %% AXIS LABELLING
                FIELD.name              = FIELD_M{plotParameterX,1}; %just to be defined
                DESVARIA.xlabelName     = FIELD_M{plotParameterX,1};
                DESVARIA.xlabelDim      = FIELD_M{plotParameterX,4};
                DESVARIA.zlabelName     = FIELD_M{PLOT.ParameterY,1};
                DESVARIA.zlabelDim      = FIELD_M{PLOT.ParameterY,4};
                DESVARIA.Task           = 'create annotation strings';
                [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,[],SWITCH,PLOT,STYLE,SAVE);
                xlabel(PLOT.xlabel)
                if SWITCH.SimplifyPlots && ~LeftmostSubplot(iPlot,1)
                    %no ylabel
                else
                    ylabel(PLOT.ylabel);
                end
                
                %% AXIS SETTINGS
                axisLimitX          = xlim; %current x-axis limits
                axisLimitY          = ylim; %current y-axis limits
                grid on
                if STYLE.MinorAxisTicks; grid minor; end
                if IN.reverseYaxis
                    DESVARIA.flipAxes = true;
                    ymin            = min(zData(:));
                    ymax            = max(zData(:));
                    if ymin<yMinOld
                        DESVARIA.axesLimits     = [axisLimitX,ymin,ymax];
                        if iPlot==nrUsedPlots; yMinOld = ymin; end
                    end
                end
                %constant axes
                if IN.constantAxis
                    DESVARIA.axesLimits     = [FIELD_M{plotParameterX,5},FIELD_M{plotParameterX,6},zmin,zmax];
                else
                    if ~FIELD.logX
                        %check for sensible x-axis values - preventing 0.00001 values
                        xMinMax         = get(gca,'xlim');
                        xminNearInt     = round(xMinMax(1,1));
                        xminDiffNearInt = abs(xMinMax(1,1)-xminNearInt);
                        if xminDiffNearInt~=0 && xminDiffNearInt<(xMinMax(1,2)-xMinMax(1,1))/100 %100 is an arbitrary value
                            set(gca,'xlim',round(get(gca,'xlim')))
                            DESVARIA.axesLimits     = [round(get(gca,'xlim')),axisLimitY];
                        end
                    end
                end
                %set manual y-axis limits
                if SWITCH.AxesLimit
                    if ~isnan(IN.YAxisMinValue)
                        axisLimitY(1,1) 	= IN.YAxisMinValue;
                    end
                    if ~isnan(IN.YAxisMaxValue)
                        axisLimitY(1,2) 	= IN.YAxisMaxValue;
                    end
                    DESVARIA.axesLimits     = [axisLimitX,axisLimitY];
                end
                DESVARIA.Task    = 'setup axes';
                [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            end
            
            %% GET MIN/MEAN/MAX VALUES
            if minmaxAvailable %separate min/max values available
                minVal      = min(mData(:,plotParameterX+1)*FIELD_M{plotParameterX,3});
                maxVal      = max(mData(:,plotParameterX+2)*FIELD_M{plotParameterX,3});
            else %no separate min/max values available
                if length(currentFieldString{1})>=4 && strcmp(currentFieldString{1}(1:4),'RefS')
                    %Combined adiabat profiles
                    minVal     = min(RS.(varString{iPlot})(:,1)*FIELD_M{plotParameterX,3});
                    maxVal     = max(RS.(varString{iPlot})(:,1)*FIELD_M{plotParameterX,3});
                else
                    minVal     = min(mData(:,plotParameterX)*FIELD_M{plotParameterX,3});
                    maxVal     = max(mData(:,plotParameterX)*FIELD_M{plotParameterX,3});
                end
            end
            
            %% DISPLAY PARAMETER INFO
            if iPlot==1; disp('Field range:'); end
            tabLength       = 20;
            lengthString    = size(fieldString{1,iPlot},2);
            whiteSpace(1,1:max(1,tabLength-lengthString)) = ' ';
            disp(['   ',STYLE.SCHAR.smallBullet,' ',fieldString{1,iPlot},whiteSpace,' ',num2str(minVal,4),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(maxVal,4),' ',FIELD_M{plotParameterX,4}])
            clearvars whiteSpace tabLength lengthString
            
            %% TEXTING
            if SWITCH.Texting
                if IN.reverseYaxis; maxYlim = min(ylim); else; maxYlim = max(ylim); end
                text(max(xlim),maxYlim,[num2str(minVal,3),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(maxVal,3),' ',FIELD_M{plotParameterX,4}],'HorizontalAlignment','right','VerticalAlignment','bottom','Color',STYLE.keyColor,'FontSize',8)
                clearvars maxYlim minVal maxVal
            end
            
            %% ADDITIONAL POST-PROCESSING
            if SWITCH.d_eta && strcmp(fieldString(iPlot),'Viscosity')
                topDepthLevel       = 500;
                botDepthLevel       = 2500;
                yyy                 = mData(:,PLOT.ParameterY)*FIELD_M{PLOT.ParameterY,3};
                xxx                 = mData(:,plotParameterX+1)*FIELD_M{plotParameterX,3};
                topValy             = min(yyy(yyy>topDepthLevel));
                topValx             = xxx(yyy==topValy);
                botValy             = max(yyy(yyy<botDepthLevel));
                botValx             = xxx(yyy==botValy);
                hold on
                indicatorColor      = [0.1 0.2 0.9];
                if strcmp(STYLE.ColorMode,'dark')
                    indicatorColor  = 1-indicatorColor;
                end
                plot([topValx,botValx],[topValy,botValy],'o','MarkerEdgeColor',indicatorColor)
                
                delta_visc          = botValx-topValx;
                delta_d             = botValy-topValy;
                gradient_eta        = delta_visc/delta_d*1000;
                gradient_etaMantle 	= delta_visc/delta_d*2890;
                d_etaFactor         = (botValx+(delta_visc/delta_d*(2890-botDepthLevel)))/topValx;
                d_etaOrdOfMag   	= (log10(botValx+(delta_visc/delta_d*(2890-botDepthLevel)))) - log10(topValx);
                
                text(topValx,botValy-(botValy-topValy)/2,{['\Delta\eta = ',num2str(gradient_eta,3),' Pas/1000km'];['\Delta\eta = ',num2str(d_etaOrdOfMag,2),' O.o.M.']},'Color',PLOT.ColorVector(1,:),'backgroundcolor',1-PLOT.ColorVector(1,:));
                disp(['(d_eta = ',num2str(gradient_eta,3),' Pas/1000km)']);
                disp(['(d_eta = ',num2str(gradient_etaMantle,3),' Pas/2890km)']);
                disp(['(d_eta = ',num2str(d_etaFactor,2),',  ',num2str(topValx,3),' <-> ',num2str(botValx+(delta_visc/delta_d*(2890-botDepthLevel)),3),')']);
                disp(['(d_eta = ',num2str(d_etaOrdOfMag,2),',  ',num2str(log10(topValx),3),' <-> ',num2str(log10(botValx+(delta_visc/delta_d*(2890-botDepthLevel))),3),' orders of magnitude)']);
            end
            
            %% REMEMBER PROFILE DATA TO SAVE
            if SAVE.Profiles
                %create array to save
                if iPlot==1; titleArray = {[FIELD_M{PLOT.ParameterY,1},'[',FIELD_M{PLOT.ParameterY,4},']']}; dataArray = zData; end
                titleArray              = [titleArray, {[FIELD_M{plotParameterX,1},'[',FIELD_M{plotParameterX,4},']']}];
                dataArray               = [dataArray, rmsData];
            end
            %pause(0.5)
        end %loop
        
        %% SAVING PROFILE DATA
        if SAVE.Profiles
            SAVE.Directory     	= FILE.directory;
            SAVE.DataName      	= [FILE.name,'_Profiles',num2str(FILE.number)];
            SAVE.DataTitle    	= strjoin(titleArray,', ');
            SAVE.data          	= dataArray;
            [SAVE.overwriteAll] = f_saveData( SAVE );
            SAVE                = rmfield(SAVE,{'data','DataTitle'});
        end
        
        %% GENERAL LABELLING
        %figure title
        if SWITCH.Title
            figure(1); hold on
            axes('position',[0 0 1 1],'visible','off')
            text(0.5,0.97,['File number = ',num2str(FILE.number-1),'; Time step = ',num2str(i_time_step),'; Time = ',num2str(PLOT.time,3),' ',PLOT.time2plotDim],'HorizontalAlignment','center')
        end
        fileNameString      = [fileNameString,'_',FILE.name];
        %pause(0.5)
        disp(' ')
    end
    
    %% LEGEND
    if ~SWITCH.QuickMode
        if isfield(PLOT,'titleString')
            if length(hRuns)>length(PLOT.titleString); PLOT.titleString{length(PLOT.titleString)+1:length(hRuns)} = '[undefined]'; end
            legendANN = PLOT.titleString';
        else
            legendANN = IN.Name';
        end
        
        for ientry=1:size(legendANN,1)
            dummy               = ['RMS',' ',legendANN{ientry,1}];
            legendANNfull(1,ientry) = cellstr(dummy);
        end
        
        try
        if PLOT.parameterRange && PLOT.rangeAsArea && hRunsArea(1,1)~=0
            for ientry=1:size(legendANN,1)
                dummy           = ['Range',' ',legendANN{ientry,1}];
                legendANNrange(1,ientry) = cellstr(dummy);
            end
            dummy               = [legendANNfull; legendANNrange];
            legendANN_combined  = dummy(:)';
            dummy               = [hRuns'; hRunsArea'];
            hRuns_combined      = dummy(:)';
            hl = legend(hRuns_combined',legendANN_combined{:,:},'Location',PLOT.LegendLocation);
        else
            hl = legend(hRuns,legendANNfull{:,:},'Location',PLOT.LegendLocation);
        end
        set(hl,'Interpreter','none');
        if SWITCH.LegendOutside
            numSubplot          = numSubplot+1;
            if SWITCH.UsePanel
                [ixPlot,izPlot] = ind2sub([nxSubplots,nzSubplots],numSubplot);
                hax(iPlot,1) 	= SAVE.P(izPlot,ixPlot).select();
            else
                hax(iPlot,1)  	= subplot(nzSubplots,nxSubplots,numSubplot);
            end
            EmptyLegendSpacePosition = get(hax(end,1),'Position');
            axis off
            LegendPosition      = get(hl,'Position');
            LegendPosition(1)   = EmptyLegendSpacePosition(1);
            LegendPosition(2)   = EmptyLegendSpacePosition(2)+EmptyLegendSpacePosition(4)-LegendPosition(4);
            set(hl,'Position',LegendPosition);
        end
        catch
            if SWITCH.Verbose
                warning off backtrace; warning('Legend could not be added; probably because of missing data.'); warning off backtrace;
            end
        end
    end
    
    %% STYLING
    if SWITCH.PlotDesign
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
        if SWITCH.LegendOutside
            axes(hax(iPlot-1,1)) %activate second last subplot (i.e., one with a graphics object)
        end
    end
    
    %% fAIO ACTION: FINISHING UP
    fAIo.task = 'finishingUp';
    [fAIo,SAVE] = f_AIo(fAIo,SWITCH,FILE,GRID,STYLE,SAVE);
    
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
        endix                       = [num2str(FILE.number),strjoin(varString,'')];
        if SWITCH.AnalysisMode;             endix = [endix,'An'];   	end
        if SWITCH.DimensionalMode;          endix = [endix,'Dim'];   	end
        if strcmp(STYLE.ColorMode,'dark');  endix = [endix,'Dark'];   	end
        SAVE.FigureName             = strcat(IN.saveString,fileNameString,'_',endix);   %file name
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
[fAIo,SAVE] = f_AIo(fAIo,SWITCH,FILE,GRID,STYLE,SAVE);

end







