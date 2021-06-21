
%%                                                  MANTLE DIAGNOSTICS 2.33
%
% calculates connected areas by accounting for MD.addFactor times the area of
% the original plot: one part added at x=0 and one other part added at x=end
%
% evaluates the residual temperature
% field indicating +1 for hot plumes and -1 for cold plumes and 0 for the residual field
% ( hot plumes are defined by e.g. Tmean+0.5*(Tmax-Tmean) )
%
%                                                Fabio Crameri, 23.03.2020

function [MANTLE] = f_DiagnosticsMantle(FILE,GRID,SETUP,SWITCH,PLOT,STYLE)

%% DEFAULTS
MD.dummy = [];
if ~isfield(MD,'plumeDefinition');      	MD.plumeDefinition      	= 2;            end     %1: solely by temperature, 2: Temperature & radial velocity
if ~isfield(MD,'plotPLUMES');               MD.plotPLUMES               = logical(0); 	end     %plot comparison of all anomalies and selected plumes
if ~isfield(MD,'addFactor');                MD.addFactor                = 0.25;         end     %percent of total x-extend that is added at x=1 and x=end
if ~isfield(MD,'upperDepthThresholdHot'); 	MD.upperDepthThresholdHot 	= 0.8;          end     %upper threshold; fraction of mantle depth (1: bottom, 0: top)
if ~isfield(MD,'lowerDepthThresholdHot'); 	MD.lowerDepthThresholdHot 	= 0.9;          end     %lower threshold; fraction of mantle depth (1: bottom, 0: top)
if ~isfield(MD,'upperDepthThresholdCold');	MD.upperDepthThresholdCold	= 0.15;         end     %upper threshold; fraction of mantle depth (1: bottom, 0: top)
if ~isfield(MD,'lowerDepthThresholdCold');	MD.lowerDepthThresholdCold	= 0.25;     	end     %lower threshold; fraction of mantle depth (1: bottom, 0: top)
if ~isfield(MD,'residualMode');             MD.residualMode             = 2;            end     %1: horizontal, 2: horizontal band, 3: global, 4: regional

%% DEFAULT OUTPUT FIELDS
MANTLE.DiagnosticsFailed        = false;
MANTLE.upWelling                = NaN;      %1 for upwelling, 0 else (defined by threshold)
MANTLE.downWelling              = NaN;      %1 for downwelling, 0 else (defined by threshold)
MANTLE.upWellingAbsolute      	= NaN;      %1 for upwelling, 0 else
MANTLE.downWellingAbsolute   	= NaN;      %1 for downwelling, 0 else
MANTLE.plumesHot                = NaN;      %1 for hot plumes, 0 else
MANTLE.plumesCold               = NaN;  	%1 for cold plumes, 0 else
MANTLE.numHotPlumes             = NaN;    	%num
MANTLE.numColdPlumes            = NaN;  	%num
MANTLE.UpwellingVolume          = NaN;   	%[nd] or [m^2]or[m^3]
MANTLE.DownwellingVolume        = NaN;  	%[nd] or [m^2]or[m^3]
MANTLE.UpwellingVolPerc         = NaN;      %in [% of model domain]
MANTLE.DownwellingVolPerc       = NaN;  	%in [% of model domain]

MANTLE.VHmaxUM                 	= NaN;    	%[plotting dim]
MANTLE.VHmaxMM                 	= NaN;    	%[plotting dim]
MANTLE.VHmaxLM                 	= NaN;    	%[plotting dim]
MANTLE.VHmeanUM               	= NaN;    	%[plotting dim]
MANTLE.VHmeanMM                	= NaN;    	%[plotting dim]
MANTLE.VHmeanLM               	= NaN;    	%[plotting dim]

MANTLE.activeUpwelling          = NaN;  	%1 for active upwelling
MANTLE.activeDownwelling        = NaN;   	%1 for active downwelling
MANTLE.numActUpwelling          = NaN;    	%num
MANTLE.numActDownwelling        = NaN;   	%num
MANTLE.ActUpwellingVolume       = NaN;     	%[nd] or [m^2]or[m^3]
MANTLE.ActDownwellingVolume    	= NaN;     	%[nd] or [m^2]or[m^3]
MANTLE.ActUpwellingVolPerc      = NaN;   	%in [% of total upwelling]
MANTLE.ActDownwellingVolPerc  	= NaN;    	%in [% of total downwelling]

MANTLE.passiveUpwelling         = NaN;    	%1 for passive upwelling
MANTLE.passiveDownwelling       = NaN;   	%1 for passive downwelling
MANTLE.numPassUpwelling         = NaN;    	%num
MANTLE.numPassDownwelling       = NaN;    	%num
MANTLE.PassUpwellingVolume   	= NaN;    	%[nd] or [m^2]or[m^3]
MANTLE.PassDownwellingVolume 	= NaN;      %[nd] or [m^2]or[m^3]
MANTLE.PassUpwellingVolPerc  	= NaN;  	%in [% of total upwelling]
MANTLE.PassDownwellingVolPerc 	= NaN;      %in [% of total downwelling]

MANTLE.plumeHotNumberUM         = NaN;    	%num
MANTLE.plumeHotVHmaxUM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHminUM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHmeanUM         = NaN;    	%[plotting dim]
MANTLE.plumeHotNumberMM         = NaN;    	%num
MANTLE.plumeHotVHmaxMM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHminMM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHmeanMM         = NaN;    	%[plotting dim]
MANTLE.plumeHotNumberLM         = NaN;    	%num
MANTLE.plumeHotVHmaxLM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHminLM          = NaN;    	%[plotting dim]
MANTLE.plumeHotVHmeanLM         = NaN;    	%[plotting dim]

MANTLE.plumeColdNumberUM      	= NaN;    	%num
MANTLE.plumeColdVHmaxUM     	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHminUM      	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHmeanUM    	= NaN;    	%[plotting dim]
MANTLE.plumeColdNumberMM    	= NaN;    	%num
MANTLE.plumeColdVHmaxMM      	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHminMM     	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHmeanMM        = NaN;    	%[plotting dim]
MANTLE.plumeColdNumberLM     	= NaN;    	%num
MANTLE.plumeColdVHmaxLM      	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHminLM     	= NaN;    	%[plotting dim]
MANTLE.plumeColdVHmeanLM        = NaN;    	%[plotting dim]

MANTLE.contLocationX        	= NaN;    	%[plotting dim]
MANTLE.contLocationZ            = NaN;      %[plotting dim]

MANTLE.contNumberUM             = NaN;    	%num
MANTLE.contVHmaxUM              = NaN;    	%[plotting dim]
MANTLE.contVHminUM              = NaN;    	%[plotting dim]
MANTLE.contVHmeanUM             = NaN;    	%[plotting dim]

MANTLE.llsvpLocationX        	= NaN;    	%[plotting dim]
MANTLE.llsvpLocationZ           = NaN;      %[plotting dim]

MANTLE.llsvpNumberLM            = NaN;    	%num
MANTLE.llsvpVHmaxLM             = NaN;    	%[plotting dim]
MANTLE.llsvpVHminLM             = NaN;    	%[plotting dim]
MANTLE.llsvpVHmeanLM            = NaN;    	%[plotting dim]

if strcmp(GRID.Type,'yinyang')
    %add here......
end

%% INPUT ADJUSTMENTS
if strcmp(GRID.Type,'yinyang') %grid might NOT be flipped in yinyang
    warning('check if these values need to be flipped upside down!');
%     MD.upperDepthThresholdHot  = 1-MD.upperDepthThresholdHot;
%     MD.lowerDepthThresholdHot  = 1-MD.lowerDepthThresholdHot;
%     MD.upperDepthThresholdCold = 1-MD.upperDepthThresholdCold;
%     MD.lowerDepthThresholdCold = 1-MD.lowerDepthThresholdCold;
end

%% SETUP SWITCHES
if strcmp(GRID.Type,'yinyang')
    nb = 2;
else
    nb = 1;
end

%% READ INPUT DATA
%temperature
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Temperature';
DATA.FieldAbbreviation      = 'T';
DATA.StopExecutionIfNotFound = true;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
interiorIndexLow            = floor( size(GRID.Z_3D,3)*1/10 );  %bottom of the mantle
interiorIndexHigh           = ceil( size(GRID.Z_3D,3)*8/10 );   %top of the mantle
if strcmp(GRID.Type,'yinyang')
    T_3D        = PLOT.T_3Dyin; T_3Dyang   = PLOT.T_3Dyang;
    Tmin        = min([PLOT.T_3Dyin(:);PLOT.T_3Dyang(:)]);
    Tmax        = max([PLOT.T_3Dyin(:);PLOT.T_3Dyang(:)]);
    dummy       = [PLOT.T_3Dyin(:,:,interiorIndexLow:interiorIndexHigh);PLOT.T_3Dyang(:,:,interiorIndexLow:interiorIndexHigh)];
    Tmedian    	= median(dummy(:)); %handle with care: median does not consider effects of an air layer or refined/spherical grid...
else %all other grid types
    T_3D        = PLOT.T_3D;
    Tmin        = min(PLOT.T_3D(:));
    Tmax        = max(PLOT.T_3D(:));
    dummy       = PLOT.T_3D(:,:,interiorIndexLow:interiorIndexHigh);
    Tmedian  	= median(dummy(:)); %handle with care: median does not consider effects of an air layer or refined/spherical grid...
end
clearvars dummy
if DATA.NotFound
    MANTLE.DiagnosticsFailed	= true;
    warning off backtrace
    disp(' '); warning('Mantle diagnostics failed due to missing temperature data.'); disp(' ');
    warning on backtrace
    return
    
end

%continental composition
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Cont. crust';
DATA.FieldAbbreviation      = 'CC';
DATA.StopExecutionIfNotFound = false;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    CC_3D    	= PLOT.CC_3Dyin; CC_3Dyang   = PLOT.CC_3Dyang;
else %all other grid types
    CC_3D  	= PLOT.CC_3D;
end
clearvars dummy

%primordial composition
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Primordial';
DATA.FieldAbbreviation      = 'PRM';
DATA.StopExecutionIfNotFound = false;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    PRM_3D    	= PLOT.PRM_3Dyin; PRM_3Dyang   = PLOT.PRM_3Dyang;
else %all other grid types
    PRM_3D  	= PLOT.PRM_3D;
end
clearvars dummy

%velocity
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Velocity';
DATA.FieldAbbreviation      = 'V';
DATA.StopExecutionIfNotFound = true;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
%     VX_3D = PLOT.VX_3Dyin; VX_3Dyang = PLOT.VX_3Dyang;
%     VY_3D = PLOT.VY_3Dyin; VY_3Dyang = PLOT.VY_3Dyang;
    VH_3D = PLOT.VH_3Dyin; VH_3Dyang = PLOT.VH_3Dyang;
    VZ_3D = PLOT.VZ_3Dyin; VZ_3Dyang = PLOT.VZ_3Dyang;
    VZmax       = max(abs([VZ_3D(:);VZ_3Dyang(:)]));   %doesn't take variable gridspacing into account....
    VZupmax     = abs(max([VZ_3D(:);VZ_3Dyang(:)]));
    VZdownmax 	= abs(min([VZ_3D(:);VZ_3Dyang(:)]));
    
else %all other grid types
%     VX_3D       = PLOT.VX_3D;
%     VY_3D       = PLOT.VY_3D;
    VH_3D       = PLOT.VH_3D;
    VZ_3D       = PLOT.VZ_3D;
    VZmax       = max(abs(VZ_3D(:)));   %doesn't take variable gridspacing into account....
    VZupmax   	= abs(max(VZ_3D(:)));
    VZdownmax  	= abs(min(VZ_3D(:)));
end
if DATA.NotFound
    MANTLE.DiagnosticsFailed	= true;
    warning off backtrace
    disp(' '); warning('Mantle diagnostics failed due to missing velocity data.'); disp(' ');
    warning on backtrace
    return
    
end

%% DIMENSIONALISATION
VH_3D           = VH_3D*SETUP.Vscale;
if strcmp(GRID.Type,'yinyang')
    VH_3Dyang 	= VH_3Dyang*SETUP.Vscale;
end

%% EXTRACT FLOW-VELOCITY CHARACTERISTICS AT CERTAIN DEPTH LEVELS
DomainDepthMax          = max(GRID.Z_3Dp(:));
DomainDepthMin          = min(GRID.Z_3Dp(GRID.Z_3Dp>0));
DomainDepth             = abs( DomainDepthMax - DomainDepthMin );
%define depth levels
UMdepth                 = 0.7;              %values between 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<
MMdepth                 = 0.5;              %values between 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<
LMdepth                 = 0.1;              %values between 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<

depthLevelDiagnostics   = [UMdepth, MMdepth, LMdepth];  %values between 1 (surface) and 0 (bottom)

%derive actual depth values and indices
depthLevelDiagnosticsP  	= DomainDepthMax - depthLevelDiagnostics.*DomainDepth;  %plotting values
depthLevelDiagnosticsIdx    = zeros(size(depthLevelDiagnosticsP)); %array indices
for id=1:length(depthLevelDiagnosticsP)
    [~,depthLevelDiagnosticsIdx(1,id)]	= min(abs(GRID.Z_3Dp(1,1,:)-depthLevelDiagnosticsP(id)));  %index of closest value to mean LAB depth further up in the plate
end
UMdepthIdx          = depthLevelDiagnosticsIdx(1,1);
MMdepthIdx          = depthLevelDiagnosticsIdx(1,2);
LMdepthIdx          = depthLevelDiagnosticsIdx(1,3);
UMdepthP            = depthLevelDiagnosticsP(1,1);
MMdepthP            = depthLevelDiagnosticsP(1,2);
LMdepthP            = depthLevelDiagnosticsP(1,3);

%extracting data
for id=1:length(depthLevelDiagnosticsP)
    %horizontal rms values
    if strcmp(GRID.Type,'yinyang')
        VHmax(1,id)   	= max(max(abs( [VH_3D(:,:,depthLevelDiagnosticsIdx(1,id)),VH_3Dyang(:,:,depthLevelDiagnosticsIdx(1,id))] )));
        VHmean(1,id)   	= mean2(abs( [VH_3D(:,:,depthLevelDiagnosticsIdx(1,id)),VH_3Dyang(:,:,depthLevelDiagnosticsIdx(1,id))] ));
    else
        VHmax(1,id)   	= max(max(abs( VH_3D(:,:,depthLevelDiagnosticsIdx(1,id)) )));
        VHmean(1,id)   	= mean2(abs( VH_3D(:,:,depthLevelDiagnosticsIdx(1,id)) ));
    end
end

%writing it to output variables
MANTLE.VHmaxUM                 	= VHmax(1,1);    	%[plotting dim]
MANTLE.VHmaxMM                 	= VHmax(1,2);    	%[plotting dim]
MANTLE.VHmaxLM                 	= VHmax(1,3);    	%[plotting dim]
MANTLE.VHmeanUM               	= VHmean(1,1);    	%[plotting dim]
MANTLE.VHmeanMM                	= VHmean(1,2);    	%[plotting dim]
MANTLE.VHmeanLM               	= VHmean(1,3);    	%[plotting dim]

%% CALCULATE RESIDUAL TEMPERATURE FIELD
%initiation
Plumes3D = zeros(size(T_3D));
Tres3D = Plumes3D;

if MD.residualMode==1 %Standard horizontal residual
    SWITCH.customFieldTask  = 'Horizontal residual';
elseif MD.residualMode==2 %Horizontal-band residual
    SWITCH.customFieldTask  = 'Horizontal-band residual';
elseif MD.residualMode==3 %Global residual
    SWITCH.customFieldTask  = 'Global residual';
elseif MD.residualMode==4 %Regional residual
    SWITCH.customFieldTask  = 'Regional residual';
end
[VAR] = f_makeField(1,FILE,GRID,SETUP,SWITCH,PLOT,1,1);
%convert variables
if strcmp(GRID.Type,'yinyang')
    Tres3D              = VAR.var2d;
    Tres3Dyang          = VAR.var2d_yang;
    %ADD horizMean etc here for yin and yang ..................
    warning('from here on not implemented for yinyang grid...')
    return
    
    horizMeanTres       = VAR.meanTres; %might be regional values
    horizMin            = VAR.minTres;
    horizMax            = VAR.maxTres;
else
    if strcmp(GRID.Dim,'2-D')
        Tres3D(:,1,:) 	= VAR.var2d;
    else %strcmp(GRID.Dim,'3-D')
        Tres3D          = VAR.var2d;
    end
    horizMeanTres       = VAR.meanTres; %might be regional values
    horizMin            = VAR.minTres;
    horizMax            = VAR.maxTres;
end
clearvars VAR

%% FIND ACTIVE PLUMES
%initialising
Plumes3D = zeros(size(Tres3D));
thrHot = Plumes3D;
thrCold = Plumes3D;
if strcmp(GRID.Type,'yinyang'); Plumes3Dyang = Plumes3D; end

%start plume diagnostics
thrHot                    	= horizMeanTres + PLOT.pHot.*(horizMax-horizMeanTres);	%after Labrosse (EPSL,2002)
thrCold                 	= horizMeanTres + PLOT.pCold.*(horizMin-horizMeanTres);	%after Labrosse (EPSL,2002)

Plumes3D(Tres3D>thrHot)   	= 1;            	%hot plumes
Plumes3D(Tres3D<thrCold)  	= -1;              	%cold plumes

if strcmp(GRID.Type,'yinyang')
    Plumes3Dyang(Tres3Dyang>thrHot_yang)    	= 1;  	%hot plumes
    Plumes3Dyang(Tres3Dyang<thrCold_yang)  	= -1;  	%cold plumes
    %XXXXX thrHot thrCold
    
end
clearvars horizMeanT horizMeanTres horizMax horizMin thrHot thrCold dummy dummyX dummyY

if logical(0)
    figure(3)
    hold on
    x(:,:) = GRID.X_3Dp(:,1,:);
    z(:,:) = GRID.Z_3Dp(:,1,:);
    p(:,:) = Plumes3D(:,1,:);
    phot = zeros(size(p));
    phot(p==1) = 1;
    phot = bwmorph(phot,'thicken',1);
    phot = bwmorph(phot,'bridge');
    
    pcold = zeros(size(p));
    pcold(p==-1) = 1;
    %    pcold = bwmorph(pcold,'thicken',2);
    pcold = bwmorph(pcold,'bridge');
    
    p2 = phot-pcold;
    contourf(x,z,p2)
    colorbar
    axis ij
    figure(1)
end

%% FIND UP- AND DOWNWELLINGS
thrVZ                           = 1/100 *VZmax;       	%define threshold for up/down welling <<<<<<<<<<<<<<<<<<
thrVZup                         = 1/100 *VZupmax;       %define threshold for upwelling <<<<<<<<<<<<<<<<<<
thrVZdown                       = 1/100 *VZdownmax;   	%define threshold for downwelling <<<<<<<<<<<<<<<<<<
upWelling                       = zeros(size(VZ_3D));  	%nothing
downWelling                     = zeros(size(VZ_3D));   %nothing
upWellingAbsolute             	= zeros(size(VZ_3D));  	%nothing
downWellingAbsolute            	= zeros(size(VZ_3D));   %nothing
upWelling(VZ_3D>thrVZup)        = 1;                  	%upwelling (exceeding threshold)
downWelling(VZ_3D<-thrVZdown)	= 1;                    %downwelling (exceeding threshold)
upWellingAbsolute(VZ_3D>0)      = 1;                  	%upwelling (absolute)
downWellingAbsolute(VZ_3D<0)   	= 1;                  	%downwelling (absolute)
if strcmp(GRID.Type,'yinyang')
    upWelling_yang                          = zeros(size(VZ_3Dyang)); 	%nothing
    downWelling_yang                        = zeros(size(VZ_3Dyang)); 	%nothing
    upWellingAbsolute_yang              	= zeros(size(VZ_3Dyang)); 	%nothing
    downWellingAbsolute_yang             	= zeros(size(VZ_3Dyang)); 	%nothing
    upWelling_yang(VZ_3Dyang>thrVZup)      = 1;                        %upwelling (exceeding threshold)
    downWelling_yang(VZ_3Dyang<-thrVZdown) = 1;                        %downwelling (exceeding threshold)
    upWellingAbsolute_yang(VZ_3Dyang>0)   	= 1;                        %upwelling (absolute)
    downWellingAbsolute_yang(VZ_3Dyang<0) 	= 1;                        %downwelling (absolute)
end

%% FIND Continents
%critical variables
thrCont                     = 0.5;  %composional threshold for continental material (0<1) <<<<<<<<<<<<<<<<<<<<<
%initialising
CONT3D = zeros(size(CC_3D));
if strcmp(GRID.Type,'yinyang'); CONT3Dyang = CONT3D; end
%start diagnostics
CONT3D(CC_3D>thrCont)   	= 1;            	%continental material

%% FIND LLSVPs
%critical variables
thrPrm                      = 0.5;  %composional threshold for primordial material (0<1) <<<<<<<<<<<<<<<<<<<<<
%initialising
LLSVP3D = zeros(size(PRM_3D));
if strcmp(GRID.Type,'yinyang'); LLSVP3Dyang = LLSVP3D; end
%start diagnostics
LLSVP3D(PRM_3D>thrPrm)   	= 1;            	%llsvp material

%% ACCOUNT FOR PERIODIC BOUNDARIES
%check for plumes crossing periodic boundaries:
add_factor          = MD.addFactor;      %percent of total x-extent that is added at x=1 and x=end
%add some slices at x=1 and x=end
PlumesOrig          = Plumes3D;
extentXgrid         = false;
extentYgrid         = false;
nx_orig             = size(PlumesOrig,1);
ny_orig             = size(PlumesOrig,2);
nx_add              = nx_orig*add_factor;
ny_add              = ny_orig*add_factor;
if strcmp(GRID.Type,'yinyang')
    extentYgrid     = true;
elseif strcmp(GRID.Type,'Cartesian')
    if strcmp(GRID.Dim,'2-D')
        if size(Plumes3D,1)>1 && size(Plumes3D,2)==1 && strcmp(GRID.BCxx,'permeable'); extentXgrid = true; end
        if size(Plumes3D,1)==1 && size(Plumes3D,2)>1 && strcmp(GRID.BCyy,'permeable'); extentYgrid = true; end
    elseif strcmp(GRID.Dim,'3-D')
        if strcmp(GRID.BCxx,'permeable'); extentXgrid = true; end
        if strcmp(GRID.BCyy,'permeable'); extentYgrid = true; end
        warning(['MANTLE DIAGNOSTICS: Grid type ',GRID.Type,' in 3-D not tested yet!'])
    end
elseif strcmp(GRID.Type,'spherical2D')
    if size(Plumes3D,1)>1 && size(Plumes3D,2)==1 && strcmp(GRID.BCxx,'permeable'); extentXgrid = true; end
    if size(Plumes3D,1)==1 && size(Plumes3D,2)>1 && strcmp(GRID.BCyy,'permeable'); extentYgrid = true; end
else
    error(['Grid type ',GRID.Type,' not found!'])
end
if extentXgrid
    X_3DpExt        = [GRID.X_3Dp((end-nx_add+1):end,:,:); GRID.X_3Dp; GRID.X_3Dp(1:nx_add,:,:)];
    Y_3DpExt        = [GRID.Y_3Dp((end-nx_add+1):end,:,:); GRID.Y_3Dp; GRID.Y_3Dp(1:nx_add,:,:)];
    Z_3DpExt        = [GRID.Z_3Dp((end-nx_add+1):end,:,:); GRID.Z_3Dp; GRID.Z_3Dp(1:nx_add,:,:)];
    Plumes3D        = [Plumes3D((end-nx_add+1):end,:,:); Plumes3D; Plumes3D(1:nx_add,:,:)];
    CONT3D          = [CONT3D((end-nx_add+1):end,:,:); CONT3D; CONT3D(1:nx_add,:,:)];
    LLSVP3D         = [LLSVP3D((end-nx_add+1):end,:,:); LLSVP3D; LLSVP3D(1:nx_add,:,:)];
    upWelling       = [upWelling((end-nx_add+1):end,:,:); upWelling; upWelling(1:nx_add,:,:)];
    downWelling     = [downWelling((end-nx_add+1):end,:,:); downWelling; downWelling(1:nx_add,:,:)];
    upWellingAbsolute       = [upWellingAbsolute((end-nx_add+1):end,:,:); upWellingAbsolute; upWellingAbsolute(1:nx_add,:,:)];
    downWellingAbsolute     = [downWellingAbsolute((end-nx_add+1):end,:,:); downWellingAbsolute; downWellingAbsolute(1:nx_add,:,:)];
    if strcmp(GRID.Type,'yinyang')
        Plumes3Dyang 	= [Plumes3Dyang((end-nx_add+1):end,:,:); Plumes3Dyang; Plumes3Dyang(1:nx_add,:,:)];
        LLSVP3Dyang 	= [LLSVP3Dyang((end-nx_add+1):end,:,:); LLSVP3Dyang; LLSVP3Dyang(1:nx_add,:,:)];
        upWelling_yang 	= [upWelling_yang((end-nx_add+1):end,:,:); upWelling_yang; upWelling_yang(1:nx_add,:,:)];
        downWelling_yang= [downWelling_yang((end-nx_add+1):end,:,:); downWelling_yang; downWelling_yang(1:nx_add,:,:)];
        upWellingAbsolute_yang 	= [upWellingAbsolute_yang((end-nx_add+1):end,:,:); upWellingAbsolute_yang; upWellingAbsolute_yang(1:nx_add,:,:)];
        downWellingAbsolute_yang= [downWellingAbsolute_yang((end-nx_add+1):end,:,:); downWellingAbsolute_yang; downWellingAbsolute_yang(1:nx_add,:,:)];
    end
else
    X_3DpExt        = GRID.X_3Dp;
    Y_3DpExt        = GRID.X_3Dp;
    Z_3DpExt        = GRID.X_3Dp;
end
if extentYgrid
    X_3DpExt     	= [X_3DpExt(:,(end-ny_add+1):end,:), X_3DpExt, X_3DpExt(:,1:ny_add,:)];
    Y_3DpExt     	= [Y_3DpExt(:,(end-ny_add+1):end,:), Y_3DpExt, Y_3DpExt(:,1:ny_add,:)];
    Z_3DpExt     	= [Z_3DpExt(:,(end-ny_add+1):end,:), Z_3DpExt, Z_3DpExt(:,1:ny_add,:)];
    Plumes3D     	= [Plumes3D(:,(end-ny_add+1):end,:), Plumes3D, Plumes3D(:,1:ny_add,:)];
    CONT3D          = [CONT3D(:,(end-ny_add+1):end,:), CONT3D, CONT3D(:,1:ny_add,:)];
    LLSVP3D     	= [LLSVP3D(:,(end-ny_add+1):end,:), LLSVP3D, LLSVP3D(:,1:ny_add,:)];
    upWelling       = [upWelling(:,(end-ny_add+1):end,:), upWelling, upWelling(:,1:ny_add,:)];
    downWelling     = [downWelling(:,(end-ny_add+1):end,:), downWelling, downWelling(:,1:ny_add,:)];
    upWellingAbsolute       = [upWellingAbsolute(:,(end-ny_add+1):end,:), upWellingAbsolute, upWellingAbsolute(:,1:ny_add,:)];
    downWellingAbsolute     = [downWellingAbsolute(:,(end-ny_add+1):end,:), downWellingAbsolute, downWellingAbsolute(:,1:ny_add,:)];
    if strcmp(GRID.Type,'yinyang')
        Plumes3Dyang 	= [Plumes3Dyang(:,(end-ny_add+1):end,:), Plumes3Dyang, Plumes3Dyang(:,1:ny_add,:)];
        LLSVP3Dyang 	= [LLSVP3Dyang(:,(end-ny_add+1):end,:), LLSVP3Dyang, LLSVP3Dyang(:,1:ny_add,:)];
        upWelling_yang	= [upWelling_yang(:,(end-ny_add+1):end,:), upWelling_yang, upWelling_yang(:,1:ny_add,:)];
        downWelling_yang= [downWelling_yang(:,(end-ny_add+1):end,:), downWelling_yang, downWelling_yang(:,1:ny_add,:)];
        upWellingAbsolute_yang	= [upWellingAbsolute_yang(:,(end-ny_add+1):end,:), upWellingAbsolute_yang, upWellingAbsolute_yang(:,1:ny_add,:)];
        downWellingAbsolute_yang= [downWellingAbsolute_yang(:,(end-ny_add+1):end,:), downWellingAbsolute_yang, downWellingAbsolute_yang(:,1:ny_add,:)];
    end
end

for ib=1:nb
    Plumes3Dh = zeros(size(Plumes3D)); Plumes3Dc = Plumes3Dh;
    CONT3Dx = zeros(size(Plumes3D));
    LLSVP3Dx = zeros(size(Plumes3D));
    if ib==1
        Plumes3Dh(Plumes3D==1)  = 1;
        Plumes3Dc(Plumes3D==-1) = 1;
        CONT3Dx                 = CONT3D;
        LLSVP3Dx                = LLSVP3D;
    elseif ib==2
        Plumes3Dh(Plumes3Dyang==1)  = 1;
        Plumes3Dc(Plumes3Dyang==-1) = 1;
        CONT3Dx                 = CONT3Dyang;
        LLSVP3Dx                = LLSVP3Dyang;
    end
    %% OPTIMISE DATA I
    if logical(1)
        if strcmp(GRID.Dim,'2-D')
            %Plumes
            dummyh(:,:)         = Plumes3Dh(:,1,:); %make 2-dimensional
            dummyc(:,:)         = Plumes3Dc(:,1,:); %make 2-dimensional
            %thicken
%             dummyh              = bwmorph(dummyh,'thicken',1);
%             dummyc              = bwmorph(dummyc,'thicken',1);
            %bridge
            dummyh              = bwmorph(dummyh,'bridge');
            dummyc              = bwmorph(dummyc,'bridge');
            Plumes3Dh(:,1,:)    = dummyh(:,:);
            Plumes3Dc(:,1,:)    = dummyc(:,:);
            clearvars dummyh dummyc
            
            for iF=1:2
                if iF==1 %Continents
                    dummy(:,:)          = CONT3Dx(:,1,:); %make 2-dimensional
                elseif iF==2 %LLSVPs
                    dummy(:,:)          = LLSVP3Dx(:,1,:); %make 2-dimensional
                else
                    error('check here')
                end
                %thicken
                %             dummy              = bwmorph(dummy,'thicken',1);
                %bridge
                dummy               = bwmorph(dummy,'bridge');
                %fill isolated pixels
                dummy               = bwmorph(dummy,'fill');
                %fill holes
                dummy               = imfill(dummy,'holes');   %<<<<<<<<<if it is filling large holes, the volume diagnostics will be inaccurate!
                if iF==1 %Continents
                    CONT3Dx(:,1,:)    	= dummy(:,:);
                elseif iF==2 %LLSVPs
                    LLSVP3Dx(:,1,:)    	= dummy(:,:);
                end
                clearvars dummy
            end
            
        elseif strcmp(GRID.Dim,'3-D')
            if SWITCH.Verbose; warning('Optimising plume data in 3-D not possible yet'); end
        end
    end
    
    if logical(0)
        figure(3)
        x(:,:) = X_3DpExt(:,1,:);
        z(:,:) = Z_3DpExt(:,1,:);
        p(:,:) = Plumes3Dh(:,1,:);
        contourf(x,z,p)
        colorbar
        axis ij
        figure(1)
    end
   
    %% CHECK CONNECTIVITY FOR AREA SIZES
    CCph            = bwconncomp(Plumes3Dh);     	%HOT PLUMES
    numPixelsPh     = cellfun(@numel,CCph.PixelIdxList);    %number of connected pixels in each connected area
    CCpc            = bwconncomp(Plumes3Dc);     	%COLD PLUMES
    numPixelsPc     = cellfun(@numel,CCpc.PixelIdxList);    %number of connected pixels in each connected area
    CCuw            = bwconncomp(upWelling);        %UPWELLINGS
    numPixelsUW     = cellfun(@numel,CCuw.PixelIdxList);    %number of connected pixels in each connected area
    CCdw            = bwconncomp(downWelling);      %DOWNWELLINGS
    numPixelsDW     = cellfun(@numel,CCdw.PixelIdxList);    %number of connected pixels in each connected area
    CCcont          = bwconncomp(CONT3Dx);       	%Continents
    numPixelsCont	= cellfun(@numel,CCcont.PixelIdxList); %number of connected pixels in each connected area
    CCllsvp         = bwconncomp(LLSVP3Dx);       	%LLSVPs
    numPixelsLLSVP	= cellfun(@numel,CCllsvp.PixelIdxList); %number of connected pixels in each connected area
    
    %% REMOVING SMALL ANOMALIES
    removeSmallAnomalies    = logical(1);
    plumesBhot  = Plumes3Dh;  	%hot plumes
    plumesBcold = Plumes3Dc;  	%cold plumes
    if removeSmallAnomalies
        %Plumes
        % *********************** small area criterion <<<<<<<<<<<<<<<<
        % [size_largest,ind_largest]	= max(numPixels);   %find single most largest area
        sizeThreshold           = 1/5* size(Plumes3Dh,3);   %find areas larger than (1/5*nz) pixels
        if strcmp(GRID.Dim,'3-D')
            sizeThreshold       = sizeThreshold *4/3*sqrt(sizeThreshold/pi); %find areas larger than x pixels
        end
        indLargePh          	= find(numPixelsPh<sizeThreshold);
        indLargePc          	= find(numPixelsPc<sizeThreshold);
        % ***********************
        for i_ind=1:size(indLargePh,2)
            plumesBhot(CCph.PixelIdxList{indLargePh(i_ind)})    = 0;
        end
        for i_ind=1:size(indLargePc,2)
            plumesBcold(CCpc.PixelIdxList{indLargePc(i_ind)})   = 0;
        end
        
        for iF=1:2
            if iF==1 %Continents
                % *********************** small area criterion <<<<<<<<<<<<<<<<
                % [size_largest,ind_largest]	= max(numPixels);   %find single most largest area
                sizeThreshold           = 10* size(CONT3Dx,3);   %find areas larger than (1/5*nz) pixels
                if strcmp(GRID.Dim,'3-D')
                    sizeThreshold       = sizeThreshold *4/3*sqrt(sizeThreshold/pi); %find areas larger than x pixels
                end
                indLargeCont          	= find(numPixelsCont<sizeThreshold);
                % ***********************
                for i_ind=1:size(indLargeCont,2)
                    CONT3Dx(CCcont.PixelIdxList{indLargeCont(i_ind)}) = 0;
                end
            elseif iF==2 %LLSVPs
                % *********************** small area criterion <<<<<<<<<<<<<<<<
                % [size_largest,ind_largest]	= max(numPixels);   %find single most largest area
                sizeThreshold           = 10* size(LLSVP3Dx,3);   %find areas larger than (1/5*nz) pixels
                if strcmp(GRID.Dim,'3-D')
                    sizeThreshold       = sizeThreshold *4/3*sqrt(sizeThreshold/pi); %find areas larger than x pixels
                end
                indLargeLlsvp          	= find(numPixelsLLSVP<sizeThreshold);
                % ***********************
                for i_ind=1:size(indLargeLlsvp,2)
                    LLSVP3Dx(CCllsvp.PixelIdxList{indLargeLlsvp(i_ind)}) = 0;
                end
            end
        end
    end
    
    %% SETTING UP HOT AND COLD ANOMALIES after criterion
    if ib==1
        plumesHot             	= plumesBhot;   %hot plumes
        plumesCold              = plumesBcold;  %cold plumes
        CONT3D                  = CONT3Dx;      %continents
        LLSVP3D                 = LLSVP3Dx;     %llsvp
    elseif ib==2
        plumesHot_yang      	= plumesBhot;   %hot plumes
        plumesCold_yang     	= plumesBcold;	%cold plumes
        CONT3Dyang         	= CONT3Dx;      %continents
        LLSVP3Dyang           	= LLSVP3Dx;     %llsvp
    end
    clearvars plumesBhot plumesBcold
end

if logical(0)
    figure(3)
    x(:,:) = X_3DpExt(:,1,:);
    z(:,:) = Z_3DpExt(:,1,:);
    p(:,:) = plumesHot(:,1,:);
%     p(:,:) = LLSVP3D(:,1,:);
    contourf(x,z,p)
    colorbar
    axis ij
    figure(1)
end

%% CONNECTIVITY CHECK WITH CORRESPONDING BOUNDARY LAYER
% ...AND CHECK FOR PLUME EXTENSION THROUGHOUT UPPER-AND LOWER-DEPTH THRESHOLDS
% This uses PLOT.Z_3Dp, which is always flipped and actual depthCCPhot   = bwconncomp(plumesHot);
CCPhot      = bwconncomp(plumesHot);
CCPcold     = bwconncomp(plumesCold);
if strcmp(GRID.Type,'yinyang')
    CCPhot_yang   = bwconncomp(plumesHot_yang);
    CCPcold_yang  = bwconncomp(plumesCold_yang);
end
Ztop    = 0;
Zbot    = max(Z_3DpExt(:));
Zdiff   = Zbot-Ztop;
if Zdiff<=0; error('Zdiff should be >0'); end
for ib=1:nb
    for plume_kind=1:2
        if plume_kind==1  %HOT
            if ib==1
                CC_kind = CCPhot; 
            else
                CC_kind = CCPhot_yang; 
            end
            ZtopThreshold  	= Ztop + Zdiff*MD.upperDepthThresholdHot;
            ZbotThreshold  	= Ztop + Zdiff*MD.lowerDepthThresholdHot;
        elseif plume_kind==2  %COLD
            if ib==1
                CC_kind = CCPcold;
            elseif ib==2
                CC_kind = CCPcold_yang;
            end
            ZtopThreshold  	= Ztop + Zdiff*MD.upperDepthThresholdCold;
            ZbotThreshold 	= Ztop + Zdiff*MD.lowerDepthThresholdCold;
        end
        for i_area=1:size(CC_kind.PixelIdxList,2)  %loop all connected areas
            surf_connection = false;
            bot_connection = false;
            for ip=1:size(CC_kind.PixelIdxList{1,i_area},1)  %loop all pixels of an area
                i_pixel = CC_kind.PixelIdxList{1,i_area}(ip,1);
                %[x,y,z] = ind2sub(size(PLUMES),i_pixel);
                if Z_3DpExt(i_pixel)<=ZtopThreshold %is connected to the top (0)
                    surf_connection = true;
                end
                if Z_3DpExt(i_pixel)>=ZbotThreshold %is connected to the bottom (1)
                    bot_connection = true;
                end
                if surf_connection && bot_connection
                    break   %this is a plume exceeding top to bottom levels! check next.
                end
            end
            
            if surf_connection && bot_connection %of defined depth thresholds
                %keep area
            else %remove area: this is no plume
                if plume_kind==1  %HOT
                    if ib==1
                        plumesHot(CC_kind.PixelIdxList{i_area})         = 0;
                    elseif ib==2
                        plumesHot_yang(CC_kind.PixelIdxList{i_area})    = 0;
                    end
                elseif plume_kind==2  %COLD
                    if ib==1
                        plumesCold(CC_kind.PixelIdxList{i_area})        = 0;
                    elseif ib==2
                        plumesCold_yang(CC_kind.PixelIdxList{i_area})   = 0;
                    end
                end
            end
        end
    end
end


%% GET CENTROIDS FOR EXTRACTED ANOMALIES
if ~strcmp(GRID.Type,'yinyang') && ~strcmp(GRID.Dim,'3-D')
    for iF=1:2
        if iF==1 %continents
            CCanom   	= bwconncomp(CONT3D);
            % CCllsvp2_yang   	= bwconncomp(CONT3Dyang);
        elseif iF==2 %llsvp
            CCanom   	= bwconncomp(LLSVP3D);
            % CCllsvp2_yang   	= bwconncomp(LLSVP3Dyang);
        end
        CentroidAnomIdx = zeros(size(CCanom.PixelIdxList,2),2);
        for i_area=1:size(CCanom.PixelIdxList,2)  %loop all connected areas
            areaCurrent(:,:)    	= false(size(LLSVP3D));
            areaCurrent(CCanom.PixelIdxList{i_area})    = true;
            %     areaCurrent     = bwareafilt(areaCurrent, 1);
            
            props                   = regionprops(areaCurrent, 'Centroid');
            dummy(i_area,:)         = round(props.Centroid); %find closest index
            CentroidAnomIdx(i_area,1)	= min(max(1,dummy(i_area,2)),size(LLSVP3D,1)); %need to flip x and z values, and limit to max/min indices
            CentroidAnomIdx(i_area,2)	= min(max(1,dummy(i_area,1)),size(LLSVP3D,3));
            
            % binaryImage = bwareafilt(binaryImage, 1); % Extract largest blob ONLY.
            % props = regionprops(binaryImage, 'Centroid');
            % centroid = props.Centroid;
        end
        %convert to actual plot position
        dummy2(1,:)        	= Z_3DpExt(1,1,CentroidAnomIdx(:,2));
        if iF==1 %continents
            CentroidCONT       = [X_3DpExt(CentroidAnomIdx(:,1),1,1), dummy2'];
            %remove dublicated entries (due to periodic side boundaries)
            CentroidCONT       = unique(CentroidCONT,'rows');
        elseif iF==2 %llsvps
            CentroidLLSVP       = [X_3DpExt(CentroidAnomIdx(:,1),1,1), dummy2'];
            %remove dublicated entries (due to periodic side boundaries)
            CentroidLLSVP       = unique(CentroidLLSVP,'rows');
        end
        clearvars dummy dummy2 areaCurrent

        if logical(0)
            figure(3)
            x(:,:) = X_3DpExt(:,1,:);
            z(:,:) = Z_3DpExt(:,1,:);
            if iF==1 %continents
                p(:,:) = CONT3D(:,1,:);
                centroidDummy = CentroidCONT;
            else
                p(:,:) = LLSVP3D(:,1,:);
                centroidDummy = CentroidLLSVP;
            end
            contourf(x,z,p)
            hold on
            for iiii=1:size(CentroidAnomIdx,1)
                plot(X_3DpExt(CentroidAnomIdx(iiii,1),1,1),Z_3DpExt(1,1,CentroidAnomIdx(iiii,2)),'r+')
            end
            plot(centroidDummy(:,1),centroidDummy(:,2),'gx')
            colorbar
            axis ij
            figure(1)
        end
    end
else
    
    %add for yinyang!
    
end


%% REMOVING GRID EXTENSION (BACK TO ORIGINAL SIZE)
if extentXgrid
    Plumes3D            = Plumes3D((nx_add+1):end-nx_add,:,:);
    plumesHot           = plumesHot((nx_add+1):end-nx_add,:,:);
    plumesCold          = plumesCold((nx_add+1):end-nx_add,:,:);
    upWelling           = upWelling((nx_add+1):end-nx_add,:,:);
    downWelling         = downWelling((nx_add+1):end-nx_add,:,:);
    upWellingAbsolute  	= upWellingAbsolute((nx_add+1):end-nx_add,:,:);
    downWellingAbsolute	= downWellingAbsolute((nx_add+1):end-nx_add,:,:);
    CONT3D              = CONT3D((nx_add+1):end-nx_add,:,:);
    LLSVP3D             = LLSVP3D((nx_add+1):end-nx_add,:,:);
    if strcmp(GRID.Type,'yinyang')
        Plumes3Dyang           = Plumes3Dyang((nx_add+1):end-nx_add,:,:);
        plumesHot_yang          = plumesHot_yang((nx_add+1):end-nx_add,:,:);
        plumesCold_yang         = plumesCold_yang((nx_add+1):end-nx_add,:,:);
        upWelling_yang          = upWelling_yang((nx_add+1):end-nx_add,:,:);
        downWelling_yang        = downWelling_yang((nx_add+1):end-nx_add,:,:);
        upWellingAbsolute_yang	= upWellingAbsolute_yang((nx_add+1):end-nx_add,:,:);
        downWellingAbsolute_yang= downWellingAbsolute_yang((nx_add+1):end-nx_add,:,:);
        CONT3Dyang             = CONT3Dyang((nx_add+1):end-nx_add,:,:);
        LLSVP3Dyang            = LLSVP3Dyang((nx_add+1):end-nx_add,:,:);
    end
end
if extentYgrid
    Plumes3D            = Plumes3D(:,(ny_add+1):end-ny_add,:);
    plumesHot           = plumesHot(:,(ny_add+1):end-ny_add,:);
    plumesCold          = plumesCold(:,(ny_add+1):end-ny_add,:);
    upWelling           = upWelling(:,(ny_add+1):end-ny_add,:);
    downWelling         = downWelling(:,(ny_add+1):end-ny_add,:);
    upWellingAbsolute  	= upWellingAbsolute(:,(ny_add+1):end-ny_add,:);
    downWellingAbsolute	= downWellingAbsolute(:,(ny_add+1):end-ny_add,:);
    CONT3D              = CONT3D(:,(ny_add+1):end-ny_add,:);
    LLSVP3D             = LLSVP3D(:,(ny_add+1):end-ny_add,:);
    if strcmp(GRID.Type,'yinyang')
        Plumes3Dyang           = Plumes3Dyang(:,(ny_add+1):end-ny_add,:);
        plumesHot_yang          = plumesHot_yang(:,(ny_add+1):end-ny_add,:);
        plumesCold_yang         = plumesCold_yang(:,(ny_add+1):end-ny_add,:);
        upWelling_yang          = upWelling_yang(:,(ny_add+1):end-ny_add,:);
        downWelling_yang        = downWelling_yang(:,(ny_add+1):end-ny_add,:);
        upWellingAbsolute_yang 	= upWellingAbsolute_yang(:,(ny_add+1):end-ny_add,:);
        downWellingAbsolute_yang= downWellingAbsolute_yang(:,(ny_add+1):end-ny_add,:);
        CONT3Dyang             = CONT3Dyang(:,(ny_add+1):end-ny_add,:);
        LLSVP3Dyang            = LLSVP3Dyang(:,(ny_add+1):end-ny_add,:);
    end
end
clearvars X_3DpExt Y_3DpExt Z_3DpExt





%% DESCRIMINATE BETWEEN ACTIVE AND PASSIVE UP-/DOWNWELLING
passiveUpwelling    = upWelling-plumesHot;      passiveUpwelling(passiveUpwelling<0) = 0;
passiveDownwelling 	= downWelling-plumesCold;   passiveDownwelling(passiveDownwelling<0) = 0;
activeUpwelling     = zeros(size(upWelling));   activeUpwelling(upWellingAbsolute & plumesHot) = 1;
activeDownwelling	= zeros(size(downWelling)); activeDownwelling(downWellingAbsolute & plumesCold) = 1;
if strcmp(GRID.Type,'yinyang')
    passiveUpwelling_yang 	= upWelling_yang-plumesHot_yang; 	passiveUpwelling_yang(passiveUpwelling_yang<0) = 0;
    passiveDownwelling_yang	= downWelling_yang-plumesCold_yang;	passiveDownwelling_yang(passiveDownwelling_yang<0) = 0;
    activeUpwelling_yang    = zeros(size(upWelling_yang));      activeUpwelling_yang(upWellingAbsolute_yang & plumesHot_yang) = 1;
    activeDownwelling_yang	= zeros(size(downWelling_yang));    activeDownwelling_yang(downWellingAbsolute_yang & plumesCold_yang) = 1;
end

%% PLUME DEFINITION
if MD.plumeDefinition==1 %defined solely by temperature
%     plumesHot       = plumesHot;
%     plumesCold      = plumesCold;
elseif MD.plumeDefinition==2 %defined by temperature and vertical velocity
    plumesHot       = activeUpwelling;
    plumesCold      = activeDownwelling;
end


if ~strcmp(GRID.Type,'yinyang')
    %% DETAILS OF FOUND ANOMALIES I: number and lateral mobility at certain depths
    %derive horizontal slices
    depthLevelCONTslices    = 0.99;                         %values between 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<<<<
    depthLevelLLSVPslices   = 0.01;                         %values between 1 (surface) and 0 (bottom) <<<<<<<<<<<<<<<<<<<<<<
    depthLevelSlicesIdx     = [UMdepthIdx,MMdepthIdx,LMdepthIdx];
    depthLevelSlicesPL      = [UMdepthP,MMdepthP,LMdepthP];
    
    %PLUMES
    %derive horizontal plume slices
    plumesHotSlice      = plumesHot(:,:,depthLevelSlicesIdx);
    plumesColdSlice   	= plumesCold(:,:,depthLevelSlicesIdx);
    plumesHotSlice      = plumesHot(:,:,[UMdepthIdx,MMdepthIdx,LMdepthIdx]);
    plumesColdSlice   	= plumesCold(:,:,[UMdepthIdx,MMdepthIdx,LMdepthIdx]);
    
    %define point sources for plumes
    numPlumesHotSlice = zeros(1,length(depthLevelSlicesPL));
    for id=1:length(depthLevelSlicesPL)
        %     CC11        = bwconncomp(plumesHotSlice(:,:,id))
        plumesHotSlice(:,:,id)      = bwmorph(plumesHotSlice(:,:,id),'shrink',Inf);
        plumesColdSlice(:,:,id) 	= bwmorph(plumesColdSlice(:,:,id),'shrink',Inf);
        try
            numPlumesHotSlice(id)    	= sum(plumesHotSlice(:,:,id),'all');
            numPlumesColdSlice(id)    	= sum(plumesColdSlice(:,:,id),'all');
        catch
            numPlumesHotSlice(id)    	= sum(sum(plumesHotSlice(:,:,id)));
            numPlumesColdSlice(id)    	= sum(sum(plumesColdSlice(:,:,id)));
        end
    end
    
    %create hot plume data arrays
    maxNumPlumesHotPerSlice         = max(numPlumesHotSlice(:));
    dummy               = zeros(size(plumesHotSlice,1),size(plumesHotSlice,2))*NaN;
    vhFullSlice         = dummy;
    vhPlumeHotSlices    = zeros(length(depthLevelSlicesPL),maxNumPlumesHotPerSlice)*NaN; %[slices,plumes]
    xPlumeHotSlices     = vhPlumeHotSlices;
    yPlumeHotSlices     = vhPlumeHotSlices;
    zPlumeHotSlices     = vhPlumeHotSlices;
    for id=1:length(depthLevelSlicesPL)
        %extract horizontal position
        if max(plumesHotSlice(:,:,id))>0
            dummyNanoms = size(dummy(plumesHotSlice(:,:,id)==1),1);
            dummy(:,:)                          = GRID.X_3Dp(:,:,depthLevelSlicesIdx(id));
            xPlumeHotSlices(id,1:dummyNanoms) 	= dummy(plumesHotSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Y_3Dp(:,:,depthLevelSlicesIdx(id));
            yPlumeHotSlices(id,1:dummyNanoms) 	= dummy(plumesHotSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Z_3Dp(:,:,depthLevelSlicesIdx(id));
            zPlumeHotSlices(id,1:dummyNanoms) 	= dummy(plumesHotSlice(:,:,id)==1);
            %extract horizontal velocity (i.e., migration magnitude) of these points in the different slices
            vhFullSlice(:,:)                    = VH_3D(:,:,depthLevelSlicesIdx(id));
            vhPlumeHotSlices(id,1:dummyNanoms)	= vhFullSlice(plumesHotSlice(:,:,id)==1);
        end
    end
    clearvars dummy dummyNanoms
    
    %create cold plume data arrays
    maxNumPlumesColdPerSlice         = max(numPlumesColdSlice(:));
    dummy               = zeros(size(plumesColdSlice,1),size(plumesColdSlice,2))*NaN;
    vhFullSlice         = dummy;
    vhPlumeColdSlices   = zeros(length(depthLevelSlicesPL),maxNumPlumesColdPerSlice)*NaN; %[slices,plumes]
    xPlumeColdSlices    = vhPlumeColdSlices;
    yPlumeColdSlices    = vhPlumeColdSlices;
    zPlumeColdSlices    = vhPlumeColdSlices;
    for id=1:length(depthLevelSlicesPL)
        %extract horizontal position
        if max(plumesColdSlice(:,:,id))>0
            dummyNanoms = size(dummy(plumesColdSlice(:,:,id)==1),1);
            dummy(:,:)                          = GRID.X_3Dp(:,:,depthLevelSlicesIdx(id));
            xPlumeColdSlices(id,1:dummyNanoms)	= dummy(plumesColdSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Y_3Dp(:,:,depthLevelSlicesIdx(id));
            yPlumeColdSlices(id,1:dummyNanoms)	= dummy(plumesColdSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Z_3Dp(:,:,depthLevelSlicesIdx(id));
            zPlumeColdSlices(id,1:dummyNanoms)	= dummy(plumesColdSlice(:,:,id)==1);
            %extract horizontal velocity (i.e., migration magnitude) of these points in the different slices
            vhFullSlice(:,:)                    = VH_3D(:,:,depthLevelSlicesIdx(id));
            vhPlumeColdSlices(id,1:dummyNanoms)= vhFullSlice(plumesColdSlice(:,:,id)==1);
        end
    end
    clearvars dummy dummyNanoms
    
    %CONTINENTs
    %derive horizontal continent slices
    depthLevelSlicesCONT 	= DomainDepthMax - depthLevelCONTslices.*DomainDepth;  %plotting values
    
    depthLevelSlicesIdx = zeros(size(depthLevelSlicesCONT)); %array indices
    for id=1:length(depthLevelSlicesCONT)
        [~,depthLevelSlicesIdx(1,id)]	= min(abs(GRID.Z_3Dp(1,1,:)-depthLevelSlicesCONT(id)));  %index of closest value to mean LAB depth further up in the plate
    end
    
    contSlice             = CONT3D(:,:,depthLevelSlicesIdx);
    
    %define point sources for continents
    numCONTSlice = zeros(1,length(depthLevelSlicesCONT));
    for id=1:length(depthLevelSlicesCONT)
        %     CC11        = bwconncomp(plumesHotSlice(:,:,id))
        contSlice(:,:,id)	= bwmorph(contSlice(:,:,id),'shrink',Inf);
        try
            numCONTSlice(id)  	= sum(contSlice(:,:,id),'all');
        catch
            numCONTSlice(id)  	= sum(sum(contSlice(:,:,id)));
        end
    end
    %create continent data arrays
    maxNumCONTsPerSlice= max(numCONTSlice(:));
    dummy           	= zeros(size(contSlice,1),size(contSlice,2))*NaN;
    vhFullSlice         = dummy;
    vhCONTsSlices       = zeros(length(depthLevelSlicesCONT),maxNumCONTsPerSlice)*NaN; %[slices,continents]
    xCONTsSlices        = vhCONTsSlices;
    yCONTsSlices        = vhCONTsSlices;
    zCONTsSlices        = vhCONTsSlices;
    for id=1:length(depthLevelSlicesCONT)
        %extract horizontal position
        if max(contSlice(:,:,id))>0
            dummyNanoms = size(dummy(contSlice(:,:,id)==1),1);
            dummy(:,:)                          = GRID.X_3Dp(:,:,depthLevelSlicesIdx(id));
            xCONTsSlices(id,1:dummyNanoms)  	= dummy(contSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Y_3Dp(:,:,depthLevelSlicesIdx(id));
            yCONTsSlices(id,1:dummyNanoms)      = dummy(contSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Z_3Dp(:,:,depthLevelSlicesIdx(id));
            zCONTsSlices(id,1:dummyNanoms)      = dummy(contSlice(:,:,id)==1);
            %extract horizontal velocity (i.e., migration magnitude) of these points in the different slices
            vhFullSlice(:,:)                    = VH_3D(:,:,depthLevelSlicesIdx(id));
            vhCONTsSlices(id,1:dummyNanoms)     = vhFullSlice(contSlice(:,:,id)==1);
        end
    end
    clearvars dummy dummyNanoms

    %LLSVPs
    %derive horizontal llsvp slices
    depthLevelSlicesLLSVP 	= DomainDepthMax - depthLevelLLSVPslices.*DomainDepth;  %plotting values
    
    depthLevelSlicesIdx = zeros(size(depthLevelSlicesLLSVP)); %array indices
    for id=1:length(depthLevelSlicesLLSVP)
        [~,depthLevelSlicesIdx(1,id)]	= min(abs(GRID.Z_3Dp(1,1,:)-depthLevelSlicesLLSVP(id)));  %index of closest value to mean LAB depth further up in the plate
    end
    
    llsvpsSlice             = LLSVP3D(:,:,depthLevelSlicesIdx);
    
    %define point sources for llsvps
    numLLSVPSlice = zeros(1,length(depthLevelSlicesLLSVP));
    for id=1:length(depthLevelSlicesLLSVP)
        %     CC11        = bwconncomp(plumesHotSlice(:,:,id))
        llsvpsSlice(:,:,id) 	= bwmorph(llsvpsSlice(:,:,id),'shrink',Inf);
        try
            numLLSVPSlice(id)    	= sum(llsvpsSlice(:,:,id),'all');
        catch
            numLLSVPSlice(id)    	= sum(sum(llsvpsSlice(:,:,id)));
        end
    end
    
    %create llsvp data arrays
    maxNumLLSVPsPerSlice= max(numLLSVPSlice(:));
    dummy           	= zeros(size(llsvpsSlice,1),size(llsvpsSlice,2))*NaN;
    vhFullSlice         = dummy;
    vhLLSVPsSlices      = zeros(length(depthLevelSlicesLLSVP),maxNumLLSVPsPerSlice)*NaN; %[slices,llsvps]
    xLLSVPsSlices       = vhLLSVPsSlices;
    yLLSVPsSlices       = vhLLSVPsSlices;
    zLLSVPsSlices       = vhLLSVPsSlices;
    for id=1:length(depthLevelSlicesLLSVP)
        %extract horizontal position
        if max(llsvpsSlice(:,:,id))>0
            dummyNanoms = size(dummy(llsvpsSlice(:,:,id)==1),1);
            dummy(:,:)                          = GRID.X_3Dp(:,:,depthLevelSlicesIdx(id));
            xLLSVPsSlices(id,1:dummyNanoms)  	= dummy(llsvpsSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Y_3Dp(:,:,depthLevelSlicesIdx(id));
            yLLSVPsSlices(id,1:dummyNanoms) 	= dummy(llsvpsSlice(:,:,id)==1);
            dummy(:,:)                          = GRID.Z_3Dp(:,:,depthLevelSlicesIdx(id));
            zLLSVPsSlices(id,1:dummyNanoms) 	= dummy(llsvpsSlice(:,:,id)==1);
            %extract horizontal velocity (i.e., migration magnitude) of these points in the different slices
            vhFullSlice(:,:)                    = VH_3D(:,:,depthLevelSlicesIdx(id));
            vhLLSVPsSlices(id,1:dummyNanoms)	= vhFullSlice(llsvpsSlice(:,:,id)==1);
        end
    end
    clearvars dummy dummyNanoms
    
    
    %output diagnostics
    %Upper Mantle
    MANTLE.plumeHotNumberUM     = numPlumesHotSlice(1);
    MANTLE.plumeHotVHmaxUM      = max(vhPlumeHotSlices(1,:));
    MANTLE.plumeHotVHminUM    	= min(vhPlumeHotSlices(1,:));
    MANTLE.plumeHotVHmeanUM   	= nanmean(vhPlumeHotSlices(1,:));
    MANTLE.plumeColdNumberUM  	= numPlumesColdSlice(1);
    MANTLE.plumeColdVHmaxUM   	= max(vhPlumeColdSlices(1,:));
    MANTLE.plumeColdVHminUM    	= min(vhPlumeColdSlices(1,:));
    MANTLE.plumeColdVHmeanUM   	= nanmean(vhPlumeColdSlices(1,:));
    %Mid Mantle
    MANTLE.plumeHotNumberMM     = numPlumesHotSlice(2);
    MANTLE.plumeHotVHmaxMM      = max(vhPlumeHotSlices(2,:));
    MANTLE.plumeHotVHminMM    	= min(vhPlumeHotSlices(2,:));
    MANTLE.plumeHotVHmeanMM   	= nanmean(vhPlumeHotSlices(2,:));
    MANTLE.plumeColdNumberMM  	= numPlumesColdSlice(2);
    MANTLE.plumeColdVHmaxMM    	= max(vhPlumeColdSlices(2,:));
    MANTLE.plumeColdVHminMM    	= min(vhPlumeColdSlices(2,:));
    MANTLE.plumeColdVHmeanMM   	= nanmean(vhPlumeColdSlices(2,:));
    %Lower Mantle
    MANTLE.plumeHotNumberLM     = numPlumesHotSlice(3);
    MANTLE.plumeHotVHmaxLM      = max(vhPlumeHotSlices(3,:));
    MANTLE.plumeHotVHminLM    	= min(vhPlumeHotSlices(3,:));
    MANTLE.plumeHotVHmeanLM   	= nanmean(vhPlumeHotSlices(3,:));
    MANTLE.plumeColdNumberLM  	= numPlumesColdSlice(3);
    MANTLE.plumeColdVHmaxLM    	= max(vhPlumeColdSlices(3,:));
    MANTLE.plumeColdVHminLM    	= min(vhPlumeColdSlices(3,:));
    MANTLE.plumeColdVHmeanLM   	= nanmean(vhPlumeColdSlices(3,:));
    
    MANTLE.contVHmaxUM          = max(vhCONTsSlices(1,:));
    MANTLE.contVHminUM          = min(vhCONTsSlices(1,:));
    MANTLE.contVHmeanUM         = nanmean(vhCONTsSlices(1,:));
    
    MANTLE.llsvpVHmaxLM         = max(vhLLSVPsSlices(1,:));
    MANTLE.llsvpVHminLM         = min(vhLLSVPsSlices(1,:));
    MANTLE.llsvpVHmeanLM        = nanmean(vhLLSVPsSlices(1,:));
    
else
    warning('This has not been implemented yet for yinyang!')
end




%% DETAILS OF FOUND ANOMALIES II
% number of hot and cold active plumes
CChot2      	= bwconncomp(plumesHot);
CCcold2         = bwconncomp(plumesCold);
if strcmp(GRID.Type,'yinyang')
    CChot2_yang         	= bwconncomp(plumesHot_yang);
    CCcold2_yang         	= bwconncomp(plumesCold_yang);
    MANTLE.numHotPlumes     = CChot2.NumObjects+CChot2_yang.NumObjects;
    MANTLE.numColdPlumes    = CCcold2.NumObjects+CCcold2_yang.NumObjects;
else
    MANTLE.numHotPlumes     = CChot2.NumObjects;
    MANTLE.numColdPlumes    = CCcold2.NumObjects;
end

%check here also for wrap-around boundaries.....................


% total volume of up- and downwelling
if strcmp(GRID.Type,'yinyang')
    MANTLE.UpwellingVolume      = sum(sum(GRID.cellVolume(upWelling(:,:,2:end)==1)))+...
        sum(sum(GRID.cellVolume(upWelling_yang(:,:,2:end)==1))); %[nd] or [m^2]or[m^3]
    MANTLE.DownwellingVolume 	= sum(sum(GRID.cellVolume(downWelling(:,:,1:end-1)==1)))+...
        sum(sum(GRID.cellVolume(downWelling_yang(:,:,1:end-1)==1))); %[nd] or [m^2]or[m^3]
else
    MANTLE.UpwellingVolume      = sum(sum(GRID.cellVolume(upWelling(:,:,2:end)==1))); %[nd] or [m^2]or[m^3]
    MANTLE.DownwellingVolume 	= sum(sum(GRID.cellVolume(downWelling(:,:,1:end-1)==1))); %[nd] or [m^2]or[m^3]
end
% percentage of total volume of up- and downwelling versus total volume model domain
MANTLE.UpwellingVolPerc     = 100 *MANTLE.UpwellingVolume /(sum(GRID.cellVolume(:))*nb); %in [% of total volume]
MANTLE.DownwellingVolPerc  	= 100 *MANTLE.DownwellingVolume /(sum(GRID.cellVolume(:))*nb); %in [% of total volume]
% area of upwellings at certain depth..................

% active up- and downwelling
MANTLE.activeUpwelling                  = activeUpwelling;
MANTLE.activeDownwelling                = activeDownwelling;
MANTLE.activeUpwelling(~plumesHot)      = 0;
MANTLE.activeDownwelling(~plumesCold)   = 0;
CCactU                                  = bwconncomp(MANTLE.activeUpwelling);
CCactD                                  = bwconncomp(MANTLE.activeDownwelling);
MANTLE.numActUpwelling                  = CCactU.NumObjects;
MANTLE.numActDownwelling                = CCactD.NumObjects;
MANTLE.ActUpwellingVolume             	= sum(sum(GRID.cellVolume(MANTLE.activeUpwelling(:,:,2:end)==1)));      %[nd] or [m^2]or[m^3]
MANTLE.ActDownwellingVolume             = sum(sum(GRID.cellVolume(MANTLE.activeDownwelling(:,:,1:end-1)==1)));  %[nd] or [m^2]or[m^3]
MANTLE.ActUpwellingVolPerc              = 100 *MANTLE.ActUpwellingVolume /MANTLE.UpwellingVolume;               %in [% of total upwelling]
MANTLE.ActDownwellingVolPerc            = 100 *MANTLE.ActDownwellingVolume /MANTLE.DownwellingVolume;           %in [% of total downwelling]

% passive up- and downwelling
MANTLE.passiveUpwelling              	= passiveUpwelling;
MANTLE.passiveDownwelling             	= passiveDownwelling;
MANTLE.passiveUpwelling(plumesHot==1)   = 0;
MANTLE.passiveDownwelling(plumesCold==1)= 0;
CCpassU                                 = bwconncomp(MANTLE.passiveUpwelling);
CCpassD                                 = bwconncomp(MANTLE.passiveDownwelling);
MANTLE.numPassUpwelling               	= CCpassU.NumObjects;
MANTLE.numPassDownwelling               = CCpassD.NumObjects;
MANTLE.PassUpwellingVolume             	= sum(sum(GRID.cellVolume(MANTLE.passiveUpwelling(:,:,2:end)==1)));      %[nd] or [m^2]or[m^3]
MANTLE.PassDownwellingVolume            = sum(sum(GRID.cellVolume(MANTLE.passiveDownwelling(:,:,1:end-1)==1)));  %[nd] or [m^2]or[m^3]
MANTLE.PassUpwellingVolPerc             = 100 *MANTLE.PassUpwellingVolume /MANTLE.UpwellingVolume;               %in [% of total upwelling]
MANTLE.PassDownwellingVolPerc           = 100 *MANTLE.PassDownwellingVolume /MANTLE.DownwellingVolume;           %in [% of total downwelling]

% continents
if exist('CentroidCONT','var')
    MANTLE.contNumberUM                    = size(CentroidCONT,1);
    MANTLE.contLocationX                   = CentroidCONT(:,1);
    MANTLE.contLocationZ                   = CentroidCONT(:,2);
end

% llsvps
if exist('CentroidLLSVP','var')
    MANTLE.llsvpNumberLM                    = size(CentroidLLSVP,1);
    MANTLE.llsvpLocationX                   = CentroidLLSVP(:,1);
    MANTLE.llsvpLocationZ                   = CentroidLLSVP(:,2);
end

%% DISPLAY INFORMATION
disp('   Horizontal Mantle Flow')
disp('     Velocity')
disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(MANTLE.VHmeanUM,3),' ',SETUP.vDim,' (mean) ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.VHmaxUM,3),' ',SETUP.vDim,' (max) ',...
    ' at ',num2str(depthLevelDiagnosticsP(1,1),4),' ',GRID.Zdim,' depth'])
disp(['     ',STYLE.SCHAR.smallBullet,' mid mantle            = ',num2str(MANTLE.VHmeanMM,3),' ',SETUP.vDim,' (mean) ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.VHmaxMM,3),' ',SETUP.vDim,' (max) ',...
    ' at ',num2str(depthLevelDiagnosticsP(1,2),4),' ',GRID.Zdim,' depth'])
disp(['     ',STYLE.SCHAR.smallBullet,' lower mantle          = ',num2str(MANTLE.VHmeanLM,3),' ',SETUP.vDim,' (mean) ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.VHmaxLM,3),' ',SETUP.vDim,' (max) ',...
    ' at ',num2str(depthLevelDiagnosticsP(1,3),4),' ',GRID.Zdim,' depth'])

disp('   Mantle Upwelling')
disp('     Volume')
disp(['     ',STYLE.SCHAR.smallBullet,' total                 = ',num2str(MANTLE.UpwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.UpwellingVolPerc,2),' %vol)'])
disp(['     ',STYLE.SCHAR.smallBullet,' active                = ',num2str(MANTLE.ActUpwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.ActUpwellingVolPerc,2),' %)'])
%     if MANTLE.numHotPlumes==0
%         disp('   No Active Hot Plume');
%     elseif MANTLE.numHotPlumes==1
%         disp('   1 Active Hot Plume:');
%     else
%         disp(['   ',num2str(MANTLE.numHotPlumes),' Active Hot Plumes:']);
%     end
disp(['     ',STYLE.SCHAR.smallBullet,' passive               = ',num2str(MANTLE.PassUpwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.PassUpwellingVolPerc,2),' %)'])

if MANTLE.plumeHotNumberLM>0
    disp('     Hot-plume mobility')
    if MANTLE.plumeHotNumberUM>0
        if MANTLE.plumeHotNumberUM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(MANTLE.plumeHotVHminUM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeHotVHmeanUM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeHotVHmaxUM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeHotNumberUM,1),' plume',pluralS,' at ',num2str(depthLevelSlicesPL(1),4),' ',GRID.Zdim,' depth)'])
    end
    if MANTLE.plumeHotNumberMM>0
        if MANTLE.plumeHotNumberMM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' mid mantle            = ',num2str(MANTLE.plumeHotVHminMM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeHotVHmeanMM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeHotVHmaxMM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeHotNumberMM,1),' plume',pluralS,' at ',num2str(depthLevelSlicesPL(2),4),' ',GRID.Zdim,' depth)'])
    end
    if MANTLE.plumeHotNumberLM>0
        if MANTLE.plumeHotNumberLM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' lower mantle          = ',num2str(MANTLE.plumeHotVHminLM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeHotVHmeanLM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeHotVHmaxLM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeHotNumberLM,1),' plume',pluralS,' at ',num2str(depthLevelSlicesPL(3),4),' ',GRID.Zdim,' depth)'])
    end
end

disp('   Mantle Downwelling')
disp('     Volume')
disp(['     ',STYLE.SCHAR.smallBullet,' total                 = ',num2str(MANTLE.DownwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.DownwellingVolPerc,2),' %vol)'])
disp(['     ',STYLE.SCHAR.smallBullet,' active                = ',num2str(MANTLE.ActDownwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.ActDownwellingVolPerc,2),' %)'])
%     if MANTLE.numColdPlumes==0
%         disp('   No Active Cold Plume');
%     elseif MANTLE.numColdPlumes==1
%         disp('   1 Active Cold Plume:');
%     else
%         disp(['   ',num2str(MANTLE.numColdPlumes),' Active Cold Plumes:']);
%     end
disp(['     ',STYLE.SCHAR.smallBullet,' passive               = ',num2str(MANTLE.PassDownwellingVolume,3),' ',GRID.VolDim,...
    ' (',num2str(MANTLE.PassDownwellingVolPerc,2),' %)'])

if MANTLE.plumeColdNumberLM>0
    disp('     Cold-plume mobility')
    if MANTLE.plumeColdNumberUM>0
        if MANTLE.plumeColdNumberUM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(MANTLE.plumeColdVHminUM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeColdVHmeanUM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeColdVHmaxUM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeColdNumberUM,1),' plume',pluralS,' at ',num2str(depthLevelSlicesPL(1),4),' ',GRID.Zdim,' depth)'])
    end
    if MANTLE.plumeColdNumberMM>0
        if MANTLE.plumeColdNumberMM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' mid mantle            = ',num2str(MANTLE.plumeColdVHminMM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeColdVHmeanMM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeColdVHmaxMM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeColdNumberMM,1),' plume',pluralS,' at ',num2str(depthLevelSlicesPL(2),4),' ',GRID.Zdim,' depth)'])
    end
    if MANTLE.plumeColdNumberLM>0
        if MANTLE.plumeColdNumberLM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' lower mantle          = ',num2str(MANTLE.plumeColdVHminLM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.plumeColdVHmeanLM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.plumeColdVHmaxLM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.plumeColdNumberLM,1),' plume',pluralS,' at ',num2str(depthLevelSlicesPL(3),4),' ',GRID.Zdim,' depth)'])
    end
end

if MANTLE.contNumberUM>0
    disp('   Continents')
    disp('     Continent mobility')
    if MANTLE.contNumberUM>0
        if MANTLE.contNumberUM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(MANTLE.contVHminUM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.contVHmeanUM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.contVHmaxUM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.contNumberUM,1),' Continent',pluralS,' at ',num2str(depthLevelSlicesCONT(1),4),' ',GRID.Zdim,' depth)'])
    end
end

if MANTLE.llsvpNumberLM>0
    disp('   LLSVPs')
    disp('     LLSVP mobility')
    if MANTLE.llsvpNumberLM>0
        if MANTLE.llsvpNumberLM>1; pluralS = 's'; else; pluralS = ''; end
        disp(['     ',STYLE.SCHAR.smallBullet,' lower mantle          = ',num2str(MANTLE.llsvpVHminLM,2),' ',STYLE.SCHAR.leftArrow,' ',num2str(MANTLE.llsvpVHmeanLM,2),...
            ' ',STYLE.SCHAR.rightArrow,' ',num2str(MANTLE.llsvpVHmaxLM,2),' ',SETUP.vDim,...
            ' (',num2str(MANTLE.llsvpNumberLM,1),' LLSVP',pluralS,' at ',num2str(depthLevelSlicesLLSVP(1),4),' ',GRID.Zdim,' depth)'])
    end
end



%% PLOTTING PLUMES SEPARATELY
if MD.plotPLUMES && strcmp(GRID.Dim,'3-D')
    figure(22),clf
    
    x = 1:size(Plumes3D,1);
    y = 1:size(Plumes3D,2);
    z = 1:size(Plumes3D,3);
    [x3d,y3d,z3d] = meshgrid(y,x,z);
    
    phot = patch(isosurface(x3d,y3d,z3d,plumesHot,+0.95)); %HOT
    isonormals(x3d,y3d,z3d,plumesHot,phot)
    set(phot,'FaceColor','red','EdgeColor','none');
    view(3);
    camlight
    lighting phong  %gouraud
    
    hold on
    pcold = patch(isosurface(x3d,y3d,z3d,plumesCold,+0.95)); %COLD
    isonormals(x3d,y3d,z3d,plumesCold,pcold)
    set(pcold,'FaceColor','blue','EdgeColor','none');
    view(3);
    axis([0 max(max(max(x3d))) 0 max(max(max(y3d))) 0 max(max(max(z3d)))])
    camlight
    lighting phong  %gouraud
    
    title('after applying criterion')
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal
    box on
    grid on
    
    figure(23),clf
    isosurface(Plumes3D,+0.95,'noshare') %HOT
    hold on
    isosurface(Plumes3D,-0.95,'noshare') %COLD
    
    title('original')
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal
    axis([0 max(max(max(x3d))) 0 max(max(max(y3d))) 0 max(max(max(z3d)))])
    box on
    grid on
    
    figure(1)  %go back to figure(1)
end


%% FUNCTION OUTPUT
MANTLE.upWelling            = upWelling;    %1 for upwelling, 0 else  (depending on threshold)
MANTLE.downWelling          = downWelling;  %1 for downwelling, 0 else  (depending on threshold)
MANTLE.upWellingAbsolute    = upWellingAbsolute;    %1 for upwelling, 0 else
MANTLE.downWellingAbsolute 	= downWellingAbsolute;  %1 for downwelling, 0 else
MANTLE.plumesHot            = plumesHot;    %1 for hot plumes, 0 else
MANTLE.plumesCold           = plumesCold;   %1 for cold plumes, 0 else
% MANTLE.numHotPlumes               %num
% MANTLE.numColdPlumes              %num
% MANTLE.UpwellingVolume            %[nd] or [m^2]or[m^3]
% MANTLE.DownwellingVolume          %[nd] or [m^2]or[m^3]
% MANTLE.UpwellingVolPerc           %in [% of model domain]
% MANTLE.DownwellingVolPerc         %in [% of model domain]

% MANTLE.activeUpwelling            %1 for active upwelling
% MANTLE.activeDownwelling        	%1 for active downwelling
% MANTLE.numActUpwelling           	%num
% MANTLE.numActDownwelling          %num
% MANTLE.ActUpwellingVolume      	%[nd] or [m^2]or[m^3]
% MANTLE.ActDownwellingVolume      	%[nd] or [m^2]or[m^3]
% MANTLE.ActUpwellingVolPerc    	%in [% of total upwelling]
% MANTLE.ActDownwellingVolPerc     	%in [% of total downwelling]

% MANTLE.passiveUpwelling           %1 for passive upwelling
% MANTLE.passiveDownwelling        	%1 for passive downwelling
% MANTLE.numPassUpwelling           %num
% MANTLE.numPassDownwelling         %num
% MANTLE.PassUpwellingVolume      	%[nd] or [m^2]or[m^3]
% MANTLE.PassDownwellingVolume      %[nd] or [m^2]or[m^3]
% MANTLE.PassUpwellingVolPerc    	%in [% of total upwelling]
% MANTLE.PassDownwellingVolPerc     %in [% of total downwelling]

if strcmp(GRID.Type,'yinyang')
    %add here......
end






















