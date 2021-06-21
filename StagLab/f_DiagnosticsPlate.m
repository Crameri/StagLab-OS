
%%                                                  PLATE DIAGNOSTICS 7.310
%
%   . calls f_Connectivity
%   . includes f_SlabTipDiagnostics
%   . includes f_BendingRadiusB
%   . includes f_slabDiagnostics
%   . includes f_SubTopoCharacteristics
%   . includes f_stagnantLidDiagnostics
%   . includes rms
%
% returns PLATE.Subduction=='noTracking' if indication not suitable/possible
% returns PLATE.StagnantLid==true if stagnant lid detected
%
%                                                Fabio Crameri, 30.10.2019

%% NOTES
% be careful with diagnostics for spherical2D geometry (i.e., use GRID.X_3Dsp for actual horizontal distance measure!)

function [PLATE,PLOT] = f_DiagnosticsPlate(FILE,GRID,SETUP,SWITCH,PLOT,STYLE)
%% SWITCHES
SlabEtaToTake               = 'min';      %'mean' or 'min' for slab-mantle viscosity difference calculation

%% DEFAULTS
%mobile lid variables
%if strcmp(GRID.Dim,'2-D')
    PLATE.SubPolarity               = 'na';
    PLATE.SubPolarityNum            = 0;
    PLATE.idxTrench                 = NaN;
    PLATE.idxUP                     = NaN;
    PLATE.idxLP                     = NaN;
    PLATE.subPlateIndices           = NaN;
    PLATE.subPlateFullFit       	= NaN;
    PLATE.subPlateFit               = NaN;
    PLATE.TrenchVelocity           	= NaN;
    PLATE.theoreticTrenchVelocity   = NaN;
    PLATE.UPvelocity                = NaN;
    PLATE.LPvelocity                = NaN;
    PLATE.PlateConvergence          = NaN;
    PLATE.PlateDivergence           = NaN;
    PLATE.leftPvelocity             = NaN;
    PLATE.rightPvelocity            = NaN;
    PLATE.PlateThicknessL           = NaN;
    PLATE.PlateThicknessR           = NaN;
    PLATE.PlateThicknessWideL       = NaN;
    PLATE.PlateThicknessWideR       = NaN;
    PLATE.PlateThicknessUP          = NaN;
    PLATE.PlateThicknessLP          = NaN;
    %topography
    PLATE.TopoOverridingPlate     	= NaN;
    PLATE.TopoSubductingPlate     	= NaN;
    PLATE.TrenchDepth               = NaN;
    PLATE.TrenchDepthX             	= NaN;
    PLATE.BackArcBasinX          	= NaN;
    PLATE.BackArcBasinZ           	= NaN;
    PLATE.BackArcBasinExtentX       = NaN;
    PLATE.BackArcBasinExtentZ       = NaN;
    PLATE.InundationX               = NaN;
    PLATE.InundationZ               = NaN;
    PLATE.InundationDistance        = NaN;
    PLATE.IslandArcX                = NaN;
    PLATE.IslandArcZ                = NaN;
    PLATE.ForeBulgeX                = NaN;
    PLATE.ForeBulgeZ                = NaN;
    PLATE.BackArcBasinArea       	= NaN;
    PLATE.UPtiltAngle             	= NaN;
    PLATE.UPtiltAngleXnear        	= NaN;
    PLATE.UPtiltAngleXfar         	= NaN;
    PLATE.UPtiltAngleZnear        	= NaN;
    PLATE.UPtiltAngleZfar         	= NaN;
    %slab
    PLATE.ShallowSlabAngle        	= NaN;
    PLATE.ShallowSlabAngleDepth    	= NaN;
    PLATE.SlabTipPosition         	= NaN;
    PLATE.SlabTipAngle              = NaN;
    PLATE.SlabTipAnglePoint1        = [NaN,NaN];
    PLATE.SlabTipAnglePoint2        = [NaN,NaN];
    PLATE.SlabTipVX                 = NaN;
    PLATE.SlabTipVZ                 = NaN;
    PLATE.SlabViscosity             = NaN;
    PLATE.SlabViscosityMin          = NaN;
    PLATE.SlabSinkingVelocity      	= NaN;
    PLATE.SlabMantleViscDiff        = NaN;
    PLATE.BendingRadius             = NaN;
    PLATE.bendingCircleCenter       = NaN;
    PLATE.ViscDissBending           = NaN;
    PLATE.ViscDissBendingRel      	= NaN;
    PLATE.RidgeVelocity             = NaN;
    PLATE.ridgeLPvelocity           = NaN;
    PLATE.ridgeRPvelocity           = NaN;
    PLATE.shallowSlabAge            = NaN;
    PLATE.waterRemainedInSlab       = NaN;
    PLATE.SlabTemperatureMin        = NaN;
    PLATE.SlabMantleTempDiff        = NaN;
    PLATE.SlabDensityMax            = NaN;
    PLATE.SlabStressMax             = NaN;
    
%elseif strcmp(GRID.Dim,'3-D')
    PLATE.SubLength                 = NaN;
    PLATE.SprLength                 = NaN;
    PLATE.SubLengthTotal         	= NaN;
    PLATE.SprLengthTotal          	= NaN;
    PLATE.SubCenterP                = NaN; 
    PLATE.SprCenterP                = NaN;  
    PLATE.UpperPlateP               = NaN;     
    PLATE.LowerPlateP               = NaN;  
    PLATE.SubNormalP                = NaN;    
    PLATE.SubStrikeP                = NaN; 
    PLATE.SprNormalP                = NaN;   
    PLATE.SprStrikeP                = NaN;
%end
%plate diagnostics
PLATE.BaseTemperature               = NaN;
PLATE.SlabContourTemperature        = NaN;
PLATE.maxStress                     = NaN;
PLATE.surfaceLevelIdx             	= NaN;
PLATE.surfXVelocity                 = NaN;
PLATE.surfYVelocity                 = NaN;
if strcmp(GRID.Type,'yinyang')
    PLATE.surfXVelocity_yang     	= NaN;
    PLATE.surfYVelocity_yang     	= NaN;
end
PLATE.coreLevelIdx                  = NaN;
PLATE.coreStress                    = NaN; %currently not saved
PLATE.coreStrainrate                = NaN; %currently not saved
PLATE.coreViscosity                 = NaN; %currently not saved
PLATE.coreStressMin                 = NaN;
PLATE.coreStressMax                 = NaN;
PLATE.coreStrainrateMin             = NaN;
PLATE.coreStrainrateMax             = NaN;
PLATE.coreViscosityMin              = NaN;
PLATE.coreViscosityMax              = NaN;
PLATE.PlateVelocity                 = NaN;
PLATE.coreXVelocity                 = NaN;
PLATE.coreYVelocity                 = NaN;
if strcmp(GRID.Type,'yinyang')
    PLATE.PlateVelocity_yang        = NaN;
    PLATE.coreXVelocity_yang     	= NaN;
    PLATE.coreYVelocity_yang     	= NaN;
end
PLATE.subductionFlowRate          	= NaN;
PLATE.ViscDissipationPlate         	= NaN;
PLATE.ViscDissipationPlateNoCrust   = NaN;
PLATE.PlateVelocityMax              = NaN;
PLATE.PlateVelocityMin              = NaN;
PLATE.PlateVelocityDiff             = NaN;
PLATE.PlateVelocityRMS              = NaN;
PLATE.PlateMobility                 = NaN;     
PLATE.PlateDrivety                  = NaN;
%mantle diagnostics
PLATE.UMantleVelocityMean           = NaN;
PLATE.UMantleVelocityMax            = NaN;
PLATE.UMantleHVelocityMean          = NaN;
PLATE.UMantleHVelocityMax          	= NaN;
PLATE.UMantleRVelocityMax           = NaN;
PLATE.UMantleViscosity              = NaN;
PLATE.UMantleTemperature            = NaN;
PLATE.UMantleDensity                = NaN;
%stagnant lid variables
PLATE.StagnantLid                   = false;
PLATE.LABdepth                      = NaN;
PLATE.maxYieldDepth                 = NaN;
PLATE.maxYieldDepthFraction         = NaN;
%plate surface variables
PLATE.SurfaceAgeMean                = NaN;
PLATE.SurfaceAgeMedian              = NaN;
PLATE.SurfaceAgeMeanCont         	= NaN;
PLATE.SurfaceAgeMedianCont        	= NaN;
PLATE.SurfaceAgeMeanOcean        	= NaN;
PLATE.SurfaceAgeMedianOcean      	= NaN;
PLATE.CrustalThicknessMean        	= NaN;
PLATE.CrustalThicknessMedian    	= NaN;
PLATE.CrustalThicknessMeanCont   	= NaN;
PLATE.CrustalThicknessMedianCont   	= NaN;
PLATE.CrustalThicknessMeanOcean   	= NaN;
PLATE.CrustalThicknessMedianOcean 	= NaN;

%% READ INPUT DATA
%temperature
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Temperature';
DATA.FieldAbbreviation      = 'T';
DATA.StopExecutionIfNotFound = true;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
interiorIndexLow    = floor( size(GRID.Z_3D,3)*1/10 );  %bottom of the mantle
interiorIndexHigh   = ceil( size(GRID.Z_3D,3)*8/10 );   %top of the mantle
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

%velocity
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Velocity';
DATA.FieldAbbreviation      = 'V';
DATA.StopExecutionIfNotFound = true;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    VX_3D = PLOT.VX_3Dyin; VX_3Dyang = PLOT.VX_3Dyang;
    VY_3D = PLOT.VY_3Dyin; VY_3Dyang = PLOT.VY_3Dyang;
    VZ_3D = PLOT.VZ_3Dyin; VZ_3Dyang = PLOT.VZ_3Dyang;
    V_3D = PLOT.V_3Dyin; V_3Dyang = PLOT.V_3Dyang;
    v_rmsApproximate    = mean2([V_3D;V_3Dyang]);   %doesn't take variable gridspacing nor yinyang obsolete grid corners into account....
else %all other grid types
    VX_3D = PLOT.VX_3D; VY_3D = PLOT.VY_3D; VZ_3D = PLOT.VZ_3D; V_3D = PLOT.V_3D;
    v_rmsApproximate    = mean2(V_3D);   %doesn't take variable gridspacing into account....
end

%viscosity
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Viscosity';
DATA.FieldAbbreviation      = 'ETA';
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    ETA_3D = PLOT.ETA_3Dyin; ETA_3Dyang = PLOT.ETA_3Dyang;
else %all other grid types
    ETA_3D = PLOT.ETA_3D;
end

%density
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Density';
DATA.FieldAbbreviation      = 'RHO';
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    RHO_3D = PLOT.RHO_3Dyin; RHO_3Dyang = PLOT.RHO_3Dyang;
else %all other grid types
    RHO_3D = PLOT.RHO_3D;
end

%stress
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Stress';
DATA.FieldAbbreviation      = 'STR';
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    STR_3D = PLOT.STR_3Dyin; STR_3Dyang = PLOT.STR_3Dyang;
else %all other grid types
    STR_3D = PLOT.STR_3D;
end

%strainrate
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Strain rate';
DATA.FieldAbbreviation      = 'EDOT';
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    EDOT_3D = PLOT.EDOT_3Dyin; EDOT_3Dyang = PLOT.EDOT_3Dyang;
else %all other grid types
    EDOT_3D = PLOT.EDOT_3D;
end

%basalt
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Basalt';
DATA.FieldAbbreviation      = 'BASALT';
[~,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    BASALT_3D = PLOT.BASALT_3Dyin; BASALT_3Dyang = PLOT.BASALT_3Dyang;
else %all other grid types
    BASALT_3D = PLOT.BASALT_3D;
end
if isnan(BASALT_3D(1))
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Composition';
    DATA.FieldAbbreviation      = 'C';
    [~,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    if strcmp(GRID.Type,'yinyang')
        C_3D = PLOT.C_3Dyin; C_3Dyang = PLOT.C_3Dyang;
    else %all other grid types
        C_3D = PLOT.C_3D;
    end
    if isnan(C_3D(1))
        %nothing to do
    else
        if max(C_3D(:))<=1
            BASALT_3D = C_3D;  %basalt should be from 0 to 1
        elseif min(C_3D(:))>=1
            BASALT_3D = min(max(0,-C_3D+2),1);  %basalt should be from 0 to 1
        else
            warning('Check here to adjust current Composition value handling!')
            BASALT_3D = min(max(0,-C_3D+2),1);  %basalt should be from 0 to 1  This might therefore need adjustment!
            BASALT_3D = BASALT_3D*NaN;
        end
    end
end

%continental crust
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Cont. crust';
DATA.FieldAbbreviation      = 'CC';
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    CC_3D = PLOT.CC_3Dyin; CC_3Dyang = PLOT.CC_3Dyang;
else %all other grid types
    CC_3D = PLOT.CC_3D;
end
if ~isnan(CC_3D(1))
    CCvalue         = CC_3D(:,:,end);       %get horizontal profile of continental crust
    CClocation      = zeros(size(CCvalue)); %apply threshold
    CClocation(CCvalue>0.2)     = 1;
    if strcmp(GRID.Type,'yinyang')
        CCvalue         = CC_3Dyang(:,:,end);       %get horizontal profile of continental crust
        CClocation_yang	= zeros(size(CCvalue)); %apply threshold
        CClocation_yang(CCvalue>0.2)	= 1;
    end
    clearvars CCvalue
end


% %air
% if strcmp(SETUP.topBC,'sticky-air')
%     DATA.Task                   = 'ImportFieldData';
%     DATA.Field2Import           = 'Air';
%     DATA.FieldAbbreviation      = 'AIR';
%     [~,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
%     if strcmp(GRID.Type,'yinyang')
%         AIR_3D = PLOT.AIR_3Dyin; AIR_3Dyang = PLOT.AIR_3Dyang;
%     else %all other grid types
%         AIR_3D = PLOT.AIR_3D;
%     end
% end

% GRAPH DATA:
%crustal thickness
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Crustal thickness';
DATA.FieldAbbreviation      = 'CR';
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    CR_3D = PLOT.CR_3Dyin; CR_3Dyang = PLOT.CR_3Dyang;
    dummy       = [CR_3D;CR_3Dyang];
    CRmedian 	= median(dummy(:)); %handle with care: median does not consider effects of an air layer or refined/spherical grid...
    CRmean      = mean2(dummy); %handle with care: median does not consider effects of an air layer or refined/spherical grid...
    if exist('CClocation','var')
        clearvars dummy
        dummy           = [CR_3D(CClocation==1);CR_3Dyang(CClocation_yang==1)];
        CRmedianCont 	= median(dummy(:));
        CRmeanCont      = mean2(dummy);
        clearvars dummy
        dummy           = [CR_3D(CClocation==0);CR_3Dyang(CClocation_yang==0)];
        CRmedianOcean 	= median(dummy(:));
        CRmeanOcean   	= mean2(dummy);
    end
    
else %all other grid types
    CR_3D = PLOT.CR_3D;
    CRmedian 	= median(CR_3D(:));
    CRmean      = mean2(CR_3D);
    if exist('CClocation','var')
        clearvars dummy
        dummy           = CR_3D(CClocation==1);
        CRmedianCont 	= median(dummy(:));
        CRmeanCont      = mean2(dummy);
        clearvars dummy
        dummy           = CR_3D(CClocation==0);
        CRmedianOcean 	= median(dummy(:));
        CRmeanOcean   	= mean2(dummy);
    end
end
if ~exist('CClocation','var')
    CRmedianCont 	= NaN;
    CRmeanCont      = NaN;
    CRmedianOcean 	= NaN;
    CRmeanOcean   	= NaN;
end

%surface age
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Surface age';
DATA.FieldAbbreviation      = 'SAGE';
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    SAGE_3D     = PLOT.SAGE_3Dyin; SAGE_3Dyang = PLOT.SAGE_3Dyang;
    dummy       = [SAGE_3D;SAGE_3Dyang];
    SAGEmedian 	= median(dummy(:)); %handle with care: median does not consider effects of an air layer or refined/spherical grid...
    SAGEmean    = mean2(dummy); %handle with care: median does not consider effects of an air layer or refined/spherical grid...
    if exist('CClocation','var')
        clearvars dummy
        dummy           = [SAGE_3D(CClocation==1);SAGE_3Dyang(CClocation_yang==1)];
        SAGEmedianCont 	= median(dummy(:));
        SAGEmeanCont  	= mean2(dummy);
        clearvars dummy
        dummy           = [SAGE_3D(CClocation==0);SAGE_3Dyang(CClocation_yang==0)];
        SAGEmedianOcean	= median(dummy(:));
        SAGEmeanOcean 	= mean2(dummy);
    end
    
else %all other grid types
    SAGE_3D     = PLOT.SAGE_3D;
    SAGEmedian 	= median(SAGE_3D(:)); 
    SAGEmean    = mean2(SAGE_3D);
    if exist('CClocation','var')
        clearvars dummy
        dummy           = SAGE_3D(CClocation==1);
        SAGEmedianCont 	= median(dummy(:));
        SAGEmeanCont    = mean2(dummy);
        clearvars dummy
        dummy           = SAGE_3D(CClocation==0);
        SAGEmedianOcean	= median(dummy(:));
        SAGEmeanOcean 	= mean2(dummy);
    end
end
if ~exist('CClocation','var')
    SAGEmedianCont 	= NaN;
    SAGEmeanCont    = NaN;
    SAGEmedianOcean	= NaN;
    SAGEmeanOcean 	= NaN;
end

%% MAKE ADDITIONAL FIELDS
%viscous dissipation (phi=sigma*edot)
if strcmp(GRID.Type,'yinyang')
%     if isfield(PLOT,'VD_3Dyin') && isfield(PLOT,'VD_3Dyang')
%         VD_3D = PLOT.VD_3Dyin; VD_3Dyang = PLOT.VD_3Dyang;
%     else
%         VD_3D       = STR_3D.*EDOT_3D;
%         VD_3Dyang  = STR_3Dyang.*EDOT_3Dyang;
%         if SWITCH.Verbose; warning('viscous dissipation for YY is not yet multiplied by grid cell volume! Adjust here.'); end
%     end
    VD_3D               = NaN; %not implemented for YY yet; see below.
    VD_3Dyang          = NaN;
    VDnoBasalt_3D       = NaN;
    VDnoBasalt_3Dyang 	= NaN;
else
    if isfield(PLOT,'VD_3D')
        VD_3D  	= PLOT.VD_3D;
    else
        if ~isnan(STR_3D(1)) && size(STR_3D,1)>1 && size(STR_3D,1)>1 && ...
                ~isnan(EDOT_3D(1)) && size(EDOT_3D,1)>1 && size(EDOT_3D,1)>1  %~NaN
            %VD_3D       = ETA_3D.*EDOT_3D.^2;   %different expression for viscous dissipation rate per unit volume [W/m^3]
            VD_3D       = STR_3D.*EDOT_3D;      %viscous dissipation rate per unit volume [W/m^3]
            %multiply by grid cell area
            VD_3D       = VD_3D.*GRID.cellVolume;  %viscous dissipation rate [W/m]
        else
            VD_3D       = NaN;
        end
        
        if size(VD_3D,1)>1 && size(VD_3D,1)>1 && size(BASALT_3D,1)>1 && size(BASALT_3D,1)>1
            VDnoBasalt_3D   = VD_3D;
            VDnoBasalt_3D(BASALT_3D>0.7) = 0;
            if max(BASALT_3D(:)>1); warning('BASALT_3D has values >1: check here!'); end
        else
            VDnoBasalt_3D   = NaN;
        end
    end
end

%% INITIAL CHECKS
%check for isothermal model
if min(T_3D(:))+0.01*min(T_3D(:))>=median(T_3D(:))
    if SWITCH.Verbose; warning('isothermal model detected: no plate diagnostics possible.'); end
    PLATE.Subduction = 'noTracking'; %flag for no indication
    return %exit function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end
%check for isoviscous model
if logical(0)
if median(ETA_3D(:))+0.01*median(ETA_3D(:))>=max(ETA_3D(:))
    if SWITCH.Verbose; warning('isoviscous model detected: no plate diagnostics possible.'); end
    PLATE.Subduction = 'noTracking'; %flag for no indication
    return %exit function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end
end

%% GRID ADJUSTMENTS
Z_3D            = GRID.Z_3Dp; %[nd],[km],[m], or [cm]
if strcmp(GRID.Type,'spherical2D')
    % RE-FLIP DEPTH VECTOR AND ACCOUNT FOR STICKY-AIR LAYER
    Z_3D        = Z_3D +SETUP.d_air*GRID.dimFactor;  %set z=0 to real surface
    Z_3D        = 1*GRID.dimFactor -Z_3D; % reverse z-axis  (1.)
    %Radius of surface for conversion to radians
    Rsurface 	= GRID.rcmb+max(GRID.Z_3D(1,1,:)); %radius [nd] or [m]
end

%% FIND APPROXIMATE MEAN PLATE THICKNESS (less accurate when strong plate-thickness variation by e.g. deep subduction)
ThorizontalMean = zeros(size(GRID.Z_3D,3),1);
dTdz = ThorizontalMean;
for iz=1:size(GRID.Z_3D,3) %scan horizontal mean T from bottom to top
    if strcmp(GRID.Type,'yinyang')
        dummy             	= [PLOT.T_3Dyin(:,:,iz);PLOT.T_3Dyang(:,:,iz)]; %not perfect as yy-patch-edges are not taken care of!!!!!!!
    else
        dummy             	= PLOT.T_3D(:,:,iz);
    end
    ThorizontalMean(iz,1)   = median(dummy(:));
    clearvars dummy
    if iz>1
        dTdz(iz,1)          = (ThorizontalMean(iz-1,1)-ThorizontalMean(iz,1))/GRID.dzp(1,1,iz);
    end
end
LABidx = NaN; LABdepth = NaN;
for iz=length(dTdz):-1:2 %scan from top to bottom
    if iz<size(dTdz,1) && GRID.Z_3Dp(1,1,iz)>1/200*max(GRID.Z_3Dp(:)) && ...
        dTdz(iz,1)<0.1*max(dTdz(end-round(length(dTdz)*2/3):end)) %no significant change of temperature with depth (neglecting lower-most part of the domain)
        LABidx              = iz;
        LABdepth            = GRID.Z_3Dp(1,1,iz);
        break
    end
end
PLATE.LABdepth            	= LABdepth;

if logical(0) %test plot
    z(:) = GRID.Z_3Dp(1,1,:);
    plot(dTdz,z)
    hold on
    plot(dTdz(LABidx,1),GRID.Z_3Dp(1,1,LABidx),'or')
    axis ij
end

%% CHECK IF PLATE EXISTS I
if isnan(LABdepth)
    warning('LABdepth could not be found!');
    PLATE.Subduction = 'noTracking'; %flag for no indication
    
    return %exit function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end


%% CRITICAL VALUES FOR PLATE THICKNESS
LithoBaseTemp1                  = 1600;                                         %[K]; for lithosphere definition - 1600
coreFractionLithoDepth          = 2/3;                                          %Fraction of plate-core depth relative to the plate base - 2/3
surfaceFractionLithoDepth       = 0.05;                                         %Fraction of plate-surface depth relative to the plate base - 0.05

%% CRITICAL VALUES FOR STAGNANT-LID CHECK
if ~SWITCH.DimensionalMode && ~SWITCH.DimensionalInput %adjustment for non-dimensional mode
    LithoBaseTemp1              = LithoBaseTemp1/SETUP.deltaT;
    LITHO.zmax                	= LABdepth;                                     %[nd]; for lithosphere definition ~130
    LITHO.zmin                	= min(GRID.Z_3D(1,1,GRID.Z_3D(1,1,:)>0));       %[nd]; for lithosphere definition ~5
else %dimensional mode
    LITHO.zmax               	= LABdepth /GRID.m2p/1e3;                       %[km]; for lithosphere definition ~130
    LITHO.zmin                 	= min(GRID.Z_3D(1,1,GRID.Z_3D(1,1,:)>0)) /1e3;  %[km]; for lithosphere definition ~5
end
% diffVal                         = 100;                                        %[km]; critical minimum plate thickness variation

% PLATE.BaseTemperature           = 1300;                                       %[K]; for lithosphere definition - 1300
% PLATE.BaseTemperature           = Tmin+ (Tmean-Tmin)*0.75;                    %mean is not correct for variable grid cell sizes & sticky-air!
PLATE.BaseTemperature           = Tmin+ (Tmedian-Tmin)*0.71;                 	%[K]; for lithosphere definition ~1300
PLATE.SlabContourTemperature    = Tmin+ (Tmedian-Tmin)*0.84;                    %[K]; for slab definition
% PLATE.BaseTemperature           = Tmin+ (Tmax-Tmin)*0.65;                 	%[K]; for lithosphere definition ~1300
% PLATE.SlabContourTemperature    = Tmin+ (Tmax-Tmin)*0.8;                      %[K]; for slab definition

if LithoBaseTemp1>max(T_3D(:))*SETUP.Tscale
    LithoBaseTemp1              = max(T_3D(:))*SETUP.Tscale-100;              	%make sure it is not exceeding maximum temperature of the model
    if SWITCH.Verbose; warning(['Lithosphere base temperature reduced to ',num2str(LithoBaseTemp1),' ',SETUP.TDim,'.']); end
end

%% SETUP
if strcmp(GRID.Type,'yinyang')
    nb = 2; %account for two patches
else
    nb = 1;
end
%variable conversions
CONAREA.xPeriodic = strcmp(GRID.BCxx,'permeable');
CONAREA.yPeriodic = strcmp(GRID.BCyy,'permeable');
if strcmp(GRID.Type,'yinyang')
    CONAREA.xPeriodic = false;
    CONAREA.yPeriodic = false;
end

%% DIMENSIONALISATION
PLATE.BaseTemperature   = PLATE.BaseTemperature *SETUP.Tscale;
PLATE.SlabContourTemperature = PLATE.SlabContourTemperature *SETUP.Tscale;
T_3D                    = T_3D .*SETUP.Tscale;                	%[K] ...or whatever
VX_3D                   = VX_3D .*SETUP.Vscale;               	%[cm/a]
VY_3D                   = VY_3D .*SETUP.Vscale;              	%[cm/a]
VZ_3D                   = VZ_3D .*SETUP.Vscale;              	%[cm/a]
V_3D                    = V_3D .*SETUP.Vscale;               	%[cm/a]
v_rmsApproximate        = v_rmsApproximate .*SETUP.Vscale;    	%[cm/a]
ETA_3D                  = ETA_3D .*SETUP.etascale;            	%[Pas]
RHO_3D                  = RHO_3D .*SETUP.rhoscale;           	%[kg/m^3]
STR_3D                  = STR_3D .*SETUP.stressscale;        	%[MPa]
EDOT_3D                 = EDOT_3D .*SETUP.edotscale;         	%[1/s]
VD_3D                   = VD_3D .*SETUP.dissscale;          	%[W/m]
VDnoBasalt_3D         	= VDnoBasalt_3D .*SETUP.dissscale;    	%[W/m]
CR_3D                   = CR_3D.*GRID.dimFactor;                %[plotting dim]
CRmean                  = CRmean*GRID.dimFactor;                %[plotting dim]
CRmedian              	= CRmedian*GRID.dimFactor;             	%[plotting dim]
CRmeanCont          	= CRmeanCont*GRID.dimFactor;         	%[plotting dim]
CRmedianCont        	= CRmedianCont*GRID.dimFactor;        	%[plotting dim]
CRmeanOcean             = CRmeanOcean*GRID.dimFactor;        	%[plotting dim]
CRmedianOcean       	= CRmedianOcean*GRID.dimFactor;       	%[plotting dim]
SAGE_3D              	= SAGE_3D.*(SETUP.timescale/SETUP.secyear/1e6);   	%[Myr]
SAGEmean            	= SAGEmean*SETUP.timescale/SETUP.secyear/1e6;    	%[Myr]
SAGEmedian             	= SAGEmedian*SETUP.timescale/SETUP.secyear/1e6;   	%[Myr]
SAGEmeanCont         	= SAGEmeanCont*SETUP.timescale/SETUP.secyear/1e6;    	%[Myr]
SAGEmedianCont       	= SAGEmedianCont*SETUP.timescale/SETUP.secyear/1e6;   	%[Myr]
SAGEmeanOcean         	= SAGEmeanOcean*SETUP.timescale/SETUP.secyear/1e6;    	%[Myr]
SAGEmedianOcean      	= SAGEmedianOcean*SETUP.timescale/SETUP.secyear/1e6;   	%[Myr]
if strcmp(GRID.Type,'yinyang')
    T_3Dyang           = T_3Dyang .*SETUP.Tscale;          	%[K]
    VX_3Dyang          = VX_3Dyang .*SETUP.Vscale;        	%[cm/a]
    VY_3Dyang          = VY_3Dyang .*SETUP.Vscale;         	%[cm/a]
    VZ_3Dyang          = VZ_3Dyang .*SETUP.Vscale;          	%[cm/a]
    V_3Dyang         	= V_3Dyang .*SETUP.Vscale;           	%[cm/a]
    ETA_3Dyang         = ETA_3Dyang .*SETUP.etascale;       	%[Pas]
    RHO_3Dyang         = RHO_3Dyang .*SETUP.rhoscale;        	%[kg/m^3]
    STR_3Dyang         = STR_3Dyang .*SETUP.stressscale;    	%[MPa]
    EDOT_3Dyang        = EDOT_3Dyang .*SETUP.edotscale;     	%[1/s]
    VD_3Dyang          = VD_3Dyang .*SETUP.dissscale;       	%[W/m]
    VDnoBasalt_3Dyang 	= VDnoBasalt_3Dyang .*SETUP.dissscale;	%[W/m]
    CR_3Dyang       	= CR_3Dyang.*GRID.dimFactor;         	%[plotting dim]
    SAGE_3Dyang      	= SAGE_3Dyang.*(SETUP.timescale/SETUP.secyear/1e6);	%[Myr]
end
if SWITCH.DimensionalMode || SWITCH.DimensionalInput
    if ~strcmp(GRID.Xdim,'km'); warning(['f_PlateThickness input should be in [',GRID.Xdim,']!']); end
end

%% SOME OUTPUT DIAGNOSTICS
PLATE.CrustalThicknessMean     	= CRmean;           %[plotting dim]
PLATE.CrustalThicknessMedian  	= CRmedian;         %[plotting dim]
PLATE.CrustalThicknessMeanCont	= CRmeanCont;      	%[plotting dim]
PLATE.CrustalThicknessMedianCont= CRmedianCont;   	%[plotting dim]
PLATE.CrustalThicknessMeanOcean	= CRmeanOcean;     	%[plotting dim]
PLATE.CrustalThicknessMedianOcean = CRmedianOcean; 	%[plotting dim]
PLATE.SurfaceAgeMean            = SAGEmean;         %[Ma]
PLATE.SurfaceAgeMedian          = SAGEmedian;       %[Ma]
PLATE.SurfaceAgeMeanCont     	= SAGEmeanCont;    	%[Ma]
PLATE.SurfaceAgeMedianCont     	= SAGEmedianCont; 	%[Ma]
PLATE.SurfaceAgeMeanOcean     	= SAGEmeanOcean;   	%[Ma]
PLATE.SurfaceAgeMedianOcean   	= SAGEmedianOcean; 	%[Ma]

%% GET DIFFERENT PLATE THICKNESSES (using various definitions)
check1 = zeros(nb,1); check2 = check1;
for iLithoDefinition=1:4
    if iLithoDefinition==1 %1600 K (i.e., plate base)
        LITHO.Tisovalue     = LithoBaseTemp1;
    elseif iLithoDefinition==2 %plate core
        LITHO.Tisovalue     = PLATE.BaseTemperature-(PLATE.BaseTemperature-Tmin)*(1-coreFractionLithoDepth); %~2/3 of the plate-base temperature
        LITHO.gridpoints    = NaN; %define here to get indices too
    elseif iLithoDefinition==3 %plate surface (i.e., close to plate surface)
        LITHO.Tisovalue     = PLATE.BaseTemperature-(PLATE.BaseTemperature-Tmin)*(1-surfaceFractionLithoDepth); %~0.05 of the plate-base temperature
        LITHO.gridpoints    = NaN; %define here to get indices too
    else %1300 K (i.e., conservative plate base for e.g., stagnant-lid check)
        LITHO.Tisovalue     = PLATE.BaseTemperature;
        LITHO.gridpoints    = NaN; %define here to get indices too
    end
    for ib=1:nb
        if ib==1
            dummy = T_3D;
        else
            dummy = T_3Dyang;
            LITHOyin = LITHO;
        end
        if SWITCH.DimensionalMode || SWITCH.DimensionalInput
            Z_3Dx = GRID.Z_3D/1e3; %to [km]
        else
            Z_3Dx = GRID.Z_3Dnd; %[nd]
        end
        LITHO.plot      = false;
        [LITHO] = f_PlateThickness(dummy,Z_3Dx,LITHO);
        if iLithoDefinition==1 %1600 K
            if ib==1
                lithoThicknessT_1600        = LITHO.thicknessFull; %scientific plate thickness definition
            else
                lithoThicknessT_1600yang    = LITHO.thicknessFull; %scientific plate thickness definition
            end
        elseif iLithoDefinition==2 %plate core
            if ib==1
                lithoThicknessCore          = LITHO.thicknessFull; %plate core thickness
                gridpointsPlateCore         = LITHO.gridpointsFull;
            else
                lithoThicknessCoreyang      = LITHO.thicknessFull; %plate core thickness
                gridpointsPlateCoreyang     = LITHO.gridpointsFull;
            end
            if ib==nb
                LITHO   = rmfield(LITHO,'gridpoints');
                LITHO   = rmfield(LITHO,'gridpointsFull');
            end
        elseif iLithoDefinition==3 %plate surface
            if ib==1
                lithoThicknessSurface     	= LITHO.thicknessFull; %plate surface thickness
                gridpointsPlateSurface   	= LITHO.gridpointsFull;
            else
                lithoThicknessSurfaceyang 	= LITHO.thicknessFull; %plate surface thickness
                gridpointsPlateSurfaceyang 	= LITHO.gridpointsFull;
            end
            if ib==nb
                LITHO   = rmfield(LITHO,'gridpoints');
                LITHO   = rmfield(LITHO,'gridpointsFull');
            end
        else %1300 K
            if ib==1
                lithoThicknessT_1300        = LITHO.thicknessFull; %necessary for some diagnostics (e.g., stagnant lid)
            else
                lithoThicknessT_1300yang    = LITHO.thicknessFull; %necessary for some diagnostics (e.g., stagnant lid)
            end
        end
        if iLithoDefinition==1 %use T=1600 for stagnant lid check
            %% STAGNANT-LID CHECK I: by plate thickness variation
            checkStagnantLid = logical(1);
            if checkStagnantLid %by plate thickness
                if max(isnan(LITHO.thickness(:)))==0 %no subduction i.e., stagnant lid
                    check1(ib) = 1;
                    if ib==nb && min(check1)==1
                        if SWITCH.Verbose; disp(['   ',STYLE.SCHAR.indicationRightArrow,' stagnant-lid detected - thickness does not exceed ',num2str(LITHO.zmax),': no boundary indication.']); end %for first field plot
                        PLATE.StagnantLid = true;
                    end
                end
                if max(LITHO.thicknessFull(:))-min(LITHO.thicknessFull(:))<mean(LITHO.thickness(:))/2 %half the plate thickness
                    check2(ib) = 1;
                    if ib==nb && min(check2)==1
                        if SWITCH.Verbose; disp(['   ',STYLE.SCHAR.indicationRightArrow,' stagnant-lid detected - low plate-thickness variation: no boundary indication.']); end %for first field plot
                        PLATE.StagnantLid = true;
                    end
                end
            end
        end
    end
end

%% CHECK IF PLATE EXISTS
if min(LITHO.thicknessFull(:))==max(LITHO.thicknessFull(:)) && max(LITHO.thicknessFull(:))<=LITHO.zmin
    PLATE.Subduction = 'noTracking'; %flag for no indication
    
    return %exit function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end

%% PROBLEM CHECK
if min(isnan(LITHO.thickness(:)))==1 %all NaN values
    if SWITCH.Verbose; warning(['LITHO.zmax = ',num2str(LITHO.zmax),' too shallow: consider increasing! - Plate diagnostics cannot be performed.']); end
    PLATE.Subduction = 'noTracking'; %flag for no indication
    
    return %exit function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end

%% CRITICAL INPUT VALUES
if strcmp(GRID.Type,'yinyang')
    depth2check         = 2.3	*mean([LITHOyin.thicknessFull(:);LITHO.thicknessFull(:)]); %adjust here the critical value for checking for subduction (i.e., cold material) in the mantle
    criticalT2check     = 4     *min([T_3D(:);T_3Dyang(:)]+0.001); %4; adjust critical temperature of the cold sinking material to find
    criticalT2checkDeep	= 4.5  	*min([T_3D(:);T_3Dyang(:)]+0.001); %4; adjust critical temperature of the cold sinking material to find
else
    % depth2check         = 3 	 *mean(LITHO.thickness(~isnan(LITHO.thickness(:)))); %3; adjust here the critical value for checking for subduction (i.e., cold material) in the mantle
    depth2check         = 2.3	*mean(LITHO.thicknessFull(:)); %1.5; adjust here the critical value for checking for subduction (i.e., cold material) in the mantle
    criticalT2check     = 4     *min(T_3D(:)+0.001); %4; adjust critical temperature of the cold sinking material to find
    criticalT2checkDeep = 4.5  	*min(T_3D(:)+0.001); %4; adjust critical temperature of the cold sinking material to find
end
deeperFac               = 4/6;                                      %deep-slab level
deeperFacUM             = 3/4;                                      %upper-mantle level

dz2deeperSlab           = deeperFac*depth2check;                    %distance to deeper checking points
dz2upperMantle          = deeperFacUM*depth2check;                  %distance to deeper checking points
depth2check2            = depth2check+ dz2deeperSlab;               %get a deeper check too (for sub. polarity check)
depth2checkUM           = depth2check+ dz2upperMantle;              %get a deeper check too (for sub. polarity check)

[~,zIdxCheck]           = min(abs(depth2check-Z_3Dx(1,1,:)));       %[nd] or [km]
[~,zIdxCheck2]          = min(abs(depth2check2-Z_3Dx(1,1,:)));      %[nd] or [km] %for deeper slab level
[~,zIdxCheckUM]         = min(abs(depth2checkUM-Z_3Dx(1,1,:)));     %[nd] or [km] %for upper mantle diagnostics level
if zIdxCheck<=zIdxCheck2+3; zIdxCheck2 = zIdxCheck-4; end       %make sure they are not the same (at least 5 grid points separated)
if zIdxCheck<=zIdxCheckUM+5; zIdxCheckUM = zIdxCheck-6; end    	%make sure they are not the same (at least 5 grid points separated)
zIdxCheckMid            = min(zIdxCheck,zIdxCheck2)+round(abs(zIdxCheck-zIdxCheck2)/2); 	%index in between measurement points
dnx2deeperSlab          = 7 *abs(zIdxCheck-zIdxCheck2);             %horizontal distance (in #gridpoints) trench<->slab to check

if logical(0) %check depth levels
    hold on
    plot(GRID.X_3Dp(:,:,zIdxCheck),GRID.Z_3Dp(:,:,zIdxCheck),'-k')
    plot(GRID.X_3Dp(:,:,zIdxCheck2),GRID.Z_3Dp(:,:,zIdxCheck2),'-b')
    plot(GRID.X_3Dp(:,:,zIdxCheckUM),GRID.Z_3Dp(:,:,zIdxCheckUM),'-r')
end

%% PROBLEM CHECKS
%check for too shallow check points
if zIdxCheck2<=0
    if SWITCH.Verbose; warning('check points too shallow: no plate diagnostics possible.'); end
    PLATE.Subduction = 'noTracking'; %flag for no indication
    
    return %exit function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
end

%% GET PLATE-CORE DEPTH THROUGHOUT MODEL DOMAIN
%indices become lower towards deeper levels!
[~,idxMeanLabDepth4core]	= min(abs(GRID.Z_3Dp(1,1,:)-LABdepth*coreFractionLithoDepth));  %index of closest value to mean LAB depth further up in the plate
[~,idxSeaLevel]             = min(abs(GRID.Z_3Dp(1,1,:)));  %index of closest value to sealevel inside the plate
if GRID.Z_3Dp(1,1,idxSeaLevel)<0; idxSeaLevel = idxSeaLevel-1; end %make sure it is below the surface
% closest = GRID.Z_3Dp(1,1,idxMeanLabDepth); %closest value
PLATE.coreLevelIdx          = min(max(gridpointsPlateCore,idxMeanLabDepth4core),idxSeaLevel);       %limit max depth to approximate plate thickness and minimum to sealevel (for plate core)
PLATE.surfaceLevelIdx       = min(max(gridpointsPlateSurface,idxMeanLabDepth4core),idxSeaLevel);    %limit max depth to approximate plate thickness and minimum to sealevel (for plate surface)
surfacePlateIdx             = min(max(LITHO.gridpointsFull,idxMeanLabDepth4core),idxSeaLevel);      %limit max depth to approximate plate thickness and minimum to sealevel (for 1300 K plate base)
surfaceAndBendingPlateIdx   = min(LITHO.gridpointsFull,idxSeaLevel);                                %limit minimum to sealevel
if GRID.Z_3D(1,1,1)<GRID.Z_3D(1,1,2); error('GRID.Z_3D is increasing with higher index. The above 3 lines need to be adjusted!'); end

%this is for indices that can be used to extract data points from 3-D arrays efficiently
[~,surfaceIdx]              = min(abs(GRID.Z_3D(1,1,:))); %find closest index to 0 depth (i.e., sea) level
coreLevelDummy              = reshape(1:numel(PLATE.coreLevelIdx),size(PLATE.coreLevelIdx))+numel(PLATE.coreLevelIdx)*(PLATE.coreLevelIdx-1);
surfaceLevelDummy        	= reshape(1:numel(PLATE.surfaceLevelIdx),size(PLATE.surfaceLevelIdx))+numel(PLATE.surfaceLevelIdx)*(PLATE.surfaceLevelIdx-1);
% surfacePlateDummy           = reshape(1:numel(surfacePlateIdx),size(surfacePlateIdx))+numel(surfacePlateIdx)*(surfacePlateIdx-1);
% surfaceAndBendingPlateDummy = reshape(1:numel(surfaceAndBendingPlateIdx),size(surfaceAndBendingPlateIdx))+numel(surfaceAndBendingPlateIdx)*(surfaceAndBendingPlateIdx-1);
%and then use:
% a               	= GRID.Z_3D(coreLevelDummy);

%this is for 0 and 1 matrices (e.g., fill plate with 1s) to extract data from 3-D arrays efficiently
surfacePlateDummyM = zeros(size(T_3D));
surfaceAndBendingPlateDummyM = surfacePlateDummyM;
for ix=1:size(PLATE.coreLevelIdx,1)
    for iy=1:size(PLATE.coreLevelIdx,2) %loop through vertical columns
        surfacePlateDummyM(ix,iy,surfacePlateIdx(ix,iy):surfaceIdx)     = 1;  %fill plate with 1 indices
        surfaceAndBendingPlateDummyM(ix,iy,surfaceAndBendingPlateIdx(ix,iy):surfaceIdx)	= 1;  %fill plate with 1 indices
    end
end

%% GET VARIOUS PARAMETERS FROM INSIDE THE PLATE AND ALONG THE PLATE CORE
dummy = zeros(size(PLATE.coreLevelIdx,1),size(PLATE.coreLevelIdx,2));
dummyNaN = nan(size(PLATE.coreLevelIdx,1),size(PLATE.coreLevelIdx,2));
if ~isnan(STR_3D(1));           PLATE.maxStress = dummy; PLATE.coreStress = dummy; end
if ~isnan(EDOT_3D(1));          PLATE.coreStrainrate = dummy; end
if ~isnan(ETA_3D(1));           PLATE.coreViscosity = dummy; end
if ~isnan(VX_3D(1));            PLATE.coreXVelocity = dummy; else; PLATE.coreXVelocity = dummyNaN; end
if ~isnan(VY_3D(1));            PLATE.coreYVelocity = dummy; else; PLATE.coreYVelocity = dummyNaN; end
if ~isnan(VX_3D(1));            PLATE.surfXVelocity = dummy; else; PLATE.surfXVelocity = dummyNaN; end
if ~isnan(VY_3D(1));            PLATE.surfYVelocity = dummy; else; PLATE.surfYVelocity = dummyNaN; end
if strcmp(GRID.Type,'yinyang')
    if ~isnan(VX_3D(1));     	PLATE.coreXVelocity_yang = dummy; else; PLATE.coreXVelocity_yang = dummyNaN; end
    if ~isnan(VY_3D(1));      	PLATE.coreYVelocity_yang = dummy; else; PLATE.coreYVelocity_yang = dummyNaN; end
    if ~isnan(VX_3D(1));    	PLATE.surfXVelocity_yang = dummy; else; PLATE.surfXVelocity_yang = dummyNaN; end
    if ~isnan(VY_3D(1));      	PLATE.surfYVelocity_yang = dummy; else; PLATE.surfYVelocity_yang = dummyNaN; end
end
if ~isnan(VD_3D(1));            sumViscDissipation = dummy; else; sumViscDissipation = NaN; end
if ~isnan(VDnoBasalt_3D(1));    sumViscDissipation2 = dummy; else; sumViscDissipation2 = NaN; end
if logical(0) %slow version
    for ix=1:size(PLATE.coreLevelIdx,1)
        for iy=1:size(PLATE.coreLevelIdx,2) %loop through vertical columns and get one value for each
            if ~isnan(STR_3D)
                PLATE.coreStress(ix,iy) 	= STR_3D(ix,iy,PLATE.coreLevelIdx(ix,iy));
                PLATE.maxStress(ix,iy)      = max(STR_3D(ix,iy,surfacePlateIdx(ix,iy):end));
            end
            if ~isnan(EDOT_3D(1));       	PLATE.coreStrainrate(ix,iy)	= EDOT_3D(ix,iy,PLATE.coreLevelIdx(ix,iy)); end
            if ~isnan(ETA_3D(1));       	PLATE.coreViscosity(ix,iy) 	= ETA_3D(ix,iy,PLATE.coreLevelIdx(ix,iy));  end
            if ~isnan(VX_3D(1));         	PLATE.coreXVelocity(ix,iy) 	= VX_3D(ix,iy,PLATE.coreLevelIdx(ix,iy));  end
            if ~isnan(VY_3D(1));         	PLATE.coreYVelocity(ix,iy) 	= VY_3D(ix,iy,PLATE.coreLevelIdx(ix,iy));  end
            if ~isnan(VX_3D(1));         	PLATE.surfXVelocity(ix,iy) 	= VX_3D(ix,iy,PLATE.surfaceLevelIdx(ix,iy));  end
            if ~isnan(VY_3D(1));         	PLATE.surfYVelocity(ix,iy) 	= VY_3D(ix,iy,PLATE.surfaceLevelIdx(ix,iy));  end
            if strcmp(GRID.Type,'yinyang')
                if ~isnan(VX_3D(1));        PLATE.coreXVelocity_yang(ix,iy) = VX_3Dyang(ix,iy,PLATE.coreLevelIdx(ix,iy));  end
                if ~isnan(VY_3D(1));     	PLATE.coreYVelocity_yang(ix,iy) = VY_3Dyang(ix,iy,PLATE.coreLevelIdx(ix,iy));  end
                if ~isnan(VX_3D(1));      	PLATE.surfXVelocity_yang(ix,iy) = VX_3Dyang(ix,iy,PLATE.surfaceLevelIdx(ix,iy));  end
                if ~isnan(VY_3D(1));     	PLATE.surfYVelocity_yang(ix,iy) = VY_3Dyang(ix,iy,PLATE.surfaceLevelIdx(ix,iy));  end
            end
            if ~isnan(VD_3D(1));        	sumViscDissipation(ix,iy) 	= sum(VD_3D(ix,iy,surfaceAndBendingPlateIdx(ix,iy):end)); end
            if ~isnan(VDnoBasalt_3D(1));  	sumViscDissipation2(ix,iy)  = sum(VDnoBasalt_3D(ix,iy,surfaceAndBendingPlateIdx(ix,iy):end)); end
        end
    end
else %faster version
    if ~isnan(STR_3D)
        PLATE.coreStress        = STR_3D(coreLevelDummy);
        PLATE.maxStress         = max(STR_3D(surfacePlateDummyM==1));
    end
    if ~isnan(EDOT_3D(1));       	PLATE.coreStrainrate	= EDOT_3D(coreLevelDummy); end
    if ~isnan(ETA_3D(1));       	PLATE.coreViscosity 	= ETA_3D(coreLevelDummy); end
    if ~isnan(VX_3D(1));         	PLATE.coreXVelocity 	= VX_3D(coreLevelDummy); end
    if ~isnan(VY_3D(1));         	PLATE.coreYVelocity 	= VY_3D(coreLevelDummy); end
    if ~isnan(VX_3D(1));         	PLATE.surfXVelocity 	= VX_3D(surfaceLevelDummy); end
    if ~isnan(VY_3D(1));         	PLATE.surfYVelocity 	= VY_3D(surfaceLevelDummy); end
    if strcmp(GRID.Type,'yinyang')
        if ~isnan(VX_3D(1));     	PLATE.coreXVelocity_yang 	= VX_3Dyang(coreLevelDummy); end
        if ~isnan(VY_3D(1));        PLATE.coreYVelocity_yang 	= VY_3Dyang(coreLevelDummy); end
        if ~isnan(VX_3D(1));        PLATE.surfXVelocity_yang 	= VX_3Dyang(surfaceLevelDummy); end
        if ~isnan(VY_3D(1));    	PLATE.surfYVelocity_yang 	= VY_3Dyang(surfaceLevelDummy); end
    end
    if ~isnan(VD_3D(1));        	sumViscDissipation      = sum(VD_3D(surfaceAndBendingPlateDummyM==1)); end
    if ~isnan(VDnoBasalt_3D(1));  	sumViscDissipation2     = sum(VDnoBasalt_3D(surfaceAndBendingPlateDummyM==1)); end
end

if strcmp(GRID.Type,'yinyang')
    warning('PLATE.coreLevelIdx and subsequent diagnostics are not yet implemented for YinYang geometry!')
end
PLATE.coreStressMin                 = min(PLATE.coreStress(:));
PLATE.coreStressMax                 = max(PLATE.coreStress(:));
PLATE.maxStress                     = max(PLATE.maxStress(:));
PLATE.coreStrainrateMin             = min(PLATE.coreStrainrate(:));
PLATE.coreStrainrateMax             = max(PLATE.coreStrainrate(:));
PLATE.coreViscosityMin              = min(PLATE.coreViscosity(:));
PLATE.coreViscosityMax              = max(PLATE.coreViscosity(:));

if strcmp(GRID.Dim,'2-D')
    PlateCoreHVelocity            	= PLATE.coreXVelocity;
else
    PlateCoreHVelocity            	= sqrt( PLATE.coreXVelocity.^2 + PLATE.coreYVelocity.^2 );  %results in absolute values
end
PLATE.PlateVelocity                 = PlateCoreHVelocity;
if strcmp(GRID.Type,'yinyang')
    PlateCoreHVelocity_yang        	= sqrt( PLATE.coreXVelocity_yang.^2 + PLATE.coreYVelocity_yang.^2 );
    PLATE.PlateVelocity_yang       	= PlateCoreHVelocity_yang;
end

PLATE.ViscDissipationPlate          = sum(sumViscDissipation(:)); 	%integral over whole surface plate plus bending area at subduction zones (if existing) including subduction channel
PLATE.ViscDissipationPlateNoCrust	= sum(sumViscDissipation2(:)); 	%same but neglecting areas with crust

if strcmp(GRID.Type,'yinyang')
    PLATE.PlateVelocityMax              = max([PlateCoreHVelocity(:);PlateCoreHVelocity_yang(:)]);
    PLATE.PlateVelocityMin              = min([PlateCoreHVelocity(:);PlateCoreHVelocity_yang(:)]);
    PLATE.PlateVelocityDiff             = sqrt((max([PLATE.coreXVelocity(:);PLATE.coreXVelocity_yang(:)])-min([PLATE.coreXVelocity(:);PLATE.coreXVelocity_yang(:)]))^2 +...
        (max([PLATE.coreYVelocity(:);PLATE.coreYVelocity_yang(:)])-min([PLATE.coreYVelocity(:);PLATE.coreYVelocity_yang(:)]))^2);
    PLATE.PlateVelocityRMS              = rms([PlateCoreHVelocity(:);PlateCoreHVelocity_yang(:)]);
    maxDiffPlateSurfVel                 = sqrt((max([PLATE.surfXVelocity(:);PLATE.surfXVelocity_yang(:)])-min([PLATE.surfXVelocity(:);PLATE.surfXVelocity_yang(:)]))^2 +...
        (max([PLATE.surfYVelocity(:);PLATE.surfYVelocity_yang(:)])-min([PLATE.surfYVelocity(:);PLATE.surfYVelocity_yang(:)]))^2);
    %DOES CURRENTLY NOT AVOID OVERLAPPING YY-EDGES.... TRY SOMETHING LIKE THIS BELOW: ........
    %         %try to avoid including overlapping yy-edges
    %         cutAway             = 50; %no. grid points to cut away
    %         maxDiffPlateVel     = max( abs(max(max(PlateVelx(:,cutAway+1:end-cutAway)))-min(min(PlateVelx(:,cutAway+1:end-cutAway)))),...
    %             abs(max(max(PlateVely(:,cutAway+1:end-cutAway)))-min(min(PlateVely(:,cutAway+1:end-cutAway)))) ); %THIS IS NOT PERFECT YET (as oblique motion is not accounted for)
else
    PLATE.PlateVelocityMax              = max(PlateCoreHVelocity(:));
    PLATE.PlateVelocityMin              = min(PlateCoreHVelocity(:));
    PLATE.PlateVelocityDiff             = sqrt((max(PLATE.coreXVelocity(:))-min(PLATE.coreXVelocity(:)))^2 + (max(PLATE.coreYVelocity(:))-min(PLATE.coreYVelocity(:)))^2);
    PLATE.PlateVelocityRMS              = rms(PlateCoreHVelocity(:));
    maxDiffPlateSurfVel                 = sqrt((max(PLATE.surfXVelocity(:))-min(PLATE.surfXVelocity(:)))^2 + (max(PLATE.surfYVelocity(:))-min(PLATE.surfYVelocity(:)))^2);
end

%% VARIABLE ADJUSTMENTS
minSurfaceLevel      	= min(Z_3D(surfaceLevelDummy));
maxSurfaceLevel      	= max(Z_3D(surfaceLevelDummy));
PlateVelx               = PLATE.coreXVelocity;
PlateSurfVelx           = PLATE.surfXVelocity;
PlateVely               = PLATE.coreYVelocity;
PlateSurfVely           = PLATE.surfYVelocity;
if strcmp(GRID.Type,'yinyang')
    PlateVelx_yang    	= PLATE.coreXVelocity_yang;
    PlateSurfVelx_yang 	= PLATE.surfXVelocity_yang;
    PlateVely_yang     	= PLATE.coreYVelocity_yang;
    PlateSurfVely_yang 	= PLATE.surfYVelocity_yang;
end
PLATE   = rmfield(PLATE,'coreStress');
PLATE   = rmfield(PLATE,'coreStrainrate');
PLATE   = rmfield(PLATE,'coreViscosity');
clearvars sumViscDissipation sumViscDissipation2

%% GET VARIOUS PARAMETER OF THE MANTLE (not adjusted for YY mode.......)
%general detection if no subduction zone is found
idxMeanLabDepth4UM          = zIdxCheckUM;
% mantleFractionLithoDepth    = 1.5;      %fraction of plate thickness to check for upper mantle parameters <<<<<<<<<<<<<<<<<<<<<<<<<<<
% [~,idxMeanLabDepth4UM]	= min(abs(GRID.Z_3Dp(1,1,:)-LABdepth*mantleFractionLithoDepth)); %index of closest value to mean LAB depth further down in the mantle
dummy                       = ETA_3D(:,:,idxMeanLabDepth4UM);
PLATE.UMantleViscosity      = median(dummy(:));
dummy                       = sqrt(VX_3D(:,:,idxMeanLabDepth4UM).^2+VY_3D(:,:,idxMeanLabDepth4UM).^2+VZ_3D(:,:,idxMeanLabDepth4UM).^2);
PLATE.UMantleVelocityMean  	= mean(dummy(:));
PLATE.UMantleVelocityMax  	= max(dummy(:));
dummy                       = sqrt(VX_3D(:,:,idxMeanLabDepth4UM).^2+VY_3D(:,:,idxMeanLabDepth4UM).^2);
PLATE.UMantleHVelocityMean 	= mean(dummy(:));
PLATE.UMantleHVelocityMax 	= max(dummy(:));
dummy                       = VZ_3D(:,:,idxMeanLabDepth4UM);
PLATE.UMantleRVelocityMax   = max(dummy(:));
dummy                       = T_3D(:,:,idxMeanLabDepth4UM);
PLATE.UMantleTemperature    = median(dummy(:));
dummy                       = RHO_3D(:,:,idxMeanLabDepth4UM);
PLATE.UMantleDensity        = median(dummy(:));


%% DIAGNOSE PLATE MOBILITY
PLATE.PlateMobility         = PLATE.PlateVelocityRMS / PLATE.UMantleVelocityMean;
PLATE.PlateDrivety          = PLATE.PlateVelocityRMS / PLATE.UMantleHVelocityMean;  %is plate driven by itself or by the mantle



%% STAGNANT-LID DIAGNOSTICS
if PLATE.StagnantLid
    LID.thickness           = lithoThicknessT_1300; %array of horizontal thickness
    LID.LABdepth            = LABdepth; %one value
    LID.LABidx              = LABidx; %one value
    if strcmp(GRID.Type,'yinyang'); LID.thicknessYang = lithoThicknessT_1300yang; end
    [LID,PLOT] = f_stagnantLidDiagnostics(FILE,GRID,LID,SWITCH,PLOT); %LID.maxYieldDepth = NaN if not successful
    
    PLATE.maxYieldDepth             = LID.maxYieldDepth;
    PLATE.maxYieldDepthFraction     = LID.maxYieldDepthFraction;
    %if you add new variables, also add them as NaN on top!

    %DISPLAY for stagnant lid
    disp('   Stagnant Lid')
    disp('     Viscosity');
    disp(['     ',STYLE.SCHAR.smallBullet,' plate core            = ',num2str(PLATE.coreViscosityMin,2),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLATE.coreViscosityMax,2),' ',SETUP.etaDim]);
    disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(PLATE.UMantleViscosity',2),' ',SETUP.etaDim,' (',num2str(GRID.Z_3Dp(1,1,zIdxCheck2),3),' ',GRID.Zdim,' depth)']);
    disp('     Stress');
    disp(['     ',STYLE.SCHAR.smallBullet,' plate max.            = ',num2str(PLATE.maxStress,3),' ',SETUP.stressDim]);
    disp(['     ',STYLE.SCHAR.smallBullet,' plate core            = ',num2str(PLATE.coreStressMin,3),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLATE.coreStressMax,3),' ',SETUP.stressDim]);
    disp('     Strainrate');
    disp(['     ',STYLE.SCHAR.smallBullet,' plate core            = ',num2str(PLATE.coreStrainrateMin,2),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLATE.coreStrainrateMax,2),' ',SETUP.edotDim2]);
    disp('     Velocity');
    disp(['     ',STYLE.SCHAR.smallBullet,' plate (RMS)           = ',num2str(PLATE.PlateVelocityRMS',2),' ',SETUP.vDim,' (',num2str(PLATE.PlateVelocityMin,2),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLATE.PlateVelocityMax,2),' ',SETUP.vDim,')']);
    disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle (mean)   = ',num2str(PLATE.UMantleVelocityMean,2),' ',SETUP.vDim,' (',num2str(GRID.Z_3Dp(1,1,zIdxCheck2),3),' ',GRID.Zdim,' depth)']);
    disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle (max)    = ',num2str(PLATE.UMantleVelocityMax,2),' (total), ',num2str(PLATE.UMantleHVelocityMax,2),' (horizontal), ',num2str(PLATE.UMantleRVelocityMax,2),' (radial) '...
        ,SETUP.vDim,' (',num2str(GRID.Z_3Dp(1,1,zIdxCheck2),3),' ',GRID.Xdim,' depth)']);
%     disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(PLATE.UMantleVelocityMax,2),' ',SETUP.vDim,' (',num2str(GRID.Z_3Dp(1,1,zIdxCheck2),3),' ',GRID.Xdim,' depth)']);
%     disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle (horiz.)	= ',num2str(PLATE.UMantleHVelocityMax,2),' ',SETUP.vDim,' (',num2str(GRID.Z_3Dp(1,1,zIdxCheck2),3),' ',GRID.Xdim,' depth)']);
%     disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle (radial)	= ',num2str(PLATE.UMantleRVelocityMax,2),' ',SETUP.vDim,' (',num2str(GRID.Z_3Dp(1,1,zIdxCheck2),3),' ',GRID.Xdim,' depth)']);
    if ~isnan(LID.maxYieldDepth)
        disp('     Other');
        disp(['     ',STYLE.SCHAR.smallBullet,' max. yield depth      = ',num2str(LID.maxYieldDepth,2),' ',GRID.Zdim,' (',num2str(LID.maxYieldDepthFraction*100,2),'% of lid)'])
    end
    PLATE.Subduction = 'noTracking'; %flag for no indication
    
    return %exit function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
end

%% FIND SUBDUCTION IN THE MANTLE  
%is there at least one slab present in the shallow mantle?
check1 = zeros(nb,1);
for ib=1:nb
    if ib==1
        dummy = T_3D;
    else
        dummy = T_3Dyang;
        Thorizon2check_yin  = Thorizon2check;
        Thorizon2check2_yin = Thorizon2check2;
    end
    Thorizon2check          = dummy(:,:,zIdxCheck);
    Thorizon2check2         = dummy(:,:,zIdxCheck2);
    minTfound               = min(Thorizon2check(:));
    if minTfound>criticalT2check
        check1(ib) = 1;
        if ib==nb && min(check1)==1
            %subduction not found
            if SWITCH.Verbose; disp(['   ',STYLE.SCHAR.indicationRightArrow,' stagnant-lid detected - no cold deep slab detected: no boundary indication.']); end %for first field plot
            PLATE.StagnantLid = true;
            PLATE.Subduction = 'noTracking'; %flag for no indication
            
            return %exit function
        end
    end
end
% .......
% to find number of slabs in domain
%use here CC = bwconncomp(BW) and numConnAreas	= cellfun(@numel,CC.NumObjects);
%and even better to account for periodic sides use f_connectivity_2D .....
% .......
if strcmp(GRID.Type,'yinyang')
    coldestSlab       	= Thorizon2check_yin<criticalT2check;
    coldestSlab_yang   	= Thorizon2check<criticalT2check;
    coldestSlab2       	= Thorizon2check2_yin<criticalT2checkDeep; %NOT USED YET 3-D (THIS IS THE EVEN DEEPER SLAB)
    coldestSlab_yang2  	= Thorizon2check2<criticalT2checkDeep;  %NOT USED YET 3-D
    coldSlabs           = coldestSlab;          %shallower
    coldSlabs_yang    	= coldestSlab_yang;   	%shallower
    coldSlabs2          = coldestSlab2;         %deeper
    coldSlabs_yang2    	= coldestSlab_yang2;  	%deeper
    
else %if strcmp(GRID.Dim,'3-D')
    coldestSlab         = Thorizon2check<criticalT2check;      %shallower
    coldestSlab2        = Thorizon2check2<criticalT2checkDeep; %deeper (NOT USED YET in 3-D)
    coldSlabs           = coldestSlab;      %shallower
    coldSlabs2          = coldestSlab2;     %deeper
    if strcmp(GRID.Dim,'2-D')
        if max(coldestSlab(:))==0 %no shallow slab found
            coldestSlab 	= NaN;
        else
            coldestSlab  	= Thorizon2check==min(Thorizon2check(coldestSlab)); %find minimum temperature
            coldestSlab    	= bwmorph(coldestSlab,'shrink',Inf);  %shrink muliple points to one center point
            coldSlabs       = bwmorph(coldSlabs,'shrink',Inf);  %shrink muliple points to one center point
        end
        if max(coldestSlab2(:))==0 %no deep slab found
            coldestSlab2 	= NaN;
        else
            coldestSlab2   	= Thorizon2check2==min(Thorizon2check2(coldestSlab2)); %find minimum temperature
            coldestSlab2  	= bwmorph(coldestSlab2,'shrink',Inf); %shrink muliple points to one center point
            coldSlabs2     	= bwmorph(coldSlabs2,'shrink',Inf);  %shrink muliple points to one center point
        end
        
        numSlabs        = sum(coldSlabs);   %shallower
        numSlabs2     	= sum(coldSlabs2);  %deeper
        
        if numSlabs~=numSlabs2
           warning('number of shallow slabs varies from number of deep slabs.'); 
        end
        
        %% FIND SLAB POLARITY 1st
        slabPolarity = zeros(size(coldestSlab));
        missingSlabDeep = zeros(size(coldestSlab));
        for ix=1:size(coldestSlab,1)
            if coldestSlab(ix,1)==1 %slab found
                if isnan(coldestSlab2) %no deep slab found at all
                    slabPolarity(ix,1)      = 0;
                else
                    if coldestSlab2(ix,1)==1
                        slabPolarity(ix,1)      = 0; %vertical slab
                    elseif max(coldestSlab2(max(1,ix-dnx2deeperSlab):ix,1))==1 && ... %DOES NOT ACCOUNT FOR PERIODIC SIDES
                            max(coldestSlab2(ix:min(size(coldestSlab,1),ix+dnx2deeperSlab),1))==0
                        slabPolarity(ix,1)      = -1; %LEFT
                    elseif max(coldestSlab2(ix:min(size(coldestSlab,1),ix+dnx2deeperSlab),1)) && ... %DOES NOT ACCOUNT FOR PERIODIC SIDES
                            max(coldestSlab2(max(1,ix-dnx2deeperSlab):ix,1))==0
                        slabPolarity(ix,1)      = +1; %RIGHT
                    else %either when nothing is found or when found on both sides
                        slabPolarity(ix,1)      = 0;
                    end
                end
            end
            
            if coldSlabs(ix,1)==1 %slab found
                if isnan(coldSlabs2) %no deep slab found at all
                    slabPolarity(ix,1)      = 0;
                    missingSlabDeep(ix,1) 	= 1;
                    coldSlabs2(ix,1)        = 1; %just temporalily, will be removed later
                else
                    if coldSlabs2(ix,1)==1
                        slabPolarity(ix,1)      = 0; %vertical slab
                    elseif max(coldSlabs2(max(1,ix-dnx2deeperSlab):ix,1))==1 && ... %DOES NOT ACCOUNT FOR PERIODIC SIDES
                            max(coldSlabs2(ix:min(size(coldSlabs,1),ix+dnx2deeperSlab),1))==0
                        slabPolarity(ix,1)      = -1; %LEFT
                    elseif max(coldSlabs2(ix:min(size(coldSlabs,1),ix+dnx2deeperSlab),1)) && ... %DOES NOT ACCOUNT FOR PERIODIC SIDES
                            max(coldSlabs2(max(1,ix-dnx2deeperSlab):ix,1))==0
                        slabPolarity(ix,1)      = +1; %RIGHT
                    else %either when nothing is found or when found on both sides
                        slabPolarity(ix,1)      = 0;
                        missingSlabDeep(ix,1) 	= 1;
                        coldSlabs2(ix,1)        = 1; %just temporarily, will be removed later
                    end
                end
            end

        end
        
        %Some diagnostics and adjustments
        if sum(coldSlabs)~=sum(coldSlabs2)
            warning('number of shallow slabs STILL varies from number of deep slabs.');
        end
        numberSlabsDetected             = sum(coldSlabs);
        if sum(coldSlabs)>1
            multipleSlabsDetected       = true;
        else
            multipleSlabsDetected       = false;
        end
    end
end

if SWITCH.Verbose; warning('REMOVE THAT HERE ONCE FINISHED IMPLEMENTING: cold deeper slabs have beeen added artificially to prevent mismatch with shallower slabs (might improve on that if more deeper slabs are found). Then go on from here to finish implementation.'); end


%% SLAB-TIP DEPTH & VELOCITY
[SLABTIP] = f_SlabTipDiagnostics(T_3D,PLATE,GRID,SWITCH);
PLATE.SlabTipPosition           = SLABTIP.SlabTipPosition;
PLATE.SlabTipAngle              = SLABTIP.SlabTipAngle;
PLATE.SlabTipAnglePoint1        = SLABTIP.SlabTipAnglePoint1;
PLATE.SlabTipAnglePoint2        = SLABTIP.SlabTipAnglePoint2;
if ~isnan(SLABTIP.SlabTipIndicesXZ(1,1)) && ~isnan(VX_3D(1))
    PLATE.SlabTipVX             = VX_3D(SLABTIP.SlabTipIndicesXZ(1,1),1,SLABTIP.SlabTipIndicesXZ(1,2));
    PLATE.SlabTipVZ             = VZ_3D(SLABTIP.SlabTipIndicesXZ(1,1),1,SLABTIP.SlabTipIndicesXZ(1,2));
end

%% SUBDUCTING-PLATE FIT, BENDING RADIUS, SLAB-TIP DEPTH, 2nd CHECK FOR SLAB POLARITY
[BENDING] = f_BendingRadiusB(T_3D,PLATE,GRID,SWITCH,LITHO);
if ~BENDING.failed
    PLATE.subPlateIndices       = BENDING.subPlateIndices;
    PLATE.subPlateFullFit       = BENDING.subPlateFullFit;
    PLATE.subPlateFit           = BENDING.subPlateFit;
    PLATE.bendingCircleCenter   = BENDING.circleCenter;
end

%% SUBDUCTING-PLATE CHARACTERISTICS ....not used yet....
if ~isnan(PLATE.subPlateIndices)
%     PLATE.SubPlateTemperature   = T_3D(PLATE.subPlateIndices);
%     PLATE.SubPlateVelocity      = V_3D(PLATE.subPlateIndices);
%     PLATE.SubPlateStress        = STR_3D(PLATE.subPlateIndices);
%     PLATE.SubPlateStrainrate    = EDOT_3D(PLATE.subPlateIndices);
%     PLATE.SubPlateViscosity     = ETA_3D(PLATE.subPlateIndices);
end

%% FINISH YINYANG DIAGNOSTICS (needs further implementation)
if strcmp(GRID.Type,'yinyang')
    PLATE.Subduction    = 'noTracking'; %flag for no indication
    
    return
end

%% SLAB POSITIONS
clearvars dummy
%shallow slab
if strcmp(GRID.Dim,'2-D')
    xColdSlabs              = GRID.X_3D(coldSlabs,zIdxCheck); 	%[nd] or [m] or [rad]
elseif strcmp(GRID.Dim,'3-D')
    dummy(:,:)              = GRID.X_3D(:,:,zIdxCheck);
    xColdSlabs              = dummy(coldSlabs);                 %[nd] or [m]
    clearvars dummy
end
zColdSlabs                  = GRID.Z_3D(1,1,zIdxCheck);         %[nd] or [m]
%deep slab
if size(coldSlabs2,1)==1 && isnan(coldSlabs2) %no deep slab found
    xColdSlabs2             = NaN;
    zColdSlabs2             = NaN;
else
    if strcmp(GRID.Dim,'2-D')
        xColdSlabs2      	= GRID.X_3D(coldSlabs2,zIdxCheck2);	%[nd] or [m] or [rad]
    elseif strcmp(GRID.Dim,'3-D')
        dummy(:,:)          = GRID.X_3D(:,:,zIdxCheck2);
        xColdSlabs2         = dummy(coldSlabs2);              	%[nd] or [m]
        clearvars dummy
    end
    zColdSlabs2             = GRID.Z_3D(1,1,zIdxCheck2);     	%[nd] or [m]
end

%% SLAB-SINKING VELOCITY
if strcmp(GRID.Dim,'2-D')
    PLATE.SlabSinkingVelocity 	= -VZ_3D(coldestSlab,1,zIdxCheck); %[plotting dimension], +ve means downward
elseif strcmp(GRID.Dim,'3-D')
    dummy                       = -VZ_3D(:,:,zIdxCheck);
    PLATE.SlabSinkingVelocity 	= median( dummy(coldestSlab) );
    clearvars dummy
end

%% SLAB-DIP ANGLE
coldSlabsMid                = round(zColdSlabs2-zColdSlabs);      %point between upper and lower slab levels
if size(xColdSlabs,1)~=size(xColdSlabs2,1) || ... %error check if same #slabs are found at both depth levels
        (size(coldestSlab2,1)==1 && isnan(coldestSlab2)) %or no deep slab found
    angleSlabs              = NaN;
    depthShallowSlabDip   	= NaN;
    PLATE.slabDipPointsX  	= [NaN, NaN];
    PLATE.slabDipPointsZ  	= [NaN, NaN];
else %normal case
    if strcmp(GRID.Type,'Cartesian')
        dx                  = abs(xColdSlabs-xColdSlabs2);      %[nd] or [m]
    elseif strcmp(GRID.Type,'spherical2D')
        dx                  = abs(xColdSlabs-xColdSlabs2)*GRID.R(1,1,zIdxCheckMid);     %L = r*alpha, [nd] or [m]
    end
    dz                      = abs(zColdSlabs-zColdSlabs2);      %[nd] or [m]
    angleSlabs              = atand(dz./dx);                    %[degrees] - tan(alpha) = dz/dx
    depthShallowSlabDip    	= zColdSlabs + (zColdSlabs2-zColdSlabs)/2; %[nd] or [m]
    
    PLATE.slabDipPointsX  	= [xColdSlabs, xColdSlabs2];
    PLATE.slabDipPointsZ   	= [zColdSlabs, zColdSlabs2];
end
if logical(0) %check measurement points
    hold on
    plot(xColdSlabs/1e3,zColdSlabs/1e3,'ok')
    plot(xColdSlabs2/1e3,zColdSlabs2/1e3,'ob')
end

%% THEORETIC TRENCH-RETREAT VELOCITY
% v_TR = v_Stokes/tan(alpha)  after Capitanio et al., 2007
PLATE.theoreticTrenchVelocity = PLATE.SlabSinkingVelocity./tand(angleSlabs);

%% SLAB DIAGNOSTICS
SLAB.angle                  = angleSlabs;
SLAB.angleDepth             = depthShallowSlabDip;
SLAB.EtaToTake              = SlabEtaToTake;
SLAB.zIdxCheck              = zIdxCheck;
SLAB.zIdxCheck2             = zIdxCheck2;
SLAB.zIdxCheckUM         	= zIdxCheckUM;
SLAB.coldSlabs              = coldSlabs;
[SLAB] = f_slabDiagnostics(SLAB,BENDING,SWITCH,GRID,T_3D,ETA_3D,RHO_3D,STR_3D,BASALT_3D);
PLATE.slabArea2check     	= SLAB.Area2check;

%% TEST FIGURE
% figure(3),clf
% plot(Tprofile2check)
% title(['taken at a depth of ',num2str(depth2check/1e3,3),' km'])


%% STAGNANT-LID CHECK II: by maximum velocity
vPlateConvMin           	= v_rmsApproximate *4/5;         %[cm/a]; critical maximum plate convergence velocity ~0.1 cm/a
checkStagnantLid = true;
if checkStagnantLid %by maximum velocity-difference
    dummy = PLATE.PlateVelocityDiff;
    if dummy<vPlateConvMin
        if SWITCH.Verbose; disp(['   ',STYLE.SCHAR.indicationRightArrow,' stagnant-lid detected - low plate velocity difference (',num2str(dummy),'): no boundary indication.']); end %for first field plot
        PLATE.StagnantLid   = true;
        PLATE.Subduction    = 'noTracking'; %flag for no indication
        
        return %exit function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
    end
end

%% SEARCH FOR horizontal velocity discontinuities (i.e., PLATE BOUNDARIES), left->right
numXpoints2check    = floor(size(Z_3Dx,3)*0.05); %normalised to nr. z-grid points, set to 1 or bigger
% checkMeanVal       = 0.6;      %critical mean value of eg. [1 0 1 1]
dvdx                = 0.2;      %0.2, fraction of the max. dv that is looked for, over a width of dx_crit
% dvdxSurf starts with first value and if no subduction zone found, takes the next
% one - the second entry is used for spreading ridge tracking:
dvdxSurf            = [0.5 0.3 0.1 0.05]; %0.3, fraction of the max. dv that is looked for, over a width of dx_crit - 

dxCrit              = 0.015;    %0.015, non-dimensional width over which velocity peak2peak is checked
dxCrit2             = 0.05;     %0.05 (152km), non-dimensional width over which adjacent plate velocities & thicknesses are checked

dxCritTrench        = 0.17;     %0.17 (490km), non-dimensional width for the trench tracking
dxCritUP         	= 0.35;     %0.17 (980km), non-dimensional wider width over which adjacent UP-thicknesses are checked

dxCritSlab          = 0.17;     %0.17 (490km), non-dimensional critical width between slab in mantle and plate core boundaries

if strcmp(GRID.Type,'spherical2D')
    dx            	= mean(GRID.dxr(:));        %[rad]
else
    dx            	= mean(GRID.dx(:));         %[nd] or [m]
end

% adjust critical values
dxCrit              = dxCrit*SETUP.D;           %[nd] or [m]
dxCrit2             = dxCrit2*SETUP.D;          %[nd] or [m]
dxCritTrench     	= dxCritTrench*SETUP.D;   	%[nd] or [m]
dxCritUP          	= dxCritUP*SETUP.D;       	%[nd] or [m]
dxCritSlab          = dxCritSlab*SETUP.D;       %[nd] or [m]
if strcmp(GRID.Type,'spherical2D')
    dxCrit          = dxCrit/Rsurface;          %conversion to [rad]
    dxCrit2         = dxCrit2/Rsurface;         %conversion to [rad]
    dxCritTrench  	= dxCritTrench/Rsurface;  	%conversion to [rad]
    dxCritUP     	= dxCritUP/Rsurface;     	%conversion to [rad]
    dxCritSlab      = dxCritSlab/Rsurface;      %conversion to [rad]
end

%check if curve fitting toolbox is installed
try
    smooth(1:5); %it is present
    smoothVelocity = logical(1);
catch me %it is not present
    warning(me.message)
    smoothVelocity = logical(0);
end

%smooth velocity
if smoothVelocity
    if strcmp(GRID.Dim,'2-D')
        PlateSurfVelx  	= smooth(PlateSurfVelx,'moving');
        PlateVelx    	= smooth(PlateVelx,'moving');
        
    elseif strcmp(GRID.Dim,'3-D')
        filterSize      = 5;
        %F = [.05 .1 .05; .1 .4 .1; .05 .1 .05]; %filter
        %F = [1/9 1/9 1/9; 1/9 1/9 1/9; 1/9 1/9 1/9]; %like 2-D filter - with one loop!
        F = fspecial('gaussian',filterSize); %filterSize is the size of matrix
        %F = fspecial('gaussian');
        
        %add ghost points to prevent side effects from smoothing
        for ib=1:nb
            if ib==1
                dummy1  = PlateSurfVelx;
                dummy2  = PlateVelx;
                dummy3  = PlateSurfVely;
                dummy4  = PlateVely;
            elseif ib==2
                dummy1  = PlateSurfVelx_yang;
                dummy2  = PlateVelx_yang;
                dummy3  = PlateSurfVely_yang;
                dummy4  = PlateVely_yang;
            end
            for ii=1:filterSize
                dummy1  = [dummy1(1,:); dummy1; dummy1(end,:)];
                dummy1  = [dummy1(:,1), dummy1, dummy1(:,end)];
                dummy2  = [dummy2(1,:); dummy2; dummy2(end,:)];
                dummy2  = [dummy2(:,1), dummy2, dummy2(:,end)];
                dummy3  = [dummy3(1,:); dummy3; dummy3(end,:)];
                dummy3  = [dummy3(:,1), dummy3, dummy3(:,end)];
                dummy4  = [dummy4(1,:); dummy4; dummy4(end,:)];
                dummy4  = [dummy4(:,1), dummy4, dummy4(:,end)];
            end
            %smoothing
            for jj=1:17
                dummy   = conv2(dummy1,F,'same');
                dummy1  = dummy;
                dummy   = conv2(dummy2,F,'same');
                dummy2  = dummy;
                dummy   = conv2(dummy3,F,'same');
                dummy3  = dummy;
                dummy   = conv2(dummy4,F,'same');
                dummy4  = dummy;
            end
            if ib==1
                %remove ghostpoints
                PlateSurfVelx       = dummy1(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
                PlateVelx           = dummy2(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
                PlateSurfVely       = dummy3(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
                PlateVely           = dummy4(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
            elseif ib==2
                PlateSurfVelx_yang	= dummy1(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
                PlateVelx_yang    	= dummy2(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
                PlateSurfVely_yang 	= dummy3(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
                PlateVely_yang     	= dummy4(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
            end
        end
    end
end

%% DERIVE PLATE BOUNDARIES
numXpoints_crit             = max(1,round(dxCrit/dx));
numXpoints_crit2            = max(1,round(dxCrit2/dx));
numXpoints_critTrench       = max(1,round(dxCritTrench/dx));
numXpoints_critUP           = max(1,round(dxCritUP/dx));
critDiffSurfVel             = maxDiffPlateSurfVel *dvdxSurf;
critDiffPlateVel            = PLATE.PlateVelocityDiff *dvdx;
if strcmp(GRID.Dim,'3-D')
    dy                      = GRID.Y_3D(1,2,1)-GRID.Y_3D(1,1,1); %[nd] or [m] or [rad]
    numYpoints_crit2        = max(1,round(dxCrit2/dy)); %dx_crit2 must be the same in y-direction
    numXpoints_critTrench	= max(1,round(dxCritTrench/dy));
    numXpoints_critUP       = max(1,round(dxCritUP/dy));
end

PlateBoundary = zeros(GRID.nx,GRID.ny);
PlateBoundaryCore = PlateBoundary; PlateBoundarySurf = PlateBoundary;
PlateThicknessL = PlateBoundary; PlateThicknessR = PlateBoundary;
PlateThicknessWideL = PlateBoundary; PlateThicknessWideR = PlateBoundary;
%PlateBoundaryPolarity = zeros(GRID.nx,GRID.ny);
PlateBoundaryPolarity = cell(GRID.nx,GRID.ny); PlateBoundaryPolarity(:,:) = {'na'};
iaccurate = 0; iaccurate2 = 0;
if strcmp(GRID.Dim,'2-D')
    for iplateSurfCriterion=1:size(critDiffSurfVel,2)
        for ix=2:GRID.nx-1
            numXpointsLeft      = min(ix-1,numXpoints_crit);
            numXpointsRight     = min(GRID.nx-ix,numXpoints_crit);
            numXpointsLeft2     = min(ix-1,numXpoints_crit2);
            numXpointsRight2    = min(GRID.nx-ix,numXpoints_crit2);
            numXpointsLeft3     = min(ix-1,numXpoints_critUP);
            numXpointsRight3    = min(GRID.nx-ix,numXpoints_critUP);
            
            %THIS IS RELATIVE TO A MEAN VALUE:
            %     for iy=1:GRID.ny
            %         if mean( PlateSurfVelx(ix-numXpointsLeft:ix-1,iy)>=mean(PlateSurfVelx(:)) )>=checkMeanVal && ...
            %                 mean( PlateSurfVelx(ix+1:ix+numXpointsRight,iy)<mean(PlateSurfVelx(:)) )>=checkMeanVal;  %rightwards moving (towards subduction trench)
            %             PlateBoundary(ix,iy) = 1; %subduction trench
            %         elseif mean( PlateSurfVelx(ix-numXpointsLeft:ix-1,iy)<mean(PlateSurfVelx(:)) )>=checkMeanVal && ...
            %                 mean( PlateSurfVelx(ix+1:ix+numXpointsRight,iy)>=mean(PlateSurfVelx(:)) )>=checkMeanVal;  %leftwards moving (towards subduction trench)
            %             PlateBoundary(ix,iy) = 2; %spreading ridge
            %         else
            %             PlateBoundary(ix,iy) = 0; %no boundary
            %         end
            %     end
            
            %THIS IS LOOKING AT THE ACTUAL CHANGE BETWEEN THE GRID POINTS
            %checks the surface and plate core velocity for subduction trenches and spreading ridges
            for iy=1:GRID.ny
                try
                    dz_surf = peak2peak( PlateSurfVelx(ix-numXpointsLeft:ix+numXpointsRight2,iy) ); %peak2peak needs the Signal Processing Toolbox
                    dz_core = peak2peak( PlateVelx(ix-numXpointsLeft:ix+numXpointsRight2,iy) );
                catch
                    dz_surf = max(PlateSurfVelx(ix-numXpointsLeft:ix+numXpointsRight2,iy))-min(PlateSurfVelx(ix-numXpointsLeft:ix+numXpointsRight2,iy));
                    dz_core = max(PlateVelx(ix-numXpointsLeft:ix+numXpointsRight2,iy))-min(PlateVelx(ix-numXpointsLeft:ix+numXpointsRight2,iy));
                end
                
                %find plate boundaries at DEEPER LEVELS
                if dz_core>critDiffPlateVel && ... %plate boundary detected in core velocity - CONVERGING BOUNDARY
                        PlateVelx(ix-1,iy) > PlateVelx(ix+1,iy) && ...
                        median(PlateVelx(ix-numXpointsLeft2:ix-1,iy)) > median(PlateVelx(ix+1:ix+numXpointsRight,iy))  %converging (towards subduction trench)
                    if min(abs(xColdSlabs-GRID.X_3D(ix,iy,1)))<dxCritSlab %check for nearby slab
                        PlateBoundaryCore(ix,iy)        = 1; %subduction trench
                    end
                elseif dz_core>critDiffPlateVel && ... %plate boundary detected in core velocity - DIVERGING BOUNDARY
                        PlateVelx(ix-1,iy) < PlateVelx(ix+1,iy) && ...
                        median(PlateVelx(ix-numXpointsLeft2:ix-1,iy)) < median(PlateVelx(ix+1:ix+numXpointsRight,iy))  %diverging (towards subduction trench)
                    PlateBoundaryCore(ix,iy)            = 2; %spreading ridge
                end
                
                %find plate boundaries AT THE SURFACE
                if dz_surf>critDiffSurfVel(1,iplateSurfCriterion) && ... %plate boundary detected in surf velocity - CONVERGING BOUNDARY
                        PlateSurfVelx(ix-1,iy) > PlateSurfVelx(ix+1,iy) && ...
                        median(PlateSurfVelx(ix-numXpointsLeft2:ix-1,iy)) > median(PlateSurfVelx(ix+1:ix+numXpointsRight,iy))  %converging (towards subduction trench)
                    PlateBoundarySurf(ix,iy)            = 1; %subduction trench
                    %get left & right plate thicknesses
                    PlateThicknessL(ix,iy)              = lithoThicknessT_1600(ix-numXpointsLeft2,iy);
                    PlateThicknessR(ix,iy)              = lithoThicknessT_1600(ix+numXpointsRight2,iy);
                    PlateThicknessWideL(ix,iy)          = lithoThicknessT_1600(ix-numXpointsLeft3,iy);  %checks thickness further away from trench (needed for UP)
                    PlateThicknessWideR(ix,iy)          = lithoThicknessT_1600(ix+numXpointsRight3,iy);
                    %indicate subduction polarity
                    dummy1 = max(1,ix-dnx2deeperSlab); %min index of range
                    dummy2 = min(size(slabPolarity,1),ix+dnx2deeperSlab); %max index of range
                    if max(slabPolarity(dummy1:dummy2,1))==+1 && ... %R
                            min(slabPolarity(dummy1:dummy2,1))==-1 %L
                        %no indication as both polarities have been found within range
                    elseif max(slabPolarity(dummy1:dummy2,1))==+1 %R
                        PlateBoundaryPolarity{ix,iy}    = 'right';
                    elseif min(slabPolarity(dummy1:dummy2,1))==-1 %L
                        PlateBoundaryPolarity{ix,iy}    = 'left';
                    else
                        %no indication of subduction polarity
                    end

                elseif dz_surf>critDiffSurfVel(1,iplateSurfCriterion) && ... %plate boundary detected in surf velocity - DIVERGING BOUNDARY
                        PlateSurfVelx(ix-1,iy) < PlateSurfVelx(ix+1,iy) && ...
                        median(PlateSurfVelx(ix-numXpointsLeft2:ix-1,iy)) < median(PlateSurfVelx(ix+1:ix+numXpointsRight,iy))  %leftwards moving (towards subduction trench)
                    PlateBoundarySurf(ix,iy)            = 2; %spreading ridge
                end
            end
        end
        
        if iplateSurfCriterion==1
           PlateBoundarySurfOriginal     	= PlateBoundarySurf; 
        end
        if ~exist('PlateBoundarySurfAccurate','var')
            PlateBoundarySurfAccurate    	= ones(size(PlateBoundarySurf))*NaN; %is also used for spreading ridges, so in case of no subduction needs to be initialised
        end
        
        
        
        
        
        
        if multipleSlabsDetected
            dummy = bwmorph(PlateBoundarySurf,'shrink',Inf);  %shrink muliple points to one center point
            if sum(dummy)>=numberSlabsDetected %as long as at least the same number of subduction trenches is found as the number of slabs
                iaccurate = iaccurate+1;
            end
        else
            if max(PlateBoundarySurf(:)==1)==1 %as long as at least one subduction trench is found
                iaccurate = iaccurate+1;
            end
        end
        
        
        
        if isnan(max(PlateBoundarySurfAccurate)) %~exist('PlateBoundarySurfAccurate','var')
            PlateBoundarySurfAccurate     	= PlateBoundarySurf;
            iaccurate2 = iaccurate2+1;
        elseif iaccurate2==1
            PlateBoundarySurfInaccurate1  	= PlateBoundarySurf;
            iaccurate2 = iaccurate2+1;
        elseif iaccurate2==2
            PlateBoundarySurfInaccurate2	= PlateBoundarySurf;
            iaccurate2 = iaccurate2+1;
        elseif iaccurate2==3
            PlateBoundarySurfInaccurate3   	= PlateBoundarySurf;  %less restrictive
            iaccurate2 = iaccurate2+1;
        end
    end
    PlateBoundarySurf           = PlateBoundarySurfOriginal;
    
elseif strcmp(GRID.Dim,'3-D')
    ddx     = numXpoints_crit2; % #gridpionts over which difference is taken FOR SUBDUCTION TRENCH
    ddy     = numYpoints_crit2;
    ddx2    = ddx*2; % #gridpionts over which difference is taken FOR SPREADING RIDGE
    ddy2    = ddy*2;
    %     ddx         = 10; % #gridpionts over which difference is taken FOR SUBDUCTION TRENCH
    %     ddy         = 10; %MUST CURRENTLY BE MULTIPLE OF 2!!!!
    %     ddx2       	= 22; % #gridpionts over which difference is taken FOR SPREADING RIDGE
    %     ddy2       	= 22; %MUST CURRENTLY BE MULTIPLE OF 2!!!!
    
    %check if multiple of 2
    is_multipleOf2 = ( 2 * round(double(ddx)/2) == ddx );
    if ~is_multipleOf2; ddx=ddx+1; end %make multiple of two
    is_multipleOf2 = ( 2 * round(double(ddy)/2) == ddy );
    if ~is_multipleOf2; ddy=ddy+1; end %make multiple of two
    is_multipleOf2 = ( 2 * round(double(ddx2)/2) == ddx2 );
    if ~is_multipleOf2; ddx2=ddx2+1; end %make multiple of two
    is_multipleOf2 = ( 2 * round(double(ddy2)/2) == ddy2 );
    if ~is_multipleOf2; ddy2=ddy2+1; end %make multiple of two
    
    %PlateSurfVel = sqrt(PlateSurfVelx.^2 + PlateSurfVely.^2); %absolute surface velocity
    PlateVel = sqrt(PlateVelx.^2 + PlateVely.^2); %absolute velocity
    if strcmp(GRID.Type,'yinyang')
        PlateVel_yang = sqrt(PlateVelx_yang.^2 + PlateVely_yang.^2); %absolute velocity
    end
    
    %add ghostpoints
    nx_add = max(ddx,ddx2)/2; ny_add = max(ddy,ddy2)/2;
    
    if CONAREA.yPeriodic %x and y swapped in 3-D
        PlateVel    = [PlateVel((end-nx_add+1):end,:);   PlateVel;   PlateVel(1:nx_add,:)];
        PlateVelx   = [PlateVelx((end-nx_add+1):end,:);   PlateVelx;   PlateVelx(1:nx_add,:)];
        PlateVely   = [PlateVely((end-nx_add+1):end,:);   PlateVely;   PlateVely(1:nx_add,:)];
        if strcmp(GRID.Type,'yinyang')
            PlateVel_yang    = [PlateVel_yang((end-nx_add+1):end,:);   PlateVel_yang;   PlateVel_yang(1:nx_add,:)];
            PlateVelx_yang   = [PlateVelx_yang((end-nx_add+1):end,:);   PlateVelx_yang;   PlateVelx_yang(1:nx_add,:)];
            PlateVely_yang   = [PlateVely_yang((end-nx_add+1):end,:);   PlateVely_yang;   PlateVely_yang(1:nx_add,:)];
        end
    end
    if CONAREA.xPeriodic %x and y swapped in 3-D
        PlateVel    = [PlateVel(:,(end-ny_add+1):end)   PlateVel   PlateVel(:,1:ny_add)];
        PlateVelx   = [PlateVelx(:,(end-ny_add+1):end)   PlateVelx   PlateVelx(:,1:ny_add)];
        PlateVely   = [PlateVely(:,(end-ny_add+1):end)   PlateVely   PlateVely(:,1:ny_add)];
        if strcmp(GRID.Type,'yinyang')
            PlateVel_yang    = [PlateVel_yang(:,(end-ny_add+1):end)   PlateVel_yang   PlateVel_yang(:,1:ny_add)];
            PlateVelx_yang   = [PlateVelx_yang(:,(end-ny_add+1):end)   PlateVelx_yang   PlateVelx_yang(:,1:ny_add)];
            PlateVely_yang   = [PlateVely_yang(:,(end-ny_add+1):end)   PlateVely_yang   PlateVely_yang(:,1:ny_add)];
        end
    end
    
    %BOUNDARY POINTS WITH NON-PERIODIC SIDES ARE INCORRECT, REALLY.......
    dvdx = zeros(size(PlateVel)); dvdx2 = dvdx;
    dvdy = dvdx; dvdy2 = dvdy;
    if strcmp(GRID.Type,'yinyang')
        dvdx_yang = dvdx; dvdx2_yang = dvdx;
        dvdy_yang = dvdx; dvdy2_yang = dvdx;
    end
    for ix=1:size(PlateVel,1)
        for iy=1:size(PlateVel,2)
            %STANDARD DERIVATIVE over ddx
            if ix<=ddx/2 || ix>size(PlateVel,1)-ddx/2 %x-boundaries
                if CONAREA.yPeriodic %x and y swapped in 3-D
                    dvdx(ix,iy) = NaN; %ghostpoints are removed anyway
                    if strcmp(GRID.Type,'yinyang')
                        dvdx_yang(ix,iy) = NaN; %ghostpoints are removed anyway
                    end
                else
                    if ix<=ddx/2
                        %dvdx(ix,iy) = (PlateVelx(ix+ddx/2,iy)-PlateVelx(1,iy))/(((ix-1)+ddx/2)*dx); %boundary points
                        dvdx(ix,iy) = PlateVelx(ix+ddx/2,iy)-PlateVelx(1,iy); %boundary points
                        if strcmp(GRID.Type,'yinyang')
                            dvdx_yang(ix,iy) = PlateVelx_yang(ix+ddx/2,iy)-PlateVelx_yang(1,iy); %boundary points
                        end
                    elseif ix>size(PlateVel,1)-ddx/2
                        %dvdx(ix,iy) = (PlateVelx(end,iy)-PlateVelx(ix-ddx/2,iy))/(((GRID.nx-ix)+ddx/2)*dx); %boundary points
                        dvdx(ix,iy) = PlateVelx(end,iy)-PlateVelx(ix-ddx/2,iy); %boundary points
                        if strcmp(GRID.Type,'yinyang')
                            dvdx_yang(ix,iy) = PlateVelx_yang(end,iy)-PlateVelx_yang(ix-ddx/2,iy); %boundary points
                        end
                    else
                        error('if statement is wrong here!')
                    end
                end
            elseif iy<=ddy/2 || iy>size(PlateVel,2)-ddy/2 %y-boundaries
                if CONAREA.xPeriodic %x and y swapped in 3-D
                    dvdy(ix,iy) = NaN; %ghostpoints are removed anyway
                    if strcmp(GRID.Type,'yinyang')
                        dvdy_yang(ix,iy) = NaN; %ghostpoints are removed anyway
                    end
                else
                    if iy<=ddy/2
                        %dvdy(ix,iy) = (PlateVely(ix,iy+ddy/2)-PlateVely(ix,1))/(((iy-1)+ddy/2)*dy); %boundary points
                        dvdy(ix,iy) = PlateVely(ix,iy+ddy/2)-PlateVely(ix,1); %boundary points
                        if strcmp(GRID.Type,'yinyang')
                            dvdy_yang(ix,iy) = PlateVely_yang(ix,iy+ddy/2)-PlateVely_yang(ix,1); %boundary points
                        end
                    elseif iy>size(PlateVel,2)-ddy/2
                        %dvdy(ix,iy) = (PlateVely(ix,end)-PlateVely(ix,iy-ddy/2))/(((GRID.ny-iy)+ddy/2)*dy); %boundary points
                        dvdy(ix,iy) = PlateVely(ix,end)-PlateVely(ix,iy-ddy/2); %boundary points
                        if strcmp(GRID.Type,'yinyang')
                            dvdy_yang(ix,iy) = PlateVely_yang(ix,end)-PlateVely_yang(ix,iy-ddy/2); %boundary points
                        end
                    else
                        error('if statement is wrong here!')
                    end
                end
            else %inside
                % dvdx(ix,iy) = (PlateVelx(ix+ddx/2,iy)-PlateVelx(ix-ddx/2,iy))/(ddx*dx);
                % dvdy(ix,iy) = (PlateVely(ix,iy+ddy/2)-PlateVely(ix,iy-ddy/2))/(ddy*dy);
                dvdx(ix,iy) = PlateVelx(ix+ddx/2,iy)-PlateVelx(ix-ddx/2,iy);
                dvdy(ix,iy) = PlateVely(ix,iy+ddy/2)-PlateVely(ix,iy-ddy/2);
                if strcmp(GRID.Type,'yinyang')
                    dvdx_yang(ix,iy) = PlateVelx_yang(ix+ddx/2,iy)-PlateVelx_yang(ix-ddx/2,iy);
                    dvdy_yang(ix,iy) = PlateVely_yang(ix,iy+ddy/2)-PlateVely_yang(ix,iy-ddy/2);
                end
            end
            %DERIVATIVE OVER LARGER DISTANCE ddx2
            if ix<=ddx2/2 || ix>size(PlateVel,1)-ddx2/2 %x-boundaries
                if CONAREA.yPeriodic %x and y swapped in 3-D
                    dvdx2(ix,iy) = NaN; %ghostpoints
                    if strcmp(GRID.Type,'yinyang')
                        dvdx2_yang(ix,iy) = NaN; %ghostpoints
                    end
                else
                    if ix<=ddx2/2
                        %dvdx2(ix,iy) = (PlateVelx(ix+ddx2/2,iy)-PlateVelx(1,iy))/(((ix-1)+ddx2/2)*dx); %boundary points
                        dvdx2(ix,iy) = PlateVelx(ix+ddx2/2,iy)-PlateVelx(1,iy); %boundary points
                        if strcmp(GRID.Type,'yinyang')
                            dvdx2_yang(ix,iy) = PlateVelx_yang(ix+ddx2/2,iy)-PlateVelx_yang(1,iy); %boundary points
                        end
                    elseif ix>size(PlateVel,1)-ddx2/2
                        %dvdx2(ix,iy) = (PlateVelx(end,iy)-PlateVelx(ix-ddx2/2,iy))/(((GRID.nx-ix)+ddx2/2)*dx); %boundary points
                        dvdx2(ix,iy) = PlateVelx(end,iy)-PlateVelx(ix-ddx2/2,iy); %boundary points
                        if strcmp(GRID.Type,'yinyang')
                            dvdx2_yang(ix,iy) = PlateVelx_yang(end,iy)-PlateVelx_yang(ix-ddx2/2,iy); %boundary points
                        end
                    else
                        error('if statement is wrong here!')
                    end
                end
            elseif iy<=ddy2/2 || iy>size(PlateVel,2)-ddy2/2 %y-boundaries
                if CONAREA.yPeriodic %x and y swapped in 3-D
                    dvdy2(ix,iy) = NaN; %ghostpoints
                    if strcmp(GRID.Type,'yinyang')
                        dvdy2_yang(ix,iy) = NaN; %ghostpoints
                    end
                else
                    if iy<=ddy2/2
                        %dvdy2(ix,iy) = (PlateVely(ix,iy+ddy2/2)-PlateVely(ix,1))/(((iy-1)+ddy2/2)*dy); %boundary points
                        dvdy2(ix,iy) = PlateVely(ix,iy+ddy2/2)-PlateVely(ix,1); %boundary points
                        if strcmp(GRID.Type,'yinyang')
                            dvdy2_yang(ix,iy) = PlateVely_yang(ix,iy+ddy2/2)-PlateVely_yang(ix,1); %boundary points
                        end
                    elseif iy>size(PlateVel,2)-ddy2/2
                        %dvdy2(ix,iy) = (PlateVely(ix,end)-PlateVely(ix,iy-ddy2/2))/(((GRID.ny-iy)+ddy2/2)*dy); %boundary points
                        dvdy2(ix,iy) = PlateVely(ix,end)-PlateVely(ix,iy-ddy2/2); %boundary points
                        if strcmp(GRID.Type,'yinyang')
                            dvdy2_yang(ix,iy) = PlateVely_yang(ix,end)-PlateVely_yang(ix,iy-ddy2/2); %boundary points
                        end
                    else
                        error('if statement is wrong here!')
                    end
                end
            else %inside
                % dvdx2(ix,iy) = (PlateVelx(ix+ddx2/2,iy)-PlateVelx(ix-ddx2/2,iy))/(ddx2*dx);
                % dvdy2(ix,iy) = (PlateVely(ix,iy+ddy2/2)-PlateVely(ix,iy-ddy2/2))/(ddy2*dy);
                dvdx2(ix,iy) = PlateVelx(ix+ddx2/2,iy)-PlateVelx(ix-ddx2/2,iy);
                dvdy2(ix,iy) = PlateVely(ix,iy+ddy2/2)-PlateVely(ix,iy-ddy2/2);
                if strcmp(GRID.Type,'yinyang')
                    dvdx2_yang(ix,iy) = PlateVelx_yang(ix+ddx2/2,iy)-PlateVelx_yang(ix-ddx2/2,iy);
                    dvdy2_yang(ix,iy) = PlateVely_yang(ix,iy+ddy2/2)-PlateVely_yang(ix,iy-ddy2/2);
                end
            end
        end
    end
    
    Array   = dvdx+dvdy;
    Array2  = dvdx2+dvdy2;
    if strcmp(GRID.Type,'yinyang')
        Array_yang   = dvdx_yang+dvdy_yang;
        Array2_yang  = dvdx2_yang+dvdy2_yang;
    end
    
    %remove ghostpoints
    if CONAREA.yPeriodic %x and y swapped in 3-D
        Array     	= Array((nx_add+1):end-nx_add,:);
        Array2     	= Array2((nx_add+1):end-nx_add,:);
        if strcmp(GRID.Type,'yinyang')
            Array_yang     	= Array_yang((nx_add+1):end-nx_add,:);
            Array2_yang     = Array2_yang((nx_add+1):end-nx_add,:);
        end
    end
    if CONAREA.xPeriodic %x and y swapped in 3-D
        Array     	= Array(:,(ny_add+1):end-ny_add);
        Array2     	= Array2(:,(ny_add+1):end-ny_add);
        if strcmp(GRID.Type,'yinyang')
            Array_yang     	= Array_yang(:,(ny_add+1):end-ny_add);
            Array2_yang     = Array2_yang(:,(ny_add+1):end-ny_add);
        end
    end
    
    if logical(0)
        figure(22),clf
        subplot(1,2,1)
        contourf(Array),colorbar,title(num2str(ddx))
        xlabel('y'); ylabel('x');
        subplot(1,2,2)
        contourf(Array2),colorbar,title(num2str(ddx2))
        xlabel('y'); ylabel('x');
    end
    
    %FIND CONNECTED PLATE BOUNDARIES (FROM DIVERGENCE, i.e., trenches and ridges)
    %SUBDUCTION TRENCHES
    %     minArray = min(Array(:));
    %     maxArray = max(Array(:));
    %     diffArray = maxArray-minArray;
    %     CONAREA.threshold       = minArray+ 0.1*diffArray;
    CONAREA.threshold       = -critDiffPlateVel;
    %     CONAREA.threshold       = -5e-6; %-4.0e-6;
    CONAREA.thresLogical	= 'smaller';
    CONAREA.nxSubgrid       = 32;
    [CONAREA] = f_Connectivity(Array,CONAREA,GRID);
    %conversion to new variables
    SubductionSurf          = CONAREA.LineP;        %[plotting dimension]
    SubSurfOutline          = CONAREA.OutlineP;     %[plotting dimension]
    SubCenterP              = CONAREA.CenterP;      %[plotting dimension]
    SubNormalP              = CONAREA.NormalP;      %unit vector
    SubStrikeP              = CONAREA.StrikeP;      %unit vector
    SubLength               = CONAREA.LengthAreas;  %[plotting dimension], or NaN
    SubLengthTotal          = CONAREA.LengthTotal;  %[plotting dimension], or NaN
    if strcmp(GRID.Type,'yinyang')
        [CONAREA_yang] = f_Connectivity(Array_yang,CONAREA,GRID);
        SubductionSurf_yang     = CONAREA_yang.LineP;       %[plotting dimension]
        SubSurfOutline_yang     = CONAREA_yang.OutlineP;    %[plotting dimension]
        SubCenterP_yang         = CONAREA_yang.CenterP;     %[plotting dimension]
        SubNormalP_yang         = CONAREA_yang.NormalP;     %unit vector
        SubStrikeP_yang         = CONAREA_yang.StrikeP;     %unit vector
        SubLength_yang          = CONAREA_yang.LengthAreas; %[plotting dimension], or NaN
        SubLengthTotal          = SubLengthTotal+CONAREA_yang.LengthTotal; %[plotting dimension], or NaN
    end
    
    %SPREADING RIDGES
    %     minArray = min(Array2(:));
    %     maxArray = max(Array2(:));
    %     diffArray = maxArray-minArray;
    %     CONAREA.threshold       = maxArray- 0.2*diffArray;
    CONAREA.threshold       = critDiffPlateVel;
    %     CONAREA.threshold       = 2e-6;
    CONAREA.thresLogical	= 'bigger';
    CONAREA.nxSubgrid       = 32;
    [CONAREA] = f_Connectivity(Array2,CONAREA,GRID);
    %conversion to new variables
    SpreadingSurf           = CONAREA.LineP;        %[plotting dimension]
    SprSurfOutline          = CONAREA.OutlineP;     %[plotting dimension]
    SprCenterP              = CONAREA.CenterP;      %[plotting dimension]
    SprNormalP              = CONAREA.NormalP;      %unit vector
    SprStrikeP              = CONAREA.StrikeP;      %unit vector
    SprLength               = CONAREA.LengthAreas;  %[plotting dimension], or NaN
    SprLengthTotal          = CONAREA.LengthTotal;  %[plotting dimension], or NaN
    if strcmp(GRID.Type,'yinyang')
        [CONAREA_yang] = f_Connectivity(Array2_yang,CONAREA,GRID);
        SpreadingSurf_yang      = CONAREA_yang.LineP;       %[plotting dimension]
        SprSurfOutline_yang     = CONAREA_yang.OutlineP;    %[plotting dimension]
        SprCenterP_yang         = CONAREA_yang.CenterP;     %[plotting dimension]
        SprNormalP_yang         = CONAREA_yang.NormalP;     %unit vector
        SprStrikeP_yang         = CONAREA_yang.StrikeP;     %unit vector
        SprLength_yang          = CONAREA_yang.LengthAreas; %[plotting dimension], or NaN
        SprLengthTotal          = SprLengthTotal+CONAREA_yang.LengthTotal; %[plotting dimension], or NaN
    end
end


%% SET KIND OF PLATE BOUNDARIES
if strcmp(GRID.Dim,'2-D')
    % subduction trenches at the surface and
    % spreading ridges at deeper levels
    xySubdCore     = [GRID.X_3D(PlateBoundaryCore==1),GRID.Y_3D(PlateBoundaryCore==1)]; %subduction core
    xySpreadCore   = [GRID.X_3D(PlateBoundaryCore==2),GRID.Y_3D(PlateBoundaryCore==2)]; %spreading core
    PlateBoundarySurf4Spreading = PlateBoundarySurfInaccurate1;
    xySpreadSurf   = [GRID.X_3D(PlateBoundarySurf4Spreading==2),GRID.Y_3D(PlateBoundarySurf4Spreading==2)]; %spreading surface
    %% SUBDUCTION TRENCH
    for iaccuracy=1:iaccurate+1
        if iaccuracy==1
            PlateBoundarySurfDummy      = PlateBoundarySurf; %original (accurate)
        elseif iaccuracy==2
            PlateBoundarySurfDummy      = PlateBoundarySurfAccurate;
        elseif iaccuracy==3
            PlateBoundarySurfDummy      = PlateBoundarySurfInaccurate1;
        elseif iaccuracy==4
            PlateBoundarySurfDummy      = PlateBoundarySurfInaccurate2;
        elseif iaccuracy==5
            PlateBoundarySurfDummy      = PlateBoundarySurfInaccurate3; %(less restrictive)
        end
        xySubdSurf     = [GRID.X_3D(PlateBoundarySurfDummy==1),GRID.Y_3D(PlateBoundarySurfDummy==1)]; %subduction surface
        
        % CHECK FOR DEEP SUBDUCTION FAULT
        % Critical Values
        dxSurfcoreCrit      = 0.15;     %non-dimensional value - Horizontal distance between deep subduction and shallow subduction ~456 km <<<<<<<<<<<<<<<<<<<<<<
        dxSurfcoreCrit      = dxSurfcoreCrit*SETUP.D;  %[nd] or [m]
        if strcmp(GRID.Type,'spherical2D')
            dxSurfcoreCrit  = dxSurfcoreCrit/Rsurface;     %conversion to [rad]
        end
        
        dummy = zeros(size(xySubdSurf,1),1);
        for i=1:size(xySubdSurf,1) %SUBDUCTION
            dx_sc           = xySubdCore(:,1)-xySubdSurf(i,1);
            if isempty(dx_sc(abs(dx_sc)<dxSurfcoreCrit)) %no core subduction close by => remove surface value
                PlateBoundarySurfDummy( GRID.X_3D(:,1,1)==xySubdSurf(i,1) ) = 0;
                dummy(i,1)  = xySubdSurf(i,1);
            end
            dx_sc = 0;
        end

        
        
        
        breakHere = false;
        if multipleSlabsDetected
            %exit if at least the same #subduction zones has been detected as #slabs
            dummy = bwmorph(PlateBoundarySurf,'shrink',Inf);  %shrink muliple points to one center point
            if sum(dummy)>=numberSlabsDetected %as long as at least the same number of subduction trenches is found as the number of slabs
                breakHere       = true;
            end
        else
            %exit if at least one subduction zone has been detected
            if max(PlateBoundarySurfDummy==1)==1
                breakHere       = true;
            end
        end
        if breakHere
            dummy(dummy==0) = [];
            if ~isempty(dummy) && SWITCH.Verbose
                disp(['   ',STYLE.SCHAR.indicationRightArrow,' no deeper fault close by: ',num2str(size(dummy,1)),'/',num2str(i),' subduction trackers removed.']);
            end
            
            break
        end
        clearvars breakHere
        
        
        
        
    end
    %% SPREADING RIDGE
    dxSurfcoreCrit      = 0.02;     %non-dimensional value - Horizontal distance between deep spreading and shallow spreading ~135 km <<<<<<<<<<<<<<<<<<<<<<
    dxSurfcoreCrit      = dxSurfcoreCrit*SETUP.D;  %[nd] or [m]
    if strcmp(GRID.Type,'spherical2D')
        dxSurfcoreCrit  = dxSurfcoreCrit/Rsurface;     %conversion to [rad]
    end
    
    dummy = zeros(size(xySpreadSurf,1),1);
    for i=1:size(xySpreadSurf,1) %SPREADING
        dx_sc               = xySpreadCore(:,1)-xySpreadSurf(i,1);
        if isempty(dx_sc(abs(dx_sc)<dxSurfcoreCrit)) %no core spreading close by => remove surface value
            PlateBoundarySurf4Spreading( GRID.X_3D(:,1,1)==xySpreadSurf(i,1) ) = 0;
            dummy(i,1)      = xySpreadSurf(i,1);
        else
            %nothing to do
        end
        dx_sc = 0;
    end
    dummy(dummy==0) = [];
    if ~isempty(dummy) && SWITCH.Verbose %for first field plot
        disp(['   ',STYLE.SCHAR.indicationRightArrow,' no deeper fault close by: some spreading trackers removed.']);
    end
    
    % ADD CHOSEN VALUE
    %subduction from plate surface
    PlateBoundary(PlateBoundarySurfDummy==1)        = 1;
    %spreading from plate core
    % PlateBoundary(PlateBoundaryCore==2) = 2;
    %spreading from plate surface
    PlateBoundary(PlateBoundarySurf4Spreading==2)  	= 2;
    
    % TEST VALUES
    if logical(0)
        figure(2),clf
        plot(PlateSurfVelx)
        hold on
        plot(smooth(PlateSurfVelx))
        % hold on
        % % plot(smooth(PlateSurfVelx,'moving'))
        % % hold on
        % % plot(smooth(PlateSurfVelx,'rlowess'))
        % % hold on
        % % plot(diff(PlateSurfVelx))
        % % hold on
        plot(PlateVelx)
        hold on
        plot(smooth(PlateVelx))
        % hold on
        % % plot(diff(PlateVelx))
        % % hold on
        plot(PlateBoundary)
        figure(1)
    end
    
elseif strcmp(GRID.Dim,'3-D')
    %find subduction polarity
    for ib=1:nb
        if ib==1
            dummySNP = SubNormalP;
            dummySCP = SubCenterP;
            dummyUPP = zeros(size(SubCenterP)).*NaN;
            dummyLPP = dummyUPP;
        else
            dummySNP = SubNormalP_yang;
            dummySCP = SubCenterP_yang;
            dummyUPP = zeros(size(SubCenterP_yang)).*NaN;
            dummyLPP = dummyUPP;
        end
        xl1 = zeros(1,2); yl1 = xl1;
        for iarea=1:size(dummySNP,1)
            for isubarea=1:size(dummySNP,2)
                if isnan(dummySCP(iarea,isubarea,1))
                    if SWITCH.Verbose; warning('NaN detected!'); end
                    continue
                end
                for iside=1:2
                    if iside==1
                        fac = 1; slabAmount = 0;
                    else
                        fac = -1;
                    end
                    %make normal line to check where slab is
                    nLinepoints = 100;
                    lineWidth1 = dxCritSlab*GRID.dimFactor;%MAYBE ADJUST dx_critSlab........
                    xl1(1,1) = dummySCP(iarea,isubarea,1);
                    xl1(1,2) = dummySCP(iarea,isubarea,1)+ fac.*dummySNP(iarea,isubarea,1).*lineWidth1;
                    xl1 = linspace(xl1(1,1),xl1(1,2),nLinepoints);
                    yl1(1,1) = dummySCP(iarea,isubarea,2);
                    yl1(1,2) = dummySCP(iarea,isubarea,2)+ fac.*dummySNP(iarea,isubarea,2).*lineWidth1;
                    yl1 = linspace(yl1(1,1),yl1(1,2),nLinepoints);
                    
                    x2d_dummy = GRID.X_3Dp(:,:,1);
                    y2d_dummy = GRID.Y_3Dp(:,:,1);
                    indexMatrix1 = 1:size(x2d_dummy,1)*size(x2d_dummy,2);
                    indexMatrix1 = reshape(indexMatrix1,size(x2d_dummy));
                    %THIS IS SLOW!********************************************************
                    lin_ind = griddata(x2d_dummy,y2d_dummy,indexMatrix1,xl1,yl1,'nearest'); % where xl1 and yl1 are the line coordinates (should be spaced not wider than slab)
                    %*********************************************************************
                    [sub_ind(1,:),sub_ind(2,:)] = ind2sub(size(x2d_dummy),lin_ind);
                    %remove duplicates
                    il=0;
                    while il<size(sub_ind,2)-1
                        il = il+1;
                        if sub_ind(1,il)==sub_ind(1,il+1) && sub_ind(1,il)==sub_ind(1,il+1)
                            sub_ind(:,il+1) = []; %erase duplicate
                        end
                    end
                    if logical(0)
                        figure(11) ; grid on ; hold on
                        plot(GRID.X_3Dp(coldestSlab),GRID.Y_3Dp(coldestSlab),'k.') %slab
                        plot(dummySCP(:,:,1),dummySCP(:,:,2),'g') %trench
                        plot(xl1,yl1); %line to check
                        plot(GRID.X_3Dp(sub_ind(1,:),1,1),GRID.Y_3Dp(1,sub_ind(2,:),1),'or')%resulting grid points to check
                        % pcolor(GRID.X_3D(:,:,1),GRID.Y_3D(:,:,1),Z) ;  % plot the domain with the illuminated blocks
                    end
                    if ib==1
                        dummy = coldestSlab(sub_ind(1,:),sub_ind(2,:));
                    else
                        dummy = coldestSlab_yang(sub_ind(1,:),sub_ind(2,:));
                    end
                    if max(dummy(:))==1 && sum(dummy(:))>slabAmount %line crosses slab location surface projection
                        %upper plate
                        dummyUPP(iarea,isubarea,:) = fac.*-dummySNP(iarea,isubarea,:);
                        dummyLPP(iarea,isubarea,:) = fac.*dummySNP(iarea,isubarea,:); %opposite
                        slabAmount = sum(dummy(:));
                    elseif iside==2 && slabAmount==0 %no slab found on both sides
                        dummyUPP(iarea,isubarea,:) = [NaN; NaN];
                        dummyLPP(iarea,isubarea,:) = [NaN; NaN]; %opposite
                    else %no slab found nearby on current side
                    end
                    clearvars sub_ind dummy lin_ind
                end
            end
        end
        if ib==1
            upperPlateP = dummyUPP;
            lowerPlateP = dummyLPP;
        else
            upperPlateP_yang = dummyUPP;
            lowerPlateP_yang = dummyLPP;
        end
    end
    clearvars dummySNP dummyUPP dummyLPP dummySCP
end


%% SOME FUNCTION OUTPUT
if strcmp(GRID.Dim,'2-D')
    %output related to specific plate boundary should be added later on!
    PLATE.PlateBoundary         = PlateBoundary;        %0: nothing, 1: subduction, 2: spreading    
    PLATE.Subduction            = [ GRID.X_3Dp(PLATE.PlateBoundary==1), GRID.Y_3Dp(PLATE.PlateBoundary==1) ];
    PLATE.SubPolarity           = PlateBoundaryPolarity(PLATE.PlateBoundary==1);
    PLATE.Spreading             = [ GRID.X_3Dp(PLATE.PlateBoundary==2), GRID.Y_3Dp(PLATE.PlateBoundary==2) ];
    
elseif strcmp(GRID.Dim,'3-D')
    PLATE.Subduction            = SubductionSurf;
    if isempty(SubductionSurf); PLATE.Subduction = 'noTracking'; end %flag for no indication
    PLATE.SubductionOut         = SubSurfOutline;
    PLATE.SubCenterP            = SubCenterP;
    PLATE.Spreading             = SpreadingSurf;
    if isempty(SpreadingSurf); PLATE.Spreading = 'noTracking'; end %flag for no indication
    PLATE.SpreadingOut          = SprSurfOutline;
    PLATE.UpperPlateP           = upperPlateP;
    PLATE.LowerPlateP           = lowerPlateP;
    PLATE.SubNormalP            = SubNormalP;
    PLATE.SubStrikeP            = SubStrikeP;
    PLATE.SprCenterP            = SprCenterP;
    PLATE.SprNormalP            = SprNormalP;
    PLATE.SprStrikeP            = SprStrikeP;
    %add here more function output.....
    if strcmp(GRID.Type,'yinyang')
        PLATE.Subduction_yang       = SubductionSurf_yang;
        if isempty(SubductionSurf_yang); PLATE.Subduction_yang = 'noTracking'; end %flag for no indication
        PLATE.SubductionOut_yang    = SubSurfOutline_yang;
        PLATE.SubCenterP_yang       = SubCenterP_yang;
        PLATE.Spreading_yang        = SpreadingSurf_yang;
        if isempty(SpreadingSurf_yang); PLATE.Spreading_yang = 'noTracking'; end %flag for no indication
        PLATE.SpreadingOut_yang     = SprSurfOutline_yang;
        PLATE.UpperPlateP_yang      = upperPlateP_yang;
        PLATE.LowerPlateP_yang  	= lowerPlateP_yang;
        PLATE.SubNormalP_yang       = SubNormalP_yang;
        PLATE.SubStrikeP_yang       = SubStrikeP_yang;
        PLATE.SprCenterP_yang       = SprCenterP_yang;
        PLATE.SprNormalP_yang       = SprNormalP_yang;
        PLATE.SprStrikeP_yang       = SprStrikeP_yang;
        %add here more function output.....
    end
end

%% CHECK FOR MULTIPLE TRACKINGS
dx_criticalSubSub       = 0.15*max(GRID.Z_3Dp(:));      %~400 km; adjust critical dx here (distance between adjacent subduction zones)
dx_critical             = 0.05*max(GRID.Z_3Dp(:));      %~140 km; adjust critical dx here (distance between adjacent pl.boundaries)
if strcmp(GRID.Type,'spherical2D')
    if ~SWITCH.DimensionalMode && ~SWITCH.DimensionalInput %non-dimensional mode
        dx_criticalSubSub   = dx_criticalSubSub/Rsurface; 	%conversion to radians
        dx_critical         = dx_critical/Rsurface;         %conversion to radians
    else
        dx_criticalSubSub   = dx_criticalSubSub/(Rsurface*GRID.m2p); 	%conversion to radians
        dx_critical         = dx_critical/(Rsurface*GRID.m2p);         %conversion to radians
    end
end
if strcmp(GRID.Dim,'2-D')
    % and combine multiple close ones to a single one
    if size(PLATE.Subduction,1)>1 || size(PLATE.Spreading,1)>1 %if multiple subduction/spreading zones
        removeSubduction = zeros(size(PLATE.Subduction,1),1);
        removeSpreading = zeros(size(PLATE.Spreading,1),1);
        removeSpreading2 = removeSpreading;
        removeSubX = []; removeSprX = []; removeSpr2X = [];
        for i=1:size(PLATE.Subduction,1) %from left to right
            if i~=1 && PLATE.Subduction(i-1,1)>PLATE.Subduction(i,1)-dx_criticalSubSub	%SUB-SUB
                keepFurthest = logical(1); %1: keeps the trench indicators furthest away from the subduction zone
                if keepFurthest
                    if strcmp(PLATE.SubPolarity(i,1),'left') %left facing subduction
                        removeSubduction(i,1) = i-1; %remove indicators on the right
                    else
                        removeSubduction(i,1) = i; %remove indicators on the left
                    end
                else
                    if strcmp(PLATE.SubPolarity(i,1),'left') %right facing subduction
                        removeSubduction(i,1) = i; %remove indicators on the right
                    else
                        removeSubduction(i,1) = i-1; %remove indicators on the left
                    end
                end
                removeSubX = [removeSubX PLATE.Subduction(i,1)];
            end
            % b(b(b>a(i)-dx_critical)<a(i)+dx_critical)
            dummy = PLATE.Spreading(PLATE.Spreading(:,1)>PLATE.Subduction(i,1)-dx_critical); %bigger than crit. value
            dummy = dummy(dummy<PLATE.Subduction(i,1)+dx_critical); %bigger and now also smaller than crit. value
            if ~isempty(dummy) %SUB-RIDGE
                dummy2 = ismember(PLATE.Spreading(:,1),dummy(:,1));
                removeSpreading2 = min(removeSpreading2+dummy2,1);
            end
        end
        
        i_local=0; idx_middle=0;
        for i=1:size(PLATE.Spreading,1) %RIDGE-RIDGE
            if i~=size(PLATE.Spreading,1) && ... %boundary condition
                    PLATE.Spreading(i+1,1)<PLATE.Spreading(i,1)+dx_critical	%there is another one to the right
                i_local = i_local+1; %no. of local trackers that need to be combined
                idx_middle = floor(i_local/2); %min is 0, max is nx/2
                %remove it
                removeSpreading(i,1) = i;
                removeSprX = [removeSprX PLATE.Spreading(i,1)];
            else
                %remove it
                removeSpreading(i,1) = i;
                removeSprX = [removeSprX PLATE.Spreading(i,1)];
                %keep middle one
                removeSpreading(i-idx_middle,1) = 0;
                removeSprX(:,end-idx_middle)    = [];
                
                i_local=0; idx_middle=0;
            end
        end
        %.
        if ~isempty(removeSubduction(removeSubduction==0))
            removeSubduction(removeSubduction==0)	= [];
        end
        if ~isempty(removeSubX) && SWITCH.Verbose %for first field plot
            disp(['   ',STYLE.SCHAR.indicationRightArrow,' dublicated subduction tracker removed at x = ',num2str(removeSubX,3)])
        end
        %.
        if ~isempty(removeSpreading(removeSpreading==0))
            removeSpreading(removeSpreading==0)   	= [];
        end
        if ~isempty(removeSprX) && SWITCH.Verbose %for first field plot
            disp(['   ',STYLE.SCHAR.indicationRightArrow,' dublicated spreading center tracker removed at x = ',num2str(removeSprX,3)]);
        end
        %.
        removeSprX2 = PLATE.Spreading(removeSpreading2==1,1);
        if ~isempty(removeSprX2) && SWITCH.Verbose %for first field plot
            disp(['   ',STYLE.SCHAR.indicationRightArrow,' to-close-to-subduction spreading center tracker removed at x = ',num2str(removeSprX2',3)]);
        end
        
        %actually remove (i.e., mark) the tracking points
        if ~size(removeSubduction,2)==0 %only if it exists
            PLATE.Subduction(removeSubduction(:,1),:)   = NaN;
        end
        if ~size(removeSpreading,2)==0 %only if it exists
            PLATE.Spreading(removeSpreading(:,1),:)     = NaN;
        end
        if ~size(removeSpreading2,2)==0 %only if it exists
            PLATE.Spreading(removeSpreading2==1,:)      = NaN;
        end
    end
    %remove NaN values
    if size(PLATE.Subduction,1)>1 %needs at least one value
        PLATE.Subduction(isnan(PLATE.Subduction(:,1)),:)    = [];
        dummy = PLATE.SubPolarity(~isnan(PLATE.Subduction(:,1)),:); PLATE.SubPolarity = dummy; %XXXXXX new check if ok
        PLATE.Spreading(isnan(PLATE.Subduction(:,1)),:)     = [];
    end
    %make sure they are the same size (needed for saving data)
    if size(PLATE.Subduction,1)~=size(PLATE.Spreading,1)
        dummy = max(size(PLATE.Subduction,1),size(PLATE.Spreading,1));
        PLATE.Subduction    = [PLATE.Subduction; zeros(dummy-size(PLATE.Subduction,1),2)*NaN];
        PLATE.SubPolarity   = [PLATE.SubPolarity; cell(dummy-size(PLATE.Subduction,1),1)]; %XXXXXX new check if ok
        PLATE.Spreading     = [PLATE.Spreading; zeros(dummy-size(PLATE.Spreading,1),2)*NaN];
    end
    % %% test output
    %     figure(2),clf
    %     plot(PlateSurfVelx)
    %     hold on
    %     plot(diff(PlateSurfVelx))
    %     hold on
    %     plot(PlateVelx)
    %     hold on
    %     plot(diff(PlateVelx))
    %     hold on
    %     plot(PlateBoundary)
    
elseif strcmp(GRID.Dim,'3-D')
    %check for multiple trackings not needed for 3-D.
end

%% HANDLE PROBLEMATIC FIRST TIME STEP
% highly influenced by isostatic adjustment
if strcmp(GRID.Dim,'2-D')
    if FILE.number==0
        if ~isempty(strcmp(PLATE.SubPolarity,'na')) || ... %size(find(PlateBoundary==1),1)>1
                isempty(PLATE.Subduction) || max(max(~isnan(PLATE.Subduction)))==0 %empty or all NaN
            % DELETE DATA: no indication in plot
            %PLATE.Subduction = PLATE.Subduction *NaN;
            PLATE.Subduction    = 'noTracking'; %flag for no indication
            PLATE.Spreading     = PLATE.Spreading *NaN;
            warning('Plate boundaries NOT indicated at the first time step.')
            return %exit function
        else
            % WARNING
            warning('Tracking plate boundaries during the first time step might be inaccurate...')
        end
    end
elseif strcmp(GRID.Dim,'3-D')
    if FILE.number==0
        PLATE.Subduction    = 'noTracking'; %flag for no indication
        if strcmp(GRID.Type,'yinyang')
            PLATE.Subduction_yang    = 'noTracking'; %flag for no indication
        end
        warning('Tracking plate boundaries during the first time step might be inaccurate...')
    end
end

%% REMOVE DELETED ENTRIES
if strcmp(GRID.Dim,'2-D')
    PLATE.Subduction(isnan(PLATE.Subduction(:,1)),:)    = []; %in GRID.Xdim
    dummy = PLATE.SubPolarity(~isnan(PLATE.Subduction(:,1)),:); PLATE.SubPolarity = dummy;
    PLATE.Spreading(isnan(PLATE.Spreading(:,1)),:)      = []; %in GRID.Xdim
    %prevent the arrays of becoming empty
    if isempty(PLATE.Subduction);   PLATE.Subduction    = [NaN, NaN];   end
    if isempty(PLATE.SubPolarity);  PLATE.SubPolarity   = {'', ''};     end
    if isempty(PLATE.Spreading);    PLATE.Spreading     = [NaN, NaN];   end
end

%% USE 2nd SLAB-POLARITY DIAGNOSTIC (if 1st failed)
if strcmp(GRID.Dim,'2-D') && size(PLATE.SubPolarity,1)==1 && (strcmp(PLATE.SubPolarity(1,1),'na') || strcmp(PLATE.SubPolarity(1,1),'')) %only if 1 subduction zone present
    if strcmp(BENDING.slabPolarity,'L')
        PLATE.SubPolarity{1,1}      = 'left';
    elseif strcmp(BENDING.slabPolarity,'R')
        PLATE.SubPolarity{1,1}      = 'right';
    end
end

%% ERROR CHECKS
couldNotBeTracked = false; couldNotBeTracked2 = false;
if strcmp(GRID.Dim,'2-D')
%     if size(PLATE.Subduction,1)>size(angleSlabs,1) %XXXXXXXXXXXX THIS NEEDS TO BE IMPROVED
%         PLATE.Subduction(isnan(PLATE.Subduction(:,1)),:) = [];
%         if size(PLATE.Subduction,1)>size(angleSlabs,1)
%             warning('#shallow-slab angles differs from #subduction trenches: some subduction trenches tracking points have been removed.')
%             PLATE.Subduction(size(angleSlabs,1)+1:end,:) = [];  %THIS JUST KEEPS THE LEFT-MOST TRACKINGS (might not be correct for multiple slabs)
%         end
%     end %XXXXXXXXXXXX THIS NEEDS TO BE IMPROVED
    
    if size(angleSlabs,1)==1 && size(angleSlabs,2)==1 && isnan(angleSlabs)
        dummy                	= size(angleSlabs,1);
        angleSlabs(dummy+1:size(PLATE.Subduction,1),:)   	= NaN; %add values
        depthShallowSlabDip(dummy+1:size(PLATE.Subduction,1),:)	= NaN;
    elseif size(PLATE.Subduction,1)<size(angleSlabs,1)
        warning('#shallow-slab angles > #subduction trenches: Shallow-shallow-slab angles could not be found.')
        angleSlabs(size(PLATE.Subduction,1)+1:size(angleSlabs,1),:)     = []; %remove values
        depthShallowSlabDip(size(PLATE.Subduction,1)+1:size(angleSlabs,1),:) 	= [];
    elseif size(PLATE.Subduction,1)>size(angleSlabs,1)
        warning('#shallow-slab angles < #subduction trenches: Some shallow-slab angles could not be found.')
        dummy                   = size(angleSlabs,1);
        angleSlabs(dummy+1:size(PLATE.Subduction,1),:)      = NaN; %add values
        depthShallowSlabDip(dummy+1:size(PLATE.Subduction,1),:)	= NaN;
    end
    
    
    
    
    
%     if size(PLATE.Subduction,1)>1 %this disables multiple subduction detections - remove later.....................
%         disp(['   ',STYLE.SCHAR.indicationRightArrow,' multiple subduction-zone detections are prevented here!']);
%         PLATE.Subduction(dummy+1:end,:)              	= [];
%         PLATE.SubPolarity(dummy+1:end,:)             	= [];
%     end




    
    %check for spreading centres inbetween trench and slab position (if found remove trench tracker)
    jumpToNextSubductionZone = false;
    for isub=1:size(PLATE.Subduction,1)
        for ispread=1:size(PLATE.Spreading,1)
            if jumpToNextSubductionZone; continue; end
            %PLATE.slabDipPointsX(1,1) ONLY CHECKS FOR ONE SLAB POSITION YET
            if ~isempty(PLATE.Subduction)
                if SWITCH.DimensionalMode || SWITCH.DimensionalInput
                    if strcmp(GRID.Type,'spherical2D')
                        dummy = 1;
                    else
                        dummy = GRID.m2p;
                    end
                    SlabXlocation   = PLATE.slabDipPointsX(1,1)*dummy;  %convert [m] to [plot]
                else
                    SlabXlocation   = PLATE.slabDipPointsX(1,1);  %[nd]
                end
                if ~multipleSlabsDetected %only apply for more than one subduction zone (messes things up with PLATE.Subduction(isub,:) = NaN;)
                    if ( PLATE.Spreading(ispread,1)>PLATE.Subduction(isub,1) && PLATE.Spreading(ispread,1)<SlabXlocation ) || ...
                            ( PLATE.Spreading(ispread,1)<PLATE.Subduction(isub,1) && PLATE.Spreading(ispread,1)>SlabXlocation )
                        if SWITCH.Verbose; disp(['   ',STYLE.SCHAR.indicationRightArrow,' ridge inbetween trench and slab: trench tracker removed at x = ',num2str(PLATE.Subduction(isub,1),3)]); end
                        PLATE.Subduction(isub,:) 	= NaN;
                        PLATE.SubPolarity{isub,1}  	= 'na';
                        jumpToNextSubductionZone  	= true;
                    end
                end
            end
        end
        if jumpToNextSubductionZone
            jumpToNextSubductionZone     	= false;
            continue
        end
    end
    
    %general check for trench tracking
    if min(isnan(PLATE.Subduction(:)))==1
        couldNotBeTracked       = true;
    end
    if min(isnan(PLATE.Spreading(:)))==1
        couldNotBeTracked2      = true;
    end
    if strcmp(PLATE.Subduction,'noTracking')
        couldNotBeTracked       = true;
    end
    %prevent empty arrays
    if isempty(PLATE.Subduction)
        PLATE.Subduction                 	= NaN;
        PLATE.SubPolarity                 	= 'na';
        couldNotBeTracked                 	= true;
    end
end

%% SOME DIAGNOSTICS - FUNCTION OUTPUT
if strcmp(GRID.Dim,'2-D')
    if couldNotBeTracked
        %keep NaN values
    else
        for isub=1:size(PLATE.Subduction,1)
            if isub==2; warning('Plate diagnostics for multiple subduction zones is still in beta and not fully functional!'); end
            %UPPER/LOWER PLATE AND TRENCH VELOCITIES
            PLATE.SubPolarityNum                = zeros(size(PLATE.Subduction,1),1);
            PLATE.UPvelocity                    = zeros(size(PLATE.Subduction,1),1);
            PLATE.LPvelocity                    = PLATE.UPvelocity;
            PLATE.PlateConvergence              = PLATE.UPvelocity;
            PLATE.leftPvelocity                 = PLATE.UPvelocity;
            PLATE.rightPvelocity                = PLATE.UPvelocity;
            PLATE.PlateThicknessL(isub,1)       = PlateThicknessL(PLATE.Subduction(isub,1)==GRID.X_3Dp(:,1,1)); %plate thickness (T=1600 K) left of trench
            PLATE.PlateThicknessR(isub,1)       = PlateThicknessR(PLATE.Subduction(isub,1)==GRID.X_3Dp(:,1,1)); %plate thickness (T=1600 K) right of trench
            PLATE.PlateThicknessWideL(isub,1) 	= PlateThicknessWideL(PLATE.Subduction(isub,1)==GRID.X_3Dp(:,1,1)); %plate thickness (T=1600 K) left of trench
            PLATE.PlateThicknessWideR(isub,1)  	= PlateThicknessWideR(PLATE.Subduction(isub,1)==GRID.X_3Dp(:,1,1)); %plate thickness (T=1600 K) right of trench
            PLATE.ShallowSlabAngle            	= angleSlabs;           	%[degrees] - (#slab x 1)
            PLATE.ShallowSlabAngleDepth       	= depthShallowSlabDip *GRID.m2p;  	%[plotting dimension]
            PLATE.SlabViscosity                 = SLAB.viscosity;         	%[Pas] or [nd] mean value
            PLATE.SlabViscosityMin              = SLAB.viscosityMin;     	%[Pas] or [nd] min value
            PLATE.SlabMantleViscDiff            = SLAB.mantleViscContrast; 	%[nd]
            PLATE.SlabTemperatureMin            = SLAB.temperatureMin;     	%[K] or [nd] min value
            PLATE.SlabMantleTempDiff            = SLAB.mantleTempContrast; 	%[nd]
            PLATE.SlabDensityMax                = SLAB.densityMax;          %[kg/m^3] or [nd] max value
            PLATE.SlabStressMax                 = SLAB.stressMax;           %[MPa] or [nd] max value
            %these are already found above for general model structure
            if ~isnan(SLAB.mantleTemperature);  PLATE.UMantleTemperature    = SLAB.mantleTemperature;   end     %[K] or [nd]
            if ~isnan(SLAB.mantleDensity);      PLATE.UMantleDensity        = SLAB.mantleDensity;       end 	%[kg/m^3] or [nd]
            if ~isnan(SLAB.mantleViscosity);    PLATE.UMantleViscosity      = SLAB.mantleViscosity;     end   	%[Pas] or [nd]

            %if you add new variables, also add them as NaN on top!
        end
        
        for isub=1:size(PLATE.Subduction,1)
            dummy = find(GRID.X_3Dp(:,1,1)==PLATE.Subduction(isub,1)); %get index of current trench
            PLATE.idxTrench(isub,1)             = dummy;
            if strcmp(PLATE.SubPolarity(isub,1),'right')
                idxTrench                      	= min(dummy+numXpoints_critTrench,size(GRID.X_3Dp,1));
                idxUP                           = min(dummy+numXpoints_critUP,size(GRID.X_3Dp,1));
                idxLP                           = max(dummy-numXpoints_crit2,1);
                idxLeftP                        = idxLP;
                idxRightP                       = idxUP;
                PLATE.PlateThicknessUP(isub,1)  = PLATE.PlateThicknessWideR(isub,1); %plate thickness (T=1600 K)
                PLATE.PlateThicknessLP(isub,1) 	= PLATE.PlateThicknessL(isub,1); %plate thickness (T=1600 K)
            elseif strcmp(PLATE.SubPolarity(isub,1),'left')
                idxTrench                    	= max(dummy-numXpoints_critTrench,1);
                idxUP                           = max(dummy-numXpoints_critUP,1);
                idxLP                           = min(dummy+numXpoints_crit2,size(GRID.X_3Dp,1));
                idxLeftP                        = idxUP;
                idxRightP                       = idxLP;
                PLATE.PlateThicknessUP(isub,1)  = PLATE.PlateThicknessWideL(isub,1); %plate thickness (T=1600 K)
                PLATE.PlateThicknessLP(isub,1)  = PLATE.PlateThicknessR(isub,1); %plate thickness (T=1600 K)
            else %unknown subduction polarity
                idxTrench                    	= NaN;
                idxUP                           = NaN;
                idxLP                           = NaN;
                idxLeftP                        = max(dummy-3*numXpoints_crit2,1);
                idxRightP                       = min(dummy+3*numXpoints_crit2,size(GRID.X_3Dp,1));
                PLATE.PlateThicknessUP(isub,1)  = NaN;
                PLATE.PlateThicknessLP(isub,1)  = NaN;
            end
            if isnan(idxUP)
                PLATE.TrenchVelocity(isub,1)    = NaN;
                PLATE.UPvelocity(isub,1)        = NaN;
                PLATE.LPvelocity(isub,1)        = NaN;
            else
                PLATE.TrenchVelocity(isub,1)    = PLATE.PlateVelocity(idxTrench,1); %nd or [cm/a]
                PLATE.UPvelocity(isub,1)        = PLATE.PlateVelocity(idxUP,1); %nd or [cm/a]
                PLATE.LPvelocity(isub,1)        = PLATE.PlateVelocity(idxLP,1); %nd or [cm/a]
            end
            PLATE.idxTR(isub,1)                 = idxTrench;
            PLATE.idxUP(isub,1)             	= idxUP;
            PLATE.idxLP(isub,1)              	= idxLP;
            PLATE.leftPvelocity(isub,1)         = PLATE.PlateVelocity(idxLeftP,1); %nd or [cm/a]
            PLATE.rightPvelocity(isub,1)        = PLATE.PlateVelocity(idxRightP,1); %nd or [cm/a]
            PLATE.PlateConvergence(isub,1)      = abs( PLATE.leftPvelocity(isub,1)-PLATE.rightPvelocity(isub,1) ); %nd or [cm/a]
            %if you add new variables, also add them as NaN on top!
        end
    end
    %'1': left, '2': right, or '0': unknown
    PLATE.SubPolarityNum(strcmp(PLATE.SubPolarity,'left'))  = -1;
    PLATE.SubPolarityNum(strcmp(PLATE.SubPolarity,'right')) = 1;
    
    %LEFT/RIGHT SPREADING PLATE AND RIDGE VELOCITIES
    PLATE.ridgeLPvelocity      	= zeros(size(PLATE.Spreading,1),1);
    PLATE.ridgeRPvelocity    	= PLATE.ridgeLPvelocity;
    PLATE.PlateDivergence     	= PLATE.ridgeLPvelocity;
    if couldNotBeTracked2
        PLATE.ridgeLPvelocity               = NaN; %nd or [cm/a]
        PLATE.ridgeRPvelocity               = NaN; %nd or [cm/a]
        PLATE.PlateDivergence               = NaN; %nd or [cm/a]
    else
        for ispr=1:size(PLATE.Spreading,1)
            dummy = find(GRID.X_3Dp(:,1,1)==PLATE.Spreading(ispr,1)); %get index of ridge(s)
            idx_rLP = max(dummy-numXpoints_crit2,1);
            idx_rRP = min(dummy+numXpoints_crit2,size(GRID.X_3Dp,1));
            
            PLATE.ridgeLPvelocity(ispr,1)  	= PLATE.PlateVelocity(idx_rLP,1); %nd or [cm/a]
            PLATE.ridgeRPvelocity(ispr,1) 	= PLATE.PlateVelocity(idx_rRP,1); %nd or [cm/a]
            PLATE.PlateDivergence(ispr,1) 	= abs( PLATE.ridgeLPvelocity(ispr,1)-PLATE.ridgeRPvelocity(ispr,1) ); %nd or [cm/a]
        end
    end
    PLATE.RidgeVelocity      	= (PLATE.ridgeLPvelocity(:,1)+PLATE.ridgeRPvelocity(:,1))./2;    %assuming perfectly symmetric spreading velocities!
    
    
    %% VISCOUS BENDING DISSIPATION
    BendingDissipationFormulation = 'Buffett2006viscoplastic';    %choose which formulation: 'Conrad1999viscous','Buffett2006viscous' or 'Buffett2006viscoplastic'
    
    % VISCOUS BENDING DISSIPATION for a viscous plate (after Conrad & Hager 1999)
    %     Cl          = 2.0;        %constant
    %     vplate      = 1;          %lower plate velocity
    %     etalitho    = 1;          %lithosphere-mantle viscosity difference
    %     hs          = 1;          %lower plate/slab thickness
    %     R           = 1;          %bending radius
    %
    %     vdLitho     = Cl * vplate^2 * etalitho * (hs/R)^3;
    Cl               	= 2.0;      %constant (after Conrad & Hager 1999)
    % depthShallowSlabDip;                 %[m] or [nd]
    if BENDING.failed || isempty(BENDING.radiusSI)
        %use less good method instead
        BENDING.radiusSI        = abs(depthShallowSlabDip./(tan(PLATE.ShallowSlabAngle./2).*sin(90-PLATE.ShallowSlabAngle))); %[m] or [nd] - see sketch Crameri et al. (2016)
        warning('Bending radius is calculated with second method! Check here.')
    end
    if SWITCH.DimensionalMode || SWITCH.DimensionalInput
        h_lp            = PLATE.PlateThicknessLP ./GRID.m2p;                %[plotting dimension]->[m]
        vplate          = abs( PLATE.LPvelocity ./(100*SETUP.secyear) );    %[cm/a]->[m/s], conversion to SI-units
    else
        h_lp            = PLATE.PlateThicknessLP;
        vplate          = abs( PLATE.LPvelocity );
    end
    if strcmp(BendingDissipationFormulation,'Conrad1999viscous')
        PLATE.ViscDissBending 	= Cl .* vplate.^2 .* PLATE.SlabMantleViscDiff .* (h_lp./BENDING.radiusSI).^3; %[m^2/s^2]  char.Value=...
    end
    
    % VISCO-PLASTIC BENDING DISSIPATION for a viscous plate (after Buffett 2006, JGR)
    %     Cl                = 1/6;      %constant
    sigma_y_ductile         = 1000e6;     %Pa; ductile yield stress (i.e., in the centre of the plate)
    if SWITCH.DimensionalMode || SWITCH.DimensionalInput
        %sigma_y_ductile         = PLATE.coreStressMax ./SETUP.stressscale;          %[MPa] -> [Pa], conversion to SI-units
        sigma_y_ductile         = PLATE.SlabStressMax ./SETUP.stressscale;          %[MPa] -> [Pa], conversion to SI-units
    else
        %sigma_y_ductile         = PLATE.coreStressMax;
        sigma_y_ductile         = PLATE.SlabStressMax;
    end
    %
    %     vdLitho           = 1/6 * vplate * sigma_y(0) * (hs^2/R);
    %     or: vdLitho       = 2/3 * vplate^2 * etaLithoEffective * (hs/R)^3
    if strcmp(BendingDissipationFormulation,'Buffett2006viscous')
        PLATE.ViscDissBending 	= 2/3 .* vplate.^2 .* PLATE.SlabViscosityMin .* (h_lp./BENDING.radiusSI).^3; %[m*kg/s^3]=[N/s]=[W/m]  char.Value=2.1649e+03  dissipation rate
    end
  
    % VISCO-PLASTIC BENDING DISSIPATION for a visco-plastic plate (after Buffett 2006, JGR - equation 42)
    if strcmp(BendingDissipationFormulation,'Buffett2006viscoplastic')
        PLATE.ViscDissBending 	= 1/6 .* vplate .* sigma_y_ductile .* (h_lp.^(2)./BENDING.radiusSI);  %[m*kg/s^3]=[N/s]=[W/m]  char.Value=874.3859  dissipation rate
    end
    
    
    %% RELATIVE BENDING DISSIPATION
    normFactorBending           = 874.3859;             %characteristic dissipation value [N/s]=[W/m]
    PLATE.ViscDissBendingRel    = PLATE.ViscDissBending./normFactorBending; %[nd]
    PLATE.BendingRadius         = BENDING.radius;       %[plotting dimension] or [nd]
    
    % PREVENT EMPTY ARRAYS
    if isempty(PLATE.ViscDissBending);  	PLATE.ViscDissBending       = NaN; end
    if isempty(PLATE.ViscDissBendingRel);  	PLATE.ViscDissBendingRel    = NaN; end
    
    %% SUBDUCTION (VOLUMETRIC) FLOW-RATE
    PLATE.subductionFlowRate    = h_lp.*vplate;       	%SI-units [m^2/s]    maybe use weight instead, i.e., give in [kg/a] or something -> present-day slab flux is ?10^15 kg/a or in volume 300 km^3/a
    PLATE.subductionFlowRate    = PLATE.subductionFlowRate*SETUP.secyear;   %[m^2/a]
    
    %% THEORETIC SLAB AGE
    PLATE.shallowSlabAge        = ((h_lp./2.32).^2)./SETUP.kappa ./SETUP.secyear./1e6; %[Ma]
    
    %% THEORETIC AMOUNT OF WATER REMAINING IN THE SLAB (WATER RETENTION)
    % W = 1.06*v_slab +0.14*age_slab -0.023*T_mantle +17    ,  [1e5*kg/m^2] = [cm/a],[Ma],[degreeC] (after Magni et al., 2014, G-cubed)
    PLATE.waterRemainedInSlab   = (1.06*PLATE.SlabSinkingVelocity +0.14*PLATE.shallowSlabAge -0.023*(PLATE.UMantleTemperature-273.15) +17)*1e5; %[kg/m^2]
    if PLATE.waterRemainedInSlab<0.0
        PLATE.waterRemainedInSlab = 0.0;
    end
    
elseif strcmp(GRID.Dim,'3-D')
    PLATE.SubLength             = SubLength;
    PLATE.SprLength             = SprLength;
    PLATE.SubLengthTotal        = SubLengthTotal;
    PLATE.SprLengthTotal        = SprLengthTotal;
    %add here more diagnostics & function output......
    if strcmp(GRID.Type,'yinyang')
        PLATE.SubLength_yang  	= SubLength_yang;
        PLATE.SprLength_yang   	= SprLength_yang;
    end
end

%% SUBDUCTION-TOPOGRAPHY CHARACTERISTICS
if ~couldNotBeTracked && strcmp(GRID.Dim,'2-D')
    [SUBTOPO] = f_SubTopoCharacteristics(PLATE,GRID,SWITCH,FILE,STYLE);
    
    PLATE.TopoOverridingPlate     	= SUBTOPO.overridingPlate;      %in [plotting dimensions]
    PLATE.TopoSubductingPlate     	= SUBTOPO.subductingPlate;      %in [plotting dimensions]
    PLATE.TrenchDepth               = SUBTOPO.TrenchZ;              %measurement points in [plotting dimensions]
    PLATE.TrenchDepthX             	= SUBTOPO.TrenchX;              %measurement points in [plotting dimensions]
    PLATE.BackArcBasinX          	= SUBTOPO.BackArcBasinX;        %measurement points in [plotting dimensions]
    PLATE.BackArcBasinZ           	= SUBTOPO.BackArcBasinZ;        %measurement points in [plotting dimensions]
    PLATE.BackArcBasinExtentX       = SUBTOPO.BackArcBasinExtentX;  %measurement points in [plotting dimensions]
    PLATE.BackArcBasinExtentZ       = SUBTOPO.BackArcBasinExtentZ;  %measurement points in [plotting dimensions]
    PLATE.InundationX               = SUBTOPO.InundationX;          %measurement points in [plotting dimensions]
    PLATE.InundationZ               = SUBTOPO.InundationZ;          %measurement points in [plotting dimensions]
    PLATE.InundationDistance        = SUBTOPO.InundationDistance;   %distance from trench in [plotting dimensions]
    PLATE.IslandArcX                = SUBTOPO.IslandArcX;          	%measurement points in [plotting dimensions]
    PLATE.IslandArcZ                = SUBTOPO.IslandArcZ;        	%measurement points in [plotting dimensions]
    PLATE.ForeBulgeX                = SUBTOPO.ForeBulgeX;        	%measurement points in [plotting dimensions]
    PLATE.ForeBulgeZ                = SUBTOPO.ForeBulgeZ;        	%measurement points in [plotting dimensions]
    PLATE.BackArcBasinArea       	= SUBTOPO.BackArcBasinArea;   	%area in [plotting dimensions^2]
    
    PLATE.UPtiltAngleXnear          = SUBTOPO.UPtiltAngleXnear;   	%measurement points in [plotting dimensions]
    PLATE.UPtiltAngleXfar           = SUBTOPO.UPtiltAngleXfar;
    PLATE.UPtiltAngleZnear          = SUBTOPO.UPtiltAngleZnear;
    PLATE.UPtiltAngleZfar           = SUBTOPO.UPtiltAngleZfar;
    PLATE.UPtiltAngle               = SUBTOPO.UPtiltAngle;        	%UP-tilt angle in [degrees]
end

%% PROBLEM CHECKS
if strcmp(GRID.Dim,'2-D')
    if size(PLATE.ShallowSlabAngle,1)~=size(PLATE.Subduction,1)
        warning('#shallow-slab angles does not match #subduction zones. Check here!');
    end
end

%% DISPLAY INFOS
if SWITCH.DimensionalMode || SWITCH.DimensionalInput
    dim_dist        = GRID.Xdim;	dimfac_dist = 1e-3; %*GRID.dimFactor
    dim_vdiss       = 'W/m';   %dissipation dimension
    dimFlowRate     = 'm^2/a'; %flux in 2-D
else %non-dimensional
    dim_dist        = ['[',GRID.Xdim,']'];       dimfac_dist = 1;
    dim_vdiss       = '[nd]';
    dimFlowRate     = '[nd]';
end
%general displays
if ~strcmp(GRID.Type,'yinyang') %not fully adjusted yet for yinyang:
    disp('   Plate diagnostics');
    disp(['     ',STYLE.SCHAR.smallBullet,' plate mobility        = ',num2str(PLATE.PlateMobility,3)]);
    disp(['     ',STYLE.SCHAR.smallBullet,' plate drivity         = ',num2str(PLATE.PlateDrivety,3)]);
end

if strcmp(GRID.Dim,'2-D')
    %subduction
    if isnan(PLATE.Subduction(1,1)) || strcmp(PLATE.Subduction,'noTracking')
        disp('   No Active Subduction Zone Detected.');
    else
        if size(PLATE.Subduction,1)==1
            disp('   1 Active Subduction Zone Detected:');
        else
            disp(['   ',num2str(size(PLATE.Subduction,1)),' Active Subduction Zones Detected:']);
        end
        if SWITCH.Verbose; disp(['     ',STYLE.SCHAR.indicationRightArrow,' plate boundaries tracked between ',num2str(minSurfaceLevel,2),' and ',num2str(maxSurfaceLevel,2),' ',dim_dist,' depth']); end
        disp('     Geometry');
        disp(['     ',STYLE.SCHAR.smallBullet,' slab polarity         = ',PLATE.SubPolarity{1,1}]);
        disp(['     ',STYLE.SCHAR.smallBullet,' shallow-slab dip      = ',num2str(PLATE.ShallowSlabAngle',2),STYLE.SCHAR.degree,' (at ',num2str(PLATE.ShallowSlabAngleDepth',3),' ',dim_dist,' depth)']);
        disp(['     ',STYLE.SCHAR.smallBullet,' slab-tip depth        = ',num2str(PLATE.SlabTipPosition(1,2)',4),' ',dim_dist]);
        disp(['     ',STYLE.SCHAR.smallBullet,' slab-tip dip          = ',num2str(PLATE.SlabTipAngle',2),STYLE.SCHAR.degree]);
        disp('     Thickness');
        if ~isnan(PLATE.PlateThicknessLP(1,1)) %only checks first subduction zone
            disp(['     ',STYLE.SCHAR.smallBullet,' lower plate           = ',num2str(PLATE.PlateThicknessLP',3),' ',dim_dist]);
            disp(['     ',STYLE.SCHAR.smallBullet,' upper plate           = ',num2str(PLATE.PlateThicknessUP',3),' ',dim_dist]);
        else %unknown polarity
            disp(['     ',STYLE.SCHAR.smallBullet,' left plate            = ',num2str(PLATE.PlateThicknessL',3),' ',dim_dist]);
            disp(['     ',STYLE.SCHAR.smallBullet,' right plate           = ',num2str(PLATE.PlateThicknessR',3),' ',dim_dist]);
        end
        disp('     Velocity');
        disp(['     ',STYLE.SCHAR.smallBullet,' plate (RMS)           = ',num2str(PLATE.PlateVelocityRMS',2),' ',SETUP.vDim,' (',num2str(PLATE.PlateVelocityMin,2),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLATE.PlateVelocityMax,2),' ',SETUP.vDim,')']);
        disp(['     ',STYLE.SCHAR.smallBullet,' upper plate           = ',num2str(PLATE.UPvelocity',2),' ',SETUP.vDim]);
        disp(['     ',STYLE.SCHAR.smallBullet,' lower plate           = ',num2str(PLATE.LPvelocity',2),' ',SETUP.vDim]);
        disp(['     ',STYLE.SCHAR.smallBullet,' plate convergence     = ',num2str(PLATE.PlateConvergence',2),' ',SETUP.vDim]);
        disp(['     ',STYLE.SCHAR.smallBullet,' slab sinking          = ',num2str(PLATE.SlabSinkingVelocity',2),' ',SETUP.vDim]);
        disp(['     ',STYLE.SCHAR.smallBullet,' slab-tip sinking      = ',num2str(-PLATE.SlabTipVZ',2),' ',SETUP.vDim]);
        disp(['     ',STYLE.SCHAR.smallBullet,' trench                = ',num2str(PLATE.TrenchVelocity',2),' ',SETUP.vDim,' (theoretic: ',num2str(PLATE.theoreticTrenchVelocity',2),' ',SETUP.vDim,')']);
        disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle (max)    = ',num2str(PLATE.UMantleVelocityMax,2),' (total), ',num2str(PLATE.UMantleHVelocityMax,2),' (horizontal), ',num2str(PLATE.UMantleRVelocityMax,2),' (radial) '...
            ,SETUP.vDim,' (',num2str(GRID.Z_3Dp(1,1,zIdxCheck2),3),' ',GRID.Xdim,' depth)']);
        disp('     Viscosity');
        disp(['     ',STYLE.SCHAR.smallBullet,' plate core            = ',num2str(PLATE.coreViscosityMin,2),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLATE.coreViscosityMax,2),' ',SETUP.etaDim]);
        disp(['     ',STYLE.SCHAR.smallBullet,' slab (mean)           = ',num2str(PLATE.SlabViscosity',2),' ',SETUP.etaDim,' (',num2str(SLAB.DepthChecked(1,1),3),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(SLAB.DepthChecked(1,2),3),' ',dim_dist,' depth)']);
        disp(['     ',STYLE.SCHAR.smallBullet,' slab (min)            = ',num2str(PLATE.SlabViscosityMin',2),' ',SETUP.etaDim,' (',num2str(SLAB.DepthChecked(1,1),3),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(SLAB.DepthChecked(1,2),3),' ',dim_dist,' depth)']);
        disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(PLATE.UMantleViscosity',2),' ',SETUP.etaDim,' (',num2str(GRID.Z_3Dp(1,1,zIdxCheckUM),3),' ',dim_dist,' depth)']);
        disp(['     ',STYLE.SCHAR.smallBullet,' slab-mantle contrast  = ',num2str(PLATE.SlabMantleViscDiff',3)]);
        disp('     Stress');
        disp(['     ',STYLE.SCHAR.smallBullet,' plate                 = ',num2str(PLATE.maxStress,3),' ',SETUP.stressDim]);
        disp(['     ',STYLE.SCHAR.smallBullet,' plate core            = ',num2str(PLATE.coreStressMin,3),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLATE.coreStressMax,3),' ',SETUP.stressDim]);
        disp('     Strainrate');
        disp(['     ',STYLE.SCHAR.smallBullet,' plate core            = ',num2str(PLATE.coreStrainrateMin,2),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLATE.coreStrainrateMax,2),' ',SETUP.edotDim2]);
        disp('     Bending');
        disp(['     ',STYLE.SCHAR.smallBullet,' radius                = ',num2str(BENDING.radiusSI'.*dimfac_dist,3),' ',dim_dist]);
        disp(['     ',STYLE.SCHAR.smallBullet,' dissipation           = ',num2str(PLATE.ViscDissBending',3),' ',dim_vdiss]);
        disp(['     ',STYLE.SCHAR.smallBullet,' rel. dissipation      = ',num2str(PLATE.ViscDissBendingRel',3),' (to ',num2str(normFactorBending,3),' W/m)']);
        if SUBTOPO.Successful
            disp('     Topography');
            disp(['     ',STYLE.SCHAR.smallBullet,' fore-bulge height     = ',num2str(PLATE.ForeBulgeZ,3),' ',dim_dist]);
            disp(['     ',STYLE.SCHAR.smallBullet,' trench depth          = ',num2str(PLATE.TrenchDepth,3),' ',dim_dist]);
            disp(['     ',STYLE.SCHAR.smallBullet,' island-arc height     = ',num2str(PLATE.IslandArcZ,3),' ',dim_dist]);
            disp(['     ',STYLE.SCHAR.smallBullet,' back-arc basin depth  = ',num2str(PLATE.BackArcBasinZ,3),' ',dim_dist]);
            disp(['     ',STYLE.SCHAR.smallBullet,' back-arc basin extent = ',num2str(PLATE.BackArcBasinExtentX,3),' ',dim_dist]);
            disp(['     ',STYLE.SCHAR.smallBullet,' back-arc basin area   = ',num2str(PLATE.BackArcBasinArea,3),' ',dim_dist,'^2']);
            disp(['     ',STYLE.SCHAR.smallBullet,' upper-plate tilt      = ',num2str(PLATE.UPtiltAngle',2),STYLE.SCHAR.degree]);
            disp(['     ',STYLE.SCHAR.smallBullet,' inundation            = ',num2str(PLATE.InundationDistance,3),' ',dim_dist]);
        end
        disp('     Other');
        disp(['     ',STYLE.SCHAR.smallBullet,' subduction flow-rate  = ',num2str(PLATE.subductionFlowRate',3),' ',dimFlowRate]);
        disp(['     ',STYLE.SCHAR.smallBullet,' shallow-slab age      = ',num2str(PLATE.shallowSlabAge',3),' Ma']);
        disp(['     ',STYLE.SCHAR.smallBullet,' visc. plate dissip.   = ',num2str(PLATE.ViscDissipationPlate',3),' ',dim_vdiss,' (',num2str(PLATE.ViscDissipationPlateNoCrust',3),' ',dim_vdiss,' without crustal areas)']);        
        disp(['     ',STYLE.SCHAR.smallBullet,' water retention       = ',num2str(PLATE.waterRemainedInSlab',3),' kg/m^2']);
    end    
    %spreading
    if isnan(PLATE.Spreading(1,1)) || strcmp(PLATE.Subduction,'noTracking')
        disp('   No Active Spreading Center Detected.');
    else
        if size(PLATE.Spreading,1)==1
            disp('   1 Active Spreading Center Detected:');
        else
            disp(['   ',num2str(size(PLATE.Spreading,1)),' Active Spreading Centers Detected:']);
        end
        disp('     Velocity');
        disp(['     ',STYLE.SCHAR.smallBullet,' plate divergence      = ',num2str(PLATE.PlateDivergence',3),' ',SETUP.vDim]);
    end
    
elseif strcmp(GRID.Dim,'3-D')
    if strcmp(GRID.Type,'Cartesian')
        dummySub = PLATE.Subduction;
        dummySpr = PLATE.Spreading;
    elseif strcmp(GRID.Type,'yinyang')
        if strcmp(PLATE.Subduction,'noTracking') && strcmp(PLATE.Subduction_yang,'noTracking')
            dummySub = 'noTracking';
        elseif strcmp(PLATE.Subduction,'noTracking') && ~strcmp(PLATE.Subduction_yang,'noTracking')
            PLATE.Subduction = [];
            dummySub = PLATE.Subduction_yang;
        elseif ~strcmp(PLATE.Subduction,'noTracking') && strcmp(PLATE.Subduction_yang,'noTracking')
            PLATE.Subduction_yang = [];
            dummySub = PLATE.Subduction;
        else
            dummySub = [PLATE.Subduction; PLATE.Subduction_yang];
        end
        dummySpr = [PLATE.Spreading; PLATE.Subduction_yang];
    end
    %subduction
    if isempty(dummySub) || strcmp(dummySub,'noTracking')
        disp('   No Active Subduction Zone Detected.');
    else
        if size(dummySub,1)==1
            disp('   1 Active Subduction Zone Detected:');
        else
            disp(['   ',num2str(size(dummySub,1)),' Active Subduction Zones Detected:']);
        end
        disp(['     ',STYLE.SCHAR.smallBullet,' trench length         = ',num2str(PLATE.SubLengthTotal,3),' ',dim_dist]);
    end
    %spreading
    if isempty(dummySpr) || strcmp(dummySub,'noTracking')
        disp('   No Active Spreading Center Detected.');
    else
        if size(dummySpr,1)==1
            disp('   1 Active Spreading Center Detected:');
        else
            disp(['   ',num2str(size(dummySpr,1)),' Active Spreading Centers Detected:']);
        end
        if SWITCH.DimensionalMode || SWITCH.DimensionalInput; dimString = GRID.Xdim; else; dimString = ['[',GRID.Xdim,']']; end
        disp(['     ',STYLE.SCHAR.smallBullet,' ridge length          = ',num2str(PLATE.SprLengthTotal,3),' ',dimString]);
    end
    
    %GET MAX. CONVERGENCE/DIVERGENCE VELOCITY
    
    disp('   Additional Plate Diagnostics');
    disp('     Viscosity');
    disp(['     ',STYLE.SCHAR.smallBullet,' plate core            = ',num2str(PLATE.coreViscosityMin,2),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLATE.coreViscosityMax,2),' ',SETUP.etaDim]);
    disp(['     ',STYLE.SCHAR.smallBullet,' upper mantle          = ',num2str(PLATE.UMantleViscosity',2),' ',SETUP.etaDim,' (',num2str(GRID.Z_3Dp(1,1,zIdxCheckUM),3),' ',dim_dist,' depth)']);
    disp('     Stress');
    disp(['     ',STYLE.SCHAR.smallBullet,' plate                 = ',num2str(PLATE.maxStress,3),' ',SETUP.stressDim]);
    disp(['     ',STYLE.SCHAR.smallBullet,' plate core            = ',num2str(PLATE.coreStressMin,3),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLATE.coreStressMax,3),' ',SETUP.stressDim]);
    disp('     Strainrate');
    disp(['     ',STYLE.SCHAR.smallBullet,' plate core            = ',num2str(PLATE.coreStrainrateMin,2),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(PLATE.coreStrainrateMax,2),' ',SETUP.edotDim2]);
    disp('     Other');
    disp(['     ',STYLE.SCHAR.smallBullet,' visc. plate dissip.   = ',num2str(PLATE.ViscDissipationPlate',3),' ',dim_vdiss,' (',num2str(PLATE.ViscDissipationPlateNoCrust',3),' ',dim_vdiss,' without crustal areas)']);

    %add here display infos......
    
end

end



% PLATE.StagnantLid         %true or false


%%for 2-D:

% PLATE.Subduction          in [GRID.Xdim], i.e., plotting values
% PLATE.SubPolarity         {'right'} or {'left'} or {'na'}
% PLATE.SubPolarityNum      %'-1': left, '1': right, or '0': unknown
% PLATE.ShallowSlabAngle   	%[degrees]
% PLATE.SlabViscosity       %[Pas] or [nd]
% PLATE.SlabViscosityMin    %[Pas] or [nd]
% PLATE.UMantleViscosity    %[Pas] or [nd]
% PLATE.SlabMantleViscDiff  %[nd]
% PLATE.ViscDissBending     %[N/s]
% PLATE.ViscDissBendingRel  %[-]
% PLATE.BendingRadius       %[plotting dimension]

% PLATE.PlateThicknessL     %nd or [plotting dimension]
% PLATE.PlateThicknessR     %nd or [plotting dimension]
% PLATE.PlateThicknessLP    %nd or [plotting dimension]
% PLATE.PlateThicknessUP    %nd or [plotting dimension]

% PLATE.PlateDivergence     %nd or [cm/a]
% PLATE.PlateConvergence    %nd or [cm/a]
% PLATE.TrenchVelocity      %nd or [cm/a]  %assuming perfect single-sided subduction
% PLATE.UPvelocity(ii,1)	%nd or [cm/a]
% PLATE.LPvelocity(ii,1) 	%nd or [cm/a]
% PLATE.leftPvelocity(ii,1)	%nd or [cm/a]
% PLATE.rightPvelocity(ii,1)%nd or [cm/a]
% PLATE.RidgeVelocity       %nd or [cm/a]  %assuming perfectly symmetric spreading velocities!
% PLATE.ridgeLPvelocity   	%nd or [cm/a]
% PLATE.ridgeRPvelocity   	%nd or [cm/a]


%%for 3-D:

% PLATE.Subduction          %[plotting dimension] - [iarea,allpoints,x_i] or 'noTracking'
% PLATE.SubductionOut   	%[plotting dimension] - [iarea,allpoints,x_i]
% PLATE.Spreading         	%[plotting dimension] - [iarea,allpoints,x_i]
% PLATE.SpreadingOut      	%[plotting dimension] - [iarea,allpoints,x_i]

% PLATE.SubLength           %nd or [m] - [iarea,d]
% PLATE.SprLength           %nd or [m] - [iarea,d]
% PLATE.SubLengthTotal   	%nd or [m] - d
% PLATE.SprLengthTotal   	%nd or [m] - d

% PLATE.SubCenterP          %[plotting dimension] - [iarea,points,x_i]
% PLATE.SprCenterP          %[plotting dimension] - [iarea,points,x_i]

% PLATE.UpperPlateP         %unit vector - [iarea,points,x_i]
% PLATE.LowerPlateP         %unit vector - [iarea,points,x_i]

% PLATE.SubNormalP          %unit vector - [iarea,points,x_i] - directed towards subducting plate
% PLATE.SubStrikeP          %unit vector - [iarea,points,x_i]
% PLATE.SprNormalP          %unit vector - [iarea,points,x_i] - random direction
% PLATE.SprStrikeP          %unit vector - [iarea,points,x_i]

%% for YinYang:
% see above plus a corresponding suite of xxx_yang variables






%% INTERNAL FUNCTIONS

%%                                         SLAB-TIP POSITION AND ANGLE 2.2
%                                                Fabio Crameri, 08.03.2019
function [SLABTIP] = f_SlabTipDiagnostics(T_3D,PLATE,GRID,SWITCH)
%% notes
% set limit of depth and width of the line fit...........
% use SLABTIP.failed = true if diagnostics is not possible

%% function input
TempThreshold                   = PLATE.SlabContourTemperature; %define threshold [K]

% defaults
SLABTIP.failed                  = false;
SLABTIP.SlabTipPosition         = [NaN,NaN];
SLABTIP.SlabTipIndicesXZ        = [NaN,NaN];    %[idx_x,idx_z]
SLABTIP.SlabTipAngle            = NaN;
SLABTIP.SlabTipAnglePoint1      = [NaN,NaN];    %[x,z]
SLABTIP.SlabTipAnglePoint2      = [NaN,NaN];    %[x,z]

%% problem checks
if strcmp(GRID.Dim,'3-D')
    return
    
end

%% get temperature contour of plate and slab
% find values below, at & above the contour value
BW = zeros(size(T_3D,1),size(T_3D,3));
BW(T_3D(:,1,:)<=TempThreshold & T_3D(:,1,:)>min(T_3D(:)))   = 1; %below or at

%% point cloud representing subducting plate
%make sure points close to each other are connected
% BW = bwmorph(BW,'thicken',2);
BW = bwmorph(BW,'bridge');
%plot(GRID.x2dp(BW),GRID.z2dp(BW),'.')

%find connected points with deepest reach
CC          = bwconncomp(BW);
maxZ = zeros(size(CC.PixelIdxList,2),1);
minZ = maxZ;
slabFound = maxZ;
for i=1:size(CC.PixelIdxList,2)
    %skip isolated small areas
    if size(CC.PixelIdxList{1,i},1)<10 %specify small area criterion <-----
        continue
    end
    maxZ(i,1) = max(GRID.z2dp(CC.PixelIdxList{1,i}));
    minZ(i,1) = min(GRID.z2dp(CC.PixelIdxList{1,i}));
    %skip areas separated from surface
    if minZ(i,1)>120e3*GRID.m2p %connection to surface criterion <---------
        continue
    end
    slabFound(i,1) = 1;
end
maxZ(~slabFound) = 0;
% minZ(~slabFound) = 0;
if max(maxZ)==0; error('Plate geometry could not be found! Check here.'); end
[~,idx] = max(maxZ);
%remove all others
for i=1:size(CC.PixelIdxList,2)
    if i~=idx
        BW(CC.PixelIdxList{1,i}) = 0;
    end
end

%% remove shallow parts
if ~isempty(PLATE.LABdepth) || ~isnan(PLATE.LABdepth)
    BW(GRID.z2dp<PLATE.LABdepth)    = 0;
end

%% get coordinates
if strcmp(GRID.Type,'Cartesian')
    x               = GRID.x2dp(BW);
elseif strcmp(GRID.Type,'spherical2D')
    x_dist        	= GRID.x2dsp(BW); %use distance
    x           	= GRID.x2dp(BW); %use degrees
end
z                   = GRID.z2dp(BW);
if isempty(x) || isempty(z)
    SLABTIP.failed                  = true;
    SLABTIP.SlabTipPosition         = [NaN,NaN];
    
    return
end

%% subducting plate contour
%find slab polarity
ApproxHorizPosShallow       = mean(x(min(z)==z));
ApproxHorizPosDeep          = mean(x(max(z)==z));
if ApproxHorizPosShallow>ApproxHorizPosDeep         %left-dipping slab
    %slabPolarity            = 'left';
    SlabTipPosition(1,1)	= min(x(max(z)==z));        %x
elseif ApproxHorizPosShallow<ApproxHorizPosDeep 	%right-dipping slab
    %slabPolarity            = 'right';
    SlabTipPosition(1,1)	= max(x(max(z)==z));        %x
else
    if SWITCH.Verbose; warning('plate polarity could not be found!'); end
    SlabTipPosition(1,1)	= mean(x(max(z)==z));       %x
end
SlabTipPosition(1,2)        = max(z);                   %z
SLABTIP.SlabTipPosition     = SlabTipPosition;              %plotting dimensions or [deg]

if logical(0)
    figure(2);clf
    plot(GRID.x2dp(BW)',GRID.z2dp(BW)','og')
    axis ij
    figure(1)
end

%% slab-tip x and z indices
[dummy,idxX]   = min(abs(GRID.X_3Dp(:,1,1)-SlabTipPosition(1,1))); %slab tip position index X
[dummy,idxZ]   = min(abs(GRID.Z_3Dp(1,1,:)-SlabTipPosition(1,2))); %slab tip position index Z
SLABTIP.SlabTipIndicesXZ    = [idxX,idxZ];

%% slab-tip angle
try
    ApproxSlabThickness         = PLATE.LABdepth;
    Radius2CheckAroundSlabTip   = 3 *ApproxSlabThickness;    %<<<<<<<<<< define checking area (i.e., length of slab segment to take angle of)
    
    % get point cloud of slab-tip portion
    SlabPoints          = [x,z];
    idx                 = rangesearch(SlabPoints,SlabTipPosition,Radius2CheckAroundSlabTip); %extract indices of close slab portions for each slab tip
    SlabTipPoints       = SlabPoints(idx{:},:)';
    
    % create min. volume ellipse around slab tip
    [A,c] = MinVolEllipse(SlabTipPoints,0.01);
    
    %FIND AREA STRIKE / NORMAL
    N = 20; % Default value for ellipse #gridpoints
    % "singular value decomposition" to extract the orientation and the
    % axes of the ellipsoid
    [~,D,V] = svd(A);
    % get the major and minor axes
    a = 1/sqrt(D(1,1)); %minor
    b = 1/sqrt(D(2,2)); %major
    theta = 0:1/N:2*pi+1/N;
    % Parametric equation of the ellipse
    state(1,:) = a*cos(theta);
    state(2,:) = b*sin(theta);
    X = V * state; %coordinate transform to obtain unit (x,y) vector of ellipse
    %radius
    R = sqrt(X(1,:).^2+X(2,:).^2);
    %unit vector minor axis
    unitVecNormal = X(:,R==min(R))./min(R); %unit pos-values of minor axis with min(R)==1
    %unitvector major axis
    unitVecStrike = X(:,R==max(R))./max(R); %unit pos-values of major axis with max(R)==1
    
    tipAnglePoint1     = [c(1,:)+b*unitVecStrike(1,:), c(2,:)+b*unitVecStrike(2,:)];
    tipAnglePoint2     = [c(1,:)-b*unitVecStrike(1,:), c(2,:)-b*unitVecStrike(2,:)];
    
    if strcmp(GRID.Type,'spherical2D')
        dx              = abs( tipAnglePoint2(:,1)-tipAnglePoint1(:,1) );
        dz              = abs( tipAnglePoint2(:,2)-tipAnglePoint1(:,2) );
        DepthOfAngle    = max(tipAnglePoint2(:,2),tipAnglePoint1(:,2)) -dz/2;
        RadiusOfAngle   = (max(GRID.Z_3Dp(:))-DepthOfAngle) +GRID.rcmb_p;
        dx              = dx.*RadiusOfAngle;        %rad to distance
        angleSlabTips   = atand(dz./dx);            %[degrees] - tan(alpha) = dz/dx
    else
        % angleSlabTips              = atand(dz./dx);                    %[degrees] - tan(alpha) = dz/dx
        angleSlabTips     	= atand(abs(unitVecStrike(2,:)./unitVecStrike(1,:))); 	%[degrees] - tan(alpha) = dz/dx
    end
    
    if logical(0)
        figure(2);clf
        plot(GRID.x2dp(BW)',GRID.z2dp(BW)','og','MarkerFaceColor','g')
        hold on
        plot(SlabTipPoints(1,:),SlabTipPoints(2,:),'or')
        axis ij
        if ~strcmp(GRID.Type,'spherical2D')
            axis equal
        end
        
        hold on
        plot(c(1),c(2),'ro')
        hold on
        Ellipse_plot(A, c);
        
        hold on
        plot([tipAnglePoint1(:,1),tipAnglePoint2(:,1)],[tipAnglePoint1(:,2),tipAnglePoint2(:,2)],'k')
        
        figure(1)
    end
    SLABTIP.SlabTipAngle            = angleSlabTips;
    SLABTIP.SlabTipAnglePoint1      = tipAnglePoint1;
    SLABTIP.SlabTipAnglePoint2      = tipAnglePoint2;
    
catch errMe
    warning(errMe.message)
    warning off backtrace
    disp(' ')
    warning(['   ',STYLE.SCHAR.indicationRightArrow,' Slab tip angle could not be determined.'])
    warning on backtrace
end

end




%%                                                      BENDING RADIUS 1.32
%                                                Fabio Crameri, 19.09.2017
function [BENDING] = f_BendingRadiusB(T_3D,PLATE,GRID,SWITCH,LITHO)
%   . calls f_plotCircle(x,y,r)

%% Notes
% use BENDING.failed = true if diagnostics is not possible

%% function input
TempThreshold               = PLATE.BaseTemperature; %define threshold [K]

%    v---------------------------------------------------- empirically derived value! (might vary depending on grid resolution!)
smoothingParam              = 1/(1+1.5*mean(GRID.dxp(:))^3*48);	%2e-6 Smoothing of fitted line: 0: total smooth (i.e., linear), 1: no smoothing
%The interesting range of p is often near 1/(1+h^3/6) where h is the average spacing of the data points

plotTest                    = logical(0);

% defaults
BENDING.failed              = false;
BENDING.subPlateIndices     = zeros(size(GRID.X_3D));
BENDING.subPlateFullFit     = [NaN,NaN];
BENDING.subPlateFit         = [NaN,NaN];
BENDING.circleCenter        = [NaN,NaN];	
BENDING.minimumPos          = [NaN,NaN];	
BENDING.radius              = NaN;	
BENDING.radiusSI            = NaN;	
BENDING.Tcontour            = NaN;	
BENDING.TcontourSlab        = NaN;
BENDING.slabPolarity        = 'na';
idxSealevel                 = find(abs(GRID.Z_3Dp(1,1,:))==min(abs(GRID.Z_3Dp(1,1,:)))); %find index close to sealevel

%% problem checks
if strcmp(GRID.Dim,'3-D')
    return
    
end

%% get temperature contour of plate and slab
% find values below, at & above the contour value
BW0 = zeros(size(T_3D,1),size(T_3D,3));
BW0(T_3D(:,1,:)<=TempThreshold & T_3D(:,1,:)>min(T_3D(:)))	= 1; %below or at

%% point cloud representing subducting plate
%thinning of data
BW              = bwmorph(BW0,'thin',2);

%delete extremely deeper than threshold
maxDepthThreshold               = 3.0* PLATE.LABdepth;  %~450km, limit depth criterion (needs to be deeper than the following depth criterion!) <-----------------
if ~isnan(PLATE.SlabTipPosition(1,1)) && ...
        3/4*PLATE.SlabTipPosition(1,2)<maxDepthThreshold
    maxDepthThreshold           = 3/4*PLATE.SlabTipPosition(1,2); %further limiting to slab depth (as tip is sometimes bent strongly)
end
BW(GRID.z2dp>maxDepthThreshold) = 0; 

%keep only lowermost values
numPoints2keep  = 2;  %should be >1 to not produce gaps or bumps
for ix=1:size(BW,1)
    for iz=2:size(BW,2)
        if BW(ix,iz)==1
            BW(ix,min(iz+numPoints2keep,size(BW,2)):end)    = 0; %delete shallow hits
            continue
        end
    end
end

%make sure points close to each other are connected
BW              = bwmorph(BW,'thicken',2);
BW              = bwmorph(BW,'bridge');
%plot(GRID.x2dp(BW),GRID.z2dp(BW),'.')

%delete points deeper than threshold
maxDepthThreshold               = 2.2* PLATE.LABdepth;  %~370km, limit depth criterion <-----------------
if ~isnan(PLATE.SlabTipPosition(1,1)) && ...
        3/4*PLATE.SlabTipPosition(1,2)<maxDepthThreshold
    maxDepthThreshold           = 3/4*PLATE.SlabTipPosition(1,2); %further limiting to slab depth (as tip is sometimes bent strongly)
end
BW(GRID.z2dp>maxDepthThreshold) = 0;

%delete points shallower than threshold
%BW(GRID.z2dp<1e3*GRID.m2p)      = 0;       %~1km, limit shallowness criterion (to prevent max curvature at spreading ridge) <-----------------
minDepthThreshold               = PLATE.LABdepth/150;   %~1km, limit shallowness criterion (to prevent max curvature at spreading ridge) <-----------------
BW(GRID.z2dp<minDepthThreshold)	= 0;

%find connected points with deepest reach
CC = bwconncomp(BW);
maxZ = zeros(size(CC.PixelIdxList,2),1);
minZ = maxZ;
slabFound = maxZ;
for i=1:size(CC.PixelIdxList,2)
    %skip isolated small areas
    if size(CC.PixelIdxList{1,i},1)<10 %specify small area criterion <-----
        continue
    end
    maxZ(i,1)       = max(GRID.z2dp(CC.PixelIdxList{1,i}));
    minZ(i,1)       = min(GRID.z2dp(CC.PixelIdxList{1,i}));
    %skip areas separated from surface
    depthTresholdSurface        = 120e3*GRID.m2p;   %COULD BE IMPROVED WITH AUTOMATISED PLATE THICKNESS VALUE
    if minZ(i,1)>depthTresholdSurface %connection to surface criterion <----------
        continue
    end
    slabFound(i,1) = 1;
end
maxZ(~slabFound) = 0;
% minZ(~slabFound) = 0;
if max(maxZ)==0
    BENDING.failed = true;
    if SWITCH.Verbose; warning('Plate geometry could not be found! Check here.'); end
    return
end
[~,idx] = max(maxZ);
%remove all others
for i=1:size(CC.PixelIdxList,2)
    if i~=idx
        BW(CC.PixelIdxList{1,i}) = 0;
    end
end

%% check subduction/plate polarity
if strcmp(GRID.Type,'Cartesian')
    x               = GRID.x2dp(BW);
elseif strcmp(GRID.Type,'spherical2D')
    clearvars dummy
    x               = GRID.X_3Dsp(BW);
    xrad            = GRID.X_3D(BW);
end
z                   = GRID.z2dp(BW);
if isempty(z)
    if SWITCH.Verbose; warning('Bending diagnostics failed.'); end
    BENDING.failed 	= true;
    return
end
if z(1,1)>z(end,1) && x(1,1)>x(end,1) || z(1,1)<z(end,1) && x(1,1)<x(end,1) %subduction on right hand side
    sPolarity       = 'R';
elseif z(1,1)>z(end,1) && x(1,1)<x(end,1) || z(1,1)<z(end,1) && x(1,1)>x(end,1) %subduction on left hand side
    sPolarity       = 'L';
else
    sPolarity       = 'na'; if SWITCH.Verbose; warning('Plate polarity could not be found! Check here.'); end
end

%% horizontal width limit
xWidthTreshold      = maxDepthThreshold *4; %width limit threshold <-------------------------
if strcmp(GRID.Type,'spherical2D') %convert limits to spherical coordinates
%     subdAngleRad    = mean(xrad-pi); %angle of subduction zone from top of spherical annulus
    xWidthTreshold 	= xWidthTreshold*2.5; %width limit threshold <-----------------------
end

% %% find subducting plate thickness and use for depth threshold to remove points
% % numXpoints_crit2              = max(1,round(dxCrit2/dx));
% % numXpointsLeft2               = min(ix-1,numXpoints_crit2);
% % PlateThicknessL(ix,iy)    	= lithoThicknessT_1600(ix-numXpointsLeft2,iy);
% dummy = x;
% if strcmp(GRID.Type,'spherical2D')
%     dummy = xrad;  %THIS NEEDS TO BE ADJUSTED FOR CYLINDRICAL GEOMETRY!!
% end
% if strcmp(sPolarity,'na')
%     subPlateThickness   = PLATE.LABdepth;
% elseif strcmp(sPolarity,'L')
%     [~,index] = min(abs( GRID.x2dp(:,idxSealevel)-(min(dummy(:))+xWidthTreshold) ));
%     subPlateThickness   = LITHO.thicknessFull(index,1);
%     
% elseif strcmp(sPolarity,'R')
%     [~,index] = min(abs( GRID.x2dp(:,idxSealevel)-(max(dummy(:))-xWidthTreshold) ));
%     subPlateThickness   = LITHO.thicknessFull(index,1);
%     
% end
% maxDepthThreshold2               = 3.1* subPlateThickness;  %~370km, limit depth criterion <-----------------
% 
% %delete points deeper than threshold
% if ~isnan(PLATE.SlabTipPosition(1,1)) && ...
%         3/4*PLATE.SlabTipPosition(1,2)<maxDepthThreshold2
%     maxDepthThreshold2           = 3/4*PLATE.SlabTipPosition(1,2); %further limiting to slab depth (as tip is sometimes bent strongly)
% end
% BW(GRID.z2dp>maxDepthThreshold2) = 0;
% dummy = z;
% x(dummy>maxDepthThreshold2)  = [];
% z(dummy>maxDepthThreshold2)  = [];
% if strcmp(GRID.Type,'spherical2D'); xrad(dummy>maxDepthThreshold2) = []; end

%% keep subducting-plate indices
LowerPlateIndices 	= bwmorph(BW,'fill',inf);
LowerPlateIndices 	= bwmorph(LowerPlateIndices,'thin',inf);

%% keep subducting-plate spline
LowerPlateFullX     = x;
LowerPlateFullZ     = z;
if strcmp(GRID.Type,'spherical2D'); xradFull = xrad; end

%% remove surface values far away from the trench (if polarity is known)
dummy = x;
if strcmp(sPolarity,'L')
    x(dummy>min(dummy(:))+xWidthTreshold)  = [];
    z(dummy>min(dummy(:))+xWidthTreshold)  = [];
    if strcmp(GRID.Type,'spherical2D'); xrad(dummy>min(dummy(:))+xWidthTreshold) = []; end
elseif strcmp(sPolarity,'R')
    x(dummy<max(dummy(:))-xWidthTreshold)  = [];
    z(dummy<max(dummy(:))-xWidthTreshold)  = [];
    if strcmp(GRID.Type,'spherical2D'); xrad(dummy<max(dummy(:))-xWidthTreshold) = []; end
end

%% remove angled endge-points to prevent artificial curvature of the spine
for ip=1:3  %remove three horizontally-outermost point rows
    dummy = x;
    if strcmp(sPolarity,'L')
        x(dummy==min(dummy(:)))  = [];
        z(dummy==min(dummy(:)))  = [];
        if strcmp(GRID.Type,'spherical2D'); xrad(dummy==min(dummy(:))) = []; end
    elseif strcmp(sPolarity,'R')
        x(dummy==max(dummy(:)))  = [];
        z(dummy==max(dummy(:)))  = [];
        if strcmp(GRID.Type,'spherical2D'); xrad(dummy==max(dummy(:))) = []; end
    end
end

%% adjust for spherical2D geometry
if strcmp(GRID.Type,'spherical2D')
    %convert to spherical coordinates FULL PLATE
    R             	= GRID.rcmb_p +max(GRID.Z_3Dp(1,1,:)) - LowerPlateFullZ;
    xFitS       	= R.*-sin(xradFull);    %.*cos(y2d_nd);
    zFitS        	= R.*-cos(xradFull);
    LowerPlateFullX	= xFitS;
    LowerPlateFullZ	= zFitS;
    %convert to spherical coordinates SUBDUCTION PART
    R             	= GRID.rcmb_p +max(GRID.Z_3Dp(1,1,:)) - z;
%     xFitDeg       	= x./R;               %[degrees]
    xFitS       	= R.*-sin(xrad);    %.*cos(y2d_nd);
    zFitS        	= R.*-cos(xrad);
    x            	= xFitS;
    z           	= zFitS;
    clearvars R xFitDeg
end

%% add spline to the slab curvature:
%add a fit (smoothingspline) to point cloud
try
    fo             	= fit(LowerPlateFullX,LowerPlateFullZ,'smoothingspline','smoothingparam',smoothingParam);
catch me
    BENDING.failed	= true;
    warning(me.message);
    return
end
numPointsFit        = 1000;
LowerPlateFitX    	= min(LowerPlateFullX):(max(LowerPlateFullX)-min(LowerPlateFullX))/numPointsFit:max(LowerPlateFullX);
LowerPlateFitZ     	= fo(LowerPlateFitX)';
%add a fit (smoothingspline) to point cloud
try
    fo              = fit(x,z,'smoothingspline','smoothingparam',smoothingParam);
catch me
    BENDING.failed 	= true;
    warning(me.message);
    return
end
numPointsFit        = 1000;
xFit                = min(x):(max(x)-min(x))/numPointsFit:max(x);
zFit                = fo(xFit)';

%first derivative
dzdx                = zeros(size(xFit));
for ix=2:size(xFit,2)-1
    %dx            	= abs( (xFit(1,ix+1)-xFit(1,ix-1))/2 );
    dx            	= (xFit(1,ix+1)-xFit(1,ix-1))/2;
    dzdx(1,ix)      = ( zFit(1,ix+1)-zFit(1,ix-1) )/(2*dx);
end
dzdx(1,1)           = dzdx(1,2);        %boundary condition
dzdx(1,end)     	= dzdx(1,end-1);    %boundary condition

%second derivative
d2zdx2 = zeros(size(xFit)); %d2zdx2x = d2zdx2;
for ix=2:size(xFit,2)-1
    dx            	= abs( (xFit(1,ix+1)-xFit(1,ix-1))/2 );
    d2zdx2(1,ix)    = ( zFit(1,ix+1)+zFit(1,ix-1)-2*zFit(1,ix) )/(dx^2);
end
d2zdx2(1,1)         = d2zdx2(1,2);        %boundary condition
d2zdx2(1,end)       = d2zdx2(1,end-1);    %boundary condition

%% ERROR CHECK
if isnan(d2zdx2(1,1)) && isnan(max(d2zdx2(:))) && isnan(min(d2zdx2(:))) %only NaNs
    warning('bending failed')
    BENDING.failed  = true;
    return
    
end
if d2zdx2(1,1)==0 && max(d2zdx2(:))==0 && min(d2zdx2(:))==0 %only zeros
    warning('bending failed')
    BENDING.failed  = true;
    return
    
end
%prevent zero values (there is a division later on)
d2zdx2(d2zdx2==0)   = NaN;

%% MINIMUM RADIUS OF CURVATURE
% Rcurvature = ( 1+(dz/dx)^2 )^(3/2)/( abs(d2z/dx^2) );
Rcurvature          = ((1+dzdx.^2).^(3/2))./(abs(d2zdx2));
RcurvatureEff     	= ((1+dzdx.^2).^(3/2))./(d2zdx2); %can be negative
minBendingRadius    = min(abs(RcurvatureEff));
minBendingRadiusEff = RcurvatureEff(Rcurvature==minBendingRadius);

xCircleContact      = xFit(Rcurvature==minBendingRadius);
zCircleContact      = zFit(Rcurvature==minBendingRadius);
dzdxContact         = dzdx(Rcurvature==minBendingRadius);

% dipAngleContact     = abs( atand(dzdxContact) );  %degrees - Winkel zur Horizontalen - always absolute
% dxRadius            = abs( sind(dipAngleContact)*minBendingRadius );
% dzRadius            = abs( cosd(dipAngleContact)*minBendingRadius );
dipAngleContact     = atand(dzdxContact);  %degrees - Winkel zur Horizontalen des Koordinatensystems
dxRadius            = -sind(dipAngleContact)*minBendingRadiusEff;
dzRadius            = cosd(dipAngleContact)*minBendingRadiusEff;
if strcmp(GRID.Type,'spherical2D')
    %angleContact    = atan(zCircleContact/xCircleContact);  %angle to the right-horizontal
    %dzRadius        = -dzRadius;   %Grid is not flipped
end
adjustSign     = 1;
adjustSignZ    = 1;
% if strcmp(sPolarity,'R')
%     adjustSign = -1;
% else
%     adjustSign = 1;
% end
% if strcmp(GRID.Type,'spherical2D')
%     if strcmp(sPolarity,'R')
%         if angleContact>pi %top half of spherical annulus
%             adjustSign      = -1;
%             adjustSignZ     = -1;       %grid is not flipped vertically
%         else
%             adjustSign      = 1;
%             adjustSignZ     = 1;       %grid is not flipped vertically
%         end
%     else
%         if angleContact>pi %top half of spherical annulus
%             adjustSign      = 1;
%             adjustSignZ     = -1;       %grid is not flipped vertically
%         else
%             adjustSign      = -1;
%             adjustSignZ     = 1;       %grid is not flipped vertically
%         end
%     end
% else
%     adjustSignZ     = 1;
% end
xCircle             = xCircleContact+dxRadius*adjustSign;
zCircle             = zCircleContact+dzRadius*adjustSignZ;

if plotTest
    figure(2),clf
    hold on
    plot(LowerPlateFullX,LowerPlateFullZ)
    plot(LowerPlateFitX,LowerPlateFitZ,'-r','LineWidth',2)
    axis equal
    axis ij
    
    figure(2),clf
    hold on
    plot(fo,x,z)
    plot(xFit,zFit,'-r','LineWidth',2)
    plot(xCircleContact,zCircleContact,'o')
%     f_plotCircle(xCircle,zCircle,minBendingRadius,[0.5 0.5 0.5],1)
    axis equal
    axis ij
        
%     figure(2),clf
%     hold on
%     plot(fo,x,z)
%     plot(xFitS',zFitS','-r','LineWidth',2)
%     plot(xCircleContactS,zCircleContactS,'o')
% %     f_plotCircle(xCircle,zCircle,minBendingRadius,[0.5 0.5 0.5],1)
%     axis equal
%     axis ij
end

%% function output
BENDING.subPlateIndices = LowerPlateIndices;                        %1 and 0s
BENDING.subPlateFullFit = [LowerPlateFitX',LowerPlateFitZ'];        %[plotting dimension] [x,z]
BENDING.subPlateFit     = [xFit',zFit'];                            %[plotting dimension] [x,z]
BENDING.circleCenter    = [xCircle,zCircle];                        %[plotting dimension] [x,z]
BENDING.minimumPos      = [xCircleContact,zCircleContact];          %[plotting dimension] [x,z]
BENDING.radius          = minBendingRadius;                         %[plotting dimension]
BENDING.radiusSI        = BENDING.radius /GRID.m2p;                 %[m] SI-units
BENDING.Tcontour        = BW0;                                      %[plotting dimension] 1 and 0  %thinned T-contour area
BENDING.TcontourSlab    = BW;                                       %[plotting dimension] 1 and 0  %thinned T-contour area
BENDING.slabPolarity    = sPolarity;                                %'L','R', or 'na'
% BENDING.failed                                                      %true or false
end



%%                     SLAB VISCOSITY & TEMPERATURE & DENSITY & STRESS 1.2
%                                                Fabio Crameri, 23.09.2018
function [SLAB] = f_slabDiagnostics(SLAB,BENDING,SWITCH,GRID,T_3D,ETA_3D,RHO_3D,STR_3D,BASALT_3D)
%% defaults
SLAB.DepthChecked       = [NaN, NaN];
SLAB.viscosity          = NaN;
SLAB.viscosityMin       = NaN;
SLAB.temperatureMin     = NaN;
SLAB.densityMax         = NaN;
SLAB.stressMax          = NaN;
SLAB.mantleTemperature  = NaN;
SLAB.mantleViscosity   	= NaN;
SLAB.mantleDensity   	= NaN;
SLAB.mantleViscContrast	= NaN;
SLAB.mantleTempContrast = NaN;
SLAB.Area2check         = NaN;

%% problem checks
if strcmp(GRID.Dim,'3-D') || isnan(BENDING.TcontourSlab(1))
    return
    
end
slabViscRoutine = 1;
if slabViscRoutine==1    
    idxSealevel  	= find(abs(GRID.Z_3Dp(1,1,:))==min(abs(GRID.Z_3Dp(1,1,:)))); %find index close to sealevel
    if strcmp(GRID.Type,'spherical2D')
        bendingMinPosR 	= sqrt(BENDING.minimumPos(1,1)^2+BENDING.minimumPos(1,2)^2);
        bendingMinPosZ 	= max(GRID.Rp(:)) - bendingMinPosR;
        bendingMinPosX  = (asin(-BENDING.minimumPos(1,1)/bendingMinPosR)+2*pi)*bendingMinPosR;  %....... shift in 2*pi needed? ......
    else
        bendingMinPosX  = BENDING.minimumPos(1,1);
        bendingMinPosZ 	= BENDING.minimumPos(1,2);
    end
    distance2check  = 9/10 *bendingMinPosZ;      %[plotting dimension]  %distance min. bending to surface <-----------
    
%     if ~SWITCH.DimensionalMode
%         distance2check = distance2check/GRID.dimFactor;
%     end
    SLAB.Area2check = BENDING.TcontourSlab;
    if strcmp(GRID.Type,'spherical2D')
        dummyR      = GRID.Rp; %radius for conversion to length (L = r*alpha)
    else
        dummyR      = ones(size(GRID.Rp)); %no adjustment needed
    end
    
%     [~,dummyInd] = min(abs(GRID.z2dp(1,:)-bendingMinPosZ(1,1))); %find z-index of min bending position
    for ix=1:size(SLAB.Area2check,1)
        if GRID.x2dp(ix,idxSealevel)*dummyR(ix,idxSealevel)>bendingMinPosX+2*distance2check || ...
                GRID.x2dp(ix,idxSealevel)*dummyR(ix,idxSealevel)<bendingMinPosX-2*distance2check
%         if GRID.x2dp(ix,dummyInd)*dummyR(ix,dummyInd)>bendingMinPosX+2*distance2check || ...
%               GRID.x2dp(ix,dummyInd)*dummyR(ix,dummyInd)<bendingMinPosX-2*distance2check
            SLAB.Area2check(ix,:) = 0; %remove contour points outside check area in x-direction
        end
    end
    clearvars dummyR
    for iz=1:size(SLAB.Area2check,2)
        if GRID.z2dp(1,iz)>bendingMinPosZ+distance2check || ...
                GRID.z2dp(1,iz)<bendingMinPosZ-distance2check
            SLAB.Area2check(:,iz) = 0; %remove contour points outside check area in z-direction
        end
    end
    %adjust area
    SLAB.Area2check     = bwmorph(SLAB.Area2check,'thicken');
    SLAB.Area2check     = bwmorph(SLAB.Area2check,'bridge');
    SLAB.Area2check     = bwmorph(SLAB.Area2check,'fill');
    %remove parts with weak crust
    if max(SLAB.Area2check(BASALT_3D>0.9)>0)
        if SWITCH.Verbose; warning('Some weak-crust regions are removed from the slab-viscosity checking area.'); end
    end
    SLAB.Area2check(BASALT_3D>0.9) = 0;
    
    if logical(0) %test resulting slab area
       figure(2)
       plot(GRID.x2dp(SLAB.Area2check),GRID.z2dp(SLAB.Area2check),'.')
    end
    SLAB.DepthChecked   = [min(GRID.z2dp(SLAB.Area2check)), max(GRID.z2dp(SLAB.Area2check))]; %min&max
    SLAB.temperatureMin	= min(min(T_3D(SLAB.Area2check)));          %shallow slab check point (1 gridpoint in each direction added)
    SLAB.viscosityMin  	= min(min(ETA_3D(SLAB.Area2check)));        %shallow slab check point (1 gridpoint in each direction added)
    SLAB.viscosity     	= mean2(ETA_3D(SLAB.Area2check));           %shallow slab check point (1 gridpoint in each direction added)
    SLAB.densityMax     = max(max(RHO_3D(SLAB.Area2check)));        %shallow slab check point (1 gridpoint in each direction added)
    SLAB.stressMax      = max(max(STR_3D(SLAB.Area2check)));        %shallow slab check point (1 gridpoint in each direction added)
    if isempty(SLAB.viscosity)
        warning('Slab viscosity could not be found: Old method is applied instead.')
        slabViscRoutine = 2; %Slab viscosity not found: apply other method
    end
end

if slabViscRoutine==2
    dummy               = bwmorph(SLAB.coldSlabs,'thicken');        %x-range to check
    dummy2              = max(SLAB.zIdxCheck-1,0):min(SLAB.zIdxCheck+1,size(ETA_3D,3)); %z-range to check
    SLAB.DepthChecked   = [min(min(GRID.z2dp(dummy,dummy2))), max(max(GRID.z2dp(dummy,dummy2)))];
    SLAB.temperatureMin	= min(min(T_3D(dummy,dummy2)));             %shallow slab check point (1 gridpoint in each direction added)
    SLAB.viscosityMin  	= min(min(ETA_3D(dummy,dummy2)));           %shallow slab check point (1 gridpoint in each direction added)
    SLAB.viscosity     	= mean2(ETA_3D(dummy,dummy2));              %shallow slab check point (1 gridpoint in each direction added)
    % etaSlab         = ETA_3D(coldSlabs,zIdxCheck);                %takes just one point into account
    SLAB.densityMax     = max(max(RHO_3D(dummy,dummy2)));           %shallow slab check point (1 gridpoint in each direction added)
    SLAB.stressMax      = max(max(STR_3D(dummy,dummy2)));           %shallow slab check point (1 gridpoint in each direction added)
end

%% PREVENT EMPTY ARRAYS
if isempty(SLAB.DepthChecked);  SLAB.DepthChecked   = [NaN, NaN]; end
if isempty(SLAB.viscosityMin);  SLAB.viscosityMin   = NaN; end
if isempty(SLAB.temperatureMin);SLAB.temperatureMin	= NaN; end
if isempty(SLAB.densityMax);    SLAB.densityMax     = NaN; end
if isempty(SLAB.stressMax);     SLAB.stressMax      = NaN; end

%% MANTLE TEMPERATURE & VISCOSITY
SLAB.mantleTemperature  = median(T_3D(:,:,SLAB.zIdxCheckUM));            %at deep slab check point level
SLAB.mantleViscosity   	= median(ETA_3D(:,:,SLAB.zIdxCheckUM));          %at deep slab check point level
SLAB.mantleDensity   	= median(RHO_3D(:,:,SLAB.zIdxCheckUM));          %at deep slab check point level

%% SLAB-MANTLE VISCOSITY- AND TEMPERATURE CONTRAST
SLAB.mantleTempContrast  = SLAB.temperatureMin./SLAB.mantleTemperature;	%slab-mantle difference
if strcmp(SLAB.EtaToTake,'min')
    SLAB.mantleViscContrast  = SLAB.viscosityMin./SLAB.mantleViscosity;	%slab-mantle difference
elseif strcmp(SLAB.EtaToTake,'mean')
    SLAB.mantleViscContrast  = SLAB.viscosity./SLAB.mantleViscosity;   	%slab-mantle difference
end
end



%%                               SUBDUCTION-TOPOGRAPHY CHARACTERISTICS 3.03
%                                                Fabio Crameri, 11.12.2018
function [SUBTOPO] = f_SubTopoCharacteristics(PLATE,GRID,SWITCH,FILE,STYLE)
%defaults
SUBTOPO.Successful          = true;
SUBTOPO.overridingPlate     = NaN;
SUBTOPO.subductingPlate     = NaN;
SUBTOPO.TrenchZ          	= NaN;
SUBTOPO.TrenchX             = NaN;
SUBTOPO.BackArcBasinX       = NaN;
SUBTOPO.BackArcBasinZ       = NaN;
SUBTOPO.BackArcBasinExtentX	= NaN;
SUBTOPO.BackArcBasinExtentZ = NaN;
SUBTOPO.InundationX      	= NaN;
SUBTOPO.InundationZ      	= NaN;
SUBTOPO.InundationDistance 	= NaN;
SUBTOPO.BackArcBasinArea   	= NaN;
SUBTOPO.IslandArcX          = NaN;
SUBTOPO.IslandArcZ          = NaN;
SUBTOPO.ForeBulgeX        	= NaN;
SUBTOPO.ForeBulgeZ        	= NaN;
SUBTOPO.UPtiltAngleXnear	= NaN;
SUBTOPO.UPtiltAngleXfar     = NaN;
SUBTOPO.UPtiltAngleZnear	= NaN;
SUBTOPO.UPtiltAngleZfar     = NaN;
SUBTOPO.UPtiltAngle         = NaN;

%% problem checks
if strcmp(GRID.Dim,'3-D')
    if SWITCH.Verbose; warning('Subduction-Topography characteristics not implemented for 3-D!'); end
    SUBTOPO.Successful      = false;
    return
    
end
%check if curve fitting toolbox is installed
try
    smooth(1:5); %it is present
catch me %it is not present
    warning(me.message)
    SUBTOPO.Successful      = false;
    return
    
end

%read topo data
[~,~,~,~,TOPO_3D,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Topography',SWITCH);
if ischar(TOPO_3D)
    if SWITCH.Verbose; warning(TOPO_3D); end
    SUBTOPO.Successful      = false;
    return
end
if size(TOPO_3D,1)==1 %exchange x and y
    dummy_3D = zeros(size(TOPO_3D,2),size(TOPO_3D,1),size(TOPO_3D,3));
    dummy_3D(:,1,:) = TOPO_3D(1,:,:);  	TOPO_3D = dummy_3D;
end
if strcmp(GRID.Dim,'2-D')
    topo(:,1)               = TOPO_3D(:,1,end);  %takes TOPO_3D(:,1,2) for e.g., StagYY
elseif strcmp(GRID.Dim,'3-D')
    if SWITCH.Verbose; warning('Subduction-Topography characteristics not implemented for 3-D!'); end
    SUBTOPO.Successful      = false;
    return
    %topo(:,:)   = TOPO_3D(:,:,2);
end
%     %actually not needed - just to be sure:  %might only be correct for 2-D?
%     topo(:,1) = topo(:,1)-mean(topo(:,1));    % remove mean value (CMB topo)
%     topo(:,2) = topo(:,2)-mean(topo(:,2));    % remove mean value (surf topo)

%adjust horizontal grid variables
%find index close to sealevel
idxSealevel = find(abs(GRID.Z_3Dp(1,1,:))==min(abs(GRID.Z_3Dp(1,1,:))));
if strcmp(GRID.Type,'Cartesian')
    xp(:,1)                 = GRID.X_3Dp(:,1,idxSealevel);
    yp(:,1)                 = GRID.Y_3Dp(1,:,idxSealevel);
elseif strcmp(GRID.Type,'spherical2D')
    xp(:,1)                 = GRID.X_3Dsp(:,1,idxSealevel);
    %yp(:,1)                 = GRID.Y_3Dsp(1,:,idxSealevel);
end
dxp                         = xp(2,1)-xp(1,1);

%dimensionalise topography (according to grid dimensions)
topo                        = topo.*GRID.dimFactor;

%smooth topography for topo characteristics
topoSmooth                  = topo;
for jj=1:2
    dummy(:,1)              = smooth(topoSmooth(:,1),5,'moving');
    topoSmooth              = dummy;
end
%smooth topography for upper-plate tilt
topoSmoothest            	= topo;
for jj=1:2
    dummy(:,1)              = smooth(topoSmoothest(:,1),20,'moving');
    topoSmoothest         	= dummy;
end

%test figure
if logical(0)
    figure(2)
    plot(topo); hold on; plot(topoSmooth,'r'); plot(topoSmoothest,'g')
    figure(1)
end

%check if subduction polarity is found
if isempty(PLATE.Subduction) || min(isnan(PLATE.Subduction(:)))==1 || strcmp(PLATE.SubPolarity{1,1},'na') || strcmp(PLATE.SubPolarity{1,1},'')
    platePolarityFound      = false;
else
    platePolarityFound      = true;
end

if platePolarityFound
    try %topography characteristics routine is specialised for high-resulution cases
        
    %find trench depth
    shiftXawayFromTrench        = 100e3;        %[m], +/- distance to check for trench minimum <<<<<<<<<<<<<<<<<<
    shiftXawayFromTrench        = shiftXawayFromTrench*GRID.m2p;     %convert to [plotting]
    shiftXawayFromTrenchIdx     = ceil(shiftXawayFromTrench/(mean(GRID.dxp(:,:,idxSealevel))))+1; %convert to index
    lowerIndex                  = PLATE.idxTrench-shiftXawayFromTrenchIdx;
    upperIndex                  = PLATE.idxTrench+shiftXawayFromTrenchIdx;
    topo2check                  = topoSmooth(max(1,lowerIndex):min(size(topoSmooth,1),upperIndex));
    [TrenchDepth,TrenchIndex] 	= min(topo2check); %in grid dimensions
    TrenchIndex                 = TrenchIndex+lowerIndex-1;
    TrenchX                     = xp(TrenchIndex,1);
    
    %add x values to topoSmooth
    topoSmooth                  = [xp, topoSmooth];
    
    %find overriding and subducting plate by tracked trench location
    if strcmp(PLATE.SubPolarity,'right')
        overridingPlateTopo   	= topoSmooth(xp(:,1)>=TrenchX(1,1),:);
        subductingPlateTopo    	= topoSmooth(xp(:,1)<=TrenchX(1,1),:);
    elseif strcmp(PLATE.SubPolarity,'left')
        overridingPlateTopo   	= topoSmooth(xp(:,1)<=TrenchX(1,1),:);
        subductingPlateTopo   	= topoSmooth(xp(:,1)>=TrenchX(1,1),:);
    end

    %find reference mean elevation of the overriding plate (at left hand side)
    xDistance4mean              = 2000e3; %m distance from trench to check for mean UP height <<<<<<<<<<<<<<<<<<<<<<<
    xWidth4mean                 = 400e3; %m distance over which mean is taken <<<<<<<<<<<<<<<<<<<<<<<
    xDistance4mean           	= xDistance4mean*GRID.m2p; %conversion to plotting values
    xWidth4mean                 = xWidth4mean*GRID.m2p; %conversion to plotting values
    
    UpperPlateDistanceFromTrench= abs( overridingPlateTopo(:,1)-TrenchX(1,1) );
    xDistance4mean              = min(xDistance4mean,max(UpperPlateDistanceFromTrench)); %adjust x distance if side boundary close by
    if xDistance4mean<200e3*GRID.m2p; warning('Upper-plate topography diagnostics are badly conditioned! Check here!'); end
    
    dummy                       = UpperPlateDistanceFromTrench(:,1)>xDistance4mean-xWidth4mean/2 & UpperPlateDistanceFromTrench(:,1)<xDistance4mean+xWidth4mean/2;
    charUPTopo                  = overridingPlateTopo(dummy,2);
    dummy                       = abs(UpperPlateDistanceFromTrench(:,1)-xDistance4mean);
    [~,idx]                     = min(dummy); %index of closest value
    charUPTopoPoint             = overridingPlateTopo(idx,:); %closest value
    meanOPelevation             = mean(charUPTopo);
    varOPelevation              = max(charUPTopo) - min(charUPTopo);
    criticalUPz                 = meanOPelevation-2.0*varOPelevation; %adjust critical z-value to define inundation's max. horizontal extent
    criticalVAR                 = 2.0*varOPelevation;
    
    %find back-arc basin minimum
    BackArcBasinZ               = overridingPlateTopo(overridingPlateTopo(:,2)==min(overridingPlateTopo(:,2)),:); %min. of overriding plate
    if BackArcBasinZ(1,1)==TrenchX(1,1) %inundation depression not found
        overridingPlateTopo2    = overridingPlateTopo;
        if strcmp(PLATE.SubPolarity,'right')
            while BackArcBasinZ(1,1)==overridingPlateTopo2(1,1)
                overridingPlateTopo2(1,:) = [];
                BackArcBasinZ     	= overridingPlateTopo2(overridingPlateTopo2(:,2)==min(overridingPlateTopo2(:,2)),:);
            end
        elseif strcmp(PLATE.SubPolarity,'left')
            while BackArcBasinZ(1,1)==overridingPlateTopo2(end,1)
                overridingPlateTopo2(end,:) = [];
                BackArcBasinZ     	= overridingPlateTopo2(overridingPlateTopo2(:,2)==min(overridingPlateTopo2(:,2)),:);
            end
        end
    end
    clearvars overridingPlateTopo2
    
    %CHECK FOR PROBLEMS
    %if varOPelevation>0.5 or so, than no back-arc basin extent can be found
    %(due to no good reference to upper plate reference height)
    maxUPdifferenceAllowed      = 700;  %[m] <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    maxUPdifferenceAllowed      = maxUPdifferenceAllowed *GRID.m2p;
    if varOPelevation>maxUPdifferenceAllowed
        preventBABxMeasurement  = true;
    else
        preventBABxMeasurement  = false;
    end
    %no back-arc depression found
    if (strcmp(PLATE.SubPolarity,'right') && BackArcBasinZ(1,1)<=TrenchX(1,1)) || ... %inundation on the subducting plate
            (strcmp(PLATE.SubPolarity,'left') && BackArcBasinZ(1,1)>=TrenchX(1,1))
        BackArcBasinZ           = [NaN NaN]; %not found
        if SWITCH.Verbose; warning('>>>> back-arc depression not found!'); end
    end
    if BackArcBasinZ(1,2)>criticalUPz %it is too shallow
        BackArcBasinZ           = [NaN NaN]; %not found
        if SWITCH.Verbose; warning off backtrace; warning('>>>> back-arc depression not found: too shallow!'); warning on backtrace; end
    end
    
    %find topo to the left of the back-arc basin
    if ~isnan(BackArcBasinZ(1,1))
        if strcmp(PLATE.SubPolarity,'right')
            topoAwayOfBAB       = overridingPlateTopo(overridingPlateTopo(:,1)>BackArcBasinZ(1,1),:);
        elseif strcmp(PLATE.SubPolarity,'left')
            topoAwayOfBAB    	= overridingPlateTopo(overridingPlateTopo(:,1)<BackArcBasinZ(1,1),:);
        end
    else %no back-arc depression
        topoAwayOfBAB           = overridingPlateTopo;
    end
    
    %find the maximum BAB-extent
    checkRange                  = 5;    %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    BackArcBasinExtent        	= [NaN NaN];
    if ~preventBABxMeasurement
        if strcmp(PLATE.SubPolarity,'right')
            while isnan(BackArcBasinExtent(1,1)) || BackArcBasinExtent(1,2)<criticalUPz-criticalVAR %still too deep
                icheck                  = 1;
                if topoAwayOfBAB(icheck,2)>criticalUPz
                    if SWITCH.Verbose; warning('>>>> back-arc depression extent not found!'); end
                    break
                end
                while topoAwayOfBAB(icheck,2)<criticalUPz && ... %check critical depth - from left to right
                        mean(topoAwayOfBAB(icheck+1:icheck+checkRange,2))>topoAwayOfBAB(icheck,2) %it's still getting shallower to the right
                    BackArcBasinExtent 	= topoAwayOfBAB(icheck,:);
                    icheck              = icheck+1;
                end
                icheck                  = 1; %reset
                checkRange              = checkRange+5; %increase checking range
            end
        elseif strcmp(PLATE.SubPolarity,'left')
            while isnan(BackArcBasinExtent(1,1)) || BackArcBasinExtent(1,2)<criticalUPz-criticalVAR %still too deep
                icheck                  = size(topoAwayOfBAB,1);
                if topoAwayOfBAB(icheck,2)>criticalUPz
                    if SWITCH.Verbose; warning('>>>> back-arc depression extent not found!'); end
                    break
                end
                while topoAwayOfBAB(icheck,2)<criticalUPz && ... %check critical depth - from right to left
                        mean(topoAwayOfBAB(icheck-checkRange:icheck-1,2))>topoAwayOfBAB(icheck,2) %it's still getting shallower to the left
                    BackArcBasinExtent 	= topoAwayOfBAB(icheck,:);
                    icheck              = icheck-1;
                end
                icheck                  = size(topoAwayOfBAB,1); %reset
                checkRange              = checkRange+5; %increase checking range
            end
        end
    end
    
    %find upper-plate inundation
    Inundation                  = [NaN, NaN];
    if max(topoAwayOfBAB(:,2))<=0 || min(topoAwayOfBAB(:,2))>0 %whole upper plate under- or above water level
        preventInundationMeasurement = true;
    else
        preventInundationMeasurement = false;
    end
    if ~preventInundationMeasurement
        if strcmp(PLATE.SubPolarity,'right')
            icheck                  = 1;
            while topoAwayOfBAB(icheck,2)<0 %check critical depth - from left to right
                Inundation          = topoAwayOfBAB(icheck,:);
                icheck              = icheck+1;
            end
        elseif strcmp(PLATE.SubPolarity,'left')
            icheck                  = size(topoAwayOfBAB,1);
            while topoAwayOfBAB(icheck,2)<0 %check critical depth - from right to left
                Inundation          = topoAwayOfBAB(icheck,:);
                icheck              = icheck-1;
            end
        end
    end
    
    %find segment between trench and BABZmin
    if strcmp(PLATE.SubPolarity,'right')
        topoBasinTrench      	= overridingPlateTopo(overridingPlateTopo(:,1)<BackArcBasinZ(1,1),:);
    elseif strcmp(PLATE.SubPolarity,'left')
        topoBasinTrench        	= overridingPlateTopo(overridingPlateTopo(:,1)>BackArcBasinZ(1,1),:);
    end
    if isempty(topoBasinTrench)
        topoBasinTrench         = [NaN NaN]; %not found
    end
    
    %find max height of volcanic arc
    if ~isnan(BackArcBasinZ(1,1))
        IslandArcHeight         = topoBasinTrench(topoBasinTrench(:,2)==max(topoBasinTrench(:,2)),:);
    else
        IslandArcHeight         = [NaN NaN]; %not found
        if SWITCH.Verbose; warning off backtrace; warning('>>>> volcanic arc not found: no back-arc basin!'); warning on backtrace; end
    end
    if ~isnan(IslandArcHeight(1,1)) && ...
            ((IslandArcHeight(1,1)<=TrenchX(1,1) && strcmp(PLATE.SubPolarity,'right')) || ...
            (IslandArcHeight(1,1)>=TrenchX(1,1) && strcmp(PLATE.SubPolarity,'left')))
        IslandArcHeight         = [NaN NaN]; %not found
        if SWITCH.Verbose; warning('>>>> volcanic arc not found!'); end
    end
    %CHECK FOR PROBLEMS
    if ~isnan(IslandArcHeight(1,1)) && ...
            ((BackArcBasinZ(1,1)<=IslandArcHeight(1,1)+2*dxp && strcmp(PLATE.SubPolarity,'right')) || ... %they are very close
            (BackArcBasinZ(1,1)>=IslandArcHeight(1,1)-2*dxp && strcmp(PLATE.SubPolarity,'left')))
        BackArcBasinZ           = [NaN NaN]; %not found
        if SWITCH.Verbose; warning('>>>> back-arc depression not found!'); end
        IslandArcHeight         = [NaN NaN]; %not found
        if SWITCH.Verbose; warning('>>>> volcanic arc not found!'); end
    end
    if ~isnan(IslandArcHeight(1,1)) && ...
            ((BackArcBasinZ(1,2)>criticalUPz && IslandArcHeight(1,2)>criticalUPz) || ... %they are both too shallow
            isnan(BackArcBasinZ(1,1))) %back-arc depression does not exist
        IslandArcHeight         = [NaN NaN]; %not found
        if SWITCH.Verbose; warning('>>>> volcanic arc not found!'); end
    end
    
    %find segment to the right of the trench deep
    x_width4forelandbulge       = 500e3; %[m]  <<<<<<<<<<<<<<<<<<<<<<<
    x_width4forelandbulge       = x_width4forelandbulge*GRID.m2p; %conversion to plotting values
    
    if strcmp(PLATE.SubPolarity,'right')
        LowerPlateTopo          = subductingPlateTopo(subductingPlateTopo(:,1)<TrenchX(1,1),:);
        topoAroundForebulge    	= LowerPlateTopo(LowerPlateTopo(:,1)>(TrenchX(1,1)-x_width4forelandbulge),:);
    elseif strcmp(PLATE.SubPolarity,'left')
        LowerPlateTopo          = subductingPlateTopo(subductingPlateTopo(:,1)>TrenchX(1,1),:);
        topoAroundForebulge   	= LowerPlateTopo(LowerPlateTopo(:,1)<(TrenchX(1,1)+x_width4forelandbulge),:);
    end
    
    %find foreland bulge
    forelandBulge               = topoAroundForebulge(topoAroundForebulge(:,2)==max(topoAroundForebulge(:,2)),:);
    %check for problems
    if isempty(forelandBulge)
        forelandBulge        	= [NaN NaN]; %not found
        if SWITCH.Verbose; warning('>>>> forelandBulge not found!'); end
    end
    
    %find minimum to the right of the foreland bulge
    if strcmp(PLATE.SubPolarity,'right')
        topoOutsideForeBulge  	= LowerPlateTopo(LowerPlateTopo(:,1)<forelandBulge(:,1),:);
    elseif strcmp(PLATE.SubPolarity,'left')
        topoOutsideForeBulge 	= LowerPlateTopo(LowerPlateTopo(:,1)>forelandBulge(:,1),:);
    end
    minOutsideForeBulge      	= topoOutsideForeBulge(topoOutsideForeBulge(:,2)==min(topoOutsideForeBulge(:,2)),:);
    
    %find back-arc depression area/volume
    if isnan(IslandArcHeight(1,1)) || isnan(BackArcBasinZ(1,1)) || isnan(BackArcBasinExtent(1,1))
        backArcArea             = NaN;
    else
        if strcmp(PLATE.SubPolarity,'right')
            topoUPoutsideVolcArc= overridingPlateTopo(overridingPlateTopo(:,1)>IslandArcHeight(1,1),:);
            topoBasin0        	= topoUPoutsideVolcArc(topoUPoutsideVolcArc(:,1)<BackArcBasinExtent(:,1),:);
        elseif strcmp(PLATE.SubPolarity,'left')
            topoUPoutsideVolcArc= overridingPlateTopo(overridingPlateTopo(:,1)<IslandArcHeight(1,1),:);
            topoBasin0        	= topoUPoutsideVolcArc(topoUPoutsideVolcArc(:,1)>BackArcBasinExtent(:,1),:);
        end
        topoBasin               = topoBasin0;
        %remove UP-mean (which is actually the z value of at the x-indundation point:
        topoBasin(:,2)          = abs( min(topoBasin0(:,2) - topoBasin0(1,2) , 0) );
        topoBasinPlot           = topoBasin;
        topoBasinPlot(:,2)      = -(topoBasin(:,2) - topoBasin0(1,2));
        backArcArea             = trapz(topoBasin(:,1),topoBasin(:,2)); %integrate topography [km^2]
    end
    
    catch erMe
        warning(['   ',STYLE.SCHAR.indicationRightArrow,' Topography characteristics could not be diagnosed due to the following error: ',erMe.message]);
        SUBTOPO.Successful      = false;
        return
        
    end
end

try %topography characteristics routine is specialised for high-resolution cases
    %find upper plate tilt
    useSmoothTopo               = true;         
    shiftXawayFromTrench        = -90e3;        %[m], FURTHER shift away from UP-measurement point & trench <<<<<<<<<<<<<<<<<<
    numPoints4mean              = 7;            %adjust the number of points to take for the mean topo value <<<<<<<<<<<<<<<<<<<
    horizontalSpacing           = 400e3;        %[m], adjust horizontal spacing between measurement points here <<<<<<<<<<<<<<<<<<<
    
    if strcmp(GRID.Type,'spherical2D')
        shiftXawayFromTrench    = shiftXawayFromTrench*1.5;    %in spherical geometry, dynamic slab-surface interaction acts further away from the trench
    end
    
    %dimension conversion
    shiftXawayFromTrench        = shiftXawayFromTrench*GRID.m2p;
    horizontalSpacing           = horizontalSpacing*GRID.m2p;
    
    shiftXawayFromTrench        = ceil(shiftXawayFromTrench/(mean(GRID.dxp(:,:,idxSealevel))))+1; %convert to index
    horizontalSpacing           = ceil(horizontalSpacing/(mean(GRID.dxp(:,:,idxSealevel))))+1; %convert to index
    idxUPnear                   = PLATE.idxUP +shiftXawayFromTrench*PLATE.SubPolarityNum;
    idxUPfar                    = PLATE.idxUP +shiftXawayFromTrench*PLATE.SubPolarityNum +horizontalSpacing*PLATE.SubPolarityNum;
    
    UPtiltAngle = ones(size(idxUPfar,2),numPoints4mean)*NaN;
    for ix=1:numPoints4mean
        %constantly move further away from the trench
        idxUPnear               = max(1,min(size(topo,1),idxUPnear+ix*PLATE.SubPolarityNum));
        idxUPfar                = max(1,min(size(topo,1),idxUPfar+ix*PLATE.SubPolarityNum));
        if useSmoothTopo
            topoUPnear              = topoSmoothest(idxUPnear,1);
            topoUPfar               = topoSmoothest(idxUPfar,1);
        else
            topoUPnear              = topo(idxUPnear,1);
            topoUPfar               = topo(idxUPfar,1);
        end
        xUPnear                 = xp(idxUPnear,1);
        xUPfar                  = xp(idxUPfar,1);
        if ix==1
            xnear = xUPnear; znear = topoUPnear;
        elseif ix==numPoints4mean
            xfar = xUPfar; zfar = topoUPfar;
        end
        UPtiltAngle(:,ix)       = atand(topoUPfar-topoUPnear)./abs(xUPnear-xUPfar); %angle to the horizon (in degrees), positive means down-dip
    end
    UPtiltAngleMean             = median(UPtiltAngle,2);

catch erMe
    warning(['   ',STYLE.SCHAR.indicationRightArrow,' Upper-plate tilt characteristics could not be diagnosed due to the following error: ',erMe.message]);
    return
    
end

%function output
if exist('overridingPlateTopo','var');      SUBTOPO.overridingPlate     = overridingPlateTopo;      end     %in [plotting dimensions]
if exist('overridingPlateTopo','var');      SUBTOPO.subductingPlate     = subductingPlateTopo;      end  	%in [plotting dimensions]
if exist('TrenchDepth','var');              SUBTOPO.TrenchZ             = TrenchDepth;              end   	%measurement points in [plotting dimensions]
if exist('TrenchX','var');                  SUBTOPO.TrenchX             = TrenchX;                  end   	%measurement points in [plotting dimensions]
if exist('BackArcBasinZ','var');            SUBTOPO.BackArcBasinX       = BackArcBasinZ(1,1);       end   	%measurement points in [plotting dimensions]
if exist('BackArcBasinZ','var');            SUBTOPO.BackArcBasinZ       = BackArcBasinZ(1,2);       end   	%measurement points in [plotting dimensions]
if exist('BackArcBasinExtent','var');       SUBTOPO.BackArcBasinExtentX	= BackArcBasinExtent(1,1); 	end   	%measurement points in [plotting dimensions]
if exist('BackArcBasinExtent','var');       SUBTOPO.BackArcBasinExtentZ	= BackArcBasinExtent(1,2);	end   	%measurement points in [plotting dimensions]
if exist('Inundation','var');               SUBTOPO.InundationX      	= Inundation(1,1);          end   	%measurement points in [plotting dimensions]
if exist('Inundation','var');               SUBTOPO.InundationZ      	= Inundation(1,2);          end   	%measurement points in [plotting dimensions]
if exist('Inundation','var');               SUBTOPO.InundationDistance 	= abs(Inundation(1,1)-TrenchX(1,1));	end   	%measurement points in [plotting dimensions]
if exist('IslandArcHeight','var');        	SUBTOPO.IslandArcX          = IslandArcHeight(1,1);  	end   	%measurement points in [plotting dimensions]
if exist('IslandArcHeight','var');          SUBTOPO.IslandArcZ          = IslandArcHeight(1,2);  	end   	%measurement points in [plotting dimensions]
if exist('forelandBulge','var');            SUBTOPO.ForeBulgeX        	= forelandBulge(1,1);       end   	%measurement points in [plotting dimensions]
if exist('forelandBulge','var');            SUBTOPO.ForeBulgeZ        	= forelandBulge(1,2);       end  	%measurement points in [plotting dimensions]
if exist('backArcArea','var');              SUBTOPO.BackArcBasinArea   	= backArcArea;              end 	%measurement points in [plotting dimensions]
SUBTOPO.UPtiltAngleXnear	= xnear;                    %measurement points in [plotting dimensions]
SUBTOPO.UPtiltAngleXfar     = xfar;
SUBTOPO.UPtiltAngleZnear	= znear;
SUBTOPO.UPtiltAngleZfar     = zfar;
SUBTOPO.UPtiltAngle         = UPtiltAngleMean;          %UP-tilt angle in [degrees]
end



%%                                            STAGNANT-LID DIAGNOSTICS 1.0
%                                                Fabio Crameri, 19.05.2017
function [LID,PLOT] = f_stagnantLidDiagnostics(FILE,GRID,LID,SWITCH,PLOT)
% outputs LID.maxYieldDepth = NaN if not successful

%defaults
defmDataExists              = true;
LID.maxYieldDepth           = NaN;
LID.maxYieldDepthFraction   = NaN;

%load data DEFORMATION MECHANISM
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Deformation mechanism';
DATA.FieldAbbreviation      = 'DEFM';
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    DEFM_3D        = PLOT.DEFM_3Dyin; DEFM_3Dyang   = PLOT.DEFM_3Dyang;
else %all other grid types
    DEFM_3D        = PLOT.DEFM_3D;
end
clearvars dummy
if defmDataExists
    %find maximum brittle yielding depth
    if strcmp(GRID.Type,'yinyang')
        LID.maxYieldDepth           = max(max(GRID.Z_3Dp(DEFM_3D==3)),max(GRID.Z_3Dp(DEFM_3Dyang==3)));
        %LID.maxYieldDepthFraction   = mean([LID.thickness,LID.thicknessYang])/LID.maxYieldDepth;
        LID.maxYieldDepthFraction   = LID.LABdepth/LID.maxYieldDepth;
    else
        LID.maxYieldDepth           = max(GRID.Z_3Dp(DEFM_3D==3));
        %LID.maxYieldDepthFraction   = mean(LID.thickness)/LID.maxYieldDepth;
        LID.maxYieldDepthFraction   = 1/LID.LABdepth*LID.maxYieldDepth;
    end
end

end



%% RMS
function y = rms(x)
%RMS Root mean squared value
y = sqrt(mean(x.*conj(x)));
end
