
%%                                                              VARIA 1.20
%
% VARIA.Task    = 'create stream line';
% [PLOT,VARIA] = f_Varia(VARIA,FILE,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE,SETUP);
% 
%                                                Fabio Crameri, 25.07.2017
                    
function [PLOT,VARIA] = f_Varia(VARIA,FILE,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE,SETUP)

%% TASKS
if strcmp(VARIA.Task,'create stream line')
    [PLOT,VARIA] = f_createStreamLine(VARIA,FILE,GRID,SWITCH,PLOT);
    
elseif strcmp(VARIA.Task,'derive plate-base topography')
    [PLOT,VARIA] = f_PlateBaseTopography(VARIA,FILE,GRID,SWITCH,PLOT,SETUP);
    
end

end


%% STREAM LINES
function [PLOT,VARIA] = f_createStreamLine(VARIA,FILE,GRID,SWITCH,PLOT)

stepsize            = 0.1;
maxNumberVertices   = PLOT.streamLength;     %streamline length in [#cells]
maxNumberVertices   = maxNumberVertices/stepsize;

sx2d                = PLOT.streamStartline(1,:);
sy2d                = PLOT.streamStartline(2,:);
sz2d                = PLOT.streamStartline(3,:);

if strcmp(GRID.Type,'Cartesian')
    DATA.Task                   = 'ImportFieldData';
    DATA.Field2Import           = 'Velocity';
    DATA.FieldAbbreviation      = 'V';
    [~,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
    
    VX_3D = PLOT.VX_3D; VY_3D = PLOT.VY_3D; VZ_3D = PLOT.VZ_3D;
    if strcmp(GRID.Dim,'3-D')
    elseif strcmp(GRID.Dim,'2-D')
        if strcmp(GRID.Type,'Cartesian')        
            vx2d(:,:)   = VX_3D(:,1,:);
            %vy2d(:,:)  = VY_3D(:,1,:);
            vz2d(:,:)   = VZ_3D(:,1,:);
            vz2d        = -vz2d;    %adjustment for flip in depth vector
        elseif strcmp(GRID.Type,'spherical2D')
%             vx2d(:,:)   = PLOT.VX_3Ds(:,1,:);
%             %vy2d(:,:)  = PLOT.VY_3Ds(:,1,:);
%             vz2d(:,:)   = PLOT.VZ_3Ds(:,1,:);
            
            vx2d(:,:)   = VX_3D(:,1,:);
            %vy2d(:,:)  = VY_3D(:,1,:);
            vz2d(:,:)   = VZ_3D(:,1,:);
        end
        %grid variables
        x2d         = GRID.x2dp;
        z2d         = GRID.z2dp;
%         if strcmp(GRID.Type,'spherical2D')
%             x2d     = GRID.x2ds;
%             z2d     = GRID.z2ds;
%         end
    end
elseif strcmp(GRID.Type,'spherical2D')
    error('fc: spherical2D streamline not yet implemented...')
end

% if strcmp(GRID.Type,'Cartesian')
    sx2d                = min(x2d(:)) +sx2d.*(max(x2d(:))-min(x2d(:)));	%convert <0-1> values to actual values of x and z
    sz2d                = min(z2d(:)) +sz2d.*(max(z2d(:))-min(z2d(:)));
% elseif strcmp(GRID.Type,'spherical2D')
%     sx2d                = min(GRID.x2d_nd(:)) +sx2d.*(max(GRID.x2d_nd(:))-min(GRID.x2d_nd(:)));	%convert <0-1> values to actual values of x and z
%     sz2d                = min(GRID.z2dp(:)) +sz2d.*(max(GRID.z2dp(:))-min(GRID.z2dp(:))); %z-values are flipped and need to be flipped back
% end

numStartingPoints   = PLOT.streamNumStartPoints;
if numStartingPoints<=1
    stepFactorX    	= sx2d(1,2);
    stepFactorZ    	= sz2d(1,2);
else
    stepFactorX    	= (sx2d(1,2)-sx2d(1,1))/(numStartingPoints-1);
    stepFactorZ   	= (sz2d(1,2)-sz2d(1,1))/(numStartingPoints-1);
end
if sx2d(1,2)-sx2d(1,1)<=0; stepFactorX = 1; end
if sz2d(1,2)-sz2d(1,1)<=0; stepFactorZ = 1; end

sxLine              = sx2d(1,1):stepFactorX:sx2d(1,2);
if size(sxLine,2)==1
    szLine          = sz2d(1,1):stepFactorZ:sz2d(1,2);
else
    szLine          = sz2d(1,1):sz2d(1,2)/size(sxLine,2):sz2d(1,2);     %must have the same size than sxline
end
if size(sxLine,2)==1; sxLine = ones(1,size(szLine,2))*sxLine; end   %in case of vertical lines
if size(szLine,2)==1; szLine = ones(1,size(sxLine,2))*szLine; end   %in case of horizontal lines
if strcmp(GRID.Dim,'3-D')
    %hs = streamline(y2d(1,:,1),x2d(:,1,1),z2d(1,1,:),vy2d',vx2d',vz2d',sy2d,sx2d,sz2d);
    streamData      = stream2(x2d(:,1),z2d(1,:),vx2d',vz2d',sxLine,szLine,[stepsize,maxNumberVertices]);
    
elseif strcmp(GRID.Dim,'2-D')
    %hs = streamline(x2d(:,1),z2d(1,:),vx2d',vz2d',sxline,szline);
        
    streamData      = stream2(x2d',z2d',vx2d',vz2d',sxLine,szLine,[stepsize,maxNumberVertices]);
%     streamData      = stream2(x2d(:,1),z2d(1,:),vx2d',vz2d',sxLine,szLine,[stepsize,maxNumberVertices]);
end
if strcmp(GRID.Type,'spherical2D')
    for id=1:length(streamData)
        %coordinate conversion
        sradius                 = abs(max(GRID.z2dp(:))-streamData{1,id}(:,2)) +GRID.rcmb_p;
        dummy                   = streamData{1,id}(:,1);
        streamData{1,id}(:,1) 	= sradius.*-sin(dummy); %.*cos(y2d_nd);
        streamData{1,id}(:,2)  	= sradius.*-cos(dummy);
    end
end

%Output
VARIA.streamData    = streamData;

VARIA.sxLine        = sxLine;
VARIA.szLine        = szLine;
VARIA.sx2d          = sx2d;
VARIA.sz2d          = sz2d;
VARIA.vx2d          = vx2d;
VARIA.vz2d          = vz2d;
VARIA.x2d           = x2d;
VARIA.z2d           = z2d;
VARIA.stepsize    	= stepsize;
VARIA.maxNumVertices= maxNumberVertices;
end



%% PLATE-BASE TOPOGRAPHY
function [PLOT,VARIA] = f_PlateBaseTopography(VARIA,FILE,GRID,SWITCH,PLOT,SETUP)

normalise2mean      = logical(1);

%% READ TEMPERATURE
%temperature
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = 'Temperature';
DATA.FieldAbbreviation      = 'T';
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    T_3D        = PLOT.T_3Dyin; T_3Dyang   = PLOT.T_3Dyang;
else %all other grid types
    T_3D        = PLOT.T_3D;
end

%% DIMENSIONALISATION
T_3D                = T_3D .*SETUP.Tscale; %[K]

%% CHECK
if isnan(T_3D(1)) || strcmp(GRID.Dim,'3-D')
    if strcmp(GRID.Dim,'3-D')
        warning('3-D plate isoline topography plot not implemented yet!')
    end
    VARIA.PlateBaseTopo     = NaN;
    return
    
end

%% CRITICAL VALUES to define isoline
LITHO.Tisovalue     = PLOT.lithoThickness(1,1);     %Temperature isovalue
LITHO.zmax          = PLOT.lithoThickness(1,2);     %depth in [km]

[LITHO] = f_PlateThickness(T_3D,GRID.Z_3Dp,LITHO);     % derive lithosphere thickness by T contour:

%% NORMALISE THE VALUES TO THE MEAN
if normalise2mean
    topo2d          = -LITHO.thickness +mean(LITHO.thickness(:));
else
    topo2d          = LITHO.thickness;
end

%% SMOOTHING
% for iloopComp = 1:nloopComp
%     if SWITCH.plotDifference && strcmp(TOPO.difference,'comp')
%         if iloopComp==1
%             PLOT.istopo2d_B = topo2d;
%             topo2d = PLOT.istopo2d_A;
%         else
%             topo2d = PLOT.istopo2d_B;
%         end
%     end
%     topoSmooth_2 = zeros(size(topo2d));
%     topoSmooth_1 = topo2d;
%     
%     if strcmp(GRID.Dim,'2-D')   % 2-D
%         for jj=1:2
%             topoSmooth_2(:,1) = smooth(topoSmooth_1(:,1),5,'moving');
%             topoSmooth_1 = topoSmooth_2;
%         end
%         topoSmooth = topoSmooth_2;
%         
%     elseif strcmp(GRID.Dim,'3-D') %3-D
%         if TOPO.ascii
%             z = topo2d(:,1); %cmb topo
%             low = min(x); high = max(x); tix = low:0.01:high;
%             low = min(y); high = max(y); tiy = low:0.01:high;
%             [xi,yi] = meshgrid(tix,tiy);
%             zi = griddata(x,y,z,xi,yi);
%         else
%             xi = y_PLOT;
%             yi = x_PLOT;
%             zi = topo2d(:,:,1); %cmb topo
%         end
%     end
% end

%% OUTPUT
VARIA.PlateBaseTopo     = topo2d;

end
