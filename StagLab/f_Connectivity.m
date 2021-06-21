
%%                                                    2-D CONNECTIVITY 1.2
%
%       . calls MinVolEllipse
%       . calls Ellipse_plot (only for figure check)
%
%                                                Fabio Crameri, 17.04.2017
%
%% NOTES
% variable grid spacing is not accounted for..........
% dx~=dy not accounted for...........
% by decreasing CONAREA.nxSubgrid you can check to remove double boundaries......
%
% calculates connected areas by accounting for f_INPUT.addFactor times the area of
% the original plot: one part added at x=0 and one other part added at x=end
%
% field indicating +1 for areas and 0 for the residual field

function [CONAREA] = f_Connectivity(Array,CONAREA,GRID)

%% DEFAULT INPUT
if ~isfield(CONAREA,'threshold');	CONAREA.threshold	= 2.0e-6; end           %threshold value defining area to check
if ~isfield(CONAREA,'thresLogical');CONAREA.thresLogical= 'bigger'; end         %'bigger' or 'smaller' than threshold
if ~isfield(CONAREA,'addFactor');	CONAREA.addFactor	= 0.25; end             %percent of total x-extend that is added at x=1 and x=end
if ~isfield(CONAREA,'xPeriodic');	CONAREA.xPeriodic	= logical(1); end       %account for periodic x-side bounday
if ~isfield(CONAREA,'yPeriodic');	CONAREA.yPeriodic	= logical(1); end       %account for periodic y-side bounday
if ~isfield(CONAREA,'connNeighbor');CONAREA.connNeighbor= 8; end                %connected neighbourhood for bwconncomp
if ~isfield(CONAREA,'nxSubgrid');   CONAREA.nxSubgrid   = 16; end           	%size of subgrid in #gridpoints (has to be power of 2!)

figureCheck = logical(0); %test figure

%% VARIABLE CONVERSIONS
X_3D        = GRID.X_3Dp;
Y_3D        = GRID.Y_3Dp;
conn        = CONAREA.connNeighbor; %connected neighbourhood
NXsubgrid   = CONAREA.nxSubgrid;
NYsubgrid   = NXsubgrid;

%% ERROR CHECKS
%subgrid cannot be bigger than main grid
if size(Array,1)<=NXsubgrid || size(Array,2)<=NYsubgrid
    error(['Subgrid (',num2str(NXsubgrid),'x',num2str(NYsubgrid),') is too big compared to the input array (',...
        num2str(size(Array,1)),'x',num2str(size(Array,2)),')! Reduce the size of CONAREA.nxSubgrid!'])
end
%subgrid size must be of power 2
[dummy,~] = log2(CONAREA.nxSubgrid);
if dummy~=0.5
    dummy = CONAREA.nxSubgrid;
    CONAREA.nxSubgrid = 2^nextpow2(dummy);
    warning(['CONAREA.nxSubgrid = ',num2str(dummy),' is not of power 2! It has now been changed to CONAREA.nxSubgrid = ',num2str(CONAREA.nxSubgrid),'.'])
end

%% APPLY CRITERION
AREASconn = zeros(size(Array));
if strcmp(CONAREA.thresLogical,'bigger')
    AREASconn = Array>CONAREA.threshold;
    
elseif strcmp(CONAREA.thresLogical,'smaller')
    AREASconn = Array<CONAREA.threshold;
end

%% ADD GHOSTPOINTS TO ACCOUNT FOR PERIODIC SIDE BOUNDARIES
%check for areas crossing periodic boundaries:
add_factor  = CONAREA.addFactor;      %percent of total x-extend that is added at x=1 and x=end
%add some slices at x=1 and x=end
AREASorig   = AREASconn;
nx_org      = size(AREASorig,1);
nx_add      = nx_org*add_factor;
ny_org      = size(AREASorig,2);
ny_add      = ny_org*add_factor;

if CONAREA.xPeriodic
    AREASconn = [AREASconn((end-nx_add+1):end,:);   AREASconn;   AREASconn(1:nx_add,:)];
    X_3D = [X_3D((end-nx_add+1):end,:,:);   X_3D;   X_3D(1:nx_add,:,:)];
    Y_3D = [Y_3D((end-nx_add+1):end,:,:);   Y_3D;   Y_3D(1:nx_add,:,:)];
end
if CONAREA.yPeriodic
    AREASconn = [AREASconn(:,(end-ny_add+1):end)   AREASconn   AREASconn(:,1:ny_add)];
    X_3D = [X_3D(:,(end-ny_add+1):end,:)   X_3D   X_3D(:,1:ny_add,:)];
    Y_3D = [Y_3D(:,(end-ny_add+1):end,:)   Y_3D   Y_3D(:,1:ny_add,:)];
end

%% FIND CONNECTED AREAS
CC = bwconncomp(AREASconn,conn);
numPixels	= cellfun(@numel,CC.PixelIdxList); %number of connected pixels in each connected area

%% REMOVE SMALL ANOMALIES
% [size_largest,ind_largest]	= max(numPixels);    	%find single most largest area
size_threshold  = 0.1; %define threshold here (non-dimensional width of square area); e.g., 0.1 = 304 km^2
dx_nd           = max(GRID.X_3Dnd(:))/size(GRID.X_3Dnd,1); %dx_nd
size_threshold  = (size_threshold/dx_nd)^2; %no. of pixels spanning defined min. area
ind_large   	= find(numPixels<size_threshold);
for i_ind=1:size(ind_large,2)
    AREASconn(CC.PixelIdxList{ind_large(i_ind)}) = 0;
end

%% REMOVE GHOSTPOINTS - REVERT TO ORIGINAL
if CONAREA.xPeriodic
    AREASconn 	= AREASconn((nx_add+1):end-nx_add,:);
    X_3D     	= X_3D((nx_add+1):end-nx_add,:,:);
    Y_3D     	= Y_3D((nx_add+1):end-nx_add,:,:);
end
if CONAREA.yPeriodic
    AREASconn 	= AREASconn(:,(ny_add+1):end-ny_add);
    X_3D     	= X_3D(:,(ny_add+1):end-ny_add,:);
    Y_3D     	= Y_3D(:,(ny_add+1):end-ny_add,:);
end

%% DERIVE CENTER LINE OF AREAS
x2d = X_3D(:,:,1); %make 2-D
y2d = Y_3D(:,:,1);

%make line out of area %APPLY THIS ONLY FOR CONNECTED AREAS!
Aline = bwmorph(AREASconn,'thin',Inf);
if figureCheck
    figure(11),clf
    subplot(1,3,1)
    plot(y2d(AREASconn),x2d(AREASconn),'.')
    axis equal
    subplot(1,3,2)
    % imshow(Aline)
    plot(y2d(Aline),x2d(Aline),'.')
    axis equal
    subplot(1,3,3)
    %axis equal
end
AREASconn = Aline;

%% UPDATE CONNECTED-AREA LIST (after removing small anomalies and creating lines)
ConA     	= bwconncomp(AREASconn,conn); %find connected area for an array of any dimension
numPixels1	= cellfun(@numel,ConA.PixelIdxList); %number of connected pixels in each connected area

%% FIND LINE POINTS AND SELECTED CHARACTERISTIC BOUNDARY POINTS
areaX = zeros(size(X_3D,1)*size(X_3D,2),size(ConA.PixelIdxList,2))*NaN; %set to maximum possible size
areaY = areaX;
LineP = zeros(size(ConA.PixelIdxList,2), max(numPixels1(:)), 2)*NaN;
CenterP = zeros(size(ConA.PixelIdxList,2), GRID.nx/NXsubgrid*GRID.ny/NYsubgrid, 2)*NaN;
StrikeP = CenterP; NormalP = CenterP;
OutlineP = zeros(size(ConA.PixelIdxList,2), GRID.nx*GRID.ny, 2)*NaN;
isubareaMaxNo = 0; kmax = 0;
for iarea=1:size(ConA.PixelIdxList,2) %loop through connected areas
    LineP(iarea,1:numPixels1(1,iarea),1) = x2d(ConA.PixelIdxList{iarea});
    LineP(iarea,1:numPixels1(1,iarea),2) = y2d(ConA.PixelIdxList{iarea});
    
    Aline(:,:) = 0; %erase
    Aline(ConA.PixelIdxList{iarea}) = 1; %only add line of current area
    AREASoutline = bwmorph(Aline,'thicken',1); %make 1 pixel thicker all around  -  plot(y2d(AREASoutline),x2d(AREASoutline),'.')
    %AREASoutline = bwmorph(AREASoutline,'remove'); %remove inside parts
    yy = y2d(AREASoutline); xx = x2d(AREASoutline);
    k = boundary(yy,xx,1.0);
    OutlineP(iarea,1:size(k,1),1:2) = [xx(k),yy(k)]; %plot(OutlineP(1,:,2),OutlineP(1,:,1),'r')
    kmax = max(kmax,size(k,1));
    
    %go into subgrid
    isubarea=0;
    for ixsubarea=1:NXsubgrid:size(x2d,1)
        for iysubarea=1:NYsubgrid:size(x2d,2)
            %CREATE SUBAREAS
            x2dsub = zeros(size(x2d))*NaN;
            y2dsub = zeros(size(y2d))*NaN;
            upperX = ixsubarea+NXsubgrid-1;
            upperY = iysubarea+NYsubgrid-1;
            
            x2dsub(ixsubarea:upperX,iysubarea:upperY) = x2d(ixsubarea:upperX,iysubarea:upperY);
            y2dsub(ixsubarea:upperX,iysubarea:upperY) = y2d(ixsubarea:upperX,iysubarea:upperY);
            
            areaX(1:size(x2d(ConA.PixelIdxList{iarea}),1),iarea) = x2dsub(ConA.PixelIdxList{iarea});
            areaY(1:size(y2d(ConA.PixelIdxList{iarea}),1),iarea) = y2dsub(ConA.PixelIdxList{iarea});
            
            area_current = [ areaX(~isnan(areaX(:,iarea)),iarea), areaY(~isnan(areaY(:,iarea)),iarea) ]; %[x,y]
            
            %CHECK IF SUBAREA IS EMPTY
            if isempty(area_current) || max(area_current(:,1))==0
                continue %skip current subarea
            else
                isubarea = isubarea+1; %count subareas used
            end
            if figureCheck
                plot(area_current(:,1),area_current(:,2),'b.')
            end
            
            if min(area_current(:,1))==max(area_current(:,1)) %data points lie on a line in x direction
                c(1,1) = area_current(1,1);
                c(2,1) = mean(area_current(:,2));
                unitVecNormal = [1;0]; %????????ADJUST (might need to swap values - CHECK THAT!)
                unitVecStrike = [0;1];
                if figureCheck
                    hold on
                    plot(c(1),c(2),'ro')
                end
                
            elseif min(area_current(:,2))==max(area_current(:,2)) %data points lie on a line in y direction
                c(1,1) = area_current(1,2);
                c(2,1) = mean(area_current(:,1));
                unitVecNormal = [0;1]; %????????ADJUST (might need to swap values - CHECK THAT!)
                unitVecStrike = [1;0];
                if figureCheck
                    hold on
                    plot(c(1),c(2),'ro')
                end
                
            else
                %APPLY MINIMUM VOLUME ELLIPSE to point cloud
                P = area_current';
                if size(P,2)>2 %does not work if only 2 points
                    % To reduce the computation time, work with the boundary points only:
                    K = convhulln(P');
                    K = unique(K(:));
                    Q = P(:,K);
                else
                    Q = P;
                end
                [A, c] = MinVolEllipse(Q, .01);
                if figureCheck
                    hold on
                    plot(c(1),c(2),'ro')
                    hold on
                    Ellipse_plot(A, c);
                end
                
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
            end
            
            %KEEP DATA VALUES of center points (c)
            CenterP(iarea,isubarea,:) = c;
            if size(unitVecStrike,2)>1 || size(unitVecNormal,2)>1
                unitVecStrike=NaN; unitVecNormal=NaN; %ellipse is circle and there is no major axis found
            end
            StrikeP(iarea,isubarea,:) = unitVecStrike;
            NormalP(iarea,isubarea,:) = unitVecNormal;
        end
    end
    isubareaMaxNo = max(isubareaMaxNo,isubarea);
end

%% FINAL ADJUSTMENTS
%centerP(iarea,...) are not necessarily the same size, hence need to remove zeros
CenterP(CenterP==0) = NaN;
CenterP(:,isubareaMaxNo+1:end,:) = [];
StrikeP(:,isubareaMaxNo+1:end,:) = [];
NormalP(:,isubareaMaxNo+1:end,:) = [];
OutlineP(:,kmax+1:end,:) = [];

%% CALCULATING ADDITIONAL PARAMETERS
% ANGLE: tan(alpha) = x/y;
StrikeAngleP = rad2deg( atan(StrikeP(:,:,1)./StrikeP(:,:,2)) );
NormalAngleP = rad2deg( atan(NormalP(:,:,1)./NormalP(:,:,2)) );
%strike points in sensible plot values
plotFac = x2d(10,1,1)-x2d(1,1,1);
% StrikePlot1 = CenterP-StrikeP*plotFac;
% StrikePlot2 = CenterP+StrikeP*plotFac;
NormalPlot1 = CenterP-NormalP*plotFac;
NormalPlot2 = CenterP+NormalP*plotFac;

if figureCheck
    figure(22); clf
    for i=1:size(CenterP,1)
        plot(LineP(i,:,1),LineP(i,:,2),'r.')
        %         hold on
        %         plot(CenterP(i,:,1),CenterP(i,:,2),'-')
        for j=1:size(CenterP,2)
            hold on
            text(CenterP(i,j,1),CenterP(i,j,2),num2str(StrikeAngleP(i,j),2))
            
            plot([NormalPlot1(i,:,1);NormalPlot2(i,:,1)],...
                [NormalPlot1(i,:,2);NormalPlot2(i,:,2)],'g-')
            %             plot([StrikePlot1(i,:,1);StrikePlot2(i,:,1)],...
            %                 [StrikePlot1(i,:,2);StrikePlot2(i,:,2)],'g-')
        end
    end
    grid on
    axis equal
end

%% FUNCTION OUTPUT
CONAREA.LineP           = LineP;            %line points [area,no.points,x_i]
CONAREA.OutlineP        = OutlineP;         %outline of line points [area,x_i,y_i]
CONAREA.CenterP         = CenterP;          %center points [area,subarea,x_i]
CONAREA.StrikeP         = StrikeP;          %unit vectors [area,subarea,x_i] at center points describing direction of strike
CONAREA.NormalP         = NormalP;          %unit vectors [area,subarea,x_i] at center points describing direction of normal
CONAREA.StrikeAngleP    = StrikeAngleP;     %strike angle [area,subarea,deg] at center points
CONAREA.NormalAngleP    = NormalAngleP;     %normal2strike angle [area,subarea,deg] at center points
if ~isempty(LineP)
    for iarea=1:size(LineP,1)
        CONAREA.LengthAreas(iarea,1) = size(~isnan(LineP(iarea,:,1)),2) *GRID.dxp(1,1,1); %[plotting dimension]  %NOT CORRECT IF DX/=DY !!!! ALSO NEGLECTS DIAGONAL POINTS... ALSO ASSUMES CONSTANT dx...
    end
    CONAREA.LengthTotal     = sum(CONAREA.LengthAreas); %[plotting dimension]
else %if nothing is found
    CONAREA.LengthAreas = NaN;
    CONAREA.LengthTotal = NaN;
end
















