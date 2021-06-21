
%%                                                     TRACKING PLUMES 1.0
%
% calculates connected areas by accounting for f_INPUT.addFactor times the area of
% the original plot: one part added at x=0 and one other part added at x=end
%
% needs as input (=PLUMES) an evaluation of the residual temperature
% field indicating +1 for hot plumes and -1 for cold plumes and 0 for the residual field
% ( hot plumes are defined by e.g. Tmean+0.5*(Tmax-Tmean) )
%
%                                                Fabio Crameri, 23.03.2016

%% NOTES
%variable grid spacing is not accounted for..........

function [PLUMEShot,PLUMEScold] = f_trackPlumes(PLUMES,f_INPUT)

%% INPUT
% f_INPUT.plotPLUMES      	= logical(0);	%plot comparison of all anomalies and selected plumes
% f_INPUT.addFactor      	= 0.25;         %percent of total x-extend that is added at x=1 and x=end  (should be exponent of 2)
% f_INPUT.upperDepthThresholdHot 	= 0.5;          %top threshold; percentage of nz (0: top, 1: bottom)   ...incl. sticky air!
% f_INPUT.lowerDepthThresholdHot	= 0.9;          %top threshold; percentage of nz (0: top, 1: bottom)
% f_INPUT.upperDepthThresholdCold	= 0.2;          %top threshold; percentage of nz (0: top, 1: bottom)   ...incl. sticky air!
% f_INPUT.lowerDepthThresholdCold	= 0.4;          %top threshold; percentage of nz (0: top, 1: bottom)

%% DEFAULTS
if ~isfield(f_INPUT,'plotPLUMES');          f_INPUT.plotPLUMES          = logical(0);   end     %plot comparison of all anomalies and selected plumes
if ~isfield(f_INPUT,'addFactor');           f_INPUT.addFactor           = 0.25;         end     %percent of total x-extend that is added at x=1 and x=end
if ~isfield(f_INPUT,'upperDepthThresholdHot'); 	f_INPUT.upperDepthThresholdHot 	= 0.5;          end     %top threshold; percentage of nz (0: top, 1: bottom)   ...incl. sticky air!
if ~isfield(f_INPUT,'lowerDepthThresholdHot'); 	f_INPUT.lowerDepthThresholdHot	= 0.9;          end     %top threshold; percentage of nz (0: top, 1: bottom)
if ~isfield(f_INPUT,'upperDepthThresholdCold');   f_INPUT.upperDepthThresholdCold	= 0.2;          end     %top threshold; percentage of nz (0: top, 1: bottom)   ...incl. sticky air!
if ~isfield(f_INPUT,'lowerDepthThresholdCold'); 	f_INPUT.lowerDepthThresholdCold	= 0.4;          end     %top threshold; percentage of nz (0: top, 1: bottom)

%% GRID DIAGNOSTICS
if size(PLUMES,1)==1 || size(PLUMES,2)==1
    gridDim = '2-D';
else
    gridDim = '3-D';
end


%% ACCOUNT FOR PERIODIC BOUNDARIES
%check for plumes crossing periodic boundaries:
add_factor  = f_INPUT.addFactor;      %percent of total x-extent that is added at x=1 and x=end
%add some slices at x=1 and x=end
PLUMESorg   = PLUMES;
nx_org      = size(PLUMESorg,2);
nx_add      = nx_org*add_factor;

PLUMES = [PLUMES(:,(end-nx_add+1):end,:)   PLUMES   PLUMES(:,1:nx_add,:)];

%% CHECK CONNECTIVITY FOR AREA SIZES
CC = bwconncomp(PLUMES);
numPixels	= cellfun(@numel,CC.PixelIdxList); %number of connected pixels in each connected area

%% CONNECTED AREA CRITERION
% *********************** area criterion
% [size_largest,ind_largest]	= max(numPixels);    	%find single most largest area
size_threshold        	= 100*size(PLUMES,3);   %find areas larger than (100*nz) pixels
ind_large               = find(numPixels<size_threshold);

top_threshold_HOT     	= f_INPUT.upperDepthThresholdHot; %top threshold; percentage of nz (0: top, 1: bottom)   ...incl. sticky air!
bottom_threshold_HOT 	= f_INPUT.lowerDepthThresholdHot; %top threshold; percentage of nz (0: top, 1: bottom)
top_threshold_COLD   	= f_INPUT.upperDepthThresholdCold; %top threshold; percentage of nz (0: top, 1: bottom)   ...incl. sticky air!
bottom_threshold_COLD 	= f_INPUT.lowerDepthThresholdCold; %top threshold; percentage of nz (0: top, 1: bottom)
% ***********************

%% REMOVING SMALL ANOMALIES
PLUMES1 = PLUMES;
for i_ind=1:size(ind_large,2)
    PLUMES1(CC.PixelIdxList{ind_large(i_ind)}) = 0;
end

%% SETTING UP HOT AND COLD ANOMALIES
PLUMEShot               = zeros(size(PLUMES1));
PLUMEScold              = zeros(size(PLUMES1));
PLUMEShot(PLUMES1==1)   = 1; %hot plumes
PLUMEScold(PLUMES1==-1) = 1; %cold plumes

%% CHECK CONNECTIVITY FOR PLUMES
% CC1     = bwconncomp(PLUMES1)
CChot   = bwconncomp(PLUMEShot);
CCcold  = bwconncomp(PLUMEScold);

% numPixels	= cellfun(@numel,CC1.PixelIdxList); %number of connected pixels in each connected area
% numPixelsHOT	= cellfun(@numel,CChot.PixelIdxList); %number of connected pixels in each connected area
% numPixelsCOLD	= cellfun(@numel,CCcold.PixelIdxList); %number of connected pixels in each connected area

%% SETTING UP DEPTH MATRIX
% nx = size(PLUMES,1);
% ny = size(PLUMES,2);
nz = size(PLUMES,3);
% x3d=ones(size(PLUMES));
% y3d=ones(size(PLUMES));
z3d=ones(size(PLUMES));

for k=1:nz
    z3d(:,:,k) = k;
end

%% CHECK CONNECTIVITY TO BOUNDARY LAYERS
for plume_kind=1:2
    if plume_kind==1  %HOT
        CC_kind = CChot; top_threshold = top_threshold_HOT; bottom_threshold = bottom_threshold_HOT;
    elseif plume_kind==2  %COLD
        CC_kind = CCcold; top_threshold = top_threshold_COLD; bottom_threshold = bottom_threshold_COLD;
    end
    
    for i_area=1:size(CC_kind.PixelIdxList,2)  %loop all connected areas
        surf_connection = false;
        bot_connection = false;
        
        for ip=1:size(CC_kind.PixelIdxList{1,i_area},1)  %loop all pixels
            i_pixel = CC_kind.PixelIdxList{1,i_area}(ip,1);
            %[x,y,z] = ind2sub(size(PLUMES),i_pixel);
            
            if z3d(i_pixel)>(nz-top_threshold*nz) %is connected to the surface
                surf_connection = true;
            end
            if z3d(i_pixel)<(nz-bottom_threshold*nz) %is connected to the bottom
                bot_connection = true;
            end
            
            if surf_connection && bot_connection
                break   %this is a plume from top to bottom! check next.
            end
        end
        
        if surf_connection && bot_connection
            
        else %remove anomaly: this is no plume   %% CHECK IF THAT MAKES SENSE! SEEMS TO REMOVE TOO MANY PLUMES?!
            if plume_kind==1  %HOT
                PLUMEShot(CC_kind.PixelIdxList{i_area}) = 0;
            elseif plume_kind==2  %COLD
                PLUMEScold(CC_kind.PixelIdxList{i_area}) = 0;
            end
        end
    end
end

%% REMOVING GRID EXTENSION (BACK TO ORIGINAL SIZE)
PLUMES     	= PLUMES(:,(nx_add+1):end-nx_add,:);
% PLUMES1  	  = PLUMES1(:,nx_add+1:end-nx_add,:);
PLUMEShot 	= PLUMEShot(:,(nx_add+1):end-nx_add,:);
PLUMEScold 	= PLUMEScold(:,(nx_add+1):end-nx_add,:);

% CChot   = bwconncomp(PLUMEShot);
% CCcold  = bwconncomp(PLUMEScold);

%% PLOTTING PLUMES SEPERATELY
if f_INPUT.plotPLUMES
    figure(22),clf
    
    x = 1:size(PLUMES,1);
    y = 1:size(PLUMES,2);
    z = 1:size(PLUMES,3);
    [x3d,y3d,z3d] = meshgrid(y,x,z);
    
    phot = patch(isosurface(x3d,y3d,z3d,PLUMEShot,+0.95)); %HOT
    isonormals(x3d,y3d,z3d,PLUMEShot,phot)
    set(phot,'FaceColor','red','EdgeColor','none');
    view(3);
    camlight
    lighting phong  %gouraud
    
    hold on
    pcold = patch(isosurface(x3d,y3d,z3d,PLUMEScold,+0.95)); %COLD
    isonormals(x3d,y3d,z3d,PLUMEScold,pcold)
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
    isosurface(PLUMES,+0.95,'noshare') %HOT
    hold on
    isosurface(PLUMES,-0.95,'noshare') %COLD
    
    title('original')
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal
    axis([0 max(max(max(x3d))) 0 max(max(max(y3d))) 0 max(max(max(z3d)))])
    box on
    grid on
end






















