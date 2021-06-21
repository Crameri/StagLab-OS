
%%                                          FIND LITHOSPHERE THICKNESS 2.12
%
% *** FOR 2-D AND 3-D ***
% and output LAB-topography (i.e., T or eta isosurface)
%
%                                                Fabio Crameri, 22.06.2018
%
%% ABOUT
% Input:
% LITHO.Tisovalue = 1400;
% LITHO.zmax      = 200;
% LITHO.plot      = false;
% z2d: in [km]; surface = 0; bottom = D;
% [LITHO] = f_PlateThickness(var2d,z2d,LITHO);
%
% Output:
% LITHO.thickness
% LITHO.thicknessFull
% LITHO.gridpoints: size(nx,ny), indicates iz,  %only produced if field LITHO.gridpoints already exists
% LITHO.gridpointsFull: size(nx,ny), indicates iz,  %only produced if field LITHO.gridpoints already exists

function [LITHO] = f_PlateThickness(var2d,z2d,LITHO)
%% Defaults
if ~isfield(LITHO,'Tisovalue');     LITHO.Tisovalue = 1600;     end %Temperature isovalue
if ~isfield(LITHO,'zmax');          LITHO.zmax = 200;           end %depth in [km]
if ~isfield(LITHO,'zmin');          LITHO.zmin = 5;             end %depth in [km]
if ~isfield(LITHO,'plot');          LITHO.plot = false;         end %plot the lithosphere thickness

%% Set minimum and maximum values to closest grid values
[~,index] = min(abs(z2d(1,1,:)-LITHO.zmin));
closestValue        = z2d(1,1,index); % Finds first one only!
if closestValue<LITHO.zmin
    LITHO.zmin      = z2d(1,1,index-1); %take the next deeper value
else
    LITHO.zmin      = z2d(1,1,index);
end
[~,index] = min(abs(z2d(1,1,:)-LITHO.zmax));
closestValue        = z2d(1,1,index); % Finds first one only!
if closestValue>LITHO.zmax
    LITHO.zmax      = z2d(1,1,index+1); %take the previous shallower value
else
    LITHO.zmax      = z2d(1,1,index);
end

%% Search Lithosphere Thickness
% var2d(:,:,end)    %is top
% var2d(:,:,1)      %is bottom
thicknessFull = zeros(size(var2d,1),size(var2d,2),1); %[ix,iy,1]
if isfield(LITHO,'gridpoints')
    getGridPoints   = true;
    gridpoints      = ones(size(var2d,1),size(var2d,2),1) .*size(var2d,3); %[ix,iy,1] %default value is surface index
    gridpointsFull  = gridpoints; %[ix,iy,1]
else
    getGridPoints   = false;
end

%optimised version:
mantlePortionsCombined = []; check = false;
for iz=size(var2d,3):-1:1 %from top to bottom
    if z2d(1,1,iz)<LITHO.zmin %make sure it has a minimum depth
        thicknessFull(:,:,1) = LITHO.zmin;
        if getGridPoints; gridpointsFull(:,:) = iz; end %fill with minimum iz
        continue
    end
    %check for hot parts and remember where mantle was found already
    mantlePortions = find( var2d(:,:,iz)>=LITHO.Tisovalue );
    mantlePortionsCombined = union(mantlePortionsCombined,mantlePortions); %combine data with no repetitions
    %find cold parts
    platePortions = find( var2d(:,:,iz)<LITHO.Tisovalue );
    platePortions = setdiff(platePortions,mantlePortionsCombined); %extract platePortions that are haven't been in mantlePortions yet
    
    if isempty(platePortions)
        check = true;
        break %no lithosphere anymore present at this depth: STOP CHECKING.
    end
    thicknessFull(platePortions)    = z2d(1,1,iz);
    if getGridPoints; gridpointsFull(platePortions) = iz; end %fill with current iz
end
if ~check; warning(['Plate thickness could not be found! - mantle temperature does not exceed critical T = ',num2str(LITHO.Tisovalue),'.']); end
%make sure it has a minimum depth (now replaced above...)
%     thicknessFull(thicknessFull<=LITHO.zmin) = LITHO.zmin;

%remove some values outside criticals
thickness                         	= thicknessFull;
thickness(thicknessFull>LITHO.zmax)	= NaN;
if getGridPoints
    gridpoints                              = gridpointsFull;
    gridpoints(thicknessFull>LITHO.zmax)    = NaN;
end

if LITHO.plot || logical(0)
    figure(22)
    plot(thickness) %2-D
    hold on
    plot(thicknessFull)
    %surf(thickness) %3-D
    xlabel('x')
    ylabel('Depth')
    %plot(gridpoints(:,:))
end

%% Output
LITHO.thickness             = thickness;            %with NaN values
LITHO.thicknessFull         = thicknessFull;        %thickness everywhere
if getGridPoints
    LITHO.gridpoints        = gridpoints;
    LITHO.gridpointsFull    = gridpointsFull;
end

end
