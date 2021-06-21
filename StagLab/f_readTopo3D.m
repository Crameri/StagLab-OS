
%%                                                READ TOPOGRAPHY DATA 2.0
%
%                                                Fabio Crameri, 20.02.2014

%% NOTES
% THIS IS CURRENTLY ALSO USED FOR OTHER ROUTINES THAN SL_FieldPlot!

% check Sergei's suggestion for the smoothing

function [TOPO] = f_readTopo3D(TOPO,SWITCH,GRID,FILE)

%% Read topo information
[~,X_3D,Y_3D,Z_3D,TOPO_3D,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Topography',SWITCH);

%% IF NOT RUN VIA SL_FieldPlot
if ~isfield(SWITCH,'dimensionalInput')
    %check if output is already dimensional (if D_dim>2) !!!
    if max(Z_3D(:))>0.5 && max(Z_3D(:))<1.5; SWITCH.DimensionalInput=false; ...
    else; SWITCH.DimensionalInput=true; SWITCH.DimensionalMode=false;
    end
    nx = size(TOPO_3D,1);
    ny = size(TOPO_3D,2);
    nz = size(TOPO_3D,3);
    if (nx>1) && (ny>1) && (nz>1); GRID.Dim='3-D'; ...
    else; GRID.Dim='2-D'; error('Grid is 2-D!!!!');
    end %check dimension
    GRID.x2dp(:,:) = X_3D(:,1,1);
    GRID.y2dp(:,:) = Y_3D(1,:,1)';
end

%% DIMENSIONALISE
xi_plot         = GRID.y2dp; %[plotting dimension]
yi_plot         = GRID.x2dp; %[plotting dimension]

if isfield(GRID,'dimFactor') %used by SL_FieldPlot
    if SWITCH.DimensionalMode; varScaleTopo = GRID.dimFactor; else; varScaleTopo = 1; end
else
    if SWITCH.DimensionalMode && ~SWITCH.DimensionalInput; varScaleTopo = TOPO.varScale; else; varScaleTopo = 1; end
end
topoP           = TOPO_3D(:,:,2)*varScaleTopo; %surf topo [m]
topoCmbP        = TOPO_3D(:,:,1)*varScaleTopo; %cmb topo [m]

%% OUTPUT TOPOGRAPHY
smoothTopo      = logical(1);
if smoothTopo
    n_start     = 2;
    n_end       = 2;
    topoSmooth  = [topoCmbP, topoP];
else
    n_start     = 1;
    n_end       = 1;
end
for i_field=n_start:n_end
    if i_field==2 %topo smooth
        %% SMOOTHING
        if strcmp(GRID.Dim,'2-D')
            smoothingMethod = 3;
            if smoothingMethod==1 %DOESN'T WORK FOR 2-D YET!!!!!!!!! TOPO BECOMES SMALL!?
                filterSize = 5; %has to be an odd number!!!!
                %F = [.05 .1 .05; .1 .4 .1; .05 .1 .05]; %filter
                %F = [1/9 1/9 1/9; 1/9 1/9 1/9; 1/9 1/9 1/9]; %like 2-D filter - with one loop!
                F = fspecial('gaussian',filterSize); %filterSize is the size of matrix
                F = F(:,ceil(end/2)); %convert for 2D
                
                %add ghost points to prevent side effects from smoothing
                dummy = topoSmooth; %surf and cmb
                for ii=1:filterSize
                    dummy = [dummy(1,:); dummy; dummy(end,:)];
                end
                topoPsmooth = zeros(size(dummy));
                %smoothing
                for jj=1:17
                    topoPsmooth(:,1) = conv(dummy(:,1),F,'same');
                    topoPsmooth(:,2) = conv(dummy(:,2),F,'same');
                    dummy = topoPsmooth;
                end
                %remove ghostpoints
                topoPsmooth = dummy(filterSize+1:end-filterSize,:);
                
            elseif smoothingMethod==2
                for jj=1:2
                    topoPsmooth(:,1) = smooth(topoSmooth(:,1),5,'loess');
                    topoPsmooth(:,2) = smooth(topoSmooth(:,2),5,'loess');
                    topoSmooth = topoPsmooth;
                end
            elseif smoothingMethod==3
                for jj=1:2
                    topoPsmooth(:,1) = smooth(topoSmooth(:,1),5,'moving');
                    topoPsmooth(:,2) = smooth(topoSmooth(:,2),5,'moving');
                    topoSmooth = topoPsmooth;
                end
            end
            topoP           = topoPsmooth(:,2);
            topoCmbP       = topoPsmooth(:,1);
            
        elseif strcmp(GRID.Dim,'3-D')
            for itopo=1:size(topoSmooth,3)
                topoSmoothS(:,:)     = topoSmooth(:,:,itopo);
                smoothingMethod = 1;
                if smoothingMethod==1
                    filterSize = 5;
                    %F = [.05 .1 .05; .1 .4 .1; .05 .1 .05]; %filter
                    %F = [1/9 1/9 1/9; 1/9 1/9 1/9; 1/9 1/9 1/9]; %like 2-D filter - with one loop!
                    F = fspecial('gaussian',filterSize); %filterSize is the size of matrix
                    %F = fspecial('gaussian');
                    
                    %add ghost points to prevent side effects from smoothing
                    dummy = topoSmoothS;
                    for ii=1:filterSize
                        dummy = [dummy(1,:); dummy; dummy(end,:)];
                        dummy = [dummy(:,1), dummy, dummy(:,end)];
                    end
                    %smoothing
                    for jj=1:17
                        topoPsmooth = conv2(dummy,F,'same');
                        dummy = topoPsmooth;
                    end
                    %remove ghostpoints
                    topoPsmooth = dummy(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
                    
                elseif smoothingMethod==2
                    %according to smooth2a
                    Nr = 3; %rows
                    Nc = 3; %colums
                    % You end up replacing element "i" by the mean of a (2*Nr+1)-by-
                    % (2*Nc+1) rectangle centered on element "i".
                    [row,col] = size(topoSmoothS);
                    eL = spdiags(ones(row,2*Nr+1),(-Nr:Nr),row,row);
                    eR = spdiags(ones(col,2*Nc+1),(-Nc:Nc),col,col);
                    % Setting all "NaN" elements of "matrixIn" to zero so that these will not
                    % affect the summation.  (If this isn't done, any sum that includes a NaN
                    % will also become NaN.)
                    A = isnan(topoSmoothS);
                    topoSmoothS(A) = 0;
                    % For each element, we have to count how many non-NaN elements went into
                    % the sums.  This is so we can divide by that number to get a mean.  We use
                    % the same matrices to do this (ie, "eL" and "eR").
                    nrmlize = eL*(~A)*eR;
                    nrmlize(A) = NaN;
                    % Actually taking the mean.
                    topoPsmooth = eL*topoSmoothS*eR;
                    topoPsmooth = topoPsmooth./nrmlize;
                end
                if itopo==1
                    topoP           = topoPsmooth;
                else
                    topoCmbP       = topoPsmooth;
                end
            end
        end
    end
end
% surf(xi_plot,yi_plot,zi_plot,'EdgeColor','none')
% contour(xi,yi,zi)

%% OUTPUT
TOPO.x          = xi_plot;      % (non-dim or dim [plotting dimension])
TOPO.y          = yi_plot;      % (non-dim or dim [plotting dimension])
TOPO.z_surf 	= topoP;    %surf topography (non-dim or dim [m])
TOPO.z_cmb      = topoCmbP;%cmb topography (non-dim or dim [m])




