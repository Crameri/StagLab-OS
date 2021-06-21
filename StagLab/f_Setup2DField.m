
%%                                          SETUP 2-D PARAMETER FIELDS 2.2
%
%                                                Fabio Crameri, 22.11.2018

function [FIELD,PLOT,SAVE] = f_Setup2DField(FIELD,GRID,SWITCH,SETUP,PLOT,FILE,SAVE)

%% DEFAULTS
FIELD2D.saveData        = logical(0);

%% PARAMETER FIELDS
if strcmp(FIELD.name,'Geoid')
    Tfieldname          = 'Geoid';
    
elseif strcmp(FIELD.name,'Heat flux')
    Tfieldname          = 'Heat flux';
    
elseif strcmp(FIELD.name,'Surface age')
    Tfieldname          = 'Surface age';

elseif strcmp(FIELD.name,'Crustal thickness')
    Tfieldname          = 'Crustal thickness';
    
else
    error('fc: unknown Tfieldname')
end

%% READ 2-D FIELD
DATA.Task                   = 'ImportFieldData';
DATA.Field2Import           = Tfieldname;
DATA.FieldAbbreviation      = 'VAR';
DATA.StopExecutionIfNotFound = false;
[DATA,PLOT] = f_Import(DATA,FILE,GRID,PLOT,SWITCH);
if strcmp(GRID.Type,'yinyang')
    VAR_3D = PLOT.VAR_3Dyin; VAR_3Dyang = PLOT.VAR_3Dyang;
else %all other grid types
    VAR_3D = PLOT.VAR_3D;
end
if strcmp(GRID.Dim,'2-D')
    field2d(:,:)   	= VAR_3D(:,1,:);
elseif strcmp(GRID.Dim,'3-D')
    field2d(:,:,:)	= VAR_3D(:,:,:);
end

%% DIAGNOSE 2-D FIELD
numFields           = size(field2d,3);

%% ADJUST 2-D FIELD
%adjust fields according to StagYY output
if strcmp(FIELD.name,'Heat flux') || strcmp(FIELD.name,'Geoid')
    %flip cmb and surface data
    field2d         = fliplr(field2d);
end
%remove obsolete dataset
if strcmp(FIELD.name,'Geoid')
    if ~PLOT.GeoidCMB
        field2d(:,2)	= [];
    elseif ~PLOT.GeoidSurf
        field2d(:,1)	= [];
    end
end

%% PLOTTING OPTIONS
%e.g., SWITCH.plotDifference

% if SWITCH.plotDifference && strcmp(FIELD2D.difference,'comp')
%     if iloopComp==1
%         PLOT.field2d_B = field2d;
%         field2d     = PLOT.field2d_A;
%     else
%         field2d     = PLOT.field2d_B;
%     end
% end

%% SMOOTHING
field_smooth_2      = zeros(size(field2d));
field_smooth_1      = field2d;
if strcmp(GRID.Dim,'2-D')
    %------------------------
    % SMOOTHING 2-D
    %------------------------
    smoothingMethod = 3;
    if smoothingMethod==1 %DOESN'T WORK FOR 2-D YET! FIELD2D BECOMES SMALL!?
        filterSize = 5; %has to be an odd number!
        %F = [.05 .1 .05; .1 .4 .1; .05 .1 .05]; %filter
        %F = [1/9 1/9 1/9; 1/9 1/9 1/9; 1/9 1/9 1/9]; %like 2-D filter - with one loop!
        F = fspecial('gaussian',filterSize); %filterSize is the size of matrix
        F = F(:,ceil(end/2)); %convert for 2D
        
        %add ghost points to prevent side effects from smoothing
        dummy = field_smooth_1; %surf and cmb
        for ii=1:filterSize
            dummy = [dummy(1,:); dummy; dummy(end,:)];
        end
        field_smooth_2 = zeros(size(dummy));
        %smoothing
        for jj=1:17
            field_smooth_2(:,1) = conv(dummy(:,1),F,'same');
            if numFields>1
                field_smooth_2(:,2) = conv(dummy(:,2),F,'same');
            end
            dummy = field_smooth_2;
        end
        %remove ghostpoints
        field_smooth_2 = dummy(filterSize+1:end-filterSize,:);
        
    elseif smoothingMethod==2
        for jj=1:2
            field_smooth_2(:,1) = smooth(field_smooth_1(:,1),5,'loess');
            if numFields>1
                field_smooth_2(:,2) = smooth(field_smooth_1(:,2),5,'loess');
            end
            field_smooth_1 = field_smooth_2;
        end
    elseif smoothingMethod==3
        for jj=1:2
            field_smooth_2(:,1) = smooth(field_smooth_1(:,1),5,'moving');
            if numFields>1
                field_smooth_2(:,2) = smooth(field_smooth_1(:,2),5,'moving');
            end
            field_smooth_1 = field_smooth_2;
        end
    end
    field_smooth = field_smooth_2;

    
elseif strcmp(GRID.Dim,'3-D')
    if numFields==1
        field_PLOT = field2d(:,:,1);
    else
        if SWITCH.Verbose; warning('Only one field plot available at once, and default is the surface field - Check here!'); end
        field_PLOT = field2d(:,:,end); %surf field
%         field_PLOT = field2d(:,:,1); %cmb field
    end
    if SWITCH.plotFIELDsmooth
        %------------------------
        % SMOOTHING 3-D
        %------------------------
        smoothingMethod = 1;
        if smoothingMethod==1
            filterSize = 5;
            %F = [.05 .1 .05; .1 .4 .1; .05 .1 .05]; %filter
            %F = [1/9 1/9 1/9; 1/9 1/9 1/9; 1/9 1/9 1/9]; %like 2-D filter - with one loop!
            F = fspecial('gaussian',filterSize); %filterSize is the size of matrix
            %F = fspecial('gaussian');
            
            %add ghost points to prevent side effects from smoothing
            dummy = field_PLOT;
            for ii=1:filterSize
                dummy = [dummy(1,:); dummy; dummy(end,:)];
                dummy = [dummy(:,1), dummy, dummy(:,end)];
            end
            %smoothing
            for jj=1:17
                field_PLOTsmooth = conv2(dummy,F,'same');
                dummy = field_PLOTsmooth;
            end
            %remove ghostpoints
            field_PLOTsmooth = dummy(filterSize+1:end-filterSize,filterSize+1:end-filterSize);
            
        elseif smoothingMethod==2
            %according to smooth2a
            Nr = 3; %rows
            Nc = 3; %colums
            % You end up replacing element "i" by the mean of a (2*Nr+1)-by-
            % (2*Nc+1) rectangle centered on element "i".
            [row,col] = size(field_PLOT);
            eL = spdiags(ones(row,2*Nr+1),(-Nr:Nr),row,row);
            eR = spdiags(ones(col,2*Nc+1),(-Nc:Nc),col,col);
            % Setting all "NaN" elements of "matrixIn" to zero so that these will not
            % affect the summation.  (If this isn't done, any sum that includes a NaN
            % will also become NaN.)
            A = isnan(field_PLOT);
            field_PLOT(A) = 0;
            % For each element, we have to count how many non-NaN elements went into
            % the sums.  This is so we can divide by that number to get a mean.  We use
            % the same matrices to do this (ie, "eL" and "eR").
            nrmlize = eL*(~A)*eR;
            nrmlize(A) = NaN;
            % Actually taking the mean.
            field_PLOTsmooth = eL*field_PLOT*eR;
            field_PLOTsmooth = field_PLOTsmooth./nrmlize;
        end
    end
end

%% PROBLEM CHECKS/PREVENTION
if ~exist('field_PLOTsmooth','var'); field_PLOTsmooth = NaN; end

%% DIMENSIONALISATION
if strcmp(GRID.Dim,'2-D')
    field_PLOT      = field2d       *FIELD.varScale;
    fieldSmooth_PLOT= field_smooth  *FIELD.varScale;
else
    field_PLOT      = field_PLOT  	*FIELD.varScale;
    field_PLOTsmooth= field_PLOTsmooth *FIELD.varScale;
end

if strcmp(GRID.Dim,'2-D')
    %% SAVING FIELD ------------------------
    if FIELD2D.saveData
        error('not adjusted yet!');
        SAVE.Directory              = FILE.directory;
        SAVE.DataName               = [FILE.name,'_heatflux',num2str(FILE.number)];
        SAVE.data                   = [x_PLOT, field_PLOT, fieldSmooth_PLOT]; %[x, field, field_smooth]
        SAVE.dat                    = logical(0);
        SAVE.txt                    = logical(0);
        SAVE.mat                    = logical(1);
        SAVE.write2existing         = logical(0);
        [SAVE.overwriteAll] = f_saveData( SAVE );
    end
    
    %% PLOTTING FIELD ----------------------------------------------
    if SWITCH.plotFIELDtrue; i_start=1; else; i_start=2; end
    if SWITCH.plotFIELDsmooth; i_end=2; else; i_end=1; end
    
    %plot field
    for i_field=i_start:i_end
        if i_field==1 %plot true field
            field2plot = field_PLOT;
        elseif i_field==2 %plot smooth field
            field2plot = fieldSmooth_PLOT;
        end
    end
    
    
elseif strcmp(GRID.Dim,'3-D')
    %% 3-D
    if SWITCH.plotFIELDtrue && SWITCH.plotFIELDsmooth
        n_start     = 1;
        n_end       = 2;
    elseif SWITCH.plotFIELDtrue
        n_start     = 1;
        n_end       = 1;
    elseif SWITCH.plotFIELDsmooth
        n_start     = 2;
        n_end       = 2;
    end
    for i_field=n_start:n_end
        if i_field==1 %plot field
            field2plot = field_PLOT;
        elseif i_field==2 %plot field smooth
            field2plot = field_PLOTsmooth;
        end
    end %3-D
    
end

%% ADJUSTMENTS FOR FUNCTION OUTPUT
clearvars field2d

%% OUTPUT
FIELD.field2d        	= field2plot;



