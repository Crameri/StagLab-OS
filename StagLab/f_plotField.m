
%%                                             PLOT PARAMETER FIELDS 10.172
%    . contains f_plotCircle
%    . contains f_cartesian2spherical
%    . calls f_DesignColourmap
%    . calls hatchfill2
%    . calls f_DesignVaria
%
% Grid is flipped (x and y) in 3-D plots.
%                                                Fabio Crameri, 25.02.2020

function [PLOT,SAVE] = f_plotField(VAR,FIELD,FIELD_M,GRID,SETUP,SWITCH,PLOT,TOPO,FILE,STYLE,SAVE,PLATE,MANTLE)

%% PROBLEM CHECKS
if strcmp(FIELD.name,'Upwelling and downwelling') && MANTLE.DiagnosticsFailed
    %% CREATE AN EMPTY PLOT
    [PLOT,SAVE] = f_EmptyPlot(SWITCH,PLOT,GRID,FIELD,STYLE,SAVE);
    return
    
end

%% FUNCTION STARTUP
if SWITCH.waitbar
    disp('   ...processing')
    PLOT.wb = waitbar(0.4,PLOT.wb,'processing...');
end

PLOT.nrSubplot = PLOT.nrSubplot+1; %count subplot

if isfield(PLOT,'quiver3Dplates') && PLOT.quiver3Dplates
    plate_velocity = logical(1);  %for 3-D velocity arrows
else
    plate_velocity = logical(0);  %for 3-D velocity arrows
end

%% VARIABLE SETUP
if strcmp(GRID.Type,'yinyang')
    VARyin  = VAR.var3d_yin; 	VARyang = VAR.var3d_yang;
else %everything else
    var2d   = VAR.var2d;
end

x2d_nd  = GRID.x2d_nd;  y2d_nd  = GRID.y2d_nd;  z2d_nd  = GRID.z2d_nd;  %[nd]
x2dp    = GRID.x2dp;    y2dp    = GRID.y2dp;    z2dp    = GRID.z2dp;    %already reversed if Cartesian, dim
nx      = GRID.nx;    	ny      = GRID.ny;   	nz      = GRID.nz;

CBAR.constant       = SWITCH.Colorbar(1,2);

grid_T              = PLOT.gridAddition(1,1);
grid_eta          	= PLOT.gridAddition(1,2);
grid_str           	= PLOT.gridAddition(1,3);
grid_edot        	= PLOT.gridAddition(1,4);

quiver_plot         = PLOT.quiverM(1,1);
quiver_T            = PLOT.quiverM(1,2);
quiver_eta          = PLOT.quiverM(1,3);
quiver_str          = PLOT.quiverM(1,4);
quiver_numDel     	= PLOT.quiverM(1,5);
scale               = PLOT.quiverM(1,6);
quiver_psi          = PLOT.quiverM(1,7);
quiver_rho          = PLOT.quiverM(1,8);
quiver_nstr         = PLOT.quiverM(1,9);
quiver_edot         = PLOT.quiverM(1,10);
quiver_v            = PLOT.quiverM(1,11);
quiver_vx        	= PLOT.quiverM(1,12);
quiver_vr        	= PLOT.quiverM(1,13);

princStressDir_plot	= PLOT.princStressM(1,1);
princStress_numDel	= PLOT.princStressM(1,2);
princStress_scale 	= PLOT.princStressM(1,3);
prstrDir_T        	= PLOT.princStressM(1,4);
prstrDir_eta     	= PLOT.princStressM(1,5);
prstrDir_str     	= PLOT.princStressM(1,6);
prstrDir_psi      	= PLOT.princStressM(1,7);
prstrDir_rho      	= PLOT.princStressM(1,8);
prstrDir_nstr    	= PLOT.princStressM(1,9);
prstrDir_edot      	= PLOT.princStressM(1,10);

% streamfun_plot  = PLOT.streamfun(1,1);
streamfun_T         = PLOT.streamfun(1,2);
streamfun_eta       = PLOT.streamfun(1,3);
streamfun_str       = PLOT.streamfun(1,4);
stream_numContours  = PLOT.streamfun(1,5);
streamfun_rho       = PLOT.streamfun(1,6);
streamfun_v         = PLOT.streamfun(1,7);
streamfun_vx     	= PLOT.streamfun(1,8);
streamfun_vr     	= PLOT.streamfun(1,9);
streamfun_edot      = PLOT.streamfun(1,10);
streamfun_udw       = PLOT.streamfun(1,11);
stream_data         = PLOT.streamdata;

streamlinePlot      = PLOT.streamline(1,1);
streamline_T        = PLOT.streamline(1,2);
streamline_eta      = PLOT.streamline(1,3);
streamline_str      = PLOT.streamline(1,4);
streamline_v        = PLOT.streamline(1,5);
streamline_vh       = PLOT.streamline(1,6);
streamline_vr       = PLOT.streamline(1,7);
startline           = PLOT.streamStartline; sx2d=startline(1,:); sy2d=startline(2,:); sz2d=startline(3,:);

defmech_T           = PLOT.defmech(1,1);
defmech_eta         = PLOT.defmech(1,2);
defmech_str         = PLOT.defmech(1,3);
defmech_edot        = PLOT.defmech(1,4);
defmech_contours    = PLOT.defmech(1,5:end);

lithoThickness_isoVal= PLOT.lithoThickness(1,1);
lithoThickness_zmax = PLOT.lithoThickness(1,2);
lithoThickness_T    = PLOT.lithoThickness(1,3);
lithoThickness_eta  = PLOT.lithoThickness(1,4);
lithoThickness_field= PLOT.LTfieldContour;
lithoThickness_nstr = PLOT.lithoThickness(1,5);

fieldContour_isoVal = PLOT.fieldContour(1,1);
fieldContour_C      = PLOT.fieldContour(1,2);
fieldContour_T      = PLOT.fieldContour(1,3);
fieldContour_eta    = PLOT.fieldContour(1,4);
fieldContour_nstr   = PLOT.fieldContour(1,5);
fieldContour_v      = PLOT.fieldContour(1,6);
fieldContour_vx     = PLOT.fieldContour(1,7);
fieldContour_vr     = PLOT.fieldContour(1,8);
fieldContour_str    = PLOT.fieldContour(1,9);
fieldContour_edot   = PLOT.fieldContour(1,10);
fieldContour_vdiss  = PLOT.fieldContour(1,11);
fieldContour_stream = PLOT.fieldContour(1,12);
fieldContour_udw    = PLOT.fieldContour(1,13);

% plume_plot          = PLOT.plume(1,1);
plume_T             = PLOT.plume(1,2);
plume_eta           = PLOT.plume(1,3);
plume_str           = PLOT.plume(1,4);
plume_rho           = PLOT.plume(1,5);
plume_v             = PLOT.plume(1,6);
plume_edot          = PLOT.plume(1,7);

%% VARIABLE ADJUSTMENTS
if strcmp(SWITCH.quiverNumArrows,'auto')
    if strcmp(GRID.Type,'Cartesian')
        quiver_numDel  	= max(0,round(size(GRID.Z_3Dp,3)/14));
    else
        quiver_numDel  	= max(0,round(size(GRID.Z_3Dp,3)/10));
    end
end
if strcmp(SWITCH.princStressNumArrows,'auto')
    if strcmp(GRID.Type,'Cartesian')
        princStress_numDel 	= max(0,round(size(GRID.Z_3Dp,3)/30));
    else
        princStress_numDel 	= max(0,round(size(GRID.Z_3Dp,3)/22));
    end
end

%% ERROR CHECKS
if isfield(VAR,'var2d') && isfield(VAR,'var3d_yin'); error('It is not clear wheter to plot var2d or var2d_yin: check here!'); end

%% FIGURE LAYOUT
%figure...
figure(1)

%subplot layout...
SPOS.task       = 'createSubplot';
[SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE);
%SP.number;    %for possible reversed layout - while PLOT.nrSubplot is still setup the old way!

%% SETUP ANNOTATION STRINGS
DESVARIA.Task 	= 'create annotation strings';
[PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);

%% CURRENT AXIS
AXcurrent       = gca; %save current axes handle

%% DIMENSIONALISATION
if strcmp(GRID.Type,'yinyang')
    VARyin      = VARyin .*FIELD.varScale;
    VARyang     = VARyang .*FIELD.varScale;
    if FIELD.logarithmic
        VARyin 	= log10(VARyin);
        VARyang	= log10(VARyang);
    end
else %everything else
    var2d       = var2d .*FIELD.varScale;
    if FIELD.logarithmic
        var2d   = log10(var2d);
    end
end

if strcmp(GRID.Type,'Cartesian')
    %% CARTESIAN 3-D PLOT ******************************************************************
    if strcmp(GRID.Dim,'3-D')   % 3-D
        AXorig  = gca;
        if ~isfield(PLOT,'isoVal'); PLOT.isoVal = [0.45 0.2 0.7 0.1 0.9]; end
        if CBAR.constant
            maxVal      = FIELD.cColorbarMAX;
            minVal      = FIELD.cColorbarMIN;
        else
            maxVal      = max(var2d(:));
            minVal      = min(var2d(:));
        end
        rangeVal        = maxVal-minVal;
        PLOT.isoVal_dim         = minVal+PLOT.isoVal.*rangeVal;
        PLOT.isoValSpecial_dim  = minVal+PLOT.isoValSpecial.*rangeVal;
        if strcmp(FIELD.name,'Upwelling and downwelling')
%             PLOT.isoVal_dim         = [1.67,0.83,-0.83,-1.67];
            PLOT.isoVal             = [1.95,0.95,-0.95,-1.95];
            PLOT.isoVal_dim         = [1.95,0.95,-0.95,-1.95];
            PLOT.isoValSpecial_dim  = [];
        end
        clearvars maxVal minVal rangeVal
        
        if SWITCH.ReduceSize
            wb = waitbar(0,'Slimming Matrix...');
            %reduce size of matrix
            maxGridpoints   = 300;   %<--------------------------=================
            nxEff           = nx;
            nyEff           = ny;
            nzEff           = nz;
            while (nxEff+nyEff+nzEff)>maxGridpoints
                waitbar(maxGridpoints/(nxEff+nyEff+nzEff),wb)
                for i=(nxEff-2):-2:2
                    x2dp(i,:)       = [];
                    var2d(i,:,:)    = [];
                end
                for j=(nyEff-2):-2:2
                    y2dp(:,j)       = [];
                    var2d(:,j,:)    = [];
                end
                for k=(nzEff-2):-2:2
                    z2dp(:,:,k) 	= [];
                    var2d(:,:,k)	= [];
                end
                x2dp(4,:)    	= [];
                var2d(4,:,:)    = [];
                y2dp(:,4)       = [];
                var2d(:,4,:)    = [];
                z2dp(:,:,4)   	= [];
                var2d(:,:,4)    = [];
                
                nxEff           = size(x2dp,1);
                nyEff           = size(y2dp,2);
                nzEff           = size(z2dp,3);
            end
            disp(['Size of plotted matrix decreased to ',num2str(size(var2d))])
            clearvars nxEff nyEff nzEff maxGridpoints
            close(wb)
        end
        
        % set up grid box...
        if SWITCH.AxesEqual; axis equal; end
        %it seems that isosurface needs the order: (y,x,z,var)...
        axis([min(y2dp) max(y2dp) min(x2dp) max(x2dp) min(z2dp) max(z2dp)])
        box on
        grid on
        %         L(1) = light; %original lighting
        %---------------------
        set(gca,'ZDir','reverse')
        %---------------------
        %set camera position
        if PLOT.cameraPosBotView
            campos(PLOT.cameraPosition.*[1 1 -1]) %flipped upside down
        else
            campos(PLOT.cameraPosition) %standard
        end
        %set perspective mode
        camproj(PLOT.camProjection)
        
        % plotting...
        %isosurface
        if isfield(PLOT,'isoVal') && ~isempty(PLOT.isoVal)
            for iplot_iso=1:size(PLOT.isoVal,2)
                %it seems that isosurface needs the order: (y,x,z,var)...
                %isosurface(x2d_plot,y2d_plot,z2d_plot,var2d_plot,PLOT.isoVal_dim(1,iplot_iso),'noshare');
                isosurface(y2dp,x2dp,z2dp,var2d,PLOT.isoVal_dim(1,iplot_iso),'noshare');
                if iplot_iso<(size(PLOT.isoVal,2)-PLOT.opacity(1,2)+1); alpha(PLOT.opacity(1,1)); end %make last isosurface opaque
                %                 camlight(PLOT.camlight); lighting(PLOT.lighting); %needed for every surface - I think
                hold on
            end
        end
        %special isosurface
        if isfield(PLOT,'isoVal_special') && ~isempty(PLOT.isoValSpecial)
            %it seems that isosurface needs the order: (y,x,z,var)...
            %             hspecial = patch(isosurface(x2d_plot,y2d_plot,z2d_plot,var2d_plot,PLOT.isoValSpecial_dim(1,1),'noshare'));
            hspecial = patch(isosurface(y2dp,x2dp,z2dp,var2d,PLOT.isoValSpecial_dim(1,1),'noshare'));
            isonormals(y2dp,x2dp,z2dp,var2d,hspecial);
            set(hspecial,'FaceAlpha',PLOT.opacitySpecial(1,1),'FaceColor',PLOT.isoColorSpecial,'EdgeColor','none');
            %             drawnow
            %             camlight(PLOT.camlight); lighting(PLOT.lighting); %needed for every surface - I think
            hold on
        end
        %saving field to file
        if FIELD.save || logical(0) %save current field
            if ~isfield(SAVE,'Directory'); SAVE.Directory = FILE.directory; end
            f_saveField(x2dp,y2dp,z2dp,var2d,FIELD,FILE,SAVE);
        end
        
        %LIGHTING
        L(1) = light; %original lighting
        set(L(1),'pos',[-1 0 -1]); %adjust lighting to reverse axis
        camlight(PLOT.camlight); lighting(PLOT.lighting); %needed for every surface - I think
        
        %warning('transparency and lighting do only work with the opengl renderer, which is a bit buggy at times... add error checks to prevent using both, ie. possible problems.')
        
        %% slices
        if PLOT.Slices
            % if PLOT.Streamfunction; PLOT.LineColor=streamColor; else;  PLOT.LineColor='none'; end
            
            %adjust the grid for slicing...
            [y2d_plot2,x2d_dummy,z2d_dummy] = meshgrid(y2dp,x2dp,z2dp);
            %PLOT MAIN FIELD SLICE...
            hold on
            sp1 = slice(y2d_plot2,x2d_dummy,z2d_dummy,var2d,PLOT.sliceY,PLOT.sliceX,PLOT.sliceZ);
            set(sp1,'edgecolor','none');
            set(sp1,'FaceColor','interp')
            alpha(sp1,PLOT.SliceOpacity);
            
            %saving slice field to file
            if logical(0)  %FIELD.save || logical(0) %save current field for a 2-D slice
                var2d_save = zeros(size(var2d,1),size(var2d,3));
                %---------------------------------------------
                var2d_save(:,:) = var2d(:,end/2,:); %<<<<adjust here which slice
                %---------------------------------------------
                z2d_save = zeros(size(z2dp,3),1); z2d_save(:,1) = z2dp(1,1,:);
                f_saveField( x2dp,1,z2d_save,var2d_save,FIELD,FILE,SAVE );
                %and time data
                if PLOT.loopField==1
                    SAVE.Directory              = ['~,',filesep,'Desktop',filesep,FILE.name,filesep];
                    SAVE.DataName               = ['time_',num2str(FILE.number)];
                    SAVE.data                   = [PLOT.time, PLOT.time_dim*SETUP.secyear,PLOT.time_dim]; %[nd, seconds, years]
                    SAVE.dat                    = logical(0);
                    SAVE.txt                    = logical(0);
                    SAVE.mat                    = logical(1);
                    SAVE.write2existing         = logical(0);
                    
                    [SAVE.overwriteAll] = f_saveData( SAVE );
                end
            end
            %PLOT SPECIAL SLICE...
            if PLOT.SliceSpecial
                hold on
                sp2 = slice(y2d_plot2,x2d_dummy,z2d_dummy,var2d,PLOT.sliceYS,PLOT.sliceXS,PLOT.sliceZS);
                set(sp2,'edgecolor','none');
                set(sp2,'FaceColor','interp')
                alpha(sp2,PLOT.SliceOpacity);
            end
            %ADD QUIVER SLICES...
            if PLOT.sliceAddQuiver  %adjust this to work also with multiple fields and only plot for certain fields
                if isfield(PLOT,'VX_3D') && isfield(PLOT,'VY_3D') && isfield(PLOT,'VZ_3D')
                    VX_3D = PLOT.VX_3D; VY_3D = PLOT.VY_3D; VZ_3D = PLOT.VZ_3D;
                else
                    % Read pressure & velocity information
                    [~,~,~,~,VX_3D,VY_3D,VZ_3D,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Velocity',SWITCH);
                    if ischar(VX_3D); error(VX_3D); end
                    if size(VX_3D,1)==1 %exchange x and y
                        dummy_3D = zeros(size(VY_3D,2),size(VY_3D,1),size(VY_3D,3));
                        dummy_3Dx = zeros(size(VY_3D,2),size(VY_3D,1),size(VY_3D,3));
                        dummy_3D(:,1,:)     = VZ_3D(1,:,:); VZ_3D = dummy_3D;
                        dummy_3Dx(:,1,:)    = VX_3D(1,:,:);
                        dummy_3D(:,1,:)     = VY_3D(1,:,:); VY_3D = dummy_3Dx; VX_3D = dummy_3D;
                    end
                    VZ_3D   = -VZ_3D; % for reverse z-axis  (1.)
                end
                X_3D_V  = zeros(size(GRID.X_3D))*NaN;
                Y_3D_V  = X_3D_V;
                Z_3D_V  = X_3D_V;
                VX_3D_V = X_3D_V;
                VY_3D_V = X_3D_V;
                VZ_3D_V = X_3D_V;
                sliceXgpV = [];
                sliceYgpV = [];
                sliceZgpV = [];
                % Extract slices close to desired values & reduce data size
                sNrWrite = PLOT.sliceNrWrite;
                %find closest x-value to desired slice
                for ii=1:size(PLOT.sliceXV,2)
                    tmp         = abs((GRID.X_3Dp(:,1,1))-PLOT.sliceXV(1,ii));
                    [minx idx]  = min(tmp); %index of closest value
                    sliceXgpV   = [sliceXgpV idx];
                end
                if ~isempty(sliceXgpV)
                    X_3D_V(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end) = GRID.X_3Dp(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end);
                    Y_3D_V(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end) = GRID.Y_3Dp(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end);
                    Z_3D_V(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end) = GRID.Z_3Dp(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end);
                    VX_3D_V(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end) = VX_3D(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end);
                    VY_3D_V(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end) = VY_3D(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end);
                    VZ_3D_V(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end) = VZ_3D(sliceXgpV,1:sNrWrite:end,1:sNrWrite:end);
                end
                %find closest y-value to desired slice
                for ii=1:size(PLOT.sliceYV,2)
                    tmp         = abs((GRID.Y_3Dp(1,:,1))-PLOT.sliceYV(1,ii));
                    [minx idx]  = min(tmp); %index of closest value
                    sliceYgpV   = [sliceYgpV idx];
                end
                if ~isempty(sliceYgpV)
                    X_3D_V(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end) = GRID.X_3Dp(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end);
                    Y_3D_V(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end) = GRID.Y_3Dp(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end);
                    Z_3D_V(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end) = GRID.Z_3Dp(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end);
                    VX_3D_V(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end) = VX_3D(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end);
                    VY_3D_V(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end) = VY_3D(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end);
                    VZ_3D_V(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end) = VZ_3D(1:sNrWrite:end,sliceYgpV,1:sNrWrite:end);
                end
                %find closest z-value to desired slice
                for ii=1:size(PLOT.sliceZV,2)
                    tmp         = abs((GRID.Z_3Dp(1,1,:))-PLOT.sliceZV(1,ii));
                    [minx idx]  = min(tmp); %index of closest value
                    sliceZgpV   = [sliceZgpV idx];
                end
                if ~isempty(sliceZgpV)
                    X_3D_V(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV) = GRID.X_3Dp(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV);
                    Y_3D_V(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV) = GRID.Y_3Dp(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV);
                    Z_3D_V(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV) = GRID.Z_3Dp(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV);
                    VX_3D_V(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV) = VX_3D(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV);
                    VY_3D_V(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV) = VY_3D(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV);
                    VZ_3D_V(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV) = VZ_3D(1:sNrWrite:end,1:sNrWrite:end,sliceZgpV);
                end
                clearvars sliceXgpV sliceYgpV sliceZgpV minx idx
                % Plot quiver
                hold on
                q1 = quiver3(Y_3D_V,X_3D_V,Z_3D_V,VX_3D_V,VY_3D_V,VZ_3D_V,scale);
                warning('check if VX_3D and VY_3D need to be switched here!')
                set(q1,'Color',PLOT.quiverColor,'LineWidth',PLOT.quiverWidth)
            end
            
            % if IN.add_streamline
            %     hold on
            %     [phi2,theta2] = meshgrid(phi_V,theta_V);
            %     %             [s1,s2] = meshgrid( phi2(:,10),theta2(:,10));     %starting points
            %     s1=phi2(5,10);
            %     s2=theta2(5,10);
            %     [sx,sy]=meshgrid(s1,s2);
            %     streamline(x2d_V,y2d_V,VX_map(:,:,SP.number),VY_map(:,:,SP.number),sx,sy)
            % end
            
            %PLOT TOPOGRAPHY SLICE ON TOP
            if PLOT.SliceAddTopo
                if TOPO.cmb; startT=1; else; startT=2; end
                for iT=startT:2
                    colorbar('hide')
                    freezeColors %freeze previous colormap
                    hold on
                    %TOPO COLORMAP
                    %if self-made colorbar:
                    if exist([SWITCH.ColormapAllT{1,1},'.mat'],'file')==2
                        load([SWITCH.ColormapAllT{1,1},'.mat']); SWITCH.ColormapAllTnum{1,1}=eval( sprintf(SWITCH.ColormapAllT{1,1}) ); %load colorbar
                    end
                    % else predefined colorbars:
                    if strcmp(SWITCH.ColormapAllT{1,2},'flip')
                        colormap(SWITCH.ColormapAllTnum{1,1})
                        colormap(flipud(colormap));
                    elseif strcmp(SWITCH.ColormapAllT{1,2},'noflip')
                        colormap(SWITCH.ColormapAllTnum{1,1})
                    else
                        error('Unknown Colormap option: Flip or no flip?')
                    end
                    TOPO.directory = FILE.directory; TOPO.fname = FILE.name; TOPO.fname_number = FILE.number;
                    TOPO.varScale = SETUP.D_dim; %km
                    
                    [TOPO] = f_readTopo3D(TOPO,SWITCH,GRID,FILE); %read topography
                    
                    %adjust plotting values--
                    TOPO.z_surf         = -TOPO.z_surf;  %flip due to depth plot
                    TOPO.z_cmb          = -TOPO.z_cmb;
                    %dimensionalise ADJUST HERE IF NECESSARY (CHECK TOPO3D FUNCTION)!
                    %                     TOPO.z_surf         = TOPO.z_surf;
                    %                     TOPO.z_cmb          = TOPO.z_cmb;
                    
                    TOPO.z_surf_plot    = TOPO.z_surf*TOPO.exagFactor;  %flip due to depth plot %TOPO.z_surf+2*(TOPO.z_surf-mean(mean(TOPO.z_surf)));
                    TOPO.z_cmb_plot     = TOPO.z_cmb*TOPO.exagFactor;
                    colormap(flipud(colormap));
                    %------------------------
                    if iT==1
                        topoZ1 = TOPO.z_cmb; topoZ2 = TOPO.z_cmb+1*D_act/conversian_factor; meanT=mean(mean(topoZ2)); %plot at the bottom;
                        topoZ2_plot = TOPO.z_cmb_plot+1*D_act/conversian_factor; meanT_plot=mean(mean(topoZ2_plot)); %plot at the bottom;
                    elseif iT==2
                        topoZ1 = TOPO.z_surf; topoZ2 = topoZ1; meanT=0;
                        topoZ2_plot = TOPO.z_surf_plot; meanT_plot=0;
                    end
                    if iT==2 && ~TOPO.surf
                    else
                        %TOPO PLOTTING
                        TOPO.slice = false;
                        if TOPO.slice
                            error('this is not yet adjusted for use in SL_FieldPlot')
                            topoDUMMY3D = ones(size(X_3D));
                            for iz=1:size(X_3D,3)
                                topoDUMMY3D(:,:,iz) = topoZ1;
                            end
                            sptopo = slice(Y_3D,Z_3D,Z_3D,topoDUMMY3D,IN.slice_yT,IN.slice_xT,IN.slice_zT);
                        else
                            sptopo = surf(GRID.y2dp,GRID.x2dp,topoZ2_plot); %actual plot
                            camlight(PLOT.camlight); lighting(PLOT.lighting); %needed for every surface - I think
                            hold on
                            sptopo2 = surf(GRID.y2dp,GRID.x2dp,topoZ2);  %for manipulating colorbar; will be deleted later
                        end
                        set(sptopo,'edgecolor','none');
                        set(sptopo,'FaceColor','interp')
                        alpha(sptopo,PLOT.SliceOpacityT);
                        hold on
                    end
                    if iT==1 && TOPO.surf && TOPO.cmb
                        topoLIM1 = max(abs(min(min(TOPO.z_surf))), abs(max(max(TOPO.z_surf))));
                        topoLIM2 = max(abs(min(min(TOPO.z_cmb))), abs(max(max(TOPO.z_cmb))));
                        topoLIM = max(topoLIM1,topoLIM2);
                        if SWITCH.Verbose
                            disp(['  surface topography: max=',num2str(max(max(TOPO.z_surf)),3),', min=',num2str(min(min(TOPO.z_surf)),3)])
                            disp(['  CMB topography:     max=',num2str(max(max(TOPO.z_cmb)),3),', min=',num2str(min(min(TOPO.z_cmb)),3)])
                        end
                    elseif iT==2 && TOPO.surf && ~TOPO.cmb; topoLIM = max(abs(min(min(TOPO.z_surf))), abs(max(max(TOPO.z_surf))));
                        if SWITCH.Verbose; disp(['  surface topography: max=',num2str(max(max(TOPO.z_surf)),3),', min=',num2str(min(min(TOPO.z_surf)),3)]); end
                    elseif iT==1 && TOPO.cmb; topoLIM = max(abs(min(min(TOPO.z_cmb))), abs(max(max(TOPO.z_cmb))));
                        if SWITCH.Verbose; disp(['  CMB topography:     max=',num2str(max(max(TOPO.z_cmb)),3),', min=',num2str(min(min(TOPO.z_cmb)),3)]); end
                    end
                    %for manipulating colorbar--
                    topoLIM_plot = topoLIM*TOPO.exagFactor;
                    if exist('sptopo2','var'); set(sptopo2,'edgecolor','none'); end
                    if SWITCH.constantColorbarT; set(gca,'clim',[meanT_plot-topoLIM_plot,meanT_plot+topoLIM_plot]); end
                    freezeColors %freeze previous colormap
                    set(sptopo,'edgecolor','none');
                    colormap(flipud(colormap));
                    %---------------------------
                    if SWITCH.constantColorbarT; set(gca,'clim',[meanT-topoLIM,meanT+topoLIM]); end
                    SWITCH.positionColorbar = 3; %1: bottom  2: right  3: south-east
                    if SWITCH.positionColorbar==1 %bottom
                        cbT = colorbar('horiz');
                        set(cbT,'Position',[0.4 0.04 0.2 0.017]);
                    elseif SWITCH.positionColorbar==2 %right
                        cbT = colorbar;
                        set(cbT,'Position',[0.9 axPos(1,2) 0.017 axPos(1,4)]);
                    elseif SWITCH.positionColorbar==3 %south-east
                        cbT = colorbar('horiz');
                        %set(cbT,'Position',[0.73 0.77 0.17 0.014]); %top
                        set(cbT,'Position',[0.73 0.12 0.17 0.014]); %bottom
                    end
                    %         set( cbT, 'XDir', 'reverse' );
                    if SWITCH.positionColorbar==2 %right
                        if SWITCH.DimensionalMode || SWITCH.DimensionalInput
                            ylabel(cbT,'Topography [km]','FontWeight','bold');
                        else
                            ylabel(cbT,'Topography [nd]','FontWeight','bold');
                        end
                    else %south-east
                        if SWITCH.DimensionalMode || SWITCH.DimensionalInput
                            title(cbT,'Topography [km]','FontWeight','bold');
                        else
                            title(cbT,'Topography [nd]','FontWeight','bold');
                        end
                    end
                    if exist('sptopo2','var') && iT==1; set(sptopo2,'Visible','off'); end
                end
                if TOPO.exagFactor~=1 && SWITCH.Verbose
                    disp(['  topo exaggeration factor: ',num2str(TOPO.exagFactor)]);
                end
                delete(sptopo2);
            end
        end
        
        legend4isosurface = logical(0);
        if legend4isosurface && ~PLOT.Slices
            DESVARIA.Task    = 'make colorbar';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            colorbar(PLOT.cb,'off');
            % legend...
            isoValAll   = [PLOT.isoVal_dim,PLOT.isoValSpecial_dim];
            for il=1:size(isoValAll,2)
                lStrings{il,1} = [num2str(isoValAll(:,il),4),' ',FIELD.dim];
            end
            if strcmp(FIELD.name,'Upwelling and downwelling')
                lStrings    = { ['active',char(11014)]; ['passive',char(11014)]; ['passive',char(11015)]; ['active',char(11015)]};
            end
            legend(lStrings)
            PLOT.cb = NaN; %no colorbar
        else
            % colorbar...
            DESVARIA.Task    = 'make colorbar';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
        end
        % title...
        DESVARIA.Task    = 'make title';
        [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
        % annotations...
        xlabel(PLOT.xlabel); ylabel(PLOT.ylabel); zlabel(PLOT.zlabel);
        if SWITCH.Texting
            text(max(xlim),min(ylim),['max val: ',num2str(max(var2d(:)),2)],'HorizontalAlignment','right','VerticalAlignment','bottom','Color',STYLE.keyColor,'FontName',STYLE.AllFontName,'FontSize',8);
            text(min(xlim),min(ylim),[num2str(nx),'x',num2str(ny),'x',num2str(nz)],'HorizontalAlignment','left','VerticalAlignment','bottom','Color',STYLE.keyColor,'FontName',STYLE.AllFontName,'FontSize',8);
            disp(['   ',FIELD.name,' variation: ',num2str(min(var2d(:)),4),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(max(var2d(:)),4),' ',FIELD.dim])
        else
            % text(min(xlim),min(ylim),[num2str(nx),'x',num2str(ny),'x',num2str(nz)],'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8)
        end
        %axis background colouring
        if ~strcmp(STYLE.ColorMode,'light')
            set(gca,'color',STYLE.ColorModeBW);
        end
        
    elseif strcmp(GRID.Dim,'2-D')
        %% CARTESIAN 2-D PLOT ******************************************************************
        %reduce size of matrix
        if SWITCH.ReduceSize
            wb = waitbar(0,'Slimming Matrix...');
            maxGridpoints       = 30000;   %<--------------------------=================
            nxEff               = nx;
            nzEff               = nz;
            while (nxEff*nzEff)>maxGridpoints
                waitbar(maxGridpoints/(nxEff+nzEff),wb)
                for i=(nxEff-2):-2:2
                    x2dp(i,:)  	= [];
                    z2dp(i,:) 	= [];
                    var2d(i,:) 	= [];
                end
                for k=(nzEff-2):-2:2
                    x2dp(:,k) 	= [];
                    z2dp(:,k)  	= [];
                    var2d(:,k) 	= [];
                end
                x2dp(4,:)       = [];
                z2dp(4,:)       = [];
                var2d(4,:)      = [];
                x2dp(:,4)       = [];
                z2dp(:,4)       = [];
                var2d(:,4)      = [];
                nxEff           = size(x2dp,1);
                nzEff           = size(x2dp,2);
            end
            disp(['Size of plotted matrix decreased to ',num2str(size(var2d))])
            clearvars nxEff nyEff nzEff maxGridpoints
            close(wb)
        end
        %setup for zoom plot
        numLoopsZoom = 1;
        if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
            numLoopsZoom = 2;
        end
        %loop for zoom plot
        AXorig  = gca;
        for iZoom=1:numLoopsZoom
            if iZoom==2  %ZOOM PLOT
                DESVARIA.Task   = 'create magnifier';
                [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
                haxZ            = DESVARIA.haxZ;
            end
            if isfield(PLOT,'saveProfile') && PLOT.saveProfile
                %skip plotting fields if plotting/saving profiles switched on
                disp(['             ',PLOT.timeString1])
            else
                %plotting
                if strcmp(PLOT.Style,'contour')
                    if strcmp(FIELD.name,'Streamfunction'); plotLineColor=PLOT.streamColor; else; plotLineColor='none'; end
                    [C,h] = contourf(x2dp,z2dp,var2d,PLOT.numContours,'LineColor',plotLineColor);
                elseif strcmp(PLOT.Style,'pcolor')
                    sp = pcolor(x2dp,z2dp,var2d); set(sp,'edgecolor','none');
                    shading interp
                else
                    error('Unknown plotting style!')
                end
                %.........
%                     drawnow  %actually only needed if topo or other graph
%                     plot is included (axis change otherwise) - now done
%                     in individual graph routines
                %.........
                %saving field to file
                if FIELD.save || logical(0) %save current field
                    if ~isfield(SAVE,'Directory'); SAVE.Directory = FILE.directory; end
                    f_saveField(x2dp,1,z2dp,var2d,FIELD,FILE,SAVE);
                end
                if iZoom==1 %main plot
                    % colorbar...
                    DESVARIA.Task    = 'make colorbar';
                    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
                    % title...
                    DESVARIA.Task    = 'make title';
                    [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
                    % annotations...
                    xlabel(PLOT.xlabel)
                    ylabel(PLOT.ylabel)
                    if SWITCH.Texting
                        %             text(max(xlim),min(ylim),['max val: ',num2str(max(var2d(:)),2)],'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',8)
                        text(max(xlim),min(ylim),['min val: ',num2str(min(var2d(:)),2),'; max val: ',num2str(max(var2d(:)),2)],'HorizontalAlignment','right','VerticalAlignment','bottom','Color',STYLE.keyColor,'FontName',STYLE.AllFontName,'FontSize',8);
                        text(min(xlim),min(ylim),[num2str(nx),'x',num2str(nz)],'HorizontalAlignment','left','VerticalAlignment','bottom','Color',STYLE.keyColor,'FontName',STYLE.AllFontName,'FontSize',8);
                        disp(['   ',FIELD.name,' variation: ',num2str(min(var2d(:)),4),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(max(var2d(:)),4),' ',FIELD.dim])
                    else
                        % text(min(xlim),min(ylim),[num2str(nx),'x',num2str(nz)],'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8);
                    end
                    
                else %zoom plot
                    if CBAR.constant
                        set(gca,'clim',[FIELD.cColorbarMIN,FIELD.cColorbarMAX])
                    else
                        cb = PLOT.cb;
                        set(gca,'clim',cb.Limits)
                    end
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    %lineColor = [0.9 0.9 0.9]; %white-ish
                    % axes...
                    DESVARIA.axesLimits     = PLOT.magnifierExtent(PLOT.loopCase,:);
                end
                % axes
                DESVARIA.flipAxes   = true;  % reverse axis (2.)
                DESVARIA.Task  	= 'setup axes';
                [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            end
            
            % EXTRACT Z-PROFILE
            if isfield(PLOT,'zProfile')
                if PLOT.zProfile
                    aspectRatio     = nx/nz;
                    factorXprof     = aspectRatio/2;
                    xProf           = PLOT.x_prof*factorXprof; %make sure it works simultaneously with different aspect ratios (standard is 2:1)
                    xProf           = round( size(z2dp,2)*xProf );
                    if isfield(PLOT,'z_prof'); zProf = round( size(x2dp,2)*PLOT.z_prof ); end %horizontal profile
                    hold on %black indication line
                    if isfield(PLOT,'z_prof')   %horizontal profile
                        plot([min(x2dp(:,1)) max(x2dp(:,1))],[z2dp(1,zProf) z2dp(1,zProf)],'w')
                        plot([min(x2dp(:,1)) max(x2dp(:,1))],[z2dp(1,zProf) z2dp(1,zProf)],'--k')
                    else                        %vertical profile
                        plot([x2dp(xProf,1) x2dp(xProf,1)],[min(z2dp(1,:)) max(z2dp(1,:))],'w')
                        plot([x2dp(xProf,1) x2dp(xProf,1)],[min(z2dp(1,:)) max(z2dp(1,:))],'--k')
                    end
                    figure(21) %profile
                    if PLOT.loopCase==1; clf; end
                    hold on
                    if isfield(PLOT,'z_prof')   %horizontal profile
                        PLOT.hprof(PLOT.loopCase,1) = plot(x2dp(:,zProf),var2d(:,zProf),...
                            'Color',PLOT.caseColor(PLOT.loopCase,:),'LineStyle',PLOT.lineStyle{1,PLOT.loopCase},...
                            'LineWidth',PLOT.lineWidth(1,PLOT.loopCase));
                        if isfield(PLOT,'saveProfile') && PLOT.saveProfile
                            %SAVE PROFILES
                            PLOT.prof_xzdata(:,PLOT.loopCase,PLOT.timestep_profile)    = x2dp(:,zProf);
                            PLOT.prof_vardata(:,PLOT.loopCase,PLOT.timestep_profile)   = var2d(:,zProf);
                        end
                    else                        %vertical profile
                        PLOT.hprof(PLOT.loopCase,1) = plot(var2d(xProf,:),z2dp(xProf,:),...
                            'Color',PLOT.caseColor(PLOT.loopCase,:),'LineStyle',PLOT.lineStyle{1,PLOT.loopCase},...
                            'LineWidth',PLOT.lineWidth(1,PLOT.loopCase));
                        if isfield(PLOT,'saveProfile') && PLOT.saveProfile
                            %SAVE PROFILES
                            PLOT.prof_xzdata(:,PLOT.loopCase,PLOT.timestep_profile)    = z2dp(xProf,:)';
                            PLOT.prof_vardata(:,PLOT.loopCase,PLOT.timestep_profile)   = var2d(xProf,:)';
                        end
                    end
                    axis tight; box on; grid on
                    if isfield(PLOT,'z_prof')   %horizontal profile
                        title(['Profile at z = ',num2str(z2dp(1,zProf),3)])
                        ylabel([num2str(FIELD.name),' [',FIELD.dim,']']); xlabel(PLOT.xlabel)
                    else
                        axis ij
                        title(['Profile at x = ',num2str(x2dp(xProf,1),3)])
                        xlabel([num2str(FIELD.name),' [',FIELD.dim,']']); ylabel(PLOT.ylabel)
                    end
                    if isfield(PLOT,'titleString'); annotProfiles = PLOT.titleString(1,1:PLOT.loopCase);
                    else annotProfiles = ['\bf{',FIELD.name,' }\rm{ ',num2str(PLOT.time2plot,3),' ',PLOT.time2plotDim,'} ']; end
                    legend(PLOT.hprof(1:PLOT.loopCase,1),annotProfiles);
                    %axis background colouring
                    if ~strcmp(STYLE.ColorMode,'light')
                        set(gca,'color',STYLE.ColorModeBW);
                    end
                    if SWITCH.PlotDesign
                        f_DesignFigure;
                    end
                    %saving profile figure
                    if SAVE.Figure
                        SAVE.FigureNr               = 21;
                        SAVE.FigureName             = ['+Profile',num2str(PLOT.number),FIELD.name];
                        if strcmp(SAVE.writeDirectory,'auto') %save to standard folder
                            SAVE.Directory          = FILE.stemSave;
                        else
                            SAVE.Directory          = SAVE.writeDirectory;
                        end
                        SAVE.keepAspect             = true;
                        f_saveFigure( SAVE );
                    end
                    figure(1)
                end
            end
            if iZoom==2; axes(AXorig); end %activate current axis
        end %zoom loop
    end %2-D/3-D
    
    
    
elseif strcmp(GRID.Type,'spherical2D')
    %% CYLINDRICAL PLOT ******************************************************************
    plotCylCart = logical(0);
    %---- spherical2D-Cartesian plot -----------------------------------
    %set up spherical2D-Cartesian figure...
    if plotCylCart
        figure(2)
        if SP.number==1; clf(figure(2)); end
        subplot(PLOT.nrPlots,1,SP.number)
        %setup grid...
        x2d_C(:,:) = GRID.X_3Dnd(:,1,:); %x2dp;
        z2d_C(:,:) = GRID.Z_3Dnd(:,1,:); %z2dp;
        %plotting...
        if strcmp(PLOT.Style,'contour')
            if strcmp(FIELD.name,'Streamfunction'); plotLineColor=PLOT.streamColor; else; plotLineColor='none'; end
            contourf(x2d_C,z2d_C,var2d,PLOT.numContours,'LineColor',plotLineColor)
            
        elseif strcmp(PLOT.Style,'pcolor')
            sp=pcolor(x2d_C,z2d_C,var2d); set(sp,'edgecolor','none');
            shading interp
        else
            error('Unknown plotting style!')
        end

        % colorbar...
        DESVARIA.Task    = 'make colorbar';
        [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
        % Texting...
        DESVARIA.Task    = 'make title';
        [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
        % annotations...
        xlabel(PLOT.xlabel)
        ylabel(PLOT.ylabel)
        % axes...
        DESVARIA.flipAxes   = true;  % reverse axis (2.)
        DESVARIA.Task  	= 'setup axes';
        [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
        figure(1) %back to original figure
    end
    
    % ---- spherical2D-spherical2D plot ---------------------------------
    var2d_S         = var2d;
    if SWITCH.closeAnnulus
        var2d_S(end+1,:)    = var2d_S(1,:);
    end
    
    axes(AXcurrent);
    
    %setup for zoom plot
    numLoopsZoom = 1;
    if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
        numLoopsZoom = 2;
    end
    %loop for zoom plot
    AXorig  = gca;
    for iZoom=1:numLoopsZoom
        if iZoom==2  %ZOOM PLOT
            DESVARIA.Task   = 'create magnifier';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            haxZ            = DESVARIA.haxZ;
        end
        %plotting...
        if strcmp(PLOT.Style,'contour')
            if strcmp(FIELD.name,'Streamfunction'); plotLineColor=PLOT.streamColor; else; plotLineColor='none'; end
            contourf(GRID.x2ds,GRID.z2ds,var2d_S,PLOT.numContours,'LineColor',plotLineColor)
        elseif strcmp(PLOT.Style,'pcolor')
            sp=pcolor(GRID.x2ds,GRID.z2ds,var2d_S); set(sp,'edgecolor','none');
            shading interp
        end
        if (FIELD.save && iZoom==1) || logical(0) %save current field
            if ~isfield(SAVE,'Directory'); SAVE.Directory = FILE.directory; end
            f_saveField(GRID.x2ds,1,GRID.z2ds,var2d_S,FIELD,FILE,SAVE);
        end
        %plot bot CMB and top air line
        hold on
        hcmb = plot(GRID.x2ds(:,end),GRID.z2ds(:,end),'Color',STYLE.keyColor,'LineWidth',0.25);
        hold on
        hair = plot(GRID.x2ds(:,1),GRID.z2ds(:,1),'Color',STYLE.keyColor,'LineWidth',0.25);
        if iZoom==1 %main plot
            % colorbar...
            DESVARIA.Task    = 'make colorbar';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            % title...
            DESVARIA.Task    = 'make title';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            % annotations...
            xlabel(PLOT.xlabel);
            ylabel(PLOT.ylabel);
            if SWITCH.spherical2DCenterText
                text(0.0,0.0,[num2str(nx),'\times',num2str(nz)],'Color',STYLE.keyColor,'FontName',STYLE.AllFontName,...
                    'HorizontalAlignment','center')
            end
            if SWITCH.Texting %display min/max values
                disp(['   ',FIELD.name,' variation: ',num2str(min(var2d_S(:)),4),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(max(var2d_S(:)),4),' ',FIELD.dim])
            end
            % axes...
            BIndicationHeight = 0; %adjust for plate boundary indications
            if PLOT.indicateTrench || PLOT.indicateRidge
                BIndicationHeight = abs(max(GRID.Z_3Dp(:))-min(GRID.Z_3Dp(:))) /10;
            end
            DESVARIA.axesLimits = [min(GRID.x2ds(:))-BIndicationHeight,max(GRID.x2ds(:))+BIndicationHeight,min(GRID.z2ds(:))-BIndicationHeight,max(GRID.z2ds(:))+BIndicationHeight];
            DESVARIA.equalAxes      = true;
            DESVARIA.axisColour     = 'none';
            DESVARIA.Task  	= 'setup axes';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            
        else %zoom plot
            if CBAR.constant
                set(gca,'clim',[FIELD.cColorbarMIN,FIELD.cColorbarMAX])
            else
                cb = PLOT.cb;
                set(gca,'clim',cb.Limits)
            end
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            % axes...
            DESVARIA.axesLimits     = PLOT.magnifierExtent(PLOT.loopCase,:);
            DESVARIA.equalAxes      = true;
            DESVARIA.Task  	= 'setup axes';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
        end
        if iZoom==2; axes(AXorig); end %activate current axis
    end  %magnifier plot
    
    % EXTRACT Z-PROFILE
    if isfield(PLOT,'zProfile')
        if PLOT.zProfile
            aspectRatio     = nx/nz;
            factorXprof     = aspectRatio/2;
            xProf           = PLOT.x_prof*factorXprof; %make sure it works simultaneously with different aspect ratios (standard is 2:1)
            xProf           = round( size(z2dp,2)*xProf );
            if isfield(PLOT,'z_prof'); zProf = round( size(x2dp,2)*PLOT.z_prof ); end %horizontal profile
            hold on %black indication line
            if isfield(PLOT,'z_prof')   %horizontal profile
                xdummy  = [min(GRID.x2dp(1,:)), max(GRID.x2dp(1,:))];
                zdummy  = [GRID.z2dp(zProf,1), GRID.z2dp(zProf,1)];
            else                        %vertical profile
                xdummy  = [GRID.x2dp(xProf,1), GRID.x2dp(xProf,1)];
                zdummy  = [min(GRID.z2dp(1,:)), max(GRID.z2dp(1,:))];
            end
            rData   = GRID.rcmb_p+max(GRID.Z_3Dp(1,1,:)) - zdummy; %calculate radii
            [xdummy,zdummy] = f_cartesian2spherical(xdummy,rData,'degrees'); %'degrees', or 'plottingDimension'
            plot([xdummy(1,1) xdummy(1,2)],[zdummy(1,1) zdummy(1,2)],'w')
            plot([xdummy(1,1) xdummy(1,2)],[zdummy(1,1) zdummy(1,2)],'--k')
            figure(21) %profile
            if PLOT.loopCase==1; clf; end
            hold on
            if isfield(PLOT,'z_prof')   %horizontal profile
                PLOT.hprof(PLOT.loopCase,1) = plot(x2dp(:,zProf),var2d(:,zProf),...
                    'Color',PLOT.caseColor(PLOT.loopCase,:),'LineStyle',PLOT.lineStyle{1,PLOT.loopCase},...
                    'LineWidth',PLOT.lineWidth(1,PLOT.loopCase));
                if isfield(PLOT,'saveProfile') && PLOT.saveProfile
                    %SAVE PROFILES
                    PLOT.prof_xzdata(:,PLOT.loopCase,PLOT.timestep_profile)    = x2dp(:,zProf);
                    PLOT.prof_vardata(:,PLOT.loopCase,PLOT.timestep_profile)   = var2d(:,zProf);
                end
            else                        %vertical profile
                PLOT.hprof(PLOT.loopCase,1) = plot(var2d(xProf,:),z2dp(xProf,:),...
                    'Color',PLOT.caseColor(PLOT.loopCase,:),'LineStyle',PLOT.lineStyle{1,PLOT.loopCase},...
                    'LineWidth',PLOT.lineWidth(1,PLOT.loopCase));
                if isfield(PLOT,'saveProfile') && PLOT.saveProfile
                    %SAVE PROFILES
                    PLOT.prof_xzdata(:,PLOT.loopCase,PLOT.timestep_profile)    = z2dp(xProf,:)';
                    PLOT.prof_vardata(:,PLOT.loopCase,PLOT.timestep_profile)   = var2d(xProf,:)';
                end
            end
            axis tight; box on; grid on
            if isfield(PLOT,'z_prof')   %horizontal profile
                title(['Profile at z = ',num2str(z2dp(1,zProf),3)])
                ylabel([num2str(FIELD.name),' [',FIELD.dim,']']); xlabel(PLOT.xlabel)
            else
                axis ij
                title(['Profile at x = ',num2str(x2dp(xProf,1),3)])
                xlabel([num2str(FIELD.name),' [',FIELD.dim,']']); ylabel(PLOT.ylabel)
            end
            if isfield(PLOT,'titleString'); annotProfiles = PLOT.titleString(1,1:PLOT.loopCase);
            else annotProfiles = ['\bf{',FIELD.name,' }\rm{ ',num2str(PLOT.time2plot,3),' ',PLOT.time2plotDim,'} ']; end
            legend(PLOT.hprof(1:PLOT.loopCase,1),annotProfiles);
            %axis background colouring
            if ~strcmp(STYLE.ColorMode,'light')
                set(gca,'color',STYLE.ColorModeBW);
            end
            if SWITCH.PlotDesign
                f_DesignFigure;
            end
            %saving profile figure
            if SAVE.Figure
                SAVE.FigureNr               = 21;
                SAVE.FigureName             = ['+Profile',num2str(PLOT.number),FIELD.name];
                if strcmp(SAVE.writeDirectory,'auto') %save to standard folder
                    SAVE.Directory          = FILE.stemSave;
                else
                    SAVE.Directory          = SAVE.writeDirectory;
                end
                SAVE.keepAspect             = true;
                f_saveFigure( SAVE );
            end
            figure(1)
        end
    end
    
    
elseif strcmp(GRID.Type,'yinyang')
    %% YINYANG MODE ******************************************************************
    add_streamline = logical(0); %not tested yet...................
    
    trackPlumes = PLOT.indicatePlumes;
    AXorig  = gca; %axis of big plot (needs to be changed later on for yy)
    
    if strcmp(SWITCH.yyPlotMode,'mapOnSphere')
        PLOT.yyNumSlices = 1; %multiple slices not possible with spherical plot
    end
    
    %% GET SLICING DEPTH-LEVELS
    gp_z    = size(GRID.X_3D,3);         % # grid points in z-direction
    znumV   = zeros(1,size(PLOT.yyZlevelV,2));
    for i_zlevel=1:size(PLOT.yyZlevelV,2)
        [~,idx]                     = min(abs(GRID.Z_3Dnd(1,1,:)-PLOT.yyZlevelV(1,i_zlevel))); %index of closest value
        PLOT.yyZlevelV(1,i_zlevel)	= GRID.Z_3Dnd(1,1,idx); %set to closest grid value
        znumV(1,i_zlevel)           = idx;  %index of z-level
    end
    if strcmp(SWITCH.yyPlotMode,'isosurface') || trackPlumes %read out complete data
        znumVall(1,1:gp_z)   = 1:gp_z;
    end
    
    %% SCALAR FIELDS
    ayy = zeros(size(VARyin,1),size(VARyin,2),2);
    nth = size(ayy,1); nph=nth*3; nthmap=2*nth; nphmap=(4*nph)/3;
    if ~strcmp(FIELD.name,'Horizontal velocity') && ~strcmp(FIELD.name,'Radial velocity') %scalar field
        if strcmp(FIELD.name,'Plate velocity') || strcmp(FIELD.name,'Topography') || strcmp(FIELD.name,'Geoid') || ...  %plot graphs
                    strcmp(FIELD.name,'Dyn. topography') || strcmp(FIELD.name,'Iso. topography') || strcmp(FIELD.name,'Res. topography') || ...
                    strcmp(FIELD.name,'Heat flux') || strcmp(FIELD.name,'Surface age') || strcmp(FIELD.name,'Crustal thickness') || ...
                    strcmp(FIELD.name,'Plate-base topography') || strcmp(FIELD.name,'Topography (self-grav)')
            if strcmp(SWITCH.yyPlotMode,'isosurface'); error('Topography plot for YY-isosurface not implemented yet!'); end
            if size(VARyin,3)>1
                if TOPO.cmb; iTopo = 1; else; iTopo = 2; end
            else
                iTopo = 1;
            end
            %disp([FIELD.name,' mean: ',num2str(mean(mean(VARyin(:,:,iTopo)))+mean(mean(VARyang(:,:,iTopo))))]) %not exact due to variable gridspacing
            ayy(:,:,1) = VARyin(:,:,iTopo);
            ayy(:,:,2) = VARyang(:,:,iTopo);
            if SWITCH.plotTOPOsmooth
                filterSize = 5;
                %F = [.05 .1 .05; .1 .4 .1; .05 .1 .05]; %filter
                %F = [1/9 1/9 1/9; 1/9 1/9 1/9; 1/9 1/9 1/9]; %like 2-D filter - with one loop!
                F = fspecial('gaussian',filterSize); %filterSize is the size of matrix
                %F = fspecial('gaussian');
                
                %add ghost points to prevent side effects from smoothing
                dummy = ayy;
                for ii=1:filterSize
                    dummy = [dummy(1,:,:); dummy; dummy(end,:,:)];
                    dummy = [dummy(:,1,:), dummy, dummy(:,end,:)];
                end
                %smoothing
                for jj=1:17
                    ayy_smooth = conv2(dummy(:,:,1),F,'same');
                    dummy(:,:,1) = ayy_smooth(:,:);
                    ayy_smooth = conv2(dummy(:,:,2),F,'same');
                    dummy(:,:,2) = ayy_smooth(:,:);
                end
                %remove ghostpoints
                ayy_smooth = dummy(filterSize+1:end-filterSize,filterSize+1:end-filterSize,:);
                ayy = ayy_smooth;
            end
            [amap(:,:,1),theta,phi] = f_YYtoMap(ayy);
            
        else %OTHER
            amap = zeros(nthmap,nphmap,size(znumV,2));
            for i_zlevel=1:size(znumV,2)  %loop through the different depth slices
                ayy(:,:,1) = VARyin(:,:,znumV(1,i_zlevel));
                ayy(:,:,2) = VARyang(:,:,znumV(1,i_zlevel));
                [amap(:,:,i_zlevel),theta,phi] = f_YYtoMap(ayy);
                clearvars ayy
            end
            if trackPlumes %convert all the T-data
                ayy = zeros(size(PLOT.T_3Dyin,1),size(PLOT.T_3Dyin,2),2);
                pmap = zeros(nthmap,nphmap,size(znumVall,2));
                for i_zlevel=1:size(znumVall,2)  %loop through the different depth slices
                    ayy(:,:,1) = PLOT.T_3Dyin(:,:,znumVall(1,i_zlevel));
                    ayy(:,:,2) = PLOT.T_3Dyang(:,:,znumVall(1,i_zlevel));
                    [pmap(:,:,i_zlevel),~,~] = f_YYtoMap(ayy);
                    clearvars ayy
                end
            end
            if strcmp(SWITCH.yyPlotMode,'isosurface') %convert all the data
                ayy = zeros(size(VARyin,1),size(VARyin,2),2);
                amap = zeros(nthmap,nphmap,size(znumVall,2));
                for i_zlevel=1:size(znumVall,2)  %loop through the different depth slices
                    ayy(:,:,1) = VARyin(:,:,znumVall(1,i_zlevel));
                    ayy(:,:,2) = VARyang(:,:,znumVall(1,i_zlevel));
                    [amap(:,:,i_zlevel),~,~] = f_YYtoMap(ayy);
                    clearvars ayy
                end
            end
        end
    end
    if SWITCH.waitbar; waitbar(0.5,PLOT.wb); end
    
    %% SETTING SWITCHES
    %quiver switch
    if strcmp(FIELD.name,'Temperature')  && quiver_T || ... %for quiver addition
            strcmp(FIELD.name,'Viscosity')    && quiver_eta || ...
            strcmp(FIELD.name,'Stress')       && quiver_str || ...
            strcmp(FIELD.name,'Strain rate')  && quiver_edot || ...
            strcmp(FIELD.name,'Density')      && quiver_rho || ...
            strcmp(FIELD.name,'Streamfunction') && quiver_psi || ...
            strcmp(FIELD.name,'Velocity') && quiver_v || ...
            strcmp(FIELD.name,'Horizontal velocity') && quiver_vx || ...
            strcmp(FIELD.name,'Radial velocity') && quiver_vr || ...
            strcmp(FIELD.name,'zz-Stress component') && quiver_nstr
        add_quiver = true;
    else
        add_quiver = false;
    end
    %streamfunction switch
    if (strcmp(FIELD.name,'Temperature')  && streamfun_T) || ... %for streamfun addition
            (strcmp(FIELD.name,'Viscosity')      && streamfun_eta) || ...
            (strcmp(FIELD.name,'Stress')         && streamfun_str) || ...
            (strcmp(FIELD.name,'Strain rate')    && streamfun_edot) || ...
            (strcmp(FIELD.name,'Density')        && streamfun_rho) || ...
            (strcmp(FIELD.name,'Velocity')       && streamfun_v) || ...
            (strcmp(FIELD.name,'Horizontal velocity') && streamfun_vx) || ...
            (strcmp(FIELD.name,'Radial velocity') && streamfun_vr)
        add_streamfun = true;
    else
        add_streamfun = false;
    end
    
    %% VECTOR FIELDS
    if strcmp(FIELD.name,'Horizontal velocity') || strcmp(FIELD.name,'Radial velocity') || strcmp(FIELD.name,'Streamfunction') || ... %VECTOR FIELD
            add_streamfun || add_quiver %or quiver/streamfunction addition
        
        % Transform coordinates for Yin & Yang grids        %<<< MAYBE DO THAT IN SL_FieldPlot.m AND SAVE TO GLOBAL GRID VARIABLE.................
        R  = GRID.Z_3Dnd+GRID.rcmb_nd;     lat = pi/4-GRID.X_3Dnd;     lon = GRID.Y_3Dnd-3*pi/4;
        XS_1 = R.*cos(lat).*cos(lon); %Yin
        YS_1 = R.*cos(lat).*sin(lon);
        ZS_1 = R.*sin(lat);
        XS_2 = -XS_1; %Yang
        YS_2 =  ZS_1;
        ZS_2 =  YS_1;
        
        % Transform velocities, if needed
        Vtheta = VX_3Dyin; Vphi = VY_3Dyin; Vr = VZ_3Dyin;            	% on Yin grid
        VX_3D_1 =  Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon);
        VY_3D_1 =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon);
        VZ_3D_1 = -Vtheta.*cos(lat)                            + Vr.*sin(lat);
        Vtheta = VX_3Dyang; Vphi = VY_3Dyang; Vr = VZ_3Dyang;        	% on Yang grid
        VX_3D_2 = -( Vtheta.*sin(lat).*cos(lon) - Vphi.*sin(lon) + Vr.*cos(lat).*cos(lon) );
        VZ_3D_2 =  Vtheta.*sin(lat).*sin(lon) + Vphi.*cos(lon) + Vr.*cos(lat).*sin(lon);
        VY_3D_2 = -Vtheta.*cos(lat)                            + Vr.*sin(lat);
        
        %conversion to spherical coordinates
        lat1    = atan2(sqrt(XS_1.^2 + YS_1.^2),ZS_1);  % cos-1 (z/r)
        long1   = atan2(YS_1,XS_1);  % tan-1 (y/x)
        lat2    = atan2(sqrt(XS_2.^2 + YS_2.^2),ZS_2);  % cos-1 (z/r)
        long2   = atan2(YS_2,XS_2);  % tan-1 (y/x)
        
        Vlat1   = VX_3D_1.*(cos(long1).*cos(lat1)) + VY_3D_1.*(sin(long1).*cos(lat1)) - VZ_3D_1.*(sin(lat1));
        Vlong1	= -VX_3D_1.*(sin(long1)) + VY_3D_1.*(cos(long1));
        Vlat2   = VX_3D_2.*(cos(long2).*cos(lat2)) + VY_3D_2.*(sin(long2).*cos(lat2)) - VZ_3D_2.*(sin(lat2));
        Vlong2  = -VX_3D_2.*(sin(long2)) + VY_3D_2.*(cos(long2));
        
        nth = size(VX_3Dyin,1); nph=nth*3; nthmap=2*nth; nphmap=(4*nph)/3;
        Vtheta_map = zeros(nthmap,nphmap,size(znumV,2));
        Vphi_map = Vtheta_map;
        Vr_map = Vtheta_map;
        for i_zlevel=1:size(znumV,2)  %loop through the different depth slices
            VXyy(:,:,1) = Vlat1(:,:,znumV(1,i_zlevel)); % Yin  grid data
            VXyy(:,:,2) = Vlat2(:,:,znumV(1,i_zlevel)); % Yang grid data
            [Vtheta_map(:,:,i_zlevel),~,~ ]   = f_YYtoMap(VXyy);
            
            VYyy(:,:,1) = Vlong1(:,:,znumV(1,i_zlevel)); % Yin  grid data
            VYyy(:,:,2) = Vlong2(:,:,znumV(1,i_zlevel)); % Yang grid data
            %             [Vphi_map(:,:,i_zlevel),~ ,~ ]   = f_YYtoMap(VYyy);
            [Vphi_map(:,:,i_zlevel),~,~ ]   = f_YYtoMap(VYyy); %use this for old matlab version
            
            VZyy(:,:,1) = VZ_3D_1(:,:,znumV(1,i_zlevel)); % Yin  grid data
            VZyy(:,:,2) = VZ_3D_2(:,:,znumV(1,i_zlevel)); % Yang grid data
            %             [Vr_map(:,:,i_zlevel),~ ,~ ]   = f_YYtoMap(VZyy);
            [Vr_map(:,:,i_zlevel),theta,phi ]   = f_YYtoMap(VZyy); %use this for old matlab version
        end
        
        %         Vphi_map is horizontal velocity
        %         Vtheta_map is vertical velocity
        
        if strcmp(FIELD.name,'Streamfunction') || add_streamfun %add_streamfunction
            %adjust latitude range
            if SWITCH.yyFlipVertical
                if exist('Vtheta_map','var')
                    Vtheta_mapSTREAM = -Vtheta_map;
                end
            end
            
            phi2d_1 = zeros(size(Vphi_map,1),size(Vphi_map,2),size(znumV,2));
            psi2d_1 = phi2d_1;
            for i_zlevel=1:size(znumV,2)
                [phi2d_1(:,:,i_zlevel),psi2d_1(:,:,i_zlevel)] = flowfun(Vphi_map(:,:,i_zlevel),Vtheta_mapSTREAM(:,:,i_zlevel));  % calculate the velocity potential (phi) and streamfunction (psi)
                %                 [phi2d_2(:,:,i_zlevel),psi2d_2(:,:,i_zlevel)] = flowfun(Vlat2(:,:,z_num_V(1,i_zlevel)),Vlong2(:,:,z_num_V(1,i_zlevel)));
            end
            if strcmp(stream_data,'phi'); stream2d(:,:,:)=phi2d_1;   %phi2d_1;
            elseif strcmp(stream_data,'psi'); stream2d(:,:,:)=psi2d_1; end
            
            if strcmp(FIELD.name,'Streamfunction')
                for i_zlevel=1:size(znumV,2)  %loop through the different depth slices
                    %                     ayy(:,:,1) = VAR_3Dyin(:,:,z_num_V(1,i_zlevel));
                    %                     ayy(:,:,2) = VAR_3Dyang(:,:,z_num_V(1,i_zlevel));
                    %                     [amap(:,:,i_zlevel),theta,phi] = f_YYtoMap(ayy);
                    
                    amap(:,:,i_zlevel) = stream2d(:,:,i_zlevel);
                    
                    figure(11),clf
                    subplot(2,2,1)
                    pcolor(Vtheta_mapSTREAM(:,:,i_zlevel))
                    subplot(2,2,2)
                    pcolor(Vphi_map(:,:,i_zlevel))
                    subplot(2,2,4)
                    pcolor(amap(:,:,i_zlevel))
                    figure(1)
                end
            end
        end
    end
    
    %for velocity plots
    if strcmp(FIELD.name,'Horizontal velocity') && ~strcmp(FIELD.name,'Radial velocity'); Vr_map = 0*Vr_map; end  %zeroing vertical component
    if strcmp(FIELD.name,'Horizontal velocity'); amap = Vtheta_map; end
    %     if plotVelocityY; amap = Vphi_map; end
    if strcmp(FIELD.name,'Radial velocity'); amap = Vr_map; end
    
    %% calculate RESIDUAL TEMPERATURE FIELD
    if strcmp(FIELD.name,'Horizontal residual') %FOR HORIZONTAL-RESIDUAL TEMPERATURE FIELD
        for z_level2=1:size(amap,3)
            %sumBoth = nansum(nansum( amap(:,:,z_level2) )); %summing up T values
            numNaN = sum(sum(isnan( amap(:,:,z_level2) )));             %count NaN's  (not really needed, just to be sure...)
            %numBoth = (size( amap(:,:,z_level2),1)*size(VAR_3Dyin(:,:,z_level2),2)) - numNaN; %number of T values
            %horizMeanBoth = sumBoth / numBoth;          %horizonatl mean
            
            horizMeanBoth       = mean2(amap(:,:,z_level2));            %horizontal mean temperature
            amap(:,:,z_level2)  = amap(:,:,z_level2) - horizMeanBoth;   %residual temperature
        end
    end
    if strcmp(FIELD.name,'Global residual') %FOR GLOBAL-RESIDUAL TEMPERATURE FIELD
        for z_level2=1:size(amap,3)
            numNaN              = sum(sum(isnan( amap(:,:,:) )));       %count NaN's  (not really needed, just to be sure...)
            horizMeanBoth       = mean2(amap(:,:,:));                   %horizontal mean temperature
            amap(:,:,z_level2)  = amap(:,:,z_level2) - horizMeanBoth;   %residual temperature
        end
    end
    if trackPlumes  %FOR PLUMES RESIDUAL TEMPERATURE FIELD
        PLUMES = zeros(size(pmap));
        for z_level2=1:size(pmap,3)
            horizMeanBoth       = mean2(pmap(:,:,z_level2));    %horizontal mean temperature
            pmap(:,:,z_level2)  = pmap(:,:,z_level2) - horizMeanBoth;  %residual temperature
            horizMean_Tres      = mean2(pmap(:,:,z_level2));    %horizontal mean Tres
            
            if trackPlumes
                horizMax    = max(max(pmap(:,:,z_level2)));     %horizontal max Tres
                horizMin    = min(min(pmap(:,:,z_level2)));     %horizontal min Tres
                
                thrHot      = horizMean_Tres + PLOT.pHot*(horizMax-horizMean_Tres);     %after Labrosse (EPSL,2002)
                thrCold     = horizMean_Tres + PLOT.pCold*(horizMin-horizMean_Tres);    %after Labrosse (EPSL,2002)
                
                plumes      = zeros(size(pmap(:,:,z_level2)));	%no plumes
                plumes(find(pmap(:,:,z_level2)>thrHot))	= 1;  	%hot plumes
                plumes(find(pmap(:,:,z_level2)<thrCold)) = -1; 	%cold plumes
                
                PLUMES(:,:,z_level2) = plumes;
            end
        end
    end
    
    %% GRID ADJUSTMENTS
    %conversions to degrees
    theta       = rad2deg(theta);
    phi         = rad2deg(phi);
    %adjust latitude range
    if SWITCH.yyFlipVertical
        theta   = -theta+90; %convert to -90?-90? instead of 0?-180?  %CORRECT NORTHENING
        if exist('Vtheta_map','var')
            Vtheta_map = -Vtheta_map;
        end
    else
        theta   = theta-90; %convert to -90?-90? instead of 0?-180?  %CORRECT NORTHENING
    end
    
    %% TRACK PLUMES
    if trackPlumes
        plot3D_plumes               = logical(0);
        f_INPUT.plotPLUMES      	= plot3D_plumes;%plot comparison of all anomalies and selected plumes
        f_INPUT.addFactor           = 0.25;         %percent of total x-extend that is added at x=1 and x=end  (should be exponent of 2)
        
        f_INPUT.upperDepthThresholdHot 	= 0.5;          %top threshold; percentage of nz (0: top, 1: bottom)   ...incl. sticky air!
        f_INPUT.lowerDepthThresholdHot	= 0.9;          %top threshold; percentage of nz (0: top, 1: bottom)
        f_INPUT.upperDepthThresholdCold	= 0.2;          %top threshold; percentage of nz (0: top, 1: bottom)   ...incl. sticky air!
        f_INPUT.lowerDepthThresholdCold	= 0.4;          %top threshold; percentage of nz (0: top, 1: bottom)
        
        [PLUMEShot,PLUMEScold] = f_trackPlumes(PLUMES,f_INPUT);
        figure(1)  %go back to figure(1)
    end
    
    %% REMOVING DATA POINTS
    if strcmp(FIELD.name,'Velocity') && add_quiver && ~add_streamline
        theta_V = theta;
        phi_V = phi;
        %remove some data points
        count=0;
        for ir=size(Vtheta_map,1):-1:1
            count = count+1;
            if count==quiver_numDel
                count = 0;
            else
                theta_V(ir,:) = [];
                Vtheta_map(ir,:,:) = [];
                Vphi_map(ir,:,:) = [];
            end
        end
        for jr=size(Vtheta_map,2)-1:-1:2
            count = count+1;
            if count==quiver_numDel
                count = 0;
            else
                phi_V(jr,:) = [];
                Vtheta_map(:,jr,:) = [];
                Vphi_map(:,jr,:) = [];
            end
        end
    end
    if SWITCH.waitbar; waitbar(0.6,PLOT.wb); end
    
    %% PLOTTING
    if strcmp(SWITCH.yyPlotMode,'map') && PLOT.yyNumSlices>1
        figure(1), clf %<<<<<<<<<  ONLY USED FOR MULTIPLE SLICES PER SUBPLOT
    end
    if SWITCH.waitbar
        disp('   ...plotting')
        PLOT.wb = waitbar(0.6,PLOT.wb,'plotting...');
    end
    if strcmp(SWITCH.yyPlotMode,'isosurface')
        %CONVERT TO CARTESIAN COORDINATES
        r               = GRID.Z_3Dndnf+GRID.rcmb_nd;
        r2d(:,1)        = r(1,1,:);
        %         phi = phi+abs(min(phi));
        %         theta = theta+abs(min(theta));
        [phiM,thetaM,rM]= meshgrid(deg2rad(phi),deg2rad(theta),r2d); %has to be in radians and from 0:360 and 0:180!!!
        [x,y,z]       	= sph2cart(phiM(:),thetaM(:),rM(:));
        
        amap2           = amap(:);
        
        %reducing data size
        isoval          = 1100;
        range           = 20;
        dummy           = amap2>(isoval+range);
        x(dummy)        = [];
        y(dummy)        = [];
        amap2(dummy)    = [];
        dummy           = amap2<(isoval-range);
        x(dummy)        = [];
        y(dummy)        = [];
        z(dummy)        = [];
        amap2(dummy)    = [];
        
        %create and interpolate onto regular meshgrid
        x_min           = min(x(:)); x_max = max(x(:));
        y_min           = min(y(:)); y_max = max(y(:));
        z_min           = min(z(:)); z_max = max(z(:));
        x_size          = size(phiM,1);
        y_size          = size(phiM,2);
        z_size          = size(phiM,3);
        dx              = (x_max-x_min)/(x_size-1);  %just make grid fine enough
        dy              = (y_max-y_min)/(y_size-1);  %just make grid fine enough
        dz              = (z_max-z_min)/(z_size-1);  %just make grid fine enough
        fac             = 1; %reducing size
        dx              = dx*fac;
        dy              = dy*fac;
        dz              = dz*fac;
        
        [Xq,Yq,Zq] = meshgrid(x_min:dx:x_max,y_min:dy:y_max,z_min:dz:z_max);
        
        disp('    [grid interpolation: this may take a while...]')
        amap_int = griddata(x,y,z,amap2,Xq,Yq,Zq); %%%slow
        
        isosurface(Xq,Yq,Zq,amap_int,isoval)
        
        axis equal off vis3d
        lighting phong
        camlight('right')
        
        
    elseif strcmp(SWITCH.yyPlotMode,'map') || strcmp(SWITCH.yyPlotMode,'mapOnSphere')
        for i_slice=1:PLOT.yyNumSlices
            if SWITCH.waitbar; waitbar(0.6+i_slice/PLOT.yyNumSlices*0.2,PLOT.wb); end
            z_num               = znumV(i_slice);           %grid number
            z_shift             = (i_slice-1)/4;
            z_level             = PLOT.yyZlevelV(i_slice);	%non-dim depth
            z_level_dim         = z_level*GRID.dimFactor; 	%dim depth [plotting dimension]
            
            %set position:
            if i_slice>1; hold on; end
            if PLOT.yyNumSlices==1
                %nothing to be done
            elseif PLOT.yyNumSlices==2;  hax = axes('Position', [0.04, 0.07+z_shift*1.6, 0.8, 0.5]);
            elseif PLOT.yyNumSlices==3;  hax = axes('Position', [0.07, 0.1+z_shift, 0.75, 0.3]);
            elseif PLOT.yyNumSlices==4;  hax = axes('Position', [0.04, 0.07+z_shift*0.85, 0.75, 0.3]);
            else; error('fc: check here!');
            end
            if i_slice>1 && PLOT.onlyBotPlumes; trackPlumes = false; end
            
            %switches
            axisVisibility          = 'off';
            if strcmp(SWITCH.yyPlotMode,'map')
                SetFrame          	= 'on';
                LabelParallel = SWITCH.yyLabelParallel;  if i_slice>1; LabelParallel='off'; end
                haxm = axesm('mollweid', 'Frame',SetFrame,'FLineWidth',0.25, ...
                    'Grid', 'on','GColor',STYLE.BlackOrWhiteColor,'GLineStyle','-','GLineWidth',0.1, ...
                    'MeridianLabel',SWITCH.yyLabelMeridian,'ParallelLabel',LabelParallel,...
                    'PLabelLocation',PLOT.yyPlabelTicks,'PLabelMeridian',PLOT.yyPlabelLocation,'Origin',[0 PLOT.yyOriginLong 0],...
                    'FontName',STYLE.AllFontName,'FontWeight','normal','FontColor',STYLE.annotationColor);
                if SWITCH.yyLimitMap
                    setm(haxm,'MapLatLimit',[PLOT.yyLimitS PLOT.yyLimitN],'MapLonLimit',[PLOT.yyLimitW PLOT.yyLimitE])
                end
                if strcmp(FIELD.name,'streamfunction'); plotLineColor = PLOT.streamColor; else; plotLineColor = 'none'; end
            end
            set(gca,'Color',STYLE.BlackOrWhiteColor,'Visible',axisVisibility)
            
            %% PLOT MAIN FIELD
            hold on
            if strcmp(SWITCH.yyPlotMode,'map')
                if strcmp(PLOT.Style,'contour')
                    %---------------------------
                    contourfm(theta,phi,amap(:,:,i_slice),PLOT.numContours,'LineColor',plotLineColor)
                    %---------------------------
                else %if strcmp(PLOT.Style,'pcolor')
                    if ~strcmp(PLOT.Style,'pcolor')
                        if SWITCH.Verbose; warning('unknown PLOT.Style!'); end
                    end
                    %---------------------------
                    sp = pcolorm(theta,phi,amap(:,:,i_slice)); set(sp,'EdgeColor','none');
                    %---------------------------                    
                end
                
            elseif strcmp(SWITCH.yyPlotMode,'mapOnSphere')
                Rmax        = max(GRID.Z_3Dp(:))+GRID.rcmb_p; %dimensional radius of planet [plotting dimension]
                Rcurrent    = Rmax-z_level; %radius of slice [plotting dimension]
                [phiM,thetaM,rM]    = meshgrid(deg2rad(phi),deg2rad(theta),Rcurrent); %has to be in radians and from 0:360 and 0:180!!!
                [x,y,z]             = sph2cart(phiM,thetaM,rM);
                %---------------------------
                haxm = gca;
                sp = surf(x,y,z,amap(:,:,i_slice),'EdgeColor','none');
                %---------------------------
                axis([-Rcurrent,Rcurrent,-Rcurrent,Rcurrent,-Rcurrent,Rcurrent]) %keep axis constant and hence, make sphere smaller
                axis equal tight vis3d
                lighting phong
                camlight('right')
            end
            drawnow %needed to prevent weird axis limits when plotting 3 maps on top each other
            
            %% QUIVER
            if SWITCH.yyOnlyTopQuiver; if i_slice~=4; add_quiver=0; else; add_quiver=1; end; end
            if SWITCH.yyOnlyBotQuiver; if i_slice==1; add_quiver=1; else; add_quiver=0; end; end
            if add_quiver
                if i_slice==4
                    arrowScale1 = arrowScale; %3
                elseif i_slice==3
                    arrowScale1 = 2;
                elseif i_slice==2
                    arrowScale1 = 1;
                elseif i_slice==1
                    arrowScale1 = 1;
                else
                    error('check this line!!! fc')
                end
                % arrowScale1 = 0; %no scaling, all same length
                hold on
                [phi2,theta2] = meshgrid(phi_V,theta_V);
                q1 = quiverm(theta2,phi2,Vtheta_map(:,:,i_slice),Vphi_map(:,:,i_slice),'k',arrowScale1);
                set(q1,'Color',PLOT.quiverColor,'LineWidth',PLOT.quiverWidth)
            end
            
            %% STREAM LINES
            if add_streamline
                hold on
                [phi2,theta2] = meshgrid(phi_V,theta_V);
                %             [s1,s2] = meshgrid( phi2(:,10),theta2(:,10));     %starting points
                s1=phi2(5,10);
                s2=theta2(5,10);
                [sx,sy]=meshgrid(s1,s2);
                streamline(theta2,phi2,Vtheta_map(:,:,i_slice),Vphi_map(:,:,i_slice),sx,sy)
            end
            
            %% TRACK PLUMES
            if trackPlumes
                if strcmp(SWITCH.yyPlotMode,'mapOnSphere')
                    error('plumes not yet implemented for mapOnSphere option!')
                end
                numberPlumes = logical(0); %labels every connected plume area (on the 2-D map area) with a number; saved in PLUMEShot_nr + PLUMEScold_nr
                if numberPlumes
                    PLUMEShot_large(:,:)        = PLUMEShot(:,:,z_num);
                    PLUMEScold_large(:,:)       = PLUMEScold(:,:,z_num);
                    add_factor                  = 1; %1 means 1 additional total field at each edge
                    nx_org                      = size(PLUMEShot,2);
                    nx_add                      = nx_org*add_factor;
                    PLUMEShot_large             = [PLUMEShot_large(:,(end-nx_add+1):end)   PLUMEShot_large   PLUMEShot_large(:,1:nx_add)]; %increase size
                    PLUMEScold_large            = [PLUMEScold_large(:,(end-nx_add+1):end)  PLUMEScold_large  PLUMEScold_large(:,1:nx_add)];
                    PLUMEShot_large_nr(:,:)     = bwlabel(PLUMEShot_large(:,:), 8); %number plumes
                    PLUMEScold_large_nr(:,:)    = bwlabel(PLUMEScold_large(:,:), 8);
                    PLUMEShot_large_nr          = PLUMEShot_large_nr(:,(nx_add+1):end-nx_add,:); %decrease size
                    PLUMEScold_large_nr         = PLUMEScold_large_nr(:,(nx_add+1):end-nx_add,:);
                    PLUMEShot_nr(:,:,i_slice)  	= PLUMEShot_large_nr; %read in matrix to save
                    PLUMEScold_nr(:,:,i_slice)	= PLUMEScold_large_nr;
                    clearvars PLUMEShot_large PLUMEScold_large PLUMEShot_large_nr PLUMEScold_large_nr
                    %                     figure(31),clf
                    %                     subplot(2,1,1); contourfm(theta,phi,PLUMEShot_nr(:,:,i_slice),numContours);
                    %                     subplot(2,1,2); contourfm(theta,phi,PLUMEScold_nr(:,:,i_slice),numContours);
                end
                hold on
                h_hot = contourm(theta,phi,PLUMEShot(:,:,z_num),PLOT.numContours,'LineColor',PLOT.plumeColorHot);
                hold on
                h_cold = contourm(theta,phi,PLUMEScold(:,:,z_num),PLOT.numContours,'LineColor',PLOT.plumeColorCold);
            end
            %     shading interp
            
            %% COLORBAR
            if (PLOT.yyNumSlices>1 && i_slice>1) && CBAR.constant
                PLOT.noColorbar = true;
            else
                PLOT.noColorbar = false;
            end
            DESVARIA.Task    = 'make colorbar';
            [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            try
                if strcmp(get(PLOT.cb,'type'),'colorbar')
                    if (PLOT.yyNumSlices>1 && i_slice==1) && ~CBAR.constant %for multiple slices
                        warning('This needs to be updated.')
                        p       = get(gca,'Position');
                        cbPos   = get(PLOT.cb,'Position');
                        if strcmp(PLOT.cbPos,'bottom')
                            %PROBABLY DOESN'T WORK........
                            set(PLOT.cb,'Location','southoutside','Position',[0.4 0.04 0.2 0.017]);
                        elseif strcmp(PLOT.cbPos,'right')
                            set(PLOT.cb,'Location','eastoutside','Position',[cbPos(1,1)+0.07 cbPos(1,2) cbPos(1,3) cbPos(1,4)]);
                        elseif strcmp(PLOT.cbPos,'southeast')
                            %PROBABLY DOESN'T WORK......
                            set(PLOT.cb,'Location','southoutside','Position',[0.73 0.12 0.17 0.014]);
                        end
                    end
                end
            catch
                %nothing to do
            end
            if CBAR.constant; set(gca,'clim',[FIELD.cColorbarMIN,FIELD.cColorbarMAX]); end
            
            %% COLORMAP
            PLOT.hax = haxm; %<<<<<<<<<<< some of it is done later on, but only for last axis...
            [SWITCH,~] = f_DesignColourmap(SWITCH,PLOT,FIELD,[],[]);
            
            %% TITLE
            if i_slice==PLOT.yyNumSlices
                DESVARIA.Task    = 'make title';
                [PLOT,DESVARIA,STYLE,SWITCH] = f_DesignVaria(DESVARIA,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE);
            end
            
            %% TEXTING
            if SWITCH.Texting %display min/max values
                disp(['   ',FIELD.name,' variation: ',num2str(min(min(amap(:,:,i_slice))),4),' ',STYLE.SCHAR.doubleSidedArrow,' ',num2str(max(max(amap(:,:,i_slice))),4),' ',FIELD.dim])
            end
            
            %% DEPTH ANNOTATION
            if PLOT.yyAnnDepth && strcmp(FIELD.plotType,'field')
                tx = text(1.75, -1.2,['z = ',num2str(z_level_dim,4),' ',GRID.Zdim]); %top right
                set(tx,'FontName',STYLE.AllFontName,'FontSize',STYLE.keyFontSize,'Color',STYLE.annotationColor);
                %create text background
                set(gca,'XLimMode', 'manual', 'YLimMode', 'manual'); % prevent the axes from resizing automatically ???
                p = get(tx,'Extent'); %Get the outer position of the text
                % create a patch around the text object
                p(1) = p(1)/2;
                p(3) = p(3)+p(1)+0.1;
                pObj = patch([p(1) p(1) p(1)+p(3) p(1)+p(3)], [p(2) p(2)+p(4) p(2)+p(4) p(2)],...
                    STYLE.annotationBackColor,'EdgeColor','none',...
                    'FaceAlpha',0.8,'Clipping','off');
                uistack(pObj,'bottom')
                clearvars tx
            end
            
            %% PATCH AREA
            if logical(0)
                % now create a  patch at the bottom
                p(1) = 0;           p(3) = 2*pi;
                p(2) = -pi/2;       p(4) = pi;
                pObj = patch([p(1) p(1) p(1)+p(3) p(1)+p(3)], [p(2) p(2)+p(4) p(2)+p(4) p(2)],...
                    [0.85 0.85 0.85],'EdgeColor','none',...
                    'FaceAlpha',0.8);
                uistack(pObj,'bottom')
            end
            
            %% SOME STYLING
            set(gca,'Color',STYLE.annotationColor,'FontSize',12)
        end
    end
    
    %% SMALL DEPTH GRAPH
    if PLOT.yyDepthGraph && PLOT.yyNumSlices==1 && strcmp(FIELD.plotType,'field')
        haxD = axes('Position', [0.05, 0.8, 0.04, 0.17]);
        plot([0 1],[z_level_dim z_level_dim],'Color',STYLE.BlackOrWhiteColor,'LineWidth',1.5)
        
        axis ij
        ylabel(['Depth [',GRID.Zdim,']'])
        axis([0 1 min(GRID.Z_3Dp(:)) max(GRID.Z_3Dp(:))])
        box on
        %         set(cla,'XTick',[])  %loescht auch den graph....
        set(gca,'XTickLabel','','XTick',[],'xColor','w','yColor',[0.5 0.5 0.5])
        axes(PLOT.hax); %re-activate main axis
    end
    
    %% SAVING PLUME DATA
    if trackPlumes && SWITCH.savePlumes  %saves output variable
        disp('   ...saving')
        
        SaveDirectory = [FILE.stemSave,FILE.name,filesep,'+plumes',filesep];
        if ~exist(SaveDirectory,'dir'); mkdir(SaveDirectory); end
        Figure_Name = ['op_',FILE.name,'_',num2str(FILE.number)];
        if SWITCH.DimensionalMode || SWITCH.DimensionalInput; Figure_Name = ['opDIM_',FILE.name,'_',num2str(FILE.number)]; end
        dirF = pwd; cd(SaveDirectory)
        save([Figure_Name,'.mat'],'PLUMES')
        cd(dirF)
    end
    
    AXorig  = gca; %axis of big plot
    
else
    error('unknown grid type!')
end %grid type





%_____________________________________________________________________________________________________________________

%% COLORMAP SPECIFICATION
if strcmp(GRID.Type,'yinyang') && PLOT.yyNumSlices>1
    %nothing to be done for multiple yy-slices
else
    AXhandles(1)        = AXorig;  %big plot axis handle
    if exist('haxZ','var')
        AXhandles(2)    = haxZ; %zoom plot axis handle
        numLoopsZoom    = 2;
    else
        numLoopsZoom    = 1;
    end
    for iZoom=1:numLoopsZoom
        if ~SWITCH.multipleColormaps && iZoom>1; continue; end
        PLOT.hax        = AXhandles(iZoom);
        
        [SWITCH,~] = f_DesignColourmap(SWITCH,PLOT,FIELD,[],[]);
    end
end

%% UNIVERSAL SWITCHES
if SWITCH.GridAlwaysOn
    grid on
    box on
end

%% GRID PLOT ADDITION
if ~strcmp(GRID.Type,'yinyang') %this is all done above for yy-cases
    if (grid_T && strcmp(FIELD.name,'Temperature')) || ...
            (grid_eta && strcmp(FIELD.name,'Viscosity')) || ...
            (grid_str && strcmp(FIELD.name,'Stress')) || ...
            (grid_edot && strcmp(FIELD.name,'Strain rate'))
                
        numLoopsZoom = 1;
        if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
            numLoopsZoom = 2;
        end
        for iZoom=1:numLoopsZoom
            if iZoom==2; axes(haxZ); end
            %add to plot
            figure(1)
            hold on
            plot(GRID.x2dp,GRID.z2dp,'k',GRID.x2dp',GRID.z2dp','k','LineWidth',0.25)
            hold off
            if logical(0) %separate plot
                figure(3)
                plot(GRID.x2dp,GRID.z2dp,'k',GRID.x2dp',GRID.z2dp','k')
                xlabel(PLOT.xlabel); ylabel(PLOT.ylabel);
                if SWITCH.AxesEqual; axis equal; end
                axis tight; axis ij
                figure(1) %return to fig 1
            end
            if numLoopsZoom==2; axes(AXcurrent); end
        end
    end
end

%% AXIS LIMIT BIG PLOT
if ~strcmp(GRID.Type,'yinyang') %this is all done above for yy-cases
    if SWITCH.AxesLimit
        if size(SWITCH.AxesLimitValues,1)>1  %individual axes values
            axis(SWITCH.AxesLimitValues(PLOT.loopCase,:));
        else %all plots the same
            axis(SWITCH.AxesLimitValues);
        end
    end
end
%_____________________________________________________________________________________________________________________





%==========================================================================
%% PLUME CONTOURS
%==========================================================================
if ~strcmp(GRID.Type,'yinyang') && ~strcmp(GRID.Dim,'3-D')
    if PLOT.indicatePlumes || ...
            (strcmp(FIELD.name,'Temperature') && plume_T) || ...
            (strcmp(FIELD.name,'Viscosity') && plume_eta) || ...
            (strcmp(FIELD.name,'Stress') && plume_str) || ...
            ((strcmp(FIELD.name,'Velocity') || strcmp(FIELD.name,'Horizontal velocity') || strcmp(FIELD.name,'Radial velocity')) ...
            && plume_v) || ...
            (strcmp(FIELD.name,'Strain rate') && plume_edot) || ...
            (strcmp(FIELD.name,'Density') && plume_rho) % plot tracked plumes
        
        % SETTING UP CORRECT GRID VARIABLES AND SIZES
        if strcmp(GRID.Type,'Cartesian')
            x2d_dummy  = x2dp; z2d_dummy = z2dp;
            if strcmp(GRID.Dim,'2-D')
                dummy      = zeros(size(MANTLE.plumesHot,1),size(MANTLE.plumesHot,3)); %adjust matrix size for later contour plot
                dummy(:,:) = MANTLE.plumesHot(:,1,:); plumesHot = dummy;
                dummy(:,:) = MANTLE.plumesCold(:,1,:); plumesCold = dummy;
            end
        elseif strcmp(GRID.Type,'spherical2D')
            x2d_dummy = GRID.x2ds; z2d_dummy = GRID.z2ds;
            if SWITCH.closeAnnulus %add a row to close spherical annulus plot
                MANTLE.plumesHot(end+1,:,:) = MANTLE.plumesHot(1,:,:);
                MANTLE.plumesCold(end+1,:,:) = MANTLE.plumesCold(1,:,:);
            end
            plumesHot(:,:)  = MANTLE.plumesHot(:,1,:);
            plumesCold(:,:) = MANTLE.plumesCold(:,1,:);
        end
        
        % PLOT ADDITION TO PLOT
        numLoopsZoom = 1;
        if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
            numLoopsZoom = 2;
        end
        for iZoom=1:numLoopsZoom
            if iZoom==2; axes(haxZ); end
            if MANTLE.numHotPlumes>0
                hold on
                h_phot = contour(x2d_dummy,z2d_dummy,plumesHot,PLOT.numContours,...
                    'LineColor',PLOT.plumeColorHot,'LineWidth',PLOT.streamWidth);
            end
            if MANTLE.numColdPlumes>0
                hold on
                h_pcold = contour(x2d_dummy,z2d_dummy,plumesCold,PLOT.numContours,...
                    'LineColor',PLOT.plumeColorCold,'LineWidth',PLOT.streamWidth);
            end
            hold off
            if numLoopsZoom==2; axes(AXcurrent); end
        end
        if strcmp(GRID.Type,'spherical2D') && SWITCH.closeAnnulus %remove row again
            MANTLE.plumesHot(end,:,:) = []; MANTLE.plumesCold(end,:,:) = [];
        end
    end
end


%==========================================================================
%% STREAMFUNCTION
%==========================================================================
if (strcmp(FIELD.name,'Temperature') && streamfun_T) || ...
        (strcmp(FIELD.name,'Viscosity') && streamfun_eta) || ...
        (strcmp(FIELD.name,'Stress') && streamfun_str) || ...
        (strcmp(FIELD.name,'Velocity') && streamfun_v) || ...
        (strcmp(FIELD.name,'Horizontal velocity') && streamfun_vx) || ...
        (strcmp(FIELD.name,'Radial velocity') && streamfun_vr) || ...
        (strcmp(FIELD.name,'Strain rate') && streamfun_edot) || ...
        (strcmp(FIELD.name,'Upwelling and downwelling') && streamfun_udw) || ...
        (strcmp(FIELD.name,'Density') && streamfun_rho) && ...
        ~strcmp(GRID.Type,'yinyang') % plot streamfunction
    
    if strcmp(GRID.Dim,'3-D')
        warning('fc: 3-D streamfunction not implemented! - Switch it off!')
    else %2-D
        if (isfield(PLOT,'VX_3D') && isfield(PLOT,'VY_3D') && isfield(PLOT,'VZ_3D'))
            VX_3D = PLOT.VX_3D; VY_3D = PLOT.VY_3D; VZ_3D = PLOT.VZ_3D;
        else
            %Read velocity information
            [~,~,~,~,VX_3D,VY_3D,VZ_3D,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Velocity',SWITCH);
            if ischar(VX_3D); error(VX_3D); end
            if size(VX_3D,1)==1 %exchange x and y
                dummy_3D = zeros(size(VX_3D,2),size(VX_3D,1),size(VX_3D,3));
                dummy_3Dx = zeros(size(VX_3D,2),size(VX_3D,1),size(VX_3D,3));
                dummy_3D(:,1,:)     = VZ_3D(1,:,:); VZ_3D = dummy_3D;
                dummy_3Dx(:,1,:)    = VX_3D(1,:,:);
                dummy_3D(:,1,:)     = VY_3D(1,:,:); VY_3D = dummy_3Dx; VX_3D = dummy_3D;
            end
        end
        if strcmp(GRID.Dim,'3-D')
            %not implemented for 3-D........
        elseif strcmp(GRID.Dim,'2-D')
            vx2d(:,:)   = VX_3D(:,1,:);
            vy2d(:,:)   = VY_3D(:,1,:);
            vz2d(:,:)   = VZ_3D(:,1,:);
        end
        [phi2d,psi2d] = flowfun(vz2d,vx2d);  % calculate the velocity potential (phi) and streamfunction (psi)
        %grid variables
        x2d_dummy = x2dp;
        z2d_dummy = z2dp;
        if strcmp(GRID.Type,'spherical2D')
            %spherical2D grid variables
            x2d_dummy   = GRID.x2ds;
            z2d_dummy   = GRID.z2ds;
            if SWITCH.closeAnnulus
                phi2d(end+1,:)	= phi2d(1,:);
                psi2d(end+1,:)	= psi2d(1,:);
            end
        end
        %     end
        if strcmp(stream_data,'phi'); var2d = phi2d; elseif strcmp(stream_data,'psi'); var2d = psi2d; else; error('fc: unknown stream data'); end
        
        numLoopsZoom = 1;
        if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
            numLoopsZoom = 2;
        end
        for iZoom=1:numLoopsZoom
            if iZoom==2; axes(haxZ); end
            hold on
            if strcmp(PLOT.streamColoured,'false')
                hstr = contour(x2d_dummy,z2d_dummy,var2d,stream_numContours,...
                    'LineColor',PLOT.streamColor,'LineWidth',PLOT.streamWidth);
                
            else %DOESN'T WORK PROPERLY ANYLONGER.. MIGHT REMOVE.
                caxis('manual');  %fix the colorbar scale (old style)
                for iiZoom=1:numLoopsZoom
                    if iiZoom==1; axes(AXcurrent); else; axes(haxZ); end
                    freezeColors;
                    if iZoom==2; axes(haxZ); end
                end
                %if self-made colorbar:
                if size(PLOT.streamColoured,1)==1 && exist([PLOT.streamColoured,'.mat'],'file')==2
                    load([PLOT.streamColoured,'.mat']); PLOT.streamColoured=eval( sprintf(PLOT.streamColoured) ); %load colorbar
                end
                %adjust colormap for streamfunction
                dummy = SWITCH.ColormapAll;
                SWITCH.ColormapAll = {PLOT.streamColoured 'noflip'};
                [SWITCH,~] = f_DesignColourmap(SWITCH,PLOT,FIELD,[],[]);
                SWITCH.ColormapAll = dummy;
                %colormap(PLOT.streamColoured);
                contour(x2d_dummy,z2d_dummy,var2d,stream_numContours);
                
                cval_max=max( abs(min(min(var2d))),abs(max(max(var2d))) );
                if CBAR.constant; set(gca,'clim',[-cval_max,cval_max]); end
                %colorbar
                for iiZoom=1:numLoopsZoom
                    if iiZoom==1; axes(AXcurrent); else; axes(haxZ); end
                    freezeColors;
                    if iZoom==2; axes(haxZ); end
                end
                %back to original:
                if strcmp(SWITCH.ColormapAll{1,2},'flip'); colormap(flipud(SWITCH.ColormapAll{1,1}));
                elseif strcmp(SWITCH.ColormapAll{1,2},'noflip'); colormap(SWITCH.ColormapAll{1,1}); end
                if SWITCH.Colorbar(1,1)
                    if iZoom==1
                        PLOT.cb = colorbar;
                        title(PLOT.cb,FIELD.dim); %on top
                        %ylabel(PLOT.cb,FIELD.dim); %on left hand side
                    end
                end
                if CBAR.constant; set(gca,'clim',[FIELD.cColorbarMIN,FIELD.cColorbarMAX]); end
            end
            hold off
            if numLoopsZoom==2; axes(AXcurrent); end
        end
    end
end



%==========================================================================
%% STREAM LINES
%==========================================================================
if streamlinePlot || ...
        (strcmp(FIELD.name,'Temperature') && streamline_T) || ...
        (strcmp(FIELD.name,'Viscosity') && streamline_eta) || ...
        (strcmp(FIELD.name,'Stress') && streamline_str) || ...
        (strcmp(FIELD.name,'Velocity') && streamline_v) || ...
        (strcmp(FIELD.name,'Horizontal velocity') && streamline_vh) || ...
        (strcmp(FIELD.name,'Radial velocity') && streamline_vr) && ...
        ~strcmp(GRID.Type,'yinyang') % plot streamlines
    
    VARIA.Task    = 'create stream line';
    [PLOT,VARIA] = f_Varia(VARIA,FILE,FIELD,GRID,SWITCH,PLOT,STYLE,SAVE,SETUP);

    if strcmp(GRID.Dim,'3-D')
        hold on
        hs = streamline(GRID.y2dp(1,:,1),GRID.x2dp(:,1,1),GRID.z2dp(1,1,:),vy2d',vx2d',vz2d',sy2d,sx2d,sz2d,[VARIA.stepsize,VARIA.maxNumVertices]);
        set(hs,'Color',PLOT.streamColor)  %,'LineWidth',1
        %needs to be adjusted..............
    elseif strcmp(GRID.Dim,'2-D')  % 2-D
        hold on
        %             hs      = streamline(VARIA.x2d(:,1),VARIA.z2d(1,:),VARIA.vx2d',VARIA.vz2d',VARIA.sxLine,VARIA.szLine,[VARIA.stepsize,VARIA.maxNumVertices]);
        for iL=1:size(VARIA.streamData,2)
            startLength     = min(length(VARIA.streamData{1,iL}(:,1)), 30/VARIA.stepsize);  %<<<<< set length
            endingLength    = min(length(VARIA.streamData{1,iL}(:,1))-1, 20/VARIA.stepsize);  %<<<<< set length
            endingLength    = min(length(VARIA.streamData{1,iL}(:,1))-startLength, endingLength);
            hL(iL)          = line(VARIA.streamData{1,iL}(1:startLength,1),VARIA.streamData{1,iL}(1:startLength,2),'Color',PLOT.streamColor,'LineWidth',0.1); %start
            hL(iL)          = line(VARIA.streamData{1,iL}(ceil(3/5*startLength):startLength,1),VARIA.streamData{1,iL}(ceil(3/5*startLength):startLength,2),'Color',PLOT.streamColor,'LineWidth',0.25); %start2
            hL(iL)          = line(VARIA.streamData{1,iL}(ceil(4/5*startLength):startLength,1),VARIA.streamData{1,iL}(ceil(4/5*startLength):startLength,2),'Color',PLOT.streamColor,'LineWidth',0.35); %start3
            hL(iL)          = line(VARIA.streamData{1,iL}(startLength:end-endingLength,1),VARIA.streamData{1,iL}(startLength:end-endingLength,2),'Color',PLOT.streamColor,'LineWidth',0.5); %middle
            hL(iL)          = line(VARIA.streamData{1,iL}(end-endingLength:end,1),VARIA.streamData{1,iL}(end-endingLength:end,2),'Color',PLOT.streamColor,'LineWidth',0.65); %end1
            hL(iL)          = line(VARIA.streamData{1,iL}(end-ceil(endingLength*6/8):end,1),VARIA.streamData{1,iL}(end-ceil(endingLength*6/8):end,2),'Color',PLOT.streamColor,'LineWidth',0.75); %end2
            hL(iL)          = line(VARIA.streamData{1,iL}(end-ceil(endingLength*5/8):end,1),VARIA.streamData{1,iL}(end-ceil(endingLength*5/8):end,2),'Color',PLOT.streamColor,'LineWidth',0.85); %end3
            hL(iL)          = line(VARIA.streamData{1,iL}(end-ceil(endingLength/2):end,1),VARIA.streamData{1,iL}(end-ceil(endingLength/2):end,2),'Color',PLOT.streamColor,'LineWidth',1.0); %end4
        end
    end
    hold off
    if streamlinePlot && SP.number==1
        figure(4); hold on;
        hs = streamline(VARIA.x2d(:,1),VARIA.z2d(1,:),VARIA.vx2d',VARIA.vz2d',VARIA.sxLine,VARIA.szLine,[VARIA.stepsize,VARIA.maxNumVertices]);
        set(hs,'Color',PLOT.streamColor)
        if SWITCH.AxesEqual; axis equal; end
        if strcmp(GRID.Type,'Cartesian')
            axis tight
            %---------------------
            axis ij  % reverse axis (2.)
            %---------------------
            xlabel(PLOT.xlabel)
            ylabel(PLOT.ylabel)
        end
        if PLOT.time_dim/1e6>100
            title(['\bf{','Streamline }\rm{ ',num2str(PLOT.time_dim/1e9,3),' Gyr'])
        else
            title(['\bf{','Streamline }\rm{ ',num2str(PLOT.time_dim/1e6,3),' Myr'])
        end
        figure(1)
    end
end


%==========================================================================
%% DEFORMATION MECHANISM
%==========================================================================
if (strcmp(FIELD.name,'Temperature') && defmech_T) || ...
        (strcmp(FIELD.name,'Viscosity') && defmech_eta) || ...
        (strcmp(FIELD.name,'Stress') && defmech_str) || ...
        (strcmp(FIELD.name,'Strain rate') && defmech_edot) && ...
        ~strcmp(GRID.Type,'yinyang') % plot deformation mechanism
    
    % Read pressure & velocity information
    [~,~,~,~,VAR_3D,~,~]     = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Deformation mechanism',SWITCH);
    if ischar(VAR_3D); error(VAR_3D); end
    if size(VAR_3D,1)==1 %exchange x and y
        dummy_3D = zeros(size(VAR_3D,2),size(VAR_3D,1),size(VAR_3D,3));
        dummy_3D(:,1,:) = VAR_3D(1,:,:);    VAR_3D = dummy_3D;
    end
    if strcmp(GRID.Dim,'3-D')
    elseif strcmp(GRID.Dim,'2-D')
        if strcmp(GRID.Type,'spherical2D') && SWITCH.closeAnnulus %add row for nicer visualisation
            VAR_3D(end+1,:,:) = VAR_3D(1,:,:);
        end
        defm2d(:,:)  = VAR_3D(:,1,:);
        %make sure values are the exact numbers they should be down to 0.00x:
        defm2d(:,:)  = round(defm2d*1e3)/1e3;
        
        if defm2d(2,2)==1 && defm2d(3,2)==4 %check for artificial values
            defm2d(2,2) = defm2d(4,2);  %and remove them
            defm2d(3,2) = defm2d(4,2);
        end
    end
    numLoopsZoom = 1;
    if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
        numLoopsZoom = 2;
    end
    for iZoom=1:numLoopsZoom
        if iZoom==2; axes(haxZ); end
        hold on
        
        caxis('manual') %fix the colorbar scale
        for i_cont = 1:size(defmech_contours,2)
            %1: diffusion creep, 2: dislocation creep, 3: plasticity, 4: visc. cutoff
            var2d_defm = zeros(size(defm2d));
            if defmech_contours(i_cont)==1 %diffusion creep
                var2d_defm(defm2d<2) = 1;
            elseif defmech_contours(i_cont)==1.1  %1/10 dislocation creep
                var2d_defm(defm2d>=1.1 & defm2d<=2) = 1;
            elseif defmech_contours(i_cont)==1.3  %1/3 dislocation creep
                var2d_defm(defm2d>=1.3 & defm2d<=2) = 1;
            else  %other
                var2d_defm(defm2d==defmech_contours(i_cont)) = 1;
            end
            if defmech_contours(i_cont)==1 || defmech_contours(i_cont)==2 || ...
                    defmech_contours(i_cont)==3 || defmech_contours(i_cont)==4
                linestyle1 = '-';
            else
                linestyle1 = ':';
            end
            if strcmp(GRID.Type,'spherical2D')
                x2d_dummy = GRID.x2ds; z2d_dummy = GRID.z2ds;
            else
                x2d_dummy = x2dp; z2d_dummy = z2dp;
            end
            contour(x2d_dummy,z2d_dummy,var2d_defm,[0.5 0.5],...
                'LineColor',PLOT.defmechColor(i_cont,:), 'LineStyle',linestyle1);
        end
        hold off
        if numLoopsZoom==2; axes(AXcurrent); end
    end
    clearvars defm2d var2d_defm
end



%==========================================================================
%% FIELD CONTOUR
%==========================================================================
%ANY FIELD CONTOUR PLOT ADDITION
if (fieldContour_C && strcmp(FIELD.name,'Composition')) || ...
        (fieldContour_T && strcmp(FIELD.name,'Temperature')) || ...
        (fieldContour_eta && strcmp(FIELD.name,'Viscosity')) || ...
        (fieldContour_nstr && strcmp(FIELD.name,'zz-Stress component')) || ...
        (fieldContour_str && strcmp(FIELD.name,'Stress')) || ...
        (fieldContour_edot && strcmp(FIELD.name,'Strain rate')) || ...
        (fieldContour_vdiss && strcmp(FIELD.name,'Viscous dissipation')) || ...
        (fieldContour_stream && strcmp(FIELD.name,'Streamfunction')) || ...
        (fieldContour_udw && strcmp(FIELD.name,'Upwelling and downwelling')) || ...
        (fieldContour_v && strcmp(FIELD.name,'Velocity')) || ...
        (fieldContour_vx && strcmp(FIELD.name,'Horizontal velocity')) || ...
        (fieldContour_vr && strcmp(FIELD.name,'Radial velocity')) && ...
        ~strcmp(GRID.Type,'yinyang')
        
    numLoopsZoom = 1;
    if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
        numLoopsZoom = 2;
    end
    for iZoom=1:numLoopsZoom
        if iZoom==2; axes(haxZ); end
        %read appropriate field
        [~,~,~,~,VAR_3D,~,~]     = f_readStagYY(FILE.directory,FILE.name,FILE.number,PLOT.fContourField,SWITCH);
        if ischar(VAR_3D); error(VAR_3D); end
        if size(VAR_3D,1)==1 %exchange x and y
            dummy_3D = zeros(size(VAR_3D,2),size(VAR_3D,1),size(VAR_3D,3));
            dummy_3D(:,1,:) = VAR_3D(1,:,:); VAR_3D = dummy_3D;
        end
        if strcmp(GRID.Dim,'2-D')
            var2dF(:,:,:)  = VAR_3D(:,:,:);
            if strcmp(GRID.Type,'spherical2D') && SWITCH.closeAnnulus %add row for nicer visualisation
                var2dF(end+1,:,:) = var2dF(1,:,:);
            end
        elseif strcmp(GRID.Dim,'3-D')
            error('field contour for 3-D not implemented yet!!!!')
        end
        var2dF = var2dF *FIELD_M{strcmp(FIELD_M(:,2),PLOT.fContourField),5}; %dimensionalize
        
        %find values below, at & above the contour value
        dummy   = var2dF;
        dummy(var2dF<fieldContour_isoVal) = 0; %below
        dummy(var2dF==fieldContour_isoVal) = 1; %at
        dummy(var2dF>fieldContour_isoVal) = 2; %above
        var2dF  = dummy;
        
        %find contour (boundary values are not checked!)
        jj=1; iii=0;
        fieldContourLoc = zeros(size(var2dF,1)*size(var2dF,2)*size(var2dF,3)/2,3)*NaN; %max. possible size/2
        for ii=2:size(var2dF,1)-1
            %             for jj=2:size(var2dF,2)-1
            for kk=2:size(var2dF,3)-1
                if var2dF(ii,jj,kk)==0 && ... %just below contour
                        ( var2dF(ii+1,jj,kk)==2 || var2dF(ii,jj,kk+1)==2 || ...
                        var2dF(ii-1,jj,kk)==2 || var2dF(ii,jj,kk-1)==2 )
                    var2dF(ii,jj,kk) = 1;
                end
                if var2dF(ii,jj,kk)==1 %create contour points array
                    iii = iii+1;
                    if strcmp(GRID.Type,'Cartesian')
                        fieldContourLoc(iii,:) = [GRID.X_3D(ii,jj,kk) GRID.Y_3D(ii,jj,kk) GRID.Z_3D(ii,jj,kk)];
                    elseif strcmp(GRID.Type,'spherical2D')
                        fieldContourLoc(iii,:) = [GRID.x2ds(ii,kk) 1 GRID.z2ds(ii,kk)];
                    end
                end
            end
            %             end
        end
        
        % FIND PLATE THICKNESS AND SLAB TIP
        findSlabTip = logical(0);
        findPlateThickness = logical(0);
        if findPlateThickness || findSlabTip
            if strcmp(GRID.Type,'spherical2D'); error('might not work for spherical2D grid, yet: Check!'); end
            fieldContourLoc(isnan(fieldContourLoc(:,1)),:) = []; %in [m]
            
            % PLATE THICKNESS
            % find depth value slightly deeper than LAB
            plateContourA = fieldContourLoc(fieldContourLoc(:,3)<=mean(fieldContourLoc(:,3)),:); %remove most of slab values
            plateContourA(plateContourA<0) = 0; %set negative values to zero to use geomean
            
            %         plateThicknessMean1 = harmmean(plateContourA(:,3)); %find most probable plate thickness harmonic mean
            plateThicknessMean2 = geomean(plateContourA(:,3));
            if PLOT.loopField==1; disp(['   > mean plate thickness is ',num2str(plateThicknessMean2/1e3),' km']); end
            %         plateThicknessMean3 = mean(plateContourA(:,3));
            
            %         figure(2)
            %         plot(fieldContourLoc(:,3)/1e3,'.')
            %         hold on
            %         plot(plateContourA(:,3)/1e3,'.r')
            %         hold on
            %         plot(plateThicknessMean1/1e3*ones(size(fieldContourLoc(:,3))),'k')
            %         hold on
            %         plot(plateThicknessMean2/1e3*ones(size(fieldContourLoc(:,3))),'g')
            %         hold on
            %         plot(plateThicknessMean3/1e3*ones(size(fieldContourLoc(:,3))),'y')
            
            % SLAB TIP
            % find depth value at slab tip
            %             slabtipContourA = fieldContourLoc(fieldContourLoc(:,3)<=mean(fieldContourLoc(:,3)),:); %derive slab values
            
            %sort for depth (deepest first):
            [slabtipContourA,sortIndex] = sort(fieldContourLoc(:,3),'descend');
            %extract only the x% points that lie deepest
            x_percent = 0.05;
            slabtipContourA = fieldContourLoc(sortIndex(1:round(size(sortIndex,1)*x_percent),1),:);
            % geometric mean of the slab tip area
            slabtipCentreA  = [geomean(slabtipContourA(:,1))    geomean(slabtipContourA(:,3))];
            slabtipBottomA  = [slabtipContourA(1,1)             slabtipContourA(1,3)];
            slabtipLeftA    = [min(slabtipContourA(:,1))        max(slabtipContourA(slabtipContourA(:,1)==(min(slabtipContourA(:,1))),3))];
            slabtipRightA   = [max(slabtipContourA(:,1))        max(slabtipContourA(slabtipContourA(:,1)==(max(slabtipContourA(:,1))),3))];
            
            %             figure(2)
            %             plot(fieldContourLoc(:,1),fieldContourLoc(:,3),'.')
            %             hold on
            %             plot(slabtipContourA(:,1),slabtipContourA(:,3),'.r')
            %             plot(slabtipCentreA(:,1),slabtipCentreA(:,2),'ro')
            %             plot(slabtipBottomA(:,1),slabtipBottomA(:,2),'go')
            %             plot(slabtipLeftA(:,1),slabtipLeftA(:,2),'go')
            %             plot(slabtipRightA(:,1),slabtipRightA(:,2),'ro')
            
            defineSlabTipAutomatically = logical(1);
            if defineSlabTipAutomatically
%                 figure(1)
                hold on
                %!! ASSUMPTION !!: if x-left is further away from x-centre than x-right - then right hand side subduction
                if abs(slabtipCentreA(:,1)-slabtipLeftA(:,1)) == abs(slabtipCentreA(:,1)-slabtipRightA(:,1)) %DOWN
                    if PLOT.loopField==1; disp('   > Vertical subduction detected'); end
                    slabtipLocA = slabtipLeftA;
                    if PLOT.loopField==1; disp(['   >    slab tip bottom: [',num2str(slabtipBottomA(:,1)/1e3,4),',',num2str(slabtipBottomA(:,2)/1e3,4),'] km']); end
                elseif abs(slabtipCentreA(:,1)-slabtipLeftA(:,1)) > abs(slabtipCentreA(:,1)-slabtipRightA(:,1)) %RIGHT
                    if PLOT.loopField==1; disp('   > Right-hand side subduction detected'); end
                    slabtipLocA = slabtipRightA;
                    if PLOT.loopField==1; disp(['   >    slab tip right:  [',num2str(slabtipRightA(:,1)/1e3,4),',',num2str(slabtipRightA(:,2)/1e3,4),'] km']); end
                else                                                                                            %LEFT
                    if PLOT.loopField==1; disp('   > Left-hand side subduction detected'); end
                    slabtipLocA = slabtipLeftA;
                    if PLOT.loopField==1; disp(['   >    slab tip left:   [',num2str(slabtipLeftA(:,1)/1e3,4),',',num2str(slabtipLeftA(:,2)/1e3,4),'] km']); end
                end
                plot(slabtipLocA(:,1)/1e3,slabtipLocA(:,2)/1e3,'x','Color',TOPO.lineColor,...
                    'MarkerSize',6,'MarkerFaceColor',TOPO.lineColor);
                hold off
                
                %write to data file
                SAVE.datafile = logical(0);
                if SAVE.datafile && PLOT.loopField==1 %only once per time step
                    SAVE.DataName               = 'slabTip_timeXZ';
                    SAVE.Directory              = [FILE.stemSave,FILE.name,'/'];
                    SAVE.data                   = [FILE.number PLOT.time slabtipLocA]; % output file - time - x - z
                    SAVE.write2existing         = logical(0);
                    SAVE.txt                    = logical(1);
                    SAVE.mat                    = logical(1);
                    
                    [SAVE.overwriteAll] = f_saveData( SAVE );
                end
                
            else
                if PLOT.loopField==1; disp(['   > slab tip centre: [',num2str(slabtipCentreA(:,1)/1e3,4),',',num2str(slabtipCentreA(:,2)/1e3,4),'] km']); end
                if PLOT.loopField==1; disp(['   > slab tip bottom: [',num2str(slabtipBottomA(:,1)/1e3,4),',',num2str(slabtipBottomA(:,2)/1e3,4),'] km']); end
                if PLOT.loopField==1; disp(['   > slab tip left:   [',num2str(slabtipLeftA(:,1)/1e3,4),',',num2str(slabtipLeftA(:,2)/1e3,4),'] km']); end
                if PLOT.loopField==1; disp(['   > slab tip right:  [',num2str(slabtipRightA(:,1)/1e3,4),',',num2str(slabtipRightA(:,2)/1e3,4),'] km']); end
%                 figure(1)
                hold on
                plot(slabtipCentreA(:,1)/1e3,slabtipCentreA(:,2)/1e3,'x','Color',PLOT.fieldContourColor,...
                    'MarkerSize',6,'MarkerFaceColor',PLOT.fieldContourColor);
                plot(slabtipBottomA(:,1)/1e3,slabtipBottomA(:,2)/1e3,'x','Color',PLOT.fieldContourColor,...
                    'MarkerSize',6,'MarkerFaceColor',PLOT.fieldContourColor);
                plot(slabtipLeftA(:,1)/1e3,slabtipLeftA(:,2)/1e3,'x','Color',PLOT.fieldContourColor,...
                    'MarkerSize',6,'MarkerFaceColor',PLOT.fieldContourColor);
                plot(slabtipRightA(:,1)/1e3,slabtipRightA(:,2)/1e3,'x','Color',PLOT.fieldContourColor,...
                    'MarkerSize',6,'MarkerFaceColor',PLOT.fieldContourColor);
                hold off
            end
        end
        
%         figure(1)
        hold on
        if strcmp(GRID.Dim,'3-D')
            error('fc: Field contour option not yet implemented in 3-D/for spherical2D!')
        elseif strcmp(GRID.Dim,'2-D')
            if strcmp(GRID.Type,'Cartesian')
                plot(GRID.X_3D(var2dF==1)/1e3,GRID.Z_3D(var2dF==1)/1e3,'o','Color',PLOT.fieldContourColor,...
                    'MarkerSize',PLOT.fieldContourSize,'MarkerFaceColor',PLOT.fieldContourColor,'LineWidth',TOPO.lineWidth);
            elseif strcmp(GRID.Type,'spherical2D')
                plot(GRID.x2ds(var2dF(:,1,:)==1),GRID.z2ds(var2dF(:,1,:)==1),'o','Color',PLOT.fieldContourColor,...
                    'MarkerSize',PLOT.fieldContourSize,'MarkerFaceColor',PLOT.fieldContourColor,'LineWidth',TOPO.lineWidth);
            end
        end
        hold off
        if numLoopsZoom==2; axes(AXcurrent); end
    end
    clearvars var2dF z2dF in
end



%==========================================================================
%% QUIVER
%==========================================================================
if quiver_plot || ...
        (strcmp(FIELD.name,'Temperature') && quiver_T) || ...
        (strcmp(FIELD.name,'Viscosity') && quiver_eta) || ...
        (strcmp(FIELD.name,'Stress') && quiver_str) || ...
        (strcmp(FIELD.name,'Streamfunction') && quiver_psi) || ...
        (strcmp(FIELD.name,'zz-Stress component') && quiver_nstr) || ...
        (strcmp(FIELD.name,'Strain rate') && quiver_edot) || ...
        (strcmp(FIELD.name,'Velocity') && quiver_v) || ...
        (strcmp(FIELD.name,'Horizontal velocity') && quiver_vx) || ...
        (strcmp(FIELD.name,'Radial velocity') && quiver_vr) || ...
        (strcmp(FIELD.name,'Density') && quiver_rho) && ...
        ~strcmp(GRID.Type,'yinyang')% plot velocity arrows
    
    if plate_velocity
        quiver_numDel	= 6;
        scale           = 20;
    end
    if isfield(PLOT,'VX_3D') && isfield(PLOT,'VY_3D') && isfield(PLOT,'VZ_3D')
        if strcmp(GRID.Type,'spherical2D')
            vx2d(:,:) = PLOT.VX_3Ds(:,1,:); vz2d(:,:) = PLOT.VZ_3Ds(:,1,:); vy2d = zeros(nx,nz); %JUST TO BE DEFINED
            VX_S(:,:) = PLOT.VX_3D(:,1,:); VZ_S(:,:) = PLOT.VZ_3D(:,1,:); %this is for spherical2D-Cartesian plot
        elseif  strcmp(GRID.Type,'Cartesian')
            VX_3D = PLOT.VX_3D; VY_3D = PLOT.VY_3D; VZ_3D = -PLOT.VZ_3D; %adjustment for flip of depth vector
        elseif  strcmp(GRID.Type,'yinyang')
            %...
        end
    else
        % Read pressure & velocity information
        [~,~,~,~,VX_3D,VY_3D,VZ_3D,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Velocity',SWITCH);
        if ischar(VX_3D); error(VX_3D); end
        if size(VX_3D,1)==1 %exchange x and y
            dummy_3D = zeros(size(Y_3D,2),size(Y_3D,1),size(Y_3D,3));
            dummy_3Dx = zeros(size(Y_3D,2),size(Y_3D,1),size(Y_3D,3));
            dummy_3D(:,1,:)     = VZ_3D(1,:,:); VZ_3D = dummy_3D;
            dummy_3Dx(:,1,:)    = VX_3D(1,:,:);
            dummy_3D(:,1,:)     = VY_3D(1,:,:); VY_3D = dummy_3Dx; VX_3D = dummy_3D;
        end
        if strcmp(GRID.Type,'spherical2D')
            vx2d        = zeros(nx,nz);
            vy2d        = zeros(nx,nz); %JUST TO BE DEFINED
            vz2d        = zeros(nx,nz);
            vx2d(:,:)   = VX_3D(1,:,:);
            vz2d(:,:)   = VZ_3D(1,:,:);
            % conversion to spherical2D grid
            VX_S = vx2d; VZ_S = vz2d; %this is for spherical2D-Cartesian plot
            vx2d        = -VX_S.*cos(x2d_nd)+VZ_S.*sin(x2d_nd);
            vz2d        = VX_S.*sin(x2d_nd)-VZ_S.*cos(x2d_nd);
        elseif  strcmp(GRID.Type,'Cartesian')
            VZ_3D       = -VZ_3D; %adjustment for flip of depth vector
        end
    end
    if strcmp(GRID.Type,'Cartesian')
        X_3D        = GRID.X_3Dp;
        Y_3D        = GRID.Y_3Dp;
        Z_3D        = GRID.Z_3Dp;
        if SWITCH.ReduceSize
            wb = waitbar(0,'Slimming Matrix (quiver)...');
            %reduce size of matrix
            maxGridpoints = 500000;  %500000;   %<--------------------------=================
            if strcmp(GRID.Dim,'3-D')
                maxGridpoints=maxGridpoints/100;
            end
            nxEff   = nx;
            nyEff   = ny;
            nzEff   = nz;
            while (nxEff*nyEff*nzEff)>maxGridpoints
                waitbar(maxGridpoints/(nxEff*nyEff*nzEff),wb)
                for i=(nxEff-2):-2:2
                    X_3D(i,:,:)     = []; Y_3D(i,:,:)     = []; Z_3D(i,:,:)     = [];
                    VX_3D(i,:,:)   	= []; VY_3D(i,:,:) 	= []; VZ_3D(i,:,:) 	= [];
                end
                for j=(nyEff-2):-2:2
                    X_3D(:,j,:)     = []; Y_3D(:,j,:)     = []; Z_3D(:,j,:)     = [];
                    VX_3D(:,j,:)   	= []; VY_3D(:,j,:) 	= []; VZ_3D(:,j,:) 	= [];
                end
                for k=(nzEff-2):-2:2
                    X_3D(:,:,k)     = []; Y_3D(:,:,k)     = []; Z_3D(:,:,k)     = [];
                    VX_3D(:,:,k)   	= []; VY_3D(:,:,k) 	= []; VZ_3D(:,:,k) 	= [];
                end
                X_3D(4,:,:) = []; X_3D(:,4,:) = []; X_3D(:,:,4) = [];
                Y_3D(4,:,:) = []; Y_3D(:,4,:) = []; Y_3D(:,:,4) = [];
                Z_3D(4,:,:) = []; Z_3D(:,4,:) = []; Z_3D(:,:,4) = [];
                VX_3D(4,:,:) = []; VX_3D(:,4,:) = []; VX_3D(:,:,4) = [];
                VY_3D(4,:,:) = []; VY_3D(:,4,:) = []; VY_3D(:,:,4) = [];
                VZ_3D(4,:,:) = []; VZ_3D(:,4,:) = []; VZ_3D(:,:,4) = [];
                nxEff   = size(X_3D,1);
                nyEff   = size(Y_3D,2);
                nzEff   = size(Z_3D,3);
            end
            disp(['Size of plotted matrix decreased to ',num2str(size(X_3D))])
            clearvars nxEff nyEff nzEff maxGridpoints
            close(wb)
        end
        if strcmp(GRID.Dim,'3-D')
            x2d2 	= X_3D; y2d2 	= Y_3D; z2d2 	= Z_3D;
            vx2d    = VX_3D; vy2d 	= VY_3D; vz2d 	= VZ_3D;
        elseif strcmp(GRID.Dim,'2-D')
            x2d2(:,:) = X_3D(:,1,:); z2d2(:,:) = Z_3D(:,1,:); y2d2 = 1; % just in order to have it defined
            vx2d(:,:) = VX_3D(:,1,:); vy2d(:,:) = VY_3D(:,1,:); vz2d(:,:) = VZ_3D(:,1,:);
        end
    end
    if strcmp(GRID.Dim,'3-D')   % reduce number of velocity arrows  (3-D)
        count=0;
        nn=size(x2d2,1)-1;
        for i=2:1:nn-1
            if count==2
                count=0.;
            else
                x2d2(i,:) = NaN; vx2d(i,:) = NaN;
            end
            count=count+1;
        end
        count=0;
        nn=size(x2d2,2)-1;
        for j=2:1:nn-1
            if count==2
                count=0.;
            else
                y2d2(j,:) = NaN; vy2d(j,:) = NaN;
            end
            count=count+1;
        end
        count=0;
        nn=size(x2d2,3)-1;
        for k=2:1:nn-1
            if count==2
                count=0.;
            else
                z2d2(k,:) = NaN; vz2d(k,:) = NaN;
            end
            count=count+1;
        end
    elseif strcmp(GRID.Dim,'2-D') % reduce number of velocity arrows  (2-D)
        if strcmp(GRID.Type,'spherical2D')
            X_Sreduced = x2d_nd; 
            x2d_Sreduced = GRID.x2ds; z2d_Sreduced = GRID.z2ds;
        elseif strcmp(GRID.Type,'Cartesian')
            %just to have them defined
            X_Sreduced = zeros(size(x2d2));
            x2d_Sreduced = zeros(size(x2d2)); z2d_Sreduced = zeros(size(x2d2));
        end
        %reduce in x-direction
        count=0;
        nn=nx-1;
        if strcmp(GRID.Type,'Cartesian')
            M2reduce=ones(size(x2d2));
        elseif strcmp(GRID.Type,'spherical2D')
            M2reduce=ones(size(X_Sreduced));
        end
        for i=2:1:nn-1
            if count==quiver_numDel
                count=0.;
            else
                M2reduce(i,:) = NaN;
            end
            count=count+1;
        end
        %remove arrows in sticky-air layer
        if strcmp(SETUP.topBC,'sticky-air') && TOPO.area && SWITCH.ActualDepth
            if strcmp(GRID.Type,'Cartesian')
                idx2remove = find(z2d2(1,:)<0);
                x2d2(:,idx2remove) = NaN; y2d2(:,idx2remove) = NaN; z2d2(:,idx2remove) = NaN;
                vx2d(:,idx2remove) = NaN; vy2d(:,idx2remove) = NaN; vz2d(:,idx2remove) = NaN;
            elseif strcmp(GRID.Type,'spherical2D')
                idx2remove = find(z2dp(1,:)<0);
                X_Sreduced(:,idx2remove) = NaN; x2d_Sreduced(:,idx2remove) = NaN; z2d_Sreduced(:,idx2remove) = NaN;
                VX_S(:,idx2remove) = NaN; VZ_S(:,idx2remove) = NaN;
                vx2d(:,idx2remove) = NaN; vy2d(:,idx2remove) = NaN; vz2d(:,idx2remove) = NaN;
            end
        end
        %reduce in z-direction
        count=0;
        nn=nz-1;
        for k=2:1:nn-1
            if count==quiver_numDel
                count=0.;
            else
                M2reduce(:,k) = NaN;
            end
            count=count+1;
        end
        %remove values (i.e., replace them by NaNs)
        if strcmp(GRID.Type,'Cartesian')
            x2d2 = x2d2.*M2reduce; z2d2 = z2d2.*M2reduce;
            vx2d = vx2d.*M2reduce; vz2d = vz2d.*M2reduce;
        elseif strcmp(GRID.Type,'spherical2D')
            if SWITCH.closeAnnulus
                x2d_Sreduced(1,:) = NaN; z2d_Sreduced(1,:) = NaN; %that is the one that gets dublicated
                vx2d(1,:) = NaN; vz2d(1,:) = NaN;
            end
            X_Sreduced = X_Sreduced.*M2reduce;
            x2d_Sreduced(1:size(M2reduce,1),:) = x2d_Sreduced(1:size(M2reduce,1),:).*M2reduce; %also works if SWITCH.closeAnnulus
            z2d_Sreduced(1:size(M2reduce,1),:) = z2d_Sreduced(1:size(M2reduce,1),:).*M2reduce;
            vx2d = vx2d.*M2reduce; vz2d = vz2d.*M2reduce;
            VX_S = VX_S.*M2reduce; VZ_S = VZ_S.*M2reduce;
        end
    end
    if plate_velocity  %CAN BE OPTIMIZED.....
        nnx=size(x2d2,1);
        nny=size(x2d2,2);
        nnz=size(x2d2,3);
        for zz=1:nnz  %loop through z direction
            count=0;
            if z2d2(1,1,zz)>0.0 && z2d2(1,1,zz)<30.0  %keep values
                for xx=2:1:nnx-1
                    if count==6
                        count=0;
                    else
                        count=count+1;
                        x2d2(xx,:,zz) = NaN; y2d2(xx,:,zz) = NaN; z2d2(xx,:,zz) = NaN;
                        vx2d(xx,:,zz) = NaN; vy2d(xx,:,zz) = NaN; vz2d(xx,:,zz) = NaN;
                    end
                end
                for yy=2:1:nny-1
                    if count==6
                        count=0;
                    else
                        count=count+1;
                        x2d2(:,yy,zz) = NaN; y2d2(:,yy,zz) = NaN; z2d2(:,yy,zz) = NaN;
                        vx2d(:,yy,zz) = NaN; vy2d(:,yy,zz) = NaN; vz2d(:,yy,zz) = NaN;
                    end
                end
            else  %remove values
                x2d2(:,:,zz) = NaN; y2d2(:,:,zz) = NaN; z2d2(:,:,zz) = NaN;
                vx2d(:,:,zz) = NaN; vy2d(:,:,zz) = NaN; vz2d(:,:,zz) = NaN;
            end
        end
    end
    if strcmp(GRID.Type,'spherical2D')  %for plotting the correct variables
        x2d2 = x2d_Sreduced; %[plotting dimension]
        z2d2 = z2d_Sreduced;
        if SWITCH.closeAnnulus
            vx2d(end+1,:)    = vx2d(1,:);
            vz2d(end+1,:)    = vz2d(1,:);
        end
    end
    if ~SWITCH.quiverScale %plot unit vectors instead
        vex2d   = vx2d./sqrt(vx2d.^2+vy2d.^2+vz2d.^2);
        vey2d   = vy2d./sqrt(vx2d.^2+vy2d.^2+vz2d.^2);
        vez2d   = vz2d./sqrt(vx2d.^2+vy2d.^2+vz2d.^2);
        vx2d    = vex2d;
        vy2d    = vey2d;
        vz2d    = vez2d;
    end
    numLoopsZoom = 1;
    if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
        numLoopsZoom = 2;
    end
    for iZoom=1:numLoopsZoom
        if quiver_T || quiver_eta || quiver_str || quiver_psi || quiver_rho || quiver_nstr || quiver_edot || ...
                    quiver_v || quiver_vx || quiver_vr
%             figure(1)
            if iZoom==2; axes(haxZ); end
            xLimits = get(gca,'XLim');  %Get the range of the x axis
            yLimits = get(gca,'YLim');  %Get the range of the y axis
            zLimits = get(gca,'ZLim');  %Get the range of the z axis
            if strcmp(GRID.Dim,'3-D')
                hold on
                quiver3(y2d2,x2d2,z2d2,vy2d,vx2d,vz2d,scale,'Color','k')
                axis([xLimits,yLimits,zLimits])  %make axes extent same as before quiver
            elseif strcmp(GRID.Dim,'2-D')
                hold on
                quiver(x2d2,z2d2,vx2d,vz2d,scale,...
                    'Color',PLOT.quiverColor,'LineWidth',PLOT.quiverWidth,'ShowArrowHead',PLOT.quiverArrowHead)
                if strcmp(PLOT.quiverArrowHead,'off')
                    scatter(x2d2(:),z2d2(:),scale/2,PLOT.quiverColor)
                end
                %save velocity field =====================================
                if logical(0)
                    %                     pause
                    T2d     = var2d;
                    time_Ga = PLOT.time_dim/1e9;
                    time_Ma = PLOT.time_dim/1e6;
                    x2d     = GRID.x2dp;
                    z2d     = GRID.z2dp;
                    save_string = ['/work/collapsing_plume/stagyy140211/+im/test/+output/plume_data',num2str(PLOT.number),'.mat'];
                    save(save_string,'time_Ga','x2d','z2d','vx2d','vz2d','T2d')
                end
                %save velocity field =====================================
                axis([xLimits,yLimits])  %make axes extent same as before quiver
            end
            hold off
            if strcmp(GRID.Type,'spherical2D') && plotCylCart  %plot also for spherical2D-Cartesian plot
                figure(2)
                %subplot(PLOT.nrPlots,1,SP.number)
                hold on
                quiver(x2d_C,z2d_C,VX_S,VZ_S,scale,'Color','k','ShowArrowHead',PLOT.quiverArrowHead)
                if strcmp(PLOT.quiverArrowHead,'off')
                    scatter(x2d_C(:),z2d_C(:),scale/2,PLOT.quiverColor)
                end
                hold off
                figure(1) %CHECK IF THAT IS ENOUGH TO GO BACK OR IF AXIS HAS TO BE SPECIFIED!!!!
            end
        end
        if quiver_plot && ~SAVE.quiverPlotExists %SEPARATE QUIVER PLOT
            SAVE.quiverPlotExists       = true;
            figure(4),clf; hold on;
            quiver(x2d2,z2d2,vx2d,vz2d,scale,'ShowArrowHead',PLOT.quiverArrowHead)
            if strcmp(PLOT.quiverArrowHead,'off')
                scatter(x2d2(:),z2d2(:),scale/2,PLOT.quiverColor)
            end
            if SWITCH.AxesEqual; axis equal; end
            if SWITCH.AxesOff
                axis off
            end
            if strcmp(GRID.Type,'Cartesian')
                axis tight
                %---------------------
                axis ij  % reverse axis (2.)
                %---------------------
                xlabel(PLOT.xlabel)
                ylabel(PLOT.ylabel);
            end
            if PLOT.time_dim/1e6>100
                title(['\bf{','Velocity }\rm{ ',num2str(PLOT.time_dim/1e9,3),' Gyr'])
            else
                title(['\bf{','Velocity }\rm{ ',num2str(PLOT.time_dim/1e6,3),' Myr'])
            end
            figure(1) %back to main figure
        end
        if numLoopsZoom==2; axes(AXcurrent); end
    end
end  %quiver plot



%==========================================================================
%% PRINCIPAL STRESS DIRECTION
%==========================================================================
if princStressDir_plot  || ...
        (strcmp(FIELD.name,'Temperature') && prstrDir_T) || ...
        (strcmp(FIELD.name,'Viscosity') && prstrDir_eta) || ...
        (strcmp(FIELD.name,'Stress') && prstrDir_str) || ...
        (strcmp(FIELD.name,'Streamfunction') && prstrDir_psi) || ...
        (strcmp(FIELD.name,'zz-Stress component') && prstrDir_nstr) || ...
        (strcmp(FIELD.name,'Strain rate') && prstrDir_edot) || ...
        (strcmp(FIELD.name,'Density') && prstrDir_rho) && ...
        ~strcmp(GRID.Type,'yinyang')% plot principal stress directions
    
%     if plate_velocity
%         quiver_numDel	= 6;
%         scale           = 20;
%     end
%     if isfield(PLOT,'VX_3D') && isfield(PLOT,'VY_3D') && isfield(PLOT,'VZ_3D')
%         if strcmp(GRID.Type,'spherical2D')
%             vx2d(:,:) = PLOT.VX_3Ds(:,1,:); vz2d(:,:) = PLOT.VZ_3Ds(:,1,:); vy2d = zeros(nx,nz); %JUST TO BE DEFINED
%             VX_S(:,:) = PLOT.VX_3D(:,1,:); VZ_S(:,:) = PLOT.VZ_3D(:,1,:); %this is for spherical2D-Cartesian plot
%         elseif  strcmp(GRID.Type,'Cartesian')
%             VX_3D = PLOT.VX_3D; VY_3D = PLOT.VY_3D; VZ_3D = -PLOT.VZ_3D; %adjustment for flip of depth vector
%         elseif  strcmp(GRID.Type,'yinyang')
%             %...
%         end
%     else
        % Read pressure & velocity information
        [~,~,~,~,PSDX_3D,PSDY_3D,PSDZ_3D,~,~] = f_readStagYY(FILE.directory,FILE.name,FILE.number,'Horizontal princ. stress',SWITCH);
        if ischar(PSDX_3D); error(PSDX_3D); end
        if size(PSDX_3D,1)==1 %exchange x and y
            dummy_3D = zeros(size(Y_3D,2),size(Y_3D,1),size(Y_3D,3));
            dummy_3Dx = zeros(size(Y_3D,2),size(Y_3D,1),size(Y_3D,3));
            dummy_3D(:,1,:)     = PSDZ_3D(1,:,:); PSDZ_3D = dummy_3D;
            dummy_3Dx(:,1,:)    = PSDX_3D(1,:,:);
            dummy_3D(:,1,:)     = PSDY_3D(1,:,:); PSDY_3D = dummy_3Dx; PSDX_3D = dummy_3D;
        end
        if strcmp(GRID.Type,'spherical2D')
            psdx2d        = zeros(nx,nz);
            psdy2d        = zeros(nx,nz); %JUST TO BE DEFINED
            psdz2d        = zeros(nx,nz);
            psdx2d(:,:)   = PSDX_3D(1,:,:);
            psdz2d(:,:)   = PSDZ_3D(1,:,:);
            % conversion to spherical2D grid
            PSDX_S = psdx2d; PSDZ_S = psdz2d; %this is for spherical2D-Cartesian plot
            psdx2d        = -PSDX_S.*cos(x2d_nd)+PSDZ_S.*sin(x2d_nd);
            psdz2d        = PSDX_S.*sin(x2d_nd)-PSDZ_S.*cos(x2d_nd);
        elseif  strcmp(GRID.Type,'Cartesian')
            PSDZ_3D = -PSDZ_3D; %adjustment for flip of depth vector
        end
%     end
    if strcmp(GRID.Type,'Cartesian')
        X_3D = GRID.X_3Dp;
        Y_3D = GRID.Y_3Dp;
        Z_3D = GRID.Z_3Dp;
        if SWITCH.ReduceSize
            wb = waitbar(0,'Slimming Matrix (quiver)...');
            %reduce size of matrix
            maxGridpoints = 500000;  %500000;   %<--------------------------=================
            if strcmp(GRID.Dim,'3-D')
                maxGridpoints=maxGridpoints/100;
            end
            nxEff   = nx;
            nyEff   = ny;
            nzEff   = nz;
            while (nxEff*nyEff*nzEff)>maxGridpoints
                waitbar(maxGridpoints/(nxEff*nyEff*nzEff),wb)
                for i=(nxEff-2):-2:2
                    X_3D(i,:,:)     = []; Y_3D(i,:,:)       = []; Z_3D(i,:,:)       = [];
                    PSDX_3D(i,:,:)	= []; PSDY_3D(i,:,:) 	= []; PSDZ_3D(i,:,:) 	= [];
                end
                for j=(nyEff-2):-2:2
                    X_3D(:,j,:)     = []; Y_3D(:,j,:)       = []; Z_3D(:,j,:)       = [];
                    PSDX_3D(:,j,:) 	= []; PSDY_3D(:,j,:)	= []; PSDZ_3D(:,j,:)	= [];
                end
                for k=(nzEff-2):-2:2
                    X_3D(:,:,k)  	= []; Y_3D(:,:,k)       = []; Z_3D(:,:,k)       = [];
                    PSDX_3D(:,:,k) 	= []; PSDY_3D(:,:,k) 	= []; PSDZ_3D(:,:,k) 	= [];
                end
                X_3D(4,:,:)     = [];   X_3D(:,4,:)     = [];   X_3D(:,:,4)     = [];
                Y_3D(4,:,:)     = [];   Y_3D(:,4,:)     = [];   Y_3D(:,:,4)     = [];
                Z_3D(4,:,:)     = [];   Z_3D(:,4,:)     = [];   Z_3D(:,:,4)     = [];
                PSDX_3D(4,:,:)  = [];   PSDX_3D(:,4,:)  = [];   PSDX_3D(:,:,4)  = [];
                PSDY_3D(4,:,:)  = [];   PSDY_3D(:,4,:)  = [];   PSDY_3D(:,:,4)  = [];
                PSDZ_3D(4,:,:)  = [];   PSDZ_3D(:,4,:)  = [];   PSDZ_3D(:,:,4)  = [];
                nxEff   	= size(X_3D,1);
                nyEff    	= size(Y_3D,2);
                nzEff     	= size(Z_3D,3);
            end
            disp(['Size of plotted matrix decreased to ',num2str(size(X_3D))])
            clearvars nxEff nyEff nzEff maxGridpoints
            close(wb)
        end
        if strcmp(GRID.Dim,'3-D')
            x2d2     = X_3D; y2d2     = Y_3D; z2d2     = Z_3D;
            psdx2d    = PSDX_3D; psdy2d    = PSDY_3D; psdz2d    = PSDZ_3D;
        elseif strcmp(GRID.Dim,'2-D')
            x2d2(:,:) = X_3D(:,1,:); z2d2(:,:) = Z_3D(:,1,:); y2d2 = 1; % just in order to have it defined
            psdx2d(:,:) = PSDX_3D(:,1,:); psdy2d(:,:) = PSDY_3D(:,1,:); psdz2d(:,:) = PSDZ_3D(:,1,:);
        end
    end
    if strcmp(GRID.Dim,'3-D')   % reduce number of velocity arrows  (3-D)
        count=0;
        nn=size(x2d2,1)-1;
        for i=2:1:nn-1
            if count==2
                count=0.;
            else
                x2d2(i,:) = NaN; psdx2d(i,:) = NaN;
            end
            count=count+1;
        end
        count=0;
        nn=size(x2d2,2)-1;
        for j=2:1:nn-1
            if count==2
                count=0.;
            else
                y2d2(j,:) = NaN; psdy2d(j,:) = NaN;
            end
            count=count+1;
        end
        count=0;
        nn=size(x2d2,3)-1;
        for k=2:1:nn-1
            if count==2
                count=0.;
            else
                z2d2(k,:) = NaN; psdz2d(k,:) = NaN;
            end
            count=count+1;
        end
    elseif strcmp(GRID.Dim,'2-D') % reduce number of velocity arrows  (2-D)
        if strcmp(GRID.Type,'spherical2D')
            X_Sreduced = x2d_nd; 
            x2d_Sreduced = GRID.x2ds; z2d_Sreduced = GRID.z2ds;
        elseif strcmp(GRID.Type,'Cartesian')
            %just to have them defined
            X_Sreduced = zeros(size(x2d2));
            x2d_Sreduced = zeros(size(x2d2)); z2d_Sreduced = zeros(size(x2d2));
        end
        %reduce in x-direction
        count=0;
        nn=nx-1;
        if strcmp(GRID.Type,'Cartesian')
            M2reduce=ones(size(x2d2));
        elseif strcmp(GRID.Type,'spherical2D')
            M2reduce=ones(size(X_Sreduced));
        end
        for i=2:1:nn-1
            if count==princStress_numDel
                count=0.;
            else
                M2reduce(i,:) = NaN;
            end
            count=count+1;
        end
        %remove arrows in sticky-air layer
        if strcmp(SETUP.topBC,'sticky-air') && TOPO.area && SWITCH.ActualDepth
            if strcmp(GRID.Type,'Cartesian')
                idx2remove = find(z2d2(1,:)<0);
                x2d2(:,idx2remove) = NaN; y2d2(:,idx2remove) = NaN; z2d2(:,idx2remove) = NaN;
                psdx2d(:,idx2remove) = NaN; psdy2d(:,idx2remove) = NaN; psdz2d(:,idx2remove) = NaN;
            elseif strcmp(GRID.Type,'spherical2D')
                idx2remove = find(z2dp(1,:)<0);
                X_Sreduced(:,idx2remove) = NaN; x2d_Sreduced(:,idx2remove) = NaN; z2d_Sreduced(:,idx2remove) = NaN;
                PSDX_S(:,idx2remove) = NaN; PSDZ_S(:,idx2remove) = NaN;
                psdx2d(:,idx2remove) = NaN; psdy2d(:,idx2remove) = NaN; psdz2d(:,idx2remove) = NaN;
            end
        end
        %reduce in z-direction
        count=0;
        nn=nz-1;
        for k=2:1:nn-1
            if count==princStress_numDel
                count=0.;
            else
                M2reduce(:,k) = NaN;
            end
            count=count+1;
        end
        %remove values (i.e., replace them by NaNs)
        if strcmp(GRID.Type,'Cartesian')
            x2d2 = x2d2.*M2reduce; z2d2 = z2d2.*M2reduce;
            psdx2d = psdx2d.*M2reduce; psdz2d = psdz2d.*M2reduce;
        elseif strcmp(GRID.Type,'spherical2D')
            if SWITCH.closeAnnulus
                x2d_Sreduced(1,:) = NaN; z2d_Sreduced(1,:) = NaN; %that is the one that gets dublicated
                psdx2d(1,:) = NaN; psdz2d(1,:) = NaN;
            end
            X_Sreduced = X_Sreduced.*M2reduce;
            x2d_Sreduced(1:size(M2reduce,1),:) = x2d_Sreduced(1:size(M2reduce,1),:).*M2reduce; %also works if SWITCH.closeAnnulus
            z2d_Sreduced(1:size(M2reduce,1),:) = z2d_Sreduced(1:size(M2reduce,1),:).*M2reduce;
            psdx2d = psdx2d.*M2reduce; psdz2d = psdz2d.*M2reduce;
            PSDX_S = PSDX_S.*M2reduce; PSDZ_S = PSDZ_S.*M2reduce;
        end
    end
    if plate_velocity  %CAN BE OPTIMIZED.....
        nnx=size(x2d2,1);
        nny=size(x2d2,2);
        nnz=size(x2d2,3);
        for zz=1:nnz  %loop through z direction
            count=0;
            if z2d2(1,1,zz)>0.0 && z2d2(1,1,zz)<30.0  %keep values
                for xx=2:1:nnx-1
                    if count==6
                        count=0;
                    else
                        count=count+1;
                        x2d2(xx,:,zz) = NaN; y2d2(xx,:,zz) = NaN; z2d2(xx,:,zz) = NaN;
                        psdx2d(xx,:,zz) = NaN; psdy2d(xx,:,zz) = NaN; psdz2d(xx,:,zz) = NaN;
                    end
                end
                for yy=2:1:nny-1
                    if count==6
                        count=0;
                    else
                        count=count+1;
                        x2d2(:,yy,zz) = NaN; y2d2(:,yy,zz) = NaN; z2d2(:,yy,zz) = NaN;
                        psdx2d(:,yy,zz) = NaN; psdy2d(:,yy,zz) = NaN; psdz2d(:,yy,zz) = NaN;
                    end
                end
            else  %remove values
                x2d2(:,:,zz) = NaN; y2d2(:,:,zz) = NaN; z2d2(:,:,zz) = NaN;
                psdx2d(:,:,zz) = NaN; psdy2d(:,:,zz) = NaN; psdz2d(:,:,zz) = NaN;
            end
        end
    end
    if strcmp(GRID.Type,'spherical2D')  %for plotting the correct variables
        x2d2 = x2d_Sreduced; %[plotting dimension]
        z2d2 = z2d_Sreduced;
        if SWITCH.closeAnnulus
            psdx2d(end+1,:)    = psdx2d(1,:);
            psdz2d(end+1,:)    = psdz2d(1,:);
        end
    end
    numLoopsZoom = 1;
    if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
        numLoopsZoom = 2;
    end
    for iZoom=1:numLoopsZoom
        if prstrDir_T || prstrDir_eta || prstrDir_str || prstrDir_psi || prstrDir_rho || prstrDir_nstr || prstrDir_edot
%             figure(1)
            if iZoom==2; axes(haxZ); end
            xLimits = get(gca,'XLim');  %Get the range of the x axis
            yLimits = get(gca,'YLim');  %Get the range of the y axis
            zLimits = get(gca,'ZLim');  %Get the range of the z axis
            if strcmp(GRID.Dim,'3-D')
                hold on
                quiver3(y2d2,x2d2,z2d2,psdy2d,psdx2d,psdz2d,princStress_scale,'Color','k','ShowArrowHead','off')
                axis([xLimits,yLimits,zLimits])  %make axes extent same as before quiver
            elseif strcmp(GRID.Dim,'2-D')
                %normalise data (i.e., remove differential length of arrows)
                psdx2d_plot = psdx2d./sqrt(psdx2d.^2+psdz2d.^2);
                psdz2d_plot = psdz2d./sqrt(psdx2d.^2+psdz2d.^2);
                %convert to match same directions but different sign
                psdx2d_plot(psdz2d_plot==abs(psdz2d_plot)) = -psdx2d_plot(psdz2d_plot==abs(psdz2d_plot));
                psdz2d_plot(psdz2d_plot==abs(psdz2d_plot)) = -psdz2d_plot(psdz2d_plot==abs(psdz2d_plot));
                hold on
%                 quiver(x2d2,z2d2,psdx2d,psdz2d,princStress_scale,...
%                     'Color',PLOT.princStressColor,'LineWidth',PLOT.princStressWidth,'ShowArrowHead','off')
                quiver(x2d2,z2d2,psdx2d_plot,psdz2d_plot,princStress_scale,...
                    'Color',PLOT.princStressColor,'LineWidth',PLOT.princStressWidth,'ShowArrowHead','off')
                %save principal stress field =====================================
                if logical(0)
                    %                     pause
                    T2d     = var2d;
                    time_Ga = PLOT.time_dim/1e9;
                    time_Ma = PLOT.time_dim/1e6;
                    x2d     = GRID.x2dp;
                    z2d     = GRID.z2dp;
                    save_string = ['/work/collapsing_plume/stagyy140211/+im/test/+output/plume_data',num2str(PLOT.number),'.mat'];
                    save(save_string,'time_Ga','x2d','z2d','psdx2d','psdz2d','T2d')
                end
                %save principal stress field =====================================
                axis([xLimits,yLimits])  %make axes extent same as before quiver
            end
            hold off
            if strcmp(GRID.Type,'spherical2D') && plotCylCart  %plot also for spherical2D-Cartesian plot
                figure(2)
                %subplot(PLOT.nrPlots,1,SP.number)
                hold on
                quiver(x2d_C,z2d_C,PSDX_S,PSDZ_S,princStress_scale,'Color','k')
                hold off
                figure(1) %CHECK IF THAT IS ENOUGH TO GO BACK OR IF AXIS HAS TO BE SPECIFIED!!!!
            end
        end
        if princStressDir_plot && ~SAVE.princStressPlotExists %SEPARATE QUIVER PLOT
            SAVE.princStressPlotExists  = true;
            figure(4),clf; hold on;
            %normalise data (i.e., remove differential length of arrows)
            psdx2d_plot = psdx2d./sqrt(psdx2d.^2+psdz2d.^2);
            psdz2d_plot = psdz2d./sqrt(psdx2d.^2+psdz2d.^2);
            %convert to match same directions but different sign
            psdx2d_plot(psdz2d_plot==abs(psdz2d_plot)) = -psdx2d_plot(psdz2d_plot==abs(psdz2d_plot));
            psdz2d_plot(psdz2d_plot==abs(psdz2d_plot)) = -psdz2d_plot(psdz2d_plot==abs(psdz2d_plot));
            %quiver(x2d2,z2d2,psdx2d_plot,psdz2d_plot,scale,'ShowArrowHead','off')
            quiver(x2d2,z2d2,psdx2d_plot,psdz2d_plot,'ShowArrowHead','off')
            if SWITCH.AxesEqual; axis equal; end
            if SWITCH.AxesOff
                axis off
            end
            if strcmp(GRID.Type,'Cartesian')
                axis tight
                %---------------------
                axis ij  % reverse axis (2.)
                %---------------------
                xlabel(PLOT.xlabel);
                ylabel(PLOT.ylabel);
            end
            if PLOT.time_dim/1e6>100
                title(['\bf{','Princ. Stress Direction }\rm{ ',num2str(PLOT.time_dim/1e9,3),' Gyr'])
            else
                title(['\bf{','Princ. Stress Direction }\rm{ ',num2str(PLOT.time_dim/1e6,3),' Myr'])
            end
            figure(1) %back to main figure
        end
        if numLoopsZoom==2; axes(AXcurrent); end
    end
end  %principal stress direction plot


%==========================================================================
%% LITHOSPHERE THICKNESS
%==========================================================================
%LITHO THICKNESS PLOT ADDITION
if (lithoThickness_T && strcmp(FIELD.name,'Temperature')) || ...
        (lithoThickness_eta && strcmp(FIELD.name,'Viscosity')) || ...
        (lithoThickness_nstr && strcmp(FIELD.name,'zz-Stress component')) && ...
        ~strcmp(GRID.Type,'yinyang')
    
    numLoopsZoom = 1;
    if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
        numLoopsZoom = 2;
    end
    for iZoom=1:numLoopsZoom
        if iZoom==2; axes(haxZ); end
        %read appropriate field
        %         x2d_plot,z2d_plot,var2d_plot
        [~,~,~,Z_3D,VAR_3D,~,~]     = f_readStagYY(FILE.directory,FILE.name,FILE.number,lithoThickness_field,SWITCH);
        if ischar(Z_3D); error(Z_3D); end
        if size(VAR_3D,1)==1 %exchange x and y
            dummy_3D = zeros(size(VAR_3D,2),size(VAR_3D,1),size(VAR_3D,3));
            dummy_3D(:,1,:)     = VAR_3D(1,:,:); VAR_3D = dummy_3D;
            dummy_3D(:,1,:)     = Z_3D(1,:,:); Z_3D = dummy_3D;
        end
        z2dF(:,:,:)     = Z_3D(:,:,:);
        var2dF(:,:,:)   = VAR_3D(:,:,:);
        var2dF          = var2dF *FIELD_M{strcmp(FIELD_M(:,2),lithoThickness_field),5}; %dimensionalize
        %z2dF = z2dF * FIELD_M.varScale(strcmp(FIELD_M.name,'topography')); %dimensionalize
        if strcmp(GRID.Dim,'2-D')
            z2dF(:,1,:) = z2dp(:,:);
        elseif strcmp(GRID.Dim,'3-D')
            error('check here: lithosphere thickness for 3-D is not implemented yet.')
        end
        
        %derive lithosphere thickness
        LITHO.Tisovalue = lithoThickness_isoVal;
        LITHO.plot      = false;
        if SWITCH.DimensionalMode || SWITCH.DimensionalInput
            LITHO.zmax      = lithoThickness_zmax; % [km]
            if ~strcmp(GRID.Xdim,'km'); warning(['f_PlateThickness input should be in [',GRID.Xdim,']!']); end
        else
            LITHO.zmax      = lithoThickness_zmax *1e3*GRID.ndFactor; %from [km] to grid dimensions
        end
        [LITHO] = f_PlateThickness(var2dF,z2dF,LITHO);
        
        lithoIsoline = LITHO.thickness;
        xdummy  = x2dp;
        zdummy 	= lithoIsoline;
        if strcmp(GRID.Type,'spherical2D')
            rData   = GRID.rcmb_p+max(GRID.Z_3Dp(1,1,:)) - zdummy; %calculate radii
            [xdummy,zdummy] = f_cartesian2spherical(xdummy,rData,'degrees'); %'degrees', or 'plottingDimension'
            clearvars rData
        end
        
        if strcmp(GRID.Dim,'3-D')
            error('fc: Litho Thickness option not yet implemented in 3-D!')
        elseif strcmp(GRID.Dim,'2-D')
            hold on
            plot(xdummy,zdummy,'Color',TOPO.lineColor,'LineWidth',TOPO.lineWidth);
            %plot(xdummy,zdummy,'Color','g');
        end
        hold off
        if numLoopsZoom==2; axes(AXcurrent); end
    end
end


%==========================================================================
%% TOPOGRAPHY
%==========================================================================
%TOPOGRAPHY PLOT ADDITION
if TOPO.field
    if strcmp(GRID.Dim,'3-D') || strcmp(GRID.Type,'yinyang')
        if SWITCH.Verbose; warning('Topography field not implemented for 3-D geometries!'); end
    else
        numLoopsZoom = 1;
        if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
            numLoopsZoom    = 2;
        end
        for iZoom=1:numLoopsZoom
            if iZoom==2; axes(haxZ); end
            
            if strcmp(SETUP.topBC,'free-slip')
                if SWITCH.Verbose && PLOT.loopField==1; disp(['   ',STYLE.SCHAR.indicationRightArrow,' Topography plot addition is set to zero due to free-slip surface']); end
                topo2d      = 0*TOPO.topo2d(:,end); %set topo (derived from zz-stress component) to zero
            else  %sticky air
                if strcmp(GRID.Type,'Cartesian')
                    topo2d  = -TOPO.topo2dp(:,end); %smooth or not - needs to be reversed due to reversed z-array
                else
                    topo2d  = TOPO.topo2dp(:,end); %smooth or not
                end
            end
            
            %figure(1)
            if strcmp(GRID.Dim,'2-D')
                hold on
                if TOPO.area %&& ~strcmp(GRID.Type,'spherical2D')
                    %     %*****************************
                    xpoints = x2dp(:,1)';
                    upper = topo2d';
                    if strcmp(SETUP.topBC,'free-slip') && SWITCH.AxesLimit %in case of free slip surface with fixed axes
                        lower = SWITCH.AxesLimitValues(1,3)*ones(size(topo2d,1),1)';
                    else
                        if strcmp(GRID.Type,'Cartesian')  %for plotting the correct variables
                            lower = min(z2dp(:))*ones(size(topo2d,1),1)';
                        elseif strcmp(GRID.Type,'spherical2D')
                            lower = min(GRID.z2ds(:))*ones(size(topo2d,1),1)';
                        end
                    end
                    %set topo-area color
                    color       = TOPO.color;
                    colorDiff   = 0.1; %0.07;
                    color2      = min([1 1 1],color+colorDiff); if strcmp(STYLE.ColorMode,'dark'); color2 = max([0 0 0],color-colorDiff); end
                    if strcmp(SETUP.topBC,'free-slip'); color = 'w'; if strcmp(STYLE.ColorMode,'dark'); color = 'k'; end; end %remove colouring
                    if strcmp(GRID.Type,'Cartesian')
                        %fill area
                        hf = fill([xpoints,fliplr(xpoints)],[upper,fliplr(lower)],color,'EdgeColor',color);
                        %add pattern
                        if PLOT.hatchAir && ~strcmp(SETUP.topBC,'free-slip')
                            hatchfill2(hf,'cross','HatchAngle',90,'HatchSpacing',5,'HatchColor',color2,'HatchLineWidth',1.75);
                        end
                    elseif strcmp(GRID.Type,'spherical2D')
                        if SWITCH.closeAnnulus
                            topo2d(end+1,:) = topo2d(1,:);
                            if iZoom==1 %only once
                                x2d_nd(end+1,:) = x2d_nd(1,:); z2d_nd(end+1,:) = z2d_nd(1,:);
                            end
                        end
                        topo2d_nd 	= (topo2d+SETUP.D_dim*GRID.m2p-SETUP.d_air*GRID.dimFactor) ./GRID.nd2dim; %non-dimensionalize topography
                        R           = topo2d_nd+GRID.rcmb_nd; %[nd]
                        x2d_ndS     = R.*-sin(x2d_nd(:,1)); %.*cos(y2d_nd(:,1));
                        z2d_ndS  	= R.*-cos(x2d_nd(:,1));
                        lower       = z2d_ndS*GRID.nd2dim; %dimensionalize
                        xpoints1    = x2d_ndS*GRID.nd2dim; %dimensionalize
                        
                        upper       = GRID.z2ds(:,end); %dimensionalize
                        xpoints2    = GRID.x2ds(:,end); %dimensionalize
                        
                        %fill area
                        hf = fill([xpoints2;flipud(xpoints1)],[upper;flipud(lower)],color,'EdgeColor',color);
                        %add pattern
                        if PLOT.hatchAir && ~strcmp(SETUP.topBC,'free-slip')
                            hatchfill2(hf,'cross','HatchAngle',90,'HatchSpacing',5,'HatchColor',color2,'HatchLineWidth',1.75);
                        end
                    end
                end
                if TOPO.line
                    if strcmp(GRID.Type,'Cartesian')
                        x2d_dummy   = x2dp;
                    elseif strcmp(GRID.Type,'spherical2D')
                        % Transform coordinates for spherical2D grid
                        % This needs to be in non-dimensional form (i.e., 0-1)
                        topo2d_nd 	= (topo2d+SETUP.D_dim*GRID.m2p-SETUP.d_air*GRID.dimFactor) ./GRID.nd2dim; %non-dimensionalize topography
                        R           = topo2d_nd+GRID.rcmb_nd; %[nd]
                        x2d_ndS     = R.*-sin(x2d_nd(:,1)); %.*cos(y2d_nd(:,1));
                        z2d_ndS  	= R.*-cos(x2d_nd(:,1));
                        topo2d      = z2d_ndS*GRID.nd2dim; %dimensionalize
                        x2d_dummy   = x2d_ndS*GRID.nd2dim; %dimensionalize
                    end
                    plot(x2d_dummy,topo2d,'Color',TOPO.lineColor,'LineWidth',TOPO.lineWidth);
                end
                
            end
            hold off
            if numLoopsZoom==2; axes(AXcurrent); end
        end
    end
end

%==========================================================================
%% INDICATE PLATE DIAGNOSTICS
%==========================================================================
if SWITCH.PlateDiagnostics && strcmp(GRID.Dim,'2-D')
    indColor1   = [0.8 0.0 0.0]; %red
    indColor2   = [0.8 0.8 0.8]; %light grey
    
    %switches
    indicateSlabViscCheckArea = false; indicateShallowSlabDip = false; indicateSlabTipDip = false; indicateSlabTip = false;
    indicateBending = false; indicatePlateFit = false; indicateUPtilt = false; indicatePlateCoreDepth = false;
    if (strcmp(FIELD.name,'Viscosity') || strcmp(FIELD.name,'Temperature'))
        indicateSlabViscCheckArea   = logical(0);
        if PLOT.indicateShallowSlabDip;     indicateShallowSlabDip = true;	end
        if PLOT.indicateSlabTipDip;         indicateSlabTipDip = true;      end
        if PLOT.indicateSlabTip;            indicateSlabTip = true;        	end
        if PLOT.indicateBending;            indicateBending = true;        	end
        if PLOT.indicatePlateFit;           indicatePlateFit = true;     	end
        if PLOT.indicateUPtilt;             indicateUPtilt = true;       	end
        if PLOT.indicatePlateCore;          indicatePlateCoreDepth = true;	end
    end
    %plot additions
    if indicateShallowSlabDip && exist('PLATE','var') && isfield(PLATE,'slabDipPointsX')
        if strcmp(GRID.Type,'Cartesian')
            xdummy = PLATE.slabDipPointsX(1,:)*GRID.m2p; zdummy = PLATE.slabDipPointsZ(1,:)*GRID.m2p;
        elseif strcmp(GRID.Type,'spherical2D')
            xdummy  = PLATE.slabDipPointsX(1,:); 
            zdummy 	= PLATE.slabDipPointsZ(1,:)*GRID.m2p;
            rData   = GRID.rcmb_p+max(GRID.Z_3Dp(1,1,:)) - zdummy; %calculate radii
            [xdummy,zdummy] = f_cartesian2spherical(xdummy,rData,'degrees'); %'degrees', or 'plottingDimension'
            clearvars rData
        end
        hold on
        plot(xdummy,zdummy,'-','Color',indColor1,'LineWidth',0.5)
        plot(xdummy,zdummy,'--','Color',indColor2,'LineWidth',0.5)
    end
    if indicateSlabTipDip && exist('PLATE','var') && isfield(PLATE,'SlabTipAnglePoint1')
        if strcmp(GRID.Type,'Cartesian')
            xdummy = [PLATE.SlabTipAnglePoint1(:,1),PLATE.SlabTipAnglePoint2(:,1)]; zdummy = [PLATE.SlabTipAnglePoint1(:,2),PLATE.SlabTipAnglePoint2(:,2)];
        elseif strcmp(GRID.Type,'spherical2D')
            xdummy  = [PLATE.SlabTipAnglePoint1(:,1),PLATE.SlabTipAnglePoint2(:,1)]; 
            zdummy 	= [PLATE.SlabTipAnglePoint1(:,2),PLATE.SlabTipAnglePoint2(:,2)];
            rData   = GRID.rcmb_p+max(GRID.Z_3Dp(1,1,:)) - zdummy; %calculate radii
            [xdummy,zdummy] = f_cartesian2spherical(xdummy,rData,'degrees'); %'degrees', or 'plottingDimension'
            clearvars rData
        end
        hold on
        plot(xdummy,zdummy,'-','Color',indColor1,'LineWidth',0.5)
        plot(xdummy,zdummy,'--','Color',indColor2,'LineWidth',0.5)
    end
    if indicateSlabTip && exist('PLATE','var') && ~isnan(PLATE.SlabTipPosition(1,1)) && ~isnan(PLATE.SlabTipPosition(1,2))
        if strcmp(GRID.Type,'Cartesian')
            xdummy = PLATE.SlabTipPosition(:,1); zdummy = PLATE.SlabTipPosition(:,2);
        elseif strcmp(GRID.Type,'spherical2D')
            xdummy = PLATE.SlabTipPosition(:,1); zdummy = PLATE.SlabTipPosition(:,2);
            rData   = GRID.rcmb_p+max(GRID.Z_3Dp(1,1,:)) - zdummy; %calculate radii
            [xdummy,zdummy] = f_cartesian2spherical(xdummy,rData,'degree'); %'degrees', or 'plottingDimension'
            clearvars rData
        end
        hold on
        plot(xdummy,zdummy,'x','Color',indColor1,'LineWidth',1.0)
    end
    if indicateUPtilt && exist('PLATE','var') && isfield(PLATE,'UPtiltAngleXnear') && ~isnan(PLATE.UPtiltAngleXnear(1,1))
        xdummy = [PLATE.UPtiltAngleXnear,PLATE.UPtiltAngleXfar];
        zdummy = [PLATE.UPtiltAngleZnear,PLATE.UPtiltAngleZfar];
        if strcmp(GRID.Type,'spherical2D')
            rData   = GRID.rcmb_p+max(GRID.Z_3Dp(1,1,:)) - zdummy; %calculate radii
            [xdummy,zdummy] = f_cartesian2spherical(xdummy,rData,'plottingDimension'); %'degrees', or 'plottingDimension'
            clearvars rData
        end
        hold on
        plot(xdummy,zdummy,'-','Color',indColor1,'LineWidth',0.5)
        plot(xdummy,zdummy,'--','Color',indColor2,'LineWidth',0.5)
    end    
    if indicatePlateFit && exist('PLATE','var') && isfield(PLATE,'subPlateFit') && ~isnan(PLATE.subPlateFit(1,1))
        hold on
        plot(PLATE.subPlateFit(:,1),PLATE.subPlateFit(:,2),'-','Color',indColor1,'LineWidth',0.5)
        hold on
        plot(PLATE.subPlateFit(:,1),PLATE.subPlateFit(:,2),'--','Color',indColor2,'LineWidth',0.5)
    end
    if indicateBending && exist('PLATE','var') && isfield(PLATE,'bendingCircleCenter') && isfield(PLATE,'BendingRadius') && ~isnan(PLATE.bendingCircleCenter(1,1))
        hold on
        f_plotCircle(PLATE.bendingCircleCenter(1,1),PLATE.bendingCircleCenter(1,2),PLATE.BendingRadius,indColor1,0.5,'-')
        hold on
        f_plotCircle(PLATE.bendingCircleCenter(1,1),PLATE.bendingCircleCenter(1,2),PLATE.BendingRadius,indColor2,0.5,'--')
    end
    if indicateSlabViscCheckArea && exist('PLATE','var') && isfield(PLATE,'slabArea2check') && ~isnan(PLATE.slabArea2check(1,1))
        if strcmp(GRID.Type,'Cartesian')
            xdummy  = GRID.x2dp;
            zdummy  = GRID.z2dp;
        elseif strcmp(GRID.Type,'spherical2D')
            xdummy  = GRID.x2ds(1:end-1,:);
            zdummy  = GRID.z2ds(1:end-1,:);
        end
        hold on         %plot addition to diagnose checking area
        contour(xdummy,zdummy,PLATE.slabArea2check,1,'LineColor','y')
    end
    if indicatePlateCoreDepth && exist('PLATE','var') && isfield(PLATE,'coreLevelIdx')
        clearvars xdummy zdummy
        if strcmp(GRID.Type,'Cartesian')
            for ix=1:size(PLATE.coreLevelIdx,1)
                xdummy(ix,1)  = GRID.x2dp(ix,PLATE.coreLevelIdx(ix,1));
                zdummy(ix,1)  = GRID.z2dp(ix,PLATE.coreLevelIdx(ix,1));
            end
        elseif strcmp(GRID.Type,'spherical2D')
            for ix=1:size(PLATE.coreLevelIdx,1)
                xdummy(ix,1)  = GRID.x2ds(ix,PLATE.coreLevelIdx(ix,1));
                zdummy(ix,1)  = GRID.z2ds(ix,PLATE.coreLevelIdx(ix,1));
            end
        end
        hold on
        plot(xdummy,zdummy,'-','Color',indColor1,'LineWidth',0.5)
        hold on
        plot(xdummy,zdummy,'--','Color',indColor2,'LineWidth',0.5)
    end
    clearvars xdummy zdummy
end


%==========================================================================
%% INDICATE MANTLE DIAGNOSTICS
%==========================================================================
if SWITCH.MantleDiagnostics && strcmp(GRID.Dim,'2-D')
    indColor1   = [0.8 0.0 0.0]; %red
    indColor2   = [0.8 0.8 0.8]; %light grey
    
    %switches
    indicateCONTposition = false;
    indicateLLSVPposition = false;
    if (strcmp(FIELD.name,'Viscosity') || strcmp(FIELD.name,'Temperature') || strcmp(FIELD.name,'Primordial') || strcmp(FIELD.name,'Cont. crust'))
        if PLOT.indicateCONTINENTposition;	indicateCONTposition = true;  	end
        if PLOT.indicateLLSVPposition;      indicateLLSVPposition = true;  	end
    end
    %plot additions
    if indicateCONTposition && exist('MANTLE','var') && isfield(MANTLE,'contLocationX') && ~isnan(MANTLE.contLocationX(1,1))
        if strcmp(GRID.Type,'Cartesian')
            xdummy = MANTLE.contLocationX(:,1); zdummy = MANTLE.contLocationZ(:,1);
        elseif strcmp(GRID.Type,'spherical2D')
            xdummy = MANTLE.contLocationX(:,1); zdummy = MANTLE.contLocationZ(:,1);
            rData   = GRID.rcmb_p+max(GRID.Z_3Dp(1,1,:)) - zdummy; %calculate radii
            [xdummy,zdummy] = f_cartesian2spherical(xdummy,rData,'degrees'); %'degrees', or 'plottingDimension'
            clearvars rData
        end
        hold on
        plot(xdummy,zdummy,'x','Color',indColor1,'LineWidth',1.0)
    end
    if indicateLLSVPposition && exist('MANTLE','var') && isfield(MANTLE,'llsvpLocationX') && ~isnan(MANTLE.llsvpLocationX(1,1))
        if strcmp(GRID.Type,'Cartesian')
            xdummy = MANTLE.llsvpLocationX(:,1); zdummy = MANTLE.llsvpLocationZ(:,1);
        elseif strcmp(GRID.Type,'spherical2D')
            xdummy = MANTLE.llsvpLocationX(:,1); zdummy = MANTLE.llsvpLocationZ(:,1);
            rData   = GRID.rcmb_p+max(GRID.Z_3Dp(1,1,:)) - zdummy; %calculate radii
            [xdummy,zdummy] = f_cartesian2spherical(xdummy,rData,'degrees'); %'degrees', or 'plottingDimension'
            clearvars rData
        end
        hold on
        plot(xdummy,zdummy,'x','Color',indColor1,'LineWidth',1.0)
    end
    
    clearvars xdummy zdummy
end


%==========================================================================
%% INDICATE PLATE BOUNDARIES
%==========================================================================
if PLOT.indicateTrench || PLOT.indicateRidge
    if strcmp(GRID.Type,'yinyang')
        warning('Plate boundary tracking not yet implemented for yinyang mode!');
    else
        
        % PLOT.SubPolarity can be larger than just one value!!!!.......
        
        
        
        if ~isfield(PLOT,'IndicationDepth');    PLOT.IndicationDepth = 0; end %[nd]  -0.0263
        if ~isfield(PLOT,'IndicationSize');     PLOT.IndicationSize  = 20;  end
        PBmarkerFaceColour = [0 0 0];
        PBmarkerEdgeColour = STYLE.keyColor;
        if strcmp(STYLE.ColorMode,'dark')
            PBmarkerFaceColour = 1-PBmarkerFaceColour; %invert colour
        end
        %derive plate boundaries
        if PLOT.loopField==1 && ~exist('PLATE','var')
            [PLATE,PLOT] = f_DiagnosticsPlate(FILE,GRID,SETUP,SWITCH,PLOT,STYLE); %only called if not previously called already
            
            PLOT.IndicationTrench = PLATE.Subduction;
            if isfield(PLATE,'SubPolarity'); PLOT.SubPolarity = PLATE.SubPolarity; end
            if isfield(PLATE,'Spreading'); PLOT.IndicationRidge = PLATE.Spreading; end
        else
            PLATE.Subduction = PLOT.IndicationTrench;
            if isfield(PLOT,'SubPolarity'); PLATE.SubPolarity = PLOT.SubPolarity; end
            if isfield(PLOT,'IndicationRidge'); PLATE.Spreading = PLOT.IndicationRidge; end
        end
        if ~strcmp(PLATE.Subduction,'noTracking')
            numLoopsZoom = 1;
            if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
                numLoopsZoom = 2;
            end
            for iZoom=1:numLoopsZoom
                if iZoom==2; axes(haxZ); end
                LPcolor             = [0.0 0.0 0.0]; %lower/subducting plate
                UPcolor             = [0.3 0.3 0.3]; %upper plate
                indicationDepth     = PLOT.IndicationDepth *GRID.nd2dim; %[nd]->[plotting dimension]
                hShift4indicationND = 0.06; %set horizontal shift for string [nd] or [rad]
                hShift4indication   = hShift4indicationND *GRID.nd2dim; %[rad]->[plotting dimension]
                if strcmp(GRID.Type,'spherical2D')
                    AxLim       = axis(gca); %[x1 x2 y1 y2 (z1) (z2)]
                    % Transform coordinates for spherical2D grid
                    R = max(GRID.Rp(:))+indicationDepth; %[plotting dimension]
                    if ~SWITCH.DimensionalMode && SWITCH.DimensionalInput
                        warning('plate boundary tracking for non-converted dimensional input not implemented; run in Dimensional Mode!')
                    end
                    XSrad_Subduction    = PLATE.Subduction(:,1);
                    XSrad_Spreading     = PLATE.Spreading(:,1);
                    
                    textLocXLeft0     	= XSrad_Subduction(:,1)-hShift4indicationND; %[nd]
                    textLocXRight0      = XSrad_Subduction(:,1)+hShift4indicationND; %[nd]
                    
                    [XS_Subduction,ZS_Subduction] = f_cartesian2spherical(XSrad_Subduction,R,'degree');
                    [XS_Spreading,ZS_Spreading] = f_cartesian2spherical(XSrad_Spreading,R,'degree');
                    textLocXLeft        = R.*-sin(textLocXLeft0); textLocZLeft = R.*-cos(textLocXLeft0);
                    textLocXRight       = R.*-sin(textLocXRight0); textLocZRight = R.*-cos(textLocXRight0);
                    
                    sinPosSubduction  	= sin(XSrad_Subduction(:,1));
                    cosPosSubduction   	= cos(XSrad_Subduction(:,1));
                    %                 angleHorizontal_sub     = rad2deg( -tan(XSnd_Subduction(:,1)) );
                    sinPosSpreading   	= sin(XSrad_Spreading(:,1));
                    cosPosSpreading    	= cos(XSrad_Spreading(:,1));
                    %                 angleHorizontal_spr     = rad2deg( -tan(XSnd_Spreading(:,1)) );
                    angleHorizontalL  	= rad2deg( -tan(textLocXLeft0) );
                    angleHorizontalR 	= rad2deg( -tan(textLocXRight0) );
                    
                    hold on
                    if PLOT.indicateTrench %TRENCH INDICATION
                        scatterPlot = logical(0);
                        if scatterPlot
                            scatter(XS_Subduction(:,1),ZS_Subduction(:,1),...
                                PLOT.IndicationSize,'v','fill','MarkerEdgeColor',PBmarkerEdgeColour,...
                                'MarkerFaceColor',PBmarkerFaceColour,'LineWidth',0.2)
                        else
                            axis manual %freezes scaling at current limits (needed to not alter axis limits when patch indicators are slightly outside domain
                            BIndicationAlpha    = 0.5;
                            % BIndicationHeight   = abs(max(GRID.Z_3Dp(:))-min(GRID.Z_3Dp(:))) /10; is defined above
                            BIndicationWidth    = sqrt(4/3*BIndicationHeight^2); %gleichschenkliges Dreieck
                            for ii=1:size(XS_Subduction,1)
                                if ~isnan(XS_Subduction(ii,1))
                                    arPosX = [XS_Subduction(ii,1)+cosPosSubduction(ii,1)*BIndicationWidth/2-sinPosSubduction(ii,1)*BIndicationHeight,...
                                        XS_Subduction(ii,1)-cosPosSubduction(ii,1)*BIndicationWidth/2-sinPosSubduction(ii,1)*BIndicationHeight,...
                                        XS_Subduction(ii,1)]; %left right bottom
                                    arPosZ = [ZS_Subduction(ii,1)-sinPosSubduction(ii,1)*BIndicationWidth/2-cosPosSubduction(ii,1)*BIndicationHeight,...
                                        ZS_Subduction(ii,1)+sinPosSubduction(ii,1)*BIndicationWidth/2-cosPosSubduction(ii,1)*BIndicationHeight,...
                                        ZS_Subduction(ii,1)]; %left right bottom
                                    if ~isnan(XS_Subduction(ii,1)) && XS_Subduction(ii,1)>AxLim(1,1) && XS_Subduction(ii,1)<AxLim(1,2)
                                        patch([arPosX(1,1) arPosX(1,2) arPosX(1,3)], ...
                                            [arPosZ(1,1) arPosZ(1,2) arPosZ(1,3)],...
                                            PBmarkerFaceColour,'FaceAlpha',BIndicationAlpha,...
                                            'EdgeColor',PBmarkerEdgeColour,'LineWidth',0.1);
                                    end
                                end
                            end
                        end
                    end
                    if PLOT.indicateRidge %RIDGE INDICATION
                        scatterPlot = logical(0);
                        if scatterPlot
                            scatter(XS_Spreading(:,1),ZS_Spreading(:,1),...
                                PLOT.IndicationSize,'d','fill','MarkerEdgeColor',PBmarkerEdgeColour,...
                                'MarkerFaceColor',PBmarkerFaceColour,'LineWidth',0.2)
                        else
                            axis manual %freezes scaling at current limits (needed to not alter axis limits when patch indicators are slightly outside domain
                            BIndicationAlpha    = 0.5;
                            % BIndicationHeight   = abs(max(GRID.Z_3Dp(:))-min(GRID.Z_3Dp(:))) /10; is defined above
                            BIndicationWidth    = sqrt(4/3*BIndicationHeight^2); %gleichschenkliges Dreieck
                            for ii=1:size(XS_Spreading,1)
                                if ~isnan(XS_Spreading(ii,1))
                                    arPosX = [XS_Spreading(ii,1)+cosPosSpreading(ii,1)*BIndicationWidth/2,...
                                        XS_Spreading(ii,1)-cosPosSpreading(ii,1)*BIndicationWidth/2,...
                                        XS_Spreading(ii,1)-sinPosSpreading(ii,1)*BIndicationHeight]; %left right bottom
                                    arPosZ = [ZS_Spreading(ii,1)-sinPosSpreading(ii,1)*BIndicationWidth/2,...
                                        ZS_Spreading(ii,1)+sinPosSpreading(ii,1)*BIndicationWidth/2,...
                                        ZS_Spreading(ii,1)-cosPosSpreading(ii,1)*BIndicationHeight]; %left right bottom
                                    if ~isnan(XS_Spreading(ii,1)) && XS_Spreading(ii,1)>AxLim(1,1) && XS_Spreading(ii,1)<AxLim(1,2)
                                        patch([arPosX(1,1) arPosX(1,2) arPosX(1,3)], ...
                                            [arPosZ(1,1) arPosZ(1,2) arPosZ(1,3)],...
                                            PBmarkerFaceColour,'FaceAlpha',BIndicationAlpha,...
                                            'EdgeColor',PBmarkerEdgeColour,'LineWidth',0.1);
                                    end
                                end
                            end
                        end
                    end
                    if PLOT.indicateUpperPlate %UPPER PLATE INDICATION
                        for ii=1:size(PLATE.Subduction,1)
                            if strcmp(PLATE.SubPolarity(ii,1),'left') %left-hand side subduction
                                text(textLocXLeft(ii,1),textLocZLeft(ii,1),'UP',...
                                    'Color',UPcolor,'HorizontalAlignment','center','VerticalAlignment','middle',...
                                    'Rotation',angleHorizontalL(1,1)','FontSize',PLOT.PlateIndSize);
                            elseif strcmp(PLATE.SubPolarity(ii,1),'right') %right-hand side subduction
                                text(textLocXRight(ii,1),textLocZRight(ii,1),'UP',...
                                    'Color',UPcolor,'HorizontalAlignment','center','VerticalAlignment','middle',...
                                    'Rotation',angleHorizontalR(1,1)','FontSize',PLOT.PlateIndSize);
                            else
                                %nothing found
                            end
                        end
                    end
                    if PLOT.indicateLowerPlate %LOWER PLATE INDICATION
                        for ii=1:size(PLATE.Subduction,1)
                            if strcmp(PLATE.SubPolarity(ii,1),'left') %left-hand side subduction
                                text(textLocXRight(ii,1),textLocZRight(ii,1),'LP',...
                                    'Color',LPcolor,'HorizontalAlignment','center','VerticalAlignment','middle',...
                                    'Rotation',angleHorizontalL(1,1),'FontSize',PLOT.PlateIndSize);
                            elseif strcmp(PLATE.SubPolarity(ii,1),'right') %right-hand side subduction
                                text(textLocXLeft(ii,1),textLocZLeft(ii,1),'LP',...
                                    'Color',LPcolor,'HorizontalAlignment','center','VerticalAlignment','middle',...
                                    'Rotation',angleHorizontalR(1,1)','FontSize',PLOT.PlateIndSize);
                            else
                                %nothing found
                            end
                        end
                    end
                    
                elseif strcmp(GRID.Type,'Cartesian')
                    hold on
                    AxLim       = axis(gca); %[x1 x2 y1 y2 (z1) (z2)]
                    if strcmp(get(gca,'XDir'),'reverse') %account for reverse axes
                        dummy=AxLim(1,2); AxLim(1,2)=AxLim(1,1); AxLim(1,1)=dummy;
                    end
                    if strcmp(get(gca,'YDir'),'reverse') %account for reverse axes
                        dummy=AxLim(1,4); AxLim(1,4)=AxLim(1,3); AxLim(1,3)=dummy;
                    end
                    
                    if strcmp(GRID.Dim,'3-D')
                        warning('fc: TRACKING PLATE BOUNDARIES IN 3-D still in beta phase!')
                        
                        if strcmp(get(gca,'ZDir'),'reverse') %account for reverse axes
                            dummy=AxLim(1,6); AxLim(1,6)=AxLim(1,5); AxLim(1,5)=dummy;
                        end
                        indZlevel = 0.1*min(zlim); %indication z-level; e.g., min(zlim)
                        PLATE.BoundariesColored     = logical(0);
                        PLATE.BoundariesEdgeColor   = logical(0);
                        BIndicationAlpha            = 1.0;
                        TriangleHeight      = 0.05 *(max(GRID.Z_3Dp(:))-min(GRID.Z_3Dp(:))); %[plotting dimension]
                        TriangleWidth       = sqrt(4/3*TriangleHeight^2); %gleichschenkliges Dreieck [plotting dimension]
                        if PLATE.BoundariesColored
                            colorBNDRY = [0 0.2 0.4]; colorSPR = [0.4 0.2 0];
                        else
                            colorBNDRY = PBmarkerFaceColour; colorSPR = colorBNDRY;
                        end
                        if PLATE.BoundariesEdgeColor
                            colorBNDRYedge = PBmarkerEdgeColour;
                        else
                            colorBNDRYedge = 'none';
                        end
                        if PLOT.indicateTrench
                            %plot3(PLATE.Subduction(:,2),PLATE.Subduction(:,1),PLATE.Subduction(:,1).*indZlevel,'b')
                            for iarea=1:size(PLATE.Subduction,1)
                                indicateTrenchTriangles = logical(1);
                                if indicateTrenchTriangles
                                    for ipoint=1:size(PLATE.SubCenterP,2) %loop only through characteristic points
                                        arPosX = [PLATE.SubCenterP(iarea,ipoint,1)-PLATE.SubStrikeP(iarea,ipoint,1)*TriangleWidth/2,...
                                            PLATE.SubCenterP(iarea,ipoint,1)+PLATE.SubStrikeP(iarea,ipoint,1)*TriangleWidth/2,...
                                            PLATE.SubCenterP(iarea,ipoint,1)+PLATE.LowerPlateP(iarea,ipoint,1)*TriangleHeight]; %left right middle
                                        arPosY = [PLATE.SubCenterP(iarea,ipoint,2)-PLATE.SubStrikeP(iarea,ipoint,2)*TriangleWidth/2,...
                                            PLATE.SubCenterP(iarea,ipoint,2)+PLATE.SubStrikeP(iarea,ipoint,2)*TriangleWidth/2,...
                                            PLATE.SubCenterP(iarea,ipoint,2)+PLATE.LowerPlateP(iarea,ipoint,2)*TriangleHeight]; %left right middle
                                        patch([arPosY(1,1) arPosY(1,2) arPosY(1,3)], ...
                                            [arPosX(1,1) arPosX(1,2) arPosX(1,3)],...
                                            [indZlevel indZlevel indZlevel],...
                                            PBmarkerFaceColour,'FaceAlpha',BIndicationAlpha,...
                                            'EdgeColor',colorBNDRYedge,'LineWidth',0.1); %trench triangles
                                        if PLOT.indicateLowerPlate %LOWER PLATE INDICATION
                                            text(PLATE.SubCenterP(iarea,ipoint,2)+PLATE.LowerPlateP(iarea,ipoint,2)*TriangleHeight*2,...
                                                PLATE.SubCenterP(iarea,ipoint,1)+PLATE.LowerPlateP(iarea,ipoint,1)*TriangleHeight*2,indZlevel,...
                                                'LP','BackgroundColor','w')
                                        end
                                        if PLOT.indicateUpperPlate %UPPER PLATE INDICATION
                                            text(PLATE.SubCenterP(iarea,ipoint,2)+PLATE.UpperPlateP(iarea,ipoint,2)*TriangleHeight*2,...
                                                PLATE.SubCenterP(iarea,ipoint,1)+PLATE.UpperPlateP(iarea,ipoint,1)*TriangleHeight*2,indZlevel,...
                                                'UP','BackgroundColor','w')
                                        end
                                    end
                                end
                                %plot(PLATE.Subduction(iarea,:,2),PLATE.Subduction(iarea,:,1),'b.')
                                patch(PLATE.SubductionOut(iarea,~isnan(PLATE.SubductionOut(iarea,:,1)),2),... %make sure to only patch non-NaN values!
                                    PLATE.SubductionOut(iarea,~isnan(PLATE.SubductionOut(iarea,:,1)),1),...
                                    ones(size(PLATE.SubductionOut(iarea,~isnan(PLATE.SubductionOut(iarea,:,1)),1))).*indZlevel,...
                                    colorBNDRY,'EdgeColor',colorBNDRYedge,'LineWidth',0.1,'FaceAlpha',BIndicationAlpha) %trench line
                            end
                        end
                        if PLOT.indicateRidge
                            TriangleHeight  = TriangleHeight/2;
                            TriangleWidth   = TriangleWidth/2;
                            %plot3(PLATE.Spreading(:,2),PLATE.Spreading(:,1),PLATE.Spreading(:,1).*indZlevel,'r')
                            for iarea=1:size(PLATE.Spreading,1)
                                indicateRidgeTriangles = logical(1);
                                if indicateRidgeTriangles
                                    for ipoint=1:size(PLATE.SprCenterP,2) %loop only through characteristic points
                                        for ipair=1:2 %loop to produce a pair of indicators
                                            if ipair==1; dummy=1; else; dummy=-1; end
                                            arPosX = [PLATE.SprCenterP(iarea,ipoint,1)-PLATE.SprStrikeP(iarea,ipoint,1)*TriangleWidth/2,...
                                                PLATE.SprCenterP(iarea,ipoint,1)+PLATE.SprStrikeP(iarea,ipoint,1)*TriangleWidth/2,...
                                                PLATE.SprCenterP(iarea,ipoint,1)+dummy*PLATE.SprNormalP(iarea,ipoint,1)*TriangleHeight]; %left right middle
                                            arPosY = [PLATE.SprCenterP(iarea,ipoint,2)-PLATE.SprStrikeP(iarea,ipoint,2)*TriangleWidth/2,...
                                                PLATE.SprCenterP(iarea,ipoint,2)+PLATE.SprStrikeP(iarea,ipoint,2)*TriangleWidth/2,...
                                                PLATE.SprCenterP(iarea,ipoint,2)+dummy*PLATE.SprNormalP(iarea,ipoint,2)*TriangleHeight]; %left right middle
                                            patch([arPosY(1,1) arPosY(1,2) arPosY(1,3)], ...
                                                [arPosX(1,1) arPosX(1,2) arPosX(1,3)],...
                                                [indZlevel indZlevel indZlevel],...
                                                PBmarkerFaceColour,'FaceAlpha',BIndicationAlpha,...
                                                'EdgeColor',colorBNDRYedge,'LineWidth',0.1); %ridge triangles
                                        end
                                    end
                                end
                                %plot(PLATE.Spreading(iarea,:,2),PLATE.Spreading(iarea,:,1),'r.')
                                patch(PLATE.SpreadingOut(iarea,~isnan(PLATE.SpreadingOut(iarea,:,1)),2),... %make sure to only patch non-NaN values!
                                    PLATE.SpreadingOut(iarea,~isnan(PLATE.SpreadingOut(iarea,:,1)),1),...
                                    ones(size(PLATE.SpreadingOut(iarea,~isnan(PLATE.SpreadingOut(iarea,:,1)),1))).*indZlevel,...
                                    colorSPR,'EdgeColor',colorBNDRYedge,'LineWidth',0.1,'FaceAlpha',BIndicationAlpha)
                            end
                        end
                        
                    elseif strcmp(GRID.Dim,'2-D')
                        axis manual %freezes scaling at current limits (needed to not alter axis limits when patch indicators are slightly outside domain
                        BIndicationAlpha    = 0.5;
                        BIndicationBot      = 0;
                        BIndicationTop      = AxLim(1,4);
                        if ~strcmp(SETUP.topBC,'sticky-air') || BIndicationTop==0 %e.g., for free-slip cases
                            BIndicationTop  = 1/20*abs(AxLim(1,3)-AxLim(1,4));  %<<<<< set indication height here <<<<<<<<
                            if strcmp(get(gca,'YDir'),'reverse') %account for reverse axes
                                BIndicationTop = -BIndicationTop;
                            end
                        end
                        BIndicationHeight   = abs(BIndicationTop-BIndicationBot);
                        BIndicationWidth    = sqrt(4/3*BIndicationHeight^2); %gleichschenkliges Dreieck
                        if PLOT.indicateTrench
                            for ii=1:size(PLATE.Subduction,1)
                                arPosX = [PLATE.Subduction(ii,1)-BIndicationWidth/2, PLATE.Subduction(ii,1)+BIndicationWidth/2, PLATE.Subduction(ii,1)]; %left right bottom
                                arPosZ = [BIndicationTop, BIndicationTop, BIndicationBot]; %left right bottom
                                if ~isnan(PLATE.Subduction(ii,1)) && PLATE.Subduction(ii,1)>AxLim(1,1) && PLATE.Subduction(ii,1)<AxLim(1,2)
                                    hptrench = patch([arPosX(1,1) arPosX(1,2) arPosX(1,3)], ...
                                        [arPosZ(1,1) arPosZ(1,2) arPosZ(1,3)],...
                                        PBmarkerFaceColour,'FaceAlpha',BIndicationAlpha,...
                                        'EdgeColor',PBmarkerEdgeColour,'LineWidth',0.1,'Clipping','off');
                                end
                            end
                        end
                        if PLOT.indicateRidge
                            for ii=1:size(PLATE.Spreading,1)
                                arPosX = [PLATE.Spreading(ii,1)-BIndicationWidth/2, PLATE.Spreading(ii,1)+BIndicationWidth/2, PLATE.Spreading(ii,1)]; %left right top
                                arPosZ = [BIndicationBot, BIndicationBot, BIndicationTop]; %left right top
                                if ~isnan(PLATE.Spreading(ii,1)) && PLATE.Spreading(ii,1)>AxLim(1,1) && PLATE.Spreading(ii,1)<AxLim(1,2)
                                    hpridge = patch([arPosX(1,1) arPosX(1,2) arPosX(1,3)], ...
                                        [arPosZ(1,1) arPosZ(1,2) arPosZ(1,3)],...
                                        PBmarkerFaceColour,'FaceAlpha',BIndicationAlpha,...
                                        'EdgeColor',PBmarkerEdgeColour,'LineWidth',0.1,'Clipping','off');
                                end
                            end
                        end
                        if PLOT.indicateUpperPlate %UPPER PLATE INDICATION
                            for ii=1:size(PLATE.Subduction,1)
                                if strcmp(PLATE.SubPolarity(ii,1),'left') %left
                                    text(PLATE.Subduction(ii,1)-hShift4indication,indicationDepth,'UP',...
                                        'Color',UPcolor,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',PLOT.PlateIndSize,'Clipping','off');
                                elseif strcmp(PLATE.SubPolarity(ii,1),'right') %right
                                    text(PLATE.Subduction(ii,1)+hShift4indication,indicationDepth,'UP',...
                                        'Color',UPcolor,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',PLOT.PlateIndSize,'Clipping','off');
                                else
                                    %nothing found
                                end
                            end
                        end
                        if PLOT.indicateLowerPlate %LOWER PLATE INDICATION
                            for ii=1:size(PLATE.Subduction,1)
                                if strcmp(PLATE.SubPolarity(ii,1),'left') %left
                                    text(PLATE.Subduction(ii,1)+hShift4indication,indicationDepth,'LP',...
                                        'Color',LPcolor,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',PLOT.PlateIndSize,'Clipping','off');
                                elseif strcmp(PLATE.SubPolarity(ii,1),'right') %right
                                    text(PLATE.Subduction(ii,1)-hShift4indication,indicationDepth,'LP',...
                                        'Color',LPcolor,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',PLOT.PlateIndSize,'Clipping','off');
                                else
                                    %nothing found
                                end
                            end
                        end
                    end
                end
                %             if strcmp(SETUP.topBC,'free-slip') && iZoom==1
                %             else  %sticky air
                %             end
                hold off
                if numLoopsZoom==2; axes(AXcurrent); end
            end
            
            
            %% SAVING SUBDUCTION-ZONE DATA ------------------------
            if SAVE.subductionZoneData && PLOT.loopField==1
                % trench x-location, subduction polarity, trench/UP vel, LP vel, Convergence vel
                % slab angle, slab viscosity, u.mantle viscosity, slab-mantle visc. difference,
                % left plate thickness, right plate thickness, lower plate thickness, upper plate thickness
                %polarity: '1': left, '2': right, or '0': unknown
                SAVE.Directory              = FILE.directory;
                SAVE.DataName               = [FILE.name,'_subZone',num2str(FILE.number)];
                SAVE.data                   = [PLATE.Subduction(:,1),PLATE.SubPolarityNum(:,1),...          %[plotting dimension]
                    PLATE.UPvelocity(:,1),PLATE.LPvelocity(:,1),PLATE.PlateConvergence(:,1),... %[cm/a]
                    PLATE.ShallowSlabAngle(:,1),... %[degree]
                    PLATE.SlabViscosity(:,1), PLATE.UMantleViscosity,PLATE.SlabMantleViscDiff(:,1),... %[Pas]
                    PLATE.PlateThicknessL(:,1),PLATE.PlateThicknessR(:,1),PLATE.PlateThicknessLP(:,1),PLATE.PlateThicknessUP(:,1),... %[plotting dimension]
                    PLATE.SlabViscosity(:,1),PLATE.UMantleViscosity(:,1),PLATE.SlabMantleViscDiff(:,1),... %[Pas]
                    PLATE.BendingRadius(:,1),PLATE.ViscDissBending(:,1),PLATE.ViscDissBendingRel(:,1)]; %[plotting dimension],[m2/s2],[nd]
                SAVE.dat                    = logical(0);
                SAVE.txt                    = logical(0);
                SAVE.mat                    = logical(1);
                SAVE.write2existing         = logical(0);
                [SAVE.overwriteAll] = f_saveData( SAVE );
            end
        end
    end
end


%% INDICATE SPECIFIC DEPTH LEVEL(S)
if PLOT.indicateDepth
    numLoopsZoom = 1;
    if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
        numLoopsZoom = 2;
    end
    for iZoom=1:numLoopsZoom
        if iZoom==2; axes(haxZ); end
        hold on
        for iline=1:length(PLOT.indicateDepthLevels)
            plot(get(gca,'xlim'),[PLOT.indicateDepthLevels(iline),PLOT.indicateDepthLevels(iline)],...
                'Color',STYLE.keyLineColor,'LineStyle','--');
        end
        hold off
        if numLoopsZoom==2; axes(AXcurrent); end
    end
end

%% INSERT RECTANGLE INDICATING ZOOM PLOT
if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
    numLoopsZoom = 1;
    if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
        numLoopsZoom = 2;
    end
    for iZoom=1:numLoopsZoom
        if iZoom==2; axes(haxZ); end
        
        for i_rect=1:2
            hold on
            %rectangle('Position',[x,y,w,h]) %white bold in the back
            h_rectangle = rectangle('Position',[PLOT.magnifierExtent(PLOT.loopCase,1),...
                PLOT.magnifierExtent(PLOT.loopCase,3),...
                PLOT.magnifierExtent(PLOT.loopCase,2)-PLOT.magnifierExtent(PLOT.loopCase,1),...
                PLOT.magnifierExtent(PLOT.loopCase,4)-PLOT.magnifierExtent(PLOT.loopCase,3)]);
            h_rectangle.LineWidth = 0.75;
            h_rectangle.EdgeColor = STYLE.BlackOrWhiteColor;
            %rectangle('Position',[x,y,w,h])
            h_rectangle = rectangle('Position',[PLOT.magnifierExtent(PLOT.loopCase,1),...
                PLOT.magnifierExtent(PLOT.loopCase,3),...
                PLOT.magnifierExtent(PLOT.loopCase,2)-PLOT.magnifierExtent(PLOT.loopCase,1),...
                PLOT.magnifierExtent(PLOT.loopCase,4)-PLOT.magnifierExtent(PLOT.loopCase,3)]);
            h_rectangle.LineWidth = 0.25;
            h_rectangle.EdgeColor = 1-STYLE.BlackOrWhiteColor;
        end
    
        if numLoopsZoom==2; axes(AXcurrent); end
    end
end

%% PLATE COMPOSITION
if SWITCH.trackC && strcmp(FIELD.name,'Composition')
    %find specified compostional field
    findCval = 1; %this is primordial material
    findRange = 0.5;
    Cfield  = var2d(var2d<findCval+findRange & var2d>findCval-findRange); %make small range instead of one exact findCvalue
    Cfieldx = x2dp(var2d<findCval+findRange & var2d>findCval-findRange);
    Cfieldz = z2dp(var2d<findCval+findRange & var2d>findCval-findRange);
    %find specified extreme of the compositional field
    %maximum height (at outermost left-hand side):
    maxValuesx  = Cfieldx(Cfieldz==min(Cfieldz(:))); %shallowest
    maxValx     = min(maxValuesx(:)); %most left
    maxValuesz  = Cfieldz(Cfieldz==min(Cfieldz(:))); %shallowest
    maxValz     = maxValuesz(maxValuesx(:)==maxValx);
    %create/add to data array
    if PLOT.number==SAVE.StartNumber %first time
        %[dimensional time in years  /  x  /  z  ]
        PLOT.topCval = [ PLOT.time maxValx maxValz ];
    else
        PLOT.topCval = [ SAVE.topCval;  PLOT.time maxValx maxValz ];
    end
    %plot indication point
    hold on
    plot(PLOT.topCval(end,2),PLOT.topCval(end,3),'or')
    %save data
    if FILE.LastFile
        SAVE.Directory              = FILE.directory;
        SAVE.DataName               = ['+topCval_',FILE.name];
        SAVE.data                   = PLOT.topCval;
        SAVE.dat                    = logical(1);
        SAVE.txt                    = logical(0);
        SAVE.mat                    = logical(1);
        SAVE.write2existng          = logical(0);
        
        [SAVE.overwriteAll] = f_saveData( SAVE );
    end
end


%==========================================================================
%% STYLE PLOTS
%==========================================================================
figure(1)
%SIMPLIFY FIGURE
if SWITCH.SimplifyPlots
    %GATHER INFORMATION ABOUT SUBPLOT LAYOUT
    PLOT.spcurrent          = SP.current;
    PLOT.sp_chax            = AXcurrent;
    if exist('haxZ','var'); PLOT.sp_chaxz = haxZ; else; PLOT.sp_chaxz = NaN; end
    PLOT.sp_ccb             = PLOT.cb;           %colorbar handle
    PLOT.title              = PLOT.titleStringCurrent;  	%title string
    
    [PLOT] = f_DesignLayout(SWITCH,PLOT,GRID,FIELD);
end
%ADJUST TICK MARKS
set(gca,'Layer','top') % make sure axes ticks are visible
set(gca,'TickDir',SWITCH.TickDirection) %set tickmarks out- or inside
if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
    axes(haxZ) %set zoom plot back on top
end
% PLOT.layout         = SP.current(1,1:2);
if isfield(PLOT,'hax'); nr_hax = size(PLOT.hax,1)+1; else; nr_hax = 1; end
PLOT.hax(nr_hax,1)  = AXcurrent;
if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
    if isfield(PLOT,'haxZ'); nr_haxZ = size(PLOT.haxZ,1)+1; PLOT.haxZ(nr_haxZ,1) = haxZ;
    else; nr_haxZ = 1; PLOT.haxZ(nr_haxZ,1) = haxZ; end
end


%% FUNCTION OUTPUT
SAVE.PlotHandles = [SAVE.PlotHandles; AXcurrent]; %remember graph handles (main subplot)
if isfield(SWITCH,'Magnifier') && SWITCH.Magnifier
    SAVE.PlotHandles = [SAVE.PlotHandles; haxZ]; %remember graph handles (zoom subplot)
end






%% INTERNAL FUNCTIONS


%% PLOT CIRCLE
function f_plotCircle(x,y,r,color,linewidth,linestyle)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle
    %0.01 is the angle step, bigger values will draw the circle faster but
    %you might notice imperfections (not very smooth)
    ang     = 0:0.01:2*pi;
    xp      = r*cos(ang);
    yp      = r*sin(ang);
    plot(x+xp,y+yp,'Color',color,'LineWidth',linewidth,'LineStyle',linestyle);
end


%% CARTESIAN-TO-SPHERICAL COORDINATES CONVERSION
function [xDataS,zDataS] = f_cartesian2spherical(xData,rData,dataDimension) %'degrees', or 'plottingDimension'
    % xData must be in radians or whatever length dimension
    % rData must be either the same size as xData or just one value!
    
    % ERROR CHECK
    if min(size(rData)~=size(xData)) && min(~size(rData)==1)
        error('rData must be either one value or have the same size as xData! Check here!')
    end
    
    % CONVERSION TO DEGREES/RADIANS (alpha = L/r)
    if strcmp(dataDimension,'degrees')
        %nothing to do
    elseif strcmp(dataDimension,'plottingDimension')
        xData       = xData ./rData; 	%[degrees/rad]
    end
    
    % COORDINATE CONVERSION
    xDataS   = rData.*-sin(xData); %.*cos(y2d_nd);
    zDataS   = rData.*-cos(xData);
end



end




