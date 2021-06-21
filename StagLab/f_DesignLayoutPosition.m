
%%                                         ADJUSTING SUBPLOT POSITIONS 2.10
% Input:
% SPOS.task = 'general'; %'specific', 'general'
% SPOS.action = 'shiftLeftSingleColumn';  %specific action
% SPOS.subplotLayout = [2 3];
% [SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE);
%
% Output:
% SP.current    %current suplot, e.g., [2,2,1]
% SP.number     %subplot number
%                                                Fabio Crameri, 21.06.2019

function [SP,SAVE] = f_DesignLayoutPosition(SPOS,PLOT,SWITCH,SAVE)

%% DEFAULTS
SP.current  = NaN;

%% CREATING SUBPLOTS
if strcmp(SPOS.task,'createSubplot')
    if isfield(PLOT,'Layout') %manual arrangement
        SP.current = [PLOT.Layout(1,1),PLOT.Layout(1,2),PLOT.nrSubplot];
    else
        if PLOT.nrPlots>=4
            if mod(PLOT.nrPlots,2)
                SP.current = [(PLOT.nrPlots+1)/2 2 PLOT.nrSubplot];
            else
                SP.current = [PLOT.nrPlots/2 2 PLOT.nrSubplot];
            end
        else
            SP.current = [PLOT.nrPlots 1 PLOT.nrSubplot];
        end
    end
    if SWITCH.ReverseLayout %top-down plotting
        if ~isfield(PLOT,'layoutNum')
            PLOT.layoutNum = zeros(SP.current(1,1),SP.current(1,2));
            idx = 0;
            for ii=1:size(PLOT.layoutNum,1)
                for jj=1:size(PLOT.layoutNum,2)
                    idx = idx+1;
                    PLOT.layoutNum(ii,jj) = idx;
                end
            end
            idx = 0; PLOT.layoutNumTransposed = PLOT.layoutNum;
            for jj=1:size(PLOT.layoutNum,2)
                for ii=1:size(PLOT.layoutNum,1)
                    idx = idx+1;
                    PLOT.layoutNumTransposed(ii,jj) = idx; %transposed layout
                end
            end
        end
        %UPDATE TO NEW LAYOUT
        %       while PLOT.nrSubplot is still setup the old way!
        SP.number       = PLOT.layoutNum(PLOT.layoutNumTransposed==PLOT.nrSubplot);
        SP.current(1,3)	= SP.number;
    else
        SP.number       = PLOT.nrSubplot;
    end
    [ixPlot,izPlot] = ind2sub([SP.current(1,2),SP.current(1,1)],SP.current(1,3));
    
    %subplot...
    if SWITCH.PlotInPlot && isfield(SPOS,'haxPlotInPlot')
        axes(SPOS.haxPlotInPlot); %axis activation; needed to prevent erasing subplot
        hold on
    else
        %drawnow
        if SWITCH.UsePanel
            hax = SAVE.P(izPlot,ixPlot).select();
            set(gcf,'CurrentAxes',hax)
        else
            subplot(SP.current(1,1),SP.current(1,2),SP.current(1,3));
        end
    end
    
    %pause...
    pause(SWITCH.spPause)


%% SPECIFIC ADJUSTMENTS
elseif strcmp(SPOS.task,'specific')
    if strcmp(SPOS.action,'shiftLeftSingleColumn')
        posAX = get(gca, 'pos');
        posAX(1,1) = 0.04;
        set(gca, 'pos', posAX);
        
    elseif strcmp(SPOS.action,'decreaseStandardGap')
        p_subplot = get(gca, 'pos');
        p_subplot(3) = p_subplot(3) + 0.05;  %add some amount to the width of the current subplot
        set(gca, 'pos', p_subplot);
        
    elseif strcmp(SPOS.action,'decreaseColorbarGap')
        p_subplot = get(gca, 'pos');
        p_subplot(3) = p_subplot(3) + 0.07;  %add some amount to the width of the current subplot
        set(gca, 'pos', p_subplot);
        
    elseif strcmp(SPOS.action,'decreaseXGap')
%         if SPOS.shiftX~=0  %KEEP UNCOMMENTED TO KEEP ALL PLOTS FIXED IN POSITION
            p_subplot = get(gca, 'pos');
            p_subplot(1) = p_subplot(1) - SPOS.shiftX;  %shift left
            set(gca, 'pos', p_subplot);
%         end
        
    elseif strcmp(SPOS.action,'decreaseYGap')
%         if SPOS.shiftY~=0  %KEEP UNCOMMENTED TO KEEP ALL PLOTS FIXED IN POSITION
            p_subplot = get(gca, 'pos');
            p_subplot(2) = p_subplot(2) + SPOS.shiftY;  %shift upwards
            set(gca, 'pos', p_subplot);
%         end
        
        %         ...
    else
        error(['SPOS.action = ',SPOS.action,' not available!'])
    end


%% GENERAL ADJUSTMENTS
elseif strcmp(SPOS.task,'general')
    if (SPOS.subplotLayout(1,1)==2 && SPOS.subplotLayout(1,2)==3) || ... %for 2x3 plot
            (SPOS.subplotLayout(1,1)==1 && SPOS.subplotLayout(1,2)==3) %for 1x3 plot
        p_subplot = get(gca, 'pos');
        p_subplot(3) = p_subplot(3) + 0.05;  %add some amount to the width of the current subplot
        set(gca, 'pos', p_subplot);
        
        p_subplot = get(gca, 'pos');
        p_subplot(1,1) = p_subplot(1,1) - 0.05;  %shift for some amount to the left
        p_subplot(1,1) = p_subplot(1,1) - 0.01*p_subplot(1,1);  %add percentage of subplot x-pos to the x-pos of the current subplot
        set(gca, 'pos', p_subplot);
    end
    if SPOS.subplotLayout(1,1)==3 && SPOS.subplotLayout(1,2)==3 %for 3x3 plot
        p_subplot = get(gca, 'pos');
        p_subplot(3) = p_subplot(3) + 0.05;  %add some amount to the width of the current subplot
        set(gca, 'pos', p_subplot);
        
        p_subplot = get(gca, 'pos');
        p_subplot(1,1) = p_subplot(1,1) - 0.05;  %shift for some amount to the left
        p_subplot(1,1) = p_subplot(1,1) - 0.02*p_subplot(1,1);  %add percentage of subplot x-pos to the x-pos of the current subplot
        set(gca, 'pos', p_subplot);
    end
    if SPOS.subplotLayout(1,1)==4 % && SPOS.subplotLayout(1,2)==3; %for 4x3 plot
        p_subplot = get(gca, 'pos');
        p_subplot(3) = p_subplot(3) + 0.05;  %add some amount to the width of the current subplot
        set(gca, 'pos', p_subplot);
        
        p_subplot = get(gca, 'pos');
        p_subplot(1,1) = p_subplot(1,1) - 0.0;  %shift for some amount to the left
        p_subplot(1,1) = p_subplot(1,1) - 0.08*p_subplot(1,1);  %add percentage of subplot x-pos to the x-pos of the current subplot
        p_subplot(1,2) = p_subplot(1,2) + 0.01*(1-p_subplot(1,2));  %add percentage of subplot y-pos to the y-pos of the current subplot
        set(gca, 'pos', p_subplot);
    end
end

%% CLEAR FIELDS
% SPOS = rmfield(SPOS,'task');
% SPOS = rmfield(SPOS,'action');

end
