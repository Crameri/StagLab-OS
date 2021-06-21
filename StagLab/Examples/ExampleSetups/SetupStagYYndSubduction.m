
%%                                                 TEST PARAMETER FILE 2-D
%                                                         for SL_FieldPlot

function [TEST] = SetupStagYYndSubduction(TEST)

SAVE.StartNumber       =   0;
SAVE.StepNumber        =   1;
SAVE.EndNumber         =   0;      %set big in order to continue to last file, set SAVE.StartNum==SAVE.EndNumberto write only one file

for iLoop=SAVE.StartNumber:SAVE.StepNumber:SAVE.EndNumber
clearvars -except iLoop SAVE; if ~exist('f_Defaults','file'); f_INSTALL; end
    [IN,PLOT,SWITCH,TOPO,STYLE,SAVE] = f_Defaults(SAVE);    %get default variables

    %% INPUT FILE(S)
    IN.Name                     =   { 'StagYYndSubduction' 'StagYYndSubduction'	'StagYYndSubduction' 'StagYYndSubduction' 'StagYYndSubduction' 'StagYYndSubduction'};	%filename; it's size controls the # of subplots
    IN.Number                   =   [   360                 390                  400                  405                  410                  422    ];            	%filenumber
    IN.Parameter                =   [  	11    ];         	%parameter dimensionalisation (see f_dimParameter)
    IN.Folder                   =   {'../'};                %file directory

    %% DIMENSIONALISATION
    SWITCH.DimensionalMode    	=   logical(1);             %converts to dimensional output as set in f_dimParameter.m  
    
    %% SUBPLOT LAYOUT
    PLOT.Layout                 =   [2 3];                  %manually set subplot layout; only valid if uncommented; [z x]
    SWITCH.TimeEvolutionModeMode=   logical(1);     
        SWITCH.timeDirection    =   'leftright';            %'topbot', 'leftright'
    
    %% POSITIONING
    PLOT.shiftUpward            =   0.05;                   %shift subplot vertically closer or further away from each other
        
    %% COLORBAR
    SWITCH.ConstantColorbar   	=   logical(1);

    %% ANNOTATIONS
    PLOT.titleString            =   {'1' '2' '3' '4' '5' '6'};	%adds manual title string if uncommented

    %% POST-PROCESSING
    PLOT.indicateTrench         =   logical(1);
    PLOT.indicateRidge          =   logical(1);
    PLOT.indicateShallowSlabDip	=   logical(1);
    PLOT.indicatePlateFit       =   logical(1);             %fit of the subducting plate portion
    PLOT.indicateBending        =   logical(1);
    
    %% SAVING FIGURE
    SAVE.Figure                 =   logical(1);             %saves plot directly to directory
  	SAVE.writeDirectory         =   '../ExampleFigures/';   %'auto': save to standard folder (.../+im/...); otherwise use e.g., '/work/stagyy/'

    %% FIELDS TO PLOT
    PLOT.Temperature            =   logical(0); 
    PLOT.Viscosity              =   logical(1);
    
    %% FIELD PLOT ADDITIONS
    IN.streamfun_v              =   logical(1);
       
    %% RUN MAIN ROUTINE >>>
    [SAVE] = SL_FieldPlot(IN,PLOT,SWITCH,TOPO,STYLE,SAVE);  %run main routine
    if SAVE.LastFile; break; end %last file reached
end

TEST.successful = true;
end

