
%%                                   SET MIN & MAX COLORBAR PARAMETERS 1.71
%
%                                                Fabio Crameri, 20.07.2020

function [CB] = f_DefaultsColorbar(IN,dimensional,loop_case)

%% DEFAULT VALUES
CB.multipleMinMax = false;
if dimensional
    CB.graph_min    = -15;          CB.graph_max    = 15;
    CB.topo_min     = -7;           CB.topo_max     = 7;
    CB.dyntopo_min  = -1;           CB.dyntopo_max  = 1;
    CB.isotopo_min  = -7;           CB.isotopo_max  = 7;
    CB.restopo_min  = -7;           CB.restopo_max  = 7;
    CB.Ph_min       = 0;         	CB.Ph_max     	= 2;
    CB.C_min        = 1;            CB.C_max        = 5;
    CB.T_min        = 300;          CB.T_max        = 2500;
    CB.Tres_min     = -400;         CB.Tres_max     = 400;
    CB.str_max      = 1000;         CB.str_min      = 0;
    CB.nstr_max     = 800;          CB.nstr_min     = -800;
    CB.edot_max     = 1e-14;        CB.edot_min     = 0;
    CB.eta_min      = 10^(18.8);    CB.eta_max      = 10^(27.8);
    CB.vel_min      = 0;            CB.vel_max      = 10;
    CB.sfun_min   	= -300;       	CB.sfun_max   	= 300;
    CB.rho_min      = 3200;         CB.rho_max      = 3450;
    CB.P_min        = 0;            CB.P_max        = 140;
    CB.geoid_min    = -50;          CB.geoid_max    = 50;
    CB.diss_min     = 0;            CB.diss_max     = 1e-5;
    CB.sfv_min      = 0;            CB.sfv_max      = 1;
    CB.water_min    = 0;            CB.water_max    = 100;
else
    CB.graph_min    = -400;       	CB.graph_max 	= 400;
    CB.topo_min     = -0.05;        CB.topo_max     = 0.05;
    CB.dyntopo_min  = -0.05;        CB.dyntopo_max  = 0.05;
    CB.isotopo_min  = -0.1;         CB.isotopo_max  = 0.1;
    CB.restopo_min  = -0.1;         CB.restopo_max  = 0.1;
    CB.Ph_min       = 0;         	CB.Ph_max     	= 2;
    CB.C_min        = 1;            CB.C_max        = 5;
    CB.T_min        = 0;            CB.T_max        = 1;
    CB.Tres_min     = -0.4;         CB.Tres_max     = 0.4;
    CB.str_max      = 1e5;          CB.str_min      = 0;
    CB.nstr_max     = 1e4;          CB.nstr_min     = -1e4;
    CB.edot_max     = 1e-14;        CB.edot_min     = 0;
    CB.eta_min      = 1e-4;         CB.eta_max      = 1e5;
    CB.vel_min      = 0;            CB.vel_max      = 400;
    CB.sfun_min   	= -12000;     	CB.sfun_max   	= 12000;
    CB.rho_min      = 0.9697;     	CB.rho_max      = 1.0455;
    CB.P_min        = 0;            CB.P_max        = 2e2;
    CB.geoid_min    = -0.05;       	CB.geoid_max    = 0.05;
    CB.diss_min     = 0;            CB.diss_max     = 1;
    CB.sfv_min      = 0;            CB.sfv_max      = 1;
    CB.water_min    = 0;            CB.water_max    = 1;
end
CB.graph_log 	= false;
CB.topo_log 	= false;
CB.dyntopo_log  = false;
CB.isotopo_log  = false;
CB.restopo_log  = false;
CB.Ph_log     	= false;
CB.C_log        = false;
CB.T_log        = false;
CB.Tres_log     = false;
CB.str_log      = false;
CB.nstr_log     = false;
CB.edot_log     = false;
CB.eta_log      = true;
CB.vel_log      = false;
CB.sfun_log     = false;
CB.rho_log      = false;
CB.P_log        = false;
CB.geoid_log    = false;
CB.diss_log     = false;
CB.sfv_log      = false;
CB.water_log    = false;

%% CHECK FOR LOG(10) & ANY INPUT VALUES
name = {'graph' 'topo' 'dyntopo' 'isotopo' 'restopo' 'Ph' 'C' 'T' 'Tres' 'str' ...
    'nstr' 'edot' 'eta' 'vel' 'sfun' 'rho' 'P' 'diss' 'sfv' 'water'};

for ind = 1:length(name)
    %check for input
    if isfield(IN,[name{ind},'_min'])
        if size(IN.([name{ind},'_min']),2)>1 %vector
            CB.([name{ind},'_min']) = IN.([name{ind},'_min'])(1,loop_case);
            CB.multipleMinMax = true;
        else %single value
            CB.([name{ind},'_min']) = IN.([name{ind},'_min']);
        end
    end
    if isfield(IN,[name{ind},'_max'])
        if size(IN.([name{ind},'_max']),2)>1 %vector
            CB.([name{ind},'_max']) = IN.([name{ind},'_max'])(1,loop_case);
            CB.multipleMinMax = true;
        else %single value
            CB.([name{ind},'_max']) = IN.([name{ind},'_max']);
        end
    end
    if isfield(IN,[name{ind},'_log'])
        CB.([name{ind},'_log']) = IN.([name{ind},'_log']);
    end
    
    %take the log10 for min/max values
    if CB.([name{ind},'_log'])
        if CB.([name{ind},'_min'])==0 %prevent min/max values of 0 if log scale is used
            CB.([name{ind},'_min']) = 1;
        end
        if CB.([name{ind},'_max'])==0 %prevent min/max values of 0 if log scale is used
            CB.([name{ind},'_max']) = 1;
        end
        CB.([name{ind},'_max']) = log10(CB.([name{ind},'_max']));
        CB.([name{ind},'_min']) = log10(CB.([name{ind},'_min']));
    end
end


end
