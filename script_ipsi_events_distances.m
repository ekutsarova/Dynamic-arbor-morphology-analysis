%% This code is based on Dynamo(version 2013), written in MATLAB (Drs. Kasper Podgorski and Kurt Haas) and trees_toolbox (Group of Dr. Hermann Cuntz). 
% I added a few functions in order to be able to calculate event pair
% distances and perform Monte Carlo simulations of the events on the axonal
% arbor (the script is associated with Kutsarova et al, 2021 bioRxiv) 
clear

addpath('.\functions\'); 
% addpath to the current stable version of trees_toolbox (download from https://www.treestoolbox.org/ )

start_trees % initialize trees_toolbox

binsheet = 'binsheet_oneandahalfhours.xlsx'; %specifies which timepoints should be taken to obtain a mean value over the conditions (dark, asynch, synch)
cellname = '20161213ipsiYellowMO_1-31Trial2'; % filename of the tracing from Dynamo without the .mat extension
group = 'MO'; % whether the cell is from any of the morpholino groups (Control, p75 or TrkB) or 'Fc' for the cells from the TrkB-Fc group
bintable = readtable(binsheet);

% find unique conditions
times_to_include = ~cellfun(@isempty, bintable.Conditions_MO);
all_conditions = unique(bintable.Conditions_MO(times_to_include));
NumConds = length(all_conditions);

%Values for input variables
FilterLength = 1.5; % minimum attained legnth (micrometers) throughout the lifetime of a newly added branch to be considered for further calculation and not just tracing error
MC_num = 100; % Monte Carlo number for randomization of the events
FilterElong = 2; % minimum change in length (micrometers) per branch per timepoint which is considered towards total elongation or retration (not used in this script but are output of filter_branches)  
%%

swc_files = dynamo_to_swc(['Matfiles/' cellname '.mat'],'ipsi_SWCs');% converts the .mat files traced in Dynamo to .swc files (one for each timepoint)

% The three rows below resample the nodes on the traced axon, so they are equidistant at 0.15 um; trees_toolbox function and save the resampled swc files (this takes a lot of time and should be performed once)
%     resam_axon = resample_tree(swc_files{i_axon},0.15,'-b''-l''-r''-v');
%     cellname_swc = [cellnames{i_axon} 'resampled_0_15.swc'];
%     swc_tree(resam_axon, cellname_swc);

resampled_swcs = strrep(swc_files, '.swc', 'resampled_0_15.swc'); % once the resampled swc files have been saved, this line can be used to call the filenames

% Load Dynamo-traced files
axon = ['Matfiles/' cellname '.mat'];
load(axon, 'savedata');
savedata = register_trees_offset(savedata); %registration based on manual landmarking in Dynamo (The new coordinates are taken straight from Dynamo)

%Get added and lost
[added,lost,Addedsum,Lostsum,~,~,~,~,~] = filter_branches(savedata, FilterLength, FilterElong);

% Binning of dark, synch, asynch based on whether the axon is from 'MO' (morpholino) or 'Fc' (TrkB-Fc) group
if group == 'Fc'
    conditions = bintable.Conditions_Fc(2:end);
else
    conditions = bintable.Conditions_MO(2:end);
end

%Clump additions and losses as events
events_sum = [Addedsum, Lostsum];
events = {added,lost};
event_names = {'added','lost'};
EventCoords = {{}, {}};

%Extract event coordinates and simulate events on the arbor; calculate pair
%distances between real or between simulated events 
for type_events=1:2
    event_all_timepoints = events_sum(:,type_events);
    EventCoords{type_events} = get_event_coords (events{type_events}, savedata, type_events);
    [rand_coords_all] = random_event_generator(event_all_timepoints, resampled_swcs, MC_num, type_events);
    
    event_coords = EventCoords{type_events};
    
    for i_cond = 1:NumConds
        cond_idx = find(strcmp(conditions, all_conditions{i_cond}));
        [~,p_dist, ~, ~,~]=find_event_dist(event_coords, cond_idx);
        events_condition = mean(event_all_timepoints(cond_idx));
        
        for i_monte_carlo = 1:MC_num
            rand_event_coords = rand_coords_all(:,i_monte_carlo);
            [~,rand_p_dist, ~,~, ~] = find_event_dist(rand_event_coords, cond_idx);
            
            MC_p_dist(i_monte_carlo,1) = rand_p_dist;
            mean_theor_p_dist = mean (MC_p_dist,1);
        end
        
        events_all(i_cond,type_events) = events_condition;
        p_dist_all(i_cond,type_events) = p_dist;
        mean_theor_p_dist_all(i_cond,type_events) = mean_theor_p_dist;
        
        R_pdist = p_dist/mean_theor_p_dist;
        R_pdist_all (i_cond,type_events) = R_pdist;
    end
end
%% Normalizes to dark  and saves tables into a .csv file

output_table = table(all_conditions);

norm_dark_pdist_all = p_dist_all(:,:)./p_dist_all(1,:);
norm_dark_theor_p_dist = mean_theor_p_dist_all(:,:)./mean_theor_p_dist_all(1,:);
norm_dark_R_pdist_all  = R_pdist_all(:,:)./R_pdist_all(1,:);
norm_events_sum = events_all(:,:)./events_all(1,:); 

variables = {norm_dark_pdist_all,norm_dark_theor_p_dist, norm_dark_R_pdist_all, norm_events_sum};
variable_names = {'norm_dark_pdist','norm_dark_theor_p_dist','norm_dark_R_pdist', 'norm_events_sum'};


for i_type = 1:size(event_names,2)
    for i_var=1:size(variables,2)
        var_to_save = variables{i_var}(:,i_type);
        output_table = addvars(output_table, var_to_save, 'NewVariableNames', [variable_names{i_var} event_names{i_type}]);
    end
end
    writetable (output_table, ['ipsi_pair_distances' cellname '_' binsheet '.csv'])