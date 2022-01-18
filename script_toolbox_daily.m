% This code is based on trees_toolbox (Group of Dr. Hermann Cuntz). 
% I added a function to be able to calculate the outer segment length based
% on Strahler analysis available in trees_toolbox
clear
addpath('.\functions\');
% addpath to the current stable version of trees_toolbox (download from https://www.treestoolbox.org/ )
% addpath to the code from Bird and Cuntz, 2019 if not included in the
% current verison of trees_toolbox, necessary for arbour spanning volume
% calculation

celltable = readtable('Conversions_all_daily.xlsx'); %conversion .csv and the proper .swc(if necessary converted to um) should be in the same folder and you should b
celltable = rmmissing(celltable(:,1:6));
cellnames = celltable.CellNames;
NumCells = size(cellnames,1);

start_trees
%%
Num_TPoints = zeros(NumCells,1);
Tightfit_Volume = zeros(NumCells,1);
All_Len_Strahler1 = zeros(NumCells,1);
All_Len_Rest = zeros(NumCells,1);

All_Histseg_Strah1 = zeros(NumCells,4);
All_Histseg_Strah1_Norm = zeros(NumCells,4);

Strahler1_Lens_Seg = {};
Strahler_All_Lens_Seg = {};

for i_axon=1:NumCells
    i_axon
    new_cellname = ['daily_SWCs\' strtrim(strrep(cellnames{i_axon},'.swc','coorUm.swc'))];   
    axon_tree = load_tree(new_cellname);
% The resampled traces were used to visualize the vertices and edges of the
% spanning volume of the axon together with the arbor. To get the resampled
% traces, use the three lines below: 
%     resam_axon = resample_tree(axon_tree,0.05,'-b''-l''-r''-v'); 
%     cellname_swc = strrep(new_cellname,'.swc','resampled.swc');
%     swc_tree(resam_axon, cellname_swc);
    
    overall_stats = stats_tree(axon_tree);
    terminal_points = T_tree(axon_tree);
    num_tpoints = sum(terminal_points);
    convex_index = convexity_tree(axon_tree);
    tight_fit = boundary_tree(axon_tree); 
    
    node_len = len_tree(axon_tree);
    node_strahler = strahler_tree(axon_tree);
    segments = find_segments_trees(axon_tree); % get the segments from the nodes
    
    num_segments = max(segments(:,1));
    Seg_Strah_Len = zeros(num_segments,2);
    
    for i_segment=1:num_segments
        i_segment;
        index = segments==i_segment;
        segment_len = sum(node_len(index));
        segment_strahler = max(node_strahler(index));
        Seg_Strah_Len(i_segment,1) = segment_len;
        Seg_Strah_Len(i_segment,2)= segment_strahler;
    end
    
    [seg_count,strah_num]= hist(Seg_Strah_Len(:,2),unique(Seg_Strah_Len(:,2)));
    
    %Number of segments of each strahler number
    num_strahler1 = seg_count(1,1);
    num_rest = sum(seg_count(2:end),2);   
    
    Strahler_Lens = zeros(max(strah_num),1);
    Strahler_Lens_Seg ={};
    
    for i_strah_num = 1:max(strah_num)
        index = Seg_Strah_Len(:,2)==i_strah_num;
        strahler_len_seg = Seg_Strah_Len(index);
        strahler_len = sum(strahler_len_seg,1);
        Strahler_Lens(i_strah_num,1) = strahler_len;
        Strahler_Lens_Seg{i_strah_num} = strahler_len_seg;
    end 
    
    %Save over all cells
    % Number of terminal points and spanning field volume
    Num_TPoints(i_axon,1) = num_tpoints;
    Tightfit_Volume(i_axon,1) = tight_fit.V;
    
    % Total length of each strahler segment type
    len_strahler1 = Strahler_Lens(1);
    len_rest = sum(Strahler_Lens(2:end),1);

    % Histograms of segment lengths
    histseg = histcounts(Seg_Strah_Len(:,1),[0,5,10,20,inf]);
    histseg_norm = histseg/sum(histseg);
   
    % Histograms of segment lengths of strahler 1 segments (outer segments)
    histseg_strah1 = histcounts(Strahler_Lens_Seg{1},[0,5,10,20,inf]);
    histseg_norm_strah1 = histseg_strah1/sum(histseg_strah1);
       
    All_Len_Strahler1(i_axon,1) = len_strahler1;
    All_Len_Rest(i_axon,1) = len_rest;
    
    All_Histseg_Strah1(i_axon,:) = histseg_strah1;
    All_Histseg_Strah1_Norm(i_axon,:) = histseg_norm_strah1;
    
    
    
    save([new_cellname 'tight_fit.mat'],'-struct','tight_fit');    close all
end
%% Saves a table for all the cells
new_celltable = addvars(celltable, Num_TPoints,...
     All_Len_Rest, All_Histseg_Strah1_Norm, Tightfit_Volume,'After','Group');
writetable(new_celltable,'Daily_Table_20220118_bigtable.xlsx') ;

%% Normalizes variables to a specific day. For example Day1 or Day4, whichever makes sense
celltable_sorted = readtable('Daily_Table_20220118_bigtable.xlsx');
group_names = unique(celltable_sorted.Group(:));
group_names = {group_names{1},group_names{3}, group_names{2}}; %ordered 
num_days = max(unique(celltable_sorted.Day(:)));
num_groups = size(group_names,2);%how many groups there are
index_day1 = celltable_sorted.Day == 1;
num_real_var = size(celltable_sorted,2)-6;
Norm_var_vec_d1 = zeros (NumCells, num_real_var);
Var_names = {};
celltable_sorted_norm_d1 = celltable_sorted;

for i_var = 7:size(celltable_sorted,2)
    var_name = celltable_sorted.Properties.VariableNames{i_var};
    variable_vec = celltable_sorted.(var_name);
    i_var
    for i_group=1:num_groups
        index_group = ismember(celltable_sorted.Group,group_names{i_group});
        index_day1_group = index_group & index_day1;
        var_vec_d1 = variable_vec(index_day1_group);
        
        for i_day = 1:num_days
            index_i_day = celltable_sorted.Day == i_day;
            index_i_day_group = index_group & index_i_day;
            var_vec_i_day =  variable_vec(index_i_day_group);
            norm_d1_var_vec_group = var_vec_i_day./var_vec_d1;
            Norm_var_vec_d1(index_i_day_group,(i_var-6))= norm_d1_var_vec_group;
        end 
    end
    
    norm_d1_var_vec = Norm_var_vec_d1(:,i_var-6);
    if ~any(isnan(norm_d1_var_vec))&& ~any(isinf(norm_d1_var_vec))
        celltable_sorted_norm_d1 = addvars(celltable_sorted_norm_d1, norm_d1_var_vec,...
            'NewVariableNames', ['Norm' var_name]);
        Var_names{i_var-6} = var_name;
    end  
end
writetable(celltable_sorted_norm_d1, 'Daily_Table_20220118_matlab_sorted_norm_d1_variables.csv')