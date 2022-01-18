% This function takes n random coordinates on the cell for each timepoint of dynamic imaging, where an event - addition or subtraction can occutr. n = the number of events at the timepoint of interes...
% events_sum is an output of filter_branches function and gives the number
% of events (either additions or losses)
%resam_trees is a list of resampled .swc created with dynamo_to_swc and
%further resampled using resample_tree down to 0.15
%MC_num specifies the number of randon draws for coordinates (Monte Carlo
%number)
% type_events specifies addition = 1; elimination = 2; For additions all of the points are considered: branch, continuation and terminal.
%For subtraction, only terminal points are considered



function [rand_coords_all] = random_event_generator(events_sum, resam_trees, MC_num, type_events)

rand_coords_all = {};

for i_time = 1:size(resam_trees,2)-1
    
    if type_events==1
        resam_tree = load_tree(resam_trees{i_time+1});
        %     index_BCpoints = find(typeN_tree(resam_tree)~=0);
        %      X = resam_tree.X(index_BCpoints);
        %     Y = resam_tree.Y(index_BCpoints);
        %     Z = resam_tree.Z(index_BCpoints);
        resam_coords = cat(2,resam_tree.X,resam_tree.Y, resam_tree.Z);
    else
        resam_tree = load_tree(resam_trees{i_time});
        index_Tpoints = find(typeN_tree(resam_tree)==0); % asks which points on the list of points are terminal == 0; and finds their index in the list
        X = resam_tree.X(index_Tpoints);
        Y = resam_tree.Y(index_Tpoints);
        Z = resam_tree.Z(index_Tpoints);
        resam_coords = cat(2, X, Y, Z);
        %  resam_coords = cat(2,resam_tree.X,resam_tree.Y, resam_tree.Z);
    end
    
    num_events = events_sum(i_time);
    
    for i_monte_carlo = 1:MC_num
        %         for i_timepoint=1:size(events_sum,1)
        %         if remove==i_timepoint
        %             fprintf('timepoint removed for Monte Carlo..' );
        %         else
        sample = randsample(size(resam_coords,1), num_events);
        rand_coords = resam_coords(sample,:);
        %             rand_coords = datasample(resam_coords, num_events, 'Replace', false); %picks n random coordinates from the resampled tree, where n is the number of events(additions or subractions for timepoint i_time)
        %         end
        rand_coords_all{i_time,i_monte_carlo} = rand_coords;
        %         end
    end
end
end

