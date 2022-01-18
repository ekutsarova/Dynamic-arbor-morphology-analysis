function [knn_dist,p_dist, median_knn_d,median_p_d,p_distances]=find_event_dist(coords, cond_idx)
if nargin>1
    
    events_coord_cond = cat(1,coords{cond_idx});
    [~,knn_distances] = knnsearch(events_coord_cond, events_coord_cond,'k', 2);
    
    if size(knn_distances,2)>1
        knn_dist = mean(knn_distances(:,2));
        median_knn_d = median(knn_distances(:,2));
    else
        knn_dist = nan;
        median_knn_d = nan;
    end
    p_distances = pdist(events_coord_cond);
    p_dist = mean(p_distances);
    median_p_d = median(p_distances);
else
    timepoints = size(coords,1);
    knn_dist = zeros(timepoints,1);
    median_knn_d = zeros(timepoints,1);
    p_dist = zeros(timepoints,1);
    median_p_d = zeros(timepoints,1);
    
    for i_time =1:timepoints
        p_distances_i_time = pdist(coords{i_time});
        p_dist_i_time = mean(p_distances_i_time);
        median_p_d_i_time = median(p_distances_i_time);
        [~,knn_distances_i_time] = knnsearch(coords{i_time}, coords{i_time},'k', 2);
        if size(knn_distances_i_time,2)>1
            knn_dist_i_tme = mean(knn_distances_i_time(:,2));
            median_knn_d_i_time = median(knn_distances_i_time(:,2));
        else
            knn_dist_i_tme = nan;
            median_knn_d_i_time = nan;
            
        end
        knn_dist(i_time) = knn_dist_i_tme;
        median_knn_d(i_time) = median_knn_d_i_time;
        p_dist(i_time) = p_dist_i_time;
        median_p_d(i_time) = median_p_d_i_time;
    end
end
end