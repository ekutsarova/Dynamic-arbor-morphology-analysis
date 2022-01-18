function [event_coords] = get_event_coords (events, cell, event_type)

event_coords = {};


if ischar(cell)
    dynamo_mat = load(cell);
    savedata = dynamo_mat.savedata;
else
    savedata = cell;
end

for i_time = 1:size(events,2)
    
    branch_idx = find(squeeze(any(events(:,i_time,:)>0,2)));
    
    event_coords_branch = zeros(length(branch_idx),3);
      
    for i_branch = 1:length(branch_idx)
        
        
        if event_type ==1
            coord = savedata.state{i_time+1}.tree{branch_idx(i_branch)}{1}(:,1) ;
            
        else
            coord = savedata.state{i_time}.tree{branch_idx(i_branch)}{1}(:,end);
        end
        event_coords_branch(i_branch,:)=coord;
    end
    event_coords{i_time} = event_coords_branch;
    
end
event_coords = event_coords';
end