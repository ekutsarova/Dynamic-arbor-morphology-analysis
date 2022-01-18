function outfiles = dynamo_to_swc(filename, outdir, is_register)
% outfiles = dynamo_to_swc(filename)
% filename: path to dynamo mat-file
% is_register: if 1, do register coordinates, default: 1
% outfiles: cell array containing paths to newly created swc-files (one per
% time point)

if nargin<3
    is_register = 1;
end

load(filename, 'savedata')

if is_register
    savedata = register_trees_offset(savedata);
end

n_time = length(savedata.state);

for i_time = 1:n_time
    
    state = savedata.state{i_time};
%     res = [state.info.xres state.info.yres state.info.zres];
    
    % Get all nodes
    
    n_branches = length(state.tree);
    
    all_nodes = zeros(0,3);
    branches = zeros(0);
    
    for i_branch = 1:n_branches
        if ~isempty(state.tree{i_branch})
            nodes = state.tree{i_branch}{1}';
%             nodes = nodes .* repmat(res, size(nodes,1), 1); % apply resolution
%             if size(nodes,1)==1
%                 asdf
%             end
            all_nodes = vertcat(all_nodes, nodes);
            branch = ones(size(nodes,1),1)*i_branch;
            branches = vertcat(branches, branch);
        end
    end
    
    % Find branchpoints and their duplicate
    
    branchpoints = find(diff(branches))+1;
    for i_bpt = 1:length(branchpoints)
        branchpoint = branchpoints(i_bpt);
        [mins,tmp] = sort(abs(sum(all_nodes - all_nodes(branchpoint,:),2)));
        matchpoint = setdiff(tmp(mins==0), branchpoints(i_bpt));
        if length(matchpoint)>1
            matchpoint = min(matchpoint);
        end
%         i_time
%         i_bpt
%         matchpoint
        branchpoints(i_bpt,2) = matchpoint;
    end
    
    % Build swc
    
    % add first branch
    first_branch = all_nodes(1:branchpoints(1,1)-1,:);
    swc = vertcat(1:size(first_branch,1), ... % index
        zeros(1,size(first_branch,1)),... % structure identifier
        first_branch', ... % coordinates
        ones(1,size(first_branch,1)),... % radius
        [-1, 1:size(first_branch,1)-1])'; % parents
    
    for i_bpt = 1:size(branchpoints,1)
        
        parent = branchpoints(i_bpt,2);
        parent = parent - sum(parent>branchpoints(:,1));
        % find nodes to add (removing branch point)
        if i_bpt==size(branchpoints,1)
            node_stop = size(all_nodes,1);
        else
            node_stop = branchpoints(i_bpt+1)-1;
        end
        nodes = all_nodes(branchpoints(i_bpt)+1:node_stop,:);
        n_nodes = size(nodes,1);
        % collect parents (branch point + continuation points if applicable)
        if n_nodes>1
            parents = [parent size(swc,1)+1:size(swc,1)+n_nodes-1];
        else
            parents = parent;
        end
        
        new_swc = vertcat(size(swc,1)+1:size(swc,1)+n_nodes,... % index
            zeros(1,n_nodes),... % ?
            nodes',... % coordinates
            ones(1,n_nodes),... % radius
            parents); % parents
        swc = vertcat(swc, new_swc');
        
    end
    
    % Save swc
    
    [~,basename,~] = fileparts(filename);
    outfile = [outdir '/' basename '_timepoint' num2str(i_time,'%02d') '.swc'];
    fileID = fopen(outfile,'w');
    for i_line = 1:size(swc,1)
        fprintf(fileID, '%d %d %.12f %.12f %.12f  %.6f %d\n', swc(i_line,:));
    end
    fclose(fileID);
    
    outfiles{i_time} = outfile;
    
end

end