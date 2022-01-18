%% This is a function from dynamo used to extract branch lengths,and subsequently branch additions/retractions.
function [tipLengths, added, subtracted] = get_added_subtracted (drawings, opts)
%A DYNAMO ANALYSIS FUNCTION
%All Dynamo analysis functions take two arguments:
%   drawing:    a Dynamo saved structure containing the fields 'states' and 'traits'
%   opts:       contains options for the function

%filo types key:
% 0 - absent
% 1 - interstitial filo
% 2 - terminal filo
% 3 - branch w int filo
% 4 - branch w term filo
% 5 - branch only

%branches_added_subtracted
%Ignores the definition of filopodia, and tracks the addition and
%subtaction of terminal projections (i.e. branchtips), as well as their
%lengths

if ~isfield(opts, 'term_dist') %distance from end to be considered terminal filopodia
    opts.term_dist = 10;
end
%accumulate
tipLengths = nan(length(drawings), size(drawings(1).state,length(drawings(1).state{end}.tree)));
for drawing = 1:length(drawings)
    for time = 1:size(drawings(drawing).state,2)
        state = drawings(drawing).state{time};
        for branchNumber = 1:size(state.tree, 2);
            branch = state.tree{branchNumber};
            
            %if the branch doesn't exist
            if isempty(branch)
                tipLengths(drawing, time, branchNumber) = nan;
                continue
            end
            
            %find the last node with children
            haschildren = find(~cellfun(@isempty, branch{2}));
            lastBP = find(haschildren, 1, 'last');
            if isempty(lastBP)
                lastBP = 1;
            end
            if lastBP==length(branch{2}) %if this isn't an endpoint
                tipLengths(drawing, time, branchNumber) = nan;
            end
            
            diffs = diff((diag([state.info.xres state.info.yres state.info.zres])*branch{1})');
            dists = sqrt(sum(diffs.^2,2)); %euclidean distance
            
            %tiplength is sum of distances from last node to end
            tipLengths(drawing, time, branchNumber) = sum(dists(lastBP:end));
        end
    end
end

%calculate added
filos = ~isnan(tipLengths);
added = filos(:, 2:end, :) & ~filos(:, 1:end-1, :);
subtracted = filos(:, 1:end-1, :) & ~filos(:, 2:end, :);
end