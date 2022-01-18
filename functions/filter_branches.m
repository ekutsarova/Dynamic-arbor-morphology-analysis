% Outputs the identity of the branches added and lost, as well as the total
% number of branches added, the total elongation
% (sum of positive change in length), total retraction (negative change in
% length), number of branches, total of the arbor per timepoint based on the function from Dynamo "get_added_subtracted". The only addition is that only branches which attained filter_length value throughout their lifetime,
% as well as elongation/branch>filter_elong are being output. 
% 
function [added,subtracted,Addedsum,Lostsum,Elongation,Retraction,...
    NumBranches,Length,tipLengths] = filter_branches(cell,filter_length, filter_elong)

if ischar(cell)
    dynamo_mat = load(cell);
    savedata = dynamo_mat.savedata;
else
    savedata = cell;
end

opts.filo_dist=10;
opts.term_dist=10;
opts.exc_basal=1;
opts.exc_axon=1;
opts.noplot=1;

state = savedata.state;

[tipLengths, ~, ~] = get_added_subtracted(savedata, opts);
tipLengths(tipLengths==0) = nan;

% Filter branches
for i_branch = 1:size(tipLengths,3)
%     i_branch
    if ~any(tipLengths(1,:,i_branch)>filter_length)
%         tipLengths(1,:,i_branch) = 0; % This was wrong
        tipLengths(1,:,i_branch) = nan;
    end
    FilteredLengths(1,:,i_branch) = tipLengths(1,:,i_branch);
end

%calculate added and subracted
filos = ~isnan(FilteredLengths);
added = filos(:, 2:end, :) & ~filos(:, 1:end-1, :);
subtracted = filos(:, 1:end-1, :) & ~filos(:, 2:end, :);

% change in length
lengths = FilteredLengths;
lengths(isnan(lengths)) = 0;
changes = lengths(:, 2:end, :) - lengths(:, 1:end-1, :);
% filter out changes smaller than filter_elong
changes(abs(changes)<filter_elong) = 0;

%sum added/timepoint
Timepoints = length(FilteredLengths(1,:,1))-1;
Addedsum = zeros(Timepoints,1);
Lostsum = zeros(Timepoints,1);
Elongation = zeros(Timepoints,1);
Retraction = zeros(Timepoints,1);
for i_time = 1:Timepoints
%     i_time
    Addedsum(i_time,1) = sum(added(:,i_time,:));
    Lostsum(i_time,1) = sum(subtracted(:,i_time,:));
    Elongation(i_time,1) = sum(changes(:,i_time,changes(:,i_time,:)>0));
    Retraction(i_time,1) = -sum(changes(:,i_time,changes(:,i_time,:)<0));
end
Length = sum(lengths,3)';
NumBranches = sum(filos,3)';

end
%%
