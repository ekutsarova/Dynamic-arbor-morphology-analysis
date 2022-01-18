function savedata_reg = register_trees_offset(savedata)

states = savedata.state;

% Grid2world
g2w = eye(4);
g2w(1,1) = states{1}.info.xres;
g2w(2,2) = states{1}.info.yres;
g2w(3,3) = states{1}.info.zres;

% Register with offset
states_reg = states;
for i_time = 1:length(states)
    aff = eye(4);
    aff(1:3,4) = -states{i_time}.info.offset;
    % here I'm switching g2w and aff, since transform to world coordinates
    % has to come after applying the offset
    states_reg{i_time}.tree = register_tree(states{i_time}.tree, g2w, aff);
    states_reg{i_time}.info.xres = 1;
    states_reg{i_time}.info.yres = 1;
    states_reg{i_time}.info.zres = 1;
end

savedata_reg = savedata;
savedata_reg.state = states_reg;

end