function segments = find_segments_trees(axon_tree)

node_type = typeN_tree(axon_tree);
segments = ones(size(node_type));

i_segment = 1;

for i_node=2:size(node_type,1)
    if node_type(i_node)==1
        assert(axon_tree.dA(i_node+1,i_node)==1)
        segments(i_node) = i_segment;
    else
        if node_type(i_node)==2
            assert(axon_tree.dA(i_node+1,i_node)==1)
        end
        segments(i_node) = i_segment;
        i_segment = i_segment + 1;
    end
end

end