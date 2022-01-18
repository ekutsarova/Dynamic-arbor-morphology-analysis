function tree_reg = register_tree(tree, aff, g2w)

tree_reg = tree;
for i_branch = 1:length(tree)
    if ~isempty(tree{i_branch})
        coord = apply_transform(tree{i_branch}{1}, g2w);
        coord_reg = apply_transform(coord, aff);
        tree_reg{i_branch}{1} = coord_reg;
    end
end

end