function coords_new = apply_transform(coords, affine)

coords(4,:) = 1;
coords_new = affine * coords;
coords_new = coords_new(1:3,:);

end