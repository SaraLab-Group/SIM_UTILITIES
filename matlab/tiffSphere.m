function S = tiffSphere(BOX_DIMS, radius, y_off, x_off, z_off)
%UNTITLED3 Summary of this function goes here
%   pretty Self Explanitory
%   generates a rasterize sphere which can be saved as a tiff stack.

    center = uint16(BOX_DIMS/2 + 1);
    S = zeros(BOX_DIMS,BOX_DIMS,BOX_DIMS);
    small_dims = radius * 2 + 1;
    lim = double(small_dims/2);
    coord_val = -lim:lim;
    sBOX = zeros(small_dims, small_dims, small_dims);
    center_small = uint16(small_dims/2);
%     x_shift = center_small + x_off;
%     y_shift = center_small + y_off;
%     z_shift = center_small + z_off;
    x_shift = center_small;
    y_shift = center_small;
    z_shift = center_small;
    for xx = 1:small_dims
        for yy = 1:small_dims
            for zz = 1:small_dims
                dist1 = (double((coord_val(xx)-coord_val(x_shift))^2+(coord_val(yy)-coord_val(y_shift))^2+(coord_val(zz)-coord_val(z_shift))^2))^(1/2);
                %fprintf("%f ", dist);
                if(dist1 <= radius)
                    sBOX(xx,yy,zz) = 1;
                    %fprintf("%d %d %d ", xx, yy, zz);
                end
            end
        end
    end
    xb_l = center + x_off - radius;
    xb_r = center + x_off + radius;
    yb_l = center + y_off - radius;
    yb_r = center + y_off + radius;
    zb_l = center + z_off - radius;
    zb_r = center + z_off + radius;
    xs_l = 1;
    xs_r = 2*radius + 1;
    ys_l = 1;
    ys_r = 2*radius + 1;
    zs_l = 1;
    zs_r = 2*radius + 1;
    if xb_l < 1
        xs_l = xs_l - xb_l + 2;
        xb_l = 1;
    end
    if xb_r > BOX_DIMS
        xs_r = xs_r - (xb_r - BOX_DIMS);
        xb_r = BOX_DIMS;
    end
    if yb_l < 1
        ys_l = ys_l - yb_l + 2;
        yb_l = 1;
    end
    if yb_r > BOX_DIMS
        ys_r = ys_r - (yb_r - BOX_DIMS);
        yb_r = BOX_DIMS;
    end
    if zb_l < 1
        zs_l = zs_l - zb_l + 2;
        zb_l = 1;
    end
    if zb_r > BOX_DIMS
        zs_r = zs_r - (zb_r - BOX_DIMS);
        zb_r = BOX_DIMS;
    end
%     fprintf("xb_l: %d xb_r: %d, yb_l: %d yb_r: %d, zb_l: %d zb_r: %d\n", xb_l, xb_r, yb_l, yb_r, zb_l, zb_r);
%     fprintf("xs_l: %d xs_r: %d, ys_l: %d ys_r: %d, zs_l: %d zs_r: %d\n", xs_l, xs_r, ys_l, ys_r, zs_l, zs_r);
    S(xb_l:xb_r,yb_l:yb_r,zb_l:zb_r) = imgaussfilt3(sBOX(xs_l:xs_r,ys_l:ys_r,zs_l:zs_r));
end