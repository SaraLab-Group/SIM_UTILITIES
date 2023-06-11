function S = sphereVolume(num_spheres,Box_Dims,radius)
%spherePlane This function will generate random spheres in a fixed central plane
%   Curently limited xy coordinates to fit 513cu.pix region
rng(0,'twister');
rnd_x = int16(randi([-128,128],1,num_spheres*10000));
rnd_y = int16(randi([-128,128],1,num_spheres*10000));
rnd_z = int16(randi([-64,64],1,num_spheres*10000));

lim = double(Box_Dims/2);
coord_val = -lim:lim;

matches = 1;
good_x = int16(zeros(num_spheres));
good_x(1) = rnd_x(1);
good_y = int16(zeros(num_spheres));
good_y(1) = rnd_y(1);
good_z = int16(zeros(num_spheres));
good_z(1) = rnd_z(1);
rnd_idx = 2;
center = floor(Box_Dims/2);
while (matches < num_spheres && rnd_idx <= num_spheres * 100)
    rndx_shift = center + rnd_x(rnd_idx);
    rndy_shift = center + rnd_y(rnd_idx);
    rndz_shift = center + rnd_z(rnd_idx);
    %z_shift = center + good_z;
    dist = 0;
    for ii = 1:(matches + 1)
        gx_shift = center + good_x(ii);
        gy_shift = center + good_y(ii);
        gz_shift = center + good_z(ii);
        dist = floor(sqrt((coord_val(rndx_shift) - coord_val(gx_shift))^2 + (coord_val(rndy_shift) - coord_val(gy_shift))^2 + (coord_val(rndz_shift) - coord_val(gz_shift))^2));
        if(floor(dist) - 3 <= radius)
           %fprintf('!');
           break 
        end
    end
    if floor(dist) - 3 > radius
       matches = matches + 1;
       good_x(matches) = rnd_x(rnd_idx);
       good_y(matches) = rnd_y(rnd_idx);
       good_z(matches) = rnd_z(rnd_idx);
       %fprintf('.');
    end
    rnd_idx = rnd_idx + 1;
   
end
fprintf("\nFound %d Good Coords\n", matches);
S = zeros(Box_Dims,Box_Dims,Box_Dims);
[ss, dims] = centered_small_Sphere(radius);
for ii = 1:matches
   S = (insertSphere(S, ss, radius,good_y(ii), good_x(ii), good_z(ii))); 
end
S(S>1) = 1;
S = imgaussfilt3(S);
end

