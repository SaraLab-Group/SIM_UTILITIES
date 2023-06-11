function S = filamentVolume(num_filaments, filament_length,Box_Dims,radius,parallel_fil, spacing)
%spherePlane This function will generate random spheres in a fixed central plane
%   Curently limited xy coordinates to fit 513cu.pix region
rng shuffle
rnd_x = int16(randi([-128,128],1,num_filaments*10000));
rnd_y = int16(randi([-128,128],1,num_filaments*10000));
rnd_z = int16(randi([-64,64],1,num_filaments*10000));

lim = double(Box_Dims/2);
coord_val = -lim:lim;

matches = 1;
good_x = int16(zeros(num_filaments));
good_x(1) = rnd_x(1);
good_y = int16(zeros(num_filaments));
good_y(1) = rnd_y(1);
good_z = int16(zeros(num_filaments));
good_z(1) = rnd_z(1);
rnd_idx = 2;
center = floor(Box_Dims/2);
while (matches < num_filaments && rnd_idx <= num_filaments * 100)
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
rng shuffle
for ii = 1:matches
   prev_x = good_x(ii);
   prev_y = good_y(ii);
   prev_z = good_z(ii);
   sign_xy = (-1)^(randi([0,100]));
   sign_z = (-1)^(randi([0,100]));
   rand_degxy = randi([0,359]);
   rand_degz = randi([-12,12]);
   if parallel_fil
       par_x1 = prev_x + sign_xy*floor(spacing*radius*cosd(rand_degxy + 90));%*cosd(rand_degz + 90));
       par_y1 = prev_y + sign_xy*floor(spacing*radius*sind(rand_degxy + 90));%*cosd(rand_degz + 90));
       %par_z1 = prev_z + floor(2*radius*sind(rand_degz + 90));
       S = (insertSphere(S, ss, radius, par_y1, par_x1, prev_z));
       par_x2 = prev_x + sign_xy*floor(spacing*radius*cosd(rand_degxy - 90));%*cosd(rand_degz - 90));
       par_y2 = prev_y + sign_xy*floor(spacing*radius*sind(rand_degxy - 90));%*cosd(rand_degz - 90));
       %par_z2 = prev_z + floor(spacing*radius*sind(rand_degz - 90));
       S = (insertSphere(S, ss, radius, par_y2, par_x2, prev_z));
   else
       S = (insertSphere(S, ss, radius,good_y(ii), good_x(ii), good_z(ii))); 
   end
   for kk = 1:filament_length
       prev_x = prev_x + sign_xy*floor(radius*cosd(rand_degxy)*cosd(rand_degz));
       prev_y = prev_y + sign_xy*floor(radius*sind(rand_degxy)*cosd(rand_degz));
       prev_z = prev_z + floor(radius*sind(rand_degz));
       
       if parallel_fil
           par_x1 = prev_x + sign_xy*floor(spacing*radius*cosd(rand_degxy + 90));%*cosd(rand_degz + 90));
           par_y1 = prev_y + sign_xy*floor(spacing*radius*sind(rand_degxy + 90));%*cosd(rand_degz + 90));
           %par_z1 = prev_z + floor(2*radius*sind(rand_degz + 90));
           S = (insertSphere(S, ss, radius, par_y1, par_x1, prev_z));
           par_x2 = prev_x + sign_xy*floor(spacing*radius*cosd(rand_degxy - 90));%*cosd(rand_degz - 90));
           par_y2 = prev_y + sign_xy*floor(spacing*radius*sind(rand_degxy - 90));%*cosd(rand_degz - 90));
           %par_z2 = prev_z + floor(2*radius*sind(rand_degz - 90));
           S = (insertSphere(S, ss, radius, par_y2, par_x2, prev_z));
       else
           S = (insertSphere(S, ss, radius, prev_y, prev_x, prev_z));
       end

       rand_degxy = randi([rand_degxy - 22,rand_degxy + 23]);% This should only allow next segment 45deg freedom
       rand_degz = randi([rand_degz - 12, rand_degz + 12]);
   end
end
S(S>1) = 1;
S = imgaussfilt3(S);
end

