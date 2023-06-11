function S = sphereFilamentPlane(num_filaments, filament_length,Box_Dims,radius,parallel_fil)
%spherePlane This function will generate random spheres in a fixed central plane
%   Curently limited xy coordinates to fit 513cu.pix region
rng(0,'twister');
rnd_x = randi([-128,128],1,num_filaments*10000);
rnd_y = randi([-128,128],1,num_filaments*10000);
rnd_z = 0;

lim = double(Box_Dims/2);
coord_val = -lim:lim;

matches = 1;
good_x = (zeros(num_filaments));
good_x(1) = rnd_x(1);
good_y = (zeros(num_filaments));
good_y(1) = rnd_y(1);
good_z = rnd_z;
rnd_idx = 2;
center = floor(Box_Dims/2);
while (matches < num_filaments && rnd_idx <= num_filaments * 100)
    rndx_shift = center + rnd_x(rnd_idx);
    rndy_shift = center + rnd_y(rnd_idx);
    z_shift = center + good_z;
    dist = 0;
    for ii = 1:(matches + 1)
        gx_shift = center + good_x(ii);
        gy_shift = center + good_y(ii);
        dist = floor(sqrt((coord_val(rndx_shift) - coord_val(gx_shift))^2+ (coord_val(rndy_shift) - coord_val(gy_shift))^2));
        if(floor(dist) - 1 <= radius)
           %fprintf('!');
           break 
        end
    end
    if floor(dist) - 1 > radius
       matches = matches + 1;
       good_x(matches) = rnd_x(rnd_idx);
       good_y(matches) = rnd_y(rnd_idx);
       %fprintf('.');
    end
    rnd_idx = rnd_idx + 1;
   
end
fprintf("\nFound %d Good Coords\n", matches);
S = zeros(Box_Dims,Box_Dims,Box_Dims);
[ss, dims] = centered_small_Sphere(radius);
for ii = 1:matches

   prev_x = good_x(ii);
   prev_y = good_y(ii);
   rand_deg = randi([0,359]);
%    rand_deg_prev = randi([rand_deg - 22,rand_deg + 23]);
%    rand_deg_next = randi([rand_deg - 22,rand_deg + 23]);
   if parallel_fil
       par_x1 = prev_x + floor(1.5*radius*cosd(rand_deg + 90));
       par_y1 = prev_y + floor(1.5*radius*sind(rand_deg + 90));
       S = (insertSphere(S, ss, radius, par_y1, par_x1, good_z));
       par_x2 = prev_x + floor(1.5*radius*cosd(rand_deg - 90));
       par_y2 = prev_y + floor(1.5*radius*sind(rand_deg - 90));
       S = (insertSphere(S, ss, radius, par_y2, par_x2, good_z));
   else
       S = (insertSphere(S, ss, radius,good_y(ii), good_x(ii), good_z));
   end
   for kk = 1:filament_length
       prev_x = prev_x + floor(radius*cosd(rand_deg));
       prev_y = prev_y + floor(radius*sind(rand_deg));
       if parallel_fil
           par_x1 = prev_x + floor(2*radius*cosd(rand_deg + 90));
           par_y1 = prev_y + floor(2*radius*sind(rand_deg + 90));
           S = (insertSphere(S, ss, radius, par_y1, par_x1, good_z));
           par_x2 = prev_x + floor(2*radius*cosd(rand_deg - 90));
           par_y2 = prev_y + floor(2*radius*sind(rand_deg - 90));
           S = (insertSphere(S, ss, radius, par_y2, par_x2, good_z));
       else
           S = (insertSphere(S, ss, radius, prev_y, prev_x, good_z));
       end
       rand_deg = randi([rand_deg - 22,rand_deg + 23]);% This should only allow next segment 45deg freedom
       
       
       
       
   end
end
S(S>1) = 1;
S = imgaussfilt3(S);
end

