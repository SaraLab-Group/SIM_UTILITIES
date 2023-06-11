%% TEST Spheres

%% Make Concentric Shells

volume = 513;
radius = 128;
Sfin = zeros(volume, volume, volume);
SCent = zeros(volume, volume, volume);
spacing = 2;
for ii = 1:floor(radius/spacing)
   S = tiffSphere(volume,radius - (ii - 1)*spacing, 0 , 0 , 0);
   Sfin = Sfin + (-1)^(ii - 1)*S;
end
Sfin(Sfin < 0) = 0;
write3Dtiff(single((Sfin)), ['conc_r',num2str(radius),'_s', num2str(spacing), '.tif']);
SCent(:,:,257) = Sfin(:,:,257);

write3Dtiff(single(abs(SCent)), ['conc_r',num2str(radius),'_s', num2str(spacing), '_center_plane.tif']);


%% make multiple small spheres (beads lawn)

% S1 = tiffSphere(513, 4, 0,0,0);
% S1 = S1 + tiffSphere(513, 4, 9,9,0);
% S1 = S1 + tiffSphere(513, 4, -9,9,0);
% S1 = S1 + tiffSphere(513, 4, 9,-9,0);
% S1 = S1 + tiffSphere(513, 4, -19,-20,0);
% S1 = S1 + tiffSphere(513, 4, -9,0,0);
% S1 = S1 + tiffSphere(513, 4, -9,-9,0);
% S1 = S1 + tiffSphere(513, 4, 20,8,0);
% S1 = S1 + tiffSphere(513, 4, -20,0,0);
% S1 = S1 + tiffSphere(513, 4, 20,-20,0);
% S1 = S1 + tiffSphere(513, 4, 19,22,0);
% S1 = S1 + tiffSphere(513, 4, 19,0,0);
% S1 = S1 + tiffSphere(513, 4, 19,31,0);
% S1 = S1 + tiffSphere(513, 4, -19,24,0);
% S1 = S1 + tiffSphere(513, 4, -19,15,0);
% S1 = S1 + tiffSphere(513, 4, -19,-10,0);
% S1 = S1 + tiffSphere(513, 4, 10,0,0);
% S1 = S1 + tiffSphere(513, 4, 0,-10,0);
% S1 = S1 + tiffSphere(513, 4, 0,11,0);
% S1 = S1 + tiffSphere(513, 4, 32,0,0);
% S1 = S1 + tiffSphere(513, 4, -33,0,0);
% S1 = S1 + tiffSphere(513, 4, -9,19,0);
% S1 = S1 + tiffSphere(513, 4, 9,21,0);
% S1 = S1 + tiffSphere(513, 4, -9,-21,0);
% S1 = S1 + tiffSphere(513, 4, 9,-19,0);
% S1 = S1 + tiffSphere(513, 4, -26,19,0);
% S1 = S1 + tiffSphere(513, 4, 20,-10,0);
% S1 = S1 + tiffSphere(513, 4, -29,-10,0);
% S1 = S1 + tiffSphere(513, 4, -26,9,0);
% S1 = S1 + tiffSphere(513, 4, 0,-19,0);
% S1 = S1 + tiffSphere(513, 4, -38,-10,0);
% S1 = S1 + tiffSphere(513, 4, -28,-18,0);
% S1 = S1 + tiffSphere(513, 4, -35,-24,0);
% S1 = S1 + tiffSphere(513, 4, -39,-32,0);
% S1 = S1 + tiffSphere(513, 4, -45,-24,0);
% S1 = S1 + tiffSphere(513, 4, -49,-12,0);
% S1 = S1 + tiffSphere(513, 4, 49,12,0);
% S1 = S1 + tiffSphere(513, 4, 35,24,0);
% S1 = S1 + tiffSphere(513, 4, 39,32,0);
% S1 = S1 + tiffSphere(513, 4, 45,24,0);
% S1 = S1 + tiffSphere(513, 4, 55,40,0);
% S1 = S1 + tiffSphere(513, 4, 48,34,0);
% S1 = S1 + tiffSphere(513, 4, -65,-50,0);
% S1 = S1 + tiffSphere(513, 4, 48,34,0);
% S1 = S1 + tiffSphere(513, 4, -75,-60,0);
% 
% 
% write3Dtiff(single(S1), "spheres in plane.tif");