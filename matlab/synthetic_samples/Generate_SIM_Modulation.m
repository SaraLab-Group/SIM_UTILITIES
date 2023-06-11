%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate a theoretical OTF for a structured illumination microscope %
% James Manton, 2019
% Modified by Antone Bajor 2022 to generate synthetic modulated psf for
% SIM Reconstruction
% Pupil calculations http://kmdouglass.github.io/posts/simple-pupil-function-calculations/
% 
% There are some issues generating an OTF/PSF with this utility as the
% Spherical shell generated is fairly rough.  Perhaps Smoothing the 
% Spherical shell will improve the OTF/PSF but I recommend Generating
% A PSF in imagej/fiji using the PSF Generator Pluggin then converting 
% PSF to OTF in Matlab or Python. I use the Richardson Wolf PSF and
% use isotropic voxel dimensions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Simulation parameters

%% Odd value gives a center point and center planes.
%% If Odd value isn't used modulation field will have a diagonal drift
field_size = 513                                                                                                                                                                                       ;
numerical_aperture_primary = 1.3;%1.3; %% Numerical Aperture of the Objective
numerical_aperture_detection = numerical_aperture_primary;
numerical_aperture_secondary = numerical_aperture_primary;
refractive_index = 1.4;%1.518;%1.4; %% Refractive Index of the immersion oil

INTERFEROMETRIC_DETECTION = 0;
SAVE_IMAGES = 1;
SAVE_3D = 0;
SAVE_OTF = 1;
SAVE_PSF = 1;
root_name = '3DSIM';

%% Values for our system

%% For Generating PSF and OTF use emmision wavelength
%% For Generating Modulation Field use excitation wavelength
wave_length = 488e-9;%525e-9;


slm_pixel_size = 9.2e-6; %% Structured Light Modulator Pixel size
pixels_per_period = 10; %% Period of pattern
f_tubelens = 180e-3;%180e-3; %% Focal Length of Tube Lens
mag_obj = 60; %% Magnification of Tube Lens
f_objective = f_tubelens/mag_obj; %% Focal Length of Objective

f_slm_lens = 610e-3;%610e-3;%500e-3;%610e-3; %% Focal Length of the SLM

rd = f_objective * numerical_aperture_primary; %% Pupil Radius

%% Difraction angle of 1st orders
theta_1st = asin(wave_length/(slm_pixel_size*pixels_per_period));
r_1st = tan(theta_1st)*f_slm_lens; %% Radius of 1st orders in pupil

%number_of_phases = 5;
%abbe = wave_length/(2*numerical_aperture_primary);
%% For best performance pixel pitch will need to be compared with the ability to laterally shift simulated
%% modulation period by desired amount, jumps in shift will cause funny seperation bands.
pix_pitch = 35e-9; %% Sampling of pixels (use at most 1/2 actual camera sensor for simulation)
wl_str = [num2str(wave_length*1e9),'nm_pix_pitch_',num2str(pix_pitch*1e9),'nm_pix field_',num2str(field_size)];
%abbe_pix = abbe/pix_pitch;

freq_NA = numerical_aperture_primary/wave_length; %High frequency supported by system 1/m
freq_NIMM = refractive_index/wave_length;
freq_samp = 1/pix_pitch; %sampling frequency
dFreq = freq_samp/(field_size); % 1/(m*pix)
fprintf('freq_NA: %f, freq_samp: %f, dFreq: %f\n', freq_NA, freq_samp, dFreq);
pupilRad = (freq_NA/dFreq); %pupil radius in pixels
sphereRad = freq_NIMM/dFreq;
fprintf('pupilRad: %f\n', pupilRad);
fprintf('sphereRad_NIMM: %f\n', sphereRad);
pupil_per_rad = pupilRad/((field_size)/2);% percent of 1/2 field radius
fprintf('pupil_per_rad: %f\n', pupil_per_rad);
prim_a_r = asin(numerical_aperture_primary/refractive_index); % half angle alpha
fprintf('half angle alpha: %f\n', prim_a_r*180/pi);
fprintf('should be same: %f\n', asin(pupilRad/sphereRad)*180/pi);
sphere_rad = pupilRad/sin(prim_a_r);
fprintf('sphere_rad: %f\n', sphere_rad);
sphere_per_rad = sphere_rad/((field_size - 1)/2); % Radius of sphere that defines spherical cap
fprintf('sphere_per_rad: %f\n', sphere_per_rad);


%%

%spot_radius = 0.95;

spot_radius = r_1st / rd; % 1st order spot radius in percentage of pupil radius
fprintf('spot_radius: %f\n', pupilRad*spot_radius);
fprintf('spot_radius per: %f\n', spot_radius);

%primary_alpha = asind(numerical_aperture_primary / refractive_index);

pupHeight = sphere_rad*cos(prim_a_r); %% Height to the bottom of spherical cap
fprintf('pupHeight: %f\n', pupHeight);
height_z = sphere_rad - pupHeight; %% Height of the spherical cap
fprintf('height_z: %f\n', height_z);
% secondary_alpha = asind(numerical_aperture_secondary / refractive_index);
maxkz_freq = height_z*dFreq;
fprintf('maxkz_freq: %f\n', maxkz_freq);
maxz_period = 1/maxkz_freq;
fprintf('maxz_period: %fum\n',maxz_period*1e6);

abbe_period = 2*wave_length/((rd/f_objective)^2);
fprintf('abbe_period: %fum\n', abbe_period*1e6);
lat_period = wave_length/(2*rd*spot_radius/f_objective);
fprintf('lat_period: %fnm\n', lat_period*1e6);



theta_spot = acos(pupilRad*spot_radius/sphere_rad);
fprintf('theta_spot deg: %f\n',theta_spot*180/pi);

spotH_kz = sphere_rad - tan(theta_spot)*pupilRad*spot_radius;
fprintf('spotH_kz: %f\n', spotH_kz);
fprintf('********** IMPORTATN **********\n'); % Meme game too stronk
f_kz = spotH_kz*dFreq;
fprintf('f_kz: %f, period_kz: %fnm\n', f_kz, 1/f_kz*1e9);

f_kxy = pupilRad*spot_radius*dFreq*2;
fprintf('f_kxy: %f, period_kxy: %fnm\n', f_kxy, 1/f_kxy*1e9);
fprintf('*******************************\n');

% beta = prim_a_r*spot_radius;%asin(rd/(f_objective*refractive_index));
% fprintf('beta: %f\n', beta);

kz_maybe = (1-cos(90 - theta_spot)*refractive_index)/wave_length;
fprintf('period_maybe: %fum\n', 1/(kz_maybe)*1e6);

%% Coordinates
lim = field_size - (field_size + 1) / 2;
x = single(-lim:lim);
y = single(-lim:lim);
z = single(-lim:lim);
[X, Y, Z] = ndgrid(x, y, z);
R = single(sqrt(X.^2 + Y.^2 + Z.^2));
phi = single(atan2(X, Y) * 180 / pi);
theta = single(atan2(X, Z) * 180 / pi);
%R = 2 * R / max(x); %% Original code coordinates scaled 100% radius to 50%
%% this gives a 3D cylindrical column representing percentage of radius
R = R / max(x);
r = sqrt(X.^2 + Y.^2);
%r = 2 * r / max(x); %% Original code coordinates scaled 100% radius to 50%
%% This gives a radius of the field as a percentage
r = r / max(x);
%write3Dtiff(single(Z),['Z','_3D.tif']);
% write3Dtiff(single(R),['R','_3D.tif']);
% write3Dtiff(single(r),['r_lil','_3D.tif']);

% Detection OTF
detection_pupil = single(zeros(field_size, field_size, field_size));
%detection_pupil = complex(detection_pupil_real,0);

detection_pupil(r < pupil_per_rad) = 1; %% generates the detection pupil
write3Dtiff(single(detection_pupil),['detection_pupil_', wl_str, num2str(pixels_per_period),'pix_3D.tif']);
detection_ctf = single(zeros(field_size, field_size, field_size,'like', 1j));

%% The Else Statement generates the sphereical cap, not using interferometric detection
if INTERFEROMETRIC_DETECTION
    detection_ctf(R < (sphere_per_rad) & R > (sphere_per_rad - 0.01) & detection_pupil > 0) = 1 + 0i;
%     detection_ctf(R < 0.5 & R > 0.49 & detection_pupil > 0) = 1 + 0i;
else
    detection_ctf(R < (sphere_per_rad + 0.01) & R > (sphere_per_rad) & detection_pupil > 0 & Z < 0) = single(1 + 0i);
%     detection_ctf(R < 0.5 & R > 0.49 & detection_pupil > 0 & Z > 0) = 1 + 0i;
end
%detection_ctf = imgaussfilt3(detection_ctf);
write3Dtiff(single(detection_ctf),['detection_ctf', wl_str, num2str(pixels_per_period),'pix_3D.tif']);

%% This generates the OTF auto corelation
detection_otf = single(fftshift(ifftn(fftn(detection_ctf) .* conj(fftn(detection_ctf)))));

% detection_atf = (ifftn(fftn(detection_ctf)));  % looks like this generates a spherical shell.
% write3Dtiff(single(abs(detection_atf)), 'detection_atf.tif');

if (SAVE_OTF)
    detection_otf(detection_otf < 0) = 0;
    write3Dtiff(single(detection_otf),['detection', '_otf_','float_', wl_str, num2str(pixels_per_period),'pix_3D.tif'])
end

fprintf("After Detection OTF\n");
%% Primary OTF

primary_pupil = zeros(field_size, field_size, field_size);

%% For some reason this part of the code uses the 1/2 (Field - 1) value for radius instead of r array
%% This first part of the code generates the 1st order spots extruded through z axis
primary_pupil(floor(abs(Y)) == round(max(x) * spot_radius * pupil_per_rad) & floor(X) == 0) =   single(1);

%% This second part generates the 0th order spot extruded thorough z axis;
primary_pupil(floor(Y) == 0 & floor(X) == 0) = 1;

write3Dtiff(single(primary_pupil),['primary_pupil','_3D.tif']);

% prim_pupX = zeros(field_size,field_size,field_size);
% prim_pupX(:,:,256:258) = primary_pupil(:,:,256:258);

%write3Dtiff(single(prim_pupX),['prim_pupX','_3D.tif']);

% special_psf = otf2psf(prim_pupX);
% write3Dtiff(single(abs(special_psf)),['special_psf','_3D.tif']);

primary_ctf = single(zeros(field_size, field_size, field_size));

%% This takes the 3 beam orders "primary_pupil" and places them only where they intersect
%% Along the spherical shell using the pupil radius is incorrect since it defines a smaller
%% radius sphere, which gives incorrect axial modulation.
primary_ctf(R < (sphere_per_rad + 0.01) & R > (sphere_per_rad) & primary_pupil > 0 & Z > 0) = single(1 + 0i);

%primary_ctf(primary_pupil > 0 & Z > 0) =   1; %Attempt to simplify OTF
%write3Dtiff(single(primary_ctf),['primary_Pre_repmat','_3D.tif']);

% Fix band weights
primary_ctf_sum = single(repmat(sum(primary_ctf, 3), 1, 1, field_size)); %% not sure this gives different value than primary_ctf in my application
write3Dtiff(single(primary_ctf_sum),['primary_repmat','_3D.tif']);
primary_ctf(primary_ctf > 0) = primary_ctf(primary_ctf > 0) ./ primary_ctf_sum(primary_ctf > 0);

write3Dtiff(single(primary_ctf),['primary_ctf_',wl_str,'_3D.tif']);

%% This step generates the 7 beam spots that can been seen in the xz center plane
%% This defines the 3D sim intensity modulation field
primary_otf = single(fftshift(ifftn(fftn(primary_ctf) .* conj(fftn(primary_ctf)))));
write3Dtiff(single(primary_otf),['modulation_otf_',wl_str,'_3D.tif']);

clear primary_pupil;
clear detection_pupil;



%% Overall OTF
%%clear primary_ctf;



%% Converts detection otf to psf
detection_psf = single(otf2psf(detection_otf));

if (SAVE_PSF)
   write3Dtiff(single(abs(detection_psf)),['Detection_psf_', wl_str,'float', '_3D.tif'])
end


%clear detection_otf;
%% Creates the realspace sim Modulation field
modulation_field = single(otf2psf(primary_otf));
%% Scales max intensity to 1, and removes checker board patterning taking abs value
modulation_field = abs(modulation_field./max(modulation_field,[],'all'));

min_modulation = min(modulation_field, [], 'all');
fprintf("min_modulation: %f\n", min_modulation); %% if this is 1 something went wrong
write3Dtiff(single((modulation_field)),['modulation_field_', wl_str,'_' ,num2str(pixels_per_period),'pix_3D.tif']);
%clear all;

