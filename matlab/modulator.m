%% Sir Mix a'lot
%% this code will modulate the sample data by chosen modulation data
%  then convert to Fourier space and multiply by the provided OTF function,
%  then collect the data plane by plane until samples is axially translated
%  through specified range. Then write the collected modulated defocused
%  image stack to disk.

clear all;
close all;

% With Great power comes great responsibility.

modPath = 'D:\SIM_Data\SIM Code Stuff\OTF_SIMULATION\working data\';
modulation = 'modulation_field_491nm_pix_pitch_23.78nm_pix field_513_10pix_3D.tif';
otf = 'PSF RW 491nm 23p78nm_otf.tif';

samplePath = 'D:\SIM_Data\SIM Code Stuff\OTF_SIMULATION\Single Bead\';
sampleName = 'single_centered_point.tif';
tic
data_out = [samplePath,sampleName(1:end - 4),'_525nm_SIMStack_psf_5_phase_491nm_nm_23.78nm vox.tif'];

threeD_mod = single(loadtiff([modPath, modulation]));
sample = single(loadtiff([samplePath, sampleName]));
main_otf = single(loadtiff([modPath, otf])); 

grab_planes = 0;
z_step = 1;

phases = 5;
% The modulation period in pixels Becareful ideally ones period should be
% evenly divisible by the phases, otherwise separation will be sub-optimal.
% make sure to verify period in pixels and rescale modulation pattern to
% something suitable by changing pixel pitch.
period = 20.5;%12.8;%15.04;

%% If grab_planes is true, writes all focal planes into w seperate files.
if(grab_planes)
    simStack = modulate_sample(sample, threeD_mod, main_otf, period, phases, 1, z_step);
    [x,y,z,w] = size(simStack);
    for ii = 1:w
        write3Dtiff(single(abs(simStack(:,:,:,ii))), [data_out(1:end-4),'_plane',num2str(ii),'.tif']);
    end
else
    simStack = modulate_sample(sample, threeD_mod, main_otf, period, phases, 0, 1);
    [x,y,z] = size(simStack);
    write3Dtiff(single(abs(simStack(:,:,:))), data_out); 
end
toc