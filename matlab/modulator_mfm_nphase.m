%% Sir Mix a'lot
%% For modulating MF-SIM Planes
clear all;
close all;

modPath = '/run/media/imaging/RawHeadRex/SIM_Data/SIM Code Stuff/OTF_SIMULATION/working data/';
modulation = 'modulation_field_488nm_pix_pitch_35nm_pix field_513_14pix_3D.tif';
otf = 'PSF RW_425nm_35nmVoxel_otf.tif';
% modulation = 'modulation_field_491nm_pix_pitch_32.25nm_pix field_513_1p42nimm_3D.tif';
% otf = 'PSF GL_515_1p43nimm_otf.tif';
samplePath = '/run/media/imaging/RawHeadRex/SIM_Data/SIM Code Stuff/OTF_SIMULATION/Single Bead/';
sampleName = 'single_centered_point.tif';

tic
data_out = [[samplePath,sampleName(1:end - 4),'_488nm_period_35nm_planes.tif']];


threeD_mod = single(loadtiff([modPath, modulation]));
sample = single(loadtiff([samplePath, sampleName]));
main_otf = single(loadtiff([modPath, otf])); 

grab_planes = 1;
z_step = 1;
angles = 3;
phases = 5;
period = 13.135135135/2;%12.8;%15.04;%

%% If grab_planes is true, writes all focal planes into w seperate files.
% if(grab_planes)
    simStack = modulate_mfm_sample(sample, threeD_mod, main_otf, period, phases, grab_planes, z_step,angles);
    [x,y,z,w] = size(simStack);
    for ii = 1:w
        write3Dtiff(single(abs(simStack(:,:,:,ii))), [data_out(1:end-4),'_plane',num2str(ii),'.tif']);
    end
% else
%     simStack = modulate_mfm_sample(sample, threeD_mod, main_otf, period, phases, 0, 1);
%     [x,y,z] = size(simStack);
%     write3Dtiff(single(abs(simStack(:,:,:))), data_out); 
% end
toc