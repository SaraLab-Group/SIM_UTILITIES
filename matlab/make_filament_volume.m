%Create a Filament Plane
clear all;

write_path = '/run/media/imaging/RawHeadRex/SIM_Data/SIM Code Stuff/OTF_SIMULATION/Filaments/';
file_name = '200_Filament_pairs_In_volume_r1_morerandL.tif';

out_file = ([write_path, file_name]);

myFVol = filamentVolume(200,80,513,1,1, 2.5);
write3Dtiff(single(myFVol(:,:,:)), out_file);