%Create a Filament Plane
clear all;

write_path = 'D:\SIM_Data\SIM Code Stuff\OTF_SIMULATION\bead volume\';
file_name = '50_Filaments_Pairs_In_Plane_r2_try.tif';

out_file = ([write_path, file_name]);

myFPlane = sphereFilamentPlane(50,20,513,2,1);
write3Dtiff(single(myFPlane(:,:,:)), out_file);