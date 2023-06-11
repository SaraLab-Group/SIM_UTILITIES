%% just gib me psf to otf
%% just a dirty little script to convert PSF into an OTF
%% using standard Matlab Function
clear all;
close all;

psfPath = '<Your Path Here>';
psfName = '<Your PSF File Name Here>.tif';


data_out = [[psfPath,psfName(1:end - 4),'_otf.tif']];

psf = loadtiff([psfPath, psfName]);

otf = fftshift(psf2otf(psf));

write3Dtiff(single(abs(otf(:,:,:))), data_out); 