%% http://courses.washington.edu/matlab1/Matlab4BS_c5.htm

clear all
close all
%digitsOld = digits(100);
%pi100 = vpa(pi);%double(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679);

phase_shift = pi/2;
pix_val = 111;%ceil(phase_shift / (2*pi)*255);
pix_val2 = 133;

vert = 1152*4;
horz = 1920*4;
x=linspace(1, horz, horz); %setting horz to sensor size.
period = 10;
period2 = 6;
shift = 5;

start_deg = 15;
rot_step = 60; %in degrees
%% sf=period; % spatial freq in cycles per image

%% Fixed Period for the HF pattern
coswave2=(1/2*cos((x*2*pi)/period2 + shift*2*pi / 5) + .5);%*(2^16-1);%%
onematrix = ones(horz, horz);
coswave2Dx =(onematrix.*coswave2);
M2 = imrotate(coswave2Dx,-62, 'bilinear');
M2 = M2(4234:6153, 4234:6153);
[Y , Yinv] = ButterworthLPF2D(M2(1:vert/4, 1:horz/4), vert/(4*2.25), 56);
[fftLPF , fftHPF] = ButterworthLPF2D(M2(1:vert/4, 1:vert/4), 1920/10.5, 100);


%fftLPF = fftLPF.*(fftLPF > 10e-30);
%fftHPF = fftHPF.*(fftHPF > 10e-30);
Yshift = fftshift(Y);
fftHPFshift = fftshift(fftHPF);
fftLPFshift = fftshift(fftLPF);
figure(1);
subplot(1,2,1);
imshow(fftHPFshift);
axis equal;
axis tight;
colormap(gray);
subplot(1,2,2)
imshow(fftLPFshift);


colormap(gray);
axis equal;
axis tight;


figure(2);
imagesc(M2);
colormap(gray);

    coswave=(1/2*cos((x*2*pi)/period) + .5);
    
    %scaled_cos = coswave;%*(2^16-1);
    %figure;
    %plot(x(1:80), scaled_cos(1:80));
    
    %%close all
% for ii = 0:179
%     
%     coswave2D =uint8(((onematrix.*coswave)*pix_val));
%     file = sprintf("sineRotation_%dpixAngle_%03d.bmp",period ,ii);
%     %coswave2D = uint16(coswave2D);
%     %figure(2);
%     %colormap(gray);
%     M = imrotate(coswave2D,ii, 'bicubic');
%     [Y_dim, X_dim] = size(M)
%     YL = floor(Y_dim/2) - vert/8;
%     YR = floor(Y_dim/2) + vert/8 - 1;
%     XL = floor(X_dim/2) - horz/8;
%     XR = floor(X_dim/2) + horz/8 - 1;
%     
%     %vert_im = coswave2D(1:vert/4, 1:horz/4);
%     imwrite(M(YL:YR,XL:XR),file,'bmp');
% end

% for ii = 0:255
%     
%     coswave2D =uint8(((onematrix.*coswave)*ii));
%     file = sprintf("sineIntensity_%03d.bmp",ii);
%     %coswave2D = uint16(coswave2D);
%     %figure(2);
%     %colormap(gray);
% %     M = imrotate(coswave2D,ii, 'bicubic');
% %     [Y_dim, X_dim] = size(M)
% %     YL = floor(Y_dim/2) - vert/8;
% %     YR = floor(Y_dim/2) + vert/8 - 1;
% %     XL = floor(X_dim/2) - horz/8;
% %     XR = floor(X_dim/2) + horz/8 - 1;
%     
%     %vert_im = coswave2D(1:vert/4, 1:horz/4);
%     imwrite(coswave2D(1:vert/4, 1:horz/4),file,'bmp');
% end


for jj = 1:3
    
%     if jj == 2
%         ang_pix_val = pix_val2;
%     else
         ang_pix_val = pix_val;
%     end
    
    for ii = 1:shift
        coswave=(1/2*cos((x*2*pi)/period + (ii-1)*2*pi / 5) + .5);
        scaled_cos = coswave;%*(2^16-1);
        %figure;
        %plot(x(1:80), scaled_cos(1:80));

        %%close all

        coswave2D =uint8((onematrix.*scaled_cos)*ang_pix_val);
        
        rot_angle = start_deg + rot_step*(jj - 1);
        file = sprintf("%dpix_%d_gv_%ddeg_phase_%d.bmp", period, ang_pix_val, rot_angle, ii);
        M = imrotate(coswave2D,rot_angle, 'bicubic');
        [Y_dim, X_dim] = size(M);
        YL = floor(Y_dim/2) - vert/8;
        YR = floor(Y_dim/2) + vert/8 - 1;
        XL = floor(X_dim/2) - horz/8;
        XR = floor(X_dim/2) + horz/8 - 1;

        %vert_im = coswave2D(1:vert/4, 1:horz/4);
        imwrite(M(YL:YR,XL:XR),file,'bmp');
        

%         %coswave2D = uint16(coswave2D);
%         %figure(2);
%         %colormap(gray);
%         %vert_im = coswave2D(1:vert/4, 1:horz/4);
%         %imagesc(vert_im);
% 
% 
%         %min_15 = L(3903:5054,2956:4875);
%         %figure(3)
%         %colormap(gray)
%         J = imrotate(coswave2D,-60, 'bilinear');
%         min_60 = J(2956:4875,2956:4875);
%         %imagesc(min_60);
% 
%         %figure(4)
%         %colormap(gray)
%         K = imrotate(coswave2D,60, 'bilinear');
%         plus_60 = K(2956:4875,2956:4875);
% 
%         L = (imrotate(coswave2D,-20, 'bicubic'));
%         %min_15 = imresize((L(2956:4875,2956:4875).*M2).*Yshift,.5);
%         %min_15 = ((L(2956:4875,2956:4875).*M2).*Yshift);
%         %ffunonn1 = fft2(min_15);
%         min_15 = ((L(3263:4416, 2956:4875)));
%         %min_15 = real(ifft2(fftLPF.*ffunonn1));
% 
%         M = imrotate(coswave2D,-75);
%         min_75nn = M(2956:4875,2956:4875);
%         M = imrotate(coswave2D,-80, 'bicubic');
%         %min_75bl = imresize((M(2956:4875,2956:4875).*M2).*Yshift,.5);
%         %min_75bl = ((M(2956:4875,2956:4875).*M2).*Yshift);
%         min_75bl = ((M(3263:4416, 2956:4875)));
%         ffunonn2 = fft2(min_75bl);
%         %ffmult = ffunonn2.*fftLPF;
%         %min_75bl2 = real(ifft2(fftLPF.*ffunonn2));
%         M = imrotate(coswave2D,-75, 'bicubic');
%         min_75bc = M(2956:4875,2956:4875);
% 
%         N = imrotate(coswave2D,40, 'bicubic');
%         %plus_45 = imresize((N(2956:4875,2956:4875).*M2).*Yshift,.5);
%         %plus_45 = (N(2956:4875,2956:4875).*M2).*Yshift;
%         plus_45 = (N(3263:4416, 2956:4875));
%         ffunonn3 = fft2(plus_45);
%         %plus_45 = real(ifft2(fftLPF.*ffunonn3));
% 
%         %imagesc(plus_60);
%     %     imwrite(vert_im,file1,'tiff');
%     %     imwrite(min_60,file2,'tiff');
%     %     imwrite(plus_60,file3,'tiff');
%          imwrite(min_15,file4,'bmp');
%          imwrite(min_75bl,file5,'bmp');
%          imwrite(plus_45,file6,'bmp');
%          shift = shift + 1;
%          min_15 = double(min_15).*Yshift;
%          min_75bl = double(min_75bl).*Yshift;
%          plus_45 = double(plus_45).*Yshift;

    end
end

