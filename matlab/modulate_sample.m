function sim_stack = modulate_sample(sample,modulation,otf,period,phases, grab_planes, z_step)
%modulate_sample This function creates a Z stack of the sample data
% translated axially through the sim modulation pattern and applies the
% defocus of the optical transfer function to the modulated stack

% sample is the source data
% modulation is the sim modulation pattern
% otf is the system optical transfer function or defocus function
% period is the sim patern latteral period for pattern translation
% phases is the number of sim pattern shifts
% grab planes is for generating multiple focus capture planes.

%%
    [x_dims, y_dims, z_dims] = size(sample);
    num_orientations = 1; % number of pattern orientations
    
    % rotational angles, since stripes start as vertical
    % add 90 degrees to these angles for equivalent xy orientation
    % -15 = 75, 45 = 135, 105 = 195 or 15 degrees
    angles=[-15 45 105]; 
    
    
    %%flags
    limit_z = 1; % to avoid generating an excessivly large stack for thin sample
    
    z_dims = x_dims*num_orientations*phases; 
    center = floor((x_dims/2) + 1);
    fprintf("center: %f\n", center);
    z_l = 1;
    z_r = x_dims;
    
    %% limit z_stack planes
    if limit_z
       if grab_planes
           %% grab 1/4 of the stack (513 - 1)
           z_dims_min = floor((x_dims-1)/(2*z_step))*num_orientations*phases;
           z_dims = floor(z_dims_min/(num_orientations*phases));
           z_l = center - floor(z_dims/2);
           z_r = center + floor(z_dims/2);          
       else
           z_dims_min = x_dims*num_orientations*phases/3;
           z_dims = floor(z_dims_min/(num_orientations*phases));
           z_l = center - floor(z_dims/2);
           z_r = center + floor(z_dims/2);
       end
    end

    
    %% For grabing multiple focus planes
    %% Plane spacing is always going to depend on your chosen voxel dimensions
    %% When you generate the PSF and Modulation Pattern
    if(grab_planes)
        slices = [245, 248, 251, 254, 257, 260, 263, 266, 269]; %108nm spacing
        %slices = [241, 245, 249, 253, 257, 261, 265, 269, 273]; %our 144nm spacing
        %slices = [233, 239, 245, 251, 257, 263, 269, 275, 281]; %216nm spacing
        %slices = [237, 242, 247, 252, 257, 262, 267, 272, 277]; %160nm spacing
        %slices = [201, 215, 229, 243, 257, 271, 285, 299, 313]; %504nm
        [r, c] = size(slices);    
        fprintf('c: %f, x_dims: %f, z_dims_min: %f\n',c, x_dims, z_dims_min);
        sim_stack = single(zeros(x_dims, x_dims, z_dims, c));
    else
        sim_stack = single(zeros(x_dims, x_dims, z_dims));
    end
    
    %gpu_count = gpuDeviceCount;
    %mod_shift = (zeros(x_dims, x_dims, x_dims));
    
    %% Uses CUDA acceleration for speeding up processing
    % a 513*513*513 3D volume needs ~8-10 GB of GPU memory
    if canUseGPU()
        gpuDevice(1);
        sample = single(gpuArray(sample));
        modulation = single(gpuArray(modulation));
        otf = single(gpuArray(otf));
        im2mod = single(gpuArray(zeros(x_dims, x_dims, x_dims)));
        mod_shift = single(gpuArray(zeros(x_dims, x_dims, x_dims)));
    else
        im2mod = (zeros(x_dims, x_dims, x_dims));
        mod_shift = (zeros(x_dims, x_dims, x_dims));
        fprintf("No GPU for YOU!!\n"); 
    end
    
    

    %% The Main Loop
    for ii = 1:num_orientations
%       fprintf("Outer Loop %d\n", ii);
        for jj = 1:phases
            %if (grab_planes)
            %    shift_amount = round(period * ((jj-1)*3*pi/7)/(2*pi));
            %else
            shift_amount = round(period * ((jj-1)/phases));
            fprintf('shift_amount: %d\n', shift_amount);
            mod_shift = gpuArray(cat(2, modulation(:,shift_amount + 1:end,:), modulation(:,1:shift_amount, :)));
            
            mod_shift = imrotate(mod_shift ,angles(ii),'bicubic','crop');
            % this is only for non rotated or shifted image
            %mod_shift = modulation;

            lead_og = z_l;
            stop_og = z_r;
              
            for kk = 1:z_dims
                tic; % show how much time it took to execute 1 loop
                lead_im = center - (z_r - z_l) + kk - 1;
                stop_im = center + kk - 1;
                
%                 fprintf("lead_im: %d, stop_im: %d\n", lead_im, stop_im);
%                 fprintf("lead_og: %d, stop_og: %d\n", lead_og, stop_og);
                
                if(grab_planes)
                    for pp = 1:c
                        %tic;
                        delta = center - slices(pp);
                        lead_og = center - (x_dims - 1)/8;
                        stop_og = center + (x_dims - 1)/8;
                        % shift sample through modulation field
                        lead_im = (x_dims - 1)/8 + 1 + (kk - 1)*z_step;
                        stop_im = (x_dims - 1)/4 + (x_dims - 1)/8 + 1 + (kk - 1)*z_step;
                        
                        im2mod = single(gpuArray(zeros(x_dims, x_dims, x_dims)));
                        im2mod(:,:,lead_im:stop_im) = sample(:,:,lead_og:stop_og); 
                        modulated = single(gpuArray((im2mod) .* mod_shift));

                        % for shifting to propper focal plane
                        if(delta > 0)
                          lead_im = delta + 1;
                          stop_im = x_dims;
                          lead_og = 1;
                          stop_og = x_dims - delta;
                        else
                          lead_im = 1;
                          stop_im = x_dims + delta; %negative delta
                          lead_og = -delta + 1; %sign makes it positive
                          stop_og = x_dims;
                        end
                        % Shift focal plane of interest to center to apply
                        % otf defocus
%                         fprintf("lead_im: %d, stop_im: %d, delta: %d\n", lead_im, stop_im, stop_im - lead_im);
%                         fprintf("lead_og: %d, stop_og: %d, delta: %d\n", lead_og, stop_og, stop_og - lead_og);
                        im2blur = single(gpuArray(zeros(x_dims, x_dims, x_dims)));
                        im2blur(:, :, lead_im:stop_im) = modulated(:, :, lead_og:stop_og);
                        clear modulated;
                        intIm = single(gpuArray(ifftn(fftn(im2blur).*fftshift(otf))));
                        sim_stack(:,:,(kk - 1)*phases*num_orientations + (jj) + (ii - 1)*phases, pp) = single(abs(intIm(:,:,center)));
                    end
                else
                    im2mod = (zeros(x_dims, x_dims, x_dims));
                    im2mod(:,:,lead_im:stop_im) = sample(:,:,lead_og:stop_og); 
                    modulated = gpuArray((im2mod) .* mod_shift);
                    %write3Dtiff(single((gather(modulated(:,:,:)))), 'modulated_sample.tif');
                    %write3Dtiff(single((gather(im2mod(:,:,:)))), 'unmodulated_sample.tif');
                    intIm = (ifftn(fftn(modulated).*fftshift(otf)));
                    %write3Dtiff(single(abs(gather(intIm))), 'blured_image.tif');
                    %write3Dtiff(single(abs(gather(intIm(:,:,center)))), 'plane_of_intrest.tif');
                    %break
                    sim_stack(:,:,(kk - 1)*phases*num_orientations + (jj) + (ii - 1)*phases) = single(abs(intIm(:,:,center)));
                end
                toc; % end of the time check
                %break
            end
            %break
        end
        %break
    end
    

end

