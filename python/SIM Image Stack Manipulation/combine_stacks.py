import tifffile
import math
import os
import numpy as np
import cv2

# This code needs to be generalized to handle n angles. Currently it only works with 3 angles. 

def combine_stacks(path,filename, angles, phases, planes, scale_z, lateral_shrink):
    
    image_stack = tifffile.imread(path + filename)
    dims = image_stack.shape

    a = np.zeros((planes, angles, phases, 1, dims[1],dims[2])).astype(image_stack.dtype)
    # a1 = np.zeros(dims[0]//angles, dims[1],dims[2],planes).astype(image_stack.dtype)
    # a2 = np.zeros(dims[0]//angles, dims[1],dims[2],planes).astype(image_stack.dtype)
    stack_size = angles*phases

    #seperate angles

    for i in range(planes):
        for j in range(phases):
            for k in range(angles):
                a[i,k,j,0,:,:] = image_stack[j + phases*k,:,:]

        #angle_stk = np.reshape(angle_stk, ((dims[0]//stack_size)*phases ,dims[1] , dims[2], angles))   
        current_file = filename[0:-5] + str(i + 1) + '.tif'
        image_stack = tifffile.imread(path + current_file)
    print(a.shape)
    print(a[:,0,0,0,:,:].shape)

    scale_by = 3
    axial_scale = ''
    new_stacks = np.zeros((planes,angles,phases,1,dims[1],dims[2])).astype(image_stack.dtype)
    for j in range(phases):
        if scale_z:

            b = np.zeros((planes*scale_by, dims[1], dims[2])).astype(a.dtype)
            c = np.zeros((planes*scale_by, dims[1], dims[2])).astype(a.dtype)
            d = np.zeros((planes*scale_by, dims[1], dims[2])).astype(a.dtype) 
            for k in range(dims[1]):
                b[:,k,:] = cv2.resize(a[:,0,j,0,k,:], dsize=(dims[1], planes*scale_by), interpolation=cv2.INTER_CUBIC)
                c[:,k,:] = cv2.resize(a[:,1,j,0,k,:], dsize=(dims[1], planes*scale_by), interpolation=cv2.INTER_CUBIC)
                d[:,k,:] = cv2.resize(a[:,2,j,0,k,:], dsize=(dims[1], planes*scale_by), interpolation=cv2.INTER_CUBIC)

            if j == 0:   
                new_stacks = np.zeros((planes*scale_by,angles,phases,1,dims[1],dims[2])).astype(image_stack.dtype)
                axial_scale = 'axial_scale'
            
            new_stacks[:,0,j,0,:,:] = b[:,:,:]
            new_stacks[:,1,j,0,:,:] = c[:,:,:]
            new_stacks[:,2,j,0,:,:] = d[:,:,:]
                
        else:
            scale_by = 1
            new_stacks[:,0,j,:,:,:] = a[:,0,j,:,:,:]
            new_stacks[:,1,j,:,:,:] = a[:,1,j,:,:,:]
            new_stacks[:,2,j,:,:,:] = a[:,2,j,:,:,:]

        

    # tifffile.imwrite(path + "atest" + str(0) + '_' + str(j) + ".tif", new_stacks[:,0,6,0,:,:])
    # tifffile.imwrite(path + "atest" + str(1) + '_' + str(j) + ".tif", new_stacks[:,1,6,0,:,:])
    # tifffile.imwrite(path + "atest" + str(2) + '_' + str(j) + ".tif", new_stacks[:,2,6,0,:,:])

    #tifffile.imwrite(path + "atest" + str(0) + ".tif", new_stacks[:,0,0,0,:,:])
    print(new_stacks.shape)
    print(np.max(new_stacks))
    stack_size = phases*angles*planes*scale_by
    final_order = np.zeros((stack_size,dims[1], dims[2])).astype(a.dtype)
    count = 0
    for i in range(planes*scale_by):
        for j in range(angles):
            for k in range(phases):
                final_order[count,:,:] = new_stacks[i,j,k,0,:,:]
                count += 1

    
    print(np.max(final_order))
    
    scale = True
    dims = final_order.shape
    if scale:
        b = np.zeros((dims[0], np.ceil(dims[1]*lateral_shrink).astype(np.uint16), np.ceil(dims[2]*lateral_shrink).astype(np.uint16))).astype(a.dtype) 
        for i in range(dims[0]):
            b[i,:,:] = cv2.resize(final_order[i,:,:], dsize=(np.ceil(dims[1]*lateral_shrink).astype(np.uint16), np.ceil(dims[2]*lateral_shrink).astype(np.uint16)), interpolation=cv2.INTER_CUBIC)
            maxb = np.max(b)
        # b = ((b/maxb)*(2**16-1)).astype(np.uint16) #convert to 16 bit    
        # b = b.astype(np.uint16)
        print(np.max(b))   
        tifffile.imwrite(path + 'final_image_stack_half_lateral_' + axial_scale + '.tif', b)
    else:
        tifffile.imwrite(path + "final_image_stack_" + axial_scale + ".tif", final_order)

        


    