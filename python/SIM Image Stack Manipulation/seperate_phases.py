import math
import os
import numpy as np
import seperate_phases

def seperate_phases(stack,phases,angles,path,scale_z):
    dims = stack.shape

    a = (np.zeros((dims[0]//phases,dims[1], dims[2], angles))).astype(stack.dtype)
    #print(a.shape)
    phase = 0

    while(phase < phases):
        # outfile = outFolder +filename[0:lend(filename) - 4] + 'phase' + str(phase) + '.tif'

        for i in range(dims[0]//phases):
            a[i,:,:] = stack[i*phases + phase,:,:]

        phase += 1