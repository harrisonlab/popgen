#!/usr/bin/env python2.7
#source hdf5py-2.5.0
#source numpy-1.9.2

import sys
import string
import re
import numpy 
import h5py


bas5=(sys.argv[1]) # the file to be edited


with h5py.File('{}'.format(bas5), 'a') as g:   #Open the bas.h5 file in append mode
    for i in range(g['/PulseData/Regions'].shape[0]):   # iterate over all the rows in the /PulseData/Regions part of the file
        if g['/PulseData/Regions'][i,1]==1:  # Check whether you are in the area that has been identified to be between the bell adapters by the software
            if int(g['/PulseData/Regions'][i,3])-int(g['/PulseData/Regions'][i,2])<500: # Ignore these regions if the length of the read is less than 500 
                continue
            else:
                g['/PulseData/Regions'][i,3]=int(g['/PulseData/Regions'][i,3])-64  # Otherwise, trim the upper end by 64 bp
                if g['/PulseData/Regions'][i,2]!=0:
                    g['/PulseData/Regions'][i,2]=int(g['/PulseData/Regions'][i,2])+64 # trim the lower end by 64 pb