import glob
import os
import fnmatch 
import numpy as np

rootdir = '/home/const/Dropbox/Documents/For DN/Very final/data'

mult = np.array([1, 1, 1, 1, 2, 2])

for dir, dirs, files in os.walk(rootdir):
    for filename in fnmatch.filter(files, 'hyper.dat'):
        fname = os.path.join(dir, filename)
        
        arr = np.loadtxt(fname, skiprows=1)
        sum = 0
        for i in arr[:,4:-1]:
            print(mult.shape, i.shape)
            sum += np.dot(mult, i)
            
        print(fname)
        print(sum)
            
            
            
        
        
