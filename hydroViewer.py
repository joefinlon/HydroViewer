########## HYDROMETEOR VIEWER ##########
# IMPORT PACKAGES
import sys, os, glob
import xarray as xr
import numpy as np
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
from matplotlib import colors 

'''
Main script that calls subroutines to plot particles meeting user-specified criteria.

Parameters
----------
imgFile: str
    Path to decompressed data generated from read_binary_* script in UIOOPS.
plotDirectory: str
    Path to save image strips to file.
campaign: str
    Name of project (e.g., 'olympex'). Allows for project-specific conditional statements to be added to plotting routines.
probeName: str
    Probe type ('Fast2DC','2DS', 'HVPS', 'CIP', 'PIP'). Used to determine decryption routines specific to the manufacturer.
annotate: int
    0 [default] - Do not annotate every 5 particles in record; 1 - Annotate
chunkSize: int
    Number of image records (frames) to plot in succession
chunkNum: int
    nth chunk of frames to process for current job (e.g., chunkNum=2 plots frames chunkSize+1 to 2*chunkSize)
    '''
imgFile = sys.argv[1]
plotDirectory = sys.argv[2]
campaign = sys.argv[3]
probeName = sys.argv[4]
annotate = int(sys.argv[5])
chunkSize = int(sys.argv[6])
chunkNum = int(sys.argv[7])

# IMAGE BUFFER & PLOTTING FUNCTIONS
def probe_defaults(probeName):
    invalidSlice = np.array([-1., -1., -1., -1., -1., -1., -1., -1.])
    if probeName=='2DS' or probeName=='HVPS':
        boundary = np.array([43690, 43690, 43690, 43690, 43690, 43690, 43690, 43690])
        boundaryTime = 0
        bufferShape = [1700, 128]
    elif probeName=='CIP' or probeName=='PIP':
        boundary = np.array([170, 170, 170, 170, 170, 170, 170, 170])
        boundaryTime = 0
        bufferShape = [1700, 64]
    elif probeName=='Fast2DC':
        boundary = np.array([170, 170, 170])
        boundaryTime = 0
        bufferShape = [512, 64]
        
    return boundary, boundaryTime, invalidSlice, bufferShape

def get_slice_endpoints(probeName, buf, boundary, boundaryTime, invalidSlice): # get indices for bounaries & start/end of particle
    numPart = 0 # number of particles in buffer
    startInd = [] # index of the start of a particle
    endInd = [] # index of the end of a particle
    boundaryInd = [] # index of the particle boundary (series of 8 consecutive '43690' values)
    
    j = 0
    while (buf[j,0] != -1) and (j+1 < buf.shape[0]):
        if (np.array_equal(buf[j,0:len(boundary)],boundary)) and ((buf[j+1,0]==boundaryTime) or (probeName=='CIP') or 
                                                                  (probeName=='PIP') or (probeName=='Fast2DC')): # boundary
            boundaryInd.append(j)
            if j>0:
                endInd.append(j-1) # index of particle end before the particle boundary (previous particle)
                numPart = numPart + 1 # for particle preceeding boundary
            startInd.append(j+2) # index of particle start after the particle boundary (next particle)
        j = j + 1
    
    boundaryInd = np.array(boundaryInd, dtype='int')
    if boundaryInd.size==0: # no particles in buffer
        pass
    else:
        if (boundaryInd.size>0) or (boundaryInd[0]>0): # first boundary after first slice in buffer
            startInd = np.insert(startInd, 0, 0)
        startInd = np.array(startInd, dtype='int')
    
        if boundaryInd[-1]+3 < buf.shape[0]: # particle occurs after last boundary if it's 3rd from last slice in record
            endInd.append(buf.shape[0]-1) # last particle ends at end of buffer
            endInd = np.array(endInd, dtype='int') 
        else:
            endInd = np.array(endInd, dtype='int') # boundary found at very end of record -- no more particles in buffer
        
    return numPart, boundaryInd, startInd, endInd
    
def buffer_integrity(partCount, boundaryInd, partStart, partEnd):
    print('There are {} particles and {} boundaries.'.format(partCount, len(boundaryInd)))
    if np.any(partEnd-partStart<0): # particle end indices NOT >= start indices [ERROR]
        print('Error with particle start/end indices!')
    else:
        print('Particle start/end indices look OK.')
        
def image_buffer(buf, bufferShape, probeName, boundaryInd): # generate matrix of 1's and 0's from buffer
    if probeName=='2DS' or probeName=='HVPS':
        boundaryData = np.tile([1,2,2,1], 32) # alternate 1's and 2's for boundary slice (white & cyan pixels)
        buf[buf==-1] = 0 # change invalid values to 0 (unshadowed segment)
        buf = 65535 - buf # 0: shadowed; 1: unshadowed
        
        # convert decimal to binary (8 image words for each slice)
        imageData = np.ones(bufferShape).astype(int) # set up image buffer (1's mean unshadowed pixels)

        for x in np.arange(buf.shape[0]):
            tempBuf = np.array([np.binary_repr(int(buf[x,0]),16), np.binary_repr(int(buf[x,1]),16),
                                np.binary_repr(int(buf[x,2]),16), np.binary_repr(int(buf[x,3]),16),
                                np.binary_repr(int(buf[x,4]),16), np.binary_repr(int(buf[x,5]),16),
                                np.binary_repr(int(buf[x,6]),16), np.binary_repr(int(buf[x,7]),16)])
            sliceBuf = []
            for y in np.arange((buf.shape[1])*16):
                sliceBuf.append(tempBuf[(np.floor(y/16)).astype(int)][np.mod(y,16)])
            sliceBuf = np.asarray(sliceBuf, dtype='int')
            imageData[x,:] = sliceBuf
    elif probeName=='CIP' or probeName=='PIP' or probeName=='Fast2DC':
        boundaryData = np.tile([1,2,2,1], 16) # alternate 1's and 2's for boundary slice (white & cyan pixels)
        buf[buf==-1] = 255 # change invalid values to 0 (unshadowed segment)
        
        # convert decimal to binary (8 image words for each slice)
        imageData = np.ones(bufferShape).astype(int) # set up image buffer (1's mean unshadowed pixels)

        for x in np.arange(buf.shape[0]):
            tempBuf = np.array([np.binary_repr(int(buf[x,0]),8), np.binary_repr(int(buf[x,1]),8),
                                np.binary_repr(int(buf[x,2]),8), np.binary_repr(int(buf[x,3]),8),
                                np.binary_repr(int(buf[x,4]),8), np.binary_repr(int(buf[x,5]),8),
                                np.binary_repr(int(buf[x,6]),8), np.binary_repr(int(buf[x,7]),8)])
            sliceBuf = []
            for y in np.arange((buf.shape[1])*8):
                sliceBuf.append(tempBuf[(np.floor(y/8)).astype(int)][np.mod(y,8)])
            sliceBuf = np.asarray(sliceBuf, dtype='int');
            imageData[x,:] = sliceBuf
            #print(imageData)
    if boundaryInd.size>0:
        imageData[boundaryInd,:] = boundaryData # write in boundary slice
        if boundaryInd[-1]+1 < buf.shape[0]:
            imageData[boundaryInd+1,:] = 1
        else:
            imageData[boundaryInd[0:-2]+1,:] = 1

    return(imageData)

def annotate_particle_incriments(boundaryInd): # get indices along buffer for every 5 particles (optional plotting)
    if boundaryInd.shape>0:
        partInds = boundaryInd[np.arange(3,(boundaryInd.size)-1, 5)]+1
    else:
        partInds = 1
    
    return partInds
#
# ==================================================
# MAIN SCRIPT
# ==================================================
cmap = colors.ListedColormap(['black', 'white', 'cyan'])
bounds = [0, 1, 2, 3]
norm = colors.BoundaryNorm(bounds, cmap.N)

[boundary, boundaryTime, invalidSlice, bufferShape] = probe_defaults(probeName)

ds = xr.open_dataset(imgFile)
numFrames = ds.dims['time']
if chunkSize*chunkNum > numFrames:
    frameStart = np.arange(chunkSize*(chunkNum-1)+1,numFrames-4,6)
else:
    frameStart = np.arange(chunkSize*(chunkNum-1)+1,chunkSize*chunkNum-4,6)
numFrames = 6*len(frameStart)
print('Now starting iteration #{}.'.format(chunkNum))

# LOOP THROUGH FRAMES FOR PLOTTING
for iter in range(len(frameStart)):
    yr = ds['year'][frameStart[iter]-1:frameStart[iter]+5].values
    mon = ds['month'][frameStart[iter]-1:frameStart[iter]+5].values
    day = ds['day'][frameStart[iter]-1:frameStart[iter]+5].values
    hr = ds['hour'][frameStart[iter]-1:frameStart[iter]+5].values
    minute = ds['minute'][frameStart[iter]-1:frameStart[iter]+5].values
    sec = ds['second'][frameStart[iter]-1:frameStart[iter]+5].values
    msec = ds['millisec'][frameStart[iter]-1:frameStart[iter]+5].values
    data = ds['data'][frameStart[iter]-1:frameStart[iter]+5].values
    if probeName=='2DS':
    	try:
    		tasRatio = ds['tasRatio'][frameStart[iter]-1:frameStart[iter]+5].values
    	except:
    		tasRatio = np.ones(6)
    else:
        tasRatio = np.ones(6)

    [partCount1, boundaryInd, partStart, partEnd] = get_slice_endpoints(
        probeName, data[0], boundary, boundaryTime, invalidSlice)
    img1 = image_buffer(data[0], bufferShape, probeName, boundaryInd)
    if annotate==1:
        partInds1 = annotate_particle_incriments(boundaryInd)

    [partCount2, boundaryInd, partStart, partEnd] = get_slice_endpoints(
        probeName, data[1], boundary, boundaryTime, invalidSlice)
    img2 = image_buffer(data[1], bufferShape, probeName, boundaryInd)
    if annotate==1:
        partInds2 = annotate_particle_incriments(boundaryInd)

    [partCount3, boundaryInd, partStart, partEnd] = get_slice_endpoints(
        probeName, data[2], boundary, boundaryTime, invalidSlice)
    img3 = image_buffer(data[2], bufferShape, probeName, boundaryInd)
    if annotate==1:
        partInds3 = annotate_particle_incriments(boundaryInd)

    [partCount4, boundaryInd, partStart, partEnd] = get_slice_endpoints(
        probeName, data[3], boundary, boundaryTime, invalidSlice)
    img4 = image_buffer(data[3], bufferShape, probeName, boundaryInd)
    if annotate==1:
        partInds4 = annotate_particle_incriments(boundaryInd)

    [partCount5, boundaryInd, partStart, partEnd] = get_slice_endpoints(
        probeName, data[4], boundary, boundaryTime, invalidSlice)
    img5 = image_buffer(data[4], bufferShape, probeName, boundaryInd)
    if annotate==1:
        partInds5 = annotate_particle_incriments(boundaryInd)

    [partCount6, boundaryInd, partStart, partEnd] = get_slice_endpoints(
        probeName, data[5], boundary, boundaryTime, invalidSlice)
    img6 = image_buffer(data[5], bufferShape, probeName, boundaryInd)
    if annotate==1:
        partInds6 = annotate_particle_incriments(boundaryInd)

    if probeName=='Fast2DC':
        fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, sharex=True, sharey=True, constrained_layout=True)
    else:
        fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, sharex=True, sharey=True, constrained_layout=True)
    #imgExtent = [0, img1.shape[0], 0, img1.shape[1]]

    ax1.imshow((img1[0:img1.shape[0]+1,:]).T, cmap=cmap, norm=norm, extent=[0, img1.shape[0], 0, tasRatio[0]*img1.shape[1]],
               interpolation='none')
    ax1.set_aspect('equal')
    fileStr = 'fr{:06d}_{:02d}{:02d}{:02d}'.format(int(frameStart[iter]), hr[0], minute[0], sec[0])
    titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(
        int(frameStart[iter]), hr[0], minute[0], sec[0], msec[0], partCount1)
    if annotate==1:
        for x in range(len(partInds1)):
            ax1.text(partInds1[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
    ax1.set_title(titleStr, size=10)
    ax1.set_xticks([]), ax1.set_yticks([]), ax1.set_adjustable('box')

    ax2.imshow((img2[0:img2.shape[0]+1,:]).T, cmap=cmap, norm=norm, extent=[0, img2.shape[0], 0, tasRatio[1]*img2.shape[1]],
               interpolation='none')
    ax2.set_aspect('equal')
    titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(
        int(frameStart[iter])+1, hr[1], minute[1], sec[1], msec[1], partCount2)
    if annotate==1:
        for x in range(len(partInds2)):
            ax2.text(partInds2[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
    ax2.set_title(titleStr, size=10)
    ax2.set_xticks([]), ax2.set_yticks([]), ax2.set_adjustable('box')

    ax3.imshow((img3[0:img3.shape[0]+1,:]).T, cmap=cmap, norm=norm, extent=[0, img3.shape[0], 0, tasRatio[2]*img3.shape[1]],
               interpolation='none')
    ax3.set_aspect('equal')
    titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(
        int(frameStart[iter])+2, hr[2], minute[2], sec[2], msec[2], partCount3)
    if annotate==1:
        for x in range(len(partInds3)):
            ax3.text(partInds3[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
    ax3.set_title(titleStr, size=10)
    ax3.set_xticks([]), ax3.set_yticks([]), ax3.set_adjustable('box')

    ax4.imshow((img4[0:img4.shape[0]+1,:]).T, cmap=cmap, norm=norm, extent=[0, img4.shape[0], 0, tasRatio[3]*img4.shape[1]],
               interpolation='none')
    ax4.set_aspect('equal')
    titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(
        int(frameStart[iter])+3, hr[3], minute[3], sec[3], msec[3], partCount4)
    if annotate==1:
        for x in range(len(partInds4)):
            ax4.text(partInds4[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
    ax4.set_title(titleStr, size=10)
    ax4.set_xticks([]), ax4.set_yticks([]), ax4.set_adjustable('box')

    ax5.imshow((img5[0:img5.shape[0]+1,:]).T, cmap=cmap, norm=norm, extent=[0, img5.shape[0], 0, tasRatio[4]*img5.shape[1]],
               interpolation='none')
    ax5.set_aspect('equal')
    titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(
        int(frameStart[iter])+4, hr[4], minute[4], sec[4], msec[4], partCount5)
    if annotate==1:
        for x in range(len(partInds5)):
            ax5.text(partInds5[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
    ax5.set_title(titleStr, size=10)
    ax5.set_xticks([]), ax5.set_yticks([]), ax5.set_adjustable('box')

    ax6.imshow((img6[0:img6.shape[0]+1,:]).T, cmap=cmap, norm=norm, extent=[0, img6.shape[0], 0, tasRatio[5]*img6.shape[1]],
               interpolation='none')
    ax6.set_aspect('equal')
    titleStr = 'Frame #{}: {:02d}{:02d}{:02d}.{:03d} | {} particles'.format(
        int(frameStart[iter])+5, hr[5], minute[5], sec[5], msec[5], partCount6)
    if annotate==1:
        for x in range(len(partInds6)):
            ax6.text(partInds6[x], 3, '*', color='red', horizontalalignment='left', verticalalignment='top')
    ax6.set_title(titleStr, size=10)
    ax6.set_xticks([]), ax6.set_yticks([]), ax6.set_adjustable('box')

    if np.remainder(iter+1,500)==0:
        print('Finished plotting frame # {} / {}'.format(iter+1, numFrames/6))

    #fig.set_tight_layout(True)
    
    outFile = '{}img{}.{}.pdf'.format(plotDirectory, probeName, fileStr)
    plt.savefig(outFile, bbox_inches='tight', pad_inches = 0, format='pdf')
    plt.clf(), plt.cla(), plt.close('all') # close figure so that it doesn't get displayed